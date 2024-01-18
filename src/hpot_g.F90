!===============================================================
!
! Copyright (C) 2005 Finite Difference Research Group
! This file is part of parsec, http://www.ices.utexas.edu/parsec/
!
! Obtains the hartree potential (vhart) on the real space
! mesh by fourier transform of the explicit hartree term 
! in G-space
!
! For compatibility the resulting potential is returned 
! via the 'rho' input parameter
!
! Written in 2015-03 by AJB, relying heavily on gw_write and ionpbc/forpbc
!
! ?This subroutine is partially threaded, but only the master PE does the transform
!
! Looking at gw_write, it looks like the charge density is best passed outside of 
! elec_st because of a bug with the "IBM compiler", I think this is not important right now.
!---------------------------------------------------------------
subroutine hpot_G(elec_st, grid, pbc, parallel, my_rho, my_vhart)

!parallel,grid,elec_st,pbc,rho

  use constants
  use electronic_struct_module
  use grid_module
  use pbc_module
  use potential_module
  use parallel_data_module
  implicit none
  !
  ! Input/Output variables:
  !
  !  electronic structure
  type (electronic_struct), intent(in) :: elec_st
  !  grid related data
  type (grid_data), intent(in) :: grid
  !  periodic boundary conditions data
  type (pbc_data), intent(in) :: pbc
  !  parallel computation related data
  !  has to be inout because we use the global funciton storage space oh bother
  type (parallel_data), intent(inout) :: parallel
  !
  real(dp), intent(in) :: my_rho(parallel%mydim,2*elec_st%nspin-1)
  ! The resulting hartree potential on my chunk
  real(dp), intent(out) :: my_vhart(parallel%ldn) !watch it! ldn and not mydim
  !
  ! Work variables:
  !
  ! the hartree energy contribution
  real(dp) :: ehart_g
  logical :: report_ehart
  !
  !  counters?
  !integer isp, jj, kk, ii, dim, g1, g2, g3
  integer ig,ii,jj, kk,ijk,ijk2,ir,cmax_index
  !
  integer :: i,j,k,iadd,js
  real(dp) :: tmpr,tmpi, invgsqr, cmax, abschd,fac
  !real(dp) :: tmpd
  complex(dpc) :: ztmp ,fi_s
  !


  ! work arrays for the master processor

  ! total charge density collected to master processor
  real(dp), allocatable :: rho_master(:)
  ! hartree potential obtained on the master processor
  real(dp), allocatable :: vhart(:)
  ! actual hartree potential to be calculated
  complex(dpc),dimension(:), allocatable :: vhart_g
  ! charge density has to be in g space beforehand
  complex(dpc), dimension(:), allocatable :: rho_tmp, rho_g

  !timers
  real(dp) t0,t1,tbase,t3
  ! allocation check
  integer alcstat
  !---------------------------------------------------------------

!setup  realspace-planewave transfer parameters
   report_ehart = .TRUE.
! already done so?

! get the entire charge density in real space
! first collect rho from all PEs
      call collect_function(parallel,my_rho(1,1))

! time for single PE work only :(
  if (parallel%iammaster) then

     allocate(rho_master(grid%nwedge),stat=alcstat)
     call alccheck('rho_master',grid%nwedge,alcstat)
     rho_master(:) = zero
          call dcopy(grid%nwedge,parallel%ftmp,1,rho_master,1)

     allocate(rho_tmp(pbc%maxdfft),stat=alcstat)
     call alccheck('rho_tmp',pbc%maxdfft,alcstat)
     rho_tmp(:) = zzero 

  ! "unfold"
  ! Remove the mapping from charge density  "unfold"
  ! in order to perform the FFT

  do k = pbc%mz, pbc%mz + pbc%n3 - 1
     do j = pbc%my, pbc%my + pbc%n2 - 1
        do i = pbc%mx, pbc%mx + pbc%n1 - 1
           iadd = ((k-pbc%mz)*pbc%n2 + j-pbc%my)*pbc%n1 + i-pbc%mx+1
           jj = grid%rindex( grid%indexg(i,j,k) )
           !consider switching with daxpy if this is slow
              rho_tmp(iadd) = rho_master(jj) !* pbc%vcell !* 2*4pi??
        enddo
     enddo
  enddo
  ! rho_master is not needed anymore
     deallocate(rho_master)

    ! transform rho_tmp  to rho_g
    ! this is "forward" and being normalized by the vector length
     call cfftw (rho_tmp,pbc%n1,pbc%n2,pbc%n3,1)

     allocate(rho_g(pbc%nstar),stat=alcstat)
     call alccheck('rho_g',pbc%nstar,alcstat)

  !but first you need "symmetrization"
  ! as written in pot_local : Must add phase from displacement of FFT origin:
  !                  coordinates of corner of FFT box are not (0,0,0).
  !

     do jj=1,pbc%ng
        js = pbc%inds(jj)
        i = 1 + pbc%kgv(1,jj)
        if (i <= 0) i=pbc%n1+i
        j = 1 + pbc%kgv(2,jj)
        if (j <= 0) j=pbc%n2+j
        k = 1 + pbc%kgv(3,jj)
        if (k <= 0) k=pbc%n3+k
        call get_address(pbc,i,j,k,iadd,fi_s)
            !iadd from get_adress
        ztmp = rho_tmp(iadd) 
        tmpr = real(ztmp,dp)
        ! preordained whether to multiply by -1
        tmpi = pbc%conj(jj)*aimag(ztmp)
        !weighted by number of vectors in each star ?
        ztmp = cmplx(tmpr,tmpi,dp)/real(pbc%mstar(js),dp)
        ! so the phase is important then.
        rho_g(js) = rho_g(js) + ztmp * pbc%phase(jj)
     enddo

  ! don't need you too anymore 
     deallocate(rho_tmp)


     write(7,*) "hpot_G will now calculate the hartree potential in G-space (DEBUG only)"
     write(7,*) "JUST SAYING: IT WILL NOT WORK"
! from here on we calculate the hartree potential and energy in G-space
     allocate(vhart_g(pbc%ng),stat=alcstat)
     call alccheck('vhart_g',pbc%ng,alcstat)
     vhart_g(:) = zzero
     ehart_g = zero
     ! 
     !fac = one/(twopi*pbc%vcell)
    ! fac = one/(twopi)
!QES:  fac = e2*fpi/2piba2 -> ehart = ehart*fac vhart=vhart*fac
      !we need e^2*4pi/(2pi/a)^2
      !fac = two*four*pi/(two*pi/ pbc%vcell**(one/three) )**2
      !fac = two*four*pi/pbc%vcell
      fac = two*four*pi
 
     write(7,*) "hpot_G: total charge from rho_g(1)*pbc%vcell = ",rho_g(1)*pbc%vcell
     write(7,*) "hpot_G: hcub,vcell,fac - ",grid%hcub,pbc%vcell,fac
    ! charge = rho_g(1) ! needed?
!$OMP PARALLEL DO &
!$OMP& SCHEDULE(RUNTIME)  &
!$OMP& DEFAULT(SHARED) &
!$OMP& PRIVATE(js,invgsqr,tmpr,tmpi) &
!$OMP& REDUCTION(+:ehart_g)
    do ig=2,pbc%ng

     !for the hartree potential
      ! "g^2" = two*pbc%ek(g_index)
      !"factor" = 1/g^2
      ! dunno about other factors right now, maybe 2 or 2pi ?

      !which star am I on?
      js = pbc%inds(ig)
      !use the kinetic energy to get G^2
      invgsqr = one/(2*pbc%ek(js))
      if (abs(invgsqr) > 1e4) then
          write(9,*) "invgsqr is large",invgsqr
      endif
      !invgsqr = fac/pbc%ek(js)
      
!      write(7,*) "DEBUGGG g^2", 1/invgsqr

     !so for the hartree energy - 
      tmpr =  real(rho_g(ig),dp)
      tmpi = aimag(rho_g(ig))
     ehart_g = ehart_g + (tmpr**2 + tmpi**2) * invgsqr
     
     ! vhart_g(ig) = invgsqr*rho_g(ig)
     !yeah, I know it's not elegant but
     !this is in the spirit of what quantum espresso does
     vhart_g(ig)=CMPLX(tmpr*invgsqr*fac,tmpi*invgsqr*fac,KIND=dp)
!      write(7,*) "DEBUGGG invgsqr (ig), vhart_g ,ig", invgsqr, vhart_g(ig), ig

     end do
!$OMP END PARALLEL DO

      write(7,*) "DEBUGG - sample 25 elements of vhart_g",vhart_g(1:25)
      write(7,*) "DEBUGG - sample 25 elements of rho_g",rho_g(1:25)

! scale the energy term?
!QES:    if gamma only ehart=ehart*omega
!QES:   if not gamma only ehart=0.5*ehart*omega
     ehart_g = ehart_g *fac/two !*pbc%vcell


    ! transform v_hartree_g to v_hart

    ! done like in pot_local
    pbc%vscr4 = zzero
    pbc%vscr2 = zero

     do jj=1,pbc%ng
        i = 1 + pbc%kgv(1,jj)
        if (i <= 0) i=pbc%n1+i
        j = 1 + pbc%kgv(2,jj)
        if (j <= 0) j=pbc%n2+j
        k = 1 + pbc%kgv(3,jj)
        if (k <= 0) k=pbc%n3+k
        call get_address(pbc,i,j,k,iadd,fi_s)
      pbc%vscr4(iadd) = pbc%vscr4(iadd)  + vhart_g(jj) * fi_s
   enddo

     allocate(vhart(grid%nwedge),stat=alcstat)
     call alccheck('vhart',grid%nwedge,alcstat)
     vhart(:) = zero

  ! Fourier transform to real space.
   call cfftw (pbc%vscr4(:),pbc%n1,pbc%n2,pbc%n3,-1)

   ! dmax = real( pbc%vscr4(1) ,dp)
   ! dmin = dmax
    cmax = zero
    cmax_index=1
   ! ierr = 0
   do k=1,pbc%n3
      do j=1,pbc%n2
         do i=1,pbc%n1
            ijk2 = ((k-1)*pbc%n2 + j-1)*pbc%n1 + i
            pbc%vscr2(ijk2) = real( pbc%vscr4(ijk2) ,dp)
            ! if (pbc%vscr2(ijk2) > dmax) dmax=pbc%vscr2(ijk2)
            ! if (pbc%vscr2(ijk2) < dmin) dmin=pbc%vscr2(ijk2)
             abschd = abs( aimag( pbc%vscr4(ijk2) ) )
             if (abschd > cmax) then
                 cmax=abschd
                 cmax_index=ijk2
             endif
            ! if (abschd > small) ierr=ierr+1
         enddo
      enddo
   enddo
   write(7,*) "hpot_G, maximum complex value for the hartree potential was",cmax
   write(7,*) "hpot_G, reported for g index",cmax_index

     do ir = 1, grid%nwedge
        i = grid%kx(ir)
        j = grid%ky(ir)
        k = grid%kz(ir)
        ijk = ((k-pbc%mz)*grid%n2 + j-pbc%my)*grid%n1 + i-pbc%mx+1
        ii = grid%indexw(i,j,k)
        vhart(ii) = pbc%vscr2(ijk) !*two !at this stage i dunno what i'm doing with all these factors
     enddo

     parallel%ftmp=vhart
     deallocate(vhart)
  endif ! done only master.

! send v_hart to all PEs 
  call export_function(parallel,my_vhart)
! recieve my v_hart part
! make sure the leftovers are zeroed out 
  my_vhart(parallel%mydim+1:parallel%ldn) = zero

  if (parallel%iammaster) then
      ! here you should report timing
      !
      if (report_ehart) then
      !report the hartree energy as calculated here
      write(7,*) "hpot_G - I registered the Hartree energy as",ehart_g
      endif
  endif

! done.


end subroutine hpot_G
!===============================================================
