!===============================================================
!
!  Copyright (C) 2005 Finite Difference Research Group
!  This file is part of parsec, http://www.ices.utexas.edu/parsec/
!
!  This subroutine calls subroutines to set up the Hartree
!  potential and solve Poisson's equation from a known charge
!  density.
!
!---------------------------------------------------------------
subroutine rho_har(grid,pot,rsymm,parallel,elec_st,pbc,lpole,ipr,exc)

  use constants
  use grid_module
  use potential_module
  use symmetry_module
  use parallel_data_module
  use electronic_struct_module
  use pbc_module
#ifdef MPI
  use mpi
#endif
  implicit none
  !
  !  Input/Output variables:
  !
  !  grid related data
  type (grid_data), intent(in) :: grid
  !  potential related data
  type (potential), intent(inout) :: pot
  !  symmetry operations in reduced group:
  type (symmetry), intent(in) :: rsymm
  !  parallel computation related data
  type (parallel_data), intent(inout) :: parallel
  !  electronic structure
  type (electronic_struct), intent(inout) :: elec_st
  !  periodic boundary conditions data
  type (pbc_data), intent(in) :: pbc
  !  order of multipole expansion
  integer, intent(in) :: lpole
  !  print flag
  integer, intent(in) :: ipr
  !  total exchange-correlation energy
  real(dp), intent(out) :: exc
  !
  !  Work variables:
  !
  !  counters
  integer :: isp, i, idx
  !  temporary variables
  real(dp) :: rhoav, vxcav, rhoav_test
  !  maximum and minimum values of charge density, in [e/bohr^3]
  real(dp), dimension(elec_st%nspin) :: dmax, dmin, dtmp
  !  work array
  real(dp) :: rho(parallel%ldn)
  !  spin descriptor
  character (len=2) :: idsp(elec_st%nspin)
  !  boundary rho from multipole expansion
  real(dp), allocatable :: brho(:)
  !  work array for Hartree CG
  real(dp), allocatable :: vect(:)
  !  allocation check
  integer alcstat
  !  printout flag
  logical :: prflag
#ifdef MPI
  integer mpinfo
#endif

  !---------------------------------------------------------------

  !  Retrieve the current charge density
  call dcopy(parallel%mydim,elec_st%rho(1,1),1,rho(1),1)
  rho(parallel%mydim+1:parallel%ldn) = zero

  !  set spin identifier
  do isp = 1, elec_st%nspin
     if (isp-1+elec_st%nspin == 1) idsp(isp) = '  '
     if (isp-1+elec_st%nspin == 2) idsp(isp) = 'up'
     if (isp-1+elec_st%nspin == 3) idsp(isp) = 'dn'
  enddo

  allocate(vect(parallel%nwedge+1),stat=alcstat)
  call alccheck('vect',parallel%nwedge+1,alcstat)
  allocate(brho(parallel%ldn),stat=alcstat)
  call alccheck('brho',parallel%ldn,alcstat)
!  brho(:) = zero !<=== AJB: brho is intent(out) in all the subroutines below

select case(pbc%per)
case (0)
        !  For confined systems, set boundary condition for Hartree
        !  potential.
        call hartset(grid, rsymm, parallel, grid%norder, lpole, &
                grid%coe2, rho, pot%clm, brho)
case (1)
        !  For wire/tube systems, set boundary condition for Hartree
        !  potential. Include the monopole contribution but ignore the
        !  logarithmic divergence (it will be cancelled by the ion-ion
        !  energy and local pseudopotentials).
        call hartset_wire(grid, pbc, rsymm, parallel, grid%norder, lpole, &
                grid%coe2, rho, brho)
case (2)
        !  For slab systems, set boundry condition for Hartree potential.
        !  Include the dipole contribution.
        call hartset_slab(grid, pbc, parallel, grid%norder, lpole, &
                grid%coe2, rho, brho)
case (3)
        if (elec_st%explicit_hartree_pbc) then
                ! for testing purposes - obtain the hartree potential using FFT
                ! do I use rho or brho?
                write(9,*) ' Oh boy, obtaining the Hartree potential in reciprocal space '

                call hpot_G(elec_st,grid,pbc,parallel,elec_st%rho,rho)
                ! rho now contains v_hart on the grid
        else
                ! Do it numerically

                !  For periodic system, find the average charge density and
                !  subtract it out to make the average zero.

                rhoav = elec_st%xele/(grid%hcub*parallel%ndim)

                ! write(9,*) ' 3D periodic system, setting the average of the charge density to zero '
                ! write(9,*) ' straightforward rhoav =',rhoav

                ! write(9,*) ' testing via averaging over rho '

                rhoav_test = sum(rho(1:parallel%mydim))
                call psum(rhoav_test,1,parallel%group_size,parallel%group_comm)
                write(9,*) ' is this the correct total charge?=',rhoav_test
                rhoav_test = rhoav_test/real(parallel%ndim)
                write(9,*) ' average charge total*vcell =',rhoav_test*pbc%vcell

                !     brho(:) = zero !since it was commented above, it should be
                !     !noted that brho is not initalized anymore,
                !     !so it is proper that we make sure the last remaining elements are zeroed out

                brho(parallel%mydim+1:parallel%ldn) = zero
                !
                !  Also, multiply the result by 8*pi in preparation for the
                !  Poisson solver (4*pi coming from the form of the equation, 2
                !  from changing Hartree to Rydberg).

                brho(1:parallel%mydim) =  &
                  (rho(1:parallel%mydim) - rhoav)*pi*eight

                ! write(9,*) ' now testing via averaging over b-rho '

                ! rhoav_test = sum(brho(1:parallel%mydim))
                ! call psum(rhoav_test,1,parallel%group_size,parallel%group_comm)
                ! write(9,*) ' is this the correct new total charge?=',rhoav_test
                ! rhoav_test = rhoav_test/real(parallel%ndim)
                ! write(9,*) ' is this the correct average?=',rhoav_test
        endif

end select

  if( .not. elec_st%explicit_hartree_pbc) then
      !  Retrieve the current Hartree potential and use it as initial
      !  guess for the Poisson solver. If there is no available Hartree
      !  potential, use the negative of the local pseudopotential.
      !
      !  AJB: this check is wasteful. should be a flag.
      if (sum(pot%vhart) == zero) then
         write(9,*) ' init Hartree guess, starting from Vlocal '
         call dcopy(parallel%mydim,pot%vion(1),1,rho(1),1)
         ! wouldn't it be better to start form the laplacian of vion?
         rho = -rho
    ! AJB: maybe - should daxpy instead, but rho is not zero since it was assigned before...
    !    call daxpy(parallel%mydim,mone,pot%vion(1),1,rho(1),1)
      else
         call dcopy(parallel%mydim,pot%vhart(1),1,rho(1),1)
      endif
      
      call hpotcg(parallel,brho,rho,grid%norder,grid%coe2, grid%lap_dir_num,vect)
  !else
      ! rho already contains hartree potential done explicitly
  endif

  deallocate(brho)

  !  Store the Hartree potential.
  call dcopy(parallel%mydim,rho(1),1,pot%vhart(1),1)

  !  Although this subroutine is called rho_har...
  !
  !  Compute the exchange-correlation potential.
  !
  ! start by adding rhoc to rho and storing as vxc.
  !
  !  "Vxc is used as core-corrected charge on input"
  !   try not to be confused by the use of vxc as charge density
  !
  !  and exchange-correlation potential on output.
  !
  ! there should be a check here like in corecd, so this won't happen
  ! all the time
  if (elec_st%nspin == 1) then
     do i = 1, parallel%mydim
        pot%vxc(i,1) = elec_st%rho(i,1) + elec_st%rhoc(i)
     enddo
  else
     do isp = 1, elec_st%nspin
        idx = isp+1
        do i = 1, parallel%mydim
           pot%vxc(i,isp) = elec_st%rho(i,idx) + half*elec_st%rhoc(i)
        enddo
     enddo
  endif

  !  Find max and min values of the charge density and report
  !  warn (but do not kill) if negative values found.

  dmax = maxval(pot%vxc,1)
  dmin = minval(pot%vxc,1)
#ifdef MPI
  call MPI_ALLREDUCE(dmax,dtmp,elec_st%nspin,MPI_DOUBLE_PRECISION, &
       MPI_MAX,parallel%group_comm,mpinfo)
  dmax = dtmp
  call MPI_ALLREDUCE(dmin,dtmp,elec_st%nspin,MPI_DOUBLE_PRECISION, &
       MPI_MIN,parallel%group_comm,mpinfo)
  dmin = dtmp
#endif

  if (parallel%iammaster) then
     do isp = 1, elec_st%nspin
        write(7,12) dmax(isp),dmin(isp),idsp(isp)
12      format(/,' Max and min values of charge density [e/bohr^3]:' &
             ,2x,e11.4,2x,e11.4,3x,a2)
     enddo
     if (minval(dmin) < zero) &
          write(7,*) 'WARNING : NEGATIVE CHARGE DENSITY FOUND!'
  endif

  if (elec_st%nspin == 1) then
     call exc_nspn(grid,parallel,elec_st%nrep,elec_st%icorr &
          ,elec_st%dioniz,pot%vxc,vect,exc)
  else
     call exc_spn(grid,parallel,elec_st%nrep,elec_st%icorr &
          ,elec_st%dioniz,pot%vxc,vect,exc)
  endif
  !
  !  If periodic boundary conditions are used, shift vxc uniformly so that
  !  its average is set to zero; do the same for vhart.
  !
  if (pbc%per == 3) then
     prflag = .false.
     if (ipr >= 0 .and. parallel%iammaster) prflag = .true.
     !if (ipr >= 1 .and. parallel%iammaster) prflag = .true.
     do isp = 1, elec_st%nspin
        call potshift(pbc%vcell,grid%hcub,elec_st%nrep &
             ,parallel%mydim,parallel%group_size &
             ,parallel%group_comm,pot%vxc(1,isp),vxcav)
        if (prflag) write(7,'(/,a,2x,f20.8,2x,a2)') &
             ' Average exchange-correlation potential [Ry]:',vxcav,idsp(isp)
     enddo
     call potshift(pbc%vcell,grid%hcub,elec_st%nrep &
          ,parallel%mydim,parallel%group_size &
          ,parallel%group_comm,pot%vhart,vxcav)
     if (prflag) write(7,'(a,2x,f20.8,/,a)') &
          ' Average Hartree potential [Ry]:',vxcav, &
          ' Attention: Average potentials reset to zero.'
  endif
  if (allocated(vect)) deallocate(vect)

end subroutine rho_har
!===============================================================
