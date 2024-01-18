!===============================================================
!
! Copyright (C) 2005 Finite Difference Research Group
! This file is part of parsec, http://www.ices.utexas.edu/parsec/
!
! Computes the local potential contribution to the Hellman-Feynman
! forces. this is done in reciprocal space taking advantage of
! its "short-ranged nature". the charge density and exchange and
! correlation potential are first transferred to reciprocal
! space by Fast Fourier Transform.
!
! WARNING: this subroutine is not fully parallelized!
! AJB 2015: URGENT - parallelize this
!
!---------------------------------------------------------------
subroutine forpbc(clust,elec_st,grid,pbc,parallel,rho_d,vxc_d)

  use constants
  use cluster_module
  use electronic_struct_module
  use grid_module
  use pbc_module
  use parallel_data_module
  implicit none
  !
  ! Input/Output variables:
  !
  ! the cluster
  type (cluster), intent(inout) :: clust
  ! electronic structure
  type (electronic_struct), intent(inout) :: elec_st
  ! grid related data
  type (grid_data), intent(in) :: grid
  ! pbc related data
  type (pbc_data), intent(inout) :: pbc
  ! parallel computation related data
  type (parallel_data), intent(inout) :: parallel

  ! distributed electron density and exchange-correlation potential
  ! (passed outside the structure to overcome a bug with the IBM compiler)
  real(dp), intent(in) :: rho_d(parallel%mydim)
  real(dp), intent(in) :: vxc_d(parallel%mydim,elec_st%nspin)
  !
  ! Work variables:
  !
  ! allocation check
  integer alcstat
  ! exchange-correlation potential collected to master processor
  real(dp), allocatable :: vxc(:)
  ! total charge density collected to master processor
  real(dp), allocatable :: rho(:)
  ! timing
  real(dp) t0,t1,tbase,t3

  integer i,isp,j,jj,k,iadd,js
  real(dp) :: tmpr,tmpi
  complex(dpc) :: ztmp,fi_s
  complex(dpc), dimension(:), allocatable :: vsc3,vsc,vsc_chg

  !---------------------------------------------------------------
!write(9,*) 'Inside forpbc local'

  ! Collect vxc over all processors and add spin components.
  if (parallel%iammaster) then
     call mysecond(tbase)

     allocate(vxc(grid%nwedge),stat=alcstat)
     call alccheck('vxc',grid%nwedge,alcstat)
     vxc(:) = zero
     allocate(rho(grid%nwedge),stat=alcstat)
     call alccheck('rho',grid%nwedge,alcstat)
     rho(:) = zero

     call mysecond(t0)
  endif

  do isp = 1, elec_st%nspin
     call collect_function(parallel,vxc_d(1,isp))
     if (parallel%iammaster) call daxpy(grid%nwedge,one,parallel%ftmp,1,vxc,1)
  enddo
  ! Collect charge density over all processors.
  ! AJB: This does not seem right? why loop over spin
  do isp = 1, elec_st%nspin
     call collect_function(parallel,rho_d)
     if (parallel%iammaster) call dcopy(grid%nwedge,parallel%ftmp,1,rho,1)
  enddo

  ! From now on, only master PE works.
  if (.not. parallel%iammaster) return

     call mysecond(t1)

     write(9,*) 
     write(9,100) t1-t0
100  format('  forpbc time to collect rho,vxc: [sec]:',1x,f10.2)

  allocate(vsc3(pbc%maxdfft),stat=alcstat)
  call alccheck('vsc3',pbc%maxdfft,alcstat)
  allocate(vsc(pbc%nstar),stat=alcstat)
  call alccheck('vsc',pbc%nstar,alcstat)
  allocate(vsc_chg(pbc%nstar),stat=alcstat)
  call alccheck('vsc_chg',pbc%nstar,alcstat)

  vsc3(:) = zzero
  pbc%vscr4(:) = zzero

     call mysecond(t0)
  ! Remove the mapping from charge density and exchange-correlation 
  ! potential in order to perform the FFTs.
  do k = pbc%mz, pbc%mz + pbc%n3 - 1
     do j = pbc%my, pbc%my + pbc%n2 - 1
        do i = pbc%mx, pbc%mx + pbc%n1 - 1
           iadd = ((k-pbc%mz)*pbc%n2 + j-pbc%my)*pbc%n1 + i-pbc%mx+1
           jj = grid%rindex( grid%indexg(i,j,k) )
           pbc%vscr4(iadd) = vxc(jj)/real(2*elec_st%nspin,dp)
           vsc3(iadd) = rho(jj) * pbc%vcell
        enddo
     enddo
  enddo

     call mysecond(t1)

     write(9,*) 
     write(9,101) t1-t0
101  format('  forpbc time to unmap to fftw[sec]:',1x,f10.2)

     call mysecond(t0)
  ! Perform the FFTs.
  call cfftw (pbc%vscr4,pbc%n1,pbc%n2,pbc%n3,1)
  call cfftw (vsc3,pbc%n1,pbc%n2,pbc%n3,1)

     call mysecond(t1)

     write(9,*) 
     write(9,102) t1-t0
102  format('  forpbc fftw time  [sec]:',1x,f10.2)
  ! Do some wraparound.

     call mysecond(t0)
  vsc(:) = cmplx(zero,zero,dp)
  vsc_chg(:) = cmplx(zero,zero,dp)
  do jj=1,pbc%ng
     js = pbc%inds(jj)
     i = 1 + pbc%kgv(1,jj)
     if (i <= 0) i=pbc%n1+i
     j = 1 + pbc%kgv(2,jj)
     if (j <= 0) j=pbc%n2+j
     k = 1 + pbc%kgv(3,jj)
     if (k <= 0) k=pbc%n3+k
     call get_address(pbc,i,j,k,iadd,fi_s)

     ztmp = pbc%vscr4(iadd) * conjg(fi_s)
     tmpr = real(ztmp,dp)
     tmpi = pbc%conj(jj)*aimag(ztmp)
     ztmp = cmplx(tmpr,tmpi,dp)/real(pbc%mstar(js),dp)
     vsc(js) = vsc(js) + ztmp

     ztmp = vsc3(iadd) * conjg(fi_s)
     tmpr = real(ztmp,dp)
     tmpi = pbc%conj(jj)*aimag(ztmp)
     ztmp = cmplx(tmpr,tmpi,dp)/real(pbc%mstar(js),dp)
     vsc_chg(js) = vsc_chg(js) + ztmp * pbc%phase(jj)
  enddo
     call mysecond(t1)

     write(9,103) t1-t0
103  format('  forpbc fftw wrap-around time  [sec]:',1x,f10.2)

     call mysecond(t3)
     write(9,104) t3-tbase
104  format('  forpbc  total time before local-force contrib.r [sec]:',1x,f10.2)
  ! Finally, calculate the local force on ions.
  call force_loc_pb(clust,pbc,vsc,vsc_chg)

     call mysecond(t1)
     write(9,105) t1-t3
105  format('  forpbc  time for local contribution [sec]:',1x,f10.2)
     write(9,*)  

  if (allocated(vsc3)) deallocate(vsc3)
  if (allocated(vsc)) deallocate(vsc)
  if (allocated(vsc_chg)) deallocate(vsc_chg)
  if (allocated(vxc)) deallocate(vxc)
  if (allocated(rho)) deallocate(rho)

end subroutine forpbc
!===============================================================
