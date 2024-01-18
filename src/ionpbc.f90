!===============================================================
!
! Copyright (C) 2005 Finite Difference Research Group
! This file is part of parsec, http://www.ices.utexas.edu/parsec/
!
! Calculates the ionic local potential (vion) on the real space
! mesh for a given ionic configuration (xatm,yatm,zatm).
!
! This subroutine is partially threaded, but not MPI
!
!---------------------------------------------------------------
subroutine ionpbc(clust,grid,pbc,parallel,vion_d,ipr)

  use constants
  use cluster_module
  use grid_module
  use pbc_module
  use parallel_data_module
  implicit none
  !
  ! Input/Output variables:
  !
  ! the cluster
  type(cluster), intent(in) :: clust
  ! grid related data
  type(grid_data), intent(in) :: grid
  ! pbc related data
  type(pbc_data), intent(inout) :: pbc
  ! parallel computation related data
  type (parallel_data), intent(inout) :: parallel
  ! local part of ionic potential
  real(dp), intent(inout) :: vion_d(parallel%mydim)
  ! printout flag
  integer, intent(in) :: ipr
  !
  ! Work variables:
  !
  ! counters
  integer i, j, k, ijk, ii, ir
  ! allocation check
  integer alcstat
  ! timing
  real(dp) t0,t1,tbase,t3
  ! work arrays for the master processor
  real(dp), allocatable :: vion(:)
  complex(dpc),dimension(:), allocatable :: vionz
  real(dp) :: tmpvec(3)

  !---------------------------------------------------------------

  if (parallel%iammaster) then

     ! Initialize potentials.
     allocate(vion(grid%nwedge),stat=alcstat)
     write(7,*) 
     call alccheck('vion',grid%nwedge,alcstat)
     vion(:) = zero
     allocate(vionz(pbc%ng),stat=alcstat)
     call alccheck('vionz',pbc%ng,alcstat)

     call mysecond(tbase)

     ! Define coordinates of atoms with respect to the corner of the
     ! periodic cell. This is needed for the FFTs and structure factor.
     ii = 0
     do i = 1, clust%type_num
        do j = 1, clust%natmi(i)
           ii = ii + 1

           tmpvec(1)=clust%xatm(ii)
           tmpvec(2)=clust%yatm(ii)
           tmpvec(3)=clust%zatm(ii)

           ! Multiplying by the bvec matrix to find 
           ! the coordinates relative to lattice vec.
           call matvec3('T',pbc%bvec,tmpvec,pbc%rat(:,j,i))

        enddo
     enddo


     ! Obtain the local ionic potential on the reciprocal space mesh.

     ! call v_first(pbc,clust%natmi,vionz)
     ! vionz(1) = cmplx(zero,zero,dp)

     ! write(9,*) 
     ! write(9,100) t1-t0
! 100  format('  ionpbc v_first time [sec]:',1x,f10.2)

!testing openmp version
     call mysecond(t0)
     call v_first_omp(pbc,clust%natmi,vionz)
     vionz(1) = cmplx(zero,zero,dp)

     call mysecond(t1)
     write(9,*) 
     write(9,101) t1-t0
101  format('  ionpbc v_first_omp time [sec]:',1x,f10.2)
     call mysecond(t0)

     ! Transfer the local ionic potential to the real space grid ...
     call pot_local(pbc,ipr,vionz)

     call mysecond(t1)
     write(9,102) t1-t0
102  format('  ionpbc pot_local time [sec]:',1x,f10.2)

     call mysecond(t0)
     ! ... and finally store it according to mapping.
     do ir = 1, grid%nwedge
        i = grid%kx(ir)
        j = grid%ky(ir)
        k = grid%kz(ir)
        ijk = ((k-pbc%mz)*grid%n2 + j-pbc%my)*grid%n1 + i-pbc%mx+1
        ii = grid%indexw(i,j,k)
        vion(ii) = pbc%vscr2(ijk)*two
     enddo
     parallel%ftmp = vion
     deallocate(vionz, vion)

     call mysecond(t1)
     write(9,103) t1-t0
103  format('  ionpbc store/index time [sec]:',1x,f10.2)

  endif                     ! parallel%iammaster
  ! Master PE distributes the ion potential.
  call export_function(parallel,vion_d)

  if (parallel%iammaster) then
     call mysecond(t3)
     write(9,104) t3-tbase
104  format('  ionpbc total [sec]:',1x,f10.2)
     write(9,*)  
  endif

end subroutine ionpbc
!===============================================================
