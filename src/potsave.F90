!========================================================================================
! (This subroutine might need a suffice .f90p)
!
!       This subroutine saves grid, hartree potential, ionic potential, 
!       exchange-correlation potential seperately
!
!       This subroutine is a part of AFM project
!
!       January 2012,   minjung@ices.utexas.edu
!========================================================================================

! modification has been made for a new AFM simulations

subroutine potsave(clust,elec_st,grid,pbc,rsymm,vion,vhart,vxc,rho,rhoc,newk,parallel,iunit)

  use constants
  use cluster_module
  use electronic_struct_module
  use grid_module
  use pbc_module
  use symmetry_module
  use parallel_data_module

  implicit none
#ifdef MPI
  include 'mpif.h'
#endif
  !
  !  Input/Output variables:
  !
  !  atomic coordinates
  type (cluster), intent(in) :: clust
  !  electronic structure
  type (electronic_struct), intent(inout) :: elec_st
  !  grid related data
  type (grid_data), intent(in) :: grid
  !  periodic boundary conditions data
  type (pbc_data), intent(in) :: pbc
  !  symmetry operations:
  type (symmetry), intent(in) :: rsymm
  !  parallel computation related data
  type (parallel_data), intent(inout) :: parallel
  !  potential related data
  real(dp),intent(in) :: vion(parallel%mydim)
  real(dp),intent(in) :: vhart(parallel%mydim)
  real(dp),intent(in) :: vxc(parallel%mydim,elec_st%nspin)
  real(dp),intent(in) :: rho(parallel%mydim,2*elec_st%nspin-1)
  real(dp),intent(in) :: rhoc(parallel%mydim)
  real(dp),intent(in) :: newk

  !  number of output unit
  integer, intent(in) :: iunit


  !  Work variables:
  !
  !  communicator
  integer comm
  !  exit code for mpi calls
  integer mpinfo
  !  counters
  integer isp, jj, kk, ii, dim
  real(dp) :: tmpd
  !
  !  flags for type of system
  !  stype(1) = nspin
  !  stype(2) = 0 : real wavefunctions
  !  stype(2) = 1 : complex wavefunctions
  !  stype(2) = 2 : spin-orbit (complex doubled wavefunctions)
  !  stype(3) = 0 : confined system (cluster)
  !  stype(3) = 1 : system periodic on all directions (bulk)
  !  stype(3) = 2 : system periodic on x direction only (wire)
  !
  integer stype(10)


#ifdef MPI
  integer status(MPI_STATUS_SIZE)
#endif

  ! iunit is 126 in parsec.f90p
  if (parallel%iammaster) open(iunit,file='pot.dat',form='unformatted')
  comm = parallel%comm
  !
  !  Build system flags.
  !
  stype = 0 ! initialization
  ! number of spin
  stype(1) = elec_st%nspin
  ! wavefunction type
  if (elec_st%cplx) then
     if(elec_st%is_so) then
         stype(2) = 2
     else
         stype(2) = 1
     endif
  endif
  ! periodicity
  select case(pbc%per) ! pbc%per: number of periodic dimension
  case(0)
     stype(3) = 0
  case(1)
     stype(3) = 2
  case(2)
     stype(3) = 3
  case(3)
     stype(3) = 1
  end select

  dim = parallel%mydim*parallel%mxwd

!====================================== pot.dat starts here

  if (parallel%iammaster) then
     rewind(iunit)  ! why is it needed?

!##### write the actual coordinates
     write(iunit) newk
     write(iunit) clust%atom_num
     do ii = 1,clust%atom_num
        write(iunit) clust%xatm(ii), clust%yatm(ii), clust%zatm(ii)!, clust%intback(:,ii)
     enddo

     write(iunit) stype

     if (pbc%is_on) then
        write(iunit) (grid%step(kk),kk=1,3)
        write(iunit) (pbc%box_size(kk),kk=1,3)
        write(iunit) ((pbc%latt_vec(jj,kk),jj=1,3),kk=1,3)
     else
        write(iunit) grid%step(1), grid%rmax
     endif


     write(iunit) grid%ndim ! total number of grid
     write(iunit) grid%shift
     write(iunit) grid%nwedge,rsymm%ntrans
     write(iunit) rsymm%rmtrx
     write(iunit) rsymm%trans
     write(iunit) rsymm%tnp
     write(iunit) rsymm%alatt, rsymm%invlat
     write(iunit) ((rsymm%chi(kk,jj),kk=1,rsymm%ntrans),jj=1,rsymm%ntrans)
     write(iunit) (grid%kx(kk),grid%ky(kk),grid%kz(kk),kk=1,grid%nwedge)
  endif

    ! note :  these potentials are not mixed data for next run
    ! ionic potential
    call collect_function(parallel,vion(1))
    if (parallel%iammaster) &
       write(iunit) (parallel%ftmp(kk),kk=1,grid%nwedge)
    ! hartree potential
    call collect_function(parallel,vhart(1))
    if (parallel%iammaster) &
       write(iunit) (parallel%ftmp(kk),kk=1,grid%nwedge)
    ! exchange-correlation potential
    call collect_function(parallel,vxc(1,1))
    if (parallel%iammaster) &
       write(iunit) (parallel%ftmp(kk),kk=1,grid%nwedge)

    ! -- If spin-orbit potential is present, 
    ! --    second spin compotent must be printed out 
!   if (elec_st%is_so) then
!      call collect_function(parallel,pot%vxc(1,2))
!      write(iunit) (parallel%ftmp(kk),kk=1,grid%nwedge)
!   endif

!---------------------------------------------------------------------------
    ! Let's save the density as well at this moment
    ! I am ignoring the spin right now
    call collect_function(parallel,rho(1,1))
    if (parallel%iammaster) &
       write(iunit) (parallel%ftmp(kk),kk=1,grid%nwedge)
    !Core charge
    call collect_function(parallel,rhoc(1))
    if (parallel%iammaster) &
       write(iunit) (parallel%ftmp(kk),kk=1,grid%nwedge)
!---------------------------------------------------------------------------


#ifdef MPI
  call MPI_Barrier(comm,mpinfo)
#endif

  if (parallel%iammaster) close(iunit) 

end subroutine
