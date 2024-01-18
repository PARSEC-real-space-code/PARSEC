subroutine calc_nscf(nscf,elec_st,pbc,pot,u_pot,solver,parallel, &
           grid, nloc_p_pot, clust, p_pot, ipr,ierr)

  use constants
  use cluster_module
  use electronic_struct_module
  use potential_module
  use non_local_psp_module
  use eigen_solver_module
  use parallel_data_module
  use pseudo_potential_module
  use nscf_module
  use pbc_module
  use grid_module
#ifdef MPI
  ! mpi definitions
  use mpi
#endif
  
  implicit none

  !
  ! Input/Output variables:
  !
  ! non selfconsistent data
  type(nscf_data), intent(inout) :: nscf
  ! electronic structure
  type (electronic_struct), intent(inout) :: elec_st
  !grid
  type (grid_data), intent(inout) :: grid
  ! pbc
  type (pbc_data), intent(inout) :: pbc
  ! potential related data
  type (potential), intent(in) :: pot
  ! on-site Coulomb interaction related data
  type (nonloc_pseudo_potential), intent(inout) :: u_pot
  ! solver related data
  type (eigen_solver), intent(inout) :: solver
  ! parallel computation related data
  type (parallel_data), intent(in) :: parallel
  ! non local pseudo potential part
  type (nonloc_pseudo_potential), intent(inout) :: nloc_p_pot
  ! pseudo_potential related data
  type (pseudo_potential), intent(in) :: p_pot
  ! cluster data
  type (cluster), intent(in) :: clust
  ! printout flag
  integer, intent(in) :: ipr
  ! error flag, 400 < ierr < 421
  integer, intent(out) :: ierr

  !
  ! Work variables:
  !
  integer :: irp,isp,kplp,scf_nkpt,i,offset,j,lastkpt,spnum,jj,nrep
  integer :: nmax,ii,nscf_nkpt
  real(dp), dimension(3)  :: vec, ktmp
  real(dp) :: min_length  
 ! Diagonalization timers
  real(dp) :: tdiag0, tdiag1, tdiag_sum 
! output info from MPI communication
  integer :: mpinfo




!until we figure out what to do with symmetry
  nrep = 1
  if (associated(grid%kecoe1)) deallocate(grid%kecoe1)
  allocate(grid%kecoe1(nscf%nkpt,6,0:grid%norder))
  grid%kecoe1(:,:,:) = zero
  do kplp = 1, nscf%nkpt
     do jj =0,grid%norder
        call matvec3('T',grid%grad_bvec_norm,nscf%kpts(1,kplp),ktmp)
        grid%kecoe1(kplp,1:3,jj) = -2*zi*grid%coe1(jj,1:3)*ktmp
     enddo
  enddo

  call destroy_eigen_solver(solver)
 
  call create_eigen_solver(solver,grid,parallel%ldn,nrep, &
       elec_st%nspin,nscf%nstate,nscf%nkpt,elec_st%cplx, &
       parallel%mxwd)



  do isp = 1, elec_st%nspin/elec_st%mxwd
    do kplp =1,max(elec_st%nkpt,1)
      do irp = 1, elec_st%nrep
        call destroy_eigenstate(elec_st%eig(irp,kplp,isp))
      enddo
    enddo
  enddo
  elec_st%nrep = 1

! Now change variables to fill nscf data in elec_st structure
  spnum = elec_st%nspin/elec_st%mxwd  
  nscf_nkpt = nscf%nkpt
  elec_st%nkpt = nscf%nkpt
  elec_st%nstate = nscf%nstate

! some array initializations...

  if (associated(elec_st%kpts)) deallocate(elec_st%kpts)
  allocate(elec_st%kpts(3,nscf%nkpt))
  if (associated(pbc%kpts)) deallocate(pbc%kpts)
  allocate(pbc%kpts(3,nscf%nkpt))
  if (associated(nloc_p_pot%kpts)) deallocate(nloc_p_pot%kpts)
  allocate(nloc_p_pot%kpts(3,nscf%nkpt))

  if (associated(elec_st%kpwt)) deallocate(elec_st%kpwt)
  allocate(elec_st%kpwt(nscf%nkpt))
  if (associated(elec_st%magmom)) deallocate(elec_st%magmom)
  allocate (elec_st%magmom (2*elec_st%nstate,nscf_nkpt))

  if (associated(elec_st%eig)) deallocate(elec_st%eig)
  allocate(elec_st%eig(elec_st%nrep,nscf_nkpt,elec_st%nspin))
  nmax = (solver%winsize + solver%nadd + nscf%nstate) * elec_st%nrep
  allocate(elec_st%irep(nmax,nscf_nkpt,elec_st%nspin))

  elec_st%irep = 1

  pbc%nkpt = elec_st%nkpt
  nloc_p_pot%nkpt = elec_st%nkpt

  elec_st%kpts(:,:)= nscf%kpts(:,:)
  pbc%kpts(:,:)=elec_st%kpts(:,:)
  nloc_p_pot%kpts(:,:)= elec_st%kpts(:,:)
  elec_st%kpwt(:) = nscf%kpwt(:)

  ii = 0
  do isp = 1, elec_st%nspin/elec_st%mxwd
    do kplp =1,max(elec_st%nkpt,1)
      do irp = 1, elec_st%nrep
        call create_eigenstate(elec_st%eig(irp,kplp,isp),parallel%ldn* &
                   parallel%mxwd,nscf%nstate+solver%winsize,elec_st%cplx)
        elec_st%eig(irp,kplp,isp)%group = ii
        ii = ii + 1
        if (mod(ii,parallel%groups_num) == 0) ii = 0
      enddo
    enddo
  enddo

  if (.not. solver%fix_neig) call initeigval(elec_st,solver%nadd)

  elec_st%eig(:,:,:)%mm = -1
  elec_st%eig(:,:,:)%nn = elec_st%nstate!+solver%winsize
  elec_st%eig(:,:,:)%nec = -1

  if (associated(elec_st%occ_in)) deallocate(elec_st%occ_in)
  allocate (elec_st%occ_in (elec_st%nstate,spnum))

  solver%maxmv = max(solver%maxmv,parallel%ldn*elec_st%nstate/elec_st%nrep)
  if (parallel%iammaster) then
     write(7,*)
     write(7,*) ' Maximum number of matrix-vector operations = ',solver%maxmv
     write(7,*)
  endif

 call nonloc(clust,grid,p_pot,nloc_p_pot,pbc,parallel,ierr)

 call mysecond(tdiag0)
 call eigval(elec_st,pot,u_pot,solver,parallel,ipr,ierr)
 ! Report diagonlization timing.
 call mysecond(tdiag1)
 write(9,*)
 write(9,42) tdiag1-tdiag0, tdiag_sum
 write(9,*)
 call myflush(9)
 tdiag_sum = tdiag_sum + (tdiag1-tdiag0)
 if (parallel%iammaster) then
    write(7,*)
    write(7,42) tdiag1-tdiag0, tdiag_sum
    write(7,*)
    call myflush(7)
 endif
42   format(' Diagonalization time [sec] :',1x,f10.2, &
          ',', 4x,'  tdiag_sum =', f12.2)
#ifdef MPI
 i = size(elec_st%irep)
 call MPI_ALLREDUCE(elec_st%irep,i,1,MPI_INTEGER,MPI_MAX, &
      parallel%comm,mpinfo)
#endif

 ! Finally call flevel for occupations
 if (parallel%iammaster) call flevel(elec_st,pbc%latt_vec,nscf,ierr)

end subroutine calc_nscf 

  
