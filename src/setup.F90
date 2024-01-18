!===============================================================
!
!  Copyright (C) 2005 Finite Difference Research Group
!  This file is part of parsec, http://www.ices.utexas.edu/parsec/
!
!  Setup a uniform, real space grid, composed of a boundary sphere
!  inside a rectangular box. 
!
!  NOTE:
!  Other geometries (e.g., an individual sphere around each atom)
!  may be attempted. They are not implemented in the current
!  version. Their addition may require adding additional tests to
!  isolat.f (see comment there).
!
!---------------------------------------------------------------
subroutine setup(grid,pbc,rsymm,parallel,nspin,kpnum,ierr)

  use constants
  use grid_module
  use pbc_module
  use symmetry_module
  use parallel_data_module

#ifdef MPI
  !  include mpi definitions
  use mpi
#endif
  implicit none
  !
  !  Input/Output variables:
  !
  !  grid related data
  type (grid_data), intent(inout) :: grid
  !  periodic boundary conditions data
  type (pbc_data), intent(inout) :: pbc
  !  symmetry operations in reduced (Abelian) group
  type (symmetry), intent(in) :: rsymm
  !  parallel computation related data
  type (parallel_data), intent(inout) :: parallel
  !  number of spins
  integer, intent(in) :: nspin
  !  number of k-points (or 1, if there are no k-points)
  integer, intent(in) :: kpnum
  !  error flag
  integer, intent(out) :: ierr
  !
  !  Work variables:
  !
  !  number of neighbors used for derivation (on one side)
  integer nord2
  !  temporary neighbor indexing
  integer neib
  !  counters
  integer ii,j,k,indx(12),numb
  !  exit code for mpi calls
  integer mpinfo
  !  number of procs available, node index
  integer nnodes,inode
  !  communicator
  integer comm
  !  rank of the master in the above communicator
  integer masterid
  integer  irows(0:parallel%procs_num)
  !  tag for the messages, it is incremented at each send/recv
  integer msgtype
 
  integer index_arr
  integer num_p,neibs_num
  integer icount,nelem,mxelem,alcstat,ishell

  !  buffer arrays
  integer, allocatable :: tmp2(:,:,:)
#ifdef MPI
  integer status(MPI_STATUS_SIZE)
#endif

  !---------------------------------------------------------------

  !
  !  Create groups of processors.
  !
#ifdef MPI
  call MPI_BCAST(parallel%groups_num, 1,MPI_INTEGER, &
       parallel%masterid,parallel%comm,mpinfo)
#endif
  call create_group_layout(parallel,rsymm%ntrans*nspin*kpnum)

  nnodes = parallel%group_size
  masterid = parallel%group_master
  comm = parallel%group_comm
  nord2 = grid%norder
  msgtype = 0
  !  Create communicator arrays.
  !
  allocate (parallel%irows(0:nnodes))

  allocate(parallel%ip1(0:nnodes))
  allocate(parallel%jp1(0:nnodes))
  allocate(parallel%irecvp(0:nnodes))
  allocate(parallel%jsendp(0:nnodes))

#ifdef MPI
  call MPI_BCAST(grid%shift, 3,MPI_DOUBLE_PRECISION, &
       parallel%masterid,parallel%comm,mpinfo)
  call MPI_BCAST(grid%n1, 3,MPI_INTEGER, &
       parallel%masterid,parallel%comm,mpinfo)
  call MPI_BCAST(grid%n2, 3,MPI_INTEGER, &
       parallel%masterid,parallel%comm,mpinfo)
  call MPI_BCAST(grid%n3, 3,MPI_INTEGER, &
       parallel%masterid,parallel%comm,mpinfo)
#endif
  !
  !  Construct regular grid, irreducible wedge, and partition.
  !
  if (parallel%iamgmaster) &
       call grid_partition(grid,pbc,rsymm,parallel,ierr)
  
#ifdef MPI
  !
  !  Broadcast Partitioning Information and PBC parameters
  !
  call MPI_BCAST(parallel%ndim, 1, MPI_INTEGER, masterid,comm,mpinfo)
  call MPI_BCAST(parallel%ldn, 1, MPI_INTEGER, masterid,comm,mpinfo)
  call MPI_BCAST(parallel%nwedge, 1, MPI_INTEGER, masterid,comm,mpinfo)
  call MPI_BCAST(parallel%irows, nnodes+1, &
       MPI_INTEGER,masterid, comm,mpinfo)
  call MPI_BCAST(pbc%mx, 1, MPI_INTEGER, masterid,comm,mpinfo)
  call MPI_BCAST(pbc%my, 1, MPI_INTEGER, masterid,comm,mpinfo)
  call MPI_BCAST(pbc%mz, 1, MPI_INTEGER, masterid,comm,mpinfo)
  call MPI_BCAST(pbc%shift, 3, MPI_DOUBLE_PRECISION, masterid,comm,mpinfo)
  call MPI_BCAST(ierr, 1, MPI_INTEGER, masterid,comm,mpinfo)
  !  Check if ierr is not zero in any group
  ii = abs(ierr)
  call MPI_ALLREDUCE(ii,ierr,1,MPI_INTEGER,MPI_MAX, &
       parallel%comm,mpinfo)
  ! broadcast grid neighbors dimensions
  call MPI_BCAST(grid%neibs_num, 1, MPI_INTEGER, masterid,comm,mpinfo)
  call MPI_BCAST(grid%nxyz, 1, MPI_INTEGER, masterid, comm, mpinfo)
#endif
  !
  !  Bail out if the irreducible wedge is incorrect
  !
  if (ierr /= 0) then
     if (parallel%iammaster) write(7,*) 'Error in grid_partition. Stop!'
     return
  endif

  irows(0:nnodes) = parallel%irows

  allocate(parallel%pint(parallel%ldn),stat=alcstat)
  call alccheck('parallel%pint',parallel%ldn,alcstat)
  allocate(parallel%senrows(parallel%nwedge),stat=alcstat)
  call alccheck('parallel%senrows',parallel%nwedge,alcstat)

  neibs_num=6+2*grid%lap_dir_num

  allocate(parallel%neibs(neibs_num*grid%norder,parallel%ldn), &
       stat=alcstat)
  call alccheck('parallel%neibs',neibs_num*grid%norder*parallel%ldn, &
       alcstat)
  parallel%neibs = 0
  allocate(parallel%tneibs(neibs_num*grid%norder,parallel%ldn), &
       stat=alcstat)
  call alccheck('parallel%tneibs', &
       neibs_num*grid%norder*parallel%ldn,alcstat)
  parallel%tneibs = 0

  if (.not. parallel%iamgmaster) then
     grid%nwedge = parallel%nwedge
     call grid_set_wedge(grid,nnodes)
     call grid_set_ist(grid,nnodes)
     if(grid%hartree_neibs_flag) &
        call grid_set_ext_neibs(grid, parallel%procs_num)
  endif

#ifdef MPI
  call MPI_BCAST(grid%kx, grid%nwedge, MPI_INTEGER, masterid,comm,mpinfo)
  call MPI_BCAST(grid%ky, grid%nwedge, MPI_INTEGER, masterid,comm,mpinfo)
  call MPI_BCAST(grid%kz, grid%nwedge, MPI_INTEGER, masterid,comm,mpinfo)
  if(grid%hartree_neibs_flag) then
  call MPI_BCAST(grid%neibs_fx, grid%neibs_num, MPI_INTEGER, masterid, comm, mpinfo) 
  call MPI_BCAST(grid%neibs_fy, grid%neibs_num, MPI_INTEGER, masterid, comm, mpinfo) 
  call MPI_BCAST(grid%neibs_fz, grid%neibs_num, MPI_INTEGER, masterid, comm, mpinfo) 
  endif
#endif
  !
  !  For each point in the sphere, make an array which points to the
  !  position of all neighbors to be used for derivation. There are
  !  3*norder of those for each point.
  !
  allocate(tmp2(neibs_num*grid%norder,parallel%ldn,2),stat=alcstat)
  call alccheck('tmp2',neibs_num*grid%norder*parallel%ldn*2,alcstat)

!  write(9,*) 'after neibs allocation'

  index_arr = 0
  do inode = 0,nnodes-1
     !  In parallel environment, send lists of neighbor points to the
     !  PEs that are going to handle them.
     msgtype = msgtype + 1
     num_p = irows(inode+1) - irows(inode)
     tmp2 = 0
     if (parallel%iamgmaster) then
        ii = 0
        do numb = irows(inode), irows(inode+1)-1
           index_arr = index_arr+1
           ii = ii + 1
           neib = 0
           do ishell = 1,nord2

  !  Doing first the main axes neighbors

              indx(1) = grid%indexg( &
                   grid%kx(numb)-ishell,grid%ky(numb),grid%kz(numb))
              indx(2) = grid%indexg( &
                   grid%kx(numb)+ishell,grid%ky(numb),grid%kz(numb))
              indx(3) = grid%indexg( &
                   grid%kx(numb),grid%ky(numb)-ishell,grid%kz(numb))
              indx(4) = grid%indexg( &
                   grid%kx(numb),grid%ky(numb)+ishell,grid%kz(numb))
              indx(5) = grid%indexg( &
                   grid%kx(numb),grid%ky(numb),grid%kz(numb)-ishell)
              indx(6) = grid%indexg( &
                   grid%kx(numb),grid%ky(numb),grid%kz(numb)+ishell)
              do j = 1, 6
                 tmp2(neib+j,ii,1) = grid%rindex(indx(j))
                 tmp2(neib+j,ii,2) = grid%rtrans(indx(j))
              enddo

              if (grid%lap_dir_num <= 0) then
                 neib = neib + neibs_num
                 cycle
              endif

              j=7 ! set the place for next copy of special directions.

  ! The following lines do the elements of laplacian
  ! coefficients that are not positioned on the main axes.

  ! Neighbors for the uv plane.
  ! changing of j is done if there is an active new laplacian
  ! direction. if the direction is not active (lap_dir(i) is 0) 
  ! then the code continues to the next direction without 
  ! advancing j!

              do k=1,3 
               if(grid%lap_dir(k) /= 0) then
                 indx(j) = grid%indexg( &
                      grid%kx(numb)-ishell*grid%lap_neig(1,k), &
                      grid%ky(numb)-ishell*grid%lap_neig(2,k), &
                      grid%kz(numb)-ishell*grid%lap_neig(3,k))
                 indx(j+1) = grid%indexg( &
                      grid%kx(numb)+ishell*grid%lap_neig(1,k), &
                      grid%ky(numb)+ishell*grid%lap_neig(2,k), &
                      grid%kz(numb)+ishell*grid%lap_neig(3,k))
                  
                 tmp2(neib+j,ii,1)=grid%rindex(indx(j))
                 tmp2(neib+j,ii,2)=grid%rtrans(indx(j))
                 j=j+1
                 tmp2(neib+j,ii,1)=grid%rindex(indx(j))
                 tmp2(neib+j,ii,2)=grid%rtrans(indx(j))
                 j=j+1
              end if
             enddo

              neib = neib + neibs_num
           enddo
        enddo
        if (parallel%group_iam == inode) then
           parallel%neibs = tmp2(:,:,1)
           parallel%tneibs = tmp2(:,:,2)
        else
#ifdef MPI
           call MPI_SEND(tmp2, neibs_num*grid%norder*parallel%ldn*2, &
                MPI_INTEGER, inode,msgtype, comm, mpinfo)
#endif
           write(9,*) ' neibs sent to node ',inode
        endif
     else
        if (parallel%group_iam == inode) then
#ifdef MPI
           call MPI_RECV(tmp2,neibs_num*grid%norder*parallel%ldn*2 &
                ,MPI_INTEGER,masterid,msgtype, comm,status,mpinfo)
           parallel%neibs = tmp2(:,:,1)
           parallel%tneibs = tmp2(:,:,2)
#endif
           write(9,*) 'RECEIVED NEIBS'
        endif
     endif
  enddo
#ifdef MPI
  call MPI_Barrier(comm,mpinfo)
#endif
  deallocate(tmp2)

  !  Master PE broadcasts indexw
#ifdef MPI
  call MPI_BCAST(grid%nxmax,1,MPI_INTEGER,masterid,comm,mpinfo)
  call MPI_BCAST(grid%nxmin,1,MPI_INTEGER,masterid,comm,mpinfo)
  call MPI_BCAST(grid%nymax,1,MPI_INTEGER,masterid,comm,mpinfo)
  call MPI_BCAST(grid%nymin,1,MPI_INTEGER,masterid,comm,mpinfo)
  call MPI_BCAST(grid%nzmax,1,MPI_INTEGER,masterid,comm,mpinfo)
  call MPI_BCAST(grid%nzmin,1,MPI_INTEGER,masterid,comm,mpinfo)
#endif
  nelem = (grid%nxmax-grid%nxmin+1)*(grid%nymax-grid%nymin+1)* &
       (grid%nzmax-grid%nzmin+1)
  if (.not. parallel%iamgmaster) then
     allocate (grid%indexw (grid%nxmin:grid%nxmax, &
          grid%nymin:grid%nymax,grid%nzmin:grid%nzmax) &
          ,stat=alcstat)
     call alccheck('indexw',nelem, alcstat)
     grid%indexw = 0
     if(grid%hartree_neibs_flag) then
      allocate (grid%neibs_index (grid%nxmin:grid%nxmax, &
                                  grid%nymin:grid%nymax, &
                                  grid%nzmin:grid%nzmax), &
                                  stat=alcstat)
      call alccheck('neibs_index',grid%nxyz,alcstat)
     endif

  endif
#ifdef MPI
  call MPI_BCAST(grid%n1,1,MPI_INTEGER,masterid,comm,mpinfo)
  call MPI_BCAST(grid%n2,1,MPI_INTEGER,masterid,comm,mpinfo)
  call MPI_BCAST(grid%n3,1,MPI_INTEGER,masterid,comm,mpinfo)
  call MPI_BCAST(grid%ndim,1,MPI_INTEGER,masterid,comm,mpinfo)
  call MPI_BCAST(grid%nwedge,1,MPI_INTEGER,masterid,comm,mpinfo)
  call MPI_BCAST(grid%indexw,nelem,MPI_INTEGER,masterid,comm,mpinfo)
  if(grid%hartree_neibs_flag) &
   call MPI_BCAST(grid%neibs_index,nelem,MPI_INTEGER,masterid,comm,mpinfo)
#endif
  if (.not. parallel%iamgmaster) &
       call grid_set_ndim(grid%ndim,grid)
#ifdef MPI
  call MPI_BCAST(grid%rindex,grid%ndim+1,MPI_INTEGER,masterid,comm,mpinfo)
  call MPI_BCAST(grid%rtrans,grid%ndim+1,MPI_INTEGER,masterid,comm,mpinfo)
  call MPI_BCAST(grid%fx,grid%ndim,MPI_INTEGER,masterid,comm,mpinfo)
  call MPI_BCAST(grid%fy,grid%ndim,MPI_INTEGER,masterid,comm,mpinfo)
  call MPI_BCAST(grid%fz,grid%ndim,MPI_INTEGER,masterid,comm,mpinfo)
#endif

  !  define local dimension mydim
  parallel%mydim = irows(parallel%group_iam+1) - irows(parallel%group_iam)
  write(9,*) 'The local dimension is :', parallel%mydim

  parallel%lap_dir_num = grid%lap_dir_num
#ifdef MPI
      call MPI_Barrier(parallel%comm,mpinfo)
#endif

  !  Find boundary information for exchange
  call comm_neigh(parallel,grid%norder)

  write(9,*) 'Information to be sent to neighbors:'
  icount = 0
  mxelem = 0
  do inode = 0, nnodes - 1
     nelem = parallel%jsendp(inode+1)-parallel%jsendp(inode)
     write(9,*) 'send PE ',inode,' #rows:',nelem
     mxelem = mxelem + nelem
     if (nelem/=0) icount = icount +1
  enddo

  parallel%maxcomm   = mxelem
  parallel%countcomm = icount
  parallel%maxcomm2  = parallel%irecvp(nnodes)-parallel%irecvp(0)
  write(9,*) 'I communicate with ',icount, &
       ' processors, a total of ',parallel%maxcomm,' elements'
  write(9,*) 'Also, I receive from these procs: ', &
       parallel%maxcomm2, ' elements'

       if (parallel%maxcomm /= parallel%maxcomm2) then
           write(9,*) ''
           write(9,*) 'WARNING, ASSYMETRIC COMM PATTERN ENCOUNTERED'
       endif
#ifdef MPI
  !  create the graph communicator and update local comm data
  write(9,*) ''
  call topo_aware(parallel) !,grid%experimental)
  !dunno where to put it
  ! this whole hack smells. need to reorder parsec structures anyhoo
#endif
  parallel%max_lap_buffers=grid%max_lap_buffers
end subroutine setup
!===============================================================
