!===============================================================
!
!  Copyright (C) 2005 Finite Difference Research Group
!  This file is part of parsec, http://www.ices.utexas.edu/parsec/
!
!  Defines:  
!      1. parallel%neibs in local form
!      2. irecvp,jsendp,senrows indices for exchanging
!         info in matvec
!      3. ip1,jp1 for exchanging info in preconditioning
!      4. pint local permutation with inter1,inter2 interior points
!
!  Rewrites parallel%neibs with a local one with the properties:
!
!              |--> parallel%neibs(neigh,row)-ioffset for all 
!              |     neighbors "neigh" local on this processor
!              |
!  parallel%neibs(neigh,row)= --> position of neigh, (> mydim) in the boundary 
!              |    values appended at the end of the current vector.
!              |    This is the "neigh" neighbor of the local "row" 
!              |   
!              |--> ndim+1 for all unused neighbors 
!
!  The boundary values needed from this node are padded at the
!  end of a local vector in a boundary-exchange information step.
!
!                     Receive and append info from 
!             local on the node    Proc1  Proc2  Proc4
!           <-------------------> <-----> <---> <------->
!
!  There is a local IP array : IP(i) = where the info from proc i
!  should be padded at. (eg IP(0) = ndim + 1, and IP(1) = ndim + 1
!  if nothing is received from proc 0). The rows that must be sent
!  to other processors are in senrows with index for each proc JP
!  Then matrix vector multiply is performed by exchanging info
!  filling the padded area and running the old mat-vec on it.
!
!  Also IP1 and JP1 are used for the communication in the 
!  preconditioning. This uses only the first shell (neighbors up
!  to norder) for every point. Since these are numbered first in
!  the boundary we put them in the same positions as in matvec.
!  Then parallel%neibs knows how to access them. Notice that IP
!  carries the positions of where they should be put. So we only
!  keep the sizes. Since we number the first shell first, the
!  senrows will contain exactly these rows if we take JP1(*) of them.
!
!  COMMUNICATION SIZE GRAPH:
!              Recv from 
!   local on the node    Proc1  Proc2  Proc4
! <-------------------> <-----> <---> <------->  index IP  (matvec)
!                       <->     <->   <-->    sizes in IP1 (precond)
!
!                  senrows
! <----------------------------------------->
! send Proc1     Proc2    Proc3      Proc4
! <-------><---------><------><------------->    index JP  (matvec)
! <-->     <-->       <->     <----->         sizes in JP1 (precond)
!
!----------------------------------------------------------------------
!  Experimentally sort first the interior rows (ie that do not 
!  require any neighbors from other nodes). Sort second the rows
!  that do not require any first shell neighbors. These will be
!  considered interior for the preconditioning. Permutation array
!  is pint and the above numbers inter1 and inter2 respectively.
!                        pint
!            <----------------------------------------->
!            <---------> inter1
!            <---------------------------> inter2
!
!	Algorithm
!  Count how many and which rows I need from each processor. (irecvp)
!  Also count the number of my rows each of my neighbors needs (jsendp)
!  (each time irow needs info from proc node, irow will be needed
!  by proc node -symmetric matrix-. If I have not received this
!  row already, add it in irecvp for reception. For sending out, a 
!  row should be sent only once to a node and not for all the 
!  neighbors of this row that reside on that node. To guarantee
!  that, for each node a flag array is updated for each row.
!  To keep the order of reception same to the one that other procs
!  will compute ---to retain locality in the access pattern of 
!     neighbors--- senrows is used as scratch to store the order.
!
!---------------------------------------------------------------
! Modified to eliminate local stack arrays.  MGR 16feb09
!---------------------------------------------------------------
subroutine comm_neigh(parallel,norder)

  use constants
  use parallel_data_module
#ifdef MPI
  !   include mpi definitions
  use mpi
#endif

  implicit none
  !
  !  Input/Output variables:
  !
  !  parallel computation related data
  type (parallel_data), intent(inout) :: parallel
  !  number of neighbors used on one side in numerical derivative
  integer, intent(in) :: norder
  !
  !  Work variables:
  !
  integer nnodes,iam,node,mpinfo
  integer ielem,jelem,istart,jstart,irow,myrow,itemp,mydim &
       ,ioffset,ncur,ish

  logical, allocatable :: sentto(:,:)
  integer, allocatable :: received(:)

  integer, dimension(:), allocatable :: irows &
       ,isendcount,irecvcount

  integer, allocatable :: recrows(:)

  integer inter1, inter2

  integer i, alcstat
  !
  !  External functions:
  !
  integer, external :: lochome

  !---------------------------------------------------------------

  !  Initializations
  iam = parallel%group_iam
  nnodes = parallel%group_size

  allocate(irows(0:nnodes))
  allocate(sentto(parallel%mydim,0:nnodes))
  allocate(isendcount(0:nnodes))
  allocate(irecvcount(0:nnodes))
  allocate(received(parallel%nwedge))
  allocate(recrows(parallel%nwedge))

  irows(:) = parallel%irows
  ioffset = irows(iam) - 1
  mydim   = parallel%mydim

  received(:) = 0
  recrows(:) = 0
  parallel%senrows(:) = 0
  parallel%irecvp(:) = 0
  parallel%jsendp(:) = 0
  parallel%ip1(:) = 0
  parallel%jp1(:) = 0
  isendcount(:) = 0
  irecvcount(:) = 0

  sentto(:,:) = .false.
  inter1 = 0
  !
  !  * Experimental *
  !  *** notice that locality on the vec won't be great in matvec
  !  *** Eventually parallel%neibs must be turned into
  !  *** parallel%neibs(rows,neighbors)
  !  *** We will have to experiment with the effects of index pint
  !
  ncur = 0
  !  for all neighbors of this row: 
  !  2*(3+lap_dir_num)*norder (first shell first)

  do ish = 1, 2*(3+parallel%lap_dir_num)*norder
     do myrow = 1, mydim
        irow = parallel%neibs(ish,myrow)
  !  find processor # of neighbor irow
        node = lochome(irow,nnodes,irows,iam)

        if (node /= iam) then
           if (.not.sentto(myrow,node)) then
              parallel%jsendp(node) = parallel%jsendp(node) + 1
              if (ish < 7) parallel%jp1(node) = parallel%jp1(node)+1
              sentto(myrow,node) = .true.
           endif
           if (received(irow) ==  0) then
              parallel%irecvp(node) = parallel%irecvp(node) + 1
              if (ish < 7) parallel%ip1(node) = parallel%ip1(node)+1
              ncur = ncur + 1
              parallel%senrows(ncur) = irow
              received(irow) = node+1
           endif
        endif
     enddo                  ! myrow = 1, mydim
  enddo                     ! ish = 1,2*(3+parallel%lap_dir_num)*norder
  !
  !  Find where info from (irecvp) and to (jsendp) each proc starts
  !
  ielem = parallel%irecvp(0)
  jelem = parallel%jsendp(0)
  parallel%irecvp(0) = 1
  parallel%jsendp(0) = 1
  do node = 1, nnodes
     istart = parallel%irecvp(node-1) + ielem
     jstart = parallel%jsendp(node-1) + jelem
     ielem = parallel%irecvp(node)
     jelem = parallel%jsendp(node)
     parallel%irecvp(node) = istart
     parallel%jsendp(node) = jstart
  enddo
  !
  !  Use received and senrows arrays to paste together the columns 
  !  I will need. Reset received to zero wherever it has something
  !  else. Rows are according to the global permutation.
  !
  do i = 1, ncur
     irow = parallel%senrows(i)
     node = received(irow)-1
     received(irow) = 0
     recrows( parallel%irecvp(node) ) = parallel%senrows(i)
     parallel%irecvp(node) = parallel%irecvp(node) + 1
  enddo
  !
  !  Restore the ip index to starting points in the appended vector
  !  for each node. This is local now (eg irecvp(0) = 1 + mydim)
  !
  do i = nnodes-1,1,-1
     parallel%irecvp(i) = parallel%irecvp(i-1) + mydim
  enddo
  parallel%irecvp(0) = 1 + mydim
  !
  !  Go through parallel%neibs array to reset the neigh global
  !  indexing to a local. Notice the order of the loops must be the
  !  same as in previous loops over ish, myrow.
  !  IMPORTANT: Since this is not 1-to-1 mapping, the same row may
  !  be accessed as a neighbor of many rows in parallel%neibs. To
  !  remember where this row resides in the new order, the received
  !  array now records the position of the first time this row
  !  (neighbor) is encountered.
  !
  do ish = 1,2*(3+parallel%lap_dir_num)*norder
  !  for this neighbor of all rows
     do myrow = 1, mydim
        irow = parallel%neibs(ish,myrow)
        node = lochome(irow,nnodes,irows,iam)
        if (node /= iam) then
           if (received(irow) ==  0) then
              ! First time encountered. Change irecvp, received
              received(irow)  = parallel%irecvp(node)
              if(parallel%irecvp(node) > parallel%irecvp(nnodes)+mydim) write(9,*) &
                   'ipn>> ',ish,myrow,irow,parallel%irecvp(node),parallel%irecvp(nnodes)
              parallel%irecvp(node) = parallel%irecvp(node) + 1
           else
              ! This row exists already in the list
              if(received(irow) > parallel%irecvp(nnodes)+mydim .or. &
                   received(irow) < 1)write(9,*) 'exs>> ',ish &
                   ,myrow,irow,received(irow)
           endif
           parallel%neibs(ish,myrow) = received(irow)
        else
           if (irow == parallel%nwedge+1) cycle
           parallel%neibs(ish,myrow) = irow - ioffset
        endif
     enddo                  ! myrow = 1, mydim
  enddo
  !
  !  Find and permute the inter1 interior points in pint
  !
  do i = 1, mydim
     parallel%pint(i) = i
  enddo
  do myrow = 1,mydim
     do ish = 1,2*(3+parallel%lap_dir_num)*norder
        irow = parallel%neibs(ish,myrow)
        if (irow > mydim) goto 10 ! interface row
     enddo
     inter1 = inter1+1      ! interior node 
     itemp  = parallel%pint(myrow)   ! permute to the begining
     parallel%pint(myrow) = parallel%pint(inter1)
     parallel%pint(inter1)= itemp
 10  continue
  enddo
  !
  !  Find interior rows (not having first shell neighbors on other
  !  procs)
  !
  inter2 = inter1
  do myrow = inter1+1, mydim
     do ish = 1,2*(3+parallel%lap_dir_num)
        irow = parallel%neibs(ish,myrow)
        if (irow > mydim) goto 11 ! interface row
     enddo
     inter2 = inter2 + 1    ! no first shell neighs
     itemp  = parallel%pint(myrow)   ! permute this row
     parallel%pint(myrow) = parallel%pint(inter2)
     parallel%pint(inter2)= itemp
 11  continue
  enddo
  !
  !  Restore again the irecvp index to starting points, but in the 
  !  non-appended vector to be used in the alltoall
  !
  do i = nnodes-1,1,-1
     parallel%irecvp(i) = parallel%irecvp(i-1)-mydim
  enddo
  parallel%irecvp(0) = 1 
  !
  !  Communicate the information from other processors' irecvp into 
  !  the appropriate jsendp locations. This is an alltoallv step 
  !  where the j-th section of recrows (from irecvp(j) to irecvp(j+1)-1) on 
  !  processor i, goes to the j-th processor senrows between its 
  !  jsendp(i) and jsendp(i+1)-1.
  !
  do i=0,nnodes-1
     isendcount(i) = parallel%irecvp(i+1)-parallel%irecvp(i)
     irecvcount(i) = parallel%jsendp(i+1)-parallel%jsendp(i)
  enddo
  do i=0,nnodes
     parallel%irecvp(i) = parallel%irecvp(i)-1
     parallel%jsendp(i) = parallel%jsendp(i)-1
  enddo

#ifdef MPI
  call MPI_Alltoallv(recrows,isendcount, parallel%irecvp, MPI_INTEGER, &
       parallel%senrows,irecvcount, parallel%jsendp, MPI_INTEGER,parallel%group_comm,mpinfo)
  !no need for barrier, this is blocking
  !call MPI_Barrier(parallel%group_comm,mpinfo) 
#endif
  !
  !  Finally, restore irecvp and jsendp indices 
  !
  do i = 0, nnodes
     parallel%irecvp(i) = parallel%irecvp(i) + 1 + mydim
     parallel%jsendp(i) = parallel%jsendp(i) + 1
  enddo

  !  printout
  write(9,*) 'INTERIOR POINTS Mvec, PreCond (', inter1,inter2, ') '
  if (inter1 == 0) then
  write(9,*) 'WARNING, NO MATVEC INTERIOR POINTS'
  write(9,*) 'CURRENT SYSTEM SIZE/PE RATIO NOT RECOMMENDED'
  endif
  parallel%inter1=inter1
  parallel%inter2=inter2
  !  Change the senrows to local numbering (subtract the offset)
  do i = 1, parallel%jsendp(nnodes)-1
     parallel%senrows(i) = parallel%senrows(i) - ioffset
  enddo

  deallocate(received)
  deallocate(recrows)
  deallocate(irows)
  deallocate(sentto)
  deallocate(isendcount)
  deallocate(irecvcount)

end subroutine comm_neigh
!===============================================================
!
!  Returns the processor that the global row irow is on.
!
!---------------------------------------------------------------
integer function lochome(irow,nnodes,irows,iam)

  implicit none
  !
  !  Input/Output variables:
  !
  integer, intent(in) :: irow, iam, nnodes
  integer, intent(in) :: irows(0:nnodes)
  !
  !  Work variables:
  !
  integer i
  !---------------------------------------------------------------
  lochome = iam
  do i = 0, nnodes-1
     if (irows(i) <= irow .and. irow < irows(i+1)) then
        lochome = i
        exit
     endif
  enddo

end function lochome
!===============================================================
!
!  Global sum of integer array
!
!---------------------------------------------------------------
subroutine pisum(ivec,nmax,nnodes,comm)

  use constants
#ifdef MPI
  !  include mpi definitions
  use mpi
#endif
  implicit none
  !
  !  Input/Output variables:
  !
  !  dimension of array vec
  integer, intent(in) :: nmax
  !  number of procs available
  integer, intent(in) :: nnodes
  !  communicator
  integer, intent(in) :: comm
  !  The original vector whose value is partial for each PE
  integer, intent(inout) :: ivec(nmax)
  !
  !  Work variables:
  !
  !  A temporary work array used to update the vectors
  integer, allocatable :: iwork(:)
  !  exit code for mpi calls
  integer mpinfo
  !  ---------------------------------------------------------------
  if (nnodes == 1) return
#ifdef MPI
  allocate(iwork(nmax))
  call MPI_allREDUCE(ivec,iwork,nmax,MPI_INTEGER, &
       MPI_SUM,comm,mpinfo)
  ivec(1:nmax) = iwork(1:nmax)
  deallocate(iwork)
#endif
end subroutine pisum
!===============================================================
