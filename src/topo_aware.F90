
!===============================================================
!
!  Copyright (C) 2005 Finite Difference Research Group
!  This file is part of parsec, http://www.ices.utexas.edu/parsec/
!
!  MPI_DIST_GRAPH_CREATE_ADJACENT is a scalable communicator
!  That would evolve to be especially good for PARSEC at the moment
!  
!  This subroutine will hopefully become something productive very soon
!  Currently does the following:
!
!  1. Creates a comm_dist_graph communicator 
!     based on local output from comm_neigh.F90, 
!  2. updates local communication data in parallel, to be used with
!     either neighborhood collectives (MPI3) or plain old isend/irecv
!
!
!
!  In the future:
!
!  1. Perhaps the send/recv buffers can be allocated here?
!
!  2. In the future this subroutine should be called before any data 
!     is distributed to the nodes, so reordering of the
!     ranks can be made based on the network topology 
!
! 
!===============================================================!
subroutine topo_aware(parallel) !,experimental)

#ifdef MPI
    use mpi
#endif

    use parallel_data_module

    implicit none
    !
    !   I/O variables
    !
    type (parallel_data), intent(inout) :: parallel
!    logical, intent(in) :: experimental !whether to actually change something in how parsec works
    !
    !   
    !
    !   Work vars
    !
    integer :: graph_count, nnodes, inode, ipack, irank0,irank1
    integer, dimension(parallel%group_size) :: nelem_in,nelem_out
    !type(MPI_Info) :: info !supposedly has some hints on the graph.
    integer :: info
    logical :: reorder   !wether to reorder ranks in new comm
    integer :: mpierr    !normal fortran mpi error int
    !
    !type(MPI_Comm) :: comm_dist_graph !the new communicator. 
    integer :: comm_dist_graph, comm_dist_graph_reorder !the new communicator. 
    !
    !   Topology data
    !
    integer :: graph_degree !size of arrays needed to describe local communication
    !   Actual local communication data
    integer :: max_in_degree,max_out_degree
    integer :: sources(0:parallel%countcomm-1)
    integer :: source_weights(0:parallel%countcomm-1)
    integer :: destinations(0:parallel%countcomm-1)
    integer :: destination_weights(0:parallel%countcomm-1)
    ! integer arrays with the j-th element containing the number
    ! of elements to send to neighbor j
    integer :: sendcounts(0:parallel%countcomm-1)
    integer :: recvcounts(0:parallel%countcomm-1)
    ! integer arrays with the j-th element containing the displacement to
    ! send/recvbuff for the j-th neighbor data
    integer :: sdispls(0:parallel%countcomm-1)
    integer :: rdispls(0:parallel%countcomm-1)
    !
    ! verification
    !
    integer :: source_check
    integer :: destination_check
    logical :: weight_check

! -------------------------------------------------------------

#ifdef MPI

    graph_degree = parallel%countcomm
    nnodes = parallel%group_size

    do inode = 0, nnodes-1
        nelem_in(inode+1) = parallel%jsendp(inode+1)-parallel%jsendp(inode)
        nelem_out(inode+1) = parallel%irecvp(inode+1)-parallel%irecvp(inode)
    enddo
    ! write(9,*) 'this is what I have right now:'
    ! write(9,*) 'nelem_in', nelem_in
    ! write(9,*) 'nelem_out', nelem_out

    graph_count=size(nelem_in)
!    write(9,*) 'size(nelems)', graph_count

    graph_count=0
    do inode = 0,nnodes-1
        if (nelem_in(inode+1) /= 0) then
           sources(graph_count)=inode
           source_weights(graph_count)=nelem_in(inode+1)
           graph_count= graph_count + 1
        endif
    enddo

    graph_count=0
    do inode = 0,nnodes-1
        if (nelem_out(inode+1) /= 0) then
           destinations(graph_count)=inode
           destination_weights(graph_count)=nelem_out(inode+1)
           graph_count= graph_count + 1
        endif
    enddo

    reorder = .FALSE.
! call MPI_INFO_CREATE(info, mpierr)
info=MPI_INFO_NULL
 call MPI_DIST_GRAPH_CREATE_ADJACENT(parallel%group_comm,graph_degree,sources,source_weights,&
                                                         graph_degree,destinations,destination_weights, &
                                                         info,reorder, comm_dist_graph, mpierr )
 call MPI_DIST_GRAPH_NEIGHBORS_COUNT(comm_dist_graph,source_check,destination_check,weight_check, mpierr)

 write(9,*) 'topo_aware: I am now aware that I communicate with',source_check, 'PEs'
 if (source_check-destination_check /= 0) then
    write(9,*) 'how come you are sending and receiving from different number of PEs?'
    write(9,*) 'This is VERY WRONG.'
    ! quit program 
 endif
 write(9,*) ''

 

 !if (experimental) then
!call MPI_TOPO_TEST(parallel%group_comm,inode,mpinfo)
!write(9,*) 'topo_test returns:',inode
!actually change stuff for communication:
!we call the graph_neighbors to make sure mpi got the message about the communication pattern
     parallel%group_comm_topo = comm_dist_graph
     write(9,*) ' CAUTION! GROUP COMMUNICATOR FOR MATVEC IS THE NEIGHBORHOOD GRAPH COMMUNICATOR!'
     max_in_degree=parallel%countcomm
     max_out_degree=parallel%countcomm
call MPI_DIST_GRAPH_NEIGHBORS(parallel%group_comm_topo, max_in_degree, sources, source_weights,&
        max_out_degree, destinations, destination_weights, mpierr)

    do inode = 0,max_in_degree-1
        sendcounts(inode) = parallel%jsendp(sources(inode)+1)-parallel%jsendp(sources(inode))
    enddo

if (parallel%countcomm /= 0) sdispls(0) = 0
    do ipack=1,max_in_degree-1
        sdispls(ipack)=sdispls(ipack-1)+sendcounts(ipack-1)
    enddo

    do inode = 0,max_out_degree-1
        recvcounts(inode) = parallel%irecvp(destinations(inode)+1)-parallel%irecvp(destinations(inode))
    enddo

if (parallel%countcomm /= 0) rdispls(0)=0
    do ipack=1,max_out_degree-1
        rdispls(ipack)=rdispls(ipack-1)+recvcounts(ipack-1)
    enddo

allocate(parallel%sendcounts(0:parallel%countcomm-1))
allocate(parallel%sources(0:parallel%countcomm-1))
allocate(parallel%sdispls(0:parallel%countcomm-1))
!allocate(parallel%source_weights(0:parallel%countcomm-1))
parallel%sendcounts=sendcounts
parallel%sources=sources
parallel%sdispls=sdispls
!parallel%source_weights=source_weights

allocate(parallel%recvcounts(0:parallel%countcomm-1))
allocate(parallel%destinations(0:parallel%countcomm-1))
allocate(parallel%rdispls(0:parallel%countcomm-1))
!allocate(parallel%destination_weights(0:parallel%countcomm-1))
parallel%recvcounts=recvcounts
parallel%destinations=destinations
parallel%rdispls=rdispls
!parallel%destination_weights=destination_weights

         write(9,*) ' NEIGHBORHOOD COMMUNICATION DATA UPDATED FOR EACH PE '
 ! tests:
 call MPI_DIST_GRAPH_CREATE_ADJACENT(parallel%group_comm,graph_degree,sources,source_weights,&
                                                         graph_degree,destinations,destination_weights, &
                                                         info,.TRUE., comm_dist_graph_reorder, mpierr )
 call MPI_COMM_RANK(comm_dist_graph, irank0, mpierr)
 call MPI_COMM_RANK(comm_dist_graph_reorder, irank1, mpierr)
 call MPI_COMM_FREE(comm_dist_graph_reorder,mpierr)
 if (irank0 /= irank1) then
    write(9,*) 'topo_ware: I am PE rank',irank0,'but MPI would like me to be rank',irank1,'!!'
 else
!    write(9,*) 'topo_ware: I am PE rank',irank0,'both in the group comm and the topo comm'
 endif

     if (parallel%iammaster) then
         write(7,*) ' GROUP COMMUNICATOR FOR MATVEC IS THE NEIGHBORHOOD GRAPH COMMUNICATOR!'
#ifdef NEIGHCOLLETIVE
         write(7,*) ' WARNING: MATVEC NEIGHBORHOOD COLLECTIVES ARE ACTIVATED '
#endif
#ifdef IALLREDUCE
         write(7,*) ' WARNING: MATVEC NON-BLOCKING ALLREDUCE IS ACTIVATED '
#endif
     endif
! endif

#endif

end subroutine topo_aware
