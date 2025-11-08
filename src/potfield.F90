!=============================================================================!
! This subroutine reads pot.dat and adds the potential field
! Part of the AFM project
!                             Jan. 2013, Minjung Kim (minjung@ices.utexas.edu)
!             Modified        Jul. 2015, Yuki Sakai (yuki@ices.utexas.edu)
!=============================================================================!

subroutine potfield(grid,pbc,parallel,v_exi,v_exh,rho0,oldk,f_name,iunit,ierr)

  use constants
  use grid_module
  use pbc_module
  use potential_module
  use parallel_data_module
  
  implicit none
#ifdef MPI
  include 'mpif.h'
#endif
  ! 1. Input/Output variables:
  !  Grid related data
  type (grid_data), intent(in) :: grid
  !  Periodic boundary conditions data
  type (pbc_data), intent(in) :: pbc
  !  Parallel computation related data
  type (parallel_data), intent(inout) :: parallel
  !  External potential (hartree)
  real(dp), intent(inout) :: v_exh(parallel%mydim)
  !  External potential (ion)
  real(dp), intent(inout) :: v_exi(parallel%mydim)
  !  External potential (xc)
  !real(dp), intent(inout) :: v_exxc(parallel%mydim)
  !  density of the sample material
  real(dp), intent(inout) :: rho0(parallel%mydim)
  !  sample kinetic energy
  real(dp), intent(out) :: oldk
  !  Input File Name
  character(len=50), intent(in)  ::  f_name
  !  Input file flag
  integer  :: iunit
  !  Error flag
  integer  :: ierr
  !parallel
  integer, allocatable :: nsize(:), displs(:)

  ! 2. Work variables:
  integer  ::  natm, ia
  real(dp),dimension(:),allocatable ::  xa, ya, za
  integer  ::  stype(3)
  !  flags for type of system
  !  stype(1) = nspin
  !  stype(2) = 0 : real wavefunctions
  !  stype(2) = 1 : complex wavefunctions
  !  stype(2) = 2 : spin-orbit (complex doubled wavefunctions)
  !  stype(3) = 0 : confined system (cluster)
  !  stype(3) = 1 : system periodic on all directions (bulk)
  !  stype(3) = 2 : system periodic on x direction only (wire)
  !  stype(3) = 3 : system periodic on x and y directions (slab)
  real(dp)  ::  step(3), box_size(3), latt_vec(3,3), rmax, shift(3)
  integer   ::  ndim, nwedge, ntrans
  integer, allocatable  ::  chi(:,:),rmtrx(:,:,:)
  real(dp),allocatable  ::  trans(:,:,:), tnp(:)
  real(dp)  ::  alatt(3,3)
  real(dp)  ::  invlat(3,3)
  integer,allocatable :: kx(:), ky(:), kz(:)
  real(dp),allocatable  ::  vion(:), vhart(:), vxc(:), rhoc(:)
  real(dp) :: ke

  ! tmp saves the electron density of the sample materials
  real(dp), allocatable  ::  tmp(:), tmph(:), tmpi(:), tmpxc(:), srho(:), tmpc(:)
  real(dp), allocatable  ::  ptmp(:), ptmph(:), ptmpi(:), ptmpxc(:), psrho(:), ptmpc(:) 
  integer,allocatable :: pkx(:), pky(:), pkz(:)
  real(dp) :: t1, t2, t3, t4, t5, t6, tmpns, tmpnw

  integer :: i, j, k, locpp,prev

  ! grid index
  integer :: nmin(3), nmax(3), nzmin
  ! allocation check
  integer alcstat

  integer mpinfo
  integer status(MPI_STATUS_SIZE)


if (parallel%iammaster) then
  allocate(tmp(grid%nwedge),stat=alcstat)
  call alccheck('tmp',grid%nwedge,alcstat)
  tmp(:) = zero
  allocate(tmph(grid%nwedge),stat=alcstat)
  call alccheck('tmph',grid%nwedge,alcstat)
  tmph(:) = zero
  allocate(tmpi(grid%nwedge),stat=alcstat)
  call alccheck('tmpi',grid%nwedge,alcstat)
  tmpi(:) = zero
  allocate(tmpxc(grid%nwedge),stat=alcstat)
  call alccheck('tmpxc',grid%nwedge,alcstat)
  tmpxc(:) = zero
  allocate(tmpc(grid%nwedge),stat=alcstat)
  call alccheck('tmpc',grid%nwedge,alcstat)
  tmpc(:) = zero
  allocate(rhoc(grid%nwedge),stat=alcstat)
  call alccheck('rhoc',grid%nwedge,alcstat)
  rhoc(:) = zero
              
  !--------  READ pot.dat  ---------!
  open(iunit,file=trim(f_name),form='unformatted',status='old')
                    
  read(iunit) oldk
  read(iunit) natm
  allocate( xa(natm), ya(natm), za(natm) )
  do ia = 1, natm
     read(iunit) xa(ia), ya(ia), za(ia)
  enddo
  read(iunit) ( stype(i), i=1,3 ) ! system flag
                   
  if ( stype(3) .eq. 0 ) then
     read(iunit) step(1), rmax
     step(2:3) = step(1)
  else
     read(iunit) ( step(i), i=1,3 )
     read(iunit) ( box_size(i), i=1,3 )
     read(iunit) ( (latt_vec(i,j), i=1,3) , j=1,3 ) 
  endif
     

  read(iunit) ndim
  read(iunit) ( shift(i), i=1,3 )
  read(iunit) nwedge, ntrans
endif

call MPI_BCAST(nwedge,1,MPI_INTEGER, &
    parallel%masterid,parallel%comm,mpinfo)
call MPI_BCAST(ntrans,1,MPI_INTEGER, &
    parallel%masterid,parallel%comm,mpinfo)
allocate( rmtrx(3,3,ntrans),trans(3,3,ntrans),chi(ntrans,ntrans),&
            kx(nwedge),ky(nwedge),kz(nwedge) )

if(parallel%iammaster) then
  read(iunit) rmtrx
  read(iunit) trans
  read(iunit) !tnp
  read(iunit) alatt, invlat
  read(iunit) ( (chi(i,j),i=1,ntrans), j=1,ntrans )
  read(iunit) ( (kx(i), ky(i), kz(i)), i=1,nwedge )
endif

allocate(vion(nwedge),vhart(nwedge),vxc(nwedge),srho(nwedge))

if(parallel%iammaster) then
  read(iunit) ( vion(i), i=1,nwedge )
  read(iunit) ( vhart(i), i=1,nwedge )
  read(iunit) ( vxc(i), i=1,nwedge )
  read(iunit) ( srho(i), i=1, nwedge)
!  read(iunit) ( rhoc(i), i=1, nwedge)
  rhoc(:) = 0.0d0

  close(iunit)

  ierr = 0 

  !1. Check if the grid spacing matches 
  !At this moment, interpolation scheme is not applied
  do i = 1, 3
     if ( grid%step(i) .ne. step(i) ) then
        write(7,*) 'Grid space does not match. Program exits'
        write(7,*) 'host step'
        write(7,*) grid%step
        write(7,*) 'field step'
        write(7,*) step
        stop 601 
     endif
  enddo

  !2. Doublecheck
  if(parallel%iammaster) then
     if (nwedge .ne. grid%nwedge) then
        write(7,*) 'The size of hamiltonian does not match'
        write(7,*) 'Current',grid%nwedge
        write(7,*) 'Read',nwedge
        stop 602
     end if
  end if

  ! 3. Check if the field is equal or smaller than the box
  if ( grid%nxmax .lt. nmax(1) ) then
     write(7,*) 'Input field is bigger in x axis. At this moment, it is not supported'
     write(7,*) 'Program exits'
     stop 603
  elseif ( grid%nymax .lt. nmax(2) ) then
     write(7,*) 'Input field is bigger in y axis. At this moment, it is not supported'
     write(7,*) 'Program exits'
     stop 603
  elseif ( grid%nzmax .lt. nmax(3) ) then
     write(7,*) 'Input field is bigger in z axis. At this moment, it is not supported'
     write(7,*) 'Program exits'
     stop 603
  endif

  nmin(1) = minval( kx ) 
  nmax(1) = maxval( kx ) 
  nmin(2) = minval( ky ) 
  nmax(2) = maxval( ky ) 
  nmin(3) = minval( kz ) 
  nmax(3) = maxval( kz )

endif

allocate(nsize(parallel%procs_num))
allocate(displs(parallel%procs_num))

!Computing array partition
if(parallel%iammaster) then
  tmpns = nwedge / parallel%procs_num
  tmpnw = nwedge
  do i = 1, parallel%procs_num
     tmpnw = tmpnw - tmpns
     if(i==parallel%procs_num .AND. tmpnw .lt. 0) then
        nsize(i) = tmpnw
     elseif(i==parallel%procs_num .AND. tmpnw .gt. 0) then
        nsize(i) = tmpns + tmpnw
     else
        nsize(i) = tmpns
     end if 
     displs(i) = tmpns * (i - 1)
  end do
endif ! MASTER Node Done

  ! 4. Assign potential values
  !Parallelization stuff
   locpp = parallel%iam + 1
#ifdef MPI
  call MPI_Barrier(parallel%comm,mpinfo)
  call MPI_BCAST(nsize,parallel%procs_num,MPI_INTEGER, &
       parallel%masterid,parallel%comm,mpinfo)
  call MPI_BCAST(displs,parallel%procs_num,MPI_INTEGER, &
       parallel%masterid,parallel%comm,mpinfo)
#endif
  allocate(ptmp(nsize(locpp)),stat=alcstat)
  call alccheck('ptmp',nsize(locpp),alcstat)
  ptmp(:) = zero
  allocate(ptmpi(nsize(locpp)),stat=alcstat)
  call alccheck('ptmpi',nsize(locpp),alcstat)
  ptmpi(:) = zero
  allocate(ptmph(nsize(locpp)),stat=alcstat)
  call alccheck('ptmph',nsize(locpp),alcstat)
  ptmph(:) = zero
  allocate(ptmpxc(nsize(locpp)),stat=alcstat)
  call alccheck('ptmpxc',nsize(locpp),alcstat)
  ptmpxc(:) = zero
  allocate(ptmpc(nsize(locpp)),stat=alcstat)
  call alccheck('ptmpc',nsize(locpp),alcstat)
  ptmpc(:) = zero
  allocate(pkx(nsize(locpp)),stat=alcstat)
  call alccheck('pkx',nsize(locpp),alcstat)
  allocate(pky(nsize(locpp)),stat=alcstat)
  call alccheck('pky',nsize(locpp),alcstat)
  allocate(pkz(nsize(locpp)),stat=alcstat)
  call alccheck('pkz',nsize(locpp),alcstat)

#ifdef MPI
  call MPI_Barrier(parallel%comm,mpinfo)
  call MPI_BCAST(vion ,nwedge,MPI_DOUBLE_PRECISION, &
       parallel%masterid,parallel%comm,mpinfo)
  call MPI_BCAST(vhart,nwedge,MPI_DOUBLE_PRECISION, &
       parallel%masterid,parallel%comm,mpinfo)
  call MPI_BCAST(vxc,  nwedge,MPI_DOUBLE_PRECISION, &
       parallel%masterid,parallel%comm,mpinfo)
  call MPI_BCAST(srho, nwedge,MPI_DOUBLE_PRECISION, &
       parallel%masterid,parallel%comm,mpinfo)
  call MPI_BCAST(kx,   nwedge,MPI_INTEGER, &
       parallel%masterid,parallel%comm,mpinfo)
  call MPI_BCAST(ky,   nwedge,MPI_INTEGER, &
       parallel%masterid,parallel%comm,mpinfo)
  call MPI_BCAST(kz,   nwedge,MPI_INTEGER, &
       parallel%masterid,parallel%comm,mpinfo)

  call MPI_Barrier(parallel%comm,mpinfo)
  CALL MPI_SCATTERV(grid%kx(1),nsize,displs,MPI_INTEGER,&
       pkx(1),nsize(parallel%iam+1),MPI_INTEGER,parallel%masterid,parallel%comm,ierr)
  CALL MPI_SCATTERV(grid%ky(1),nsize,displs,MPI_INTEGER,&
       pky(1),nsize(parallel%iam+1),MPI_INTEGER,parallel%masterid,parallel%comm,ierr)
  call MPI_Barrier(parallel%comm,mpinfo)
  CALL MPI_SCATTERV(grid%kz(1),nsize,displs,MPI_INTEGER,&
       pkz(1),nsize(parallel%iam+1),MPI_INTEGER,parallel%masterid,parallel%comm,ierr)
#endif

  prev = 0

  do j = 1, nsize(locpp)
     if ( pkx(j)==kx(prev+1) .and. pky(j)==ky(prev+1) &
           .and. pkz(j)==kz(prev+1) ) then
!          .and. pkz(j)==kz(i)-(nmin(3)-nzmin) ) then 
        ptmpi(j) = vion(prev+1)
        ptmph(j) = vhart(prev+1)
 !       ptmpxc(j) = vxc(i)
        ptmp(j) = srho(prev+1)
!        ptmpc(j) = rhoc(prev+1)
        prev = prev+1
        goto 243
     endif

     do i = 1, nwedge
        if ( pkx(j)==kx(i) .and. pky(j)==ky(i) &
           .and. pkz(j)==kz(i) ) then 
!          .and. pkz(j)==kz(i)-(nmin(3)-nzmin) ) then 
           ptmpi(j) = vion(i)
           ptmph(j) = vhart(i)
 !          ptmpxc(j) = vxc(i)
           ptmp(j) = srho(i)
!           ptmpc(j) = rhoc(i)
           prev = i
           goto 243
        endif
     enddo
243 continue
  enddo
  
#ifdef MPI
  call MPI_Barrier(parallel%comm,mpinfo)
  CALL MPI_GATHERV(ptmpi(1), nsize(parallel%iam+1), MPI_DOUBLE_PRECISION, & 
       tmpi(1), nsize,displs, MPI_DOUBLE_PRECISION, parallel%masterid,parallel%comm, ierr)
  CALL MPI_GATHERV(ptmph(1), nsize(parallel%iam+1), MPI_DOUBLE_PRECISION, &
       tmph(1), nsize,displs, MPI_DOUBLE_PRECISION, parallel%masterid,parallel%comm, ierr)
  CALL MPI_GATHERV(ptmp(1),  nsize(parallel%iam+1), MPI_DOUBLE_PRECISION, & 
       tmp(1),  nsize,displs, MPI_DOUBLE_PRECISION, parallel%masterid,parallel%comm, ierr)
! CALL MPI_GATHERV(ptmpc(1),  nsize(parallel%iam+1), MPI_DOUBLE_PRECISION, & 
!      tmpc(1),  nsize,displs, MPI_DOUBLE_PRECISION, parallel%masterid,parallel%comm, ierr)
#endif

if (parallel%iammaster) then
! do i = 1, grid%nwedge
!    tmp(i) = tmp(i) + tmpc(i)
! end do
  parallel%ftmp = tmp
endif ! MASTER Node Done

! Let all processors know the external potential so that it can be added 
!     if potential_field is true
#ifdef MPI
  call MPI_Barrier(parallel%comm,mpinfo)
#endif
  call export_function(parallel,rho0)

if (parallel%iammaster) then
  parallel%ftmp = tmpxc
endif ! MASTER Node Done

!#ifdef MPI
!  call MPI_Barrier(parallel%comm,mpinfo)
!#endif
!  call export_function(parallel,v_exxc)

if (parallel%iammaster) then
  parallel%ftmp = tmpi
endif ! MASTER Node Done

#ifdef MPI
  call MPI_Barrier(parallel%comm,mpinfo)
#endif
  call export_function(parallel,v_exi)

if (parallel%iammaster) then
  parallel%ftmp = tmph
  ! Deallocate
  deallocate( rmtrx, trans, chi, kx, ky, kz, vion, vhart, vxc, tmp, tmpi, tmph, tmpxc, srho)  
endif ! MASTER Node Done

#ifdef MPI
  call MPI_Barrier(parallel%comm,mpinfo)
#endif
  call export_function(parallel,v_exh)

  deallocate(ptmp, ptmpi, ptmph, ptmpxc, pkx, pky, pkz) 
end subroutine potfield
