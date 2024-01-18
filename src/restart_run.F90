!===============================================================
!
!  Copyright (C) 2005 Finite Difference Research Group
!  This file is part of parsec, http://www.ices.utexas.edu/parsec/
!
!  This subroutine reads grid, potential, wave function, energy
!  levels, etc, from a previous calculation, as stored in parsec.dat,
!  for restart purposes.
!  Extrapolation (i.e., reshaping of charge density, and wave
!  functions) is always done, although it is only relevant if
!  periodic boundary conditions are not used and the previous
!  parsec.dat refer to a calculation with smaller enclosing sphere.
!
!---------------------------------------------------------------
subroutine restart_run(elec_st,grid,pbc,parallel,solver &
     ,vold,rho,iunit,ixtrpflg,ierr)

  use constants
  use electronic_struct_module
  use grid_module
  use pbc_module
  use parallel_data_module
  use eigen_solver_module
#ifdef MPI
  use mpi
#endif
  implicit none
  !  include mpi definitions
  !
  !  Input/Output variables:
  !
  !  electronic structure
  type (electronic_struct), intent(inout) :: elec_st
  !  grid related data
  type (grid_data), intent(in) :: grid
  !  periodic boundary conditions data
  type (pbc_data), intent(in) :: pbc
  !  parallel computation related data
  type (parallel_data), intent(inout) :: parallel
  !  solver related data
  type (eigen_solver), intent(inout) :: solver

  !  distributed self-consistent potential and electron density
  !  (passed outside the structure to overcome a bug with the IBM compiler)
  real(dp), intent(out) :: vold(parallel%mydim,elec_st%nspin)
  real(dp), intent(out) :: rho(parallel%mydim,2*elec_st%nspin-1)

  !  number of output unit
  integer, intent(in) :: iunit
  !  flag indicating whether or not extrapolation to a larger sphere
  !  size is performed. It is used in the main code for deciding
  !  whether or not the Hartree and exchange-correlation potentials
  !  need to be calculated upon a restart.
  integer, intent(out) :: ixtrpflg
  !  compatibility flag: read output to parsec.dat if false; read
  !  error flag, 360 < ierr < 381
  integer, intent(out) :: ierr
  !
  !  Work variables:
  !
  !  rank of the master in the above communicator
  integer masterid
  !  communicator
  integer comm
  !  exit code for mpi calls
  integer mpinfo
  !  variables for old sphere size, old grid spacing, and
  !  old Hamiltonian size from a previous run
  real(dp) :: rmaxold,hold(3),avec_old(3,3)
  integer ndimold,nwold
  !  the following are actually used only when extrapolating from a
  !  grid pertaining to a smaller sphere
  !  arrays containing 1d->3d grid mapping for the stored data
  integer, dimension (:), allocatable :: kxold, kyold, kzold
  !  array containing the mapping from the 1d vector of the STORED
  !  grid to the 1d vector of the NEW grid
  integer, allocatable :: map(:)
  !  array for temporary storage of vectors that need to be mapped
  !  from the old 1d ordering to the new 1d ordering
  real(dp), allocatable :: fnold(:)
  complex(dpc), allocatable :: zfnold(:)
  !  counters
  integer i, ii, jj, isp, igrp
  real(dp) :: rtmp
  !  array for temporary storage of eigenvectors
  complex(dpc), allocatable :: zwf(:)
  real(dp), allocatable :: wf(:)
  !  temporary integer for ibuff
  integer ibuff,ibuff_send
  !  actual number of spins used
  integer nspin
  !  sphere size
  real(dp) :: rmax
  !  current Hamiltonian size, wedge size
  integer ndim,nwedge
  !  allocation check
  integer alcstat
  !  date label on parsec.dat file
  character (len=26) :: datelabel
  !  counters for representations
  integer nrep,nrold,irp
  !  temporary arrays for eigenvalues/occupations, before reshaping
  !  according to symmetry operations
  real(dp), dimension(:), allocatable :: en_tmp
  !  jrep keeps track of how many eigenstates are already in each
  !  representation
  integer :: jrep(elec_st%nrep)
  !  old character table
  integer :: chiold(elec_st%nrep,elec_st%nrep)
  !
  !  flags for type of system
  !  stype(1) = number of spins
  !  stype(2) = 0 : real wavefunctions
  !  stype(2) = 1 : complex wavefunctions
  !  stype(3) = 0 : confined system
  !  stype(3) = 1 : system periodic on all directions
  !
  integer stype(10)
  ! tolerance threshold
  real(dp), parameter :: tol = 1.d-10
  ! selected eigenvectors found in parsec.dat
  integer :: nsave
  integer, allocatable :: indxsave(:)
#ifdef MPI
  integer status(MPI_STATUS_SIZE)
#endif

  !  kpoint variables
  integer kpnum, kpold, kplp
  !  coordinates of each kpoint
  real(dp), dimension(:,:), allocatable :: kpts
  real(dp) :: kp_diff(3)
  !  true if there is mismatch between kpoints in parsec.dat and kpoints
  !  in memory
  logical :: kpt_mismatch

  !---------------------------------------------------------------

  kpnum = max(elec_st%nkpt,1)

  nspin = elec_st%nspin
  rmax = grid%rmax
  ndim = grid%ndim
  nwedge = grid%nwedge
  masterid = parallel%masterid
  comm = parallel%comm
  nrep = elec_st%nrep
  ierr = 0

  if (parallel%iammaster) then
     !  Open unit.
     open(iunit,file='parsec.dat',form='unformatted',status='old',iostat=ii)
     if (ii /= 0) then
        write(7,*) 'ERROR: parsec.dat file not found.'
        write(7,*) 'Stop in restart_run.'
        ierr = 371
        goto 14
     endif

     !  Read the old grid spacing, sphere size, Hamiltonian size.
     read(iunit) datelabel
     write(7,*)
     write(7,*) 'Restarting from previous run'
     write(7,*) '----------------------------'
     write(7,*)
     write(7,*) ' parsec.dat file created on   ',datelabel
     write(7,*)
     read(iunit) stype
     if (stype(1) /= nspin) then
        write(7,*)
        write(7,*) 'ERROR: stored number of spins is incorrect, ', &
             stype(1),nspin
        write(7,*) 'STOP in restart_run.'
        ierr = 361
        goto 14
     endif
     if ((elec_st%cplx .and. stype(2) == 0) .or. &
          (.not. elec_st%cplx .and. stype(2) == 1)) then
        write(7,*)
        write(7,*) 'ERROR: complex algebra incorrect, ', &
             stype(2),elec_st%cplx
        write(7,*) 'STOP in restart_run.'
        ierr = 362
        goto 14
     endif
     if (stype(2) == 2) elec_st%is_so = .true.
     select case(pbc%per)
     case(0)
        jj = 0
     case(1)
        jj = 2
     case(2)
        jj = 3
     case(3)
        jj = 1
     end select
     if (stype(3) /= jj) then
        write(7,*)
        write(7,*) 'ERROR: boundary conditions are incorrect, ', &
             stype(3),pbc%is_on
        write(7,*) 'STOP in restart_run.'
        ierr = 363
        goto 14
     endif
     if (pbc%is_on) then
        read(iunit) (hold(ii),ii=1,3),(rmaxold,ii=1,4)
        read(iunit) ((avec_old(jj,ii),jj=1,3),ii=1,3)
        read(iunit)
        read(iunit)
        read(iunit)
        read(iunit) kpold
        read(iunit)
        read(iunit)
        read(iunit)
        allocate(kpts(3,kpold))
        read(iunit) ((kpts(jj,ii),jj=1,3),ii=1,kpold)
        read(iunit)
        ! Sanity tests.
        hold = hold - grid%step
        avec_old = avec_old - pbc%latt_vec
        if (maxval(abs(hold)) > tol) then
           write(7,*)
           write(7,*) 'ERROR: stored grid spacing conflicts with'
           write(7,*) 'current grid spacing'
           write(7,*) hold + grid%step
           write(7,*) grid%step(1:3)
           write(7,*) 'STOP in restart_run.'
           ierr = 364
        goto 14
        else if (maxval(abs(avec_old)) > tol) then
           write(7,*)
           write(7,*) 'ERROR: the size of the supercell conflicts'
           write(7,*) 'with the stored one'
           write(7,*) avec_old + pbc%latt_vec
           write(7,*) pbc%latt_vec
           write(7,*) 'STOP in restart_run.'
           ierr = 365
        goto 14
        endif

     else
        read(iunit) hold(1), rmaxold

        kpold = 1

        if (hold(1) /= grid%step(1)) then
           write(7,*)
           write(7,*) 'ERROR: stored grid spacing conflicts with'
           write(7,*) 'current grid spacing'
           write(7,*) hold(1),grid%step(1)
           write(7,*) 'STOP in restart_run.'
           ierr = 366
        goto 14
        endif

     endif
     read(iunit) ndimold
     read(iunit)
     read(iunit) nwold,nrold

  !
  !  If the trend of the Hamiltonian size is inconsistent with the
  !  trend of the sphere size - complain and quit.
  !
     if (pbc%per == 3) rmax = rmaxold
     if ((rmax == rmaxold) .and. ((ndim /= ndimold) .or. &
          (nwedge /= nwold))) then
        write(7,*) 'ERROR: Despite using the same grid spacing'
        write(7,*) 'and sphere size, the stored Hamiltonian size is'
        write(7,*) 'incompatible with the current one. Something'
        write(7,*) 'is seriously wrong. Good luck!'
        write(7,*) 'STOP in restart_run.'
        write(7,*) rmax,rmaxold,ndim,ndimold,nwedge,nwold
        ierr = 367
        goto 14
     else if ((rmax > rmaxold) .and. ((ndim < ndimold) .or. &
          (nwedge < nwold))) then
        write(7,*) 'ERROR: Despite using the same grid spacing'
        write(7,*) 'and a larger sphere size, the stored Hamiltonian'
        write(7,*) 'size is smaller than the current one. Something'
        write(7,*) 'is seriously wrong. Good luck!'
        write(7,*) 'STOP in restart_run.'
        write(7,*) rmax,rmaxold,ndim,ndimold,nwedge,nwold
        ierr = 368
        goto 14
     endif
     !  Check if number of representations matches.
     if (nrold /= nrep) then
        write(7,*)
        write(7,*) 'ERROR: number of irreducible representations'
        write(7,*) 'does not match. ',nrold,nrep
        write(7,*) 'STOP in restart_run.'
        ierr = 369
        goto 14
     endif

     !  Read previous character table and compares it with the current
     !  one; skip additional information about symmetries.
     do ii = 1, 4
        read(iunit)
     enddo
     read(iunit) ((chiold(ii,jj),ii=1,nrold),jj=1,nrold)
     chiold = chiold - elec_st%chi
     ii = minval(chiold)
     jj = maxval(chiold)
     if (ii /= 0 .or. jj /= 0) then
        write(7,*)
        write(7,*) 'ERROR: character table does not match'
        write(7,*) chiold
        write(7,*) elec_st%chi
        ierr = 370
        goto 14
     endif
  endif ! parallel%iammaster
14 continue
  !  Done with sanity checks, send flag to other PEs.
#ifdef MPI
  call MPI_BCAST(ierr,1,MPI_INTEGER,masterid,comm,mpinfo)
#endif
  if (ierr /= 0) return
  !
  ! If restarting from spin-orbit calculation, update flags and array sizes.
  !
#ifdef MPI
  call MPI_BCAST(elec_st%is_so,1,MPI_LOGICAL,masterid,comm,mpinfo)
#endif
  if (elec_st%is_so) then
     elec_st%cplx = .true.
     elec_st%so_pert = .false.
     elec_st%mxwd = 2
     parallel%mxwd = 2
     nspin = 1
     if (parallel%iammaster) then
        if (associated(parallel%zftmp)) deallocate(parallel%zftmp)
        allocate(parallel%zftmp &
             (grid%nwedge*parallel%mxwd),stat=alcstat)
        call alccheck('parallel%zftmp', &
             grid%nwedge*parallel%mxwd,alcstat)
     endif
     elec_st%nstate = 2 * elec_st%nstate
     jj = (solver%winsize + solver%nadd + elec_st%nstate)*elec_st%nrep
     deallocate(elec_st%irep)
     allocate(elec_st%irep(jj,kpnum,elec_st%nspin))
  endif

     if (elec_st%cplx) then
        allocate(zwf(parallel%mydim*parallel%mxwd),stat=alcstat)
        call alccheck('zwf',parallel%mydim*parallel%mxwd,alcstat)
     else
        allocate(wf(nwedge),stat=alcstat)
        call alccheck('wf',nwedge,alcstat)
     endif
  if (parallel%iammaster) then
     allocate(kxold(nwold),stat=alcstat)
     call alccheck('kxold',nwold,alcstat)
     allocate(kyold(nwold),stat=alcstat)
     call alccheck('kyold',nwold,alcstat)
     allocate(kzold(nwold),stat=alcstat)
     call alccheck('kzold',nwold,alcstat)
     allocate(map(nwold),stat=alcstat)
     call alccheck('map',nwold,alcstat)
     map(:) = 0
     allocate(fnold(nwold),stat=alcstat)
     call alccheck('fnold',nwold,alcstat)
     if (elec_st%cplx) then
        allocate(zfnold(nwold*parallel%mxwd),stat=alcstat)
        call alccheck('zfnold',nwold*parallel%mxwd,alcstat)
     endif

     rho(:,:) = zero
     vold(:,:) = zero

     !  Read the irreducible wedge.
     read(iunit) (kxold(i),kyold(i),kzold(i),i=1,nwold)

     !  Set extrapolation flag to zero.
     ixtrpflg = 0

     if (rmax /= rmaxold) then
        !  If current sphere larger or smaller than stored one, extrapolation
        !  is necessary. Write a warning to parsec.out, and set the
        !  extrapolation flag to one.
        write(7,*)
        write(7,*) 'WARNING: Extrapolating an initial guess based on'
        write(7,*) 'stored data for a sphere of radius ',rmaxold,'.'
        write(7,*)
        ixtrpflg = 1
     endif
     !
     !  Compute the map for converting the order in the stored vectors
     !  to the order needed for the present vectors, by converting from
     !  1d to 3d "old style" and converting back from 3d to 1d using
     !  the current indexg map. Must skip grid points outside the boundary
     !  sphere if old radius is bigger than current radius.
     !
     do ii = 1, nwold
        hold(1) = ( kxold(ii) + grid%shift(1) )*grid%step(1)
        hold(2) = ( kyold(ii) + grid%shift(2) )*grid%step(2)
        hold(3) = ( kzold(ii) + grid%shift(3) )*grid%step(3)
        hold(1:pbc%per) = zero
        rtmp = sqrt(dot_product(hold,hold))
        if (rtmp > rmax) cycle
        jj = grid%indexw(kxold(ii),kyold(ii),kzold(ii))
        if (jj > 0 .and. jj <= nwedge) map(ii) = jj
     enddo
  endif !  parallel%iammaster

#ifdef MPI
  call MPI_BCAST(ixtrpflg,1,MPI_INTEGER,parallel%masterid &
       ,parallel%comm,mpinfo)
  call MPI_BCAST(kpold,1,MPI_INTEGER,parallel%masterid &
       ,parallel%comm,mpinfo)
#endif
  !
  !  Start the main loop over spin, where eigenvalues, eigenvectors,
  !  representations and occupancy factors are read from each state.
  !  If eigenstate structures are not defined, skip eigenvalues/vectors
  !  and keep only potentials and charge density.
  !
  do isp = 1, nspin
     jj = isp - 1 + nspin
     do kplp = 1, kpold
        if (parallel%iammaster) read(iunit) ibuff
#ifdef MPI
        call MPI_BCAST(ibuff,1,MPI_INTEGER,masterid,comm,mpinfo)
#endif
        !
        !  If this k-point is not in the current set, skip eigenvalues etc.
        !
        kpt_mismatch = .false.
        if (parallel%iammaster .and. elec_st%nkpt > 0 ) then
           kp_diff = kpts(:,kplp) - elec_st%kpts(:,kplp)
           if (maxval(abs(kp_diff)) > tol) then
              kpt_mismatch = .true.
              read(iunit)
              read(iunit)
              read(iunit)
           else
              kpt_mismatch = .false.
           endif
        endif
#ifdef MPI
        call MPI_BCAST(kpt_mismatch,1,MPI_LOGICAL,masterid,comm,mpinfo)
#endif
        if (kpt_mismatch) then
           if (parallel%iammaster) then
              do ii = 1, 3
                 read(iunit)
              enddo
           endif
           cycle
        endif

        ibuff_send = min(ibuff,elec_st%nstate)
        elec_st%irep(:,kplp,isp) = 0
        if (parallel%iammaster) &
             read(iunit) (elec_st%irep(i,kplp,isp),i=1,ibuff_send)
#ifdef MPI
        call MPI_BCAST(elec_st%irep(1,kplp,isp),ibuff_send,MPI_INTEGER &
             ,masterid,comm,mpinfo)
#endif
        !
        !  Update the number of computed states per representation.
        !
        elec_st%eig(:,kplp,isp)%nec = 0
        do i = 1, ibuff_send
           irp = elec_st%irep(i,kplp,isp)
           elec_st%eig(irp,kplp,isp)%nec = elec_st%eig(irp,kplp,isp)%nec + 1
        enddo
        do irp = 1, nrep
           jrep(irp) = elec_st%eig(irp,kplp,isp)%nec
        enddo
        do irp = 1, elec_st%nrep
           if (elec_st%eig(irp,kplp,isp)%nec > elec_st%eig(irp,kplp,isp)%nn) &
                elec_st%eig(irp,kplp,isp)%nn = elec_st%eig(irp,kplp,isp)%nec
           elec_st%eig(irp,kplp,isp)%mm = elec_st%eig(irp,kplp,isp)%nec + &
                2*solver%nadd + solver%winsize
           !
           !  If this representation/spin is assigned to the group of this
           !  processor, update size of wave-function array.
           !
           if ( elec_st%eig(irp,kplp,isp)%group == parallel%mygroup ) then
              call eig_adj_size(elec_st%eig(irp,kplp,isp),parallel%ldn, &
                   elec_st%eig(irp,kplp,isp)%nec,elec_st%cplx)
           else
              call eig_e_adj_size(elec_st%eig(irp,kplp,isp), &
                   elec_st%eig(irp,kplp,isp)%nec)
           endif
        enddo
        elec_st%ntotal(isp) = sum(elec_st%eig(:,kplp,isp)%nec)
        !
        !  Arrays that store eigenvalues and occupancy factors should be
        !  split among the various representations; use jrep as counter.
        !
        if (parallel%iammaster) then
           allocate(en_tmp(ibuff_send))
           jrep = 0
           read(iunit) (en_tmp(i),i=1,ibuff_send)
           read(iunit)
           do i = 1, ibuff_send
              irp = elec_st%irep(i,kplp,isp)
              jrep(irp) = jrep(irp) + 1
              elec_st%eig(irp,kplp,isp)%en(jrep(irp)) = en_tmp(i)
           enddo
           deallocate(en_tmp)
        endif
#ifdef MPI
        do irp = 1, elec_st%nrep
           call MPI_BCAST(elec_st%eig(irp,kplp,isp)%en,elec_st%eig(irp,kplp,isp)%nec &
                ,MPI_DOUBLE_PRECISION,masterid,comm,mpinfo)
        enddo
#endif
     enddo
     !
     !  For all the vectors that follow: First, place them in the
     !  temporary array fnold. Second, initialize the vector to zero,
     !  for the CURRENT Hamiltonian size (i.e., loop to ndim and not to
     !  ndimold). Third, map the old vector (as stored in fnold)
     !  to the needed new order.
     !
     !  NOTE: The vectors that follow are only the charge density and
     !  the wave functions. When extrapolating, there is no need to
     !  extrapolate the old potential, because it needs to be
     !  recalculated anyway (V_hartree is not zero outside the
     !  old sphere).
     !
     !  Store old SCF potential.
     if (parallel%iammaster) then
        read(iunit) (fnold(i),i=1,nwold)
        parallel%ftmp(:) = zero
        do i = 1, nwold
           if (map(i) /= 0) parallel%ftmp(map(i)) = fnold(i)
        enddo
     endif
     call export_function(parallel,vold(1,isp))
     !  Store charge density.
     if (parallel%iammaster) then
        read(iunit) (fnold(i),i=1,nwold)
        parallel%ftmp(:) = zero
        do i = 1, nwold
           if (map(i) /= 0) parallel%ftmp(map(i)) = fnold(i)
        enddo
     endif
     call export_function(parallel,rho(1,jj))
  enddo
  !
  !  If spin-orbit potential is present, then the potential and density
  !  for the second spin component must be read separately. The loop above
  !  goes only over the first spin component.
  !
  if (elec_st%is_so) then
     if (parallel%iammaster) then
        read(iunit) (fnold(i),i=1,nwold)
        parallel%ftmp(:) = zero
        do i = 1, nwold
           if (map(i) /= 0) parallel%ftmp(map(i)) = fnold(i)
        enddo
     endif
     call export_function(parallel,vold(1,2))
     if (parallel%iammaster) then
        read(iunit) (fnold(i),i=1,nwold)
        parallel%ftmp(:) = zero
        do i = 1, nwold
           if (map(i) /= 0) parallel%ftmp(map(i)) = fnold(i)
        enddo
     endif
     call export_function(parallel,rho(1,3))
  endif

  elec_st%eig(:,:,:)%nec = 0
  elec_st%ntotal(:) = 0
  !
  !  If the computation is spin-polarized, construct the total
  !  charge density as the sum of the stored up and down charge
  !  densities.
  !
  if (elec_st%nspin == 2) then
     do i = 1, parallel%mydim
        rho(i,1) = rho(i,2) + rho(i,3)
     enddo
  endif

  if (parallel%iammaster) then
     deallocate(kxold, kyold, kzold)
     deallocate(map)
     deallocate(fnold)
     if (elec_st%cplx) then
        deallocate(zfnold)
        deallocate(zwf)
     else
        deallocate(wf)
     endif
  endif

#ifdef MPI
  call MPI_BARRIER(comm,mpinfo)
#endif
end subroutine restart_run
!===============================================================
