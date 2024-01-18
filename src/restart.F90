!===============================================================
!
! Copyright (C) 2005 Finite Difference Research Group
! This file is part of parsec, http://www.ices.utexas.edu/parsec/
!
! This subroutine reads grid, potential, wave function, energy
! levels, etc, from a previous calculation, as stored in wfn.dat,
! for restart purposes.
! Extrapolation (i.e., reshaping of charge density, and wave
! functions) is always done, although it is only relevant if
! periodic boundary conditions are not used and the previous
! wfn.dat refer to a calculation with smaller enclosing sphere.
!
!---------------------------------------------------------------
subroutine restart(elec_st,grid,pbc,parallel,solver, &
     vold,rho,ixtrpflg,olddat,ierr)

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
  ! include mpi definitions
  !
  ! Input/Output variables:
  !
  ! electronic structure
  type (electronic_struct), intent(inout) :: elec_st
  ! grid related data
  type (grid_data), intent(in) :: grid
  ! periodic boundary conditions data
  type (pbc_data), intent(in) :: pbc
  ! parallel computation related data
  type (parallel_data), intent(inout) :: parallel
  ! solver related data
  type (eigen_solver), intent(inout) :: solver

  ! distributed self-consistent potential and electron density
  ! (passed outside the structure to overcome a bug with the IBM compiler)
  real(dp), intent(out) :: vold(parallel%mydim,elec_st%nspin)
  real(dp), intent(out) :: rho(parallel%mydim,2*elec_st%nspin-1)

  ! flag indicating whether or not extrapolation to a larger sphere
  ! size is performed. It is used in the main code for deciding
  ! whether or not the Hartree and exchange-correlation potentials
  ! need to be calculated upon a restart.
  integer, intent(out) :: ixtrpflg
  ! compatibility flag: read output to wfn.dat if false; read
  ! output to old.dat if true
  logical, intent(in) :: olddat
  ! error flag, 320 < ierr < 341
  integer, intent(out) :: ierr
  !
  ! Work variables:
  !
  ! rank of the master in the above communicator
  integer masterid
  ! communicator
  integer comm
  ! exit code for mpi calls
  integer mpinfo
  ! read in variables for old sphere size, old grid spacing, and
  ! old Hamiltonian size from a previous run
  real(dp) :: rmaxold,hold
  real(dp) :: alold_x,alold_y,alold_z,hold_x,hold_y,hold_z
  integer ndimold,nwold
  ! the following are actually used only when extrapolating from a
  ! grid pertaining to a smaller sphere
  ! arrays containing 1d->3d grid mapping for the stored data
  integer, dimension (:), allocatable :: kxold, kyold, kzold
  ! array containing the mapping from the 1d vector of the STORED
  ! grid to the 1d vector of the NEW grid
  integer, allocatable :: map(:)
  ! array for temporary storage of vectors that need to be mapped
  ! from the old 1d ordering to the new 1d ordering
  real(dp), allocatable :: fnold(:)
  complex(dpc), allocatable :: zfnold(:)
  ! counters
  integer i, jj, kk, isp
  ! array for temporary storage of eigenvectors
  complex(dpc), allocatable :: zwf(:)
  real(dp), allocatable :: wf(:)
  ! temporary integer for ibuff
  integer ibuff,ibuff_send
  integer lstart,node,msgtype,ldim
  ! actual number of spins used
  integer nspin
  ! sphere size
  real(dp) :: rmax
  ! current Hamiltonian size, wedge size
  integer ndim,nwedge
  ! allocation check
  integer alcstat
  ! date label on wfn.dat file
  character (len=26) :: datelabel
  ! counters for representations
  integer nrep,nrold,irp
  ! temporary arrays for eigenvalues/occupations, before reshaping
  ! according to symmetry operations
  real(dp), dimension(:), allocatable :: en_tmp,occ_tmp
  ! jrep keeps track of how many eigenstates are already in each
  ! representation
  integer :: jrep(elec_st%nrep)
  ! old character table
  integer :: chiold(elec_st%nrep,elec_st%nrep)
  integer :: kpnum

  integer :: igrp, gnode

  !---------------------------------------------------------------


  ! Assume no k-points since wfn.dat has no Brillouin zone.
  kpnum = 1

  nspin = elec_st%nspin
  rmax = grid%rmax
  ndim = grid%ndim
  nwedge = grid%nwedge
  masterid = parallel%masterid
  comm = parallel%comm
  nrep = elec_st%nrep

  allocate(kxold(nwedge),stat=alcstat)
  call alccheck('kxold',nwedge,alcstat)
  allocate(kyold(nwedge),stat=alcstat)
  call alccheck('kyold',nwedge,alcstat)
  allocate(kzold(nwedge),stat=alcstat)
  call alccheck('kzold',nwedge,alcstat)
  allocate(map(nwedge),stat=alcstat)
  call alccheck('map',nwedge,alcstat)
  allocate(fnold(nwedge),stat=alcstat)
  call alccheck('fnold',nwedge,alcstat)
  if (elec_st%cplx) then
     allocate(zfnold(nwedge),stat=alcstat)
     call alccheck('zfnold',nwedge,alcstat)
     allocate(zwf(nwedge),stat=alcstat)
     call alccheck('zwf',nwedge,alcstat)
  else
     allocate(wf(nwedge),stat=alcstat)
     call alccheck('wf',nwedge,alcstat)
  endif

  rho(:,:) = zero
  vold(:,:) = zero

  ! initialize the msgtype counter
  msgtype = 0    

  ! open wfn.dat - this is a binary file for reading/writing
  ! final results on eigenvalues, occupations, charge density,
  ! potentials etc. (see detailed list in the output variables list)
  if (olddat) then
     open(4,file='old.dat',form='unformatted',status='unknown')
  else
     open(4,file='wfn.dat',form='unformatted',status='unknown')
  endif

  ! read the old grid spacing, sphere size, Hamiltonian size
  if (olddat) then
     if (pbc%is_on) then
        read(4) hold_x,hold_y,hold_z,alold_x,alold_y &
             ,alold_z,ndimold !, nspin, datelabel
     else
        read(4) hold, rmaxold, ndimold !, nspin, datelabel
     endif
     nwold = ndimold
     nrold = nrep
  else
     read(4) datelabel
     write(7,*)
     write(7,*) 'Restarting from previous run'
     write(7,*) '----------------------------'
     write(7,*)
     write(7,*) ' wfn.dat file created on   ',datelabel
     write(7,*)
     if (pbc%is_on) then
        read(4) hold_x,hold_y,hold_z,alold_x,alold_y &
             ,alold_z,ndimold,nwold,nrold !, nspin
     else
        read(4) hold, rmaxold, ndimold,nwold,nrold !, nspin
     endif
  endif

  ! if the grid spacing is incompatible with the present one,
  ! complain and quit - no interpolation scheme was implemented!

  if (pbc%is_on) then
     if ((hold_x /= grid%step(1)) .or. (hold_y /= grid%step(2 &
          )) .or. (hold_z /= grid%step(3))) then
        write(7,*)
        write(7,*) 'ERROR: stored grid spacing conflicts with'
        write(7,*) 'current grid spacing'
        write(7,*) hold_x,hold_y,hold_z
        write(7,*) grid%step(1:3)
        write(7,*) 'STOP in restart'
        ierr = 321
     else if ((alold_x /= pbc%box_size(1)) .or. &
          (alold_y /= pbc%box_size(2)) .or. &
          (alold_z /= pbc%box_size(3))) then
        write(7,*)
        write(7,*) 'ERROR: the size of the supercell conflicts'
        write(7,*) 'with the stored one'
        write(7,*) alold_x,alold_y,alold_z
        write(7,*) pbc%box_size(1:3)
        write(7,*) 'STOP in restart'
        ierr = 322
     endif

  else
     if (hold /= grid%step(1)) then
        write(7,*)
        write(7,*) 'ERROR: stored grid spacing conflicts with'
        write(7,*) 'current grid spacing'
        write(7,*) hold,grid%step(1)
        write(7,*) 'STOP in restart'
        ierr = 323
        ! If the stored sphere size is larger than the present one,
        ! complain and quit - no "sphere squeezing" is implemented.
     else if (rmaxold > rmax) then
        write(7,*) 'ERROR: stored sphere size larger than current one.'
        write(7,*) 'This case is not implemented.'
        write(7,*) rmaxold,rmax
        write(7,*) 'STOP in restart'
        ierr = 324
     endif
  endif
  !
  ! If the trend of the Hamiltonian size is inconsistent with the
  ! trend of the sphere size - complain and quit.
  !
  if (pbc%is_on) then
     if ((ndim /= ndimold) .or. (nwold /= nwedge)) then
        write(7,*) 'ERROR: Despite using the same grid spacing'
        write(7,*) 'and cell size, the stored Hamiltonian size is'
        write(7,*) 'incompatible with the current /= Something'
        write(7,*) 'is seriously wrong. Good luck!'
        write(7,*) 'STOP in restart'
        write(7,*) ndim,ndimold,nwedge,nwold
        ierr = 325
     endif
  else
     if ((rmax == rmaxold) .and. ((ndim /= ndimold) .or. &
          (nwedge /= nwold))) then
        write(7,*) 'ERROR: Despite using the same grid spacing'
        write(7,*) 'and sphere size, the stored Hamiltonian size is'
        write(7,*) 'incompatible with the current /= Something'
        write(7,*) 'is seriously wrong. Good luck!'
        write(7,*) 'STOP in restart'
        write(7,*) rmax,rmaxold,ndim,ndimold,nwedge,nwold
        ierr = 326
     else if ((rmax > rmaxold) .and. ((ndim < ndimold) .or. &
          (nwedge < nwold))) then
        write(7,*) 'ERROR: Despite using the same grid spacing'
        write(7,*) 'and a larger sphere size, the stored Hamiltonian'
        write(7,*) 'size is smaller than the current /= Something'
        write(7,*) 'is seriously wrong. Good luck!'
        write(7,*) 'STOP in restart'
        write(7,*) rmax,rmaxold,ndim,ndimold,nwedge,nwold
        ierr = 327
     endif
  endif
  ! check if number of representations matches
  if (nrold /= nrep) then
     write(7,*)
     write(7,*) 'ERROR: number of irreducible representations'
     write(7,*) 'does not match. ',nrold,nrep
     write(7,*) 'STOP in restart'
     ierr = 328
  endif

  ! read previous character table and compares it with the current
  ! one, skip additional information about symmetries
  if (olddat) then
     chiold = elec_st%chi
  else
     do kk = 1, 5
        read(4)
     enddo
     read(4) ((chiold(kk,jj),kk=1,nrold),jj=1,nrold)
  endif
  chiold = chiold - elec_st%chi
  kk = minval(chiold)
  jj = maxval(chiold)
  if (kk /= 0 .or. jj /= 0) then
     write(7,*)
     write(7,*) 'ERROR: character table does not match'
     write(7,*) chiold
     write(7,*) elec_st%chi
     ierr = 329
  endif
  ! done with sanity checks, send flag to other PEs
#ifdef MPI
  call MPI_BCAST(ierr,1,MPI_INTEGER,masterid,comm,mpinfo)
#endif
  if (ierr /= 0) return

  ! read the irreducible wedge
  read(4) (kxold(i),kyold(i),kzold(i),i=1,nwold)

  ! set extrapolation flag to zero
  ixtrpflg = 0

  if (rmax > rmaxold) then
     ! If current sphere larger than stored /= extrapolation is
     ! necessary. Write a warning to parsec.out, and set the
     ! extrapolation flag to one.
     write(7,*)
     write(7,*) 'WARNING: Extrapolating an initial guess based on'
     write(7,*) 'stored data for a smaller sphere.'
     write(7,*)
     ixtrpflg = 1
  endif
  !
  ! Compute the map for converting the order in the stored vectors
  ! to the order needed for the present vectors, by converting from
  ! 1d to 3d "old style" and converting back from 3d to 1d using
  ! the current indexg map.
  !
  do i=1,nwold
     map(i)=grid%indexw(kxold(i),kyold(i),kzold(i))
  enddo
  !
  ! Start the main loop over spin, where eigenvalues, eigenvectors,
  ! representations and occupancy factors are read from each state.
  ! If eigenstate structures are not defined, skip eigenvalues/vectors
  ! and keep only potentials and charge density.
  !
  do isp=1,nspin
     jj = isp-1+nspin
     read(4) ibuff
     ibuff_send = min(ibuff,elec_st%nstate)
#ifdef MPI
     call MPI_BCAST(ibuff_send,1,MPI_INTEGER,masterid,comm,mpinfo)
#endif
     if (olddat) then
        elec_st%irep(:,kpnum,isp) = 0
        elec_st%irep(1:ibuff_send,kpnum,isp) = 1
     else
        elec_st%irep(:,kpnum,isp) = 0
        read(4) (elec_st%irep(i,kpnum,isp),i=1,ibuff_send)
     endif
     !
     ! Update the number of computed states per representation.
     !
     elec_st%eig(:,kpnum,isp)%nec = 0
     do i = 1, ibuff_send
        irp = elec_st%irep(i,kpnum,isp)
        elec_st%eig(irp,kpnum,isp)%nec= &
             elec_st%eig(irp,kpnum,isp)%nec+1
     enddo
#ifdef MPI
     call MPI_BCAST(elec_st%irep(1,kpnum,isp),ibuff_send,MPI_INTEGER &
          ,masterid,comm,mpinfo)
     do irp = 1, nrep
        jrep(irp) = elec_st%eig(irp,kpnum,isp)%nec
     enddo
     call MPI_BCAST(jrep,nrep,MPI_INTEGER &
          ,masterid,comm,mpinfo)
     do irp = 1, nrep
        elec_st%eig(irp,kpnum,isp)%nec = jrep(irp)
     enddo
#endif
     do irp = 1, elec_st%nrep
        if (elec_st%eig(irp,kpnum,isp)%nec >  &
             elec_st%eig(irp,kpnum,isp)%nn) &
             elec_st%eig(irp,kpnum,isp)%nn =  &
             elec_st%eig(irp,kpnum,isp)%nec
        elec_st%eig(irp,kpnum,isp)%mm =  &
             elec_st%eig(irp,kpnum,isp)%nec + &
             2*solver%nadd + solver%winsize
        !
        ! If this representation/spin is assigned to the group of this
        ! processor, update size of wave-function array.
        !
        if ( elec_st%eig(irp,kpnum,isp)%group == &
             parallel%mygroup ) then
           call eig_adj_size(elec_st%eig(irp,kpnum,isp),parallel%ldn, &
                elec_st%eig(irp,kpnum,isp)%nec,elec_st%cplx)
        else
           call eig_e_adj_size(elec_st%eig(irp,kpnum,isp), &
                elec_st%eig(irp,kpnum,isp)%nec)
        endif
     enddo
     elec_st%ntotal(isp) = sum(elec_st%eig(:,kpnum,isp)%nec)
     !
     ! Arrays that store eigenvalues and occupancy factors should be
     ! split among the various representations; use jrep as counter.
     !
     allocate(en_tmp(ibuff_send))
     allocate(occ_tmp(ibuff_send))
     jrep = 0
     read(4) (en_tmp(i),i=1,ibuff_send)
     read(4) (occ_tmp(i),i=1,ibuff_send)
     do i = 1, ibuff_send
        irp = elec_st%irep(i,kpnum,isp)
        jrep(irp) = jrep(irp) + 1
        elec_st%eig(irp,kpnum,isp)%en(jrep(irp)) = en_tmp(i)
        elec_st%eig(irp,kpnum,isp)%occ(jrep(irp)) = occ_tmp(i)
     enddo
     elec_st%ifmax(isp) = 0
     do i = ibuff_send, 1, -1
        if (occ_tmp(i) /= zero) then
           elec_st%ifmax(isp) = i
           exit
        endif
     enddo
#ifdef MPI
     do irp = 1, elec_st%nrep
        call MPI_BCAST(elec_st%eig(irp,kpnum,isp)%en &
             ,elec_st%eig(irp,kpnum,isp)%nec &
             ,MPI_DOUBLE_PRECISION,masterid,comm,mpinfo)
     enddo
#endif
     !
     ! Update solver.
     !
     if (solver%do_subsp) then
        do irp = 1, elec_st%nrep
           do kk = 1, elec_st%eig(irp,kpnum,isp)%nec
              solver%eval_loc(kk,irp,isp) = &
                   elec_st%eig(irp,kpnum,isp)%en(kk)
           enddo
           solver%eig_init(irp,kpnum,isp) = .true.
        enddo
     endif

     deallocate(en_tmp)
     deallocate(occ_tmp)
     ! store old SCF potential
     read(4) (fnold(i),i=1,nwold)
     parallel%ftmp(:) = zero
     do i=1,nwold
        parallel%ftmp(map(i)) = fnold(i)
     enddo
     call export_function(parallel,vold(1,isp))
     !
     ! For all the vectors that follow: First, place them in the
     ! temporary array fnold. Second, initialize the vector to zero,
     ! for the CURRENT Hamiltonian size (i.e., loop to ndim and not to
     ! ndimold). Third, map the old vector (as stored in fnold)
     ! to the needed new order.
     !
     ! NOTE: The vectors that follow are only the charge density and
     ! the wave functions. When extrapolating, there is no need to
     ! extrapolate the old potential, because it needs to be
     ! recalculated anyway, because V_hartree is not zero outside the
     ! old sphere.
     !

     ! store charge density
     read(4) (fnold(i),i=1,nwold)
     parallel%ftmp(:) = zero
     do i=1,nwold
        parallel%ftmp(map(i)) = fnold(i)
     enddo
     call export_function(parallel,rho(1,jj))
     !
     ! Store wave functions (again, split according to symmetry
     ! representation) and distribute them among processors.
     ! jrep is reinitialized and used to reshape wfn arrays.
     !
     jrep = 0
     do kk=1,ibuff
        if (elec_st%cplx) then
           read(4) (zfnold(i),i=1,nwold)
        else
           read(4) (fnold(i),i=1,nwold)
        endif
        if (kk <= ibuff_send) then
           irp = elec_st%irep(kk,kpnum,isp)
           jrep(irp) = jrep(irp) + 1
           if (elec_st%cplx) then
              zwf(:) = zzero
              do i=1,nwold
                 zwf(map(i)) = zfnold(i)
              enddo
           else
              wf(:) = zero
              do i=1,nwold
                 wf(map(i)) = fnold(i)
              enddo
           endif
           igrp = elec_st%eig(irp,kpnum,isp)%group
           do gnode = 0, parallel%group_size-1
              node = parallel%gmap(gnode + 1, igrp + 1)
              msgtype = msgtype + 1
              lstart = parallel%irows(gnode)
              ldim = parallel%irows(gnode+1) - lstart
              if (elec_st%cplx) then
                 if (parallel%iam == node) then
                    call zcopy(ldim,zwf(lstart),1, &
                         elec_st%eig(irp,kpnum,isp)%zwf(1,jrep(irp)),1)
                 else
#ifdef MPI
                    call MPI_SEND(zwf(lstart), ldim, &
                         MPI_DOUBLE_COMPLEX, node,msgtype, comm, mpinfo)
#endif
                 endif
              else
                 if (parallel%iam == node) then
                    call dcopy(ldim,wf(lstart),1, &
                         elec_st%eig(irp,kpnum,isp)%wf(1,jrep(irp)),1)
                 else
#ifdef MPI
                    call MPI_SEND(wf(lstart), ldim, &
                         MPI_DOUBLE_PRECISION, node,msgtype, comm, mpinfo)
#endif
                 endif
              endif
           enddo
           write(9,*) 'unit 4:wf for level= ',kk, &
                ' has been spread'
        endif
     enddo
  enddo                     ! do isp=1,nspin
  !
  ! If the computation is spin-polarized, construct the total
  ! charge density as the sum of the stored up and down charge
  ! densities.
  !
  if (nspin == 2) then
     do i=1,parallel%mydim
        rho(i,1) = rho(i,2) + rho(i,3)
     enddo
  endif

  close(4)

  deallocate(kxold, kyold, kzold)
  deallocate(map)
  deallocate(fnold)
  if (elec_st%cplx) then
     deallocate(zfnold)
     deallocate(zwf)
  else
     deallocate(wf)
  endif

#ifdef MPI
  call MPI_BARRIER(comm,mpinfo)
#endif
end subroutine restart
!===============================================================
!
! Processors without access to i/o receive restart information
! from master PE.
!
!---------------------------------------------------------------
subroutine restart_mpi(elec_st,parallel,solver,vold,rho,ierr)

  use constants
  use electronic_struct_module
  use parallel_data_module
  use eigen_solver_module
#ifdef MPI
  use mpi
#endif
  implicit none
  !
  ! Input/Output variables:
  !
  ! electronic structure
  type (electronic_struct), intent(inout) :: elec_st
  ! parallel computation related data
  type (parallel_data), intent(inout) :: parallel
  ! solver related data
  type (eigen_solver), intent(inout) :: solver

  ! distributed self-consistent potential and electron density
  ! (passed outside the structure to overcome a bug with the IBM compiler)
  real(dp), intent(out) :: vold(parallel%mydim,elec_st%nspin)
  real(dp), intent(out) :: rho(parallel%mydim,2*elec_st%nspin-1)
  ! error flag
  integer, intent(out) :: ierr
  !
  ! Work variables:
  !
  integer isp,kk,ibuff,msgtype,mpinfo,nrep,irp
  integer jrep(elec_st%nrep)
#ifdef MPI
  integer status(MPI_STATUS_SIZE)
#endif
  integer :: kpnum
  !---------------------------------------------------------------

  ! Assume no k-points since wfn.dat has no Brillouin zone.
  kpnum = 1

  msgtype = parallel%iam + 1
  nrep = elec_st%nrep
  ierr = 0

  rho(:,:) = zero
  vold(:,:) = zero

#ifdef MPI
  ! Bail out if wfn.dat file did not pass sanity checks
  call MPI_BCAST(ierr,1,MPI_INTEGER, &
       parallel%masterid,parallel%comm,mpinfo)
  if (ierr /= 0) return

  do isp = 1,elec_st%nspin

     call MPI_BCAST(ibuff,1,MPI_INTEGER, &
          parallel%masterid,parallel%comm,mpinfo)
     elec_st%irep(:,kpnum,isp) = 0
     call MPI_BCAST(elec_st%irep(1,kpnum,isp),ibuff,MPI_INTEGER, &
          parallel%masterid,parallel%comm,mpinfo)
     do irp = 1, nrep
        jrep(irp) = elec_st%eig(irp,kpnum,isp)%nec
     enddo
     call MPI_BCAST(jrep,nrep,MPI_INTEGER, &
          parallel%masterid,parallel%comm,mpinfo)
     do irp = 1, nrep
        elec_st%eig(irp,kpnum,isp)%nec = jrep(irp)
     enddo
     do irp = 1, elec_st%nrep
        if (elec_st%eig(irp,kpnum,isp)%nec > &
             elec_st%eig(irp,kpnum,isp)%nn) &
             elec_st%eig(irp,kpnum,isp)%nn =  &
             elec_st%eig(irp,kpnum,isp)%nec
        elec_st%eig(irp,kpnum,isp)%mm =  &
             elec_st%eig(irp,kpnum,isp)%nec + &
             2*solver%nadd + solver%winsize
        !
        ! If this representation/spin is assigned to the group of this
        ! processor, update size of wave-function array.
        !
        if ( elec_st%eig(irp,kpnum,isp)%group == &
             parallel%mygroup ) then
           call eig_adj_size(elec_st%eig(irp,kpnum,isp),parallel%ldn, &
                elec_st%eig(irp,kpnum,isp)%nec,elec_st%cplx)
        else
           call eig_e_adj_size(elec_st%eig(irp,kpnum,isp), &
                elec_st%eig(irp,kpnum,isp)%nec)
        endif
     enddo
     elec_st%ntotal(isp) = sum(elec_st%eig(:,kpnum,isp)%nec)

#ifdef MPI
     do irp = 1, elec_st%nrep
        call MPI_BCAST(elec_st%eig(irp,kpnum,isp)%en, &
             elec_st%eig(irp,kpnum,isp)%nec,MPI_DOUBLE_PRECISION, &
             parallel%masterid,parallel%comm,mpinfo)
     enddo
#endif
     !
     ! Update solver.
     !
     if (solver%do_subsp) then
        do irp = 1, elec_st%nrep
           do kk = 1, elec_st%eig(irp,kpnum,isp)%nec
              solver%eval_loc(kk,irp,isp) = &
                   elec_st%eig(irp,kpnum,isp)%en(kk)
           enddo
           solver%eig_init(irp,kpnum,isp) = .true.
        enddo
     endif

     call export_function(parallel,vold(1,isp))
     kk = isp-1+elec_st%nspin
     call export_function(parallel,rho(1,kk))

     jrep = 0
     do kk = 1,ibuff
        irp = elec_st%irep(kk,kpnum,isp)
        jrep(irp) = jrep(irp) + 1
        if ( elec_st%eig(irp,kpnum,isp)%group == &
             parallel%mygroup ) then
           if (elec_st%cplx) then
              call MPI_RECV(elec_st%eig(irp,kpnum,isp)%wf(1,jrep(irp)), &
                   parallel%mydim,MPI_DOUBLE_COMPLEX, &
                   parallel%masterid,msgtype,parallel%comm,status, &
                   mpinfo)
           else
              call MPI_RECV(elec_st%eig(irp,kpnum,isp)%wf(1,jrep(irp)), &
                   parallel%mydim,MPI_DOUBLE_PRECISION, &
                   parallel%masterid,msgtype,parallel%comm,status, &
                   mpinfo)
           endif
        endif
        msgtype = msgtype + parallel%group_size
     enddo
  enddo
  call MPI_BARRIER(parallel%comm,mpinfo)
#endif

  if (elec_st%nspin == 2) then
     do kk=1,parallel%mydim
        rho(kk,1) = rho(kk,2) + rho(kk,3)
     enddo
  endif

end subroutine restart_mpi
!===============================================================
