!===============================================================================
!
!                          +-+-+-+-+-+-+-+-+-+-+-+-+
!                            P A R S E C  -  1.4.X
!                          +-+-+-+-+-+-+-+-+-+-+-+-+
!
! PARESC, the real-space, high-order finite difference electronic structure code.
!
! PARSEC solves the Kohn-Sham equations by expressing electron wave-functions
! directly in real space, without the use of explicit basis sets.  It uses
! norm-conserving pseudopotentials (Troullier-Martins and other varieties). It
! is designed for ab-initio quantum-mechanical calculations of the electronic
! structure of matter, within density-functional theory.  PARSEC is optimized
! for massively parallel computing environment, but it is also compatible with
! serial machines. A finite-difference approach is used for the calculation of
! spatial derivatives.  Owing to the sparsity of the Hamiltonian matrix, the
! Kohn-Sham equations are solved by direct diagonalization, with the use of
! extremely efficient sparse-matrix eigensolvers. 
! Some of its features are:
! * Choice of boundary conditions: periodic (any direction), or confined.
! * Structural relaxation.
! * Simulated annealing.
! * Langevin molecular dynamics.
! * Polarizability calculations (confined-system boundary conditions only).
! * Spin-orbit coupling.
! * Non-collinear magnetism.
!
! Copyright (C) 2005 Finite Difference Research Group
! This file is part of parsec, http://www.ices.utexas.edu/parsec/
!
! THIS VERSION IS UNDER DEVELOPMENT It is based on 1.3.6_2D_NCL version branch
! 
! Basic roadmap:
!
! 1. Prepare sets of simple tests, examples and benchmarks suitable For modern
! hardware.
!
! 2. Find more bugs and fix them
!
! 3. The spin-orbit part should be throughly checked and the integration should
! be tighter. Right now the scaling enhancements and openmp stuff is not
! implemented for the SO matvec.
!
! 4. All short names for variables and functions should be changed to more
! informative strings. 
!
! 5. Update manual and documentation, move to dOxygen
!
!
! This program is free software; you can redistribute it and/or modify it
! under the terms of the GNU General Public License as published by the 
! Free Software Foundation; either version 3 of the License, or (at your
! option) any later version.
!
! This program is distributed in the hope that it will be useful, but 
! WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY 
! or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License 
! for more details.
!
! You should have received a copy of the GNU General Public License along 
! with this program.  If not, see <https://www.gnu.org/licenses/>.
!
! The original version of this code was written in 1993-94 by J. R. Chelikowsky,
! N. Troullier, and X. Jing at the University of Minnesota Mathematical
! subroutines were written by Y. Saad and A. Stathopoulos.
!
! Subsequent significant changes and upgrades were made by the following people
!
! (INCOMPLETE list):
!
! H. Kim - core correction, multipole expansion
! I. Vasiliev - Polarizability, asymptotically-corrected
! functionals
! M. Jain and L. Kronik - GGA, generalized Broyden mixer, massive refactoring
! K. Wu - improvements in iterative diagonalizer and Hartree solver
! R. Burdick - major rewrite and improvement of multipole expansion
! M. Alemany and M. Jain - periodic boundary conditions, parallel implementation
! M. Alemany and E. Lorin - ARPACK solver (simplified matvec interface by AJB)
! H. Kwak and M. Tiago - LDA+U
! M. Tiago - Many other things
! Y. Zhou - More Eigensolvers
! A. Natan - PBC and more
! D. Nave - Spin orbit portions
! A. J. Biller -  Matvec rehaul, threading support (openMP), passive aggressive
! code-commenting, etc.
!
! And many others...
!-------------------------------------------------------------------------------
program parsec
! use self_restraint_when_refactoring
  use constants
  use matvec_module
  use cluster_module
  use electronic_struct_module
  use grid_module
  use pseudo_potential_module
  use molecular_dynamic_module
  use pbc_module
  use mixer_module
  use movement_module
  use symmetry_module
  use parsec_global_data
  use bandstruc_module
  use nscf_module
#ifdef MPI
  use mpi
#endif
#ifdef USEHDF5
  use hdf5
#endif
#ifdef OMPFUN
  use omp_lib
#endif

  implicit none
#ifdef ITAC
  ! trace-analyzer definitions
    include 'VT.inc'
  ! sadly one needs global variables for this
    include 'vtcommon.inc'
#endif 
  ! other variable definitions
  include 'def.h'

  ! the cluster
  type (cluster) :: clust
  ! electronic structure
  type (electronic_struct) :: elec_st
  ! bandstructure structure
  type (bandstruc) :: band_st
  ! nscf structure
  type (nscf_data) :: nscf
  ! grid related data
  type (grid_data) :: grid
  ! pseudopotential related data
  type (pseudo_potential) :: p_pot 
  ! molecular dynamic related data
  type (molecular_dynamic) :: mol_dynamic
  ! periodic boundary conditions data
  type (pbc_data) :: pbc
  ! mixer related data
  type (mixer_data) :: mixer
  ! movement data
  type (movement) :: move
  ! symmetry operations:
  ! symm: full Point group; rsymm: Abelian subgroup
  type (symmetry) :: symm, rsymm
  ! hdf5 file 
#ifdef USEHDF5 
  CHARACTER(LEN=9) :: filename = "parsec.h5" ! File name
  INTEGER(HID_T) :: file_id                  ! File identifier
#else
  ! let FORTRAN deal with it
  integer :: file_id
#endif
#ifdef OLDDEBUG
  integer :: kk, pp
#endif
#ifdef ITAC
  integer :: vt_ierr
#endif
  ! to be moved to their own subroutine:
  integer :: kgrid_tmp(3), idx
  real(dp) :: kgrid_shift_tmp(3),tmp
  ! Transitional object (equivalent to a stuffed toy)
  real(dp), allocatable :: vect(:)

  ! Initialize error flag.
  ierr = 0

#ifdef ITAC
  call VTTRACEOFF(vt_ierr)
#endif
  ! ===============================================================
  ! Define parallel environment.
  ! ===============================================================
  ! Create_parallel_data collectively calls MPI_INIT_THREAD,
  ! A modern MPI_INIT that allows us to use MPI_THREAD_FUNNELED if
  ! openMP is used (one unique master thread for MPI)
  ! Thus it also checks if this mode is available in the MPI runtime
  ! and aborts if not.
  !
  ! Next, the global MPI_COMM_WORLD is assigned, though the comments
  ! discuss something about a group, probably relics from some other times.
  ! The master is found, and iam, iammaster are assigned.
  ! If MPI is not used, a "sterile" parallel object is created.
  !
  ! This, along with all the other object instantiation should actually
  ! be done by type-bound procedures
  ! so, something like this: ierr = parallel%init_parallel
  !                          call exit_err etc...
  !
  ! stuff to test?
  ! that both openmp and regular modes work, with and without mpi. 
  call create_parallel_data (parallel)
  


#ifdef ITAC
  ! ===============================================================
  ! Initialize and define Intel Traceanalyzer API hooks
  ! ===============================================================
  ! This is a bookeeping procedure for low level instrumentation 
  ! only that the name does not fit the actual function.
  ! can't think of anything to test here.
  call init_trace_data (ierr)
#endif

  ! ===============================================================
  ! Init stuff - Nullify's all possible pointers so associated()
  !     is false
  ! ===============================================================
  ! Naturally this should be wrapped into a single procedure
  ! and all of this init should be done using type-bound procedures
  ! also, init's and nullification should be differentiated.
  call init_parallel_data(parallel) 
  call init_cluster(clust)  
  call init_electronic_struct(elec_st) !nullifies pointess, AND sets mxwd,ndim,nrep to -1
  call init_grid(grid) 
  call init_pseudo_potential(p_pot) 
  call init_molecular_dynamic(mol_dynamic)
  call init_mixer(mixer) 
  call init_pbc(pbc) !also sets pbc%is_on=.false.
  call init_nonloc_pseudo_pot(nloc_p_pot) 
  call init_movement(move) 
  call init_symmetry(symm)
  call init_symmetry(rsymm)
  call init_matvec() !notice natmi,dwvec,zwvec are to be global, so no data struct argument. Also, why?
  call init_nscf(nscf) !also sets nscf_on=.false. and 5 other assignments
  call init_potential(pot) 

  ! ===============================================================
  ! Start the clock
  ! ===============================================================
  if (parallel%iammaster) then
  call mysecond(tstrt)
  endif
  ! set the setup state on. why? see later.
  is_setup = .true.

  ! ===============================================================
  ! Initialize standard output:
  ! out.* (unit 9) = scratch output from all PEs, including master
  ! ===============================================================
  ! This should be encapsulated within a function, since
  ! eventually we would like this to be done only if a certain verbosity level was stated,
  ! since with thousands of PEs this is just messy.
  ! also, we should haev one of three _mutually exclusive_ output types - 
  ! mpi i/o, hdf5 and std "write".
  !
  ! Set the extension string for my node.
  write(idstring,'(I4.4)') parallel%iam !limited for 10k processes
  open(9,file='out.'//idstring,form='formatted',status='unknown')
  write(9,*) 'Processor No ',parallel%iam,'  has started'
  write(9,*) 'Number of processors: ',parallel%procs_num
#ifdef OMPFUN
  write(9,*) 'Number of threads/processor: ',OMP_GET_MAX_THREADS()
#else
  write(9,*) 'Number of threads/processor: *** No openMP support'
#endif

  ! Initialize error flag.
  ierr = 0
#ifdef USEHDF5 
!    Initialize FORTRAN interface.
     CALL h5open_f (ierr)
     ! Create a new file using default properties.
     CALL h5fcreate_f(filename, H5F_ACC_TRUNC_F, &
       file_id, ierr)
     ! Terminate access to the file.
!     CALL h5fclose_f(file_id, ierr)
!    Close FORTRAN interface.
!     CALL h5close_f(ierr)
#endif

  ! ===============================================================
  ! Read in user defined parameters. Performed by master PE only.
  ! ===============================================================
  ! This is bad. usrinput should only read usrinput, not print stuff on screen, 
  ! and certianly not modify the main datastructures.
  ! This is how it's going to be:
  ! 0. printed out data to parsec.out has to begin before this stage.
  ! 1. read and parse user input into unique data structure
  ! 2. sanity checks for weird parameters and memory 
  ! 3. if master -> then print to screen what's going to happen
  ! 4. send the user input structure to all PEs
  ! 5. only now make the modification to the other data structures based on the input
  
  if (parallel%iammaster) call usrinput(clust,band_st,elec_st,grid, &
       pot,p_pot,mol_dynamic,pbc,mixer,solver,move,istart, &
       mxiter,vconv,vconv_approach,npolflg,field,outflag,ipr, &
       export_griddata_flag,chsym, &
       parallel%groups_num,parallel%procs_num,wstrt,nscf, &
       oldinpformat,outevflag,readwfndat,ignoresym,outputgw,enable_data_out,file_id,ierr)

#ifdef MPI
  ! so why is it so important to do this now?
  call MPI_BCAST(clust%type_num,1,MPI_INTEGER, &
       parallel%masterid,parallel%comm,mpinfo)

  if (.not. parallel%iammaster)  then
     allocate(p_pot%rws(clust%type_num))
     p_pot%rws(:) = zero
  endif
#endif

  ! so why is it so important to do this now?
  ! BTW this is a clear code-smell. 
  ! why does this number need to be updated twice?
  parallel%mxwd = elec_st%mxwd


  ! "oh by the way, let's check if the input file was read properly"
  call exit_err (clust,elec_st,grid,pot,p_pot,nloc_p_pot,u_pot, &
       move, mol_dynamic,pbc,mixer,solver,symm,rsymm,parallel,&
       file_id,ierr)

  ! this should be moved to the init section
  call create_movement (move)

  ! I wonder what this does?
  if (parallel%iammaster .and. ignoresym) solver%nadd = 1

  ! ===============================================================
  ! Read in pseudopotential files and compute associated quantities.
  ! Performed by master PE only.
  ! ===============================================================
  ! pseudo is a mess, since not all pseudopotential formats were created equal
  ! 1. each different pseudo parsing should be encapsulated in a different function
  ! 2.1 support for UPF should be added, obsolete formats should be removed
  ! 2.2 the extra input flags for fhi are cumbersome
  ! 2.3 should prepare support for PAW 
  ! 3. this subroutine calls "fornberg" - a subroutine from 1994 
  ! that does 1st and 2nd order (and more) derivatives, it needs to be checked.
  ! 4. This part might benefit from openmp enhancement
  ! 5. The pseudo- wavefunctions and potentials are on a radial log grid, and are
  ! being transferred to a uniform (regular?) grid, the new entities are not guaranteed 
  ! to have "support" on the grid spacing (h1,h2,h3). This is a major source of error,
  ! and it should be controlled.
  ! 6. For this purpose (and some others, i guess) there is an option to do a cubic_spline
  ! interpolation, but I know nothing of this, also - it is a 1992 code from NR, which needs
  ! a special license to actually reuse. since this is a silly subroutine, we can replace it
  ! with something from GSL instead, which will probably be even safer.

  if (parallel%iammaster) then
      if ( solver%name /= TEST ) then
     call pseudo(clust,p_pot,elec_st%nstate,elval,elec_st%icorr,ierr)
      endif
  endif
  call exit_err (clust,elec_st,grid,pot,p_pot,nloc_p_pot,u_pot, &
       move, mol_dynamic,pbc,mixer,solver,symm,rsymm,parallel,&
       file_id,ierr)

  ! ===============================================================
   ! fourier_f was once here ! moved to inipbc
  ! ===============================================================

  ! ===============================================================
  ! Initialize various variables and files.
  ! Performed by master PE only.
  ! ===============================================================
  ! Finally initing some stuff. Do not be confused, this isn't init_var
  ! it is rather yet another super-function that should be dissolved - 
  !
  ! who is with intent inout?
  !   clust, elec_st,mol_dynamic, move.
  ! who is with intent out?
  !  ifield (polarization), ierr
  !
  !   Deals with restarting molecular dynamics,
  ! also:
  !   centers the cluster around origin (if not restart, and not turned off)
  ! also:
  !   sets xele in elec_st to be (elval - elec_st%ncharge)
  ! also:
  !   deals with polarization stuff and opens polar.dat if "on" 
  ! also:
  !   if atoms move, opens atom.cor and deals with some stuff for that as well
  ! also:
  !   if the minimization is restarted deals with that
  ! also:
  !   if molecular dynamics are done, deals with those files 
  !      sets the required parameters for the MD run
  ! also:
  !   deals with artifical occupations from occup.in
  ! also:
  !   line 231: AMIR - have to check this change with Doron !!!!!
  !   what was done? back_to_cell was deactivated.


  if (parallel%iammaster) call initial(clust,elec_st,pbc, &
       mol_dynamic,move,elval,npolflg,field,ifield, &
       istart,outevflag,grid%domain_shape,ierr)
  call exit_err (clust,elec_st,grid,pot,p_pot,nloc_p_pot,u_pot, &
       move, mol_dynamic,pbc,mixer,solver,symm,rsymm,parallel,&
       file_id,ierr)
! WE ALL SHOULD BE INSIDE INIT_VAR, ~200 lines below!
#ifdef MPI
  call MPI_BCAST(readwfndat,1,MPI_LOGICAL, &
       parallel%masterid,parallel%comm,mpinfo)
  call MPI_BCAST(chsym,1,MPI_LOGICAL, &
       parallel%masterid,parallel%comm,mpinfo)
  call MPI_BCAST(istart,1,MPI_INTEGER, &
       parallel%masterid,parallel%comm,mpinfo)
  call MPI_BCAST(mxiter,1,MPI_INTEGER, &
       parallel%masterid,parallel%comm,mpinfo)
  call MPI_BCAST(outflag,1,MPI_INTEGER, &
       parallel%masterid,parallel%comm,mpinfo)
  call MPI_BCAST(outevflag,1,MPI_INTEGER, &
       parallel%masterid,parallel%comm,mpinfo)
  call MPI_BCAST(npolflg,1,MPI_INTEGER, &
       parallel%masterid,parallel%comm,mpinfo)
  call MPI_BCAST(ifield,1,MPI_INTEGER, &
       parallel%masterid,parallel%comm,mpinfo)
  call MPI_BCAST(vconv,1,MPI_DOUBLE_PRECISION, &
       parallel%masterid,parallel%comm,mpinfo)
  call MPI_BCAST(vconv_approach,1,MPI_DOUBLE_PRECISION, &
       parallel%masterid,parallel%comm,mpinfo)
  call MPI_BCAST(field,1,MPI_DOUBLE_PRECISION, &
       parallel%masterid,parallel%comm,mpinfo)
  call MPI_BCAST(move%is_on,1,MPI_LOGICAL, &
       parallel%masterid,parallel%comm,mpinfo)
  call MPI_BCAST(move%mxmove,1,MPI_INTEGER, &
       parallel%masterid,parallel%comm,mpinfo)
  call MPI_BCAST(mol_dynamic%is_on,1,MPI_LOGICAL, &
       parallel%masterid,parallel%comm,mpinfo)
  call MPI_BCAST(mol_dynamic%step_num,1,MPI_INTEGER, &
       parallel%masterid,parallel%comm,mpinfo)
  call MPI_BCAST(elec_st%mxwd, &
       1,MPI_INTEGER,parallel%masterid,parallel%comm,mpinfo)
  call MPI_BCAST(parallel%mxwd, &
       1,MPI_INTEGER,parallel%masterid,parallel%comm,mpinfo)
  call MPI_BCAST(elec_st%nspin, &
       1,MPI_INTEGER,parallel%masterid,parallel%comm,mpinfo)
  call MPI_BCAST(elec_st%ncl, &
       1,MPI_LOGICAL,parallel%masterid,parallel%comm,mpinfo)
  call MPI_BCAST(p_pot%rws,clust%type_num,MPI_DOUBLE_PRECISION, &
       parallel%masterid,parallel%comm,mpinfo)
 ! if (p_pot%is_so) then
 !  call MPI_BCAST(p_pot%so, p_pot%type_num,MPI_LOGICAL &
 !          ,parallel%masterid,parallel%comm,mpinfo)
 ! endif

#endif

  ! ===============================================================
  ! Initialize symmetry operations. Performed by master PE only.
  ! ===============================================================
  !
  ! If molecular dynamics is being performed, do not impose
  ! symmetry operations. The random motion of atoms will break any
  ! symmetry. Also, do not use symmetries if polarizability is
  ! calculated, since the applied electric field also breaks
  ! symmetry. For manual movement of atoms, ignoring symmetry
  ! operations prevents a "dumb-user" mistake of using symmetry
  ! operations inconsistent with user-provided atomic coordinates.
  ! if a band structure is calculated symmetry is switched off
  !
  ! ierr is used to define how much symmetry is used
  ! ierr=0 -> full symmetry
  ! ierr=2 -> no symmetry
  ! 
  ! this text should be writted from within the function, even if symm is suppresed
  if (parallel%iammaster) then
     write(7,'(/,a)') ' Symmetry Properties:'
     write(7,'(a)') ' --------------------'

     if(ignoresym) then
        ierr = 2
        write(7,*) "Symmetry suppressed by the Ignore_Symmetry flag."
     else
        ierr = 0

        if (elec_st%nkpt /= 0) ierr = 1
        if (.not. all(grid%shift == half)) ierr = 1
        if (pbc%per==2) ierr = 2
        if (p_pot%is_so .or. elec_st%ncl) ierr = 1
        if (band_st%bands_on) ierr = 1
        if (nscf%nscf_on) ierr = 1

        if (mol_dynamic%is_on) ierr = 2
        if (npolflg == 1) ierr = 2
        if (move%is_on .and. move%name == MANUAL) ierr = 2

        if ((.not. pbc%is_on) .and. grid%domain_shape > 0) then
           write(7,*)
           write(7,*) "WARNING: Symmetry is not currently supported with"
           write(7,*) "non-spherical domain cluster calculations."
           write(7,*) "Symmetry will not be used."
           ierr = 2
        endif
     endif
     call symmetries(clust,grid,pbc,symm,rsymm,ipr,ierr)
  endif
  call exit_err (clust,elec_st,grid,pot,p_pot,nloc_p_pot,u_pot, &
       move, mol_dynamic,pbc,mixer,solver,symm,rsymm,parallel,&
       file_id,ierr)

  ! ===============================================================
  ! If periodic system - initialize Fourier space quantities.
  ! performed by master PE only.
  ! ===============================================================
  ! Actually, changes both grid and pbc.
  ! Most importantly, makes sure that the grid spacings are commensurate with the lattice vectors,
  ! as well as any frac. translation in non-symmorphic operations (adjust_step)
  !
  if (parallel%iammaster .and. pbc%is_on) &
       call inipbc(clust,p_pot,grid,pbc,symm,ipr,ierr)
  call exit_err (clust,elec_st,grid,pot,p_pot,nloc_p_pot,u_pot, &
       move, mol_dynamic,pbc,mixer,solver,symm,rsymm,parallel,&
       file_id,ierr)
  !
  ! SO stuff, notice that this is being setup only for master, then 
  ! sent to others, while it should be done for everyone, no?
  if (parallel%iammaster) then
      if (elec_st%ncl) then
          solver%ncl = .true.
          parallel%mxwd = 2
          if(nloc_p_pot%is_so) elec_st%is_so = .true.
      elseif (nloc_p_pot%is_so) then
          if (elec_st%mxwd == 1) then
              elec_st%is_so = .false.
              parallel%mxwd = 1
          elseif(elec_st%mxwd == 2) then
              elec_st%is_so = .true.
              elec_st%mxwd = 2
              parallel%mxwd = 2
          endif
      endif
  endif
  ! sending of parameters from above ^
#ifdef MPI
  call MPI_BCAST(parallel%mxwd, &
       1,MPI_LOGICAL,parallel%masterid,parallel%comm,mpinfo)
  call MPI_BCAST(elec_st%mxwd, &
       1,MPI_LOGICAL,parallel%masterid,parallel%comm,mpinfo)
  call MPI_BCAST(elec_st%is_so, &
       1,MPI_LOGICAL,parallel%masterid,parallel%comm,mpinfo)
  call MPI_BCAST(nloc_p_pot%is_so, &
       1,MPI_LOGICAL,parallel%masterid,parallel%comm,mpinfo)
  call MPI_BCAST(solver%ncl, &
       1,MPI_LOGICAL,parallel%masterid,parallel%comm,mpinfo)
#endif

#ifdef AJB_DEBUG
  if (parallel%iammaster) then
     write(7,*)
     write(7,*) ' BCATSED spin-orbit flags'
     write(7,*)
  endif
     write(9,*) 'i am here bfore init_var !'
#endif
  ! ===============================================================
  ! Send some basic info over to the procs.
  ! ===============================================================
  ! "some" - an understatement if I ever saw one.
  ! Another superfunction that should be dissolved.
  ! Semi-comp. list of things that are done:
  ! 0. Does something like 180 different bcasts, who can tell if all the data is sent?
  ! 1. Rescaling the grid
  ! 2. Moves information from rsymm to elec_st
  ! 3. Creates elec_st on other PEs than master (which created it on usrinput FCOL)
  ! 4. Same for symmetry
  ! 5. Constructing the M-P kpoint grid if needed
  ! 6. Allocating magmom just because
  ! 7. Calling fornberg (see above) to get the derivative coefficients in grid%coe
  ! 8. Calculate laplacian coefficients if needed 
  ! 9. Creating nonloc_pseudo struct
  ! 10.Creating eigenstate structures etc. 

  call init_var(clust,band_st,elec_st,mixer,solver,symm,rsymm,parallel, &
       nloc_p_pot,grid,p_pot,pbc,ipr,export_griddata_flag,nscf,outputgw,enable_data_out,ierr)
  call exit_err (clust,elec_st,grid,pot,p_pot,nloc_p_pot,u_pot, &
       move, mol_dynamic,pbc,mixer,solver,symm,rsymm,parallel,&
       file_id,ierr)
#ifdef AJB_DEBUG
  if (parallel%iammaster) then
     write(7,*)
     write(7,*) ' survived init var'
     write(7,*)
  endif
#endif

#ifdef ITAC
  call VTTRACEON(vt_ierr)  !AJB: I decided that tracing should begin here 
  call VTBEGIN( vt_grid_setup_state , vt_ierr)  !set the VT state to grid_setup
  if (parallel%iammaster) then
     write(7,*)
     write(7,*) ' VT TRACE BEGIN'
     write(7,*)
  endif
#endif
  ! ===============================================================
  ! Setup the real-space grid and associated index arrays.
  ! ===============================================================
  ! but first, let's do some unrelated stuff:
  if (p_pot%is_so .or. elec_st%ncl) then
     isp = 1
  else
     isp = elec_st%nspin
  endif
  kpnum = max(elec_st%nkpt,1)
  ! setup should really be called grid_setup or something 
  ! also does communication stuff. also is a super function, but here it is
  ! probably justified
  !
  ! 1. create groups of processors if needed, see document in SVN about it
  ! 2. bcasts the grid shift and n1,n2,n3 integers to everyone (finally?)
  ! 3. calls grid partition, this is too large to comment upon here
  ! 4. sends the result of grid partition, also calls aux. grid functions
  ! 5. allocates communication arrays
  ! 6. sets up and looks for neighbor data for communication - so slow
  ! 7. more BCASTs, including data that was already sent
  ! 8. defines the important parallel%mydim
  ! 9. calls comm_neigh which finds out what the communication pattern is
  ! 10. creates the new topologically aware communicators for matvec 

  call setup(grid,pbc,rsymm,parallel,isp,kpnum,ierr)


  call exit_err (clust,elec_st,grid,pot,p_pot,nloc_p_pot,u_pot, &
       move, mol_dynamic,pbc,mixer,solver,symm,rsymm,parallel,&
       file_id,ierr)
  ! hello, I am here because I needed to know about "grid" and ldn.
  ! should also be type-bound to the solver.
  call create_eigen_solver(solver,grid,parallel%ldn,elec_st%nrep, &
       elec_st%nspin,elec_st%nstate,elec_st%nkpt,elec_st%cplx, &
       parallel%mxwd)

  ! ===============================================================
  ! Assign representations/spins round-robin accross the blocks.
  ! If spin-orbital potentials are present, both spin components must
  ! be assigned the same block.
  ! ===============================================================
  ! Blocks? These are groups.
  ! This should be encapsulated, also why here of all places?
  ii = 0
  if (nloc_p_pot%is_so .or. elec_st%ncl) then
     do kplp = 1, max(elec_st%nkpt,1)

        elec_st%eig(1,kplp,1)%group = ii
        ii = ii + 1
        if (mod(ii,parallel%groups_num) == 0) ii = 0

! this can probably be removed because of the 'else' clause
!       AMIR - have to check this with DORON.
!        do irp = 1, elec_st%nrep
!           elec_st%eig(irp,kplp,:)%group = ii
!           ii = ii + 1
!           if (mod(ii,parallel%groups_num) == 0) ii = 0
!        enddo

     enddo
  else
     do isp = 1, elec_st%nspin
        do kplp = 1, max(elec_st%nkpt,1)
           do irp = 1, elec_st%nrep
              elec_st%eig(irp,kplp,isp)%group = ii
              ii = ii + 1
              if (mod(ii,parallel%groups_num) == 0) ii = 0
           enddo
        enddo
     enddo
  endif

  ! ===============================================================
  ! Allocate important variables.
  ! ===============================================================
  ! There you go, no more information needed:
  !
  ! Actually allocate memory for charge
  call electronic_struct_set_charge(parallel%mydim,elec_st)

  if (elec_st%do_vdw) then
  ! same for vdw stuff. 
  call electronic_struct_set_charge_vdw(parallel%mydim,elec_st,clust%atom_num)
    allocate (clust%dist_atm (parallel%mydim,clust%atom_num))
  endif

  ! Actually allocate mixer, the mixer type was predetermined 
  if (.not. nscf%nscf_on) call set_mixer (parallel%mydim,elec_st%nspin,mixer)

  ! Actually allocate potential arrays
  call potential_set_ndim (parallel%mydim,elec_st%nspin,pot,elec_st%is_cur)
  
  ! This is THE best place to put this.
  solver%maxmv = max(solver%maxmv,parallel%ldn*elec_st%nstate/elec_st%nrep)
  if (parallel%iammaster) then
     write(7,*)
     write(7,*) ' Maximum number of matrix-vector operations = ',solver%maxmv
     write(7,*)
  endif
  
  ! the "group" master allocates place to put a variable whose size
  ! is the entire grid. This is probably going to be impossible in the future
  if (parallel%iamgmaster) then
     allocate(parallel%ftmp(grid%nwedge),stat=alcstat)
     call alccheck('parallel%ftmp',grid%nwedge,alcstat)
     parallel%ftmp = zero
     if (elec_st%cplx .or. elec_st%scf_so .or. elec_st%ncl) then
        allocate(parallel%zftmp(grid%nwedge*parallel%mxwd),stat=alcstat)
        call alccheck('parallel%zftmp',grid%nwedge*parallel%mxwd,alcstat)
        parallel%zftmp = zzero
     endif
  endif

  ! ===============================================================
  ! Define factorial parameters and check for isolated atoms.
  ! ===============================================================
  ! Because these two should be together forever.
  ! notice this allocation without check:
  allocate(pot%clm(0:solver%lpole,0:solver%lpole))
  ! setclm is a wonder. it sets the clm. what is clm?
  ! These are the coefficients for the multipole expansion. so yeah,
  ! more than factorial parameters.
  ! btw, these only go up to order 9, and are supposedly accurate.
  call setclm (solver%lpole,pot%clm)
  ! isolat: checks for "insulated atoms"
  ! for PBC - if an atom is too far away from the origin (rmax := one)
  ! for CLUSTER - uses "outside_domain" to check if any atom is outside grid.
  ! This needs some distinction for 2D and 1D PBCs!!!! Amir knows what I'm talking about 
  if (parallel%iammaster) call isolat(clust,grid,pbc%per,ierr)
  call exit_err (clust,elec_st,grid,pot,p_pot,nloc_p_pot,u_pot, &
       move, mol_dynamic,pbc,mixer,solver,symm,rsymm,parallel,&
       file_id,ierr)

#ifdef ITAC
  call VTEND( vt_grid_setup_state, vt_ierr)
  call VTBEGIN( vt_nonloc, vt_ierr)
#endif
  ! ===============================================================
  ! Calculate the non-local contribution to the Hamiltonian
  ! Master PE sets up variables and then export to all procs.
  ! ===============================================================
  ! Yes. Worth noting that this means |Vl-Vloc|phi_lm>
  ! The export is dangerous since it uses allreduces hoping that the arrays
  ! are zeroed out beforehand. also, many BCASTs.
  ! Should mention: the "global" non local array that was sent to every process
  ! is then trimmed and restructured by every process so it is local
  if(solver%name /= TEST) then
  if (parallel%iammaster) then
     write(7,*)
     write(7,*) ' Calculating the non-local contribution to the Hamiltonian...'
     write(7,*)
  endif
  call nonloc(clust,grid,p_pot,nloc_p_pot,pbc,parallel,ierr)
  call exit_err (clust,elec_st,grid,pot,p_pot,nloc_p_pot,u_pot, &
       move, mol_dynamic,pbc,mixer,solver,symm,rsymm,parallel,&
       file_id,ierr)
  else
    write(9,*) 'TEST MODE - Skipping V_non-local'
  endif
#ifdef ITAC
  call VTEND( vt_nonloc, vt_ierr)
  call VTBEGIN( vt_grid_setup_state, vt_ierr)
#endif
  ! ===============================================================
  ! Calculate the on-site Coulomb interaction contribution
  ! to the Hamiltonian.
  ! Master PE sets up variables and then export to all procs.
  ! ===============================================================
  ! HERE BE DRAGONS
  if(solver%name /= TEST) then
  call upot(clust,grid,p_pot,u_pot,pbc,parallel,elec_st%nspin, &
       elec_st%cplx,ierr)
  call exit_err (clust,elec_st,grid,pot,p_pot,nloc_p_pot,u_pot, &
       move, mol_dynamic,pbc,mixer,solver,symm,rsymm,parallel,&
       file_id,ierr)
  else
    write(9,*) ' TEST MODE - Skipping on-site Coulomb interaction'
  endif

  ! ===============================================================
  ! Report memory usage and minimum memory needed for the
  ! calculation.
  ! ===============================================================
  ! Yes, this is outdated. Does not think about chebff, and certainly does not include buffers
  if (parallel%iammaster) call chkmemory(elec_st,grid, &
       pbc,mixer,solver,parallel)

  ! ===============================================================
  ! Compute core-correction charge density. If... maxval(p_pot%icore)>0
  ! ===============================================================
  ! sets up superposition of core-correction charge density 
  if(solver%name /= TEST) then
  call corecd(clust,grid,p_pot,pbc,parallel,elec_st%rhoc,ierr)
  call exit_err (clust,elec_st,grid,pot,p_pot,nloc_p_pot,u_pot, &
       move, mol_dynamic,pbc,mixer,solver,symm,rsymm,parallel,&
       file_id,ierr)
  else
      write(9,*) ' TEST MODE - Skipping core-correction'
  endif

#ifdef ITAC
  call VTEND( vt_grid_setup_state, vt_ierr)
  call VTBEGIN( vt_pbc_ion, vt_ierr)
#endif




! ===============================================================
! Set ionic potential - superposition of local ionic potential
! ===============================================================
! 
if (solver%name /= TEST) then
        select case(pbc%per)
        case(0)
                ! sets pot%vion, for grid points away from ions, the ionic
                ! potential Zion/r  is used.  If any gridpoint is within rc of an ion,
                ! the local pseudopotential is interpolated on the spot 
                call ionpot(clust, grid, p_pot, parallel, pot%vion, oldinpformat)
        case(1)
                ! constructs the local potential using a superposition of the
                ! ions and compensating gaussian charges
                ! using the domain decomposition each PE eventually does:
                !   pot%vion(1:mydim)=vion_local(1:mydim)-rho_gauss(1:mydim)
                ! 
                ! rho_gauss is a potential compensating for the fake charge added to vion_local, see below
                ! 
                ! vion_local is a superimposed contribution from each atom, 
                ! although with Gaussian charge dist. cancellation to reduce the need for summing over
                ! the cell replicas. The Zion/r form is used here as in ionpot, 
                ! else, an on-the-spot interpolation is made for the local part of the
                ! pseudopotential
                !
                ! WARNING: The interpolation is taken from ionpot, but only the part that
                ! assumes a logarithmic gird is found, so "old style" is actually not supported
                ! 
                call ionpot_wire(clust, grid, pot, p_pot, pbc, rsymm, parallel, &
                        solver%lpole)
        case(2)
                ! pretty much the same as ionpot_wire only in 2D, including the warning.
                call ionpot_slab(clust, grid, pot, p_pot, pbc, rsymm, parallel, &
                        solver%lpole)
        case(3)
                ! Completely serial starts in g-space,
                ! with the v_first subroutine, which calculates 4 things
                !
                ! 1. vql(type,star) ionic potential per type and star (also used for forces)
                ! 2. vionz , i.e. sum(vql*structure_factor): the total ionic potential
                ! 3. dnc(type,star) core charge per type and g-star (for forces) , can be
                !     done elsewhere
                ! 4. FOR NO APPARENT REASON, now commented out: dvqj - approx for ionic
                !      pot. derivative, 
                !    
                ! transfers vionz to real space using cfftw:dfft
                ! and then changes the mapping to fit the real space grid and uses an
                ! outdated method to export the potential pot%vion to all other PEs
                call ionpbc(clust, grid, pbc, parallel, pot%vion, ipr)
        end select

else
        write(9,*) 'TEST MODE - Skipping ionic potential'
endif

  ! this is the place to write out the potential and check what's going on
  ! Write pot%vion to external file if requested.
  ! something like..
  !    call write_cube(pbc,grid,pot%vion,origin_xyz,outfile_name)

  ! ===============================================================
  ! Calculate ion-ion interaction energy and force on a ion from
  ! all other ions.
  ! ===============================================================
  if(solver%name /= TEST) then
      if (parallel%iammaster) then
         select case(pbc%per)
         case(0)
            call forceion(clust,p_pot%zion,ipr,enuc)
         case(1)
            call forceion_wire(clust,pbc%latt_vec(1,1), enuc,ipr,p_pot%zion,ierr)
         case(2)
            call ewald_slab(clust,pbc,enuc,1,p_pot%zion)
         case(3)
            call ewald_sum(clust,pbc,enuc,ipr,p_pot%zion)
         end select
      endif
      call exit_err (clust,elec_st,grid,pot,p_pot,nloc_p_pot,u_pot, &
           move, mol_dynamic,pbc,mixer,solver,symm,rsymm,parallel,&
           file_id,ierr)
  else
    write(9,*) 'TEST MODE - Skipping ion-ion, ewald etc. '
  endif

#ifdef ITAC
  call VTEND( vt_pbc_ion, vt_ierr)
  call VTBEGIN( vt_pre_scf_state, vt_ierr)
#endif
  ! ===============================================================
  ! IF INITIAL RUN
  ixtrpflg = 0
  if (istart == 0) then

     ! ===============================================================
     ! Create an initial guess for the charge density based on a
     ! superposition of the atomic charge densities. This charge
     ! density can be spin-polarized.
     ! ===============================================================
     
     if (parallel%iammaster) then
     call mysecond(ttimer0)
     endif

     call initchrg(clust,elec_st,grid,p_pot,pbc,parallel,elec_st%rho)

     if (parallel%iammaster) then
     call mysecond(ttimer1)
        write(7,*)
        write(7,20) ttimer1-ttimer0
        write(7,*)
20      format(' Charge init time [sec]:',1x,f10.2)
     endif

     ! ===============================================================
     ! IF RESTARTED RUN
  else if (istart == 1) then
     ! ===============================================================
     ! Load potential, charge density, wave functions, etc., from file
     ! ===============================================================
     !
     ! For the moment, keep interface with old wfn.dat file.
     if (readwfndat) then
        if (parallel%iammaster) then
           call restart(elec_st,grid,pbc,parallel,solver, &
                pot%vold,elec_st%rho,ixtrpflg,.false.,ierr)
        else
           call restart_mpi(elec_st,parallel,solver, &
                pot%vold,elec_st%rho,ierr)
        endif
        call exit_err (clust,elec_st,grid,pot,p_pot,nloc_p_pot,u_pot, &
             move, mol_dynamic,pbc,mixer,solver,symm,rsymm,parallel,&
       file_id,ierr)
#ifdef MPI
        call MPI_BCAST(ixtrpflg,1,MPI_INTEGER,parallel%masterid, &
             parallel%comm,mpinfo)
#endif
     else
        call restart_run(elec_st,grid,pbc,parallel,solver, &
             pot%vold,elec_st%rho,4,ixtrpflg,ierr)
        call exit_err (clust,elec_st,grid,pot,p_pot,nloc_p_pot,u_pot, &
             move, mol_dynamic,pbc,mixer,solver,symm,rsymm,parallel,&
       file_id,ierr)
        ! If restarting from previous spin-orbit calculation, update
        ! data.
        if (elec_st%is_so) then
           if (.not. p_pot%is_so) then
              if (parallel%iammaster) then
                 write(7,*) 'ERROR! restart file has spin-orbit interactions '
                 write(7,*) 'but spin-orbit potential is not present'
                 write(7,*) ' STOP '
              endif
           endif
           if (associated(u_pot%denmat)) then
              irp = size(u_pot%denmat,dim=1)
              isp = size(u_pot%denmat,dim=3)
              deallocate(u_pot%denmat)
              allocate(u_pot%zdenmat(irp,irp,isp,elec_st%nspin))
           endif
        endif
     endif
  endif
#ifdef MPI
  call MPI_Barrier(parallel%comm,mpinfo)
#endif

  thart_sum = zero
  tset_sum  = zero
  tforce_sum = zero

  if (elec_st%do_vdw) then
  call initchrg_f(clust,elec_st,grid,p_pot,pbc,parallel,elec_st%rhof,clust%dist_atm) 
  end if

  ! ===============================================================
  ! If initial run or if extrapolating from a smaller sphere size,
  ! must compute the Hartree and exchange-correlation potentials.
  ! ===============================================================

  ! istart == 0 => initial run
  if ((istart == 0) .or. (ixtrpflg == 1)) then

     ! Compute and report total number of electrons in the system.
     tdummy = sum(elec_st%rho(:,1))
     call psum(tdummy,1,parallel%group_size,parallel%group_comm)
     if (parallel%iammaster) then
     call mysecond(thart0)
        write(7,23) tdummy*grid%hcub * real(rsymm%ntrans,dp)
23      format(' Tot. electron charge from atomic density ', &
             '[e] = ',g16.5)
     endif

#ifdef ITAC
  call VTEND( vt_pre_scf_state,vt_ierr )
  call VTBEGIN( vt_hart_xc, vt_ierr)
#endif 
     call rho_har(grid,pot,rsymm,parallel,elec_st,pbc,solver%lpole,ipr,exc)

     if (parallel%iammaster) then
     ! Report timing.
     call mysecond(thart1)
        write(7,*)
        write(7,24) thart1-thart0
        thart_sum = thart_sum + thart1-thart0
        write(7,*)
24      format(' Hartree potential time [sec]:',1x,f10.2)
     endif
     ! ===============================================================
     ! Compute initial guess for the screened potential.
     ! ===============================================================
#ifdef ITAC
  call VTEND( vt_hart_xc, vt_ierr)
  call VTBEGIN( vt_pre_scf_state,vt_ierr )
#endif 
     do isp = 1, elec_st%nspin
        do ii = 1, parallel%mydim
           pot%vold(ii,isp) = pot%vion(ii)+pot%vxc(ii,isp)+pot%vhart(ii)
        enddo
     enddo
  endif
  !
  ! First run - new and old potentials (with respect to mixing) are
  ! equal.
  !
  pot%vnew(:,:) = pot%vold(:,:)
  ! ===============================================================
  ! If non-collinear calculation, compute initial guess for B_xc
  ! ===============================================================
  if (elec_st%ncl .and. istart == 0) &
     call bxc(pot,parallel,elec_st,solver,.true.,nloc_p_pot)


  tdiag_sum = zero

  ! ===============================================================
  ! Define the vector potential in case of external magnetic field.
  ! ===============================================================
!  pot%vecpot = zero
!  if(elec_st%is_mag) &
!     call magvec(rsymm,elec_st,grid,pot,parallel,solver)


  ! ===============================================================
  ! Define the number of desired eigenvalues in elec_st structure.
  ! ===============================================================
  if (.not. solver%fix_neig) call initeigval(elec_st,solver%nadd)

  ! ===============================================================
  ! ===============================================================
  ! ATOM MOVEMENT AND POLARIZABILITY LOOPS START HERE !!
  ! ===============================================================
  ! ===============================================================
  ! Code returns to this point everytime atoms move (applicable
  ! when molecular dynamics or structural relaxation are applied).
  ! If (and only if) calculating polarizability, the code cycles
  ! over seven directions of the electric field:
  ! zero, x, -x, y, -y , z, -z
  ! The code returns to this point for each of these seven
  ! calculations.
  ! ===============================================================

#ifdef ITAC
  call VTEND ( vt_pre_scf_state, vt_ierr)
#endif

30 continue

  ! Set timers for self-consistent field and movement timings.
  if (parallel%iammaster) then
  call mysecond(tmove0)
     if (is_setup) then
        ! Report setup time.
        write(7,*)
        tset_sum = tmove0-tstrt
        write(7,32) tset_sum
        is_setup = .false.
        write(7,*)
        write(7,*) ('-',ii=1,65)
        write(7,*)
32      format('Setup time [sec] :',1x,f10.2)
        call myflush(7)
     endif
  endif

  if(solver%name == TEST) then
      if (parallel%iammaster) write(7,*) "TEST MODE COMPLETE"
      write(9,*) "TEST MODE COMPLETE"
      goto 80
  endif
  ! ===============================================================
  ! Initialize temporary arrays in matvec module.
  ! ===============================================================
  call matvec_init(clust,parallel%nwedge)

  ! ===============================================================
  ! If non self consistent calculation skip to 80
  ! ===============================================================
  if (nscf%nscf_on) goto 80

  ! ===============================================================
  ! ===============================================================
  ! ======= SELF-CONSISTENT ITERATION LOOP STARTS HERE !! =========
  ! ===============================================================
  ! ===============================================================

  ! if dynamic tolerance is used, tolerance of first
  ! iteration is fixed to 0.1 here
  if (solver%dyntol) then
     solver%dyntoler = solver%toler
     solver%toler = 1.d-1
  endif

!#ifdef AJB_DEBUG
  if (parallel%iammaster) then
     write(7,*)
     write(7,*) ' Starting SCF... '
     write(7,*)
  endif
!#endif
99 continue
  do iter = 1, mxiter
       
#ifdef ITAC
  call VTBEGIN( vt_scf_state, vt_ierr)
  ! sync time stamp
  call VTTIMESYNC(mpinfo)
#endif 
     ierr = 0
     ! ===============================================================
     ! Calculate the density matrix for LDA+U.
     ! WARNING: density matrix should be calculated using a sum over
     ! all k-points! Must symmetrize (z)denmat over the irreducible
     ! Brillouin zone.
     ! ===============================================================
     if (u_pot%maxnloc > 0 .and. iter > 1) then
        if (elec_st%cplx) then
           call zcalculate_den(elec_st,u_pot,parallel)
        else
           call calculate_den(elec_st,u_pot,parallel)
        endif
        call exit_err (clust,elec_st,grid,pot,p_pot,nloc_p_pot, &
             u_pot,move,mol_dynamic,pbc,mixer,solver,symm,rsymm,parallel,&
       file_id,ierr)
     endif

     if (parallel%iammaster) then
     call mysecond(tdiag0)
     endif
#ifdef ITAC
  call VTEND( vt_scf_state, vt_ierr)
#endif 
     call eigval(elec_st,pot,u_pot,solver,parallel,ipr,ierr)
     ! Report diagonlization timing.

#ifdef ITAC
  call VTBEGIN( vt_scf_state, vt_ierr)
#endif 

     if (parallel%iammaster) then
     call mysecond(tdiag1)
     tdiag_sum = tdiag_sum + (tdiag1-tdiag0)
     endif
     ! write(9,*)
     ! write(9,42) tdiag1-tdiag0, tdiag_sum
     ! write(9,*)
     ! call myflush(9)
     if (parallel%iammaster) then
        write(7,*)
        write(7,42) tdiag1-tdiag0, tdiag_sum
        write(7,*)
        call myflush(7)
     endif
42   format(' Diagonalization time [sec] :',1x,f10.2, &
          ',', 4x,'  tdiag_sum =', f12.2)
     call exit_err (clust,elec_st,grid,pot,p_pot,nloc_p_pot, &
          u_pot,move,mol_dynamic,pbc,mixer,solver,symm,rsymm,parallel,&
       file_id,ierr)
#ifdef MPI
     i = size(elec_st%irep)
     call MPI_ALLREDUCE(elec_st%irep,i,1,MPI_INTEGER,MPI_MAX, &
          parallel%comm,mpinfo)
#endif
     ! ===============================================================
     ! Determine state occupations.
     ! ===============================================================
88   continue
     if (parallel%iammaster) call flevel(elec_st,pbc%latt_vec,nscf,ierr)
     call exit_err (clust,elec_st,grid,pot,p_pot,nloc_p_pot,u_pot, &
          move, mol_dynamic,pbc,mixer,solver,symm,rsymm,parallel,&
       file_id,ierr)


     ! patch for restarted run with spin-orbit
 !    if (istart*iter == 1 .and. elec_st%mxwd == 2) then
 !       do isp = 1, elec_st%nspin
 !          elec_st%totel(isp) = elec_st%xele*half + &
 !             real(((-1)**isp),dp)*elec_st%net_magmom
 !       enddo
 !    endif
     ! end of patch
     ! ===============================================================
     ! Compute new noncollinear magnetization from wave-functions
     ! ===============================================================
     if(elec_st%ncl) call bxc(pot,parallel,elec_st,solver,.false., &
       nloc_p_pot)


     ! ===============================================================
     ! Calculate new charge density from the wave functions.
     ! ===============================================================
     if (parallel%iammaster) then
     call mysecond(thart0)
     endif
#ifdef ITAC
  call VTEND( vt_scf_state,vt_ierr )
  call VTBEGIN( vt_hart_xc, vt_ierr)
#endif 
     call newrho(elec_st,parallel,grid%hcub2)

     ! ===============================================================
     ! Symmetrize the charge density.
     ! ===============================================================
     if (chsym) then
        do isp = 1, 2*elec_st%nspin - 1
           if (parallel%iammaster) then
              write(7,*)
              write(7,*) ' Symmetrizing charge isp = ',isp
           endif
           call symm_scalar(grid,rsymm,symm,parallel, &
                elec_st%rho(:,isp),pbc%is_on)
        enddo
     endif
     ! ===============================================================
     ! Compute new Hartree potential.
     ! ===============================================================
     call rho_har(grid,pot,rsymm,parallel,elec_st,pbc, &
          solver%lpole,ipr,exc)

     ! Report timing.
     if (parallel%iammaster) then
     call mysecond(thart1)
        write(7,*)
        write(7,24) thart1-thart0
        thart_sum = thart_sum + thart1-thart0
        write(7,*)
        call myflush(7)
     endif
#ifdef ITAC
  call VTEND( vt_hart_xc, vt_ierr)
  call VTBEGIN( vt_scf_state,vt_ierr )
#endif 

     ! ===============================================================
     ! Compute noncollinear stencile Bxc*S
     ! ===============================================================
     if (elec_st%ncl) then !{{{
       solver%bxc = zzero
       do  i = 1,parallel%mydim
          if (abs(elec_st%spin3dn(i)) > zero) then

             solver%bxc(i,1) = half * (pot%vxc(i,1)-pot%vxc(i,2)) * &
                cmplx(elec_st%spin3d(i,3),zero) / elec_st%spin3dn(i)

             solver%bxc(i,2) = half * (pot%vxc(i,1)-pot%vxc(i,2)) * &
                cmplx(elec_st%spin3d(i,1),-elec_st%spin3d(i,2)) / &
                elec_st%spin3dn(i)
          endif
       enddo
     endif !}}}

     ! ===============================================================
     ! Save the old Hartree-exchange-correlation potential.
     ! ===============================================================
     ! The sum of the Hartee and exchange-correlation potentials from
     ! the previous iteration is used in the computation of the total
     ! energy. We save it into an array here, just before vnew is
     ! updated, so old values are still reflected.

     pot%vold = pot%vnew

     do isp = 1, elec_st%nspin
        do i = 1, parallel%mydim
           pot%vhxcold(i,isp) = pot%vnew(i,isp) - pot%vion(i)
           pot%vnew(i,isp) = pot%vhart(i) + pot%vxc(i,isp) + pot%vion(i)
        enddo
     enddo

777 format(3(f15.9,1x))

     ! ===============================================================
     ! Find total energy.
     ! ===============================================================
     call totnrg(elec_st,pot,pbc,parallel,exc,enuc,grid%hcub, &
         clust%atom_num,bdev,.FALSE.)  !<=last argument is hack for vdw_flag
     ! ===============================================================
     ! Compute the self-consistent residual error (SRE).
     ! ===============================================================
     call getsre(elec_st,pot,parallel,grid%hcub,move%num,iter)
     ! if dynamic tolerance is used, update the tolerance
     ! according to the new SRE
     if (solver%dyntol) then
        solver%toler = min(solver%toler,minval(elec_st%sre)/real(10,dp))
        solver%toler = max(solver%toler,solver%dyntoler)
        if (parallel%iammaster) then
           write(7,'(a,/)') 'Dynamic tolerance is used in eigensolver'
           write(7,'(a,f11.5)') 'New tolerance:',solver%toler
        endif
     endif
     ! ===============================================================
     ! Mix old and new potential. 
     ! ===============================================================
     ! Prepare sum and difference of up and down as inputs for the
     ! mixing, only if spin-polarized run.
     ! In case of non self consistent calculation, do not update the potential.
     select case (mixer%name)
     case (ANDERSON)
        if (elec_st%nspin == 2) then
           call spinmix(pot,mixer)
           call anderson_mix(mixer,parallel,iter,elec_st%nrep, &
                mixer%xin,mixer%xout,ierr)
           call spinunmix(pot,mixer)
        else
           call anderson_mix(mixer,parallel,iter,elec_st%nrep, &
                pot%vold,pot%vnew,ierr)
        endif
     case (BROYDEN) 
        if (elec_st%nspin == 2) then
           call spinmix(pot,mixer)
           call genbro(mixer,parallel,iter,elec_st%nrep, &
                mixer%xin,mixer%xout,ierr)
           call spinunmix(pot,mixer)
        else
           call genbro(mixer,parallel,iter,elec_st%nrep, &
                pot%vold,pot%vnew,ierr)
        endif
     case (MSECANT1)
        if (elec_st%nspin == 2) then
           call spinmix(pot,mixer)
           call msecant1_mix(mixer,parallel,iter,elec_st%nrep, &
                mixer%xin,mixer%xout,ierr)
           call spinunmix(pot,mixer)
        else
           call msecant1_mix(mixer,parallel,iter,elec_st%nrep, &
                pot%vold,pot%vnew,ierr)
        endif
     case (MSECANT2)
        if (elec_st%nspin == 2) then
           call spinmix(pot,mixer)
           call msecant2_mix(mixer,parallel,iter,elec_st%nrep, &
                mixer%xin,mixer%xout,ierr)
           call spinunmix(pot,mixer)
        else
           call msecant2_mix(mixer,parallel,iter,elec_st%nrep, &
                pot%vold,pot%vnew,ierr)
        endif
     case (MSECANT3)
        if (elec_st%nspin == 2) then
           call spinmix(pot,mixer)
           call msecant3_mix(mixer,parallel,iter,elec_st%nrep, &
                mixer%xin,mixer%xout,ierr)
           call spinunmix(pot,mixer)
        else
           call msecant3_mix(mixer,parallel,iter,elec_st%nrep, &
                pot%vold,pot%vnew,ierr)
        endif
     end select
     call exit_err (clust,elec_st,grid,pot,p_pot,nloc_p_pot,u_pot, &
     move, mol_dynamic,pbc,mixer,solver,symm,rsymm,parallel,&
       file_id,ierr)
     !if(elec_st%ncl) call nclmom(pot,parallel,elec_st,solver,nloc_p_pot,.false.)
#ifdef ITAC
  call VTEND( vt_scf_state, vt_ierr)
  call VTBEGIN( vt_post_scf_state,vt_ierr )
#endif
     ! ===============================================================
     ! Write out results for charge plots and restarts.
     ! ===============================================================
    if(enable_data_out) then
     if (mod(outflag/2,2) == 1 .and. outflag>=0) then
#ifdef AJB_DEBUG
         write(9,*) ' doing wfnsave for charge plots and restarts, ncl ', elec_st%ncl
#endif
         if (elec_st%ncl) then 
             call wfnsave(elec_st,grid, &
                 pbc,rsymm,parallel,pot%vnew,elec_st%rho,4,elec_st%spin3d)
         else
             call wfnsave(elec_st,grid, &
                 pbc,rsymm,parallel,pot%vnew,elec_st%rho,4)
         endif
     endif

     if (mod((outevflag-1)/2,2) == 1) then
         if (parallel%iammaster) then
#ifdef AJB_DEBUG
             write(7,*) ' calling eigensave for charge plots and restarts '
#endif
             call eigensave(elec_st,move%num,61,outevflag)
         endif
     endif
   endif
#ifdef ITAC
  call VTEND( vt_post_scf_state,vt_ierr )
#endif
!HERE BE DRAGONS!

     ! ===============================================================
     ! If spin-orbit potentials are used and if they were already
     ! included as perturbation, add the spin-orbit potentials to the
     ! Hamiltonian and start a fresh SCF run.
     ! ===============================================================
     if(elec_st%so_pert) then !{{{
        elec_st%so_pert = .false.
        if(elec_st%scf_so) then !{{{
           if (parallel%iammaster) then !{{{
              write(7,*) ('*', i= 1, 51)
              write(7,*) '  Spin orbit correction as a perturbation is done'
              write(7,*) ' Calculating spin orbit correction self consistently'
              write(7,*) ('*', i= 1, 51)
#ifdef OLDDEBUG
              do kk = 1, elec_st%nstate
                 do pp = 1, kpnum
                    write(87,*)elec_st%eig(1,pp,1)%en(kk)-elec_st%efermi
                 enddo
              enddo
#endif
           endif !}}}
           call destroy_eigen_solver(solver)
           !if(elec_st%is_mag) &
              !call magvec(rsymm,elec_st,grid,pot,parallel,solver)
           call create_eigen_solver(solver,grid,parallel%ldn,1,1, &
                elec_st%nstate,elec_st%nkpt,.true.,parallel%mxwd)
           if (solver%do_subsp) then
              solver%eig_init(:,:,:) = .false.
              solver%firstfilt(:,:) = .true.
           endif
           go to 99
        else
           mixer%scf_fail = .false. 
           go to 60
        endif !}}}
     endif !}}}

     ! Below 'iter > 5', condition is imposed just to provide 
     ! protection from very early decrement of solver%polym. 
     ! 'solver%polym .gt. 10' makes sure that this decrement of polym
     ! mechanism, we only use when polym is large enough. These conditions
     ! are just protections and practicalities and not parameters. 

     if ( (iter > 5) .and. (solver%polym > 10) .and.  &
          (maxval(elec_st%sre) < vconv_approach)) &
          solver%polym = solver%polym-1

     ! ===============================================================
     ! Escape iterative loop if self-consistency obtained.
     ! ===============================================================
     if ( elec_st%use_plain_sre ) then
         current_sre = maxval(elec_st%plain_sre)
     else
         current_sre = maxval(elec_st%sre)
     endif

     if (current_sre < vconv) then
        solver%polym=solver%polym_t

         if (elec_st%do_vdw) then
             call efvdw(clust,elec_st,grid,parallel,pbc,clust%vdw_forces)
             !do you mean to say that only now you want the vdw corrections? seems logical
            if (parallel%iammaster) then
              write(7,*) '    TOTAL ENERGY REPORT WITH TS-VDW CORRECTIONS'
            endif

             call totnrg(elec_st,pot,pbc,parallel,exc,enuc,grid%hcub, &
                 clust%atom_num,bdev,elec_st%do_vdw) 
         end if

        if (nloc_p_pot%is_so) then 
        ! ===============================================================
        ! If spin-orbit potential is included, go back to SCF with the 
        ! spin-orbit potential. 
        ! ===============================================================
#ifdef AJB_DEBUG
        write(9,*) ' my value for nloc_p_pot%is_so is ', nloc_p_pot%is_so
#endif
           if (.not. elec_st%is_so) then !{{{
               if (parallel%iammaster) then !{{{


#ifdef OLDDEBUG
                 do kk = 1, elec_st%nstate
                    do pp = 1, kpnum
                       write(88,*)elec_st%eig(1,pp,1)%en(kk)-elec_st%efermi, &
                            elec_st%eig(1,pp,1)%en(kk)-elec_st%efermi
                    enddo
                 enddo
#endif

                 write(7,*) ('*',i= 1, 51)
                 write(7,*) '     The scalar relativistic SCF has converged' 
                 write(7,*) ' Taking spin orbit correction as a perturbation now'
                 write(7,*) ('*',i= 1, 51)
               endif !}}}
              elec_st%is_so = .true.
              elec_st%mxwd = 2
              parallel%mxwd = 2
              !
              ! Must update the size of communication arrays.
              !
              if (parallel%iammaster) then
                 if (associated(parallel%zftmp)) deallocate(parallel%zftmp)
                 allocate(parallel%zftmp &
                      (grid%nwedge*parallel%mxwd),stat=alcstat)
                 call alccheck('parallel%zftmp', &
                      grid%nwedge*parallel%mxwd,alcstat)
              endif

              if (associated(u_pot%denmat)) then
                 irp = size(u_pot%denmat,dim=1)
                 isp = size(u_pot%denmat,dim=3)
                 deallocate(u_pot%denmat)
                 allocate(u_pot%zdenmat(irp,irp,isp,elec_st%nspin))
              endif
              !
              ! Must flush eigenvalues if chebdav/subspace is used.
              if (solver%do_subsp) then
                 solver%eig_init(:,:,:) = .false.
                 solver%firstfilt(:,:) = .true.
                 solver%ncompcnt(:,:) = 0
              endif

              elec_st%ndim = grid%ndim
              elec_st%so_pert = .true.
              call pls(elec_st,solver,nloc_p_pot,parallel)
              go to 88 
           endif !}}}

           if(elec_st%is_so .and. elec_st%scf_so) then !{{{
              if (parallel%iammaster) then
#ifdef OLDDEBUG
                 do kk = 1, elec_st%nstate
                    do pp = 1, kpnum
                       write(86,*) elec_st%eig(1,pp,1)%en(kk)- &
                            elec_st%efermi
                    enddo
                 enddo
#endif
                 write(7,*) ('*',i= 1, 51)
                 write(7,*) ' Convergence of Spin Orbit Corrected Hamiltonian'
                 write(7,*) ('*',i= 1, 51)
              endif
              mixer%scf_fail = .false.
              ierr=0
              go to 60
           endif !}}}
        else                !if (nloc_p_pot%is_so)
!           if (parallel%iammaster) write(7,*)'No Spin Orbit Correction'
        !=== OR JUST GO TO 60 if no SO is done ===!
           mixer%scf_fail = .false.
           ierr = 0
           go to 60 
        endif 
    endif

     ! AJB: 
     ! this is not so useful in its current form. commenting out for the moment
     ! if (parallel%iammaster) then
     !    inquire(file='stop_scf',exist=ioflag)
     !    if (ioflag) then
     !       ierr = -10
     !       write(7,*)
     !       write(7,*) 'WARNING: abort SCF '
     !       write(7,*)
     !    endif
     ! endif
! #ifdef MPI
     ! call MPI_BCAST(ierr,1,MPI_INTEGER, &
     !      parallel%masterid,parallel%comm,mpinfo)
! #endif
     ! if (ierr /= 0) exit

  enddo                     ! iter = 1, mxiter

  ! solver%polym reset to its original value below:
  solver%polym=solver%polym_t

  ierr = 0
  if (iter > mxiter) iter = mxiter

  !if(elec_st%ncl) call nclmom(pot,parallel,elec_st,solver,nloc_p_pot,.false.)

  ! ===============================================================
  ! ===============================================================
  ! ======== SELF-CONSISTENT ITERATION LOOP ENDS HERE !! ==========
  ! ===============================================================
  ! ===============================================================


#ifdef ITAC
  call VTBEGIN( vt_post_scf_state, vt_ierr)
#endif
  ! Declare failure of convergence if iterative loop concluded
  ! normally (maximum iterations exceeded).
  if (parallel%iammaster) then
     write(7,*)
     write(7,*) ' WARNING !!! Self-consistent field not achieved'
     write(7,'(a,i4,a)') '  even after',iter,' iterations'
     write(7,*)
  endif

  mixer%scf_fail = .true.

  ! Report self-consistent field timing
60 continue
  if (parallel%iammaster) then
  call mysecond(tscf)
     write(7,*)
     write(7,62) tscf-tmove0, tdiag_sum
     write(7,*)
62   format('Time for self-consistent field [sec] :',1x,f10.2, &
          ',',4x, 'tdiag_sum =', f12.2 )
  endif                     ! parallel%iammaster
  !
  ! Delete bro.* files after a successful SCF calculation
  do ii = 1,(iter-1)/(mixer%block-1) + 1
     write(fnam,'(a4,i1,a1)') 'bro.',ii,'.'
     open(10,file=fnam//idstring,form='unformatted')
     close(10,status='delete')
  enddo
  open(11,file='bro.df.'//idstring,form='unformatted')
  close(11,status='delete')
  if (solver%do_subsp) then
     solver%firstfilt(:,:) = .true.
     solver%ncompcnt(:,:) = 0
  endif

  ! ===============================================================
  ! Calculate the total charge and dipole moment of the system.
  ! ===============================================================
  if (pbc%per<=2) &
       call dipole(clust,rsymm,elec_st,grid,parallel,p_pot%zion)

  ! ===============================================================
  ! Write out results for charge plots and restarts.
  ! ===============================================================
!HERE BE DRAGONS!
if(enable_data_out) then
  if(elec_st%ncl) then
#ifdef AJB_DEBUG
      write(9,*) ' doing wfnsave for charge plots and restarts, with spin3d, after dipole noncl '
#endif
     call wfnsave(elec_st,grid,pbc,rsymm,parallel, &
          pot%vnew,elec_st%rho,4,elec_st%spin3d)
#ifdef USEHDF5
     call h5_wfnsave(elec_st,grid,pbc,rsymm,parallel, &
          pot%vnew,elec_st%rho,4,file_id,elec_st%spin3d)
#endif
  else
#ifdef AJB_DEBUG
      write(9,*) ' doing wfnsave for charge plots and restarts, no spin3d, after dipole noncl '
#endif
     call wfnsave(elec_st,grid,pbc,rsymm,parallel, &
          pot%vnew,elec_st%rho,4)
#ifdef USEHDF5
     call h5_wfnsave(elec_st,grid,pbc,rsymm,parallel, &
          pot%vnew,elec_st%rho,4,file_id)
#endif
  endif


  if ((mod((outevflag-1)/2,2) == 0).and.(outevflag /= 0)) then
     if (parallel%iammaster) then
#ifdef AJB_DEBUG
      write(7,*) ' calling eigensave'
#endif
          call eigensave(elec_st,move%num,61,outevflag)
#ifdef AJB_DEBUG
      write(7,*) ' finished eigensave'
#endif
      endif
  endif
endif


  ! ===============================================================
  ! Compute polarizability. If solving with an electric field, use
  ! the solution without an electric field as an initial guess.
  ! ===============================================================
  if (npolflg == 1) then

     call polar(rsymm,elec_st,grid,pot,parallel,clust%atom_num, &
          field,ifield)
     ! AJB: 
     ! this is not so useful in its current form. commenting out for the moment
         ! if (parallel%iammaster) then
         !    inquire(file='stop_scf',exist=ioflag)
         !    if (ioflag) then
         !       ierr = -10
         !       write(7,*)
         !       write(7,*) 'WARNING: abort SCF '
         !       write(7,*)
         !    endif
         ! endif
! #ifdef MPI
         ! call MPI_BCAST(ierr,1,MPI_INTEGER &
         !      ,parallel%masterid,parallel%comm,mpinfo)
! #endif

         if (ifield < 8 .and. ierr  ==  0) goto 30
     if (ifield < 8) goto 30

  endif
#ifdef ITAC
  call VTEND( vt_post_scf_state, vt_ierr)
#endif

  ! ===============================================================
  ! POLARIZABILITY LOOP ENDS HERE !!
  ! ===============================================================

  ! ===============================================================
  ! Calculate forces on each nucleus. Start with local contribution.
  ! ===============================================================
#ifdef ITAC
  call VTBEGIN( vt_forces_state, vt_ierr)
#endif
  if(parallel%iammaster) then
     call mysecond(tforce0)
  endif

  select case(pbc%per)
  case(3)
     call forpbc(clust,elec_st,grid,pbc,parallel,elec_st%rho(1,1),pot%vxc)
  case(2)
     ! A temporary call, until forloc_slab will be written (ayelet)
     call forcloc_slab(clust,grid,pbc,pot,p_pot,rsymm,parallel, &
          elec_st%rho(1,1),ipr)
  case(1)
     call forcloc_wire(clust,grid,pot,p_pot,rsymm,parallel, &
          pbc%latt_vec(1,1),elec_st%rho(1,1),ipr)
  case(0)
     call forcloc(clust,grid,pot,p_pot,rsymm,parallel, &
          elec_st%rho(1,1),ipr,oldinpformat)
  end select
  if(parallel%iammaster) then
     call mysecond(tforce1)
     tforce_sum = tforce_sum + tforce1 - tforce0
     write(7,*)
     write(7,63) tforce1-tforce0, tforce_sum
     write(7,*)
63   format('Time for local forces [sec] :',1x,f10.2, &
          ',',4x, 'tforce_sum =', f12.2 )
  endif
  ! ===============================================================
  ! Add non-local pseudopotential contribution to the force to
  ! obtain the total force.
  ! ===============================================================

#ifdef AJB_DEBUG
write(9,*) ' Adding non-local pseudopotential contribution to the force'
#endif
  if(parallel%iammaster) then
     call mysecond(tforce0)
  endif

  if ( maxval(abs(p_pot%uu)) > zero .or. maxval(abs(p_pot%jj)) > zero ) &
       call forcnloc_u(clust,elec_st,p_pot,u_pot, &
       symm,rsymm,parallel,ipr,ierr)
  call forcnloc(clust,elec_st,p_pot,nloc_p_pot, &
      symm,rsymm,parallel,pbc%is_on,ipr,ierr)
#ifdef AJB_DEBUG
write(9,*) '    ....Done with non-local force contribution!'
#endif

  if(parallel%iammaster) then
     call mysecond(tforce1)
     tforce_sum = tforce_sum + tforce1 - tforce0
     write(7,*)
     write(7,64) tforce1-tforce0, tforce_sum
     write(7,*)
64   format('Time for non-local forces [sec] :',1x,f10.2, &
          ',',4x, 'tforce_sum =', f12.2 )
  endif
  call exit_err (clust,elec_st,grid,pot,p_pot,nloc_p_pot,u_pot, &
       move, mol_dynamic,pbc,mixer,solver,symm,rsymm,parallel,&
       file_id,ierr)

#ifdef ITAC
  call VTEND( vt_forces_state, vt_ierr)
#endif
  ! Report movement timing. - why here?
  if (parallel%iammaster) then
  call mysecond(tmove1)
     write(7,*)
     write(7,65) move%num, tmove1-tmove0
     write(7,*)
     write(7,*) ('=',ii=1,65)
     write(7,*)
65   format('Time for movement',1x,i4,1x,'is [sec]:',1x,f10.2)
  endif
  !
  ! Exit if no atom movements are needed.
  !
  if ((move%is_on) .and. (move%mxmove == 0)) goto 80
  if ( (mol_dynamic%is_on) .and. (mol_dynamic%step_num == 0) ) goto 80
  ! ===============================================================
  ! Relax atomic coordinates to minimize energy.
  ! ===============================================================
#ifdef ITAC
  call VTBEGIN( vt_move_state, vt_ierr)
#endif
  if (move%is_on) then
     if (parallel%iammaster) then
!ifdef AJB_DEBUG
     write(7,*)
     write(7,*) ' calling domove'
     write(7,*)
!#endif
        call domove(clust,move,pbc,elec_st%etot,bdev,ipr,ierr)
!        call myflush(66)
!ifdef AJB_DEBUG
     write(7,*)
     write(7,*) ' finished domove'
     write(7,*)
!        call myflush(7)
!endif
     endif
     call exit_err (clust,elec_st,grid,pot,p_pot,nloc_p_pot,u_pot, &
          move, mol_dynamic,pbc,mixer,solver,symm,rsymm,parallel,&
       file_id,ierr)
#ifdef MPI
     call MPI_BCAST(move%done,1,MPI_LOGICAL, &
          parallel%masterid,parallel%comm,mpinfo)
     call MPI_BCAST(move%min_frc_found,1,MPI_INTEGER, &
          parallel%masterid,parallel%comm,mpinfo)
     call MPI_BCAST(clust%xatm,clust%atom_num,MPI_DOUBLE_PRECISION, &
          parallel%masterid,parallel%comm,mpinfo)
     call MPI_BCAST(clust%yatm,clust%atom_num,MPI_DOUBLE_PRECISION, &
          parallel%masterid,parallel%comm,mpinfo)
     call MPI_BCAST(clust%zatm,clust%atom_num,MPI_DOUBLE_PRECISION, &
          parallel%masterid,parallel%comm,mpinfo)
#endif

#ifdef ITAC
  call VTEND( vt_move_state, vt_ierr)
#endif
     ! AJB: 
     ! this is not so useful in its current form. commenting out for the moment
     ! if (parallel%iammaster) then
     !    inquire(file='stop_scf',exist=ioflag)
     !    if (ioflag) then
     !       ierr = -10
     !       write(7,*)
     !       write(7,*) 'WARNING: SCF aborted by user!'
     !       write(7,*)
     !    endif
     ! endif
     ! #ifdef MPI
     !      call MPI_BCAST(ierr,1,MPI_INTEGER, &
     !           parallel%masterid,parallel%comm,mpinfo)
     ! #endif
     !      if (ierr /= 0) goto 80

     ! if minimum force criterion achieved, leave atom movement loop
     if (move%min_frc_found) goto 80
     ! if optimization algorithm terminated, leave atom movement loop
     if (move%done) goto 80
     ! update movement counter
     move%num = move%num + 1
  endif
  ! ===============================================================
  ! Perform molecular dynamics simulation.
  ! ===============================================================
#ifdef ITAC
       call VTBEGIN( vt_move_state, vt_ierr)
#endif
  if (mol_dynamic%is_on .and. parallel%iammaster) then
       call moldyn(clust,mol_dynamic,pbc,move%num,bdev)
  endif
#ifdef MPI
  call MPI_BCAST(move%num,1,MPI_INTEGER, &
       parallel%masterid,parallel%comm,mpinfo)
  call MPI_BCAST(clust%xatm,clust%atom_num,MPI_DOUBLE_PRECISION, &
       parallel%masterid,parallel%comm,mpinfo)
  call MPI_BCAST(clust%yatm,clust%atom_num,MPI_DOUBLE_PRECISION, &
       parallel%masterid,parallel%comm,mpinfo)
  call MPI_BCAST(clust%zatm,clust%atom_num,MPI_DOUBLE_PRECISION, &
       parallel%masterid,parallel%comm,mpinfo)
#endif
#ifdef ITAC
       call VTEND( vt_move_state, vt_ierr)
#endif
     ! AJB: 
     ! this is not so useful in its current form. commenting out for the moment
      ! if (parallel%iammaster) then
      !    inquire(file='stop_scf',exist=ioflag)
      !    if (ioflag) then
      !       ierr = -10
      !       write(7,*)
      !       write(7,*) 'WARNING: abort SCF '
      !       write(7,*)
      !    endif
      ! endif
! #ifdef MPI
      ! call MPI_BCAST(ierr,1,MPI_INTEGER, &
      !      parallel%masterid,parallel%comm,mpinfo)
! #endif
      ! if (ierr /= 0) goto 80
  !
  ! If no more atom movements are needed - leave the atom
  ! movemement loop.
  if ((.not. move%is_on) .and. (.not.mol_dynamic%is_on)) goto 80
  if ((move%is_on) .and. (move%num > move%mxmove)) goto 80
  if ( (mol_dynamic%is_on) .and. &
       (move%num > mol_dynamic%step_num) ) goto 80

  ! ===============================================================
  ! ===============================================================
  ! ===============================================================
  ! ===============================================================

  ! In preparation for the next atom movement, do things that need
  ! to be done again for the new atomic coordinates, before
  ! repeating the self-consistent cycle for the new coordinates.

  ! ===============================================================
  ! Check for isolated atoms.
  ! ===============================================================
#ifdef AJB_DEBUG
     if (parallel%iammaster) then
     write(7,*)
     write(7,*) ' now calling isolat'
     write(7,*)
     endif
#endif
#ifdef ITAC
       call VTBEGIN( vt_post_scf_state, vt_ierr )
#endif
  if (parallel%iammaster) then 
      call mysecond(tset0)
      call isolat(clust,grid,pbc%per,ierr)
  endif
  call exit_err (clust,elec_st,grid,pot,p_pot,nloc_p_pot,u_pot, &
       move, mol_dynamic,pbc,mixer,solver,symm,rsymm,parallel,&
       file_id,ierr)

#ifdef ITAC
       call VTEND( vt_post_scf_state, vt_ierr )
#endif
  ! ===============================================================
  ! Calculate the non-local contribution to the Hamiltonian.
  ! ===============================================================
#ifdef AJB_DEBUG
     if (parallel%iammaster) then
     write(7,*)
     write(7,*) ' now calling nonloc'
     write(7,*)
     endif
#endif

#ifdef ITAC
  call VTBEGIN( vt_nonloc, vt_ierr)
#endif
  call nonloc(clust,grid,p_pot,nloc_p_pot,pbc,parallel,ierr)
  call exit_err (clust,elec_st,grid,pot,p_pot,nloc_p_pot,u_pot, &
       move, mol_dynamic,pbc,mixer,solver,symm,rsymm,parallel,&
       file_id,ierr)
#ifdef ITAC
  call VTEND( vt_nonloc, vt_ierr)
  call VTBEGIN( vt_grid_setup_state, vt_ierr)
#endif
#ifdef AJB_DEBUG
     if (parallel%iammaster) then
     write(7,*)
     write(7,*) ' now calling upot'
     write(7,*)
     endif
#endif
  call upot(clust,grid,p_pot,u_pot,pbc,parallel,elec_st%nspin, &
       elec_st%cplx,ierr)
  call exit_err (clust,elec_st,grid,pot,p_pot,nloc_p_pot,u_pot, &
       move, mol_dynamic,pbc,mixer,solver,symm,rsymm,parallel,&
       file_id,ierr)

  ! Update matvec module
  call matvec_init(clust,parallel%nwedge)

  ! ===============================================================
  ! Compute core-correction charge density.
  ! ===============================================================
  call corecd(clust,grid,p_pot,pbc,parallel,elec_st%rhoc,ierr)
  call exit_err (clust,elec_st,grid,pot,p_pot,nloc_p_pot,u_pot, &
       move, mol_dynamic,pbc,mixer,solver,symm,rsymm,parallel,&
       file_id,ierr)

  ! ===============================================================
  ! Update the potential.
  ! ===============================================================

  ! Remove old ionic part from potential.
  do isp=1,elec_st%nspin
     do i = 1, parallel%mydim
        pot%vnew(i,isp) = pot%vnew(i,isp)-pot%vion(i)
     enddo
  enddo
#ifdef ITAC
  call VTEND( vt_grid_setup_state, vt_ierr)
  call VTBEGIN( vt_pbc_ion, vt_ierr)
#endif
  ! Calculate the new local ionic pseudopotential.
  select case(pbc%per)
  case (0)
     call ionpot(clust,grid,p_pot,parallel,pot%vion,oldinpformat)
  case(1)
     call ionpot_wire(clust,grid,pot,p_pot,pbc,rsymm,parallel, &
          solver%lpole)
  case(2)
     call ionpot_slab(clust,grid,pot,p_pot,pbc,rsymm,parallel, &
          solver%lpole)
  case(3)
     call ionpbc(clust,grid,pbc,parallel,pot%vion,ipr)
  end select

  ! Add the new ionic potential component.
  do isp = 1, elec_st%nspin
     do ii = 1, parallel%mydim
        pot%vnew(ii,isp) = pot%vnew(ii,isp) + pot%vion(ii)
     enddo
  enddo
  !
  ! Calculate ion-ion interaction energy and force on a ion from
  ! all other ions.
  if (parallel%iammaster) then
     select case(pbc%per)
     case(0)
        call forceion(clust,p_pot%zion,ipr,enuc)
     case(1)
        call forceion_wire(clust,pbc%latt_vec(1,1), enuc,ipr,p_pot%zion,ierr)
     case(2)
        call ewald_slab(clust,pbc,enuc,ipr,p_pot%zion)
     case(3)
        call ewald_sum(clust,pbc,enuc,ipr,p_pot%zion)
     end select
     call mysecond(tset1)
     tset_sum = tset_sum + (tset1 -tset0)
  endif                     ! parallel%iammaster
#ifdef ITAC
  call VTEND( vt_pbc_ion, vt_ierr)
#endif
  call exit_err (clust,elec_st,grid,pot,p_pot,nloc_p_pot,u_pot, &
       move, mol_dynamic,pbc,mixer,solver,symm,rsymm,parallel,&
       file_id,ierr)

  if ((move%is_on) .and. (move%num <= move%mxmove)) goto 30
  if( (mol_dynamic%is_on) .and.  &
       (move%num <= mol_dynamic%step_num) )  goto 30

  ! ===============================================================
  ! ===============================================================
  ! ATOM MOVEMENT LOOP ENDS HERE !!
  ! ===============================================================
  ! ===============================================================



80 continue

#ifdef ITAC
  call VTBEGIN( vt_post_processing_state, vt_ierr)
#endif
! band structure calculation

 if (band_st%bands_on) then
       if (parallel%iammaster) write(7,*) 'calling calc_bands '
       call calc_bands(band_st, elec_st, pbc, &
       pot, u_pot, solver, parallel, grid, nloc_p_pot, &
       clust, p_pot, ipr, ierr)
       if (parallel%iammaster) write(7,*) 'band structure calculation done'
 endif

! non self consistent calculation 

 if (nscf%nscf_on) then
       if (parallel%iammaster) write(7,*) 'calling calc_nscf '
       call calc_nscf(nscf, elec_st, pbc, &
       pot, u_pot, solver, parallel, grid, nloc_p_pot, &
       clust, p_pot, ipr, ierr)
       if (parallel%iammaster) write(7,*) '...NSCF calculation done'
 endif

#ifdef GW
  if (outputgw) then

! Make sure vxc has been calculated and is stored in pot%vxc. This is especially
! true if one does non-scf calclations but even in the scf case, dont shift the 
! vxc potential to zero

     !
     !  Compute the exchange-correlation potential.
     !
     !  Vxc is used as core-corrected charge on input
     !  and exchange-correlation potential on output.
     !
     allocate(vect(parallel%nwedge+1),stat=alcstat)
     call alccheck('vect',parallel%nwedge+1,alcstat)

     if (elec_st%nspin == 1) then
        do i = 1, parallel%mydim
           pot%vxc(i,1) = elec_st%rho(i,1) + elec_st%rhoc(i)
        enddo
        call exc_nspn(grid,parallel,elec_st%nrep,elec_st%icorr &
          ,elec_st%dioniz,pot%vxc,vect,tmp)
     else
        do isp = 1, elec_st%nspin
          idx = isp+1
          do i = 1, parallel%mydim
             pot%vxc(i,isp) = elec_st%rho(i,idx) + half*elec_st%rhoc(i)
          enddo
        enddo
        call exc_spn(grid,parallel,elec_st%nrep,elec_st%icorr &
          ,elec_st%dioniz,pot%vxc,vect,tmp)
     endif
   
     deallocate(vect)

     kgrid_tmp(:) = elec_st%mpgrid(:)
     kgrid_shift_tmp(:) = elec_st%mpshift(:)
     
     elec_st%mpgrid(:) = nscf%kgrid(:)
     elec_st%mpshift(:) = nscf%kshift(:)

     if (parallel%iammaster) write(7,*) 'calling gw_write '
     if (parallel%iammaster) write(7,*) 

     call gw_write(elec_st,clust,grid,pbc,symm,parallel,elec_st%rho,pot%vxc)

     if (parallel%iammaster) write(7,*) 'GW output is written'
     if (parallel%iammaster) write(7,*) 

     elec_st%mpgrid(:) = kgrid_tmp(:)
     elec_st%mpshift(:) = kgrid_shift_tmp(:)
  endif !outputgw
#endif 

  if (parallel%iammaster) then
     ! If minimization was performed:
     if (move%is_on) then
        ! Report maximum number of movements used if applicable.
        if (move%num > move%mxmove) then
           write(7,*) 'NOTICE: Maximum number of movements used.'
           write(66,*) 'NOTICE: Maximum number of movements used.'
           write(7,*)
        endif
        ! Report whether it was successful.
        if (move%min_frc_found .and. move%name /= MANUAL) then
           write(7,*) ' Minimum force is satisfied '
           write(7,*)
        elseif (move%name /= MANUAL) then
           write(7,*) ' WARNING: Minimum force is NOT satisfied!!'
           write(7,*)
        endif
        ! Write warning if SCF loop failure.
        if (mixer%scf_fail) then
           write(7,*) ' WARNING: SCF convergence not achieved' 
           write(7,*)
        endif

        ! Write last set of coordinates for next run.
        write(66,*)
        write(66,*) ' Final coordinates for next run '
        do j = 1,clust%atom_num
           write(66,84) clust%xatm(j),clust%yatm(j),clust%zatm(j)
        enddo
     endif
84   format(3(2x,f11.6))
     ! writing the density of states file
     if (pbc%create_dos .and. parallel%iammaster) then
        if (elec_st%mpgrid(1) /= elec_st%dos_mpgrid(1) .or. elec_st%mpgrid(2) /= elec_st%dos_mpgrid(2) &
        .or. elec_st%mpgrid(3) /= elec_st%dos_mpgrid(3)) then 
           write(7,*) 'Warning ! You are creating the DOS is created using the original MP Grid'
           write(7,*) 'In order to use the DOS grid perfrom a restart run with'
           write(7,*) 'the same input'
        endif
        call dos(clust,elec_st,grid,pbc,parallel,ierr)
     endif
  endif !parallel%iammaster

  ! Write grid data to external files if requested.
  if (export_griddata_flag(1) /= 0) then
     call export_grid_data(export_griddata_flag,parallel,pbc,grid,elec_st,pot)
  endif

  if (parallel%iammaster) then
     ! Report timing
     call mysecond(tfinish)
     write(7,*)
     write(7,*) ('=',ii=1,65)
     write(7,86) ' Setup time [sec] : ',tset_sum
     write(7,86) ' Time spent on Hartree potential [sec] : ',thart_sum
     write(7,86) ' Time spent on diagonalization [sec] : ',tdiag_sum
     write(7,86) ' Time spent on force calc. [sec] : ',tforce_sum
     write(7,*)
     write(7,86) ' CPU time [sec] :',tfinish-tstrt
     call custom_date_time(datelabel,wfinish)
     write(7,86) ' Wall-clock time [sec] :',wfinish-wstrt
     write(7,*) 'Current date/time: ',datelabel,' UTC'
     write(7,*)
     write(7,*) ('=',ii=1,65)
     write(7,*)
86   format(a,1x,f10.2)

  endif !parallel%iammaster

  write(9,*)
  write(9,*) 'Closing file on PE #',parallel%iam
  write(9,*)

#ifdef ITAC
  call VTEND(vt_post_processing_state, vt_ierr)
#endif
 call destroy_eigen_solver(solver)
           call create_eigen_solver(solver,grid,parallel%ldn, &
                elec_st%nrep,1,elec_st%nstate,elec_st%nkpt, &
                elec_st%cplx,parallel%mxwd)
           if (solver%do_subsp) then
              solver%eig_init(:,:,:) = .false.
              solver%firstfilt(:,:) = .true.
           endif
!AJB: This is already above, to be removed:
! if (band_st%bands_on) then
!	write(7,*) 'calling calc_bands ' 
!     call calc_bands(band_st, elec_st, pbc, & 
!      pot, u_pot, solver, parallel, grid, nloc_p_pot, & 
!     clust, p_pot, ipr, ierr) 
! endif
  call exit_err (clust,elec_st,grid,pot,p_pot,nloc_p_pot,u_pot, &
       move,mol_dynamic,pbc,mixer,solver,symm,rsymm,parallel,file_id,0)

! AJB: FAIL SAFE -  will this work?
#ifdef MPI
     call MPI_FINALIZE(mpinfo)
#endif
end program parsec
!===============================================================
