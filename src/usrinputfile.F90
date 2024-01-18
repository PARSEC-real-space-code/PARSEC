!===============================================================
!
! Copyright (C) 2005 Finite Difference Research Group
! This file is part of parsec, http://www.ices.utexas.edu/parsec/
!
! Read in user parameters from parsec.in.
!
!---------------------------------------------------------------
subroutine usrinput(clust,band_st,elec_st,grid,pot,p_pot,mol_dynamic,pbc,mixer, &
     solver,move,istart,mxiter,vconv,vconv_approach,npolflg,field,outflag, &
     ipr,export_griddata_flag,chsym,ngroups,nnodes,wstrt,nscf,oldinpformat,outevflag,readwfndat, &
     ignoresym,outputgw,enable_data_out,file_id,ierr)

  use constants
  use cluster_module
  use electronic_struct_module
  use grid_module
  use potential_module
  use pseudo_potential_module
  use molecular_dynamic_module
  use pbc_module
  use mixer_module
  use eigen_solver_module
  use movement_module
  use bandstruc_module
  use nscf_module
  use esdf
#ifdef USEHDF5
  use hdf5
#endif
#ifdef OMPFUN
  use omp_lib
#endif


  implicit none
  !
  ! Input/Output variables:
  !
  ! the cluster
  type(cluster), intent(inout) :: clust
  ! band structure
  type(bandstruc), intent(inout) :: band_st
  ! electronic structure
  type(electronic_struct), intent(inout) :: elec_st
  ! grid related data
  type(grid_data), intent(inout) :: grid
  ! potential related data
  type (potential), intent(inout) :: pot
  ! pseudopotential related data
  type(pseudo_potential), intent (inout) :: p_pot
  ! molecular dynamic related data
  type(molecular_dynamic), intent (inout) :: mol_dynamic
  ! periodic boundary conditions data
  type(pbc_data), intent (inout) :: pbc
  ! mixer related data
  type(mixer_data), intent (inout) :: mixer
  ! solver related data
  type(eigen_solver), intent (inout) :: solver
  ! movement data
  type(movement), intent (inout) :: move
  ! non selfconsistent data
  type(nscf_data), intent(inout) :: nscf

  ! restart flag
  integer, intent(out) :: istart
  ! maximal # of iterations
  integer, intent(out) :: mxiter
  ! convergence criterion
  real(dp), intent(out) :: vconv
  ! convergence criterion approaching
  real(dp), intent(out) ::  vconv_approach
  ! polarizability flag
  integer, intent(out) :: npolflg
  ! polarizability field
  real(dp), intent(out) :: field
  !
  ! output flag for writing parsec.dat file:
  ! mod(outflag,2) = 1 : write all calculated wave functions 
  ! mod(outflag,2) = 0 : write only the ones for occupied states
  ! mod(outflag/2,2) = 1 : write parsec.dat only after SCF is finished
  ! mod(outflag/2,2) = 0 : write parsec.dat between steps of SCF loop
  !
  integer, intent(out) :: outflag
  !
  ! output flag for additional info in parsec.out
  ! default is ipr = 0 (only basic output); higher values will
  ! cause additional information to be printed out
  ! ipr >= 1 : print out local/non-local components of force and
  !              initial eigenvalues
  ! ipr >= 2 : print out all symmetry data and additional arpack output
  ! ipr >= 3 : print out more BFGS output
  ! ipr >= 4 : print out all BFGS output
  ! ipr >= 5 : print out all g-stars (warning: _lots_ of data)
  ! ipr >= 6 : print out local potentials (warning: _lots_ of data)
  !
  integer, intent(out) :: ipr
  ! output flag for exporting of grid data to external files
  ! (size matches definition in def.h)
  integer, intent(inout) :: export_griddata_flag(MAX_EXPORT_OPTS)
  ! flag for symmetrization of charge density
  logical, intent(out) :: chsym

  ! number of groups of processors
  integer, intent(out) :: ngroups

  ! number of processors
  integer, intent(in) :: nnodes

  ! starting time
  real(dp), intent(inout) :: wstrt

  logical, intent(out) :: oldinpformat
  logical, intent(out) :: readwfndat
  !
  ! output flag for writing eigen.dat file:
  ! mod(outevflag,2) = 0 : write all calculated eigenvalues in eigen.dat
  ! mod(outevflag,2) = 1 : write only the ones for occupied states
  ! mod((outevflag-1)/2,2) = 0 : write eigen.dat only after SCF is finished
  ! mod(outevflag/2,2) = 1 : write eigen.dat between steps of SCF loop
  !
  integer, intent(out) :: outevflag
  ! true if all symmetry should be ignored
  logical, intent(out) :: ignoresym

  ! flag for GW
  logical, intent(out) :: outputgw
  ! flag for data output
  logical, intent(out) :: enable_data_out

  ! error flag; 100 < ierr < 201
  integer, intent(out) :: ierr
  !
  ! Work variables:
  !
  ! counters
  integer ity, i, itmp, ilow, iup, jj, i1, i2, i3
  real(dp) :: dtmp, tmp
  ! dummy number
  integer ntmp

  ! spin number
  integer nspin
  ! correlation type
  character (len=2) :: icorr
  ! ionization energy difference, used with cc correlation only
  real(dp) :: dioniz
  ! electron temperature (in Kelvin), tfermi
  ! [fixed occupations flag if negative]
  real(dp) :: tfermi
  ! string flag
  character (len=60) ::  strflag, shapeflag
  ! string with local components of pseudopotential
  character (len=1) :: locstr
  ! flag for selecting which atoms are allows to move
  ! flag is used for constructing mvat only.
  ! flag is NOT passed to the main code
  integer nslctflg
  ! default values for the different types
  integer integer_default
  real(dp) ::  double_default
  character string_default
  logical boolean_default

  ! boundary sphere size (atomic units) - for confined systems
  ! (zero boundary condition)
  real(dp) :: rmax
  ! parameters for describing the domain shape
  ! in zero boundary conditions
  ! these parameters are copied to grid at the end
  integer :: domain_shape
  real(dp), dimension(3) :: d_shape_param
  integer, dimension(1) :: i_shape_param

  ! scale used for atomic coordinates
  real(dp) :: cscale(3,3),vtmp(3),lattvec(3,3)
  ! grid spacing (atomic units)
  real(dp) :: h
  ! net charge in the system (in units of electron charge)
  real(dp) ncharge
  ! number of different atom types
  integer naty
  ! total actual number of atoms (from all types combined)
  ! used as counter during atom coordinate read-in 
  integer natom,natm
  ! number of replicas of unit cell in a super-cell calculation
  integer :: supercell(3)
  ! amount of vacuum in super-cell
  integer :: cell_vac(3)

  ! number of lines in block data
  integer nlines
  ! number of states
  integer nstate
  ! order of finite difference expansion
  integer norder
  integer, parameter :: mxorder = 20

  logical istart_tmp
  logical npolflg_tmp
  logical so_start

  logical do_vdw
  integer :: hirsh_cell_one , hirsh_cell_two , hirsh_cell_three
  integer :: hirsh_cell_total(1,3)
  real(dp) :: stpmax

  logical lflag, outputallstates

  ! temporary arrays for atom information
  integer, parameter :: maxatm = 100000
  integer, allocatable :: tmove(:)
  integer, allocatable :: atype(:)
  real(dp), allocatable :: acoord(:,:)
  real(dp), allocatable :: initmag(:,:)

  ! point charge variables
  ! number of point charge types
  integer :: npttyp  
  ! number of point charges
  integer :: nptchrg
  ! counter
  integer :: ipttyp
  ! charged sheet variables
  ! number of charged sheet types
  integer :: number_charged_sheet_type
  ! number of charged sheets
  integer :: number_charged_sheets
  ! counter
  integer :: icharged_sheet_type
  ! name of main output file
  character (len=200) :: parsec_out_filename
  ! date label
  character (len=26) :: datelabel
#ifdef USEHDF5
  integer(hid_t), intent(in) :: file_id
  integer(hid_t) :: filetype, space, dset
  integer, parameter :: lenline = 200
  character(len=lenline), dimension(:), allocatable :: cml
  INTEGER(HSIZE_T), DIMENSION(1:1) :: dims
  INTEGER(HSIZE_T), DIMENSION(2) :: data_dims
  INTEGER(HSIZE_T), DIMENSION(:), ALLOCATABLE :: str_len
  CHARACTER(LEN=20) :: dataset
  integer :: hdferr, linenumber, fileflag
#else
  integer, intent(in) :: file_id 
#endif
  ! Array for the esdf_split function (llength defined in esdf module)
  character(len=llength) :: flag_words(MAX_EXPORT_OPTS)
  character(32) :: format_buffer

  !---------------------------------------------------------------

  ! initializing default values
  integer_default = 0
  double_default = zero
  string_default = ' '
  boolean_default = .false.

  allocate(tmove(maxatm))
  allocate(acoord(3,maxatm))
  allocate(atype(maxatm))

#ifdef USEHDF5
  ! read in entire file
  open(unit=3,file='parsec.in',status='old')
  linenumber=1
  fileflag = 0
  do while (fileflag == 0)
          read(3,'(a)',iostat=fileflag)
          write(*,'(a)')
          if(fileflag==-1) exit
          linenumber = linenumber + 1
  enddo
  rewind(unit=3)
  dims = linenumber
  data_dims = (/lenline, linenumber/)
  allocate(cml(linenumber))
  allocate(str_len(linenumber))
  str_len = lenline
  linenumber=1
  fileflag = 0
  do while (fileflag == 0)
          read(3,'(a)',iostat=fileflag) cml(linenumber)
          write(*,'(a)') cml(linenumber)
          if(fileflag==-1) exit
          linenumber = linenumber + 1
  enddo
  close(3)
  ! creating our own kind of string here, I suppose?
  CALL H5Tcopy_f(H5T_STRING, filetype, hdferr)
  ! modifying it to pad a certain kind of way
  CALL H5Tset_strpad_f(filetype, H5T_STR_NULLPAD_F, hdferr)
  CALL h5screate_simple_f(1, dims, space, hdferr)
  CALL h5dcreate_f(file_id, "parsec.in", filetype, space, dset, hdferr)
  CALL h5dwrite_vl_f(dset, filetype, cml, data_dims, str_len, hdferr, space)
  CALL h5dclose_f(dset , hdferr)
  CALL H5Tclose_f(filetype, hdferr)
  CALL h5sclose_f(space, hdferr)
#endif
  

  ! initialize error flag
  ierr = 0

  ! initialization of domain_shape parameters for cluster calculations
  domain_shape = 0
  d_shape_param = zero
  i_shape_param = 0

  ! User parameter input from file 'parsec.in'	
  ! Progress output written into 'parsec.out' (this name can be
  ! modified by the user)
  call esdf_init ('parsec.in')

  ! Initialize output:
  parsec_out_filename = esdf_string('Output_File_Name','parsec.out')
  open(7,file=trim(parsec_out_filename),form='formatted')  
  write(7,*)
  write(7,*) ('=',i=1,65)
  write(7,*)
  write(7,*) ' PARSEC 1.4_DEV - Real-space DFT Program'
  write(7,*) ' under development - use at your own risk !!'
  write(7,*) ' '

 ! Pack in the build details to the parsec.out file, so have some idea in the
 ! future of what was used to generate the files. In between the GIT number
 ! mach file, and defines, we should have a pretty good idea of what was
 ! underpinning this
 write(7,*) ' begin BUILD_DETAILS '
 write(format_buffer,'("(a17,a",i03,")")') len(trim(build_gitcommit()))
 write(7,format_buffer) ' git commit  : ',trim(build_gitcommit())
 write(format_buffer,'("(a17,a",i03,")")') len(trim(build_config()))
 write(7,format_buffer) '  config file  : ',trim(build_config())
 write(7,*) ' end BUILD_DETAILS '

#ifdef MPI
  write(7,*) ' parallel run - using ', nnodes,' PEs '
#else
  write(7,*) ' serial run, no MPI interface'
#endif
#ifdef OMPFUN
  write(7,*) ' OMP enabled. Max. Threads/PE: ',OMP_GET_MAX_THREADS()
  write(7,*) ' DEBUG: Runtime OMP thread affinity set to ',OMP_GET_PROC_BIND()
  write(7,*) ' where...'
  write(7,*)  "  false state stat ="  , omp_proc_bind_false
  write(7,*)  "  true  state stat ="   , omp_proc_bind_true
  write(7,*) ' alternatively...'
  write(7,*)  "  master   binding =" , omp_proc_bind_master
  write(7,*)  "  close    binding ="  , omp_proc_bind_close
  write(7,*)  "  spread   binding =" , omp_proc_bind_spread
#endif
  call custom_date_time(datelabel,wstrt)
  write(7,'(/,a,a,a,/)') ' starting run on ',datelabel,' UTC'
  write(7,*) ('=',i=1,65)
  write(7,*)
  call myflush(7)

  ! read in restart flag: 
  ! .false. for initial run, .true. for restaring from parsec.dat
  istart_tmp = esdf_boolean ('Restart_run', boolean_default)
  if (istart_tmp) then
     istart = 1
     write(7,*) ' Restarted Run - from parsec.dat'
  else
     istart = 0
     write(7,*) 'Initial Run - starting from atomic potentials'
  endif

  ! specify wfn.dat as restarting file
  readwfndat = esdf_boolean ('Restart_from_wfndat', boolean_default)
  if (readwfndat) then
     write(7,*) ' WARNING! Restarting from file parsec.dat '
     write(7,*) ' You should think seriously about ','upgrading source code!'
  endif

  ! should symmetry be ignored
  ignoresym = esdf_boolean ('Ignore_Symmetry', boolean_default)

  write(7,*) 'ignoresym=',ignoresym

  ! read in Langevin molecular dynamics flag 
  mol_dynamic%is_on = esdf_boolean ('Molecular_Dynamics',boolean_default)

  write(7,*)
  ! read in boundary condition type 
  write(7,*) 'Grid data:'
  write(7,*) '~~------~~'
  write(7,*)

  ! read in (boundary sphere size / side of box), and grid spacing
  strflag = esdf_reduce(esdf_string('Boundary_Conditions','cluster'))
  select case(trim(strflag))
  case ('cluster','0d')
     pbc%per = 0
     pbc%is_on = .false.
     elec_st%nkpt=0
  case ('wire','tube','1d')
     pbc%per = 1
     pbc%is_on = .true.
  case ('slab','2d')
     pbc%per = 2
     pbc%is_on = .true.
  case ('bulk','3d','pbc')
     pbc%per = 3
     pbc%is_on = .true.
  case default
     write(7,'(/,a,a)') ' ERROR IN BOUNDARY CONDITIONS! STOP. ',strflag
     ierr = 101
     return
  end select

  select case(pbc%per)
  case (3)
#if defined USEFFTW2 || USEFFTW3
     if (esdf_block('Cell_Shape',nlines)) then
        if ( nlines == 1) then
           write(7,*) ' Assume orthorhombic periodic cell'
           pbc%latt_vec = zero
           read (block_data(1),*,iostat=ierr) (pbc%latt_vec(i,i),i=1,3)
           goto 30
           !that's so mature guys.
        endif
        if ( nlines /= 3) then
           write(7,*) 'ERROR: there should be 1 or 3 lines instead of ', nlines
           write(7,*) 'in the "Cell_Shape" block for 3D PBC'
           write(7,*) 'STOP in usrinput'
           ierr = 102
           return
        endif
        do i = 1, nlines
           read (block_data(i),*,iostat=ierr) pbc%latt_vec(:,i)
        enddo
30      continue
        if (ierr /= 0) then
           write(7,*) 'ERROR: bad specification of "Cell_Shape" '
           write(7,*) ' block for periodic system! '
           write(7,'(10(a,/))') block_data(1:nlines)
           write(7,*) 'STOP in usrinput'
           ierr = 103
           return
        endif
     else
        write(7,*) 'ERROR: specification of "Cell_Shape" block'
        write(7,*) ' for periodic system not found! '
        write(7,*) 'STOP in usrinput'
        ierr = 104
        return
     endif
     dtmp = esdf_physical ('Lattice_Vector_Scale',one,'bohr')
     pbc%latt_vec = dtmp * pbc%latt_vec
#else
     write(7,*) 'ERROR! Periodic boundary conditions are input '
     write(7,*) 'in parsec.in but this executable does not seem '
     write(7,*) 'to have support for FFTW. Stop. Hammer time.'
     ierr = 105
     return
#endif
  case (2)
     if (esdf_block('Cell_Shape',nlines)) then
        if ( nlines /= 2) then
           write(7,*) 'ERROR: there should be 2 lines instead of ',nlines
           write(7,*) 'in the "Cell_Shape" block'
           write(7,*) 'STOP in usrinput'
           ierr = 106
           return
        endif
        pbc%latt_vec = zero
        do i = 1, nlines
                read (block_data(i),*,iostat=ierr) pbc%latt_vec(:,i)
        enddo
        if (ierr /= 0) then
           write(7,*) 'ERROR: bad specification of "Cell_Shape" '
           write(7,*) ' block for periodic system! '
           write(7,'(10(a,/))') block_data(1:nlines)
           write(7,*) 'STOP in usrinput'
           ierr = 107
           return
        endif
     else
        write(7,*) 'ERROR: specification of "Cell_Shape" block'
        write(7,*) ' for slab system not found! '
        write(7,*) 'STOP in usrinput'
        ierr = 108
        return
     endif

     dtmp = esdf_physical ('Lattice_Vector_Scale',one,'bohr')
     pbc%latt_vec = dtmp * pbc%latt_vec

     rmax = esdf_physical ('Boundary_Sphere_Radius',double_default,'bohr')
     if (rmax == double_default) then
        write(7,*) 'ERROR: specification of boundary length for'
        write(7,*) ' non-periodic direction not found! '
        write(7,*) 'STOP in usrinput'
        ierr = 109
        return
     endif

  case (1)
     if (esdf_block('Cell_Shape',nlines)) then
        if ( nlines /= 1) then
           write(7,*) 'ERROR: there should be 1 line instead of ',nlines
           write(7,*) 'in the "Cell_Shape" block'
           write(7,*) 'STOP in usrinput'
           ierr = 106
           return
        endif
        pbc%latt_vec = zero
        read (block_data(1),*,iostat=ierr) pbc%latt_vec(1,1)
        if (ierr /= 0) then
           write(7,*) 'ERROR: bad specification of "Cell_Shape" '
           write(7,*) ' block for periodic system! '
           write(7,'(10(a,/))') block_data(1:nlines)
           write(7,*) 'STOP in usrinput'
           ierr = 107
           return
        endif
     else
        write(7,*) 'ERROR: specification of "Cell_Shape" block'
        write(7,*) ' for slab system not found! '
        write(7,*) 'STOP in usrinput'
        ierr = 108
        return
     endif
     dtmp = esdf_physical ('Lattice_Vector_Scale',one,'bohr')
     pbc%latt_vec = dtmp * pbc%latt_vec

     rmax = esdf_physical ('Boundary_Sphere_Radius',double_default,'bohr')
     if (rmax == double_default) then
        write(7,*) 'ERROR: specification of boundary length for'
        write(7,*) ' non-periodic direction not found! '
        write(7,*) 'STOP in usrinput'
        ierr = 109
        return
     endif
  case (0)
     shapeflag = esdf_reduce(esdf_string('Cluster_Domain_Shape','undef'))
     select case(trim(shapeflag))
     case ('sphere')
        domain_shape = 0
        rmax = esdf_physical('Boundary_Sphere_Radius',double_default,'bohr')
        d_shape_param(1) = rmax

        if (rmax == double_default) then
           write(7,*) 'ERROR: specification of boundary radius for'
           write(7,*) ' finite system not found! '
           write(7,*) 'STOP in usrinput'
           ierr = 110
           return
        endif

        if (rmax <= zero) then
           write(7,*) 'ERROR: spherical radius must be positive'
           write(7,*) 'STOP in usrinput'
           ierr = 110
           return
        endif
     case ('ellipsoid')
        domain_shape = 1
        if(esdf_block('Domain_Shape_Parameters', nlines)) then
           if(nlines /= 1) then
              write(7,*) "ERROR: ellipsoidal radii should be specified"
              write(7,*) " in one line: x_radius y_radius z_radius"
              write(7,*) 'STOP in usrinput'
              ierr = 110
              return
           endif

           read(block_data(1),*,iostat=ierr) d_shape_param(1), &
              d_shape_param(2), d_shape_param(3)

           if(ierr /= 0) then
              write(7,*) "ERROR: bad specification of ellipsoidal"
              write(7,*) " boundary conditions. Should be:"
              write(7,*) " x_radius y_radius z_radius"
              write(7,*) 'STOP in usrinput'
              ierr = 110
              return
           endif

           if(d_shape_param(1) <= zero .or. d_shape_param(2) <= zero &
                 .or. d_shape_param(3) <= zero) then
              write(7,*) "ERROR: bad specification of ellipsoidal radii"
              write(7,*) 'STOP in usrinput'
              ierr = 110
              return
           endif
        else
           write(7,*) "ERROR: A Domain_Shape_Parameters block is"
           write(7,*) " required for ellipsoidal cluster domains."
           write(7,*) 'STOP in usrinput'
           ierr = 110
           return
        endif
        rmax = maxval(d_shape_param)
     case ('cylinder')
        domain_shape = 2
        if(esdf_block('Domain_Shape_Parameters', nlines)) then
           if(nlines /= 1) then
              write(7,*) "ERROR: cylindrical boundaries should be specified"
              write(7,*) " in one line: radius length orientation"
              write(7,*) " -- orientation is 1 for x axis, 2 for y, 3 for z"
              write(7,*) 'STOP in usrinput'
              ierr = 110
              return
           endif

           read(block_data(1),*,iostat=ierr) d_shape_param(1), &
              d_shape_param(2), i_shape_param(1)

           if(ierr /= 0) then
              write(7,*) "ERROR: bad specification of cylindrical"
              write(7,*) " boundary conditions. Should be:"
              write(7,*) " radius length orientation"
              write(7,*) 'STOP in usrinput'
              ierr = 110
              return
           endif

           if(d_shape_param(1) <= zero) then
              write(7,*) "ERROR: bad specification of cylindrical"
              write(7,*) " radius (should be positive)"
              write(7,*) 'STOP in usrinput'
              ierr = 110
              return
           endif

           if(d_shape_param(2) <= zero) then
              write(7,*) "ERROR: bad specification of cylindrical"
              write(7,*) " length (should be positive)"
              write(7,*) 'STOP in usrinput'
              ierr = 110
              return
           endif

           if(i_shape_param(1) < 1 .or. i_shape_param(1) > 3) then
              write(7,*) "ERROR: bad specification of cylindrical"
              write(7,*) " orientation. Should be 1 for x-axis,"
              write(7,*) " 2 for y-axis, 3 for z-axis"
              write(7,*) 'STOP in usrinput'
              ierr = 110
              return
           endif
        else
           write(7,*) "ERROR: A Domain_Shape_Parameters block is"
           write(7,*) " required for cylindrical cluster domains."
           write(7,*) 'STOP in usrinput'
           ierr = 110
           return
        endif
        rmax = sqrt((d_shape_param(2)/2.0d0)**2 + d_shape_param(1)**2)
     case ('box')
        domain_shape = 3
        if(esdf_block('Domain_Shape_Parameters', nlines)) then
           if(nlines /= 1) then
              write(7,*) "ERROR: box boundaries should be specified"
              write(7,*) " in one line: x_length y_length z_length"
              write(7,*) 'STOP in usrinput'
              ierr = 110
              return
           endif

           read(block_data(1),*,iostat=ierr) d_shape_param(1), &
              d_shape_param(2), d_shape_param(3)

           if(ierr /= 0) then
              write(7,*) "ERROR: bad specification of box"
              write(7,*) " boundary conditions. Should be:"
              write(7,*) " x_length y_length z_length"
              write(7,*) 'STOP in usrinput'
              ierr = 110
              return
           endif

           if(d_shape_param(1) <= zero .or. d_shape_param(2) <= zero &
                 .or. d_shape_param(3) <= zero) then
              write(7,*) "ERROR: bad specification of box"
              write(7,*) " dimensions.  All lengths should be positive."
              write(7,*) 'STOP in usrinput'
              ierr = 110
              return
           endif
        else
           write(7,*) "ERROR: A Domain_Shape_Parameters block is"
           write(7,*) " required for box cluster domains."
           write(7,*) 'STOP in usrinput'
           ierr = 110
           return
        endif
        rmax = sqrt((d_shape_param(1)/2.0d0)**2 + (d_shape_param(2)/2.0d0)**2 &
               + (d_shape_param(3)/2.0d0)**2)
     case default
        write(7,*) "WARNING: Cluster Calculation Domain Shape Unspecified"
        write(7,*) " Assuming Spherical Shape"
        domain_shape = 0
        rmax = esdf_physical('Boundary_Sphere_Radius',double_default,'bohr')
        if (rmax == double_default) then
           write(7,*) 'ERROR: specification of boundary radius for'
           write(7,*) ' finite system not found! '
           write(7,*) 'STOP in usrinput'
           ierr = 110
           return
        endif
        d_shape_param(1) = rmax
     end select ! end of domain_shape

  case default
     write(7,*) 'ERROR: unspecified boundary conditions! '
     write(7,*) 'STOP in usrinput ',pbc%per
     ierr = 111
     return
  end select

  ! copy the domain_shape variables to the grid structure for use later
  grid%domain_shape = domain_shape
  grid%d_shape_param = d_shape_param
  grid%i_shape_param = i_shape_param

  !
  ! Add parameters for super-cell.
  !
  supercell = 1
  if (esdf_block('Super_Cell',nlines)) then
     if ( nlines /= 1) then
        write(7,*) 'ERROR: there should be 1 line instead of ',nlines
        write(7,*) 'in the "Super_Cell" block'
        write(7,*) 'STOP in usrinput'
        ierr = 158
        return
     endif
     read (block_data(1),*,iostat=ierr) supercell(1:pbc%per)
     if (ierr /= 0) then
        write(7,*) 'ERROR: bad specification of "Super_Cell" '
        write(7,*) ' block for periodic system! '
        write(7,'(10(a,/))') block_data(1:nlines)
        write(7,*) 'STOP in usrinput'
        ierr = 159
        return
     endif
     write(7,*) ' Periodic cell is being replicated according to super-cell'
     write(7,*) ' specification. Super-cell size = ',supercell(1:pbc%per)
     do i = 1, pbc%per
        pbc%latt_vec(:,i) = pbc%latt_vec(:,i) * real(supercell(i),dp)
     enddo
  endif
  !
  ! Add vacuum in super-cell.
  !
  cell_vac = 0
  if (esdf_block('Super_Cell_Vac',nlines)) then
     if ( nlines /= 1) then
        write(7,*) 'ERROR: there should be 1 line instead of ',nlines
        write(7,*) 'in the "Super_Cell_Vac" block'
        write(7,*) 'STOP in usrinput'
        ierr = 160
        return
     endif
     read (block_data(1),*,iostat=ierr) cell_vac(1:pbc%per)
     if (ierr /= 0) then
        write(7,*) 'ERROR: bad specification of "Super_Cell_Vac" '
        write(7,*) ' block for periodic system! '
        write(7,'(10(a,/))') block_data(1:nlines)
        write(7,*) 'STOP in usrinput'
        ierr = 161
        return
     endif
     write(7,*) ' Vacuum is added to super-cell'
     write(7,*) ' amount of vacuum space = ',cell_vac(1:pbc%per)
     do i = 1, pbc%per
        pbc%latt_vec(:,i) = pbc%latt_vec(:,i) * ( one + real(cell_vac(i),dp) )
     enddo
  endif
  !
  ! Add 1 to the diagonal part of pbc%latt_vec for non-periodic directions.
  !
  do i = pbc%per + 1, 3
     pbc%latt_vec(i,i) = one
  enddo
  !
  ! Define grid shift:
  ! Non-periodic directions : shift = 0.5;
  ! Periodic directions orthogonal to all others: shift = 0.5;
  ! Periodic directions non-orthogonal to some other direction: shift = zero.
  !
  grid%shift = half
  if(ignoresym) grid%shift = zero ! AMIR - this should be checked !!!!!!!
  lattvec = matmul(transpose(pbc%latt_vec),pbc%latt_vec)
  do i1 = 1, pbc%per
     do i2 = 1, pbc%per
        if ( (i1 /= i2) .and. (abs(lattvec(i1,i2)) > zero) ) &
             grid%shift(i1) = zero
     enddo
  enddo

  do jj = 1, pbc%per
     lattvec(:,jj) = pbc%latt_vec(:,jj)/ &
          (one + real(cell_vac(jj),dp))/real(supercell(jj),dp)
  enddo

  if (pbc%per == 3) rmax = one
  h = esdf_physical ('Grid_Spacing',double_default,'bohr')
  if (h == double_default) then
     write(7,*) 'ERROR: specification of grid spacing not found! '
     write(7,*) 'STOP in usrinput'
     ierr = 112
     return
  endif
  call create_grid (rmax,h,grid)

  select case(pbc%per)
  case (0)
     write(7,*) 'Confined system (cluster) with zero boundary condition!'
     select case(grid%domain_shape)
     case (0)
        write(7,*) ' Spherical domain shape:'
        write(7,'(a,f10.6,a)') ' --- Radius is ',rmax, ' bohrs'
     case (1)
        write(7,*) ' Ellipsoidal domain shape:'
        write(7,'(a,f10.6,a)') ' --- x radius is ',grid%d_shape_param(1), &
           ' bohrs'
        write(7,'(a,f10.6,a)') ' --- y radius is ',grid%d_shape_param(2), &
           ' bohrs'
        write(7,'(a,f10.6,a)') ' --- z radius is ',grid%d_shape_param(3), &
           ' bohrs'
     case (2)
        write(7,*) ' Cylindrical domain shape:'
        select case(grid%i_shape_param(1))
        case (1)
           write(7,*) ' --- Oriented along the x axis'
        case (2)
           write(7,*) ' --- Oriented along the y axis'
        case (3)
           write(7,*) ' --- Oriented along the z axis'
        end select
        write(7,'(a,f10.6,a)') ' --- Radius is ',grid%d_shape_param(1), &
           ' bohrs'
        write(7,'(a,f10.6,a)') ' --- Length is ',grid%d_shape_param(2), &
           ' bohrs'
     case (3)
        write(7,*) ' Box domain shape:'
        write(7,'(a,f10.6,a)') ' --- x length is ',grid%d_shape_param(1), &
           ' bohrs'
        write(7,'(a,f10.6,a)') ' --- y length is ',grid%d_shape_param(2), &
           ' bohrs'
        write(7,'(a,f10.6,a)') ' --- z length is ',grid%d_shape_param(3), &
           ' bohrs'
     end select
     write(7,'(a,f9.6,a)') ' Grid spacing is ', h, ' bohrs'
     itmp=0
  case (1)
     pbc%box_size = zero
     pbc%box_size(1) = abs(pbc%latt_vec(1,1))
         
     write(7,*) 'Tube boundary conditions!' 
     write(7,*) 'System is non-periodic on the yz plane'
     write(7,*) 'and periodic along x direction.'
     write(7,*)
     write(7,*) 'Periodicity length [bohr] : '
     write(7,'(3(f15.10, 1x))') pbc%latt_vec(1,1)
     write(7,*)
     write(7,'(a,f10.6,a)') ' Boundary cylinder radius is ', rmax, ' bohrs'

     write(7,'(a,f9.6,a)') ' User-provided grid spacing is ', h,' bohrs'
     write(7,*) 'WARNING: grid spacing may be rescaled below!!'
     write(7,*)

  case (2)
     pbc%box_size = zero
     do i = 1, 2
        pbc%box_size(i) = sqrt(sum(pbc%latt_vec(:,i)**2))
     enddo

     write(7,*) 'Periodic boundary conditions!'
     write(7,*) 'System is periodic on the xy plane'
     write(7,*) 'and non-periodic on the z direction.'
     write(7,*)
     write(7,*) 'Unit lattice vectors [bohr] : '
     do i = 1, 3
        write(7,'(3(f15.10, 1x))') pbc%latt_vec(:,i)
     enddo
     write(7,*)
     write(7,'(a,f10.6,a)') ' Boundary slab height is ', rmax, ' bohrs'

     write(7,'(a,f9.6,a)') ' User-provided grid spacing is ', h,' bohrs'
     write(7,*) 'WARNING: grid spacing may be rescaled below!!'
     write(7,*)

  case (3)
     do i = 1, 3
        pbc%box_size(i) = sqrt(sum(pbc%latt_vec(:,i)**2))
     enddo

     write(7,*) 'Periodic boundary conditions!' 
     write(7,*) 'System is fully periodic '
     write(7,*) 'Unit lattice vectors [bohr] : '
     do i = 1, 3
        write(7,'(3(f15.10, 1x))') pbc%latt_vec(:,i)
     enddo
     write(7,*)

     write(7,'(a,f9.6,a)') ' User-provided grid spacing is ',h,' bohrs'
     write(7,*) 'WARNING: grid spacing may be rescaled below!!'
     write(7,*)

      !Are we explicitly calculating the hartree potential?
      elec_st%explicit_hartree_pbc = esdf_boolean('Calc_Hartree_G',boolean_default)
      if (elec_st%explicit_hartree_pbc) then
         write(7,*)
         write(7,*) 'CAUTION: The Hartree potential will be calculated in reciprocal space'
         write(7,*)
      endif

  end select

  itmp = 1
  grid%ndouble = esdf_integer('Double_grid_order',itmp)
  if (grid%ndouble > 1) then
     write(7,'(a,i4)') ' Order of double grid is ',grid%ndouble
     write(7,'(a,/)') ' WARNING! DOUBLE GRID WILL MODIFY FORCES AND ', &
     'TOTAL ENERGY.'
  endif
  if (grid%ndouble < 1) then
     write(7,*) 'ERROR! double grid order must be greater than zero! ', &
         grid%ndouble
     write(7,*) 'STOP in usrinput'
     ierr = 113
     return
  endif
  
  write(7,*) 'Grid points are shifted from origin! '
  write(7,'(a,3f8.4,a)') ' shift vector = ',grid%shift, &
       '   [units of grid spacing]'

  ! Read information about kpoint sampling.

    band_st%bands_on = .false. ! this flag is here to make sure that 
                               ! it is off in cluster runs.  

    nscf%nscf_on = .false.     ! this flag is here to make sure that 
                               ! it is off in cluster runs -- although
                               ! later it might make sense to reseruct it
  if(pbc%is_on) then

     strflag = esdf_reduce(esdf_string ('kpoint_method', 'none'))

     select case(trim(strflag))
     case ('none')
        write(7,*) ' Kpoint sampling is off.'
        write(7,*) ' Calculation at gamma point'
        elec_st%nkpt = 0
     case ('mp')
        elec_st%kptmethod = MONKHORST_PACK
        elec_st%nkpt = 1
        write(7,*) ' kpoints selected by Monkhort-Pack algorithm'
        if(esdf_block('Monkhorst_Pack_Grid',nlines))then
           if(nlines /= 1) then
              write(7,*) ' error in specifying Monkhorst-Pack Grid: ',nlines
              write(7,*) ' STOP in usrinput'
              ierr = 114
              return
           endif
           elec_st%mpgrid = 1
           read(block_data(nlines),*,iostat=ierr) &
                (elec_st%mpgrid(i), i = 1, pbc%per)
           if((ierr /= 0) .or. (elec_st%mpgrid(1) <= 0) .or. &
                (elec_st%mpgrid(2) <= 0) .or. (elec_st%mpgrid(3) <= 0)) then
              write(7,*) ' error in specifying Monkhorst-Pack Grid'
              write(7,*) ' STOP in usrinput'
              ierr = 115
           endif
        endif
        write(7,*) ' Monkhorst-Pack grid dimensions specified'
            write(7,37) (elec_st%mpgrid(i), i = 1, pbc%per)
37      format(1x,i3,1x,'X',i3,1x,'X',i3)
        if(esdf_block('Monkhorst_Pack_Shift',nlines))then
           if(nlines /= 1) then
              write(7,*) ' error in specifying Monkhorst-Pack Shift'
              write(7,*) ' STOP in usrinput'
              ierr = 116
              return
           endif
           elec_st%mpshift = zero
           read(block_data(nlines),*,iostat=ierr) &
                (elec_st%mpshift(i), i = 1, pbc%per)
           if((ierr /= 0) .or. (elec_st%mpshift(1) < 0) .or. &
                (elec_st%mpshift(2) < 0) .or. (elec_st%mpshift(3) < 0)) then
              write(7,*) ' error in specifying Monkhorst-Pack Shift'
              write(7,*) ' STOP in usrinput'
              ierr = 117
           endif
        endif
        write(7,*)' Monkhorst-Pack shift specified'
        write(7,38) (elec_st%mpshift(i), i = 1, pbc%per)
38      format(1x,f7.3,1x,'X',f7.3,1x,'X',f7.3)
     case ('manual')
        elec_st%kptmethod = KP_MANUAL
        write(7,*)' kpoints selected manually'
        if(esdf_block('Kpoint_List',nlines))then
           elec_st%nkpt = nlines
           allocate(elec_st%kpts(3,elec_st%nkpt))
           allocate(elec_st%kpwt(elec_st%nkpt))
           elec_st%kpts = zero
           elec_st%kpwt = zero
           do i=1,nlines
              read(block_data(i),*) &
                   (elec_st%kpts(itmp,i), itmp = 1, 3),elec_st%kpwt(i)
           enddo
        else
           write(7,*) 'ERROR: Could not find list of kpoint coordiantes'
           write(7,*) 'STOP in usrinput.'
           ierr = 118
        endif
        strflag = esdf_reduce(esdf_string ('Kpoint_Unit', &
             'cartesian_inverse_bohr'))
        ! cscale : reciprocal unit lattice vectors; necessary only if
        ! kpoint_unit is unit_reciprocal_lattice_vectors.
        cscale(:,:) = zero
        cscale = transpose(lattvec)
        call mtrxin(cscale,vtmp(1),vtmp(2))
        cscale = cscale * twopi
        if (trim(strflag) == 'reciprocal_lattice_vectors') then
           do i = 1, elec_st%nkpt
              call matvec3('N',cscale,elec_st%kpts(1,i),vtmp)
              elec_st%kpts(:,i) = vtmp
           enddo
        endif
        write(7,*)' Kpoints (units 1/Bohr) and sampling weights'
        do i = 1, elec_st%nkpt
           write(7,'(3x,f8.5,1x,f8.5,1x,f8.5,1x,f8.5)') &
                elec_st%kpts(:,i), elec_st%kpwt(i)
           if (elec_st%kpwt(i) < 0) then
              write(7,*) 'Sampling weight must be positive'
              write(7,*) 'STOP in usrinput.'
              ierr = 119
              return
           endif
        enddo
     end select

  

    !read band structure parameters from input file   
    if(esdf_block('bandstruc',nlines))then
      band_st%bands_on = .true. 
      band_st%nlines = nlines
      allocate(band_st%blines(nlines))  
      do i=1,nlines
              read(block_data(i),*) &
                  band_st%blines(i)%line_num, band_st%blines(i)%start(1:3), &
                  band_st%blines(i)%end(1:3), band_st%blines(i)%line_name
      enddo
      itmp=20
      band_st%npoints=esdf_integer('bandstruc_points',itmp)
      write(7,*) 'Calculating band structure'
    endif

    ! Read non self consistent calculation kpoints from file
    if(esdf_block('nscf_kpoints',nlines))then
      nscf%nscf_on = .true. 
      nscf%nkpt = nlines
      allocate(nscf%kpts(3,nscf%nkpt))
      allocate(nscf%kpwt(nscf%nkpt))
      nscf%kpts = zero
      nscf%kpwt = zero
      do i=1,nlines
         read(block_data(i),*) &
               (nscf%kpts(itmp,i), itmp = 1, 3),nscf%kpwt(i)
         write(*,*) & ! this is a debug line
               (nscf%kpts(itmp,i), itmp = 1, 3),nscf%kpwt(i)
      enddo
      write(7,*) ' Calculating non self consistent electronic structure'
    endif
    if (nscf%nscf_on) then
    strflag = esdf_reduce(esdf_string ('nscf_kpoint_unit', &
             'cartesian_inverse_bohr'))
    ! cscale : reciprocal unit lattice vectors; necessary only if
    ! kpoint_unit is unit_reciprocal_lattice_vectors.
    cscale(:,:) = zero
    cscale = transpose(lattvec)
    call mtrxin(cscale,vtmp(1),vtmp(2))
    cscale = cscale * twopi
    if (trim(strflag) == 'reciprocal_lattice_vectors') then
       do i = 1, nscf%nkpt
          call matvec3('N',cscale,nscf%kpts(1:3,i),vtmp)
          nscf%kpts(:,i) = vtmp
       enddo
    endif
    ! Renormalize the weights to 1
    dtmp = sum(nscf%kpwt,1) 
    nscf%kpwt = nscf%kpwt / dtmp

    write(7,*)' Kpoints (units 1/Bohr) and sampling weights'
    do i = 1, nscf%nkpt
       write(7,'(3x,f8.5,1x,f8.5,1x,f8.5,1x,f8.5)') &
             nscf%kpts(:,i), nscf%kpwt(i)
       if (nscf%kpwt(i) < 0) then
          write(7,*) 'Sampling weight must be positive'
          write(7,*) 'STOP in usrinput.'
          ierr = 119
          return
       endif
    enddo
    endif

    if (nscf%nscf_on) then
      nscf%nstate = esdf_integer ('nscf_States_Num',elec_st%nstate)

      if (istart /= 1) then
         write(7,*) ' In case of non self consistent calculations'
         write(7,*) ' Please specify restart flag to read potential'
         write(7,*) ' from a previous run'
         write(7,*) ' STOP in usrinput'
         ierr = 114
         return
      endif
    
   endif !nscf%nscf_on

  endif ! pbc%is_on
 
  pbc%create_dos = esdf_boolean ('create_dos', .false.)
  itmp = 1000
  pbc%dos_pnum = esdf_integer('dos_pnum',itmp)

if (pbc%create_dos .and. pbc%is_on) then 
        write(7,*)' Full DOS creation is on'

        lflag = esdf_block('DOS_MP_Grid',nlines)

        if (lflag .and. nlines==1) then
                elec_st%dos_mpgrid = 1

                read(block_data(nlines),*,iostat=ierr) &
                        (elec_st%dos_mpgrid(i), i = 1, pbc%per)

                if ( (ierr /= 0) .or. &
                     (elec_st%dos_mpgrid(1) <= 0) .or. &
                     (elec_st%dos_mpgrid(2) <= 0) .or. &
                     (elec_st%dos_mpgrid(3) <= 0) ) then

                        write(7,*)' error in specifying DOS MP Grid'
                        write(7,*)' STOP in usrinput'
                        ierr=1
                endif

                write(7,*)' DOS MP grid dimensions specified'
                write(7,37) (elec_st%dos_mpgrid(i), i = 1, pbc%per)

                if(esdf_block('DOS_MP_Shift',nlines))then
                        if(nlines /= 1) then
                                write(7,*)' error in specifying DOS MP Shift'
                                write(7,*)' STOP in usrinput'
                                ierr=1
                        endif

                        elec_st%dos_mpshift = zero

                        read(block_data(nlines),*,iostat=ierr) &
                                (elec_st%dos_mpshift(i), i = 1, pbc%per)

                        if ( (ierr /= 0) .or. &
                             (elec_st%dos_mpshift(1) < 0) .or. &
                             (elec_st%dos_mpshift(2) < 0) .or. &
                             (elec_st%dos_mpshift(3) < 0)) then

                                write(7,*)' error in specifying DOS MP Shift'
                                write(7,*)' STOP in usrinput'
                                ierr=1
                        endif
                endif

                write(7,*)' DOS MP shift specified'
                write(7,38) (elec_st%dos_mpshift(i), i = 1, pbc%per)
        else
                if (elec_st%kptmethod == KP_MANUAL) then
                        write(7,*)' WARNING: DOS cannot be created based on a manual MP grid,'
                        write(7,*)'          please specify a DOS MP grid. DOS creation aborted !'
                        pbc%create_dos = .false.
                else
                        write(7,*)' DOS MP Grid not specified, using one used for density calculation'
                        elec_st%dos_mpgrid=elec_st%mpgrid
                        elec_st%dos_mpshift=elec_st%mpshift   
                endif
        endif 

        if (istart == 1) then
                elec_st%kptmethod = MONKHORST_PACK
                elec_st%mpgrid=elec_st%dos_mpgrid
                elec_st%mpshift=elec_st%dos_mpshift
                outflag = -1
        endif

        itmp = -1
        pbc%ylmdos_l = esdf_integer('ylm_dos_l',itmp)

endif

  
  ! read request for orbital magnetism
   elec_st%orb_mag = esdf_boolean ('orb_mag',.false.)  
   if(elec_st%orb_mag)then
     pot%is_cur = .true.
     elec_st%is_cur = .true.
   else
     pot%is_cur = .false.
     elec_st%is_cur = .false.
   endif

  ! read value of external magnetic field
  elec_st%is_mag = .false.
  elec_st%mag_field = esdf_double('External_Mag_Field',double_default)
  if (abs(elec_st%mag_field)>zero) then
    elec_st%mag_field = elec_st%mag_field / 2.5d+9
    elec_st%is_mag = .true.
    write(7,'(/,a)') ' Applying external uniform magnetic field'
    write(7,'(a,f11.5,a)') ' Magnetic field: H =',elec_st%mag_field, &
        ' Gauss'
  else
    elec_st%is_mag = .false.
    !write(7,'(/,a)') ' No external uniform magnetic field H=0 !!'
  endif

  ! Read user defined complex or real wfn value.
  elec_st%cplx = esdf_boolean ('complex_wfn', .false.)
  if (elec_st%nkpt > 0) elec_st%cplx = .true.
  if (elec_st%is_mag) elec_st%cplx = .true.
  !if (elec_st%mxwd == 2) elec_st%cplx = .true. !AJB? bugfix UNINITIALIZED
  if (elec_st%cplx) then
     write(7,*) 'WAVEFUNCTIONS ARE COMPLEX!'
  else
     write(7,*) 'WAVEFUNCTIONS ARE REAL!'
     if (elec_st%nkpt > 0) then
        write(7,*) '...but kpoint sampling is ON!'
        write(7,*) 'STOP in usrinput.'
        ierr = 120
        return
     endif
  endif

  ! Read order of finite difference expansion.
  itmp = 12
  norder = esdf_integer ('Expansion_Order',itmp)

  write(7,'(a,i2,/)') ' The Finite-difference expansion is of order ', norder
  if ((norder < 1) .or. (norder > mxorder)) then
     write(7,*)
     write(7,*) 'ERROR: Expansion_Order is the number of'
     write(7,*) 'neighbors used on two sides in numerical'
     write(7,*) 'derivative. It should be an integer' 
     write(7,*) 'between 2 and ', mxorder,' Input as ',norder
     write(7,*) 'STOP in usrinput.'
     ierr = 121
     return
  endif

  if (mod(norder,2) /= 0) then
     write(7,*)
     write(7,*) 'ERROR: Expansion_Order should be an even number'
     write(7,*) 'between 2 and ', mxorder,' Input as ',norder
     write(7,*) 'STOP in usrinput.'
     ierr = 122
     return
  endif

  grid%norder = norder/2

  ! Read in the number of eigenstates to be calculated, nstate;
  ! net cluster charge, ncharge.
  nstate = esdf_integer ('States_Num',integer_default)
  if (nstate == integer_default) then
     write(7,*) 'ERROR: unknown number of states. '
     write(7,*) ' Input line "States_Num" not found!'
     write(7,*) 'STOP in usrinput'
     ierr = 123
     return
  endif

  if(nstate < 0) then
      write(7,*) 'States_Num is negative. Preparing for Special TEST MODE'
      nstate = 0
  endif

  ncharge = esdf_physical ('Net_Charges',double_default,'e')
  write(7,*) 'Eigenvalue data: '
  write(7,*) '---------------- '
  write(7,*) 'Number of states: ', nstate
  write(7,'(a,g12.4,a)') ' Net cluster charge = ', ncharge,'  [e] '


  ! read in electron temperature (in Kelvin)
  dtmp = 80.0d0
  tfermi = esdf_physical ('Fermi_Temp',dtmp,'k')
  write(7,*)
  if (tfermi >= zero) then
     write(7,'(a,f9.2,a)') ' Fermi temperature = ',tfermi, ' [K]'
  else 
     write(7,*) 'WARNING: reading occupations from occup.in file'
     ! occupations read from file in initial.f
  endif

  write(7,*)
  write(7,*) 'Self-consistency data: '
  write(7,*) '----------------------'
  ! read in maximal # of iterations before giving up on self-consistency
  itmp = 50
  mxiter = esdf_integer ('Max_Iter',itmp)
  write(7,*) 'Maximum number of iterations is ', mxiter

  ! read in eigensolver
  strflag = esdf_reduce(esdf_string ('Eigensolver', 'chebdav'))
  !!-- chebdav was the default for previous version, we may set the default to chebff soon
  !strflag = esdf_reduce(esdf_string ('Eigensolver', 'chebff'))

#ifndef USEARPACK
  if (trim(strflag) == 'arpack') then
     write(7,*) ' ERROR: arpack is selected as eigensolver but'
     write(7,*) ' this executable seems to have missing arpack.'
     write(7,*) ' Stop in usrinput.'
     ierr = 124
     return
  endif
#endif

  select case (trim(strflag))
  case ('arpack')
     write(7,*) 'Using ARPACK (implicitly restarted Arnoldi eigensolver)'
     solver%name = ARPACK
     grid%experimental = esdf_boolean('Exp_Comm',boolean_default)
     solver%experimental = grid%experimental
     solver%mv_blksize = min(nstate/2,4)
     solver%mv_blksize  = esdf_integer ('Matvec_Blocksize', solver%mv_blksize)
         write(7,*) 'The matvec operations are (mostly) called to with... '
         write(7,*) '... block size: ',solver%mv_blksize
     grid%max_lap_buffers = esdf_integer('max_lap_buffers',2)
     if (grid%max_lap_buffers<2) then
         write(7,*) 'Maximum number of communication buffers in matvec: (default)',2
         grid%max_lap_buffers=2
     else
         write(7,*) 'Maximum number of communication buffers in matvec:',grid%max_lap_buffers
     endif

  case ('chebff') 
     write(7,*) 'Using first-filter solver'
     solver%name = CHEBFF
     itmp = 20
     !continue to use 'Chebdav_degree' without introducing a new variable/parameter
     solver%polym0 = esdf_integer('Chebdav_degree',itmp) 
     write(7,'(a,1x,i5)') ' Polynomial degree for First-filter is ', &
          solver%polym0
     !if the degree is between 10 and 15, can still try without erring out
     if (solver%polym0 < 10) then  
        write(7,*) ' warning: for chebff, the chebdav_degree better be at least 10'
        write(7,*) ' your chebdav_degree = ',solver%polym0
        write(7,*) ' automatically update chebdav_degree to 15 '
        solver%polym0 = 15
     endif
     solver%mv_blksize = min(nstate/2,6)
     solver%mv_blksize  = esdf_integer ('Matvec_Blocksize', solver%mv_blksize)
         write(7,*) 'The matvec operations are (mostly) called to with... '
         write(7,*) '... block size: ',solver%mv_blksize
     grid%max_lap_buffers = esdf_integer('max_lap_buffers',2)
     if (grid%max_lap_buffers<2) then
         write(7,*) 'Maximum number of communication buffers in matvec: (default)',2
         grid%max_lap_buffers=2
     else
         write(7,*) 'Maximum number of communication buffers in matvec:',grid%max_lap_buffers
     endif
     itmp = 2
     solver%ff_maxiter = esdf_integer('FF_MaxIter',itmp)
     write(7,'(a,1x,i5)') ' Maximum iteration number for First-filter is ', &
          solver%ff_maxiter
     if (solver%ff_maxiter < 1 .or. solver%ff_maxiter >= 10) then  
        write(7,*) ' warning: input FF_MaxIter is not within [1, 9] '
        write(7,*) ' your FF_MaxIter = ',solver%ff_maxiter
        write(7,*) ' automatically update ff_maxiter to 2 '
        solver%ff_maxiter = 2
     endif     

     grid%experimental = esdf_boolean('Exp_Comm',boolean_default)
     solver%experimental = grid%experimental
     !if (grid%experimental) then
     !    write(7,*) ' CAUTION! EXPERIMENTAL COMMUNICATION MODE CHOSEN!'
     ! endif

  case ('chebdav') 
     write(7,*) 'Using Chebyshev-Davidson eigensolver'
     solver%name = CHEBDAV
     itmp = 20
     solver%polym0 = esdf_integer('Chebdav_degree',itmp)
     write(7,'(a,1x,i5)') ' Polynomial degree for Chebyshev-Davidson is ', &
          solver%polym0
     if (solver%polym0 < 15) then
        write(7,*) ' ERROR: chebdav_degree should be at least 15!'
        write(7,*) ' chebdav_degree = ',solver%polym0
        write(7,*) ' STOP in usrinput.'
        ierr = 125
        return
     endif
     solver%mv_blksize = min(nstate/2,6)
     solver%mv_blksize  = esdf_integer ('Matvec_Blocksize', solver%mv_blksize)
         write(7,*) 'The matvec operations are (mostly) called to with... '
         write(7,*) '... block size: ',solver%mv_blksize
     grid%max_lap_buffers = esdf_integer('max_lap_buffers',0)
     if (grid%max_lap_buffers<2) then
         write(7,*) 'Maximum number of communication buffers in matvec: (default)',2
         grid%max_lap_buffers=2
     else
         write(7,*) 'Maximum number of communication buffers in matvec:',grid%max_lap_buffers
     endif
     grid%experimental = esdf_boolean('Exp_Comm',boolean_default)
     solver%experimental = grid%experimental
     ! if (grid%experimental) then
     !    write(7,*) ' CAUTION! EXPERIMENTAL COMMUNICATION MODE CHOSEN!'
     ! endif
  case ('diagla')
     write(7,*) 'Using DIAGLA (block preconditioned ', &
          'Lanczos eigen-problem solver)'
     solver%name = DIAGLA
        write(7,*) 'WARNING: DIAGLA eigensolver is to be removed from PARSEC'
        write(7,*) 'Please choose a different one...'
        ierr = 126
        return
  case ('trlanc','trlan')
     ! write(7,*) 'Using Thick Restarted Lanczos eigensolver'
     ! solver%name = TRLANC
     ! if (elec_st%cplx) then
        write(7,*) 'ERROR: eigensolver TRLANC removed from PARSEC'
        write(7,*) 'Please choose a different one...'
        ierr = 126
        return
     !endif
  case('none')
      write(7,*) "Eigensolver is 'none'. Preparing for Special TEST MODE"
     solver%name = TEST
  case default
     write(7,*) 'wrong input eigen solver name: ',trim(strflag)
     write(7,*) 'STOP in usrinput.'
     ierr = 127
     return
  end select
  write(7,*)

  solver%fix_neig = esdf_boolean ('Fix_Number_Eigenvalues',boolean_default)
  if (solver%fix_neig) then
     write(7,*) 'Fixing the number of eigenvalues per representation'
     write(7,*) 'WARNING!! This option should only be used when ', &
          'restarting from previous SCF.'
  endif
  write(7,*)

  solver%do_subsp = esdf_boolean ('Chebyshev_filtering',boolean_default)

  ! always do subspace filtering at latter steps
  if (solver%name == CHEBFF .or. solver%name == CHEBDAV) solver%do_subsp = .true.

  if (solver%do_subsp) then
     itmp = 15
     solver%polym = esdf_integer('Chebyshev_degree',itmp)
     solver%polym_t = solver%polym
     write(7,*) 'Performing Chebyshev subspace filtering'
     write(7,'(a,1x,i5)') ' Polynomial degree for Chebyshev filtering is ', &
          solver%polym
    if (solver%polym < 10) then
      itmp = 1
    else
      itmp = 3
    endif
    solver%dpm = esdf_integer('Chebyshev_degree_delta',itmp)
    write(7,'(a,1x,i5)') ' Change in polynomial degree (dpm) for Chebyshev filtering is ', &
          solver%dpm 
  write(7,*) 'Performing Rayleigh-Ritz rotation on all SCF iterations'
  else
     write(7,*) 'No Chebyshev subspace filtering'
  endif
  write(7,*)

  ! read in value for convergence criterion and diagonalization tolerance
  dtmp = 2.d-4
  vconv = esdf_physical ('Convergence_criterion', dtmp,'ry')
  dtmp = 100.d0 * vconv 
  vconv_approach = esdf_physical ('Convergence_criterion_approach', dtmp,'ry')

  elec_st%use_plain_sre = esdf_boolean('use_plain_sre',boolean_default) 

  dtmp = 1.d-4
  solver%toler = esdf_double('Diag_Tolerance', dtmp)

  ! read in value for dynamic convergence tolerance
  solver%dyntol = esdf_boolean('dynamic_diag_tol',boolean_default)
  ! insert dynamic tolerance for faster convergence
  if (solver%dyntol) then
     if (istart_tmp) then
        solver%dyntol = .false.
        write(7,*) 'dynamic tolerance is not working with restart run!'
     else
        solver%toler = 0.1d0
        write(7,*) 'YOU ARE USING DYNAMIC TOLERANCE FOR CONVERGENCE'
     endif
  endif

  itmp = 6
  solver%nadd = esdf_integer('Subspace_buffer_size',itmp)
  ! write(7,'(a,f12.7,a)') ' Self-consistency convergence criterion is ', &
  !      vconv, ' Ry'
  write(7,*) ' Self-consistency convergence criterion is ', &
       vconv, ' Ry'
  ! write(7,'(a,f12.7,a)') ' Self-consistency convergence indicator criterion is ', &
  !      vconv_approach, ' Ry'
  write(7,*) ' Self-consistency convergence indicator is at', &
       vconv_approach, ' Ry'
  write(7,'(a)') ' If applicable, Chebyshev polynomial order will decrease to '
  write(7,'(a)') ' 10 after SRE becomes less than the self-consistency ', &
       ' convergence indicator.'
   !write(7,'(a,f12.7)') ' Diagonalization tolerance is ',solver%toler
   if (elec_st%use_plain_sre) then
  write(7,'(a)') ''
  write(7,'(a)') '  Using the unweighted SRE for convergence checking'
  write(7,'(a)') ''
   endif
  write(7,*) ' Diagonalization tolerance is ',solver%toler
  write(7,'(a,i5,/)') ' Buffer size in subspace is  ',solver%nadd

  itmp = 30000
  solver%maxmv = esdf_integer('Maximum_matvec',itmp)

  if (solver%do_subsp .and. solver%nadd < 6) then
     solver%nadd = 6
     write(7,*) ' With Chebyshev subspace, buffer size  must be at least 6'
     write(7,'(a,1x,i5)') ' New buffer size is ',solver%nadd
  endif

  write(7,*)
  write(7,*) 'Mixer data:'
  write(7,*) '-----------'

  mixer%scf_fail = .false.  !set flag to false

  ! read in mixer type
  strflag = esdf_reduce(esdf_string ('Mixing_Method', 'anderson'))

  select case (trim(strflag))
  case ('anderson')
     mixer%name = ANDERSON
  case ('broyden')
     mixer%name = BROYDEN
  case ('msecant1')
     mixer%name = MSECANT1
  case ('msecant2')
     mixer%name = MSECANT2
  case ('msecant3')
     mixer%name = MSECANT3
  case default
     write(7,*)
     write(7,*) 'ERROR: Unknown mixer type: ',trim(strflag)
     write(7,*) 'STOP in usrinput'
     ierr = 128
     return
  end select

  itmp = 4
  mixer%memory =  esdf_integer ('Memory_Param', itmp)

  itmp = 20
  mixer%restart = esdf_integer ('Restart_Mixing', itmp)

  dtmp = 0.30d0
  mixer%param = esdf_double ('Mixing_Param', dtmp)

  itmp = 9
  solver%lpole = esdf_integer ('Solver_Lpole',itmp)
  write(7,*) 'solver lpole is : ', solver%lpole

  solver%full_hartree_flag = esdf_boolean ('Full_Hartree',.false.)
  if(solver%full_hartree_flag) then
    write(7,*) 'Full hartree boundary conditions are used'
    grid%hartree_neibs_flag=.true.
  else
    grid%hartree_neibs_flag=.false.
  endif

  

  itmp = 1
  mixer%group_size = esdf_integer ('Mixing_Group_Size', itmp)

  itmp = 0
  mixer%en_stage = esdf_integer ('Mixing_EN_Like', itmp)

  dtmp = 0.3
  mixer%restart_factor = esdf_double ('Mixing_Restart_Factor', dtmp)

  dtmp = 2.0
  mixer%expand_factor = esdf_double ('Mixing_Memory_Expand_Factor', dtmp)

  ! Preferred type of update, used by hybrid methods implemented in
  ! msecant3.f90.
  itmp = 1
  mixer%preferred_type = esdf_integer ('Mixing_Preferred_Type', itmp)

  ! Maximum memory parameter, mixer%block > mixer%memory+1
  mixer%block = 20

  if (mixer%name == ANDERSON) then
     ! read mixing parameter (Anderson mixer)
     write(7,*) 'Anderson mixer'
     write(7,'(a,f6.3,a,i2)') ' Initial Jacobian: ', mixer%param, &
          '  Mixing memory is ', mixer%memory
     write(7,*) 'Mixing restarted after ',mixer%restart,' iterations'

  else if (mixer%name == BROYDEN) then
     write(7,*) ' Generalized Broyden mixer'
     !!if (mixer%memory+1 == mixer%block) then  !!a bug, it won't error out when mixer%block is too small
     if (mixer%memory+1 > mixer%block) then
        write(7,*)
        write(7,*) 'Warning: insufficient memory for mixing.'
        !!write(7,*) 'Increase mixer%block in structures.F to at least', &
        write(7,*) '  now reset mixer%block to ', &
             mixer%memory+1
        !!write(7,*) 'STOP in usrinput'
        !!ierr = 129
        !!return  !!in fact there is no need to error out in this case, simply reset mixer%block
        mixer%block = mixer%memory+1
     endif
     write(7,'(a,f6.3,a,i2)') ' Initial Jacobian: ', mixer%param, &
          '  Mixing memory is ', mixer%memory
     write(7,*) 'Mixing restarted after ',mixer%restart,' iterations'

  else if (mixer%name == MSECANT1 .or. mixer%name == MSECANT2) then
     if (mixer%name == MSECANT1) then
        if (mixer%en_stage == 0) then
           write(7,*) 'Multisecant method, Broyden-like Type-I update'
        else
           write(7,*) 'Multisecant method, EN-like Type-I update'
        endif
     else if (mixer%en_stage == 0) then
        write(7,*) 'Multisecant method, Broyden-like Type-II update'
     else
        write(7,*) 'Multisecant method, EN-like Type-II update'
     endif
     if (mixer%group_size == 0) then
        write(7,'(a,f6.3,a,i2)') ' Mixing parameter: ', mixer%param, &
          '  All secant equations taken into account in one group'
     else
        write(7,'(a,f6.3,a,i2)') ' Mixing parameter: ', mixer%param, &
          '  Group size of secant equations: ', mixer%group_size
     endif
     write(7,'(a,f6.3)') ' Restarting factor is: ', mixer%restart_factor

  else if (mixer%name == MSECANT3) then
     if (mixer%en_stage == 0) then
        write(7,*) 'Multisecant method, Broyden-like hybrid update'
     else
        write(7,*) 'Multisecant method, EN-like hybrid update'
     endif
     if (mixer%preferred_type == 1) then
        write(7,*) 'Preferred type of update is Type-I'
     else
        write(7,*) 'Preferred type of update is Type-II'
     endif
     write(7,'(a,f6.3,a,i2)') ' Mixing parameter: ', mixer%param, &
          '  Group size of secant equations: ', mixer%group_size
     write(7,'(a,f6.3)') ' Restarting factor is: ', mixer%restart_factor
  endif
  write(7,*)

  write(7,*) 'Atom data:'
  write(7,*) '--~--~--~-'
  ! read in and write number of atom types
  naty = esdf_integer('Atom_Types_Num',integer_default)
  if (solver%name == TEST) then
      naty = 0
      ! what happens in fortran when you allocate (0)?
  else
      if (naty == integer_default) then
         write(7,*) 'ERROR: unknown number of atom types. '
         write(7,*) ' Input line "Atom_Types_Num" not found!'
         write(7,*) 'STOP in usrinput'
         ierr = 130
         return
      endif
      write(7,'(/,a,i5,/)') ' Tot. # of atom types is ', naty
  endif

  strflag = esdf_reduce(esdf_string ('Coordinate_Unit','cartesian_bohr'))
  ! cscale : for periodic systems, hold unit lattice vectors in
  ! column-wise form
  cscale(:,:) = zero

  select case (trim(strflag))
  case ('cartesian_bohr')
     write(7,*) 'Atom coordinates input in Bohr radii'
     do i = 1, 3
        cscale(i,i) = one
     enddo
  case ('cartesian_ang')
     write(7,*) 'Atom coordinates input in angstroms'
     do i = 1, 3
        cscale(i,i) = one/angs
     enddo
  case ('lattice_vectors')
     write(7,*) 'Atom coordinates input in lattice vector units'
     if (.not. pbc%is_on) then
        write(7,*) ' ERROR: this is not a periodic system'
        write(7,*) ' Lattice vectors not defined. STOP.'
        ierr = 131
        return
     endif
     cscale = lattvec
     do i = pbc%per + 1, 3
        cscale(i,i) = one
     enddo

  case default
     write(7,*)
     write(7,*) 'ERROR: unknown coordinate_unit: ',trim(strflag)
     write(7,*) 'STOP in usrinput'
     ierr = 132
     return
  end select
  write(7,*)

  clust%type_num = naty
  allocate (clust%natmi (naty))
  allocate (clust%name (naty))
  allocate (clust%ylmdos_cutoff (naty))

  ! creating the pseudo_potential module
  i = norder/2
  call create_pseudo_potential (naty,i,p_pot)
  p_pot%nlocp = 0
  p_pot%eleatm(:,:) = zero

  ! read in minimization data
  strflag = esdf_reduce(esdf_string ('Minimization', 'none'))

  select case (trim(strflag))
  case ('none')
     move%name = NONE
  case ('simple')
     move%name = STEEPDESC
  case ('bfgs')
     move%name = BFGS
  case ('manual')
     move%name = MANUAL
  case default
     write(7,*)
     write(7,*) 'ERROR: unknown minimization method: ',trim(strflag)
     write(7,*) 'STOP in usrinput'
     ierr = 133
     return
  end select

  ! For each atom type:
  natom = 0
  do ity = 1,naty
     ! read name and core cutoff
     clust%name(ity) = esdf_string ('Atom_Type', string_default)
     write(7,'(/,a,a)') ' Chemical element : ',trim(clust%name(ity))

     strflag = esdf_reduce(esdf_string ('Pseudopotential_Format','martins_new'))
!ccm
     write(7,*) strflag
!ccm

     select case (trim(strflag))
     case ('martins_old')
        p_pot%format(ity) = MARTINS_OLD
        write(7,*) ' pseudopotential format : old Martins'
     case ('martins_new')
        p_pot%format(ity) = MARTINS_NEW
        write(7,*) ' pseudopotential format : new Martins'
     case ('martins_wang')
        p_pot%format(ity) = MARTINS_WANG
        write(7,*) ' pseudopotential format : Martins-Wang'
     case ('fhipp')
        p_pot%format(ity) = FHIPP
        write(7,*) ' pseudopotential format : FHIPP'
     case default
        write(7,*) ' Unknown pseudopotential format :',trim(strflag)
        ierr = 134
        return
     end select

     ! does this atom type have spin-orbit psp?
     p_pot%so(ity) = esdf_boolean ('SO_PSP',.false.)
     if (p_pot%format(ity) == FHIPP .and. p_pot%so(ity)) then
        write(7,*) ' ERROR: FHIPP pseudopotential is input for '
        write(7,*) ' atom type ',clust%name(ity),' and spin-orbit'
        write(7,*) ' is enabled for this atom. But FHIPP does not'
        write(7,*) ' have support for spin-orbit potentials! STOP'
        ierr = 135
        return
     endif

     ! should the valence charge be read from psp for this atom type?
     p_pot%rvcd(ity) = esdf_boolean ('Read_VCD',.false.)
     if (p_pot%format(ity) /= MARTINS_NEW .and. p_pot%rvcd(ity)) then
        write(7,*) ' ERROR: check the type of pseudopotential    '
        write(7,*) ' atom type ',clust%name(ity)
        write(7,*) ' Can not read valence charge density from file'
        ierr = 234
        return
     endif

     ! Local component, read for all formats except Martins-Wang
     if ( p_pot%format(ity) /= MARTINS_WANG ) then
        locstr = esdf_reduce(esdf_string('Local_Component',string_default))

        if(locstr == 's') then
           p_pot%loc(ity) = 1
        else if (locstr == 'p') then
           p_pot%loc(ity) = 2
        else if (locstr == 'd') then
           p_pot%loc(ity) = 3
        else if (locstr == 'f') then
           p_pot%loc(ity) = 4
        else
           write(7,*)
           write(7,*) 'ERROR: problem with atom type', ity
           write(7,*) 'Impossible local pseudopotential'
           write(7,*) 'Must be s, p, d, or f'
           write(7,*) 'STOP in usrinput'
           ierr = 136
           return
        endif
     endif
     !
     ! Other parameters for this atom type: core cut-off radius,
     ! number of psp channels, number of electrons per pseudo-orbital.
     ! Read for formats Martins_old and FHIPP only.
     !
     p_pot%nlocp(ity) =  esdf_integer ('Potential_Num',integer_default)
     if ( p_pot%format(ity) == MARTINS_OLD .or. &
          p_pot%format(ity) == MARTINS_WANG .or. &
          p_pot%format(ity) == FHIPP ) then
        p_pot%rcore(ity) = esdf_physical ('Core_Cutoff_Radius', &
             double_default,'bohr')
     endif
     if ( p_pot%format(ity) == MARTINS_OLD .or. &
          p_pot%format(ity) == FHIPP ) then
        if (esdf_block('Electron_Per_Orbital',nlines)) then
           if ( nlines /= 1) then
              write(7,*) 'ERROR: there should be only one line '
              write(7,*) 'in the "Atom_Per_Orbital" block'
              write(7,*) 'STOP in usrinput'
              ierr = 137
              return
           endif
           read (block_data(nlines),*) &
                (p_pot%eleatm(ity,itmp),itmp=1,p_pot%nlocp(ity))
        else
           write(7,*) 'Cannot find #electrons/orbital block'
        endif
     endif
     !
     ! FHIPP format also need the non-linear core correction flag to
     ! be specified.
     !
     if ( p_pot%format(ity) == FHIPP ) then
        lflag = esdf_boolean('Nonlinear_Core_Correction' ,boolean_default)
        if (lflag) then
           p_pot%icore(ity) = 1
        else
           p_pot%icore(ity) = 0
        endif
     endif
     !
     ! Test sanity of pseudopotential parameters.
     ! Martins_old and FHIPP formats only.
     !
     if ( p_pot%format(ity) == MARTINS_OLD .or. &
          p_pot%format(ity) == FHIPP ) then
        if ((p_pot%nlocp(ity) <= 0) .or. (p_pot%nlocp(ity) > 4)) then
           write(7,*)
           write(7,*) 'ERROR: problem with atom type', ity
           write(7,*) 'impossible number of potentials! ',p_pot%nlocp(ity)
           write(7,*) 'STOP in usrinput.'
           ierr = 138
           return
        endif
        if (p_pot%loc(ity) > p_pot%nlocp(ity)) then
           write(7,*)
           write(7,*) 'ERROR: problem with atom type', ity
           write(7,*) 'local pseudopot. beyond # of components'
           write(7,*) 'STOP in usrinput.'
           ierr = 139
           return
        endif
     endif

     ! read initial spin polarization; this is used only in spin-polarized
     ! calculations and when not restating from previous charge density.
     dtmp = 0.1d0
     p_pot%spol(ity) = esdf_double('Initial_Spin_Polarization',dtmp)

     !this is now only available for pbc
     if(pbc%is_on) then
     ! read parameters for Fourier filtering of pseudopotentials and
     ! core charge (if it exists)
     dtmp = zero
     p_pot%alpha(ity) = esdf_double ('Alpha_filter', dtmp)
     dtmp = zero
     p_pot%beta1(ity) = esdf_double ('Beta1_filter', dtmp)
     dtmp = p_pot%alpha(ity)
     p_pot%acore(ity) = esdf_double ('Core_filter', dtmp)
     endif

     ! Use cubic spline interpolation for pseudopotentials?
     p_pot%uspline(ity) = esdf_boolean ('Cubic_spline', boolean_default)
     ! read number of atoms of this type
     nslctflg  = esdf_integer ('Move_Flag', integer_default)
     clust%ylmdos_cutoff(ity) = esdf_physical ('ylm_dos_cutoff',one,'bohr')
     ! read their coordinates
     if (esdf_block ('Atom_Coord',clust%natmi(ity) )) then   
        if(nslctflg > clust%natmi(ity)) then
           write(7,*)
           write(7,*) 'ERROR: invalid move_flag value: ',nslctflg
           write(7,*) 'STOP in usrinput.'
           ierr = 140
           return
        end if
        i = supercell(1) * supercell(2) * supercell(3)
        write(7,12) clust%natmi(ity)*i, clust%name(ity)
        write(7,*) 'and their initial coordinates are:'
        write(7,*)
        write(7,16)

        do i = 1, clust%natmi(ity)
           natom = natom + 1
           atype(natom) = ity
           if (natom > maxatm) then
              write(7,*) ' ERROR: too many atoms in input! ', &
                   natom,maxatm,' Increase maxatm in usrinput. '
              ierr = 141
              return
           endif
           if (nslctflg == 0) then
              read(block_data(i),*) acoord(:,natom)
              tmove(natom) = 1
           else if (nslctflg > 0) then
              read(block_data(i),*) acoord(:,natom)
              if (i <= nslctflg) then 
                 tmove(natom)=1
              else
                 tmove(natom)=0
              endif
           else
              if ( (move%name == NONE) .and. &
                   ( .not.(mol_dynamic%is_on) ) )  then
                 read (block_data(i),*) acoord(:,natom)
                 tmove(natom)=1
              else
                 read (block_data(i),*) acoord(:,natom),tmove(natom)
              endif
              if((tmove(natom) /= 0) .and. (tmove(natom) /= 1))then
                 write(7,*)
                 write(7,*) 'ERROR: invalid atom movement index: ',tmove(natom)
                 write(7,*) 'STOP in usrinput'
                 ierr = 142
                 return
              endif
           endif
           call matvec3('N',cscale,acoord(1,natom),vtmp)
           ! If multiple periodic cells are defined in a super-cell, augment
           ! number of atoms.
           if (supercell(1) * supercell(2) * supercell(3) /= 1) then
              jj = tmove(natom)
              natom = natom - 1
              select case(pbc%per)
              case(1)
                 do i1 = 0, supercell(1) - 1
                    natom = natom + 1
                    acoord(:,natom) = vtmp + real(i1,dp)*lattvec(:,1)
                    tmove(natom) = jj
                    write(7,18) acoord(:,natom),tmove(natom)
                 enddo
              case(2)
                 do i1 = 0, supercell(1) - 1
                    do i2 = 0, supercell(2) - 1
                       natom = natom + 1
                       acoord(:,natom) = vtmp + &
                            real(i1,dp)*lattvec(:,1) + &
                            real(i2,dp)*lattvec(:,2)
                       tmove(natom) = jj
                       write(7,18) acoord(:,natom),tmove(natom)
                    enddo
                 enddo
              case(3)
                 do i1 = 0, supercell(1) - 1
                    do i2 = 0, supercell(2) - 1
                       do i3 = 0, supercell(3) - 1
                          natom = natom + 1
                          acoord(:,natom) = vtmp + &
                               real(i1,dp)*lattvec(:,1) + &
                               real(i2,dp)*lattvec(:,2) + &
                               real(i3,dp)*lattvec(:,3)
                          tmove(natom) = jj
                          write(7,18) acoord(:,natom),tmove(natom)
                       enddo
                    enddo
                 enddo
              end select
           else
              acoord(:,natom) = vtmp
              write(7,18) acoord(:,natom),tmove(natom)
           endif
        enddo
        clust%natmi(ity) = clust%natmi(ity) * supercell(1) * &
             supercell(2) * supercell(3)
     else
        write(7,*) 'ERROR: no atom coordinates were found'
        write(7,*) 'for atom type ', clust%name(ity)
        write(7,*) 'STOP in usrinput'
        ierr = 143
        return
     endif
     ! read parameters for the LDA+U (if applicable)
     ! assume LDA+U only for d channel (lp1 = 3)
     p_pot%uu(:,ity) = zero
     p_pot%jj(:,ity) = zero
     p_pot%uu(3,ity) = esdf_physical ('ldaplusu_u',double_default,'ry')
     p_pot%jj(3,ity) = esdf_physical ('ldaplusu_j',double_default,'ry')
     if (p_pot%uu(3,ity) /= zero .or. p_pot%jj(3,ity) /= zero) &
         write(7,'(/,a,f10.3,a,f10.3,a)') ' LDA+U parameters : U = ', &
         p_pot%uu(3,ity)*rydberg,' eV   J = ',p_pot%jj(3,ity)*rydberg,' eV'

  enddo                     !do ity = 1,naty
12 format(1x,'There are ', i6, 1x, a2, 1x, 'atoms')
16 format(4x,'x [bohr]',10x,'y [bohr]',10x,'z [bohr]',5x,'movable?')
18 format(3(f15.9,3x),3x,i1)

  ! Does any chemical element have spin-orbit?
  p_pot%is_so = .false.
  do ity = 1, naty
     if (p_pot%so(ity)) p_pot%is_so = .true.
  enddo

if (p_pot%is_so) then
        !AJB: please explain why hcub==sqrt(hcub)
        p_pot%so_hcub = sqrt(grid%hcub)

        write(7,'(a,/)')'  Spin-Orbit corrected computation!'

        if (pbc%per > 0) then
                write(7,'(a,/)') '  WARNING: I have used the input grid &
                        &parameter to setup parameters for SO, but the &
                        &actual grid parameter COULD be modified &
                        &by PARSEC below'
        endif

        !     elec_st%mxwd = 1
endif

  so_start = .false.
  ! read spin-orbit flag for scf (default is perturbation theory)
  elec_st%scf_so = esdf_boolean ('SCFSO',boolean_default)
  if (elec_st%scf_so) write(7,'(a,/)') &
       ' Self-Consistent Spin-Orbit Correction!'
  allocate (p_pot%rws(naty))
  p_pot%rws(:) = zero
  do ity = 1, naty
     p_pot%rws(ity) = esdf_double('Atomic_Radius',double_default)
  enddo

  ! read spin-orbit flag for spinor representation from scratch
    elec_st%ncl =  esdf_boolean ('Non_Collinear_magnetism',.false.)
    if (elec_st%ncl) then
       elec_st%nspin = 2
       elec_st%mxwd = 2
       elec_st%cplx = .true.
       write(7,'(a,/)') 'Performing Non-Collinear Clculation!'
       if(elec_st%scf_so .or. p_pot%is_so) so_start = .true.
       allocate(initmag(3,maxatm))
       initmag = zero
       natm = 0
       allocate(elec_st%init_mag(3,maxatm))
       elec_st%init_mag(:,:) = zero
       do ity = 1,naty
          if (esdf_block ('Initial_NCL_Moment',clust%natmi(ity) )) then
             do i = 1, clust%natmi(ity)
                natm = natm + 1
                read(block_data(i),*) initmag(:,natm)
                tmp = sqrt(dot_product(initmag(:,natm),initmag(:,natm)))
                elec_st%init_mag(:,natm) = initmag(:,natm)/tmp
             enddo
          endif
       enddo
    endif
  ! read spin-orbit flag for spinor representation from scratch
  so_start = esdf_boolean ('SO_from_scratch',so_start)
  if (so_start) then
    elec_st%cplx = .true.
    elec_st%nspin = 2
    elec_st%scf_so = .true.
    elec_st%mxwd = 2
    elec_st%is_so = .true.
    write(7,'(a,/)') 'Starting Self-Consistent Spin-Orbit Clculation!'
  else
    if(.not. elec_st%ncl) elec_st%mxwd = 1
    !AJB: Nowhere in the code is elec_st%is_so set to false so I put it here for now
    elec_st%is_so= .false.
  endif

  nstate = nstate*elec_st%mxwd

  ! Make sure there are atoms with spin-orbit if p_pot%is_so = .true.
  if (p_pot%is_so) then
     lflag = .false.
     do ity = 1, naty
        if (p_pot%so(ity)) lflag = .true.
     enddo
     if (.not. lflag) then
        write(7,*) ' WARNING: spin-orbit interaction is enabled but'
        write(7,*) ' there are no atom types with spin-orbit potential'
        p_pot%is_so = .false.
     endif
  endif

  ! if (solver%name == TRLANC) then
  !    if (p_pot%is_so) then
  !       write(7,*) 'ERROR: eigensolver TRLANC not available with spin-orbit.'
  !       write(7,*) 'Chose a different one...'
  !       ierr = 144
  !       return
  !    endif
  ! endif

  ! creating the 'empty' cluster
  if( naty > 0  ) then
      call create_cluster (natom,clust)
      clust%mvat(1:natom) = tmove(1:natom)
      clust%xatm(1:natom) = acoord(1,1:natom)
      clust%yatm(1:natom) = acoord(2,1:natom)
      clust%zatm(1:natom) = acoord(3,1:natom)
      clust%atype(1:natom) = atype(1:natom)

      deallocate(tmove)
      deallocate(acoord)
      if (allocated(initmag)) deallocate(initmag)
      write(7,'(/,a,i7,/)') ' Tot. number of atoms = ', clust%atom_num
  endif
  ! read in point charge flag
  clust%has_ptchrg = esdf_boolean ('Add_Point_Charges',boolean_default)

  if (clust%has_ptchrg) then
     if (pbc%is_on) then
        write(7,*) 'ERROR: support for point charges is only available in'
        write(7,*) 'non-periodic systems. STOP in usrinput'
        ierr = 145
        return
     endif
     allocate(acoord(3,maxatm))
     write(7,*) 'Point Charge data:'
     write(7,*) '------------------'
     !  read in and write number of atom types
     npttyp = esdf_integer('Point_Typ_Num',integer_default)
     if (npttyp == integer_default) then
        write(7,*) 'ERROR: unknown number of point charge types. '
        write(7,*) ' Input line "Point_Typ_Num" not found!'
        write(7,*) 'STOP in usrinput'
        ierr = 146
        return
     endif
     write(7,*) 'Tot. # of point charge types is ', npttyp
     clust%npttyp = npttyp
     allocate (clust%nptt (npttyp))
     allocate (clust%qpt (npttyp))

     !   For each point charge type:
     nptchrg = 0
     do ipttyp = 1,npttyp
        ! read amount of charge per type
        clust%qpt(ipttyp) = esdf_physical('Pt_Chg',double_default,'e')
        ! read coordinates
        if (esdf_block ('Point_Coord',clust%nptt(ipttyp) )) then
           write(7,22) clust%nptt(ipttyp), clust%qpt(ipttyp)
           write(7,*) 'and their coordinates are:'
           write(7,*)
           write(7,26)
           do i=1,clust%nptt(ipttyp)
              nptchrg = nptchrg + 1
              if (nptchrg > maxatm) then
                 write(7,*) 'ERROR: too many point charges in input! ', &
                      'Increase maxatm in usrinputfile. '
                 ierr = 147
                 return
              endif
              read(block_data(i),*) acoord(:,nptchrg)
              call matvec3('N',cscale,acoord(1,nptchrg),vtmp)
              acoord(:,nptchrg) = vtmp
              write(7,18) acoord(:,nptchrg)
           enddo
        else
           write(7,*) 'ERROR: no point charge coordinates were found'
           write(7,*) 'for point charge type ', ipttyp
           write(7,*) 'STOP in usrinput'
           ierr = 148
           return
        endif
     enddo

22   format(1x,'There are',1x,i4,1x,'point charges with charge ',f5.2,' e')
26   format(4x,'x [bohr]',10x,'y [bohr]',10x,'z [bohr]')

     ! create appropriate arrays in cluster structure
     allocate (clust%xpt (nptchrg))
     allocate (clust%ypt (nptchrg))
     allocate (clust%zpt (nptchrg))

     clust%xpt(1:nptchrg) = acoord(1,1:nptchrg)
     clust%ypt(1:nptchrg) = acoord(2,1:nptchrg)
     clust%zpt(1:nptchrg) = acoord(3,1:nptchrg)

     clust%nptchrg = nptchrg

     write(7,*)
     write(7,'(a,i4,/)') ' Tot. number of point charges = ',nptchrg
     deallocate(acoord)
  endif ! end point charges

  ! read in charged sheet flag
  clust%has_charged_sheet = esdf_boolean ('add_charged_sheets',boolean_default)

  if (clust%has_charged_sheet) then
     if (pbc%per /= 2 ) then
        write(7,*) 'ERROR: support for charged sheets is only available in'
        write(7,*) '2-D (slab) systems. STOP in usrinput'
        ierr = 145 !TODO: choose a different error number
        return
     endif
     allocate(acoord(3,maxatm))
     write(7,*) 'Charged sheet data:'
     write(7,*) '------------------'
     !  read in and write number of charged sheet types
     number_charged_sheet_type = esdf_integer('num_charged_sheets',integer_default)
     if (number_charged_sheet_type == integer_default) then
        write(7,*) 'ERROR: unknown number of charged sheets. '
        write(7,*) ' Input line "num_charged_sheets" not found!'
        write(7,*) 'STOP in usrinput'
        ierr = 146 !TODO: choose a different error number
        return
     endif
     write(7,*) 'Tot. # of charged sheets is ', number_charged_sheet_type
     clust%number_charged_sheet_type = number_charged_sheet_type
     allocate (clust%number_charged_sheets_per_type (number_charged_sheet_type))
     allocate (clust%sheet_charge_of_type (number_charged_sheet_type))

     !   For each charged sheet type:
     number_charged_sheets = 0
     do icharged_sheet_type = 1,number_charged_sheet_type
        ! read amount of charge per type
        clust%sheet_charge_of_type(icharged_sheet_type) = esdf_physical('sheet_charge',double_default,'e')
        ! read coordinates
        if (esdf_block ('Sheet_Coord',clust%number_charged_sheets_per_type(icharged_sheet_type) )) then
           write(7,10022) clust%number_charged_sheets_per_type(icharged_sheet_type), clust%sheet_charge_of_type(icharged_sheet_type)
           write(7,*) 'and their coordinates are:'
           write(7,*)
           write(7,10026)
           do i=1,clust%number_charged_sheets_per_type(icharged_sheet_type)
              number_charged_sheets = number_charged_sheets + 1
              if (number_charged_sheets > maxatm) then
                 write(7,*) 'ERROR: too many charged sheets in input! ', &
                      'Increase maxatm in usrinputfile. '
                      ierr = 147 !TODO choose a different error number 
                 return
              endif
              read(block_data(i),*) acoord(:,number_charged_sheets)
              call matvec3('N',cscale,acoord(1,number_charged_sheets),vtmp)
              acoord(:,number_charged_sheets) = vtmp
              write(7,18) acoord(:,number_charged_sheets)
           enddo
        else
           write(7,*) 'ERROR: no charged sheet coordinates were found'
           write(7,*) 'for sheet charge number: ', icharged_sheet_type
           write(7,*) 'STOP in usrinput'
           ierr = 148 !TODO choose a different error number
           return
        endif
     enddo

10022   format(1x,'There are',1x,i4,1x,'charged sheets with charge ',f5.2,' e')
10026   format(4x,'x [bohr]',10x,'y [bohr]',10x,'z [bohr]')

     ! create appropriate arrays in cluster structure
!     allocate (clust%x_charged_sheets (number_charged_sheets))
!     allocate (clust%y_charged_sheets (number_charged_sheets))
     allocate (clust%z_charged_sheets (number_charged_sheets))

!     clust%x_charged_sheets(1:number_charged_sheets) = acoord(1,1:number_charged_sheets)
!     clust%y_charged_sheets(1:number_charged_sheets) = acoord(2,1:number_charged_sheets)
     clust%z_charged_sheets(1:number_charged_sheets) = acoord(3,1:number_charged_sheets)

     clust%number_charged_sheets = number_charged_sheets

     write(7,*)
     write(7,'(a,i4,/)') ' Tot. number of charged sheets = ',number_charged_sheets
     deallocate(acoord)
  endif ! end charged sheets

  ! Read in spin polarization flag.
  lflag = esdf_boolean ('Spin_Polarization', boolean_default)
  if (lflag) then
     nspin = 2
  else
     nspin = 1
  endif

  ! Read in correlation type. Default spin-polarized functional is
  ! LSDA Perdew-Wang. Default non-spin-polarized functional is
  ! LDA Perdew-Zunger.

  write(7,*) 'Correlation data:'
  write(7,*) '-----------------'
  strflag = esdf_reduce(esdf_string ('Correlation_Type','xx'))
  if (trim(strflag) == 'xx') then
     if (nspin == 2) then
        strflag = 'pl'
     else
        strflag = 'ca'
     endif
  endif

  select case (trim(strflag))
  case ('xa')
     icorr = 'xa'
     write(7,20) icorr,'LDA, Slater x-alpha'
  case ('wi')
     icorr = 'wi'
     write(7,20) icorr,'LDA, Wigner interpolation'
  case ('hl')
     icorr = 'hl'
     write(7,20) icorr,'LDA, Hedin-Lundqvist'
  case ('ca','pz','lda')
     icorr = 'ca'
     write(7,20) icorr,'LDA, Ceperley-Alder, Perdew-Zunger parametrization'
  case ('pl','pw92','pwlda')
     icorr = 'pl'
     write(7,20) icorr,'LDA, Perdew-Wang parametrization'
  case ('pw','pw91','pwgga')
     icorr = 'pw'
     write(7,20) icorr,'GGA, Perdew-Wang parametrization'
  case ('pb','pbe')
     icorr = 'pb'
     write(7,20) icorr,'GGA, Perdew-Burke-Ernzerhof parametrization'
  case ('lb')
     icorr = 'lb'
     write(7,20) icorr,'ACLDA, Leeuwen-Baerends correction'
  case ('cs')
     icorr = 'cs'
     write(7,20) icorr,'ACLDA, Casida-Salahub correction'
     ! read-in the ionization energy difference (actually used 
     ! only with 'cs' corr.)
     dtmp = 0.d0
     dioniz = esdf_physical('Ion_Energy_Diff', dtmp,'ry')
     write(7,'(a,f7.4,a)') ' Ionization energy difference is: ',dioniz, ' Ry'
  case default
     write(7,*)
     write(7,*) 'ERROR: unknown correlation type:',trim(strflag)
     write(7,*) 'STOP in usrinput.'
     ierr = 149
     return
  end select
20 format(1x,'Exchange-Correlation functional is',1x,a2,/,1x,a)
  write(7,*) 
!!! VDW PART
  do_vdw=esdf_boolean('Van_Der_Waals' , boolean_default)
  if (do_vdw) then
  write(7,*)'including Tkatchenko-Scheffler dispersion correction!'
      if (icorr /='pb') then
          write(7,*)'Error, the TS correction implemented only for GGA-PBE functional!'
          ierr = 158
          return
      end if
      if (.not. ignoresym) then
          write(7,*)'Error, for the TS correction symmetry has to be switched off'
          ierr = 159
          return
      end if
      if (pbc%is_on) then
          if (pbc%per/=3) then
              write(7,*)'Error, in this version the TS correction implemented only for 3D periodic boundary conditions (sorry)'
              ierr = 160
              return
          end if      
          write(7,*)'TS-VDW: with periodic calculation'
          itmp = 0
          hirsh_cell_one = esdf_integer('Hf_One' , itmp)
          hirsh_cell_two = esdf_integer('Hf_Two' , itmp)
          hirsh_cell_three = esdf_integer('Hf_Three' , itmp)

          hirsh_cell_total(1,1)=hirsh_cell_one
          hirsh_cell_total(1,2)=hirsh_cell_two
          hirsh_cell_total(1,3)=hirsh_cell_three
          write(7,*)'Got the following Hf parameters:'
          write(7,*) hirsh_cell_total
          write(7,*)'The numbers for the Hf cells are: 2*(Hf_One, Hf_Two, Hf_Three)+1 '
          write(7,*) hirsh_cell_one*2+1, hirsh_cell_two*2+1,hirsh_cell_three*2+1
          write(7,*)'Note: The thresold for the TS energy component is 10^-6 [eV]'
      else
          hirsh_cell_total(:,:)=0
      end if
  else
  !    write(7,*)'Note: The calculation does not include dispersion corrections.'
  end if  

  itmp = 0
!!! END VDW PART
  write(7,*) 'Other input data:'
  write(7,*) '-----------------'

  !check to see that data output is enabled
  enable_data_out = esdf_boolean ( 'Output_Data', .TRUE. )

  if( .NOT. enable_data_out) then
     write(7,*) 'NOTICE: Raw Data output has been completely disabled'
     write(7,*) 'NOTICE: All other output flags but GW will be ingnored'
  endif

  ! read in output flag: 
  outflag = 0
  outputallstates = esdf_boolean ('Output_All_States', boolean_default)
  if (outputallstates)  then
     ! save all eigenvectors to disk
     outflag = outflag + 1
     elec_st%nsave = nstate*elec_st%mxwd
     allocate(elec_st%indxsave(elec_st%nsave))
     do i = 1, elec_st%nsave
        elec_st%indxsave(i) = i
     enddo
     write(7,*) 'Writing all calculated ','wave functions to parsec.dat file'
  else
     ! read in set of eigenvectors to be written to disk
     if (esdf_block('save_wfn_bands',nlines)) then
        if (nlines == 0) then
           elec_st%nsave = 0
        elseif (nlines == 1) then
           read(block_data(1),*,iostat = ierr) ilow, iup
           if (ierr == 0) then
              elec_st%nsave = iup - ilow + 1
              allocate(elec_st%indxsave(elec_st%nsave))
              do i = ilow, iup
                 elec_st%indxsave(i - ilow + 1) = i
              enddo
           else
              elec_st%nsave = nlines
              allocate(elec_st%indxsave(nlines))
              do i = 1, nlines
                 read(block_data(i),*,iostat = ierr) elec_st%indxsave(i)
              enddo
           endif
        else
           elec_st%nsave = nlines
           allocate(elec_st%indxsave(nlines))
           do i = 1, nlines
              read(block_data(i),*,iostat = ierr) elec_st%indxsave(i)
           enddo
        endif
        if (elec_st%nsave > 0) then
           write(7,*) 'Writing ',elec_st%nsave, &
                ' selected wave functions to parsec.dat file:'
           write(7,'(10i6)') elec_st%indxsave(1:elec_st%nsave)
        else
           write(7,*) 'No wave functions are written to  parsec.dat file'
        endif
        if (ierr /= 0) then
           write(7,*) 'ERROR: bad specification of "Save_Wfn_Bands" block!'
           write(7,*) 'STOP in usrinput.'
           ierr = 150
           return
        endif
     else
        ! save only eigenvectors of occupied levels to disk
        elec_st%nsave = -1
        write(7,*) 'Writing only wave functions of occupied levels ', &
             'to parsec.dat file'
     endif
  endif
  write(7,*)

  outputgw = esdf_boolean ('Output_GW', boolean_default)
  if (outputgw) then
#ifdef GW
     write(7,*) 'Writing files for GW calculation using BerkeleyGW'
#else
     write(7,*) 'Please recompile the code with -DGW flag to enable this'
     write(7,*) 'and link it to appropriate libraries from BerkeleyGW'
     write(7,*) 'STOP in usrinput.'
     ierr = 150
     return
#endif
     if (nscf%nscf_on) then
      if(esdf_block('nscf_kgrid',nlines))then
       if(nlines /= 1) then
          write(7,*) ' error in specifying kgrid for nscf calc: ',nlines
          write(7,*) ' STOP in usrinput'
          ierr = 114
          return
        endif
      nscf%kgrid = 1
     read(block_data(nlines),*,iostat=ierr) &
          (nscf%kgrid(i), i = 1, 3)
      write(7,*) ' Non-self consistent kgrid dimensions'
      write(7,37) (nscf%kgrid(i), i = 1, 3)
     else
      write(7,*) ' Have to specify nscf_kgrid'
      write(7,*) ' STOP in usrinput'
      ierr = 114
      return
     endif
     if(esdf_block('nscf_kgrid_shift',nlines))then
       if(nlines /= 1) then
          write(7,*) ' error in specifying Monkhorst-Pack Shift'
          write(7,*) ' STOP in usrinput'
          ierr = 116
          return
       endif
     nscf%kshift = zero
     read(block_data(nlines),*,iostat=ierr) &
          (nscf%kshift(i), i = 1, 3)
     if((ierr /= 0) .or. (nscf%kshift(1) < 0) .or. &
        (nscf%kshift(2) < 0) .or. (nscf%kshift(3) < 0)) then
        write(7,*) ' error in specifying Monkhorst-Pack Shift'
        write(7,*) ' STOP in usrinput'
        ierr = 117
     endif
     write(7,*)' Non-selfconsistent kgrid shift'
     write(7,38) (nscf%kshift(i), i = 1, 3)
     else
      write(7,*) ' Have to specify nscf_kgrid_shift'
      write(7,*) ' STOP in usrinput'
      ierr = 114
      return
    endif
    elec_st%efermi = esdf_physical ('nscf_fermi_level',double_default,'Ry')
    if (elec_st%efermi == double_default) then
        write(7,*) 'ERROR: specification of Fermi level for'
        write(7,*) ' non self consistent calculation not found! '
        write(7,*) ' Please specify nscf_efermi taken from scf calculation '
        write(7,*) 'STOP in usrinput'
        ierr = 114
        return
     endif
   endif

  endif

  lflag = boolean_default
  if (elec_st%nkpt > 0) lflag = .true.
  chsym = esdf_boolean('Symmetrize_Charge_Density',lflag)
  if (chsym) write(7,*) 'Symmetrizing charge density according to symm.', &
       ' operations in the full point group'

  lflag = esdf_boolean ('Save_Intermediate_Charge_Density',boolean_default)
  if (lflag)  then
     write(7,*) 'Writing charge density to parsec.dat file after each ', &
          'step of SCF'
     write(7,*)
     outflag = outflag + 2
  endif


  if (nspin == 2) then

     write(7,'(a,/)') ' Spin-polarized computation!'
     write(7,*) 'Initial spin polarization: '
     do ity = 1, naty
        write(7,'(2x,a2,5x,f11.4)') clust%name(ity),p_pot%spol(ity)
     enddo
     write(7,*)
  else
     nspin = 1
     write(7,'(a,/)') ' No spin effects!'
     if (p_pot%is_so) then
        write(7,*) 'ERROR: spin-orbit requires spin polarized functional.'
        write(7,*) 'Must declare "spin_polarization .true." in parsec.in.'
        write(7,*) 'STOP'
        ierr = 151
        return
     endif
  endif

  oldinpformat = esdf_boolean ('Old_Interpolation_Format',boolean_default)
  write(7,*)
  if (oldinpformat) then
     write(7,*) 'Using old (direct) interpolation of pseudopot. !'
  else
     write(7,*) 'Using new (r*V_ps) interpolation of pseudopot. !'
  endif

  !
  ! Output level flag
  ipr = esdf_integer('Output_Level',integer_default)
  write(7,'(/,a,i2,/)') ' output level [1 - 6] = ',ipr

  !
  ! Grid data output flag
  flag_words = esdf_split(esdf_reduce(esdf_string ('output_grid_data', 'none')), &
         MAX_EXPORT_OPTS, ntmp)
  export_griddata_flag = NONE
  i = 0
  do while (i < ntmp)
     i = i+1
     strflag = flag_words(i)
     select case (trim(strflag))
     case ('none')
        i = ntmp
     case ('charge_density')
        write(7,*) 'Writing total charge density to external file'
        write(7,*)
        export_griddata_flag(i) = CHGDENS
     case ('v_external')
        write(7,*) 'Writing external potential to external file'
        write(7,*)
        export_griddata_flag(i) = VEXT
     case ('v_hartree')
        write(7,*) 'Writing hartree potential to external file'
        write(7,*)
        export_griddata_flag(i) = VHAR
     case ('v_external+v_hartree')
        write(7,*) 'Writing sum of external and hartree potentials to external file'
        write(7,*)
        export_griddata_flag(i) = VEXTHAR
     case ('all')
        write(7,*) 'Writing charge density, external potential, hartree potentials and the sum of the two to external files'
        write(7,*)
        export_griddata_flag(1) = CHGDENS
        export_griddata_flag(2) = VEXT
        export_griddata_flag(3) = VHAR
        export_griddata_flag(4) = VEXTHAR
        i = ntmp
     case default
        write(7,*) 'Wrong input to output_grid_data: ',trim(strflag)
        write(7,*) 'STOP in usrinput.'
        ierr = 1001 !TODO: choose a different error number
        return
     end select
  enddo

  ! read in  eigenstate to separated file output  flag:
  ! .false. / .true.
  outevflag=0
  lflag = esdf_boolean ('Output_Eigenvalues', boolean_default)

  if (elec_st%nkpt > 0) lflag = .true.

  if (lflag) then
     outevflag = outevflag + 1
     write(7,*) 'Generating Eigenvalue Output File!' 
     lflag = esdf_boolean ('Output_All_States', boolean_default)
     if (outputallstates)  then
        write(7,*) 'Writing all calculated eigenvalues to eigen.dat'
        write(7,*)
        outevflag = outevflag + 1
     endif
     lflag = esdf_boolean ('Save_Intermediate_Eigenvalues',boolean_default)
     if (lflag)  then
        write(7,*) 'Writing eigenvalues to eigen.dat after each step of SCF'
        write(7,*)
        outevflag = outevflag + 2
     endif
  else
     write(7,*) 'No eigenvalue output File!'
  endif
  write(7,*) 'outevflag = ', outevflag
  write(7,*)

  if (move%name == NONE) then
     move%is_on = .false.
     write(7,*) 'No minimization!'
  else
     move%is_on = .true.
  endif

  if (move%is_on) then
     itmp = 500
     move%mxmove = esdf_integer ('Movement_Num', itmp)
     dtmp = 1.d-2
     move%fmin = esdf_physical ('Force_Min', dtmp,'ry/bohr')
     stpmax = 0.05d0
     move%stepmax = esdf_physical ('Max_Step', stpmax,'bohr')
     dtmp = 1.d-4
     move%stepmin = esdf_physical ('Min_Step', dtmp,'bohr')
     move%is_restart = esdf_boolean ('relax_restart',boolean_default)
     write(7,*) 'Minimization data:'
     write(7,*) '------------------'
     if (move%fmin < 1.d-5) then
        write(7,*) 'ERROR: minimum force below numerical accuracy'
        write(7,*) 'increase it to at least 1E-5'
        write(7,*) 'STOP in usrinput'
        ierr = 152
        return
     endif
     if (move%name == STEEPDESC) then
        write(7,*) 'Simple minimization!'
        write(7,*) 'This scheme is not robust!'
        write(7,*) 'make sure the construction of displacements ', &
             'in domove.f fits your purpose.'
        if (move%stepmax < 1.d-5) then
           write(7,*) 'ERROR: maximum step size too small!'
           write(7,*) 'increase it to at least 1E-5'
           write(7,*) 'STOP in usrinput'
           ierr = 153
           return
        endif
     else if (move%name == BFGS) then
        write(7,*) 'Limited-memory BFGS minimization!'
        dtmp = 0.d0
        itmp = 7
        move%mmax = esdf_integer ('BFGS_number_corr', itmp)
        if (move%stepmax == stpmax) then
           ! unrestricted BFGS, no bounds on atom coordinates
           move%stepmax  = -1.d0
        else
           ! restricted BFGS, move%stepmax defines the maximum displacement
           write(7,*) 'Restricted BFGS'
           write(7,*) 'Maximum displacement of atoms from initial ', &
                'position is ',move%stepmax
           write(7,*) ' bohrs along each cartesian direction.'
        endif
     else if (move%name == MANUAL) then
        write(7,*) 'atoms move according to info in manual.dat'
        ! make all atoms movable
        clust%mvat = 1
        ! Set force limit to zero (this avoids undesirable outputs from
        ! subroutine domove).
        move%fmin = zero
     end if
     write(7,'(a,i4)') ' Maximum number of movements  = ',move%mxmove
     if (move%name == STEEPDESC) then
        write(7,91) ' Maximum step size  = ', move%stepmax
        write(7,91) ' Mininum step size  = ', move%stepmin
     endif
     if (move%name /= MANUAL) write(7,'(a,f8.5,a)') &
          ' Minimum force  = ',move%fmin, ' [Ry/au]'
     if (move%is_restart) then
        write(7,*) 'Restart relaxation, information on previous'
        write(7,*) 'runs read from relax_restart.dat file.'
     endif
     if (move%name == BFGS) write(7,'(a,a,i4)') ' Number of ', &
          'corrections used in LBFGS update = ',move%mmax
  end if

91 format(a,f8.5,' [au]')
  write(7,*)

  ! read in LMD parameters
  ! initial and final temperatures (in kelvin), 
  ! step temperature for stair cooling
  ! # of step, cooling scheme (1 - stair, 2- linear, 3 - log)
  ! time step, friction coefficient, # of time steps
  strflag = esdf_reduce(esdf_string ('Cooling_Method','stair'))

  select case (trim(strflag))
  case ('stair')
     mol_dynamic%cool_type = MD_STAIR
  case ('linear')
     mol_dynamic%cool_type = MD_LINEAR
  case ('log')
     mol_dynamic%cool_type = MD_LOG
  case default
     write(7,*)  
     write(7,*)'ERROR: unkown molecular dynamics cooling scheme: ', &
          trim(strflag)
     write(7,*)'STOP in usrinput'
     ierr = 154
     return
  end select

  dtmp = 500.0d0
  mol_dynamic%tempi = esdf_physical ('Tinit', dtmp,'k')
  dtmp = 500.0d0
  mol_dynamic%tempf = esdf_physical ('T_Final', dtmp,'k')
  dtmp = 150.0d0
  mol_dynamic%step_tmp = esdf_physical ('T_Step', dtmp,'k')

  dtmp = 150.0
  mol_dynamic%time_step =  esdf_physical ('Time_Step', dtmp,'ps')
  mol_dynamic%friction_coef = esdf_double ('Friction_Coefficient', &
       double_default)
  itmp = 250
  mol_dynamic%step_num = esdf_integer ('Step_num', itmp)
  ! For debug purposes - fix the rng seed
  itmp = 0
  mol_dynamic%rng_seed = esdf_integer('rng_seed', itmp)
  mol_dynamic%use_test_rng = esdf_boolean ('test_rng',boolean_default)

  if (.not.  mol_dynamic%is_on) then
     write(7,*) 'No molecular dynamics!'
  else
     call molecular_dynamic_turn_on (mol_dynamic,clust%atom_num)
     mol_dynamic%is_restart = esdf_boolean('Restart_mode',boolean_default)

     write(7,*) 'Molecular dynamics data:'
     write(7,*) '-------------------------'
     if (mol_dynamic%is_restart) then
        write(7,*) 'WARNING: restarted molecular dynamics run'
        write(7,*) 'restart data assumed present in mdinit.dat'
     endif
     if (mol_dynamic%rng_seed > 0) then
        write(7,*) 'WARNING: running with user defined RNG seed: ' , mol_dynamic%rng_seed
     endif
     if (mol_dynamic%use_test_rng ) then
        write(7,*) 'WARNING: running with simple RNG engine ' 
     endif

     select case (mol_dynamic%cool_type)
     case (MD_STAIR)
        write(7,*) '"Stairstep" cooling scheme '
        write(7,'(a,f6.1,a)') ' Temperature step = ',mol_dynamic%step_tmp, &
        ' [K]'
        i = (mol_dynamic%tempi-mol_dynamic%tempf) /mol_dynamic%step_tmp + 1
        mol_dynamic%stride = mol_dynamic%step_num/i
        if (mol_dynamic%stride * i < mol_dynamic%step_num) &
             mol_dynamic%stride = mol_dynamic%stride + 1
        write(7,*) 'Temperature modified after ',mol_dynamic%stride, &
             ' iterations'
     case (MD_LINEAR)
        write(7,*) 'Linear cooling scheme '
     case (MD_LOG)
        write(7,*) 'Logarithmic cooling scheme '
     end select

     write(7,52) ' Initial Temperature = ', mol_dynamic%tempi
     write(7,52) ' Final Temperature = ', mol_dynamic%tempf
     write(7,'(a,f8.2,a,f6.3,a)') ' Time step = ' ,mol_dynamic%time_step, &
          ' [au] = ',mol_dynamic%time_step*timeaups,' [ps]'
     write(7,'(a,f8.5,a)') ' Friction Coefficient = ', &
          mol_dynamic%friction_coef,' [au]'
      if (mol_dynamic%friction_coef < 10d-6 ) then
          write(7,*) 'Small coefficient, doing MICROcanonical ensemble'
          write(7,*) ' '
      else
          write(7,*) 'Large coefficient, doing (MACRO)canonical ensemble'
          write(7,*) ' '
      endif
     write(7,'(a,i5)') ' Number of MD iterations = ',mol_dynamic%step_num
     write(7,*) ' '
  endif
52 format(a,f6.1,' [K]')

  ! read in polarizability flag and electric field
  npolflg_tmp = esdf_boolean ('Polarizability', boolean_default)
  if (npolflg_tmp) then
     npolflg = 1
  else
     npolflg = 0
  endif
  dtmp = 0.001
  field = esdf_physical ('Electric_Field', dtmp,'ry/bohr/e')
  if (.not. npolflg_tmp) then
     npolflg = 0
     write(7,*) 'No polarizability!'
  else
     npolflg = 1
     write(7,*) 'Polarizability data:'
     write(7,*) '--------------------'
     write(7,*) 'Using seven-fields polarizability approach'
     write(7,'(a,f8.5,a)') ' Electric field  = ', field, ' [Ry/au]'
  endif
  write(7,*)
  !AJB: IMPORTANT
  ngroups = esdf_integer ('MPI_Groups', integer_default)
  !let's keep things simple, since the graph communicator is messed up
  ! if there are groups - RIGHT NOW
  !ngroups = esdf_integer ('MPI_Groups', 1)

  ! flag compatibility tests
  if((mol_dynamic%is_on) .or. (move%name /= NONE)) then
     if (npolflg_tmp) then
        write(7,*)  
        write(7,*) 'ERROR: Cannot do polarizability and move atoms'
        write(7,*) 'set MD and Minimization Flags = 0'
        write(7,*) 'STOP in usrinput'
        ierr = 155
        return
     endif
  endif

  if (pbc%per>2 .and. npolflg_tmp) then
     write(7,*)
     write(7,*) 'ERROR: Cannot do polarizability with periodic'
     write(7,*) 'boundary conditions!'
     write(7,*) 'STOP in usrinput'
     ierr = 156
     return
  endif

  if((mol_dynamic%is_on) .and. (move%name /= NONE)) then
     write(7,*)  
     write(7,*) 'ERROR: Cannot do minimization and MD'
     write(7,*) 'set either MD or minimization flag to zero'
     write(7,*) 'STOP in usrinput'
     ierr = 157
     return
  endif

  call esdf_close

  ! creating the 'empty' electronic_structure
  nstate = nstate 
  call create_electronic_struct (nspin,nstate,elec_st,elec_st%mxwd)
  !VDW:
  elec_st%do_vdw= do_vdw
  elec_st%hirsh_3d_cell = hirsh_cell_total
  elec_st%ncharge = ncharge 
  !
  elec_st%icorr = icorr
  if (icorr == 'cs') elec_st%dioniz = dioniz
  elec_st%tfermi = tfermi

  pot%nspin = nspin

  if (pbc%is_on) then
     call pbc_turn_on (natom,naty,pbc)
  endif

contains

        character(16) function build_gitcommit()
        implicit none
#ifdef GIT_COMMIT
        build_gitcommit = trim( GIT_COMMIT )
#else
        build_gitcommit = "Unknown"
#endif
        end function build_gitcommit

        character(1024) function build_config()
        implicit none
#ifdef MACH
        build_config = trim( MACH )
#else
        build_config = "Unknown"
#endif
        end function build_config

 
end subroutine usrinput
!===============================================================
