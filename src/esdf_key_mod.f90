!===============================================================
!
! Copyright (C) 2005 Finite Difference Research Group
! This file is part of parsec, http://www.ices.utexas.edu/parsec/
!
! Module to hold keyword list. this must be updated as
! new keywords are brought into existence.
!
! The 'label' is the label as used in calling the esdf routines
! 'typ' defines the type, with the following syntax. it is 3
! characters long.
! the first indicates:
!      i - integer
!      s - single
!      d - double
!      p - physical
!      t - string (text)
!      e - defined (exists)
!      l - boolean (logical)
!      b - block
! the second is always a colon (:)
! the third indicates the "level" of the keyword
!      b - basic
!      i - intermediate
!      e - expert
!      d - dummy
!
! 'Dscrpt' is a description of the variable. it should contain a
! (short) title enclosed between *! ... !*, and then a more detailed
! description of the variable.
!
!---------------------------------------------------------------
module esdf_key

  implicit none
  ! maximum number of kewords
  integer, parameter :: numkw = 200

  character (len=80)   :: kw_label(numkw)
  character (len=3)    :: kw_typ(numkw)
  character (len=3000) :: kw_dscrpt(numkw)

  integer :: kw_index(numkw)

  ! now define the keywords

  data kw_label(1)    / 'restart_run' /
  data kw_typ(1)      / 'L:B' /
  data kw_dscrpt(1)   / '*! Inital run flag !*' /
  data kw_label(2)    / 'periodic_system' /
  data kw_typ(2)      / 'L:B' /
  data kw_dscrpt(2)   / '*!System type (confined or periodic)!*'/
  data kw_label(3)    / 'boundary_sphere_radius' /
  data kw_typ(3)      / 'P:B' /
  data kw_dscrpt(3)   / '*! Radius of boundary sphere !*' /
  data kw_label(4)    / 'grid_spacing' /
  data kw_typ(4)      / 'P:E' /
  data kw_dscrpt(4)   / '*! Grid spacing (h) !*' /
  data kw_label(5)    / 'expansion_order' /
  data kw_typ(5)      / 'I:B' /
  data kw_dscrpt(5)   / '*!finite difference expansion order!*'/
  data kw_label(6)    / 'states_num' /
  data kw_typ(6)      / 'I:B' /
  data kw_dscrpt(6)   / '*! number of states !*' /
  data kw_label(7)    / 'net_charges' /
  data kw_typ(7)      / 'P:E' /
  data kw_dscrpt(7)   / '*! net charge !*' /
  data kw_label(8)    / 'fermi_temp' /
  data kw_typ(8)      / 'P:I' /
  data kw_dscrpt(8)   / '*! Fermi temp. !*' /
  data kw_label(9)    / 'max_iter' /
  data kw_typ(9)      / 'I:E' /
  data kw_dscrpt(9)   / '*! max iter. for self-consistency!*'/
  data kw_label(10)   / 'convergence_criterion' /
  data kw_typ(10)     / 'P:E' /
  data kw_dscrpt(10)  / '*! Convergence criterion !*' /
  data kw_label(11)   / 'diag_tolerance' /
  data kw_typ(11)     / 'D:E' /
  data kw_dscrpt(11)  / '*!  Diagonalization tolerance !*' /
  data kw_label(12)   / 'eigensolver' /
  data kw_typ(12)     / 'T:E' /
  data kw_dscrpt(12)  / '*! eigensolver - chebdav or arpack or else !*' /
  data kw_label(13)   / 'mixing_method' /
  data kw_typ(13)     / 'T:B' /
  data kw_dscrpt(13)  / '*! Mixer type (Anderson,Broyden,msecant1,msecant2,msecant3) !*' /
  data kw_label(14)   / 'mixing_param' /
  data kw_typ(14)     / 'D:B' /
  data kw_dscrpt(14)  / '*! Mixing parameter !*' /
  data kw_label(15)   / 'memory_param' /
  data kw_typ(15)     / 'I:B' /
  data kw_dscrpt(15)  / '*! Memory parameter (Broyden only) !*' /
  data kw_label(16)   / 'atom_types_num' /
  data kw_typ(16)     / 'I:E' /
  data kw_dscrpt(16)  / '*! Number of atom types !*' /
  data kw_label(17)   / 'atom_type' /
  data kw_typ(17)     / 'T:E' /
  data kw_dscrpt(17)  / '*! Type of the atom !*' /
  data kw_label(18)   / 'core_cutoff_radius' /
  data kw_typ(18)     / 'P:E' /
  data kw_dscrpt(18)  / '*!pseudopotential core cutoff radius!*'/
  data kw_label(19)   / 'potential_num' /
  data kw_typ(19)     / 'I:E' /
  data kw_dscrpt(19)  / '*! Number of potentials !*' /
  data kw_label(20)   / 'local_component' /
  data kw_typ(20)     / 'T:B' /
  data kw_dscrpt(20)  / '*! The local component (s,p or d) !*' /
  data kw_label(21)   / 'move_flag' /
  data kw_typ(21)     / 'I:I' /
  data kw_dscrpt(21)  / '*!Which atom to move(all,some,first n)!*'/
  data kw_label(22)   / 'lattice_vector_scale' /
  data kw_typ(22)     / 'P:D' /
  data kw_dscrpt(22)  / '*! Unit for lattice vectors !*' /
  data kw_label(23)   / 'correlation_type' /
  data kw_typ(23)     / 'T:D' /
  data kw_dscrpt(23)  / '*! Correlation type !*' /
  data kw_label(24)   / 'ion_energy_diff' /
  data kw_typ(24)     / 'P:I' /
  data kw_dscrpt(24)  / '*! Ion energy diff. (only for "cc") !*'/
  data kw_label(25)   / 'pseudopotential_format' /
  data kw_typ(25)     / 'T:I' /
  data kw_dscrpt(25)  / '*! Format of pseudopotential files !*'/
  data kw_label(26)   / 'minimization' /
  data kw_typ(26)     / 'T:E' /
  data kw_dscrpt(26)  / '*! Minimization type (none,simple,BFGS) !*'/
  data kw_label(27)   / 'movement_num' /
  data kw_typ(27)     / 'I:E' /
  data kw_dscrpt(27)  / '*! Number of movements !*' /
  data kw_label(28)   / 'force_min' /
  data kw_typ(28)     / 'P:I' /
  data kw_dscrpt(28)  / '*! Force minimum !*' /
  data kw_label(29)   / 'max_step' /
  data kw_typ(29)     / 'P:B' /
  data kw_dscrpt(29)  / '*! Max step!*' /
  data kw_label(30)   / 'min_step' /
  data kw_typ(30)     / 'P:E' /
  data kw_dscrpt(30)  / '*! Min step!*' /
  data kw_label(31)   / 'relax_restart' /
  data kw_typ(31)     / 'L:B' / 
  data kw_dscrpt(31)  / '*! Flag for relaxation restarting from previous run!*' /
  data kw_label(32)   / 'molecular_dynamics' /
  data kw_typ(32)     / 'L:B' /
  data kw_dscrpt(32)  / '*! Flag for calculating molecular dynamics !*' /
  data kw_label(33)   / 'restart_mode' /
  data kw_typ(33)     / 'L:I' /
  data kw_dscrpt(33)  / '*! Flag for restart mode !*' /
  data kw_label(34)   / 'cooling_method' /
  data kw_typ(34)     / 'T:I' /
  data kw_dscrpt(34)  / '*!Cooling type (log,linear,stair)!*'/
  data kw_label(35)   / 'tinit' /
  data kw_typ(35)     / 'P:I' / 
  data kw_dscrpt(35)  / '*! Tinit !*' /
  data kw_label(36)   / 't_final' /
  data kw_typ(36)     / 'P:E' / 
  data kw_dscrpt(36)  / '*! tfinal !*' /
  data kw_label(37)   / 't_step' /
  data kw_typ(37)     / 'P:I' / 
  data kw_dscrpt(37)  / '*! tstep (for "stairstep" method)!*' /
  data kw_label(38)   / 'time_step' /
  data kw_typ(38)     / 'P:I' /
  data kw_dscrpt(38)  / '*! Time step !*' /
  data kw_label(39)   / 'friction_coefficient' /
  data kw_typ(39)     / 'D:I' /
  data kw_dscrpt(39)  / '*! Friction coefficient !*' /
  data kw_label(40)   / 'step_num' /
  data kw_typ(40)     / 'I:E' /
  data kw_dscrpt(40)  / '*! Number of steps !*' /
  data kw_label(41)   / 'polarizability' /
  data kw_typ(41)     / 'L:B' /
  data kw_dscrpt(41)  / '*! Flag for polarizability !*' /
  data kw_label(42)   / 'electric_field' /
  data kw_typ(42)     / 'P:E' /
  data kw_dscrpt(42)  / '*! Applied electric field !*' /
  data kw_label(43)   / 'spin_polarization' /
  data kw_typ(43)     / 'L:E' /
  data kw_dscrpt(43)  / '*! Flag for spin polarization !*' /
  data kw_label(44)   / 'electron_per_orbital' /
  data kw_typ(44)     / 'B:E' /
  data kw_dscrpt(44)  / '*! # electrons per each orbital!*' /
  data kw_label(45)   / 'atom_coord' /
  data kw_typ(45)     / 'B:E' /
  data kw_dscrpt(45)  / '*! coordinates of atom!*'/
  data kw_label(46)   / 'coordinate_unit' /
  data kw_typ(46)     / 'T:D' /
  data kw_dscrpt(46)  / '*! Unit for atomic coordinates !*' /

  data kw_label(47)   / 'output_eigenvalues' /
  data kw_typ(47)     / 'L:B' /
  data kw_dscrpt(47)  / '*! Output eigenstate flag !*' /

  data kw_label(48)   / 'save_intermediate_eigenvalues' /
  data kw_typ(48)     / 'L:B' /
  data kw_dscrpt(48)  / '*! Save info on eigen.dat between SCF interations! Good for Debug Mode !*'/

  data kw_label(49)   / 'cell_shape' /
  data kw_typ(49)     / 'B:E' /
  data kw_dscrpt(49)  / '*! length of sides of periodic box !*' /

  data kw_label(50)   / 'output_all_states' /
  data kw_typ(50)     / 'L:B' /
  data kw_dscrpt(50)  / '*! output flag for parsec.dat !*' /

  data kw_label(51)    / 'boundary_conditions' /
  data kw_typ(51)      / 'T:B' /
  data kw_dscrpt(51)   / '*!Type of boundary conditions: cluster, tube, slab, crystal!*'/

  data kw_label(52)   / 'initial_spin_polarization' /
  data kw_typ(52)     / 'D:I' /
  data kw_dscrpt(52)  / '*!Initial spin polarization for each chemical element!*' /

  data kw_label(53)   / 'save_intermediate_charge_density' /
  data kw_typ(53)     / 'L:B' /
  data kw_dscrpt(53)  / '*!Save info on parsec.dat between SCF interations!*'/

  data kw_label(54)   / 'bfgs_number_corr' /
  data kw_typ(54)     / 'I:I' /
  data kw_dscrpt(54)  / '*!number of BFGS corrections!*' /

  data kw_label(55)   / 'restart_mixing' /
  data kw_typ(55)     / 'I:I' /
  data kw_dscrpt(55)  / '*!Number of iterations before mixing is restarted!*' /

  data kw_label(56)   / 'old_interpolation_format' /
  data kw_typ(56)     / 'L:B' /
  data kw_dscrpt(56)  / '*choose between interpolating Vps (old) and r*Vps (new) ' /

  data kw_label(57)   / 'output_level' /
  data kw_typ(57)     / 'I:I' /
  data kw_dscrpt(57)  / '*!level of output data in parsec.out!*' /

  data kw_label(58)   / 'maximum_matvec' /
  data kw_typ(58)     / 'I:I' /
  data kw_dscrpt(58)  / '*!Maximum number of matrix-vector operations allowed!*' /

  data kw_label(59)   / 'subspace_buffer_size' /
  data kw_typ(59)     / 'I:I' /
  data kw_dscrpt(59)  / '*! Number of additional states for each representation !*' /

  data kw_label(60)   / 'van_der_waals' /
  data kw_typ(60)     / 'L:B' /
  data kw_dscrpt(60)  / '*! Calculate the vdw contribution to forces !*' /

  data kw_label(61)   / 'alpha_filter' /
  data kw_typ(61)     / 'D:I' /
  data kw_dscrpt(61)  / '*!alpha parameter in Fourier filtering!*' /

  data kw_label(62)   / 'beta1_filter' /
  data kw_typ(62)     / 'D:I' /
  data kw_dscrpt(62)  / '*!beta_1 parameter in Fourier filtering!*' /

  data kw_label(63)   / 'core_filter' /
  data kw_typ(63)     / 'D:I' /
  data kw_dscrpt(63)  / '*!alpha parameter in core charge Fourier filtering!*' /

  data kw_label(64)   / 'chebyshev_degree' /
  data kw_typ(64)     / 'I:E' /
  data kw_dscrpt(64)  / '*! Chebyshev polynomial degree, usually within 7 to 25 !*' /

  data kw_label(65)   / 'save_wfn_bands' /
  data kw_typ(65)     / 'B:E' /
  data kw_dscrpt(65)  / '*! List of selected bands/levels to be saved!*' /

  data kw_label(66)   / 'mpi_groups' /
  data kw_typ(66)     / 'I:E' /
  data kw_dscrpt(66)  / '*! number of groups of PEs!*'/

  data kw_label(67)   / 'fix_number_eigenvalues' /
  data kw_typ(67)     / 'L:E' /
  data kw_dscrpt(67)  / '*! Skip dynamical adjustment of neig !*' /

  data kw_label(68)   / 'symmetrize_charge_density' /
  data kw_typ(68)     / 'L:E' /
  data kw_dscrpt(68)  / '*! Symmetrize charge density during SCF!*' /

  data kw_label(69)   / 'chebyshev_filtering' /
  data kw_typ(69)     / 'L:E' /
  data kw_dscrpt(69)  / '*! Flag for Chebyshev subspace filtering !*' /

  data kw_label(70)   / 'chebdav_degree' /
  data kw_typ(70)     / 'I:E' /
  data kw_dscrpt(70)  / '*! Chebyshev-Davidson degree, usually 15 or more !*' /

  data kw_label(71)   / 'add_point_charges' /
  data kw_typ(71)     / 'L:B' /
  data kw_dscrpt(71)  / '*! Point charges flag !*' /

  data kw_label(72)   / 'point_typ_num' /
  data kw_typ(72)     / 'I:E' /
  data kw_dscrpt(72)  / '*! Number of point charge types !*' /

  data kw_label(73)   / 'pt_chg' /
  data kw_typ(73)     / 'P:I' /
  data kw_dscrpt(73)  / '*! charge of point charge of given type !*' /

  data kw_label(74)   / 'point_coord' /
  data kw_typ(74)     / 'B:E' /
  data kw_dscrpt(74)  / '*! coordinates of point charges!*'/

  data kw_label(76)   / 'kpoint_method' /
  data kw_typ(76)     / 'T:D' /
  data kw_dscrpt(76)  / '*! method to determine kpts !*' /

  data kw_label(77)   / 'kpoint_list' /
  data kw_typ(77)     / 'B:E' /
  data kw_dscrpt(77)  / '*! kpoints entered manually!*'/

  data kw_label(78)   / 'complex_wfn' /
  data kw_typ(78)     / 'L:B' /
  data kw_dscrpt(78)/ '*! Flag to identify wfn cmplx or real!*' / 

  data kw_label(79)   / 'monkhorst_pack_grid' /
  data kw_typ(79)     / 'B:E' /
  data kw_dscrpt(79)  / '*! MP grid !*' /

  data kw_label(80)   / 'monkhorst_pack_shift' /
  data kw_typ(80)     / 'B:E' /
  data kw_dscrpt(80)  / '*! shift in MP grid !*' /

  data kw_label(81)   / 'restart_from_wfndat' /
  data kw_typ(81)     / 'L:E' /
  data kw_dscrpt(81)  / '*! Restart from file wfn.dat !*' /

  data kw_label(82)   / 'output_file_name' /
  data kw_typ(82)     / 'T:E' /
  data kw_dscrpt(82)  / '*! Name of major output file !*' /

  data kw_label(83)   / 'double_grid_order' /
  data kw_typ(83)     / 'I:E' /
  data kw_dscrpt(83)  / '*! order of double grid in pseudopotentials!*'/

  data kw_label(84)   / 'ldaplusu_u' /
  data kw_typ(84)     / 'P:E' /
  data kw_dscrpt(84)  / 'The Hubbard U parameter for LDA+U' /

  data kw_label(85)   / 'ldaplusu_j' /
  data kw_typ(85)     / 'P:E' /
  data kw_dscrpt(85)  / 'The Hubbard J parameter for LDA+U' /

  data kw_label(86)   / 'nonlinear_core_correction' /
  data kw_typ(86)     / 'L:E' /
  data kw_dscrpt(86)  / 'Non-linear core correction ' /

  data kw_label(87)   / 'ignore_symmetry' /
  data kw_typ(87)     / 'L:B' /
  data kw_dscrpt(87)  / '*! turn all symmetry operations off !*' /

  data kw_label(88)    / 'so_psp' /
  data kw_typ(88)      / 'L:B' /
  data kw_dscrpt(88)   / '*! Apply SO correction to specific atom!*'/

  data kw_label(89)    / 'scf_so' /
  data kw_typ(89)      / 'L:B' /
  data kw_dscrpt(89)   / '*!Apply SO correction self-consistently!*'/

  data kw_label(90)    / 'dynamic_diag_tol' /
  data kw_typ(90)      / 'L:B' /
  data kw_dscrpt(90)   /'*!Apply dynamic diagonalization tolerance!*'/

  data kw_label(91)   / 'kpoint_unit' /
  data kw_typ(91)     / 'T:D' /
  data kw_dscrpt(91)  / '*! Unit for coordinates of k-points !*' /

  data kw_label(92)   / 'chebyshev_degree_delta' /
  data kw_typ(92)     / 'I:E' /
  data kw_dscrpt(92)  / '*! Chebyshev polynomial degree (delta), usually within 1 to 3 !*' /

  data kw_label(93)   / 'convergence_criterion_approach' /
  data kw_typ(93)     / 'P:E' /
  data kw_dscrpt(93)  / '*! Indicator of convergence approach, usually 0.01 !*' /

  data kw_label(94)   / 'mixing_group_size' /
  data kw_typ(94)     / 'I:E' /
  data kw_dscrpt(94)  / '*! Size of groups of secant equations !*' /

  data kw_label(95)   / 'mixing_en_like' /
  data kw_typ(95)     / 'I:E' /
  data kw_dscrpt(95)  / '*! Flag to enable the Broyden-like update !*' /

  data kw_label(96)   / 'mixing_preferred_type' /
  data kw_typ(96)     / 'I:E' /
  data kw_dscrpt(96)  / '*! Preferred type of hybrid method !*' /

  data kw_label(97)   / 'mixing_restart_factor' /
  data kw_typ(97)     / 'D:E' /
  data kw_dscrpt(97)  / '*! Restart in multisecant mixing methods !*' /

  data kw_label(98)   / 'mixing_memory_expand_factor' /
  data kw_typ(98)     / 'D:E' /
  data kw_dscrpt(98)  / '*! Factor used to expand memory for mixing !*' /

  data kw_label(99)   / 'super_cell' /
  data kw_typ(99)     / 'B:E' /
  data kw_dscrpt(99)  / '*! Duplicate of the super-cell size !*' /

  data kw_label(100)   / 'super_cell_vac' /
  data kw_typ(100)     / 'B:E' /
  data kw_dscrpt(100)  / '*! Amount of vacuum to be put in the system !*' /

  data kw_label(101)   / 'external_mag_field'/
  data kw_typ(101)     / 'D:E' /
  data kw_dscrpt(101)  / '*! External magnetic field in Gauss !*' /

  data kw_label(102)   / 'solver_lpole' /
  data kw_typ(102)     / 'I:E' /
  data kw_dscrpt(102)  / '*! lpole order for the multipole expansion !*' /

  data kw_label(103)   / 'cluster_domain_shape' /
  data kw_typ(103)     / 'T:E' /
  data kw_dscrpt(103)  / '*! choice of domain shape for cluster calculations !*' /

  data kw_label(104)   / 'domain_shape_parameters' /
  data kw_typ(104)     / 'B:E' /
  data kw_dscrpt(104)  / '*! parameters for the cluster calculation domain !*' /

  data kw_label(105)   / 'bandstruc'/
  data kw_typ(105)     / 'B:E' /
  data kw_dscrpt(105)  / '*! lines in kspace for band calculation  !*' /

  data kw_label(106)   / 'bandstruc_points'/
  data kw_typ(106)     / 'I:E' /
  data kw_dscrpt(106)  / '*! num of pts for shortest line in band structure !*' /

  data kw_label(107)   / 'rng_seed'/
  data kw_typ(107)     / 'I:E' /
  data kw_dscrpt(107)  / '*! Override randomly generated seed for MD !*' /

  data kw_label(108)   / 'test_rng'/
  data kw_typ(108)     / 'L:E' /
  data kw_dscrpt(108)  / '*! Use simple RNG for MD (debug only!) !*' /


  data kw_label(110)    / 'create_dos' /
  data kw_typ(110)      / 'L:B' /
  data kw_dscrpt(110)   /'Whether to create a Density of States file'/


  data kw_label(111)   / 'dos_mp_grid' /
  data kw_typ(111)     / 'B:E' /
  data kw_dscrpt(111)  / '*! MP grid for creating DOS !*' /

  data kw_label(112)   / 'dos_mp_shift' /
  data kw_typ(112)     / 'B:E' /
  data kw_dscrpt(112)  / '*! shift in MP grid for creating DOS !*' /

  data kw_label(113)    / 'ylm_dos_l' /
  data kw_typ(113)      / 'I:E' /
  data kw_dscrpt(113)   /'The l of the Ylm used for Angular DOS (l<=ylm_dos)'/

  data kw_label(114)    / 'ylm_dos_cutoff' /
  data kw_typ(114)      / 'P:E' /
  data kw_dscrpt(114)   /'Cuttoff for Ylm caluclation in atomic unit'/  

  data kw_label(115)    / 'dos_pnum' /
  data kw_typ(115)      / 'I:E' /
  data kw_dscrpt(115)   /'Number of points in DOS graph'/  

  data kw_label(116)    / 'orb_mag' /
  data kw_typ(116)      / 'L:B' /
  data kw_dscrpt(116)   / '*!Apply orbital magnetism!*'/

  data kw_label(117)    / 'Non_Collinear_magnetism' /
  data kw_typ(117)      / 'L:B' /
  data kw_dscrpt(117)   / '*!Apply non-collinear magnetism!*'/

  data kw_label(118)    / 'SO_from_scratch' /
  data kw_typ(118)      / 'L:B' /
  data kw_dscrpt(118)   / '*!Use spinor representation from scratch!*'/

  data kw_label(119)   / 'Initial_NCL_Moment' /
  data kw_typ(119)     / 'B:E' /
  data kw_dscrpt(119)  / '*! initial 3D magnetic moment of atom!*'/

  data kw_label(120)   / 'Atomic_Radius' /
  data kw_typ(120)     / 'D:E' /
  data kw_dscrpt(120)  / '*!Radius for each type of atoms. Default in rtable.f90!*'/

  data kw_label(121)   / 'cubic_spline' /
  data kw_typ(121)     / 'L:I' /
  data kw_dscrpt(121)  / '*!Use cubic spline interpolation for pseudopotentials!*'/

  data kw_label(122)   / 'Read_VCD' /
  data kw_typ(122)     / 'L:I' /
  data kw_dscrpt(122)  / '*!Read the vallence charge density from file, for ionic semicore PPs!*'/

  data kw_label(123)   / 'Full_Hartree' /
  data kw_typ(123)     / 'L:I' /
  data kw_dscrpt(123)  / '*!use full hartree boundary conditions!*'/

  data kw_label(124)   / 'matvec_blocksize' /
  data kw_typ(124)     / 'I:E' /
  data kw_dscrpt(124)  / '*! Block size for everything matvec related !*' /

  data kw_label(125)   / 'output_gw' /
  data kw_typ(125)     / 'L:B' /
  data kw_dscrpt(125)  / '*!Whether we have GW output for BerkeleyGW!*' /

  data kw_label(126)   / 'vxc_matrix_elements' /
  data kw_typ(126)     / 'T:E' /
  data kw_dscrpt(126)  / '*!vxc matrix elements needed bof BerkeleyGW!*' /

  data kw_label(127)   / 'nscf_kpoints' /
  data kw_typ(127)     / 'B:E' /
  data kw_dscrpt(127)  / '*!K-points on which we perform the non-self consistent calculation!*' /

  data kw_label(128)   / 'nscf_kgrid' /
  data kw_typ(128)     / 'B:E' /
  data kw_dscrpt(128)  / '*! MP grid for nscf kpoints (needed for BerkeleyGW)!*' /

  data kw_label(129)   / 'nscf_kgrid_shift' /
  data kw_typ(129)     / 'B:E' /
  data kw_dscrpt(129)  / '*! shift in MP grid for nscf kpoints (needed for BerkeleyGW) !*' /

  data kw_label(130)   / 'nscf_kpoint_unit' /
  data kw_typ(130)     / 'T:D' /
  data kw_dscrpt(130)  / '*! Unit for coordinates of nscf k-points !*' /

  data kw_label(131)    / 'nscf_states_num' /
  data kw_typ(131)      / 'I:B' /
  data kw_dscrpt(131)   / '*! number of states for nscf calculation !*' /

  data kw_label(132)    / 'nscf_fermi_level' /
  data kw_typ(132)      / 'P:B' /
  data kw_dscrpt(132)   / '*! Fermi level in Ry for nscf calculation !*' /

  data kw_label(133)    / 'add_charged_sheets' /
  data kw_typ(133)      / 'L:E' /
  data kw_dscrpt(133)   / '*! Flag for addition of one (or more) fixed sheet(s) of charge!*' /

  data kw_label(134)    / 'num_charged_sheets' /
  data kw_typ(134)      / 'I:E' /
  data kw_dscrpt(134)   / '*! Number of charged sheets to be inserted !*' /
  
  data kw_label(135)   / 'sheet_charge' /
  data kw_typ(135)     / 'P:E' /
  data kw_dscrpt(135)  / '*! charge of sheet charge of given type !*' /

  data kw_label(136)   / 'sheet_coord' /
  data kw_typ(136)     / 'B:E' /
  data kw_dscrpt(136)  / '*! coordinates of sheet of given type !*'/

  data kw_label(137)   / 'Exp_Comm' /
  data kw_typ(137)     / 'L:E' /
  data kw_dscrpt(137)  / '*! Experimental flag to allow new hidden communication mode !*'/

  data kw_label(138)   / 'output_grid_data' /
  data kw_typ(138)     / 'T:I' /
  data kw_dscrpt(138)  / '*! Output grid data to external file. &
        &Options: charge_density, v_external, v_hartree, v_external+v_hartree (sum) !*' /

  data kw_label(139)    / 'ff_maxiter' /
  data kw_typ(139)      / 'I:E' /
  data kw_dscrpt(139)   / '*! number of max iterations for the first-filter eigen-solver !*' /

  data kw_label(140)    / 'max_lap_buffers' /
  data kw_typ(140)      / 'I:E' /
  data kw_dscrpt(140)   / '*! maximum number of buffers to store comm neighbors on during matvec !*' /

  data kw_label(141)   / 'hf_one' /
  data kw_typ(141)     / 'I:B' /
  data kw_dscrpt(141)  / '*! (not Hartree-Fock) 2*hf_d+1 is the number of cells for hf in 3D per!*' /

  data kw_label(142)   / 'hf_two' /
  data kw_typ(142)     / 'I:B' /
  data kw_dscrpt(142)  / '*! (not Hartree-Fock)  2*hf_d+1 is the number of cells for hf in 3D per!*' /

  data kw_label(143)   / 'hf_three' /
  data kw_typ(143)     / 'I:B' /
  data kw_dscrpt(143)  / '*! (not Hartree-Fock)  2*hf_d+1 is the number of cells for hf in 3D per!*' /

  data kw_label(144)   / 'Use_Plain_SRE' /
  data kw_typ(144)     / 'L:E' /
  data kw_dscrpt(144)  / '*! do not weight the SRE by the grid charge - try if there is lots of vacuum!*'/
  
  data kw_label(145)   / 'Calc_Hartree_G' /
  data kw_typ(145)     / 'L:E' /
  data kw_dscrpt(145)  / '*! For 3D PBC, can calculate the hartree potential in G space!*'/

  data kw_label(146)   / 'Output_Data' /
  data kw_typ(146)     / 'L:E' /
  data kw_dscrpt(146)  / '*! Whether to write raw data at all or not!*'/
end module esdf_key
!===============================================================
