! ===============================================================
! Copyright (C) 2005 Finite Difference Research Group
! This file is part of parsec, http://www.ices.utexas.edu/parsec/
!
! author: A. Makmal (2005).
!
! ==================== Cluster structure ========================

module cluster_module
  use constants
  implicit none

  type cluster    

     real(dp), dimension (:), pointer :: xatm
     real(dp), dimension (:), pointer :: yatm
     real(dp), dimension (:), pointer :: zatm
     integer, dimension (:), pointer :: mvat 

     ! amass - atomic mass of each ion
     real(dp),  dimension (:), pointer :: amass

     ! number of atoms for each atom type
     integer, dimension (:), pointer :: natmi

     ! force in each axis for each atom
     real(dp), dimension (:,:), pointer :: force

     ! vdw force for each atom in each axis
     real(dp), dimension (:,:), pointer :: vdw_forces

     ! distance of atoms from grid points
     real(dp), dimension (:,:), pointer :: dist_atm

     ! contribution of spin-orbit to non-local force
     complex(dpc), dimension (:,:), pointer :: force_so

     character(len=2),  dimension (:), pointer :: name

     ! total number of atom types
     integer :: type_num
     ! total number of atoms
     integer :: atom_num

     ! flag for the existence of point charges
     logical :: has_ptchrg
     ! number of point charge types
     integer :: npttyp
     ! charge of each point charge type
     real(dp), dimension (:), pointer :: qpt
     ! number of point charges of each type
     integer, dimension (:), pointer :: nptt
     ! total number of point charges
     integer :: nptchrg
     ! type for each atom  !AJB: this should move up!
     integer, dimension (:), pointer :: atype
     ! coordinates of point charges
     real(dp), dimension (:), pointer :: xpt
     real(dp), dimension (:), pointer :: ypt
     real(dp), dimension (:), pointer :: zpt

     ! flag for the existence of charged sheets
     logical :: has_charged_sheet
     ! number of charged sheet types
     integer :: number_charged_sheet_type
     ! number of charged sheets 
     integer :: number_charged_sheets
     ! charge of each sheet charge type
     real(dp), dimension (:), pointer :: sheet_charge_of_type
     ! number of charged sheets of each type
     integer, dimension (:), pointer :: number_charged_sheets_per_type
     ! coordinates of point charges
     ! real(dp), dimension (:), pointer :: x_charged_sheets
     ! real(dp), dimension (:), pointer :: y_charged_sheets
     real(dp), dimension (:), pointer :: z_charged_sheets

     ! the cutoff radius used to calculate the angular dos
     ! specified for each type
     real(dp), dimension (:), pointer :: ylmdos_cutoff

  end type cluster

contains

  subroutine init_cluster(clust)
    implicit none
    type (cluster), intent (inout) :: clust

    nullify(clust%xatm , clust%yatm, clust%zatm, clust%mvat, clust%amass, clust%natmi, clust%force)
    nullify(clust%force_so, clust%name, clust%ylmdos_cutoff, clust%xpt, clust%ypt)
    nullify(clust%zpt, clust%qpt, clust%nptt, clust%atype)
    nullify(clust%number_charged_sheets_per_type,clust%sheet_charge_of_type,clust%z_charged_sheets)

  end subroutine init_cluster

  subroutine destroy_cluster (clust)
    implicit none
    type (cluster), intent (inout) :: clust

    if (associated (clust%xatm)) deallocate (clust%xatm) 
    if (associated (clust%yatm)) deallocate (clust%yatm) 
    if (associated (clust%zatm)) deallocate (clust%zatm)      
    if (associated (clust%mvat)) deallocate (clust%mvat)    
    if (associated (clust%amass)) deallocate (clust%amass)  
    if (associated (clust%natmi)) deallocate (clust%natmi)
    if (associated (clust%force)) deallocate (clust%force)
    if (associated (clust%vdw_forces)) deallocate (clust%vdw_forces)
    if (associated (clust%force_so)) deallocate (clust%force_so)
    if (associated (clust%name)) deallocate (clust%name)
    if (associated (clust%ylmdos_cutoff)) deallocate (clust%ylmdos_cutoff)      

    if (associated (clust%xpt)) deallocate (clust%xpt)
    if (associated (clust%ypt)) deallocate (clust%ypt)
    if (associated (clust%zpt)) deallocate (clust%zpt)
    if (associated (clust%qpt)) deallocate (clust%qpt)
    if (associated (clust%nptt)) deallocate (clust%nptt)
    if (associated (clust%atype)) deallocate (clust%atype)

    if (associated (clust%dist_atm)) deallocate(clust%dist_atm)
    if (associated (clust%sheet_charge_of_type)) deallocate (clust%sheet_charge_of_type)
    if (associated (clust%number_charged_sheets_per_type)) deallocate (clust%number_charged_sheets_per_type)
    if (associated (clust%z_charged_sheets)) deallocate (clust%z_charged_sheets)

  end subroutine destroy_cluster

  subroutine create_cluster (atom_num, clust)
    implicit none
    integer, intent (inout) :: atom_num
    type (cluster), intent (inout) :: clust

    clust%atom_num = atom_num
    allocate (clust%mvat (atom_num))
    allocate (clust%zatm (atom_num))
    allocate (clust%yatm (atom_num))
    allocate (clust%xatm (atom_num))
    allocate (clust%amass (atom_num))
    allocate (clust%atype (atom_num))
    allocate (clust%force (3,atom_num))
! note that this array is transposed 10110!
    allocate (clust%vdw_forces (atom_num,3))
    allocate (clust%force_so (3,atom_num))
    clust%force_so(:,:) = zzero
    clust%vdw_forces(:,:) = zero

  end subroutine create_cluster

end module cluster_module

! ==================== electronic_struct structure ==============
module electronic_struct_module
  use constants
  implicit none

  type eigenstates

     ! label of the PE group in charge of this representation
     integer :: group
     ! number of computed eigenstates for this representation
     integer :: mm
     ! target number of eigenstates (may differ from mm after a call
     ! to eigval.F because of dynamical adjustment)
     integer :: nn
     ! number of converged eigenstates for this representation
     integer :: nec
     ! eigenvalues
     real(dp), dimension (:), pointer :: en
     ! occupation of each state
     real(dp), dimension (:), pointer :: occ
     ! wave function array
     real(dp), dimension (:,:), pointer :: wf
     complex(dpc), dimension (:,:), pointer :: zwf

  end type eigenstates

  type electronic_struct

     ! number of kpoints (<=0 means gamma only)
     integer :: nkpt
     ! kpoint selection method
     integer :: kptmethod
     ! Monkhorts-Pack grid
     integer :: mpgrid(3)
     ! Monkhorts-Pack shift
     real(dp) :: mpshift(3)
     ! DOS Monkhorts-Pack grid
     integer :: dos_mpgrid(3)
     ! DOS Monkhorts-Pack shift
     real(dp) :: dos_mpshift(3)
            
                    
     ! kpoint list (Cartesian coordinates, units of inverse Bohr radius)
     real(dp), dimension (:,:), pointer :: kpts
     ! kpoint weights 
     real(dp), dimension (:), pointer :: kpwt
     ! number of kpoints in full Brillouin zone
     integer :: nkptf
     ! kpoints in full Brillouin zone
     real(dp), dimension (:,:), pointer :: kptf
     ! index of kpoint in IBZ equivalent to each kpoint in full BZ
     integer, dimension (:), pointer :: kptindr
     ! symmetry operation which brings a kpoint in full BZ to the IBZ
     integer, dimension (:), pointer :: kptitran
     ! index of kpoint in IBZ for every kpoint in full BZ, ordered in the MP grid
     integer, dimension (:,:,:), pointer :: kptgrid

     ! number of spins (1 or 2)
     integer :: nspin
     ! number of requested eigenstates
     integer :: nstate
     ! total number of computed eigenstates for each spin orientation,
     ! after droping all states above the maximum computed for any repr.
     integer, dimension(:), pointer :: ntotal
     ! net charge in the system (in units of electron charge)
     real(dp) :: ncharge
     ! actual number of valence electrons in structure
     real(dp) :: xele
     ! electron temperature (in Kelvin)
     real(dp) :: tfermi
     ! fermi level
     real(dp) :: efermi
     ! LDA+U energy
     real(dp) :: etot_plusu
     ! total energy in Rydberg, total energy per atom in eV
     real(dp) :: etot, bdev
     ! exchange-correlation energy
     real(dp) :: exc
     ! ionization energy difference (used with cc correlation flag only)
     real(dp) :: dioniz
     ! choice of exchange-correlation functional (see exc_spn.f and
     ! exc_nspn.f for possible choices)
     character (len=2) :: icorr
     ! flag for real/complex eigenvectors:
     ! .true. if eigenvectors are complex
     ! .false. otherwise
     logical :: cplx

     ! dipole along the three axes
     real(dp) :: dipx, dipy, dipz
     ! charge of cluster
     real(dp) :: charge

     ! total number of electrons in each spin channel
     real(dp), dimension (:), pointer :: totel

     ! number of occupied states in a given iteration
     ! (including partially occupied states and all representations)
     integer, dimension (:), pointer :: ifmax
     ! user-provided occupation of each state (including all
     ! representations)
     real(dp), dimension (:,:), pointer :: occ_in
     ! charge density, global and distributed
     real(dp), dimension (:,:), pointer :: rho
     ! core charge density constructed from superimposing the atomic
     ! core-correction charge densities
     real(dp), dimension (:), pointer :: rhoc
!!! VDW 
     logical :: do_vdw
     ! vdw energy contribution to total energy
     real(dp) :: evdw
     ! charge density for the free atoms
     real(dp), dimension (:,:,:), pointer :: rhof
     ! free_volume for periodic structures
     real(dp), dimension (:,:), pointer :: v_free_periodic
     ! hirshfeld function for periodic structures
     real(dp), dimension (:,:,:), pointer :: hirshfeld_periodic
     ! r for hirshfeld function periodic
     real(dp), dimension (:,:,:), pointer :: dist_periodic
     ! 2*hirsh_3d_cell+1 is the number of cells in each direction for hirshfeld function
     integer, dimension(:,:), pointer :: hirsh_3d_cell
!!! END VDW PART
     ! charge-weighted SRE
     real(dp), dimension (:), pointer :: sre
     ! non wighted SRE (just (\delta v)^2)
     real(dp), dimension (:), pointer :: plain_sre
     ! which sre to actually use?
     logical :: use_plain_sre

     ! number of symmetry representations
     integer :: nrep

     ! representation of each state, after eigenvalues are arranged in
     ! ascending order
     integer, dimension(:,:,:), pointer :: irep

     ! table of characters for the reduced symmetry group
     integer, dimension(:,:), pointer :: chi

     ! set of eigenvectors to be saved to parsec.dat
     integer :: nsave
     integer, dimension(:), pointer :: indxsave

     ! eigenvalues/eigenvectors for each symmetry representation and spin
     type(eigenstates), dimension(:,:,:), pointer :: eig

     ! flag for perturbative spin-orbit correction
     logical :: is_so
     logical :: scf_so
     logical :: so_pert
     ! flag for non-collinear calculation
     logical :: ncl
     !  non-collinear contribution to total energy
     real(dp) :: totnrg_ncl
     ! initial atomic 3D magnetic moment
     real(dp), dimension (:,:), pointer :: init_mag
     ! integrated 3d magnetization if non-collinear calculation
     real(dp) :: s3dtot(3)
     ! 3d magnetization if non-collinear calculation
     real(dp), dimension (:,:), pointer ::  spin3d
     ! norm of 3d magnetization if non-collinear calculation
     real(dp), dimension (:), pointer ::  spin3dn
     ! additional spin dimension of wave-function due to spin-orbit
     ! and\or non-collinear states
     integer :: mxwd
     integer :: ndim
     ! In case of external magnetic field, it's magnitude
     real(dp) :: mag_field
     ! Flag for external magnetic field
     logical :: is_mag
     ! In case of spin-orbit: Net magnetic moment of all states
     real(dp) :: net_magmom
     ! In case of spin-orbit: magnetic moment of each state
     real(dp), dimension (:,:), pointer ::  magmom
     ! flag for orbital magnetism
     logical :: orb_mag
     ! flag for current dft
     logical :: is_cur
     ! flag for performing an explicit evaluation of the hartree potential for pbc%per=3
     logical :: explicit_hartree_pbc
     ! current density array
     real(dp), dimension (:,:,:), pointer :: jop
     ! exchange current density array
     real(dp), dimension (:,:,:), pointer :: jex

  end type  electronic_struct

contains

  subroutine init_electronic_struct(elec_st)
    implicit none
    type (electronic_struct), intent (inout) :: elec_st

    nullify(elec_st%ntotal)
    nullify(elec_st%totel)
    nullify(elec_st%ifmax)
    nullify(elec_st%occ_in)
    nullify(elec_st%rho)
    nullify(elec_st%rhoc)
    nullify(elec_st%sre)
    nullify(elec_st%plain_sre)
    nullify(elec_st%irep)
    nullify(elec_st%chi)
    nullify(elec_st%kpts)
    nullify(elec_st%kpwt)
    nullify(elec_st%kptf)
    nullify(elec_st%kptindr)
    nullify(elec_st%kptitran)
    nullify(elec_st%indxsave)
    nullify(elec_st%magmom)
    nullify(elec_st%spin3d)
    nullify(elec_st%spin3dn)
    nullify(elec_st%jop)
    nullify(elec_st%jex)
    nullify(elec_st%eig)

    elec_st%mxwd = -1
    elec_st%ndim = -1
    elec_st%nrep = -1
    elec_st%use_plain_sre = .FALSE.
    elec_st%explicit_hartree_pbc = .FALSE.

  end subroutine init_electronic_struct

  subroutine destroy_electronic_struct (elec_st)
    implicit none
    type (electronic_struct), intent (inout) :: elec_st
    integer :: irp,isp,kplp

    if (associated (elec_st%ntotal)) deallocate (elec_st%ntotal)
    if (associated (elec_st%totel)) deallocate (elec_st%totel)
    if (associated (elec_st%ifmax)) deallocate (elec_st%ifmax)
    if (associated(elec_st%occ_in)) deallocate(elec_st%occ_in)
    if (associated (elec_st%rho)) deallocate (elec_st%rho)
    if (associated (elec_st%rhoc)) deallocate (elec_st%rhoc)
!!! VDW PART
    if (associated (elec_st%rhof)) deallocate (elec_st%rhof)
    if (associated (elec_st%v_free_periodic)) deallocate (elec_st%v_free_periodic)
    if (associated (elec_st%hirshfeld_periodic)) deallocate (elec_st%hirshfeld_periodic)
    if (associated (elec_st%dist_periodic)) deallocate (elec_st%dist_periodic)
    if (associated (elec_st%hirsh_3d_cell)) deallocate (elec_st%hirsh_3d_cell)
!!! END VDW PART
    if (associated (elec_st%sre)) deallocate (elec_st%sre)
    if (associated (elec_st%plain_sre)) deallocate (elec_st%plain_sre)
    if (associated (elec_st%irep)) deallocate (elec_st%irep)
    if (associated (elec_st%chi)) deallocate (elec_st%chi)
    if (associated (elec_st%kpts)) deallocate (elec_st%kpts)
    if (associated (elec_st%kpwt)) deallocate (elec_st%kpwt)
    if (associated (elec_st%kptf)) deallocate (elec_st%kptf)
    if (associated (elec_st%kptindr)) deallocate (elec_st%kptindr)
    if (associated (elec_st%kptitran)) deallocate (elec_st%kptitran)
    if (associated (elec_st%indxsave)) deallocate(elec_st%indxsave)
    if (associated (elec_st%magmom)) deallocate(elec_st%magmom)
    if (associated (elec_st%spin3d)) deallocate(elec_st%spin3d)
    if (associated (elec_st%spin3dn)) deallocate(elec_st%spin3dn)
    if (associated(elec_st%jop)) deallocate(elec_st%jop)
    if (associated(elec_st%jex)) deallocate(elec_st%jex)
      
    if (associated(elec_st%eig)) then
        do isp = 1, elec_st%nspin/elec_st%mxwd
           do kplp =1,max(elec_st%nkpt,1)
          do irp = 1, elec_st%nrep
             call destroy_eigenstate(elec_st%eig(irp,kplp,isp))
          enddo
           enddo
        enddo
        deallocate(elec_st%eig)
    endif

  end subroutine destroy_electronic_struct
  subroutine destroy_eigenstate (eigen)
    implicit none
    type (eigenstates), intent (inout) :: eigen

    if (associated(eigen%en)) deallocate(eigen%en)
    if (associated(eigen%occ)) deallocate(eigen%occ)
    if (associated(eigen%wf)) deallocate(eigen%wf)
    if (associated(eigen%zwf)) deallocate(eigen%zwf)

  end subroutine destroy_eigenstate

  subroutine create_eigenstate (eigen,ldn,nn,cplx)
    implicit none
    type (eigenstates), intent (inout) :: eigen
    integer, intent(in) :: ldn, nn
    logical, intent(in) :: cplx
    integer :: alcstat

    allocate(eigen%en(nn))
    allocate(eigen%occ(nn))
    eigen%en(:) = zero
    eigen%occ(:) = zero
    if (cplx) then
       allocate(eigen%zwf(ldn,nn),stat=alcstat)
       nullify(eigen%wf) !bugfix
       call alccheck('eigen%zwf',ldn*nn,alcstat)
       eigen%zwf(:,:) = zzero
    else
       allocate(eigen%wf(ldn,nn),stat=alcstat)
       nullify(eigen%zwf) !bugfix
       call alccheck('eigen%wf',ldn*nn,alcstat)
       eigen%wf(:,:) = zero
    endif

  end subroutine create_eigenstate

  subroutine create_electronic_struct(nspin,nstate,elec_st,mxwd)
    implicit none
    integer, intent (inout) :: nspin
    integer, intent (inout) :: nstate
    integer, intent (in) :: mxwd
    
    ! default values
    integer, parameter :: int_def = -1
    real(dp), parameter :: real_def = zero
    type (electronic_struct), intent (inout) :: elec_st  

    elec_st%nspin = nspin
    elec_st%nstate = nstate
    elec_st%mxwd = mxwd

    allocate (elec_st%ntotal (nspin/mxwd))
    allocate (elec_st%totel (nspin))
    allocate (elec_st%ifmax (nspin/mxwd))
    allocate (elec_st%occ_in (nstate,nspin/mxwd))
    allocate (elec_st%sre (nspin))
    allocate (elec_st%plain_sre (nspin))
    !VDW:
    allocate (elec_st%hirsh_3d_cell(1,3))

    elec_st%ntotal(:) = int_def
    elec_st%ncharge = real_def
    elec_st%xele = real_def
    elec_st%tfermi = real_def
    elec_st%etot = real_def
    elec_st%etot_plusu = real_def
    elec_st%bdev = real_def
    elec_st%exc = real_def

    ! elec_st%dioniz = real_def - this var is initialized elsewhere.

    elec_st%dipx = real_def
    elec_st%dipy = real_def
    elec_st%dipz = real_def

    elec_st%so_pert = .false.
    !VDW:
    elec_st%do_vdw = .false.
    elec_st%evdw = real_def
    elec_st%hirsh_3d_cell(:,:)=0

  end subroutine create_electronic_struct

  subroutine electronic_struct_set_charge (ndim,elec_st)
    implicit none 
    integer, intent (in) :: ndim
    type (electronic_struct), intent (inout) :: elec_st 

    integer :: nspin
    integer :: alcstat
    nspin = elec_st%nspin
 
    allocate (elec_st%rho (ndim,2*nspin-1),stat=alcstat)
    call alccheck('rho',ndim*(2*nspin-1), alcstat)
    allocate (elec_st%rhoc (ndim),stat=alcstat)
    call alccheck('rhoc',ndim, alcstat)

    if (elec_st%ncl) then
       allocate(elec_st%spin3d(ndim,3))
       call alccheck('spin3d',ndim*3, alcstat)
       allocate(elec_st%spin3dn(ndim))
       call alccheck('spin3dn',ndim, alcstat)
    endif

  end subroutine electronic_struct_set_charge

  subroutine electronic_struct_set_charge_vdw (ndim,elec_st,atom_num)
    implicit none 
    integer, intent (in) :: ndim
    type (electronic_struct), intent (inout) :: elec_st 
    integer, intent (in) :: atom_num

    integer :: nspin
    integer :: alcstat
    integer :: num_hirsh_super 

    nspin = elec_st%nspin


    allocate (elec_st%rhof (ndim,2*nspin-1,atom_num))
    allocate (elec_st%v_free_periodic (1,atom_num))
    num_hirsh_super = (elec_st%hirsh_3d_cell(1,1)*2+1)*(elec_st%hirsh_3d_cell(1,2)*2+1)*(elec_st%hirsh_3d_cell(1,3)*2+1)

    allocate (elec_st%hirshfeld_periodic (ndim,num_hirsh_super,atom_num))
    allocate (elec_st%dist_periodic (ndim,num_hirsh_super,atom_num))


  end subroutine electronic_struct_set_charge_vdw

end module electronic_struct_module

! ==================== grid structure ========================
module grid_module
  use constants
  implicit none

  type grid_data

     ! boundary sphere size (atomic units) - for confined systems
     real(dp) :: rmax
     ! user provided  grid spacing (atomic units) - h 
     real(dp) :: stepin
     ! step = program corrected h in each direction
     ! (different from user given one only for PBC)
     ! WARNING: if PBC are not used, step must be the same on all 3
     ! directions, otherwise the boundary charge in hartset will be wrong.
     real(dp) :: step(3)
     ! hcub = h**3, hcub2 = 2/(h**3),h_2 = h**2
     real(dp) :: hcub, hcub2, h_2
     ! total number of effective grid points
     integer :: ndim
     ! total number of grid points in irreducible wedge
     integer :: nwedge
     ! number of neighbors used on one side in numerical derivative
     integer :: norder

     ! number of external neighbors
     integer :: neibs_num
     !
     ! order of double grid
     ! the number of dense-grid points between two neighbors in the
     ! original (coarser) grid is grid%ndouble - 1
     integer :: ndouble
     ! inverse order
     real(dp) :: invd
     ! coefficients in double-grid interpolation
     real(dp), dimension(:), pointer :: fdouble

     ! nxmax-nxmin+1 -- number of grid points along x-axis (including
     ! points for derivative); same for y and z directions
     integer :: nxmin, nxmax
     integer :: nymin, nymax
     integer :: nzmin, nzmax
     integer :: nxyz !
     ! n1,n2,n3 = number of grid points along x-axis, y-axis, and
     ! z-axis respectively (the ones that are actually used)
     integer :: n1, n2, n3

     ! 3d to 1d conversion table for grid point indexing - gives the
     ! 1d index based on the 3-d position (in h units, relative to the
     ! origin). Returns value of ndim+1 if point is outside the boundary sphere
     ! indexg(i,j,k) = ndim+1 for points outside the current boundary sphere
     ! indexg(i,j,k) = 1 to ndim for points inside the current boundary sphere
     integer, dimension (:,:,:), pointer :: indexg
     ! indexw(i,j,k) : has the same meaning as indexg, but it is
     ! defined only in the irreducible wedge, after considering
     ! symmetry operations in the reduced, Abelian subgroup
     integer, dimension (:,:,:), pointer :: indexw
     ! 1d to 3d conversion tables for grid point indexing - giving INDEX
     ! These arrays retrieve the three 3d indices for the ith 1d grid point
     ! 1d to 3d conversion tables for grid point indexing - giving
     ! cartesian coordinates (for irreducible wedge, replace fx with kx etc):
     ! xx = (shift(1) + kx)*h, yy=(shift(2) + ky)*h, zz=(shift(3) + kz)*h
     integer, dimension (:), pointer :: kx, ky, kz
     integer, dimension (:), pointer :: fx, fy, fz
     ! hartree_neibs_flag is a logical flag that indicates
     ! whether the arrays neibs_index, neibs_fx, neibs_fy and neibs_fz
     ! are needed for the full hartree boundary conditions
     logical :: hartree_neibs_flag
     ! neibs_index is an array with an arbitrary index of all external 
     ! neighbors , this array is used to calculate some of the boundary
     ! conditions
     integer, dimension (:,:,:), pointer :: neibs_index
     ! neib_kx, neib_ky, neib_kz - has the same meanning as kx, ky and kz
     ! but for external grid points
     integer, dimension (:), pointer :: neibs_fx, neibs_fy, neibs_fz
     !
     ! shift of grid points
     real(dp), dimension(3) :: shift
     ! irreducible wedge parameters
     integer, dimension (:), pointer :: rindex
     integer, dimension (:), pointer :: rtrans

     ! data and arrays for the calculation of laplacian in
     ! non-orthogonal grid (u,v,w).

     ! number of additional directions for calculating derivatives:
     integer lap_dir_num

     ! direction value stored in lap_dir(3). Index meaning is 1 - uv 
     ! direction, 2 uw direction, 3 vw direction. Value meaning:
     ! 0 direction is not used, 1 u+v direction is used, -1 u-v
     ! direction is used.
     integer, dimension(3) :: lap_dir

     ! step size for each direction is stored in lap_dir_step(3)
     real(dp), dimension(3) :: lap_dir_step

     ! an array with the pointers of the 3 nearest neighbors that
     ! are used
     integer, dimension(3,3) :: lap_neig

     ! laplacian coefficients values for different directories
     ! this values are calculated by the routine pbc_grid_coefs
     real(dp), dimension(6) :: b_lap

     ! inverse of the normalized lattice vectors matrix.
     real(dp), dimension(3,3) :: grad_bvec_norm

     ! coefficients for the finite difference expansion
     real(dp), dimension (:,:), pointer :: coe1, coe2
     ! coefficients for the gradient term of kinetic energy
     complex(dpc), dimension (:,:,:), pointer :: kecoe1
     ! work array for grid partitioning
     integer, dimension(:), pointer :: ist

     ! flag for domain shape in cluster BCs
     ! possible shapes:
     ! 0 - sphere, this is the default
     ! 1 - ellipsoid
     ! 2 - cylindrical
     ! 3 - box
     ! more options can be added.

     integer :: domain_shape

     ! Parameters relevant to the cluster shape
     ! The size of the array and the meaning of its parameters
     ! depend on the cluster shape.
     !
     ! Sphere: d_shape_param(1) = radius
     ! Ellipse: d_shape_param(1) = x radius
     !          d_shape_param(2) = y radius
     !          d_shape_param(3) = z radius
     ! Cylinder: d_shape_param(1) = radius
     !           d_shape_param(2) = length (centered at 0)
     !           i_shape_param(1) = orientation
     !                              (x = 1, y = 2, z = 3)
     ! Box: d_shape_param(1) = x length (centered at 0)
     !      d_shape_param(2) = y length (centered at 0)
     !      d_shape_param(3) = z length (centered at 0)
 
     real(dp), dimension(3)  :: d_shape_param
     integer, dimension(1) :: i_shape_param
     ! IMPORTANT: if the dimesions of these shape_param arrays
     ! are changed, update the broadcast statements in init_var
     !
     ! tells the different subroutines that experimental mode is on :)
     logical :: experimental
     !number of buffers to use in matvec, read in here, fed to parallel later
     integer :: max_lap_buffers


  end type grid_data

contains

  subroutine init_grid (grid)
    implicit none
    type (grid_data), intent (inout) :: grid


    nullify(grid%fdouble) 
    nullify(grid%indexg) 
    nullify(grid%indexw) 
    nullify(grid%kx) 
    nullify(grid%ky) 
    nullify(grid%kz) 
    nullify(grid%fx) 
    nullify(grid%fy) 
    nullify(grid%fz) 
    nullify(grid%rindex) 
    nullify(grid%rtrans) 
    nullify(grid%neibs_index) 
    nullify(grid%neibs_fx) 
    nullify(grid%neibs_fy) 
    nullify(grid%neibs_fz) 
    nullify(grid%coe1) 
    nullify(grid%kecoe1) 
    nullify(grid%coe2) 
    nullify(grid%ist) 

  end subroutine init_grid

  subroutine destroy_grid (grid)
    implicit none
    type (grid_data), intent (inout) :: grid

    if (associated (grid%fdouble)) deallocate (grid%fdouble) 
    if (associated (grid%indexg)) deallocate (grid%indexg) 
    if (associated (grid%indexw)) deallocate (grid%indexw)

    if (associated (grid%kx)) deallocate (grid%kx) 
    if (associated (grid%ky)) deallocate (grid%ky) 
    if (associated (grid%kz)) deallocate (grid%kz)
    
    if (associated (grid%fx)) deallocate (grid%fx)
    if (associated (grid%fy)) deallocate (grid%fy)
    if (associated (grid%fz)) deallocate (grid%fz)

    if (associated (grid%rindex)) deallocate (grid%rindex)
    if (associated (grid%rtrans)) deallocate (grid%rtrans)

    if (associated (grid%neibs_index)) deallocate (grid%neibs_index)
    if (associated (grid%neibs_fx)) deallocate (grid%neibs_fx)
    if (associated (grid%neibs_fy)) deallocate (grid%neibs_fy)
    if (associated (grid%neibs_fz)) deallocate (grid%neibs_fz)
    if (associated (grid%coe1)) deallocate (grid%coe1)
    if (associated (grid%kecoe1)) deallocate (grid%kecoe1)
    if (associated (grid%coe2)) deallocate (grid%coe2)

    if (associated(grid%ist)) deallocate(grid%ist)

  end subroutine destroy_grid

  subroutine create_grid(rmax,step,grid)
      use constants
    implicit none
    real(dp), intent (in) :: rmax
    real(dp), intent (in) :: step

    type (grid_data), intent (out) :: grid
    ! default value
    integer, parameter :: int_def = -1

    grid%rmax = rmax
    grid%stepin = step
    grid%hcub = step**3
    grid%hcub2 = two/(step**3)
    !AJB: this would seem natural to do as well:
    grid%step(:)=zero

    grid%nxmin = int_def
    grid%nxmax = int_def
    grid%nymin = int_def
    grid%nymax = int_def
    grid%nzmin = int_def
    grid%nzmax = int_def

  end subroutine create_grid

  subroutine grid_set_index(nxmin,nxmax,nymin,nymax,nzmin,nzmax,grid)
    implicit none
    integer, intent (in) :: nxmin,nxmax,nymin,nymax,nzmin,nzmax
    type (grid_data), intent (inout) :: grid

    integer :: alcstat, nxyz

    grid%nxmax = nxmax
    grid%nxmin = nxmin
    grid%nymax = nymax
    grid%nymin = nymin
    grid%nzmax = nzmax
    grid%nzmin = nzmin

    nxyz = (nxmax-nxmin+1)*(nymax-nymin+1)*(nzmax-nzmin+1)
    allocate (grid%indexg (nxmin:nxmax,nymin:nymax,nzmin:nzmax),stat=alcstat)
    call alccheck('indexg',nxyz, alcstat)
    grid%indexg(:,:,:) = 0
    allocate (grid%indexw (nxmin:nxmax,nymin:nymax,nzmin:nzmax),stat=alcstat)
    call alccheck('indexw',nxyz, alcstat)
    grid%indexw(:,:,:) = 0

    if(grid%hartree_neibs_flag) then
     allocate (grid%neibs_index (nxmin:nxmax,nymin:nymax,nzmin:nzmax), &
               stat=alcstat)
     call alccheck('neibs_index',nxyz,alcstat)
    end if
  
    grid%nxyz = nxyz

  end subroutine grid_set_index

  subroutine grid_set_ndim(ndim,grid)
    implicit none
    integer, intent (in) :: ndim
    type (grid_data), intent (inout) :: grid

    integer :: alcstat

    grid%ndim = ndim
    allocate (grid%rindex (ndim+1),stat=alcstat)
    call alccheck('grid%rindex',ndim+1, alcstat)
    allocate (grid%rtrans (ndim+1),stat=alcstat)
    call alccheck('grid%rtrans',ndim+1, alcstat)
    allocate (grid%fx(ndim),stat=alcstat)
    call alccheck('grid%fx',ndim,alcstat)
    allocate (grid%fy(ndim),stat=alcstat)
    call alccheck('grid%fy',ndim,alcstat)
    allocate (grid%fz(ndim),stat=alcstat)
    call alccheck('grid%fz',ndim,alcstat)

  end subroutine grid_set_ndim

  subroutine grid_set_ext_neibs(grid, nnodes)
    implicit none
    ! half order of finite difference
    type (grid_data), intent (inout) :: grid
    integer, intent(in) :: nnodes

    integer :: alcstat

    allocate(grid%neibs_fx (grid%neibs_num), stat=alcstat)
    call alccheck('grid%neibs_fx',grid%neibs_num,alcstat)
    allocate(grid%neibs_fy (grid%neibs_num), stat=alcstat)
    call alccheck('grid%neibs_fy',grid%neibs_num,alcstat)
    allocate(grid%neibs_fz (grid%neibs_num), stat=alcstat)
    call alccheck('grid%neibs_fz',grid%neibs_num,alcstat)


  end subroutine grid_set_ext_neibs

  subroutine grid_destroy_ext_neibs(grid,nnodes)
    implicit none
    ! half order of finite difference
    type (grid_data), intent (inout) :: grid
    integer, intent(in) :: nnodes

    if(associated(grid%neibs_fx)) deallocate(grid%neibs_fx)
    if(associated(grid%neibs_fy)) deallocate(grid%neibs_fy)
    if(associated(grid%neibs_fz)) deallocate(grid%neibs_fz)
    if(associated(grid%neibs_index)) deallocate(grid%neibs_index)

  end subroutine grid_destroy_ext_neibs


  subroutine grid_set_wedge(grid,nnodes)
    implicit none
    ! half order of finite difference
    type (grid_data), intent (inout) :: grid
    integer, intent(in) :: nnodes

    integer :: alcstat

    allocate (grid%kx (grid%nwedge),stat=alcstat) 
    call alccheck('grid%kx',grid%nwedge, alcstat)
    allocate (grid%ky (grid%nwedge),stat=alcstat)
    call alccheck('grid%ky',grid%nwedge, alcstat)
    allocate (grid%kz (grid%nwedge),stat=alcstat)
    call alccheck('grid%kz',grid%nwedge, alcstat)

    allocate(grid%ist(0:nnodes))

  end subroutine grid_set_wedge

  subroutine grid_set_ist(grid,nnodes)
    implicit none
    ! half order of finite difference
    type (grid_data), intent (inout) :: grid
    integer, intent(in) :: nnodes

    integer ii, isize, remainder

    grid%ist(:) = 0

    isize = grid%nwedge / nnodes
    remainder = grid%nwedge - isize*nnodes
    ! distribute (grid%nwedge) points so that the first PEs receive
    ! (isize+1) points and the last ones receive (isize)
    grid%ist(0) = 1
    do ii = 1, remainder
       grid%ist(ii) = grid%ist(ii-1) + isize + 1
    enddo
    do ii = remainder+1,nnodes-1
       grid%ist(ii) = grid%ist(ii-1) + isize
    enddo
    grid%ist(nnodes) = grid%nwedge + 1

  end subroutine grid_set_ist

end module grid_module

! ==================== pseudo structure ========================
module pseudo_potential_module
  use constants
  implicit none

  type  pseudo_potential

     integer, dimension (:), pointer :: format
     ! pseudopotential core cutoff for each type of atom
     real(dp), dimension (:), pointer :: rcore

     ! For each atom type:
     ! # of potentials, local component, # electron/orbit L
     integer, dimension (:), pointer :: nlocp
     integer, dimension (:), pointer :: loc
     real(dp), dimension (:,:), pointer :: eleatm
     ! initial spin polarization: ( Z_up - Z_down )/( Z_up + Z_down )
     real(dp), dimension (:), pointer :: spol
     !
     ! parameters read from ptable
     ! zion - positive charge of pseudo ion (equals # of valence electrons)
     real(dp), dimension (:), pointer :: zion
     !
     ! data read directly from the pseudopotential file(s):
     ! ns - number of points on the radial pseduopotential grid
     ! icore - 0: no core correction, 1 - core correction
     ! rs - coordinates of radial pseudopotential grid
     ! vion - local pseudopotential
     ! denc - fixed core correction charge density
     integer, dimension (:), pointer :: ns
     integer, dimension (:), pointer :: icore
     real(dp), dimension (:,:), pointer :: rs, vion, denc
     !
     ! data calculated based on info in pseudopotential files(s):
     ! dvion - radial derivative of vion
     ! d2vion - second radial derivative of vion
     ! ddenc - radial derivative of denc 
     ! d2denc - second radial derivative of denc 
     ! rho_r, drhodr, d2rhodr- radial distribution of (volume) atomic charge
     ! density and its radial derivative, respectively
     real(dp), dimension (:,:), pointer :: dvion, d2vion, ddenc, d2denc, rho_r, drhodr, d2rhodr
     ! wfspd - pseudo wavefunction (up to four angular channels)
     real(dp), dimension (:,:,:), pointer :: wfspd
     ! dwfspd - derivative of pseudo wavefunction
     real(dp), dimension (:,:,:), pointer :: dwfspd
     !
     ! quantities for Kleinman-Bylander transformation
     ! vw - wfspd*vspd/r
     ! dvw - radial derivative of vw
     ! d2vw - second radial derivative of vw
     ! ekbi - Kleinman-Bylander normalization integral 
     real(dp), dimension (:,:,:), pointer :: vw
     real(dp), dimension (:,:,:), pointer :: dvw
     real(dp), dimension (:,:,:), pointer :: d2vw
     real(dp), dimension (:,:), pointer :: ekbi

     ! The pseudopotential is given on a log grid !
     ! r(i) = aaa*(exp(bbb*(i-1))) - ccc.
     ! aaa, bbb, ccc : parameters of logarithmic pseudopotential grid
     real(dp), dimension (:), pointer :: par_a
     real(dp), dimension (:), pointer :: par_b
     real(dp), dimension (:), pointer :: par_c

     ! parameters for Fourier filtering, as defined by Briggs et al.,
     ! PRB 54, 14362 (1996). Core filter parameter is the reciprocal
     ! space cutoff (equivalent to alpha) for core charge
     real(dp), dimension (:), pointer :: alpha, beta1, acore

     ! order for spline interpolation
     logical, dimension (:), pointer :: uspline

     ! parameters for LDA+U, as in Anisimov et al., PRB 44, 943 (1991)
     ! uu - the Hubbard parameter U
     ! jj - the Hund's-rule exchange parameter J
     real(dp), dimension (:,:), pointer :: uu, jj
     ! number of atom types
     integer :: type_num
     ! mxpot = max(ns)
     integer :: mxpot
     ! number of neighbors used on one side in numerical derivative
     integer :: norder

     ! flag for spin-orbit correction
     logical :: is_so
     ! flag for spin-orbit for each atom type
     logical, dimension (:), pointer :: so

     ! flag for reading valence charge density of each atom type
     logical, dimension (:), pointer :: rvcd

     ! normalization factor so_hcub=sqrt(grid%hcub)
     real (dp) :: so_hcub

     ! spin-orbit normalization coefficients
     real(dp),  dimension (:,:), pointer :: cc

     ! radial part of (spin-orbit potential)*(radial pseudo wave function)
     real(dp), dimension (:,:,:), pointer :: vsor
     real(dp), dimension (:,:,:), pointer :: vionr

     ! derivative of radial parts
     real(dp), dimension (:,:,:), pointer :: dvr_so
     real(dp), dimension (:,:,:), pointer :: dvr_ion

     ! Wigner-Seitz radii of atoms
     real(dp), dimension (:), pointer :: rws

     ! second derivative of radial parts
     real(dp), dimension (:,:,:), pointer :: d2vr_so
     real(dp), dimension (:,:,:), pointer :: d2vr_ion

  end type pseudo_potential

contains
  subroutine init_pseudo_potential (p_pot)
    implicit none
    type (pseudo_potential), intent (inout) :: p_pot

    nullify(p_pot%format)
    nullify(p_pot%rcore)
    nullify(p_pot%nlocp)
    nullify(p_pot%loc)
    nullify(p_pot%eleatm)
    nullify(p_pot%spol)
    nullify(p_pot%zion)

    nullify(p_pot%ns)
    nullify(p_pot%icore)
    nullify(p_pot%rs)
    nullify(p_pot%denc)
    nullify(p_pot%rho_r)
    nullify(p_pot%drhodr)
    nullify(p_pot%d2rhodr)
    nullify(p_pot%vion)
    nullify(p_pot%dvion)
    nullify(p_pot%d2vion)
    nullify(p_pot%vw)
    nullify(p_pot%dvw)
    nullify(p_pot%d2vw)
    nullify(p_pot%ddenc)
    nullify(p_pot%d2denc)
    nullify(p_pot%ekbi)
    nullify(p_pot%par_a)
    nullify(p_pot%par_b)
    nullify(p_pot%par_c)
    nullify(p_pot%alpha)
    nullify(p_pot%beta1)
    nullify(p_pot%acore)
    nullify(p_pot%uspline)

    nullify(p_pot%wfspd)
    nullify(p_pot%dwfspd)
    nullify(p_pot%uu)
    nullify(p_pot%jj)

    nullify(p_pot%so)
    nullify(p_pot%rvcd)
    nullify(p_pot%cc)
    nullify(p_pot%vsor)
    nullify(p_pot%vionr)
    nullify(p_pot%dvr_so)
    nullify(p_pot%dvr_ion)
    nullify(p_pot%rws)
    nullify(p_pot%d2vr_so)
    nullify(p_pot%d2vr_ion)

  end subroutine init_pseudo_potential

  subroutine destroy_pseudo_potential (p_pot)
    implicit none
    type (pseudo_potential), intent (inout) :: p_pot

    if (associated (p_pot%format)) deallocate (p_pot%format)
    if (associated (p_pot%rcore)) deallocate (p_pot%rcore)
    if (associated (p_pot%nlocp)) deallocate (p_pot%nlocp)
    if (associated (p_pot%loc)) deallocate (p_pot%loc)
    if (associated (p_pot%eleatm)) deallocate (p_pot%eleatm)
    if (associated (p_pot%spol)) deallocate (p_pot%spol)
    if (associated (p_pot%zion)) deallocate (p_pot%zion)

    if (associated (p_pot%ns)) deallocate (p_pot%ns)
    if (associated (p_pot%icore)) deallocate (p_pot%icore)
    if (associated (p_pot%rs)) deallocate (p_pot%rs)
    if (associated (p_pot%denc)) deallocate (p_pot%denc)
    if (associated (p_pot%rho_r)) deallocate (p_pot%rho_r)
    if (associated (p_pot%drhodr)) deallocate (p_pot%drhodr)
    if (associated (p_pot%d2rhodr)) deallocate (p_pot%d2rhodr)
    if (associated (p_pot%vion)) deallocate (p_pot%vion)
    if (associated (p_pot%dvion)) deallocate (p_pot%dvion)
    if (associated (p_pot%d2vion)) deallocate (p_pot%d2vion)
    if (associated (p_pot%vw)) deallocate (p_pot%vw)
    if (associated (p_pot%dvw)) deallocate (p_pot%dvw)
    if (associated (p_pot%d2vw)) deallocate (p_pot%d2vw)
    if (associated (p_pot%ddenc)) deallocate (p_pot%ddenc)
    if (associated (p_pot%d2denc)) deallocate (p_pot%d2denc)
    if (associated (p_pot%ekbi)) deallocate (p_pot%ekbi)
    if (associated (p_pot%par_a)) deallocate (p_pot%par_a)
    if (associated (p_pot%par_b)) deallocate (p_pot%par_b)
    if (associated (p_pot%par_c)) deallocate (p_pot%par_c)
    if (associated (p_pot%alpha)) deallocate (p_pot%alpha)
    if (associated (p_pot%beta1)) deallocate (p_pot%beta1)
    if (associated (p_pot%acore)) deallocate (p_pot%acore)
    if (associated (p_pot%uspline)) deallocate (p_pot%uspline)

    if (associated (p_pot%wfspd)) deallocate (p_pot%wfspd)
    if (associated (p_pot%dwfspd)) deallocate (p_pot%dwfspd)
    if (associated (p_pot%uu)) deallocate (p_pot%uu)
    if (associated (p_pot%jj)) deallocate (p_pot%jj)

    if (associated (p_pot%so)) deallocate (p_pot%so)
    if (associated (p_pot%rvcd)) deallocate (p_pot%rvcd)
    if (associated (p_pot%cc)) deallocate (p_pot%cc)
    if (associated (p_pot%vsor)) deallocate (p_pot%vsor)
    if (associated (p_pot%vionr)) deallocate (p_pot%vionr)
    if (associated (p_pot%dvr_so)) deallocate (p_pot%dvr_so)
    if (associated (p_pot%dvr_ion)) deallocate (p_pot%dvr_ion)
    if (associated (p_pot%rws)) deallocate (p_pot%rws)
    if (associated (p_pot%d2vr_so)) deallocate (p_pot%d2vr_so)
    if (associated (p_pot%d2vr_ion)) deallocate (p_pot%d2vr_ion)

  end subroutine destroy_pseudo_potential

  subroutine create_pseudo_potential (type_num,norder,p_pot)
    implicit none
    integer, intent (inout) :: type_num
    integer, intent (inout) :: norder

    type (pseudo_potential), intent (inout) :: p_pot
    ! default values
    integer, parameter :: int_def = -1
    
    if (type_num == 0) then
    p_pot%type_num = int_def
    p_pot%norder = int_def
    p_pot%mxpot = int_def
      return
    endif

    p_pot%type_num = type_num
    p_pot%norder = norder
    p_pot%mxpot = int_def

    allocate (p_pot%format (type_num))
    allocate (p_pot%rcore (type_num))
    allocate (p_pot%nlocp (type_num))
    allocate (p_pot%loc (type_num))
    allocate (p_pot%eleatm (type_num,4))
    allocate (p_pot%spol (type_num))
    allocate (p_pot%ns (type_num))
    allocate (p_pot%icore (type_num))
    allocate (p_pot%ekbi (4,type_num))
    allocate (p_pot%par_a (type_num))
    allocate (p_pot%par_b (type_num))
    allocate (p_pot%par_c (type_num))
    allocate (p_pot%alpha (type_num))
    allocate (p_pot%beta1 (type_num))
    allocate (p_pot%acore (type_num))
    allocate (p_pot%uspline (type_num))
    allocate (p_pot%zion (type_num))
    allocate (p_pot%uu (4,type_num))
    allocate (p_pot%jj (4,type_num))
    allocate (p_pot%so (type_num))
    allocate (p_pot%rvcd (type_num))
    allocate (p_pot%cc (2,type_num))

  end subroutine create_pseudo_potential

  subroutine pseudo_potential_set_mxpot (mxpot,p_pot)
    implicit none
    integer, intent (in) :: mxpot
    type (pseudo_potential), intent (inout) :: p_pot

    ! work variables
    integer :: norder
    integer :: type_num

    p_pot%mxpot = mxpot

    norder = p_pot%norder
    type_num = p_pot%type_num

    allocate (p_pot%rs (-norder:mxpot,type_num))

    allocate (p_pot%denc (-norder:mxpot,type_num))
    allocate (p_pot%rho_r (-norder:mxpot,type_num))
    allocate (p_pot%drhodr (-norder:mxpot,type_num))
    allocate (p_pot%d2rhodr (-norder:mxpot,type_num))
    allocate (p_pot%vion (-norder:mxpot,type_num))
    allocate (p_pot%dvion (-norder:mxpot,type_num))
    allocate (p_pot%d2vion (-norder:mxpot,type_num))
    allocate (p_pot%vw (-norder:mxpot,type_num,4))
    allocate (p_pot%dvw (-norder:mxpot,type_num,4))
    allocate (p_pot%d2vw (-norder:mxpot,type_num,4))
    allocate (p_pot%ddenc (-norder:mxpot,type_num))
    allocate (p_pot%d2denc (-norder:mxpot,type_num))
    allocate (p_pot%wfspd (-norder:mxpot,type_num,4))
    allocate (p_pot%dwfspd (-norder:mxpot,type_num,4))

    allocate (p_pot%vsor (-norder:mxpot, type_num, 2))
    allocate (p_pot%vionr (-norder:mxpot, type_num, 2))
    allocate (p_pot%dvr_so (-norder:mxpot, type_num, 2))
    allocate (p_pot%dvr_ion (-norder:mxpot, type_num, 2))
    allocate (p_pot%d2vr_so (-norder:mxpot, type_num, 2))
    allocate (p_pot%d2vr_ion (-norder:mxpot, type_num, 2))

  end subroutine pseudo_potential_set_mxpot

end module pseudo_potential_module

! ==================== molecular dynamic structure ==============
module molecular_dynamic_module
  use constants
  implicit none

  type  molecular_dynamic

     ! cooling scheme (1 - stair, 2- linear, 3 - log)
     integer :: cool_type
     ! how many step to apply
     integer :: step_num      
     ! initial temperature (in kelvins)
     real(dp) :: tempi
     ! final temperature (in kelvins)
     real(dp) :: tempf
     ! step temperature provided by user and adjusted by code for stair cooling
     real(dp) :: step_tmp
     ! size of temperature step in K, calculated according to the cooling
     ! scheme
     real(dp) :: tscale
     ! time step
     real(dp) :: time_step
     ! friction coefficient
     real(dp) :: friction_coef

     ! Counter of atom movement cycles in molecular dynamics.
     ! Different from imove if a MD simulation is restarted.
     integer ::  iframe

     ! cooling stride (number of MD steps between two consecutive
     ! coolings); used for stair-step cooling only
     integer :: stride

     ! flag indicates if molecular dynamic is on, in current
     ! computation (provided by the user)
     logical :: is_on

     ! flag to turn on limited movement MD (some atoms fixed and some free)
     logical :: limited

     ! flag that indicates if this is a restart run
     logical :: is_restart

     ! Langevin trajectory variables

     ! current positions
     real(dp), dimension (:), pointer :: xcur, ycur, zcur
     ! old positions (from previous step)      
     real(dp), dimension (:), pointer :: xold, yold, zold

     ! old velocities (from previous step)
     real(dp), dimension (:), pointer :: vxold, vyold, vzold

     ! old and older accelerations (from previous two steps)
     real(dp), dimension (:), pointer ::accxold,accyold,acczold

     ! For testing/debugging: Use a simple RNG instead of system's
     logical :: use_test_rng 
     ! For testing/debugging: Override the reinitalizaion of the RNG with your own seed
     integer :: rng_seed
     ! Size of random-number array
     integer :: nrandom
     ! Position counter of random-number array. The last used random
     ! vector in vrandom is at position iset. Initialized to 0
     integer :: iset
     ! Random-number array. This array is reconstructed every time all
     ! numbers are used. See subroutine getrandom
     real(dp), dimension(:), pointer :: vrandom

  end type molecular_dynamic

contains

  subroutine init_molecular_dynamic (mol_dynamic)
    implicit none
    type (molecular_dynamic), intent (inout) :: mol_dynamic

    nullify(mol_dynamic%xcur)
    nullify(mol_dynamic%ycur)
    nullify(mol_dynamic%zcur)
    nullify(mol_dynamic%xold)
    nullify(mol_dynamic%yold)
    nullify(mol_dynamic%zold)
    nullify(mol_dynamic%vxold)
    nullify(mol_dynamic%vyold)
    nullify(mol_dynamic%vzold)
    nullify(mol_dynamic%accxold)
    nullify(mol_dynamic%accyold)
    nullify(mol_dynamic%acczold)
    nullify(mol_dynamic%vrandom)

  end subroutine init_molecular_dynamic

  subroutine destroy_molecular_dynamic (mol_dynamic)
    implicit none
    type (molecular_dynamic), intent (inout) :: mol_dynamic

    if (associated (mol_dynamic%xcur)) deallocate (mol_dynamic%xcur)
    if (associated (mol_dynamic%ycur)) deallocate (mol_dynamic%ycur)
    if (associated (mol_dynamic%zcur)) deallocate (mol_dynamic%zcur)
    if (associated (mol_dynamic%xold)) deallocate (mol_dynamic%xold)
    if (associated (mol_dynamic%yold)) deallocate (mol_dynamic%yold)
    if (associated (mol_dynamic%zold)) deallocate (mol_dynamic%zold)
    if (associated (mol_dynamic%vxold)) deallocate (mol_dynamic%vxold)
    if (associated (mol_dynamic%vyold)) deallocate (mol_dynamic%vyold)
    if (associated (mol_dynamic%vzold)) deallocate (mol_dynamic%vzold)
    if (associated (mol_dynamic%accxold)) deallocate (mol_dynamic%accxold)
    if (associated (mol_dynamic%accyold)) deallocate (mol_dynamic%accyold)
    if (associated (mol_dynamic%acczold)) deallocate (mol_dynamic%acczold)
    if (associated (mol_dynamic%vrandom)) deallocate (mol_dynamic%vrandom)

  end subroutine destroy_molecular_dynamic
      
  subroutine molecular_dynamic_turn_on (mol_dynamic,atom_num)
    implicit none

    type (molecular_dynamic), intent (inout) :: mol_dynamic
    integer, intent(in) :: atom_num

    mol_dynamic%is_restart = .false.

    allocate (mol_dynamic%xcur (atom_num))
    allocate (mol_dynamic%ycur (atom_num))
    allocate (mol_dynamic%zcur (atom_num))
    allocate (mol_dynamic%xold (atom_num))
    allocate (mol_dynamic%yold (atom_num))
    allocate (mol_dynamic%zold (atom_num))
    allocate (mol_dynamic%vxold (atom_num))
    allocate (mol_dynamic%vyold (atom_num))
    allocate (mol_dynamic%vzold (atom_num))
    allocate (mol_dynamic%accxold (atom_num))
    allocate (mol_dynamic%accyold (atom_num))
    allocate (mol_dynamic%acczold (atom_num))

    mol_dynamic%nrandom = 10000
    allocate (mol_dynamic%vrandom(mol_dynamic%nrandom))
    mol_dynamic%iset = 0

  end subroutine molecular_dynamic_turn_on

end module molecular_dynamic_module

! ==================== non local pseudo potential structure =====
module non_local_psp_module
  use constants
  implicit none

  type  nonloc_pseudo_potential

     ! maxnloc,maxnlm,atom_num
     integer :: maxnloc 
     ! number of non-local projectors over all atoms in the system
     integer :: nlm
     integer :: atom_num
     ! maxnloc in Wigner-Seits radii
     integer :: wsmaxnloc

     ! The following arrays hold information for each atom:

     ! array holding the values of the ket V_l-V_local|psi_lm>. The
     ! array only stores the value of the ket at non-local points -
     ! elsewhere it is certainly zero! There is a different ket for each 
     ! atom and for each lm value.
     ! if l=local the array just stores zeroes, as appropriate.
     real(dp), dimension (:,:), pointer ::  anloc
     ! flag for non-collinear calculation
     logical :: ncl
     ! initial intensity non-collinear magnetization 
     real(dp), dimension (:,:), pointer ::  gmag

     ! Array holding the derivates of the Kleinman-Bylander projectors.
     ! Needed in the calculation of non-local parts of the forces
     ! (this array used to be local to forcnloc.F).
     real(dp), dimension (:,:,:), pointer :: vylmd

     ! actual number of non-local points around each atom,
     ! based on the core-cutoff of that atom type and the atom
     ! coordinates.
     integer, dimension (:), pointer ::  nlatom
     integer, dimension (:), pointer ::  wsnlatom

     ! position inside the 1-d vector of each non-local point, defined
     ! in the irreducible wedge
     integer, dimension (:,:), pointer :: indw

     ! position inside the 1-d vector of each Wigner-Seitz non-local point, defined
     ! in the irreducible wedge
     integer, dimension (:,:), pointer :: wsindw

     ! symmetry operation associated to each non-local point
     ! (operation needed to bring it back to the irreducible wedge)
     integer, dimension (:,:), pointer :: tran

     ! symmetry operation associated to each non-local point inside Wigner-Seits radii
     ! (operation needed to bring it back to the irreducible wedge)
     integer, dimension (:,:), pointer :: wstran

     ! number of lm components
     integer, dimension (:), pointer :: nlmatm

     ! skbi - same as ekbi, but for each m seprately (although
     ! independent of m) and for each atom instead of each atom type.
     real(dp), dimension (:), pointer :: skbi

     ! density matrix for the on-site Coulomb interaction
     ! only used with u_pot structure.
     real(dp), dimension (:,:,:,:), pointer :: denmat
     complex(dpc), dimension (:,:,:,:), pointer :: zdenmat
     ! Constants derived from the 3-j symbols, used in the LDA+U calculation.
     ! The relation between U and 3-j symbols can be found e.g. in PRB
     ! 52, R5467 (1995)
     complex(dpc), dimension(:,:,:,:), pointer :: p3,p4
     ! spin index
     integer :: isp

     ! left and right side corrections for Vnonloc in matvec
     complex(dpc), dimension (:,:,:), pointer :: right
     complex(dpc), dimension (:,:,:), pointer :: left

     ! Hubbard U and Hund's J for LDA+U
     ! p_pot%uu/jj(:,ity) -> u_pot%uu/jj(:,ja)
     real(dp), dimension (:), pointer :: uu
     real(dp), dimension (:), pointer :: jj

     ! array holding the spin-orbit potentials
     complex(dpc), dimension (:,:,:), pointer :: v_so
     complex(dpc), dimension (:,:,:), pointer :: v_ion

     ! array holding the derivatives of spin-orbit potentials
     complex(dpc), dimension (:,:,:,:), pointer :: dv_so
     complex(dpc), dimension (:,:,:,:), pointer :: dv_ion

     ! index of atoms with spin orbit correction
     integer, dimension (:), pointer :: so_indx

     ! flag for spin-orbit psp for each atom
     logical,dimension (:), pointer :: so

     ! mxpot: length of vectors read for radial part of spin-orbit potentials
     integer :: mxpot
     integer :: so_num,so_type_num
     ! true if spin-orbit potentials are used
     logical :: is_so
     ! sign factor
     real (dp), dimension (:,:), pointer :: cc
     ! Flag for external magnetic field
     logical :: is_mag
     ! repeat of elec_st%kpts, %nkpts
     integer :: nkpt
     real(dp), dimension (:,:), pointer :: kpts

  end type nonloc_pseudo_potential

contains

  subroutine   init_nonloc_pseudo_pot (nloc_p_pot)
    implicit none
    type (nonloc_pseudo_potential), intent (inout) :: nloc_p_pot

    nullify(nloc_p_pot%anloc)
    nullify(nloc_p_pot%vylmd)
    nullify(nloc_p_pot%nlatom)
    nullify(nloc_p_pot%wsnlatom)
    nullify(nloc_p_pot%indw)
    nullify(nloc_p_pot%tran)
    nullify(nloc_p_pot%wsindw)
    nullify(nloc_p_pot%wstran)
    nullify(nloc_p_pot%nlmatm)
    nullify(nloc_p_pot%skbi)
    nullify(nloc_p_pot%gmag)

    nullify(nloc_p_pot%uu)
    nullify(nloc_p_pot%jj)
    nullify(nloc_p_pot%denmat)
    nullify(nloc_p_pot%zdenmat)
    nullify(nloc_p_pot%p3)
    nullify(nloc_p_pot%p4)

    nullify(nloc_p_pot%v_so)
    nullify(nloc_p_pot%v_ion)
    nullify(nloc_p_pot%dv_so)
    nullify(nloc_p_pot%dv_ion)
    nullify(nloc_p_pot%so_indx)
    nullify(nloc_p_pot%so)
    nullify(nloc_p_pot%cc)

  end subroutine init_nonloc_pseudo_pot

  subroutine destroy_nonloc_pseudo_pot (nloc_p_pot)
    implicit none
    type (nonloc_pseudo_potential), intent (inout) :: nloc_p_pot

    if (associated (nloc_p_pot%anloc)) deallocate (nloc_p_pot%anloc)
    if (associated (nloc_p_pot%vylmd)) deallocate (nloc_p_pot%vylmd)
    if (associated (nloc_p_pot%nlatom)) deallocate (nloc_p_pot%nlatom)
    if (associated (nloc_p_pot%wsnlatom)) deallocate (nloc_p_pot%wsnlatom)
    if (associated (nloc_p_pot%indw)) deallocate (nloc_p_pot%indw)
    if (associated (nloc_p_pot%tran)) deallocate (nloc_p_pot%tran)
    if (associated (nloc_p_pot%wsindw)) deallocate (nloc_p_pot%wsindw)
    if (associated (nloc_p_pot%wstran)) deallocate (nloc_p_pot%wstran)
    if (associated (nloc_p_pot%nlmatm)) deallocate (nloc_p_pot%nlmatm)
    if (associated (nloc_p_pot%skbi)) deallocate (nloc_p_pot%skbi)
    if (associated (nloc_p_pot%gmag)) deallocate (nloc_p_pot%gmag)

    if (associated (nloc_p_pot%uu)) deallocate (nloc_p_pot%uu)
    if (associated (nloc_p_pot%jj)) deallocate (nloc_p_pot%jj)
    if (associated (nloc_p_pot%denmat)) deallocate (nloc_p_pot%denmat)
    if (associated (nloc_p_pot%zdenmat)) deallocate (nloc_p_pot%zdenmat)
    if (associated (nloc_p_pot%p3)) deallocate (nloc_p_pot%p3)
    if (associated (nloc_p_pot%p4)) deallocate (nloc_p_pot%p4)

    if (associated (nloc_p_pot%v_so)) deallocate (nloc_p_pot%v_so)
    if (associated (nloc_p_pot%v_ion)) deallocate (nloc_p_pot%v_ion)
    if (associated (nloc_p_pot%dv_so)) deallocate (nloc_p_pot%dv_so)
    if (associated (nloc_p_pot%dv_ion)) deallocate (nloc_p_pot%dv_ion)
    if (associated (nloc_p_pot%so_indx)) deallocate (nloc_p_pot%so_indx)
    if (associated (nloc_p_pot%so)) deallocate (nloc_p_pot%so)
    if (associated (nloc_p_pot%cc)) deallocate (nloc_p_pot%cc)
  end subroutine destroy_nonloc_pseudo_pot

  subroutine create_nonloc_pseudo_pot (nlm,atom_num,nloc_p_pot)
    implicit none
    integer, intent (in) :: nlm
    integer, intent (in) :: atom_num

    type (nonloc_pseudo_potential), intent (out) :: nloc_p_pot

    nloc_p_pot%nlm = nlm

    nloc_p_pot%atom_num = atom_num

    allocate (nloc_p_pot%nlatom (atom_num))
    nloc_p_pot%nlatom(:) = 0
    allocate (nloc_p_pot%wsnlatom (atom_num))
    nloc_p_pot%wsnlatom(:) = 0
    allocate (nloc_p_pot%nlmatm (atom_num))
    nloc_p_pot%nlmatm(:) = 0
    allocate (nloc_p_pot%skbi (nloc_p_pot%nlm))
    nloc_p_pot%skbi(:) = zero
    allocate (nloc_p_pot%uu (nloc_p_pot%nlm))
    nloc_p_pot%uu(:) = zero
    allocate (nloc_p_pot%jj (nloc_p_pot%nlm))
    nloc_p_pot%jj(:) = zero

    allocate (nloc_p_pot%so (atom_num))

  end subroutine create_nonloc_pseudo_pot

  subroutine nonloc_pseudo_pot_set_maxnloc (nloc_p_pot)
    implicit none
    type (nonloc_pseudo_potential), intent (out) :: nloc_p_pot
    ! work variables
    integer :: maxnloc, wsmaxnloc
    integer :: atom_num
    integer :: nlm

    maxnloc = maxval (nloc_p_pot%nlatom)
    wsmaxnloc = maxval (nloc_p_pot%wsnlatom)
    nlm = nloc_p_pot%nlm
    atom_num = nloc_p_pot%atom_num
    nloc_p_pot%maxnloc = maxnloc
    nloc_p_pot%wsmaxnloc = wsmaxnloc

    if (associated (nloc_p_pot%anloc)) deallocate (nloc_p_pot%anloc)
    if (associated (nloc_p_pot%vylmd)) deallocate (nloc_p_pot%vylmd)
    if (associated (nloc_p_pot%indw)) deallocate (nloc_p_pot%indw)
    if (associated (nloc_p_pot%tran)) deallocate (nloc_p_pot%tran)
    if (associated (nloc_p_pot%wsindw)) deallocate (nloc_p_pot%wsindw)
    if (associated (nloc_p_pot%wstran)) deallocate (nloc_p_pot%wstran)
    if (associated (nloc_p_pot%left)) deallocate (nloc_p_pot%left)
    if (associated (nloc_p_pot%right)) deallocate (nloc_p_pot%right)
    if (associated (nloc_p_pot%v_so)) deallocate (nloc_p_pot%v_so)
    if (associated (nloc_p_pot%v_ion)) deallocate (nloc_p_pot%v_ion)
    if (associated (nloc_p_pot%dv_so)) deallocate (nloc_p_pot%dv_so)
    if (associated (nloc_p_pot%dv_ion)) deallocate (nloc_p_pot%dv_ion)
    if (associated (nloc_p_pot%so_indx)) deallocate (nloc_p_pot%so_indx)
    if (associated (nloc_p_pot%cc)) deallocate (nloc_p_pot%cc)

    if (maxnloc > 0) then
#ifdef AJB_DEBUG
     write(9,*) 'nonloac_psuedo...maxnloc:allocating anloc'
#endif
       allocate (nloc_p_pot%anloc (maxnloc,nlm))
       nloc_p_pot%anloc(:,:) = zero
#ifdef AJB_DEBUG
     write(9,*) 'nonloac_psuedo...maxnloc:allocating vylmd'
#endif
       allocate (nloc_p_pot%vylmd (3,maxnloc,nlm))
       nloc_p_pot%vylmd(:,:,:) = zero
#ifdef AJB_DEBUG
     write(9,*) 'nonloac_psuedo...maxnloc:allocating indw,tran'
#endif
       allocate (nloc_p_pot%indw (maxnloc,atom_num))
       nloc_p_pot%indw(:,:) = 0
       allocate (nloc_p_pot%tran (maxnloc,atom_num))
       nloc_p_pot%tran(:,:) = 0
       allocate (nloc_p_pot%wsindw (wsmaxnloc,atom_num))
       nloc_p_pot%wsindw(:,:) = 0
       allocate (nloc_p_pot%wstran (wsmaxnloc,atom_num))
       nloc_p_pot%wstran(:,:) = 0

       if(nloc_p_pot%nkpt > 0) then 
#ifdef AJB_DEBUG
     write(9,*) 'nonloac_psuedo...maxnloc:allocating left/right'
#endif
          allocate(nloc_p_pot%left(maxnloc,nloc_p_pot%nkpt,atom_num))
          nloc_p_pot%left = zzero
          allocate(nloc_p_pot%right(maxnloc,nloc_p_pot%nkpt,atom_num))
          nloc_p_pot%right = zzero
       endif
       
       if (nloc_p_pot%so_num > 0) then
#ifdef AJB_DEBUG
     write(9,*) 'nonloac_psuedo...maxnloc:allocating SO stuff'
#endif
          allocate (nloc_p_pot%v_so (maxnloc, 8,nloc_p_pot%so_num))
          nloc_p_pot%v_so(:,:,:) = zzero
          allocate (nloc_p_pot%v_ion (maxnloc, 8,nloc_p_pot%so_num))
          nloc_p_pot%v_ion(:,:,:) = zzero
          allocate (nloc_p_pot%dv_so (3,maxnloc, 8,nloc_p_pot%so_num))
          nloc_p_pot%dv_so(:,:,:,:) = zzero
          allocate (nloc_p_pot%dv_ion (3,maxnloc, 8,nloc_p_pot%so_num))
          nloc_p_pot%dv_ion(:,:,:,:) = zzero
          allocate (nloc_p_pot%cc (2,nloc_p_pot%so_num))
          nloc_p_pot%cc(:,:) = zero
          allocate (nloc_p_pot%so_indx (nloc_p_pot%so_num))
          nloc_p_pot%so_indx(:) = 0
       endif
    endif
#ifdef AJB_DEBUG
     write(9,*) 'finished nonloac_psuedo...maxnloc'
#endif

  end subroutine nonloc_pseudo_pot_set_maxnloc


end module non_local_psp_module

! ==================== periodic boundary conditions structure ===
module pbc_module
  use constants
  implicit none

  type  pbc_data

     ! true if periodic boundary conditions are used
     logical :: is_on
     ! number of g-vectors
     integer :: ng
     ! number of stars of g-vectors (sets of g-vectors related by symmetry)
     integer :: nstar
     ! aliases for clust%atom_num, clust%type_num
     integer :: atom_num, type_num
     ! total size of FFT mesh
     integer :: maxdfft
     ! coordinates of the first corner of FFT mesh, needed for define
     ! the map between real space and reciprocal space
     ! In slabs: pbc%mz = 0
     ! In wires: pbc%my = pbc%mz = 0
     ! In clusters: pbc%mx = pbc%my = pc%mz = 0
     integer :: mx, my, mz
     ! In bulk (fully periodic) systems: pbc%n1,n2,n3 = grid%n1,n2,n3
     ! In slabs: pbc%n1,n2 = grid%n1,n2, pbc%n3 = 1
     ! In wires: pbc%n1 = grid%n1, pbc%n2 = pbc%n3 = 1
     integer :: n1, n2, n3
     ! coordinates of g-vectors, reciprocal lattice coordinates,
     ! in ascending order of length
     integer, dimension(:,:), pointer :: kgv
     ! points to which star a g-vector belongs to
     integer, dimension(:), pointer :: inds
     ! number of g-vectors in each star
     integer, dimension(:), pointer :: mstar
     ! alpha energy (used for total energy calculation)
     real(dp) :: ealpha
     ! number of periodic dimensions (0, 1, 2, 3)
     integer :: per
     
     ! lattice vectors as read from the input file.
     !
     ! In bulk (fully periodic) systems :
     !            -                           -
     !            |   a_1_x   a_2_x   a_3_x   |
     ! latt_vec = |   a_1_y   a_2_y   a_3_y   |  vectors a_1, a_2, a_3
     !            |   a_1_z   a_2_z   a_3_z   |
     !            -                           -
     !
     ! In slabs (periodic on the xy plane only) :
     !
     !            -                           -
     !            |   a_1_x   a_2_x     0     |
     ! latt_vec = |   a_1_y   a_2_y     0     |  vectors a_1, a_2
     !            |     0       0       1     |
     !            -                           -
     !
     ! In wires (periodic on the x axis only) :
     !
     !            -                           -
     !            |   a_1_x     0       0     |
     ! latt_vec = |     0       1       0     |  vector a_1
     !            |     0       0       1     |
     !            -                           -
     !
     ! In any case, the only relevant portion of latt_vec is
     ! latt_vec(1:per , 1:per), where per = 3 (bulk), 2 (slab) or 1
     ! (wire).
     !
     real(dp) :: latt_vec(3,3)
     ! this matrix is arranged as latt_vec but the lattice
     ! vectors are normalized to unity.
     real(dp) :: avec_norm(3,3)

     ! repeat of elec_st%kpts, %nkpts
     integer :: nkpt
     real(dp), dimension (:,:), pointer :: kpts

     ! box_size = rmax
     ! box_size(1 + per:3) = 0
     real(dp) :: box_size(3)
     ! vcell = box_size(1)*box_size(2)*box_size(3) !! AJB, O.SINAI: not sure this is true for non-orthorhombic
     real(dp) :: vcell
     ! a_surface_cell = abs(det[v1 v2]) in slab geometry with lattice vectors v1,v2
     real(dp) :: a_surface_cell

     ! adot = transpose(latt_vec) * latt_vec
     real(dp) :: adot(3,3)
     ! reciprocal lattice vectors, transpose(bvec) * latt_vec = 2*pi
     real(dp) :: bvec(3,3)
     ! bdot = transpose(bvec) * bvec
     real(dp) :: bdot(3,3)

     ! shift of grid points, in units of lattice vector
     ! In slabs: pbc%shift(3) = 0
     ! In wires: pbc%shift(2:3) = 0
     ! In clusters: pbc%shift(1:3) = 0
     real(dp) :: shift(3)
     ! atomic coordinates with respect to the corner of FFT box
     real(dp), dimension(:,:,:), pointer :: rat

     ! phase factors
     complex(dpc), dimension(:), pointer :: phase
     ! conj(n) = -1 if one must take the complex conjugate of x*phase
     real(dp), dimension(:), pointer :: conj
     ! kinetic energy, |G|^2, of stars of g-vectors
     real(dp), dimension(:), pointer :: ek

     ! number of points in 1-D Fourier space mesh (pseudopotentials)
     integer, dimension(:),pointer :: nq
     ! maxdlqp = maxval(nq)
     integer :: maxdlqp
     ! Fourier-transformed local component of pseudopotentials
     real(dp), dimension(:,:), pointer :: vloc
     ! Fourier-transformed core charge
     real(dp), dimension(:,:), pointer :: dcor
     ! spacing of 1-D Fourier space mesh  (for pseudopotentials)
     real(dp), dimension(:), pointer :: delq

     ! true if the user asked to generate full dos
     logical :: create_dos
     ! number of points in DOS graph
     integer :: dos_pnum
     ! -1 - don't create LDOS, otherwise signifies the l of the
     ! angular functions used for LDOS (Ylm for l <= ylmdos_l)
     integer :: ylmdos_l
     !
     real(dp), dimension(:,:), pointer :: vql
     !
     real(dp), dimension(:,:), pointer :: dnc
     !
     real(dp), dimension(:), pointer :: vscr2
     !
     complex(dpc), dimension(:), pointer :: vscr4

  end type pbc_data

contains

  subroutine    init_pbc (pbc)
    implicit none
    type (pbc_data), intent (inout) :: pbc

    nullify(pbc%nq)
    nullify(pbc%kgv)
    nullify(pbc%inds)
    nullify(pbc%mstar)
    nullify(pbc%delq)
    nullify(pbc%phase)
    nullify(pbc%conj)
    nullify(pbc%ek)
    nullify(pbc%vql)
    nullify(pbc%dnc)
    nullify(pbc%vscr2)
    nullify(pbc%vscr4)
    nullify(pbc%rat)
    nullify(pbc%vloc)
    nullify(pbc%dcor)
    nullify(pbc%kpts)
    
    pbc%is_on = .false.

  end subroutine    init_pbc
  subroutine destroy_pbc (pbc)
    implicit none
    type (pbc_data), intent (inout) :: pbc

    if (associated (pbc%nq)) deallocate (pbc%nq)
    if (associated (pbc%kgv)) deallocate (pbc%kgv)
    if (associated (pbc%inds)) deallocate (pbc%inds)
    if (associated (pbc%mstar)) deallocate (pbc%mstar)
    if (associated (pbc%delq)) deallocate (pbc%delq)
    if (associated (pbc%phase)) deallocate (pbc%phase)
    if (associated (pbc%conj)) deallocate (pbc%conj)
    if (associated (pbc%ek)) deallocate (pbc%ek)
    if (associated (pbc%vql)) deallocate (pbc%vql)
    if (associated (pbc%dnc)) deallocate (pbc%dnc)
    if (associated (pbc%vscr2)) deallocate (pbc%vscr2)
    if (associated (pbc%vscr4)) deallocate (pbc%vscr4)
    if (associated (pbc%rat)) deallocate (pbc%rat)
    if (associated (pbc%vloc)) deallocate (pbc%vloc)
    if (associated (pbc%dcor)) deallocate (pbc%dcor)
    if (associated (pbc%kpts)) deallocate (pbc%kpts)

  end subroutine destroy_pbc

  subroutine pbc_turn_on (atom_num,type_num,pbc)
    implicit none
    integer, intent (in) :: atom_num
    integer, intent (in) :: type_num

    type (pbc_data), intent (out) :: pbc

    pbc%atom_num = atom_num
    pbc%type_num = type_num
    allocate (pbc%nq (type_num))
    allocate (pbc%delq (type_num))
    allocate (pbc%rat (3,atom_num,type_num))
  end subroutine pbc_turn_on

  subroutine pbc_set_maxdfft (maxdfft, pbc)

    implicit none
    integer, intent (in) :: maxdfft
    type (pbc_data), intent (inout) :: pbc
    integer :: alcstat

    pbc%maxdfft = maxdfft

    allocate (pbc%vscr2 (maxdfft),stat=alcstat)
    call alccheck('vscr2',maxdfft, alcstat)
    allocate (pbc%vscr4 (maxdfft),stat=alcstat)
    call alccheck('vscr4',maxdfft, alcstat)

  end subroutine pbc_set_maxdfft

  subroutine pbc_set_ng (ng, pbc)
    implicit none
    integer, intent (in) :: ng
    type (pbc_data), intent (inout) :: pbc

    pbc%ng = ng
    allocate (pbc%kgv (3,ng))
    allocate (pbc%inds (ng))
    allocate (pbc%phase (ng))
    allocate (pbc%conj (ng))
  end subroutine pbc_set_ng

  subroutine pbc_set_nstar (pbc,ns)
    implicit none
    type (pbc_data), intent (inout) :: pbc
    integer, intent (in) :: ns

    pbc%nstar = ns
    allocate (pbc%mstar (pbc%nstar))
    allocate (pbc%ek (pbc%nstar))
    allocate (pbc%vql (pbc%type_num,pbc%nstar))
    allocate (pbc%dnc (pbc%type_num,pbc%nstar))

  end subroutine pbc_set_nstar

  subroutine pbc_set_maxdlqp (maxdlqp,pbc)
    implicit none
    type (pbc_data), intent (inout) :: pbc

    integer, intent (in) :: maxdlqp

    ! work variables
    integer :: type_num

    type_num = pbc%type_num
    pbc%maxdlqp = maxdlqp

    allocate (pbc%vloc (maxdlqp, type_num))
    allocate (pbc%dcor (maxdlqp, type_num))

  end subroutine pbc_set_maxdlqp

end module pbc_module

! ==================== potential  structure =====================
module potential_module
  use constants
  implicit none

  type  potential

     ! aliases for parallel%mydim and elec_st%nspin
     integer :: ndim, nspin
     ! total new potential (from current iteration)
     real(dp), dimension (:,:), pointer :: vnew
     ! total old potential (from previous iteration)
     real(dp), dimension (:,:), pointer :: vold
     ! Hartree potential, exchange-correlation potential, ionic potential
     real(dp), dimension (:), pointer :: vhart
     real(dp), dimension (:,:), pointer :: vxc
     real(dp), dimension (:), pointer :: vion
     ! Gaussian charge distribution, for 1-dimensional (wire) systems
     real(dp), dimension (:), pointer :: rho_gauss
     ! sum of Hartree & XC potentials from previous iteration
     real(dp), dimension (:,:), pointer :: vhxcold
     ! array of coefficients for Hartree multipole expansion
     real(dp),  dimension (:,:), pointer ::  clm

     ! AMIR - changes of Doron - for current dft?

     ! current scalar potential
     real(dp), dimension (:,:), pointer :: vxcj
     ! hartree-like current vector potential
     real(dp), dimension (:,:,:), pointer :: axch
     ! exchange-like current vector potential
     real(dp), dimension (:,:,:), pointer :: axc
     ! flag for current potential
     logical :: is_cur
     ! Vector potential squared (vec{A} \cdot vec{A})
     real(dp), dimension (:), pointer :: vecpot



  end type potential

contains
  subroutine init_potential (pot)
    implicit none
    type (potential), intent (inout) :: pot 

    nullify (pot%vnew)
    nullify (pot%vold)
    nullify (pot%vhart)
    nullify (pot%vxc)
    nullify (pot%vion)
    nullify (pot%rho_gauss)
    nullify (pot%vhxcold)
    nullify (pot%clm)

    ! AMIR - changes of Doron for current dft?

    nullify (pot%vxcj)
    nullify (pot%axch)
    nullify (pot%axc)
    nullify (pot%vecpot)

  end subroutine init_potential

  subroutine destroy_potential (pot)
    implicit none
    type (potential), intent (inout) :: pot 

    if (associated (pot%vnew)) deallocate (pot%vnew)
    if (associated (pot%vold)) deallocate (pot%vold)
    if (associated (pot%vhart)) deallocate (pot%vhart)
    if (associated (pot%vxc)) deallocate (pot%vxc)
    if (associated (pot%vion)) deallocate (pot%vion)
    if (associated (pot%rho_gauss)) deallocate (pot%rho_gauss)
    if (associated (pot%vhxcold)) deallocate (pot%vhxcold)
    if (associated (pot%clm)) deallocate (pot%clm)

    ! AMIR - changes of Doron for current dft?

    if (associated (pot%vxc)) deallocate (pot%vxcj)
    if (associated (pot%vxc)) deallocate (pot%axch)
    if (associated (pot%vxc)) deallocate (pot%axc)
    if (associated (pot%vecpot)) deallocate (pot%vecpot)

  end subroutine destroy_potential

  subroutine potential_set_ndim (ndim,nspin,pot,is_cur)
    implicit none
    integer, intent (in) :: ndim, nspin

    logical, intent(in)  :: is_cur

    type (potential), intent (inout) :: pot

    integer :: alcstat

    pot%ndim = ndim
    pot%nspin = nspin

    allocate (pot%vnew (ndim,nspin),stat=alcstat)
    call alccheck('vnew',ndim*nspin, alcstat)
    allocate (pot%vold (ndim,nspin),stat=alcstat)
    call alccheck('vold',ndim*nspin, alcstat)
    allocate (pot%vhart (ndim),stat=alcstat)
    call alccheck('vhart',ndim, alcstat)
    allocate (pot%vxc (ndim,nspin),stat=alcstat)
    call alccheck('vxc',ndim*nspin, alcstat)
    allocate (pot%vion (ndim),stat=alcstat)
    call alccheck('vion',ndim, alcstat)
    allocate (pot%rho_gauss (ndim),stat=alcstat)
    call alccheck('rho_gauss',ndim, alcstat)
    allocate (pot%vhxcold (ndim,nspin),stat=alcstat)
    call alccheck('vhxcold',ndim*nspin, alcstat)
    allocate (pot%vecpot (ndim),stat=alcstat)
    call alccheck('vecpot',ndim, alcstat)
!    if (is_cur) then
!        allocate (pot%vxcj (ndim,nspin),stat=alcstat)
!        call alccheck('vxcj',ndim*nspin, alcstat)
!        allocate (pot%axch (ndim,2,3),stat=alcstat)
!        call alccheck('axch',3*ndim*nspin, alcstat)
!        allocate (pot%axc (ndim,2,3),stat=alcstat)
!        call alccheck('axc',3*ndim*nspin, alcstat)
!    endif

    pot%vhart(:) = zero
!    pot%vxc(:,:) = zero
!    pot%vxcj(:,:) = zero
!    pot%axch(:,:,:) = zero
!    pot%axc(:,:,:) = zero
    pot%vxc(:,:) = zero

  end subroutine potential_set_ndim

end module potential_module

! ==================== mixer  structure ========================
module mixer_module
  use constants
  implicit none

  type  mixer_data

     integer :: name
     real(dp)  :: param
     ! number of previous iterations used for mixing
     integer :: memory
     ! number of iterations before mixing is restarted
     integer :: restart
     ! flag set to true if scf loop fails
     logical :: scf_fail

     ! aliases for parallel%mydim and elec_st%nspin
     integer :: ndim, nspin

     ! size parameter for memory block allocated for the Broyden
     ! update arrays;  used to be kss0 in param.h
     integer :: block

     ! In Broyden, resid1 used as Initial guess for the Jacobian
     ! (assuming it is diagonal) [used to be gmix]
     real(dp), dimension (:), pointer :: resid1
     ! previous vold and second-previous vold 
     ! (updated for next run by the subroutine)
     real(dp), dimension (:), pointer :: old1 
     real(dp), dimension (:), pointer :: old2

     real(dp), dimension (:), pointer :: xin
     real(dp), dimension (:), pointer :: xout
     real(dp), dimension (:,:), pointer :: t
     real(dp), dimension (:,:), pointer :: tinv

     ! saved for Anderson mixing
     real(dp), dimension (:,:), pointer :: xinold
     real(dp), dimension (:,:), pointer :: xoutold

     ! Variables defined for multisecant methods.
     integer :: en_stage
     integer :: nsecant
     integer :: group_size
     real(dp) :: restart_factor
     real(dp) :: expand_factor
     real(dp), dimension (:,:), pointer :: dx
     real(dp), dimension (:,:), pointer :: df
     real(dp), dimension (:), pointer :: x0
     real(dp), dimension (:), pointer :: f0
     ! Above variables are used by msecant1/2/3.
     ! N is used by msecant1/3.
     real(dp), dimension (:,:), pointer :: n
     ! preferred_type, update_type, dx1, dx2, df2 are used only by msecant3.
     integer :: preferred_type, update_type
     real(dp), dimension (:,:), pointer :: dx1
     real(dp), dimension (:,:), pointer :: dx2
     real(dp), dimension (:,:), pointer :: df2

  end type mixer_data
contains

  subroutine    init_mixer (mixer)
    implicit none
    type (mixer_data), intent (inout) :: mixer 

    nullify(mixer%resid1)
    nullify(mixer%old1)
    nullify(mixer%old2)
    nullify(mixer%xin)
    nullify(mixer%xout)
    nullify(mixer%t)
    nullify(mixer%tinv)
    nullify(mixer%xinold)
    nullify(mixer%xoutold)

    ! Free memory used by multisecant methods
    nullify(mixer%dx)
    nullify(mixer%df)
    nullify(mixer%n)
    nullify(mixer%dx1)
    nullify(mixer%dx2)
    nullify(mixer%df2)
    nullify(mixer%x0)
    nullify(mixer%f0)

  end subroutine    init_mixer
  subroutine destroy_mixer (mixer)
    implicit none
    type (mixer_data), intent (inout) :: mixer 

    if (associated (mixer%resid1)) deallocate (mixer%resid1)
    if (associated (mixer%old1)) deallocate (mixer%old1)
    if (associated (mixer%old2)) deallocate (mixer%old2)
    if (associated (mixer%xin)) deallocate (mixer%xin)
    if (associated (mixer%xout)) deallocate (mixer%xout)
    if (associated (mixer%t)) deallocate (mixer%t)
    if (associated (mixer%tinv)) deallocate (mixer%tinv)
    if (associated (mixer%xinold)) deallocate (mixer%xinold)
    if (associated (mixer%xoutold)) deallocate (mixer%xoutold)

    ! Free memory used by multisecant methods
    if (associated (mixer%dx)) deallocate (mixer%dx)
    if (associated (mixer%df)) deallocate (mixer%df)
    if (associated (mixer%n)) deallocate (mixer%n)
    if (associated (mixer%dx1)) deallocate (mixer%dx1)
    if (associated (mixer%dx2)) deallocate (mixer%dx2)
    if (associated (mixer%df2)) deallocate (mixer%df2)
    if (associated (mixer%x0)) deallocate (mixer%x0)
    if (associated (mixer%f0)) deallocate (mixer%f0)

  end subroutine destroy_mixer

  subroutine set_mixer(ndim,nspin,mixer)

    use constants
    implicit none
    !
    ! Input/Output variables:
    !
    ! mixer related data
    type (mixer_data), intent(inout) :: mixer

    integer, intent (in) :: ndim
    integer, intent (in) :: nspin
    !
    ! Work variables:
    !
    integer :: dim
    integer :: alcstat
    integer :: sz

    mixer%ndim = ndim
    mixer%nspin = nspin
    dim = ndim*nspin

    allocate (mixer%xin (dim), stat=alcstat)
    call alccheck('mixer%xin',ndim, alcstat)
    allocate (mixer%xout (dim), stat=alcstat)
    call alccheck('mixer%xout',ndim, alcstat)
      
    if (mixer%name == ANDERSON) then
       allocate (mixer%xinold (1:dim,1:mixer%memory), stat=alcstat)
       call alccheck('mixer%xinold',dim*mixer%memory, alcstat)
       allocate (mixer%xoutold (1:dim,1:mixer%memory), stat=alcstat)
       call alccheck('mixer%xoutold',dim*mixer%memory, alcstat)

    elseif (mixer%name == BROYDEN) then
       allocate (mixer%resid1 (1:dim), stat=alcstat)
       call alccheck('mixer%resid1',ndim, alcstat)
       allocate (mixer%old1 (1:dim), stat=alcstat)
       call alccheck('mixer%old1',ndim, alcstat)
       allocate (mixer%old2 (1:dim), stat=alcstat)
       call alccheck('mixer%old2',ndim, alcstat)
       allocate(mixer%t(mixer%memory,mixer%memory))
       allocate(mixer%tinv(mixer%memory,mixer%memory))

    elseif (mixer%name == MSECANT1) then
       allocate (mixer%x0 (dim), stat=alcstat)
       call alccheck('mixer%x0', ndim,alcstat)
       allocate (mixer%f0 (dim), stat=alcstat)
       call alccheck('mixer%f0', ndim,alcstat)
       allocate (mixer%dx (1:dim,1:mixer%memory), stat=alcstat)
       call alccheck('mixer%dx', dim*mixer%memory, alcstat)
       allocate (mixer%df (1:dim,1:mixer%memory), stat=alcstat)
       call alccheck('mixer%df', dim*mixer%memory, alcstat)
       if (mixer%group_size == 0) then
           sz = mixer%memory
       else
           sz = mixer%group_size
       endif
       allocate (mixer%n (1:dim,1:sz), stat=alcstat)
       call alccheck('mixer%n', dim*sz, alcstat)

    elseif (mixer%name == MSECANT2) then
       allocate (mixer%x0 (dim), stat=alcstat)
       call alccheck('mixer%x0', ndim,alcstat)
       allocate (mixer%f0 (dim), stat=alcstat)
       call alccheck('mixer%f0', ndim,alcstat)
       allocate (mixer%dx (1:dim,1:mixer%memory), stat=alcstat)
       call alccheck('mixer%dx', dim*mixer%memory, alcstat)
       allocate (mixer%df (1:dim,1:mixer%memory), stat=alcstat)
       call alccheck('mixer%df', dim*mixer%memory, alcstat)

    elseif (mixer%name == MSECANT3) then
       allocate (mixer%x0 (dim), stat=alcstat)
       call alccheck('mixer%x0', ndim,alcstat)
       allocate (mixer%f0 (dim), stat=alcstat)
       call alccheck('mixer%f0', ndim,alcstat)
       allocate (mixer%dx (1:dim,1:mixer%memory), stat=alcstat)
       call alccheck('mixer%dx', dim*mixer%memory, alcstat)
       allocate (mixer%df (1:dim,1:mixer%memory), stat=alcstat)
       call alccheck('mixer%df', dim*mixer%memory, alcstat)
       if (mixer%group_size == 0) then
           if (mixer%preferred_type == 1) then
               allocate (mixer%n (1:dim,1:mixer%memory), stat=alcstat)
               call alccheck('mixer%n', dim*mixer%memory, alcstat)
           endif
       else
           allocate (mixer%n (1:dim,1:mixer%group_size), stat=alcstat)
           call alccheck('mixer%dx1', dim*mixer%group_size, alcstat)
           allocate (mixer%dx1 (1:dim,1:mixer%group_size), stat=alcstat)
           call alccheck('mixer%dx2', dim*mixer%group_size, alcstat)
           allocate (mixer%dx2 (1:dim,1:mixer%group_size), stat=alcstat)
           call alccheck('mixer%df2', dim*mixer%group_size, alcstat)
           allocate (mixer%df2 (1:dim,1:mixer%group_size), stat=alcstat)
       endif

    endif

  end subroutine set_mixer

end module mixer_module

! ==================== solver  structure ========================
module eigen_solver_module
  use constants
  implicit none

  type  eigen_solver

     ! eigensolver flag
     integer :: name

     ! diagonalization tolerance
     real(dp) :: toler

     ! dynamic diagonalization tolerance
     logical :: dyntol

     ! dynamic diagonalization tolerance
     real(dp) :: dtoler(2)
     real(dp) :: dyntoler

     ! order of multipole expansion used in computation of
     ! the Hartree potential (MAXIMUM VALUE CURRENTLY SUPPORTED - 9)
     integer :: lpole

     ! full_hartree_flag - if true - a full hartree calculation is 
     ! done for the boundaries instead of the multipole expansion. this
     ! is needed in some cases where the multipole expansion is not appropriate
     ! (for example - long cylinder).

     logical :: full_hartree_flag

     ! maximum number of matrix-vector multiplications allowed.
     integer :: maxmv
     ! number of matrix-vector multiplications performed (for
     ! statistics purposes only)
     integer :: totalmv

     ! number of additional eigenstates to be computed for each
     ! representation (important if first guess of eigenvalues is too
     ! different from the converged ones)
     integer :: nadd

     integer :: kss0

     ! window size (used in diagla and chebdav)
     integer :: winsize
     ! polynomial degree for chebff and chebdav (Chebyshev-Davidson), first diagonalization only
     integer :: polym0, ff_maxiter
     ! polynomial degree for Chebyshev filtering
     integer :: polym
     integer :: polym_t
     ! delta_degree: difference between the polynomial degree used at low
     ! energy versus high energy.
     integer :: dpm

     ! input/output flag (input is used by ARPACK only, as flag for
     ! residue restart vector)
     integer, dimension (:,:,:), pointer :: info
     ! residual restart vector, used only for arpack
     complex(dpc), dimension (:,:,:,:), pointer :: zres_res
     real(dp), dimension (:,:,:,:), pointer :: res_res
     ! residual norms, used only for diagla
     real(dp), dimension (:,:,:,:), pointer :: resn
     ! distributed diagonal part of Hamiltonian
     real(dp), dimension (:), pointer :: adiag
     ! distributed noncollinear diagonal part of Hamiltonian
     complex(dpc), dimension (:,:), pointer :: bxc
     ! alias for grid%norder, order of finite order expansion
     integer :: norder
     ! alias for grid%coe2
     real(dp), dimension (:,:), pointer :: coe2
     ! alias for grid%kecoe1
     complex(dpc), dimension (:,:,:), pointer :: kecoe1
     ! alias for rsymm%nrep
     integer :: nrep
     ! characters for current irreducible representation
     real(dp), dimension (:), pointer :: chi
     ! true if noncollinear calculation
     logical :: ncl
     ! true if Chebyshev filtering is being done
     logical :: do_subsp
     ! true if the number of eigenvalues per representation is kept fixed
     logical :: fix_neig
     ! set of Ritz values, or Kohn-Sham eigenvalues, used for
     ! Chebyshev filtering
     real(dp), dimension (:,:,:), pointer :: eval_loc
     ! lower bound of spectrum in Chebyshev filtering
     real(dp), dimension (:,:), pointer :: lowerb
     ! blocksize for everything matvec related
     integer :: mv_blksize

     integer, dimension (:,:), pointer :: ncompcnt

     ! true if a first guess of eigenvalues and eigenvectors is
     ! available for this representation and spin
     logical, dimension (:,:,:), pointer :: eig_init
     ! true if this is the first iteration in Chebyshev filtering
     logical, dimension (:,:), pointer :: firstfilt

     ! number of converged eigenvalues to be used in eigval 
     ! with chebdav routine
     integer,  dimension (:,:,:), pointer :: nconvt
     ! flag for external manetic field
     logical :: is_mag
     ! Vector potential: rot{A}=H_ext
     real(dp), dimension (:,:), pointer :: vecpot
     !experimental flag to use experimental communication mode in matvec
     logical :: experimental

  end type eigen_solver
contains

  subroutine destroy_eigen_solver (solver)
    implicit none
    type (eigen_solver), intent (inout) :: solver

    if (associated(solver%info)) deallocate(solver%info)
    if (associated(solver%nconvt)) deallocate(solver%nconvt)
    if (associated(solver%res_res)) deallocate(solver%res_res)
    if (associated(solver%zres_res)) deallocate(solver%zres_res)
    if (associated(solver%resn)) deallocate(solver%resn)
    if (associated(solver%adiag)) deallocate(solver%adiag)
    if (associated(solver%bxc)) deallocate(solver%bxc)
    if (associated(solver%coe2)) deallocate(solver%coe2)
    if (associated(solver%kecoe1)) deallocate(solver%kecoe1)
    if (associated(solver%chi)) deallocate(solver%chi)
    if (associated(solver%eval_loc)) deallocate(solver%eval_loc)
    if (associated(solver%lowerb)) deallocate(solver%lowerb)
    if (associated(solver%ncompcnt)) deallocate(solver%ncompcnt)
    if (associated(solver%eig_init)) deallocate(solver%eig_init)
    if (associated(solver%firstfilt)) deallocate(solver%firstfilt)
    if (associated(solver%vecpot)) deallocate(solver%vecpot)

  end subroutine destroy_eigen_solver

  subroutine create_eigen_solver (solver,grid,ldn,nrep,nspin,nstate, &
       nkpt_in,cplx,mxwd)
    use grid_module
    implicit none
    type (eigen_solver), intent (out) :: solver
    type (grid_data), intent(in) :: grid
    integer, intent(in) :: ldn,nrep,nspin,nstate,nkpt_in,mxwd
    logical, intent(in) :: cplx

    integer :: alcstat,nmax,ntmp,nkpt
    
    if (solver%name == TEST) then
        !we dont use the eigensolver
         return !?
    endif
    !
    ! nkpt = 0 if kpoints are not used; must increase it to one
    nkpt = max(nkpt_in,1)
    !
    ! should be inserted by the user inside parsec.in, maximum is 9
    solver%lpole = 9
    ! should also be inserted by the user inside parsec.in
    solver%kss0 = 20
    ! diagla/chebdav need extra window
    if (solver%name == DIAGLA) then
       solver%winsize = 5
    elseif (solver%name == CHEBDAV .or. solver%name == CHEBFF) then
       solver%winsize = 12
    else
       solver%winsize = 0
    endif
    solver%totalmv = 0

    solver%nrep = nrep

    allocate(solver%chi(nrep))
    allocate(solver%adiag(ldn),stat =alcstat)
    call alccheck('solver%adiag',ldn, alcstat)
    if (mxwd == 2) then
       allocate(solver%bxc(ldn,2),stat =alcstat)
       call alccheck('solver%bxc',2*ldn, alcstat)
    endif
    allocate(solver%info(nrep,nkpt,nspin))
    solver%info = 0

    solver%norder = grid%norder
    ntmp=3+grid%lap_dir_num
    allocate(solver%coe2(-solver%norder:solver%norder,ntmp))
    solver%coe2(:,:) = grid%coe2(:,:)

    if (nkpt_in /= 0) then
       allocate(solver%kecoe1(nkpt_in,6,0:solver%norder))
       solver%kecoe1(:,:,:)=grid%kecoe1(:,:,:)
    endif

    allocate(solver%eig_init(nrep,nkpt,nspin))

    if (solver%do_subsp) then
       solver%eig_init(:,:,:) = .false.
       allocate(solver%firstfilt(nrep,nspin))
       solver%firstfilt(:,:) = .true.
       allocate(solver%lowerb(nrep,nspin))
       allocate(solver%ncompcnt(nrep,nspin))
       solver%ncompcnt(:,:) = 0
       nmax = nstate + solver%nadd + solver%winsize
       allocate(solver%eval_loc(nmax,nrep,nspin))
       solver%eval_loc(:,:,:) = zero
    endif

    if (solver%name == ARPACK) then
       if (cplx) then
          allocate(solver%zres_res(ldn*mxwd,nrep,nkpt,nspin),stat=alcstat)
          call alccheck('solver%zres_res',ldn*mxwd*nrep*nkpt*nspin,alcstat)
          solver%zres_res(:,:,:,:) = zzero
       else
          allocate(solver%res_res(ldn,nrep,nkpt,nspin),stat=alcstat)
          call alccheck('solver%res_res',ldn*nrep*nkpt*nspin,alcstat)
          solver%res_res(:,:,:,:) = zero
       endif
    endif

    ! initialize nconvt to zero
    allocate(solver%nconvt(nrep,nkpt,nspin))
    solver%nconvt(:,:,:) = 0
    if(solver%is_mag) then
      allocate(solver%vecpot(ldn,3),stat =alcstat)
      call alccheck('solver%vecpot',3*ldn, alcstat)
      solver%vecpot = zero
    endif
!    solver%experimental = .false.
  end subroutine create_eigen_solver
  
end module eigen_solver_module

! ==================== movement data structure ==================
module movement_module
  use constants
  implicit none

  type movement

     ! maximum number of movements allowed
     integer :: mxmove
     ! counter for number of movements
     integer :: num
     ! type of movements: 1 - Simple, 2 - BFGS, 3 - manual (see const.f)
     integer :: name
     ! max step(au), min step (au),  force minimum
     real(dp)  :: stepmax, stepmin, fmin
     ! true if movement of atoms is being done
     logical :: is_on
     ! 'was minimization achieved' flag
     logical :: min_frc_found
     ! 'routine is stopping' flag
     logical :: done
     ! coordinates of the movable atoms
     ! in systems without periodic boundary conditions, this is
     ! equivalent to coord%xatm,coord%yatm,coord%zatm; in systems with
     ! periodic boundary conditions, they may differ by a unit lattice vector
     real(dp), dimension (:,:), pointer :: rcur
     ! true if structural relaxation is done restarting from previous run
     logical :: is_restart
     ! number of corrections used in the BFGS update
     integer :: mmax
     ! lower and upper bounds for the positions of all atoms
     ! (restricted BFGS only)
     real(dp), dimension (:), pointer :: lower, upper
     ! workspace arrays for BFGS
     character (len=120) :: cbfgs
     logical, dimension(:), pointer :: lbfgs
     integer, dimension (:), pointer :: ibfgs
     real(dp), dimension (:), pointer :: dbfgs, wbfgs

  end type movement
contains

  subroutine create_movement (move)
    implicit none
    type (movement), intent (inout) :: move

    move%min_frc_found = .false.
    move%done = .false.
    move%num = 0
  end subroutine create_movement

  subroutine init_movement (move)
    implicit none
    type (movement), intent (inout) :: move

    nullify(move%rcur)
    nullify(move%lower)
    nullify(move%upper)
    nullify(move%lbfgs)
    nullify(move%ibfgs)
    nullify(move%dbfgs)
    nullify(move%wbfgs)
  end subroutine init_movement

  subroutine destroy_movement (move)
    implicit none
    type (movement), intent (inout) :: move
    if(associated(move%rcur))  deallocate(move%rcur)
    if(associated(move%lower))  deallocate(move%lower)
    if(associated(move%upper))  deallocate(move%upper)
    if(associated(move%lbfgs))  deallocate(move%lbfgs)
    if(associated(move%ibfgs))  deallocate(move%ibfgs)
    if(associated(move%dbfgs))  deallocate(move%dbfgs)
    if(associated(move%wbfgs))  deallocate(move%wbfgs)

  end subroutine destroy_movement

end module movement_module

! ==================== symmetry data structure ==================
module symmetry_module
  use constants
  implicit none

  type symmetry

     ! true if symmetry operations are used
     logical :: use_symm
     ! number of point group operations
     integer :: ntrans
     ! integer indicating whether hexagonal or cubic symmetry
     integer :: cell_symmetry
     ! matrix of rotations in g-lattice coordinates
     integer, dimension(:,:,:), pointer :: gmtrx
     ! matrix of rotations in r-lattice coordinates
     integer, dimension(:,:,:), pointer :: rmtrx
     ! matrix of rotations in real space
     real(dp), dimension(:,:,:), pointer :: trans
     ! fractional translation in lattice coordinates
     real(dp), dimension(:,:), pointer :: tnp
     ! unit lattice vectors for the supercell (if no PBC are used, the
     ! supercell will be a cubic box that encloses the system)
     real(dp) :: alatt(3,3)
     ! inverse of alatt matrix (needed for the transformation of
     ! coordinates from cartesian to lattice vector units)
     real(dp) :: invlat(3,3)
     ! character table (used only for Abelian subgroup, so it is always
     ! integer)
     integer, dimension(:,:), pointer :: chi

  end type symmetry

contains
  subroutine init_symmetry (symm)
    implicit none
    type (symmetry), intent (inout) :: symm
    nullify(symm%gmtrx)
    nullify(symm%rmtrx)
    nullify(symm%trans)
    nullify(symm%tnp)
    nullify(symm%chi)
  end subroutine init_symmetry

  subroutine destroy_symmetry (symm)
    implicit none
    type (symmetry), intent (inout) :: symm
    if (associated (symm%gmtrx)) deallocate (symm%gmtrx)
    if (associated (symm%rmtrx)) deallocate (symm%rmtrx)
    if (associated (symm%trans)) deallocate (symm%trans)
    if (associated (symm%tnp)) deallocate (symm%tnp)
    if (associated (symm%chi)) deallocate (symm%chi)
  end subroutine destroy_symmetry

  subroutine create_symmetry (ntrans,symm)
    implicit none
    integer, intent (inout) :: ntrans
    type (symmetry), intent (inout) :: symm

    symm%ntrans = ntrans
    allocate (symm%gmtrx (3,3,ntrans))
    allocate (symm%rmtrx (3,3,ntrans))
    allocate (symm%trans (3,3,ntrans))
    allocate (symm%tnp (3,ntrans))
    allocate (symm%chi (ntrans,ntrans))

  end subroutine create_symmetry

end module symmetry_module

! ==================== parallel data structure ==================

module parallel_data_module
  use constants
  implicit none

  type parallel_data

     ! total number of processors
     integer :: procs_num
     ! ID of master
     integer :: masterid
     ! true for master processor, false for all other ones
     logical :: iammaster
     ! rank of all the processors in the groupms
     integer :: iam
     ! communicator for inter-processor communication
     integer :: comm

     ! handle of the world group (containing all processors)
     integer :: world_handle
     ! total number of groups
     integer :: groups_num
     ! number of processors per group ( group_size * groups_num = procs_num)
     integer :: group_size
     ! order of group to which each processor belongs to
     integer :: mygroup
     ! true for the master processors of this group (each group has one
     ! and only one master)
     logical :: iamgmaster
     ! rank of the master processor of this group
     integer :: group_master
     ! rank of processors within a group
     integer :: group_iam
     ! handle of a group
     integer :: group_handle
     ! communicator within a group
     integer :: group_comm
     ! optimized comm for matvec
     integer :: group_comm_topo
     ! plan of the groups: all processors belonging to j-th group
     ! have parallel%iam = parallel%gmap(parallel%group_iam,j)
     integer, dimension(:,:), pointer :: gmap
     ! rank of processors in the group of masters
     integer :: gmaster_iam
     ! handle of the group of masters
     integer :: gmaster_handle
     ! communicator in the group of masters
     integer :: gmaster_comm

     ! actual size of hamiltonian
     integer ndim, nwedge
     ! local size of hamiltonian (equal to number of grid points
     ! held by each processor)
     integer mydim
     ! dimension of distributed grid arrays; in order to account
     ! for points outside the domain, choose ldn > max(numb)
     integer ldn
     ! alias for grid%lap_dir_num
     integer :: lap_dir_num

     ! work array for functions on the irreducible wedge,
     ! allocated by master PE only!
     real(dp), pointer :: ftmp(:)
     complex(dpc), pointer :: zftmp(:)

     ! Position of the first block of rows (=points in irreducible
     ! wedge) that each processor holds. The number of rows in
     ! processor ipe is irows(ipe+1)-irows(ipe). Global
     integer, dimension (:), pointer :: irows
     ! arrays for inter-processor communication
     ! they are all local (i.e., different values on different procs)
     ! see more about the meaning of these variables in comm_neigh.F
     integer, dimension(:), pointer :: ip1, jp1, irecvp, jsendp
     integer, dimension (:), pointer :: senrows
     integer, dimension (:), pointer :: pint
     integer :: inter1,inter2
     ! neighbor points to each point belonging to a computing PE;
     ! this index array indicates the position, in the 1-d array, 
     ! of the ith neighbor of point j (i from 1 to norder*3, because
     ! there are three directions, j from 1 to ndim). used for
     ! calculating derivatives; notice that neibs points to the
     ! index in irreducible wedge, not full grid!
     integer, dimension (:,:), pointer :: neibs
     ! index array indicating the symmetry operation needed to
     ! bring this neighor point back to the irreducible wedge
     ! NOTE: tneibs should be used whenever generic functions are
     ! evaluated at a neighbor point, with the exception of totally
     ! symmetric functions (like charge density, potentials etc.);
     ! that is because totally symmetric functions have character one
     ! (i.e., they do not gain phase factors upon a symmetry operation)
     integer, dimension (:,:), pointer :: tneibs
     ! local size for communication buffers (edges,size of send,size of receive)
     integer :: countcomm, maxcomm, maxcomm2
     ! local information about communication:
     integer, allocatable, dimension (:) :: &
         sources,destinations,sendcounts,recvcounts,rdispls,sdispls
     !info about laplacian buffers used for communication:
     integer :: max_lap_buffers
     ! additional spin dimension of wave-function due to spin-orbit
     ! and\or non-collinear states
     integer :: mxwd,nshift

  end type parallel_data
contains

  subroutine init_parallel_data(parallel)
    implicit none
    type (parallel_data), intent (inout) :: parallel
     nullify (parallel%irows)
     nullify (parallel%ftmp)
     nullify (parallel%zftmp)
     nullify (parallel%ip1)
     nullify (parallel%jp1)
     nullify (parallel%irecvp)
     nullify (parallel%jsendp)
     nullify (parallel%senrows)
     nullify (parallel%pint)
     nullify (parallel%neibs)
     nullify (parallel%tneibs)
     nullify (parallel%gmap)

  end subroutine init_parallel_data

  subroutine destroy_parallel_data (parallel)
    implicit none
    type (parallel_data), intent (inout) :: parallel

    if (associated (parallel%irows)) deallocate (parallel%irows)
    if (associated (parallel%ftmp)) deallocate (parallel%ftmp)
    if (associated (parallel%zftmp)) deallocate (parallel%zftmp)
    if (associated (parallel%ip1)) deallocate (parallel%ip1)
    if (associated (parallel%jp1)) deallocate (parallel%jp1)
    if (associated (parallel%irecvp)) deallocate (parallel%irecvp)
    if (associated (parallel%jsendp)) deallocate (parallel%jsendp)
    if (associated (parallel%senrows)) deallocate (parallel%senrows)
    if (associated (parallel%pint)) deallocate (parallel%pint)
    if (associated (parallel%neibs)) deallocate (parallel%neibs)
    if (associated (parallel%tneibs)) deallocate (parallel%tneibs)
    if (associated (parallel%gmap)) deallocate (parallel%gmap)

  end subroutine destroy_parallel_data
!===================================================================!
  subroutine create_parallel_data (parallel)
#ifdef MPI
! include the mpi parameters file
! 
    use mpi
#endif
    implicit none
    type (parallel_data), intent (out) :: parallel
#ifdef MPI
    integer mpinfo
#endif
  ! check for required thread support
    integer thread_provided
    integer thread_requested

#if OMPFUN && MPI
    thread_requested = MPI_THREAD_FUNNELED
#elif MPI
    thread_requested = MPI_THREAD_SINGLE
#endif

#ifdef MPI
    ! Initialise MPI, get size and my id, create communicator
    call MPI_INIT_THREAD(thread_requested,thread_provided,mpinfo)
    !    call MPI_INIT(mpinfo)
    if (thread_provided < thread_requested) then
        write(7,*) 'Error: MPI init does not support requested thread level '
         call MPI_ABORT(MPI_COMM_WORLD,-1,mpinfo)
    endif
    ! Start a communicator for the entire bunch
    parallel%comm = MPI_COMM_WORLD
    ! Determine the size of this group. This group should have total
    ! number of processors as its size
    call MPI_COMM_SIZE(parallel%comm, parallel%procs_num, mpinfo)
    ! Among this group define the processor with rank 0 as the master
    ! processor
    parallel%masterid = 0
    ! Now determine the ranks
    call MPI_COMM_RANK(parallel%comm, parallel%iam, mpinfo)
#else
    ! For a non-MPI calculation, define degenerate environment
    parallel%procs_num = 1
    parallel%masterid = 0
    parallel%comm = 0
    parallel%iam = 0
#endif
    if (parallel%masterid == parallel%iam) then
       parallel%iammaster = .true.
    else
       parallel%iammaster = .false.
    endif

    parallel%mxwd = 1
    parallel%nshift = 1

  end subroutine create_parallel_data

!=============================Intel TraceAnalyzer API ======================================
#ifdef ITAC
  subroutine init_trace_data (vtierr)
    include 'VT.inc'
    include 'vtcommon.inc'
      integer vtierr
      vtierr = 0
!define the classes of state handles
      call VTCLASSDEF( '"Main"',       calc_class ,vtierr )
      call VTCLASSDEF( 'Eigval',     eig_class ,vtierr )
      call VTCLASSDEF( 'Chebdav',      chebdav_class ,vtierr )
      call VTCLASSDEF( 'ChebFF',       chebff_class ,vtierr )
      call VTCLASSDEF( 'Subspace',     subspace_class ,vtierr )
      call VTCLASSDEF( 'Matvec',       mvec_class ,vtierr )
      call VTCLASSDEF( 'Setup',        setup_class ,vtierr )
      call VTCLASSDEF( 'Post-Calc',    post_class ,vtierr )

!define basic states handles
      call VTFUNCDEF( 'Grid_Setup',    setup_class,  vt_grid_setup_state     , vtierr)
      call VTFUNCDEF( 'Other_Setup',  setup_class,  vt_pre_scf_state        , vtierr)

      call VTFUNCDEF( 'HXC',           calc_class,  vt_hart_xc        , vtierr)
      call VTFUNCDEF( 'Other_lap',           calc_class,  vt_lapmvs        , vtierr)
      call VTFUNCDEF( 'NonLocal',      calc_class,  vt_nonloc        , vtierr)
      call VTFUNCDEF( 'SCF',           calc_class,  vt_scf_state             , vtierr)
      call VTFUNCDEF( 'PostSCF',       calc_class,  vt_post_scf_state        , vtierr)
      call VTFUNCDEF( 'Ions',          calc_class,  vt_pbc_ion               , vtierr)
      call VTFUNCDEF( 'Forces',        calc_class,  vt_forces_state          , vtierr)

      call VTFUNCDEF( 'ARPACK', eig_class,  vt_arpack    , vtierr)
      call VTFUNCDEF( 'ARPACK:work', eig_class,  vt_naupd    , vtierr)
      call VTFUNCDEF( 'ARPACK:post_processing', eig_class,  vt_seupd    , vtierr)
      call VTFUNCDEF( 'Spin-Orbit?',   eig_class,  vt_so            , vtierr)
      call VTFUNCDEF( 'SortEVAL',      eig_class,  vt_eigensort_state        , vtierr)
      call VTFUNCDEF( 'SendEIG',       eig_class,  vt_eigencomm_state        , vtierr)

      call VTFUNCDEF( 'Orthnormal',    eig_class,  vt_orthnormal , vtierr)
      call VTFUNCDEF( 'DGKS',          eig_class,  vt_dgks , vtierr)
      call VTFUNCDEF( 'Orth-2ndpass',  eig_class,  vt_gemv    , vtierr)
      call VTFUNCDEF( 'Lanczos',       eig_class,  vt_lancz          , vtierr)
!
      call VTFUNCDEF( 'H_nl.p',     mvec_class,  vt_matvec_ps             , vtierr)
      call VTFUNCDEF( 'H.p',    mvec_class,  vt_matvec          , vtierr)
      call VTFUNCDEF( 'preBDXC',       mvec_class,  vt_bdxc_leftover         , vtierr)
      call VTFUNCDEF( 'BDXC',          mvec_class,  vt_bdxc                  , vtierr)
      call VTFUNCDEF( 'Buffering',          mvec_class,  vt_buffer                  , vtierr)
!

      call VTFUNCDEF( 'Chebff_Main',   chebff_class,  vt_chebff, vtierr)
      call VTFUNCDEF( 'Update',        chebff_class,  vt_chebff_update    , vtierr)
      call VTFUNCDEF( 'Decomp',        chebff_class,  vt_chebff_decomp     , vtierr)
      call VTFUNCDEF( 'Filter',        chebff_class,  vt_chebff_filter     , vtierr)

      call VTFUNCDEF( 'Chebdav_Main',  chebdav_class,  vt_chebdav            , vtierr)
      call VTFUNCDEF( 'VHV'    ,       chebdav_class,  vt_chebdav_MM         , vtierr)
      call VTFUNCDEF( 'Decomp',        chebdav_class,  vt_chebdav_decomp     , vtierr)
      call VTFUNCDEF( 'Update_I',      chebdav_class,  vt_chebdav_Iupdate    , vtierr)
      call VTFUNCDEF( 'Update_II',     chebdav_class,  vt_chebdav_IIupdate   , vtierr)
      call VTFUNCDEF( 'Reorder',       chebdav_class,  vt_chebdav_reorder    , vtierr)
!
      call VTFUNCDEF( 'Subspace_Main',      subspace_class,  vt_subspace_main          , vtierr)
      call VTFUNCDEF( 'Filter',      subspace_class,  vt_subspace_filter          , vtierr)
      call VTFUNCDEF( 'Update',      subspace_class,  vt_subspace_update          , vtierr)

      call VTFUNCDEF( 'Move',          post_class,  vt_move_state            , vtierr)
      call VTFUNCDEF( 'Postprocessing',post_class,  vt_post_processing_state , vtierr)
!


  end subroutine init_trace_data
#endif
!========================================================================================
  subroutine create_group_layout(parallel,groups_num_in)
#ifdef MPI
! include the mpi parameters file
! 
    use mpi
#endif
    implicit none
    type (parallel_data), intent (inout) :: parallel
    integer, intent(in) :: groups_num_in

    integer groups_num
    integer ii, mykey
#ifdef MPI
    integer  mpinfo, group_comm, group_handle, igrp
#endif
    integer, allocatable :: master_group(:)

    if (parallel%groups_num == 0) then
       groups_num = groups_num_in
    else
       groups_num = parallel%groups_num
    endif

    do ii = groups_num, 1, -1
       if (mod(parallel%procs_num,ii) == 0 .and. &
            mod(groups_num ,ii) == 0) then
          parallel%groups_num = ii
          parallel%group_size = parallel%procs_num / parallel%groups_num
          exit
       endif
    enddo
    parallel%mygroup = parallel%iam/parallel%group_size
    mykey = mod(parallel%iam,parallel%group_size)
    allocate(parallel%gmap(parallel%group_size,parallel%groups_num))
    parallel%gmap = 0
    parallel%gmap(mykey + 1,parallel%mygroup + 1) = parallel%iam
    call pisum(parallel%gmap,parallel%procs_num, &
         parallel%procs_num,parallel%comm)

#ifdef MPI
    call MPI_COMM_GROUP(parallel%comm,parallel%world_handle,mpinfo)
    allocate(master_group(0:parallel%group_size-1))
    do igrp = 0, parallel%groups_num - 1
       do ii = 0, parallel%group_size - 1
          master_group(ii) = parallel%gmap(ii + 1,igrp + 1)
       enddo
       call MPI_GROUP_INCL(parallel%world_handle,parallel%group_size, &
            master_group,group_handle,mpinfo)
       call MPI_COMM_CREATE(parallel%comm,group_handle, group_comm,mpinfo)
       if (igrp == parallel%mygroup) then
          parallel%group_handle = group_handle
          parallel%group_comm = group_comm
          call MPI_GROUP_RANK(parallel%group_handle, parallel%group_iam,mpinfo)
       endif
    enddo
    deallocate(master_group)
#else
    parallel%group_iam = 0
    parallel%group_comm = 0
    parallel%mygroup = 0
#endif
    !
    ! Define masters within each group.
    ! Warning: the master PE *must* be master in its group!
    !
    parallel%group_master = 0

    if (parallel%group_iam == parallel%group_master) then
       parallel%iamgmaster = .true.
    else
       parallel%iamgmaster = .false.
    endif
    allocate(master_group(0:parallel%groups_num-1))
    do ii = 0, parallel%groups_num - 1
       master_group(ii) = ii*parallel%group_size + parallel%group_master
    enddo
#ifdef MPI
    call MPI_GROUP_INCL(parallel%world_handle ,parallel%groups_num, &
         master_group ,parallel%gmaster_handle,mpinfo)
    call MPI_COMM_CREATE(parallel%comm,parallel%gmaster_handle, &
         parallel%gmaster_comm,mpinfo)
    if (parallel%iamgmaster) call MPI_GROUP_RANK(parallel%gmaster_handle, &
         parallel%gmaster_iam,mpinfo)
#else
    parallel%gmaster_iam = 0
    parallel%gmaster_comm = 0
#endif
    deallocate(master_group)

    if (parallel%iammaster) then
       write(7,*) (' ',ii=1,20),('*',ii=1,20),(' ',ii=1,20)
       write(7,*)
       write(7,*) ' Distributing ',parallel%procs_num, &
            ' processors among ',parallel%groups_num,' groups.'
       write(7,*) ' Each group contains ',parallel%group_size ,' processors.'
       write(7,*)
       write(7,*) (' ',ii=1,20),('*',ii=1,20),(' ',ii=1,20)
       call myflush(7)
    endif
    write(9,*) (' ',ii=1,20),('*',ii=1,20),(' ',ii=1,20)
    write(9,*)
    write(9,*) ' Group partition information: '
    write(9,*) ' I am processor rank ',parallel%group_iam , &
         ' in group ',parallel%mygroup
    if (parallel%iamgmaster) then
       write(9,*) ' I am master in my group, with master rank ' , &
            parallel%gmaster_iam
    else
       write(9,*) ' I am not master in my group.'
    endif
    write(9,*)
    write(9,*) (' ',ii=1,20),('*',ii=1,20),(' ',ii=1,20)
    call myflush(9)

  end subroutine create_group_layout
!=================================================================================
  subroutine export_function(parallel,f_distrib)
    !
    ! Distributes array parallel%ftmp from master PE to the other
    ! processors. The distribution follows grid partition defined
    ! in parallel%irows. Ordinary processors keep
    ! the distributed data in array f_distrib.
    !
#ifdef MPI
    use mpi
#endif
    implicit none
    type (parallel_data), intent (inout) :: parallel
    real(dp), intent(inout) :: f_distrib(parallel%mydim)

    ! communication variables
    integer node, nelems, mpinfo
#ifdef MPI
    integer status(MPI_STATUS_SIZE)
#endif

    f_distrib(:) = zero
    !
    ! First, master processor broadcasts copies of parallel%ftmp to
    ! the masters within each group.
#ifdef MPI
    if (parallel%iamgmaster) then
       call MPI_BCAST(parallel%ftmp,parallel%nwedge,MPI_DOUBLE_PRECISION, &
            parallel%masterid,parallel%gmaster_comm,mpinfo)
    endif
#endif

    do node = 1, parallel%group_size
       nelems = parallel%irows(node) - parallel%irows(node-1)
       if (parallel%iamgmaster) then
          if (node == parallel%group_iam + 1) then
             call dcopy(nelems,parallel%ftmp(parallel%irows(node-1)),1, &
                  f_distrib,1)
          else
#ifdef MPI
             call MPI_SEND(parallel%ftmp(parallel%irows(node-1)), &
                  nelems,MPI_DOUBLE_PRECISION, &
                  node-1,node,parallel%group_comm,mpinfo)
#endif
          endif
       else
          if (node == parallel%group_iam + 1) then
#ifdef MPI
             call MPI_RECV(f_distrib,nelems,MPI_DOUBLE_PRECISION, &
                  parallel%group_master,node,parallel%group_comm,status,mpinfo)
#endif
          endif
       endif
#ifdef MPI
       call MPI_Barrier(parallel%group_comm,mpinfo)
#endif
    enddo

  end subroutine export_function

!=================================================================================
  subroutine export_zfunction(parallel,f_distrib)
    !
    ! Similar to export_function, but for complex arrays with size
    ! parallel%mydim * parallel%mxwd.
    !
#ifdef MPI
    use mpi
#endif
    implicit none
    type (parallel_data), intent (in) :: parallel
    complex(dpc), intent(inout) :: f_distrib(parallel%mydim*parallel%mxwd)

    ! counters
    integer :: ii, istart
    ! communication variables
    integer node, nelems, mpinfo
    integer :: offset
#ifdef MPI
    integer status(MPI_STATUS_SIZE)
#endif

    f_distrib(:) = zzero
    !
    ! First, master processor broadcasts copies of parallel%zftmp to
    ! the masters within each group.
    !
#ifdef MPI
    if (parallel%iamgmaster) then
       call MPI_BCAST(parallel%zftmp,parallel%nwedge*parallel%mxwd, &
            MPI_DOUBLE_COMPLEX, &
            parallel%masterid,parallel%gmaster_comm,mpinfo)
    endif
#endif

    do ii = 1, parallel%mxwd
       istart = 1 + (ii - 1)*parallel%mydim
       do node = 1, parallel%group_size
          nelems = parallel%irows(node) - parallel%irows(node-1)
          offset = parallel%irows(node-1) + (ii - 1)*parallel%nwedge
          if (parallel%iamgmaster) then
             if (node == parallel%group_iam + 1) then
                call zcopy(nelems,parallel%zftmp(offset),1,f_distrib(istart),1)
             else
#ifdef MPI
                call MPI_SEND(parallel%zftmp(offset), &
                     nelems,MPI_DOUBLE_COMPLEX, &
                     node-1,node,parallel%group_comm,mpinfo)
#endif
             endif
          else
             if (node == parallel%group_iam + 1) then
#ifdef MPI
                call MPI_RECV(f_distrib(istart),nelems,MPI_DOUBLE_COMPLEX, &
                     parallel%group_master,node,parallel%group_comm,status,mpinfo)
#endif
             endif
          endif
#ifdef MPI
          call MPI_Barrier(parallel%group_comm,mpinfo)
#endif
       enddo
    enddo
  end subroutine export_zfunction
!=================================================================================

  subroutine collect_function(parallel,f_distrib)
    !
    ! Collects array parallel%ftmp to master PE. The
    ! distribution follows grid partition defined in parallel%irows.
    ! Input dat is in array f_distrib
    !
#ifdef MPI
    use mpi
#endif
    implicit none
    type (parallel_data), intent (inout) :: parallel
    real(dp), intent(in) :: f_distrib(parallel%mydim)

    ! communication variables
    integer node, nelems
#ifdef MPI
    integer status(MPI_STATUS_SIZE), mpinfo
#endif
#ifdef AJB_DEBUG
!      write(9,*) ' inside collect_function now I'
#endif

    if (parallel%iammaster) parallel%ftmp(:) = zero
    do node = 1, parallel%group_size
       nelems = parallel%irows(node) - parallel%irows(node-1)
       if (.not. parallel%iamgmaster) then
#ifdef AJB_DEBUG
!      write(9,*) ' inside collect_function now II'
#endif
          if (node == parallel%group_iam + 1) then
#ifdef MPI
             call MPI_SEND(f_distrib,nelems,MPI_DOUBLE_PRECISION, &
                  parallel%group_master,node,parallel%group_comm,mpinfo)
#endif
          endif
       else
          if (node == parallel%group_iam + 1) then
             call dcopy(nelems,f_distrib,1, &
                  parallel%ftmp(parallel%irows(node-1)),1)
          else
#ifdef MPI
             call MPI_RECV(parallel%ftmp(parallel%irows(node-1)), &
                  nelems,MPI_DOUBLE_PRECISION, &
                  node-1,node,parallel%group_comm,status,mpinfo)
#endif
          endif
       endif
#ifdef MPI
#ifdef AJB_DEBUG
      !write(9,*) ' collect function: hit the barrier'
#endif
       call MPI_Barrier(parallel%group_comm,mpinfo)
#endif
    enddo

  end subroutine collect_function

!=================================================================================
  subroutine collect_zfunction(parallel,f_distrib)
    !
    ! Collects array parallel%zftmp to master PE. The
    ! distribution follows grid partition defined in parallel%irows.
    ! Input dat is in array f_distrib
    !
#ifdef MPI
    use mpi
#endif
    implicit none
    type (parallel_data), intent (inout) :: parallel
    complex(dpc), intent(in) :: f_distrib(parallel%mydim*parallel%mxwd)

    ! counters
    integer :: ii, istart
    ! communication variables
    integer :: node, nelems, mpinfo
    integer :: offset
#ifdef MPI
    integer status(MPI_STATUS_SIZE)
#endif

    if (parallel%iammaster) parallel%zftmp(:) = zzero
    do ii = 1, parallel%mxwd
       istart = 1 + (ii - 1)*parallel%mydim
       do node = 1, parallel%group_size
          nelems = parallel%irows(node) - parallel%irows(node-1)
          offset = parallel%irows(node-1) + (ii - 1)*parallel%nwedge
          if (.not. parallel%iamgmaster) then
             if (node == parallel%group_iam + 1) then
#ifdef MPI
                call MPI_SEND(f_distrib(istart),nelems,MPI_DOUBLE_COMPLEX, &
                     parallel%group_master,node,parallel%group_comm,mpinfo)
#endif
             endif
          else
             if (node == parallel%group_iam + 1) then
                call zcopy(nelems,f_distrib(istart),1,parallel%zftmp(offset),1)
             else
#ifdef MPI
                call MPI_RECV(parallel%zftmp(offset), &
                     nelems,MPI_DOUBLE_COMPLEX, &
                     node-1,node,parallel%group_comm,status,mpinfo)
#endif
             endif
          endif
#ifdef MPI
          call MPI_Barrier(parallel%group_comm,mpinfo)
#endif
       enddo
    enddo

  end subroutine collect_zfunction

!=================================================================================
  subroutine group_reduction(parallel,f_distr)
    !
    ! Performs a reduction of function f_distr distributed across
    ! groups. The reduction is done over all processors which share
    ! the same set of grid points.
    !
    ! AJB: oh well... maybe something better than recv send?
#ifdef MPI
    use mpi
#endif
    implicit none
    type (parallel_data), intent (in) :: parallel
    real(dp), intent(inout) :: f_distr(parallel%mydim)

    real(dp) :: f_tmp(parallel%mydim)
    integer :: masternode, inode, msgtype, mpinfo
#ifdef MPI
    integer status(MPI_STATUS_SIZE)
#endif

#ifdef MPI
    masternode = mod(parallel%iam,parallel%group_size) !AJB: ?
    if (parallel%iam == masternode) f_tmp = f_distr
    do inode = masternode, parallel%procs_num - 1, parallel%group_size
       if (inode == masternode) cycle
       msgtype = inode
       if (parallel%iam == masternode) then
          call MPI_RECV(f_distr,parallel%mydim,MPI_DOUBLE_PRECISION, &
               inode,msgtype,parallel%comm,status,mpinfo)
          f_tmp = f_tmp + f_distr
       else
          if (parallel%iam == inode) call MPI_SEND(f_distr,parallel%mydim, &
               MPI_DOUBLE_PRECISION,masternode,msgtype,parallel%comm,mpinfo)
       endif
    enddo

    call MPI_Barrier(parallel%comm,mpinfo)

    do inode = masternode, parallel%procs_num - 1, parallel%group_size
       if (inode == masternode) cycle
       msgtype = inode + parallel%procs_num
       if (parallel%iam == masternode) then
          call MPI_SEND(f_tmp,parallel%mydim,MPI_DOUBLE_PRECISION, &
               inode,msgtype,parallel%comm,mpinfo)
       else
          if (parallel%iam == inode) &
          call MPI_RECV(f_tmp,parallel%mydim,MPI_DOUBLE_PRECISION, &
               masternode,msgtype,parallel%comm,status,mpinfo)
       endif
    enddo
    f_distr = f_tmp
    call MPI_BARRIER(parallel%comm,mpinfo)
#endif

  end subroutine group_reduction

!=================================================================================
end module parallel_data_module

! ==================== parallel data structure ==================
!
! This is a superstructure that holds variables used in the
! matvec subroutine
!
module parsec_global_data
  !
  ! type declarations
  !
  use constants
  use potential_module
  use non_local_psp_module
  use eigen_solver_module
  use parallel_data_module
  implicit none
  !
  ! define the actural data being used
  !
  ! potential related data
  type (potential) :: pot 

  ! non local pseudopotential related data
  type (nonloc_pseudo_potential) :: nloc_p_pot 

  ! on-site Coulomb interaction related data
  type (nonloc_pseudo_potential) :: u_pot

  ! solver related data
  type (eigen_solver) :: solver

  ! parallel computation related data
  type (parallel_data) :: parallel

end module parsec_global_data

! ==================== matvec structure ========================
module matvec_module
  use constants
  implicit none
  integer, dimension(:), pointer::natmi
  real(dp), dimension(:), pointer :: dwvec
  complex(dpc), dimension(:), pointer :: zwvec
  integer :: type_num
  ! flag: true if Zmatvec2 is to be used, false otherwise
  ! what is matvec2??!
  !logical :: use_matvec2

contains

  subroutine matvec_init(clust, ndim)

    use constants
    use cluster_module
    implicit none

    type (cluster), intent(in) :: clust
    integer, intent(in) :: ndim
    !logical, intent(in) :: flag_in

    integer alcstat

    !use_matvec2 = flag_in
    type_num = clust%type_num

    call matvec_destroy()

    allocate(natmi(clust%type_num))
    natmi(:) = clust%natmi(:)
    !okay, what? why both
    write(9,*) "allocating both dwvec and zwvec of size (why!):",ndim+1
    allocate(dwvec(ndim+1), stat=alcstat)
    call alccheck('dwvec',ndim+1,alcstat)
    allocate(zwvec(ndim+1), stat=alcstat)
    call alccheck('zwvec', ndim+1,alcstat)
  end subroutine matvec_init

  subroutine matvec_destroy()
    implicit none
    if (associated(natmi)) deallocate(natmi)
    if (associated(dwvec)) deallocate(dwvec)
    if (associated(zwvec)) deallocate(zwvec)
  end subroutine matvec_destroy

  subroutine init_matvec()
    nullify(natmi)
    nullify(dwvec)
    nullify(zwvec)
  end subroutine init_matvec
  
end module matvec_module

! ==================== band structure ========================
module bandstruc_module
  use constants
  implicit none

  type bandline
    
    ! line number
    integer :: line_num
    ! coordinates of start and end points of each line in direct coordinates
    real(dp), dimension(3)  :: start , end
    ! coordinates of start and end points of each line in
    ! cartesian coordinates of inverse Bohr radius
    real(dp), dimension(3)  :: startb , endb
    ! line name
    character(len=30) :: line_name
    ! line length 
    real(dp) :: length
    ! number of k-points along the line
    integer :: nkpt
    ! kpoint list (Cartesian coordinates, units of inverse Bohr radius)
    real(dp), dimension (:,:), pointer :: kpts
    ! eigenvlues by kpoint number, spin index, state number
    real(dp), dimension (:,:,:), pointer :: eigs

  end type bandline

  type bandstruc

    !flag for whether band structure calculation is on
    logical :: bands_on
    !number of points for the shortest line for band structure calculation
    integer :: npoints
    !defines the resolution of the band structure
    integer :: nlines
    type(bandline), dimension(:), pointer :: blines    

  end type bandstruc

end module bandstruc_module

!!!!!!!!!!!!! Manish
! This data structure are from PARATEC. They are used
! here with heavy modification. 
!
!     ========== gspace structure ==========================
!

module gspace_module  

  use constants
  implicit none  
  !
  !     first some public data which is global to ALL GSPACES
  !

  type gspace  
     !
     !     ------ true information --------------
     !
     real(dp) :: gmax             ! cutoff length
     real(dp) :: rk(3)            ! the k-vector by which it is shifted
     integer :: length            ! number of g-space vectors 
     integer :: nstar             ! total number of stars in gspace
     integer :: ig0               ! index of the g=0. If not on this proc, = -1
     logical :: imap              ! if indexg mapping is available
     logical :: igvec             ! if gvectors are available
     integer :: kmax(3)
     integer, pointer :: &
          gvec(:,:)         ! the gspace vectors themselves
     integer, pointer :: &
          indexg(:,:,:)         ! Given the gvec, this array gives the 
                              ! index in the gspace vector array of length
     !
     !     ------ FFT information ---------------
     !
     integer :: fftsize(3)        ! FFT grid boundaries for gspace
     character(len=16) :: name    ! a string identifying the gspace

  end type gspace

contains  

  subroutine nullify_gspace(gs)  

    implicit none
    type(gspace), intent(inout) :: gs  

    nullify(gs%gvec)  
    nullify(gs%indexg)  

  end subroutine nullify_gspace

  subroutine destroy_gspace(gs)  

    implicit none
    type(gspace), intent(inout) :: gs  

    if (associated(gs%gvec)) deallocate(gs%gvec)  
    if (associated(gs%indexg)) deallocate(gs%indexg)  

    call nullify_gspace(gs)

  end subroutine destroy_gspace

end module gspace_module

module nscf_module
  use constants
  implicit none

  type nscf_data
    
    !flag for whether non-selfconsistent calculation is on
    logical :: nscf_on
    ! Number of eigenstates to be calculated
    integer :: nstate
    ! number of kpoints
    integer :: nkpt
    ! kgrid
    integer :: kgrid(3)
    ! shift for the kgrid
    real(dp) :: kshift(3)
    ! kpoint list (Cartesian coordinates, units of inverse Bohr radius)
    real(dp), dimension (:,:), pointer :: kpts
    ! kpoint weights
    real(dp), dimension (:), pointer :: kpwt

  end type nscf_data
contains
  subroutine init_nscf(nscf)
  use constants
  type(nscf_data), intent(inout) :: nscf

  nscf%nscf_on = .false.
  nscf%nstate = -1
  nscf%nkpt   = -1
  nscf%kgrid = 0
  nscf%kshift = 0

  nullify(nscf%kpts)
  nullify(nscf%kpwt)
  end subroutine init_nscf

  subroutine destroy_nscf(nscf)
  use constants
  type(nscf_data), intent(inout) :: nscf

  if (associated(nscf%kpts)) deallocate(nscf%kpts)
  if (associated(nscf%kpwt)) deallocate(nscf%kpwt)

  end subroutine destroy_nscf

end module nscf_module

! =================== End of structures ======================
! ===============================================================
