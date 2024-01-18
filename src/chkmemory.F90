!===============================================================
!
! Copyright (C) 2005 Finite Difference Research Group
! This file is part of parsec, http://www.ices.utexas.edu/parsec/
!
! Report memory requirements.
! WARNING: Future modifications of this code should include an
! update of the list of "big arrays" contained in this subroutine.
! Big arrays are defined as arrays with size proportional to the
! number of grid-points (or reciprocal space points, equivalently)
! or with size proportional to the square of number of computed
! states.
!
! Currently, these are the classes of big arrays:
! (1) global arrays, allocated by all PEs during most of the
!     calculation;
! (2) 0PE arrays, allocated only in single-PE mode (no MPI
!     interface), during most of the calculation;
! (3) M+0PE arrays, allocated by master PE or allocated in 
!     single-PE mode, during most of the calculation;
! (4) temp M+0PE arrays, allocated by master PE or allocated in 
!     single-PE mode, but used only in short run-time subroutines;
! (5) temp arrays, allocated by all PEs,
!     but used only in short run-time subroutines;
! (7) temp S arrays,  allocated only in MPI environment and used
!     only in short run-time subroutines;
!
! For each class i of arrays, there is a counter memc(i) which
! stores the amount of bytes used to store those arrays. At the 
! end of the subroutine, memory counts are converted to MB and
! printed out.
!
! Author: M. Tiago (2004)
!
! AJB: Terribly out of date and inaccurate!
!---------------------------------------------------------------
subroutine chkmemory(elec_st,grid,pbc,mixer,solver,parallel)

  use constants
  use electronic_struct_module
  use grid_module
  use pbc_module
  use mixer_module
  use eigen_solver_module
  use parallel_data_module

  implicit none
  !
  ! Input/Output variables:
  !
  ! electronic structure
  type(electronic_struct), intent(in) :: elec_st
  ! grid related data
  type(grid_data), intent(in) :: grid
  ! periodic boundary conditions data
  type(pbc_data), intent (in) :: pbc
  ! mixer related data
  type(mixer_data), intent (in) :: mixer
  ! solver related data
  type(eigen_solver), intent (in) :: solver
  ! parallel computation related data
  type (parallel_data), intent(in) :: parallel
  !
  ! Work variables:
  !
  ! counters
  integer :: ii
  real(dp) :: tmp,memc(20)
  !
  ! aliases for parameters contained in the structures:
  !
  ! parallel%procs_num : number of computing processors
  integer :: npe
  ! grid%ndim : number of points in real-space grid
  integer :: ndim
  ! grid%nwedge : number of points in irreducible wedge
  integer :: nwedge
  ! ldn : estimated number of grid points owned by each computing processor
  integer :: ldn
  ! grid%nxmax    : number of grid points
  ! grid%nymax    : along each direction x,y,z
  ! grid%nzmax    : 
  integer :: nxyzmax
  ! grid%norder : number of neighbors in numerical derivative
  integer :: norder
  ! elec_st%nstate : number of computed Kohn-Sham eigenstates
  integer :: nstate
  ! elec_st%nspin : number of spins
  integer :: nspin
  ! kss0 : iteration parameter in hpotcg
  integer :: kss0
  ! number of kpoints
  integer :: kpnum
  !---------------------------------------------------------------

  npe = parallel%group_size
  ndim = grid%ndim
  nwedge = grid%nwedge
#ifdef MPI
  ldn = nwedge/npe + npe
#else
  ldn = nwedge
#endif
  nxyzmax = (2*grid%nxmax+1)*(2*grid%nymax+1)*(2*grid%nzmax+1)
  norder = grid%norder
  nstate = elec_st%nstate
  nspin = elec_st%nspin
  kpnum = max(elec_st%nkpt,1)
  memc = 0
  !
  ! memc(1) = size of Global arrays
  !
  ! elec_st%eig%wf, double precision : wave function arrays
  ! WARNING: very approximate estimate for memory (this is a lower
  ! bound only, ignores the additional states from solver%nadd)
  ii = max(elec_st%nstate*kpnum,2*elec_st%nstate*kpnum/elec_st%nrep)
  if (elec_st%cplx) then
     ii = 2*ii + solver%nadd + solver%winsize
  else
     ii = ii + solver%nadd + solver%winsize
  endif
  memc(1) = memc(1) + eight*ldn*nspin*ii

  if (solver%name == ARPACK) then
     !
     ! Vbasis, workd, resid, %resid, workl, double precision : work arrays
     ! used in arpack
     ! WARNING: very approximate estimate for memory (now, this is
     ! overestimate; it assumes that nstate levels are computed per reprs.)
     memc(1) = memc(1) + 8*ldn*(3*nstate + 10)/elec_st%nrep
     memc(1) = memc(1) + 8*(2*nstate + 5)*(2*nstate + 13)
  elseif (solver%name == TRLANC) then
     !
     ! Vbasis, double precision : work array for thick-restart Lanczos
     ! WARNING: very approximate estimate for memory (now, this is
     ! overestimate; it assumes that nstate levels are computed per reprs.)
     memc(1) = memc(1) + 8*ldn*(2*nstate + 10)/elec_st%nrep
  elseif (solver%name == DIAGLA) then
     !
     ! ww, double precision : temporary array used in diagla
     memc(1) = memc(1) + 8*ldn*(solver%kss0 + 1)
  endif

  ! add the extra memory used in Chebyshev subspace
  if (solver%do_subsp ) then
     memc(1)=memc(1) + 8*ldn
     ! count the extra work arrays (can be overestimate)
     memc(1)=memc(1)+ldn*(nstate/elec_st%nrep+solver%nadd)+ldn*30
  endif
  !
  ! brho, double precision : temporary array used in diagla and hpotcg
  memc(1) = memc(1) + 8*ldn*(solver%kss0 + 1)
  !
  ! adiag, double precision : diagonal part of Hamiltonian
  memc(1) = memc(1) + 8*ldn
  !
  ! rho,vhart, double precision : electronic structure arrays
  memc(1) = memc(1) + 8*ldn*2*nspin
  !
  ! grid%fx,fy,fz, integer : grid layout (full grid)
  memc(1) = memc(1) + 4*ndim*3
  ! grid%rindex, rtrans : index of points and symmetry transformation
  ! needed to bring them to IW
  memc(1) = memc(1) + 4*ndim*2
  ! grid%kx, ky, kz, parallel%neibs, tneibs, senrows, pint, integer :
  ! indices of grid points
  memc(1) = memc(1) + 4*ldn*( 2 + 2*6*norder) + 4*3*nwedge
  !
  ! pot%vold,vnew,vhxcold,vxc,vhart,vion, double precision : 
  ! potentials in "pot" structure
  memc(1) = memc(1) + 8*ldn*nspin*4 + 8*ldn*2
  !
  ! vect : array used in parallel Poisson solver (hpotcg.F)
  memc(1) = memc(1) + 8*nwedge
  !
  ! elec_st%rho,rhoc, double precision : charge density
  memc(1) = memc(1) + 8*ldn*nspin*2
  !
  ! mixer%resid1,old1,old2,xin,xout, double precision :
  ! mixer arrays
  memc(1) = memc(1) + 8*ldn*5
  ! mixer%xinold,xoutold are also used by Anderson mixing
  if (mixer%name == ANDERSON) memc(1) = memc(1) + 8*ldn*2*nspin*mixer%memory
  !
  ! grid%indexw, integer : index of points in 3-d grid
  memc(1) = memc(1) + 4*nxyzmax
  !
  ! memc(2) = size of 0PE arrays
  ! 
  ! wkm wk2, double precision : work arrays used in hpotcg (Poisson solver)
  ! WARNING: the value below for kss0 must be consistent with hpotcg.F
  kss0 = 6
  memc(2) = memc(2) + 8*nwedge*( kss0 + 1 )
  !
  ! memc(3) = size of arrays stored by master and in single
  ! processor mode
  !
  ! grid%indexg, integer : index of points in 3-d grid
  memc(3) = memc(3) + 4*nxyzmax
  !
  ! memc(4) = size of temporary arrays stored by master and in
  ! single processor mode; not stored at all times, take just
  ! maximum value
  !
  ! dvhcc,vxcav, double precision : arrays for local forces (forcloc)
  if (.not. pbc%is_on) then
     tmp = 8*nwedge*2
     memc(4) = max(memc(4),tmp)
  endif
  !
  ! double precision, temporary arrays in restart.f
  tmp = 8*nwedge*10
  memc(4) = max(memc(4),tmp)
  !
  ! tmp2, integer : temporary array for neighbors, in setup.F
  tmp = 4*nwedge*norder*2*3*2
  memc(4) = max(memc(4),tmp)
  !
  ! ninv,rinv, integer : temporary array for irreducible wedge, in
  ! grid_partition
  tmp = 4*ndim*2
  memc(4) = max(memc(4),tmp)
  !
  ! expg, pbc%vscr*,vsci* and all the double precision/complex arrays used
  ! for FFT
  if ( pbc%is_on ) then
     tmp = 8*nxyzmax*4 + 8*pbc%ng + 8*pbc%nstar*4
     memc(4) = max(memc(4),tmp)
  endif
  !
  ! memc(5) = size of temporary arrays; not stored at all times,
  ! take just maximum value
  !
  ! df,u,v, double precision : mixer arrays used in Broyden mixing
  if (mixer%name == BROYDEN) then
     tmp = 8*ldn*nspin*(mixer%block*2 + 1)*2 + 8*ldn*nspin*mixer%memory
     memc(5) = max(memc(5),tmp)
  endif
  !
  ! delta1, delta2, double precision : mixer arrays used in Anderson
  if (mixer%name == ANDERSON) then
     tmp = 8*ldn*nspin*2
     memc(5) = max(memc(5),tmp)
  endif
  !
  ! gradients for exchange-correlation (exc_nspn.f exc_spn.f)
  memc(5) = memc(5) + 8*ldn*18
#ifdef MPI
  !
  ! memc(7) = temporary arrays stored by computing PEs but only
  ! with MPI interface
  !
  ! wk,wk2,jake,pint, double precision : temporary arrays used in
  ! parallel Poisson solver (hpotcg.F)
  tmp = 8*ldn*( kss0 + 2 + 3*2*norder )
  !
  ! senrows, integer : array for MPI communication in hpotcg.F
  tmp = tmp + 4*nwedge
  memc(7) = max(memc(7),tmp)
  !
  ! senrows/recvrow temporary array for setup
  tmp = 8*nwedge*3 + 4*nwedge*3
  memc(7) = max(memc(7),tmp)
#endif
  !
  ! convert all memory counts to megabytes
  !
  memc = memc/1024/1024
  !
  write(7,*)
  write(7,*) 'Memory usage:'
  write(7,*) '-------------'
#ifdef MPI
  write(7,20) ' estimated memory usage per PE    : ', memc(1)+memc(7),' MB'
  write(7,20) ' additional memory in master PE   : ', memc(3)+memc(4),' MB'
  write(7,20) ' maximum estimated memory per PE  : ', &
       memc(1)+memc(3)+memc(4)+memc(7),' MB'
#else 
  write(7,20) ' estimated memory usage (no MPI)  : ', sum(memc(1:4)),' MB'
#endif
  write(7,*)
  call myflush(7)

20 format(a,f8.2,a)

end subroutine chkmemory
!===============================================================
