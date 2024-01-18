!===============================================================
!
!  Copyright (C) 2005 Finite Difference Research Group
!  This file is part of parsec, http://www.ices.utexas.edu/parsec/
!
!  Version of ionpot for two-dimensional (slab) systems. Sets up the
!  ionic potential with two parts:
!     vion_local : ionic potential created by the local pseudopotential
!        and by a set of Gaussian charges centered at the atom sites; and
!     v_gauss    : negative of the potential created by the set of
!        Gaussian charges.
!
!  By construction, vion_local is short-ranged because the potential of
!  Gaussian charges compensates the local pseudopotential at large distances
!  (both of them are Coulombic). We calculated this part by doing straight
!  summation over atom sites on the unit cell and a finite number of cell
!  images.
!  The v_gauss potential is calculated by setting up the total charge
!  and solving Poisson's equation for it.
!
!  Some definitions:
!
!  rho_gauss (r) = Sum_atoms { 1/( sqrt(pi)*sigma )^3 * Z_ion *
!                              exp( -(rr/sigma)^2 ) }
!
!  V_gauss (r) = Sum_atoms { e^2 * Z_ion * ( 1 - erfc(rr/sigma) )/rr }
!
!  erfc(x) = complementary error function.
!  rr = distance between point "r" and each atom in the sum above.
!
! -------------------------------------------------------------
! O. Sinai (Weizmann), 2015: Added section that adds effects
! of fixed charged sheet (electrostatic) to Vion. Block marked
! below with CHARGED-SHEET.
!
!---------------------------------------------------------------
subroutine ionpot_slab(clust,grid,pot,p_pot,pbc,rsymm,parallel,lpole)

  use constants
  use cluster_module
  use grid_module
  use potential_module
  use pseudo_potential_module
  use pbc_module
  use symmetry_module
  use parallel_data_module
  implicit none

  !
  !  Input/Output variables:
  !
  !  the cluster
  type (cluster), intent(in) :: clust
  !  grid related data
  type (grid_data), intent(in) :: grid
  !  potential related data
  type (potential), intent(inout) :: pot
  !  pseudo_potential related data
  type (pseudo_potential), intent(in) :: p_pot
  !  periodic boundary conditions data
  type (pbc_data), intent(in) :: pbc
  !  symmetry operations in reduced group:
  type (symmetry), intent(in) :: rsymm
  !  parallel computation related data
  type (parallel_data), intent(in) :: parallel
  !  order of multipole expansion
  integer, intent(in) :: lpole

  !
  !  Work variables:
  !
  real(dp), dimension(1) :: psumvec

  !  Number of cell images used to calculate vion_local, on each (u,v) axis
  integer nreplicau, nreplicav
  
  ! Gaussian charge density related variables
  real (dp) :: &
       sigma, &           ! the Gaussian's standart deviation
       charge, &          ! the charge at each atom, the total area of the Gaussian
       fdens, &           ! the charge density
       netchrg            ! total electric charge in the unit cell

  ! Indexing and possition variables
  real(dp) :: &
       rindex, &          ! potential radius
       delr, &            ! potential deriviative
       rr(3), &           ! non-cartesian grid coordinates
       rr1(3), &          ! cartesian grid coordinates
       dist, &            ! the distance from each grid point
       distsigma          ! = dist/sigma


  ! Potential calculations variables
  real(dp) :: &
       vtemp, &           ! the contribution of each grid point to the local potential
       erfzero            ! the value of the (1-erfc(r)) deriviative at zero

  ! Potential and Charge Density arrays
  real(dp), allocatable :: &
       vion_local(:), &   ! Local pseudopotential
       rho_gauss(:), &    ! Charge density of the Gaussians
       vect(:), &         ! Work array for Poisson's solver
       brho(:)            ! Boundary charge density

  ! Counters
  integer ja, ity, iat, igrid, jgrid, jok, icell, jcell, i, j, ii, ics
  ! Distances and charges for charged-sheet part
  real(dp) :: dzabs, dzz, qsarea

  !  Allocation check
  integer alcstat

  !  External function
  real(dp), intrinsic :: erf

  !  Tolerance in compensated ionic charge
  real(dp), parameter :: tol_c = 1.d-2


  !---------------------------------------------------------------
  !
  !  Start by calculating vion_local and setting the total Gaussian charge.
  !

  ! Allocate and Initiate vectors

  allocate(vion_local(parallel%mydim),stat=alcstat)
  call alccheck('vion_local',parallel%mydim, alcstat)
  vion_local(:) = zero
  allocate(rho_gauss(parallel%ldn),stat=alcstat)
  call alccheck('rho_gauss',parallel%ldn, alcstat)
  rho_gauss(:) = zero


  ! Define the deriviative of (1-erfc(r)) at r-->0
  ! for future use in L'hospital's law.
  erfzero = two/sqrt(pi)

  
  ! Initaiate variables
  ja = 0
  netchrg = zero


  ! Go over all the atoms, in a large radius around the center cell, 
  ! and calculate their contribution to the local potential.
  ! Since we are adding a fake gaussian density, the sum of the two
  ! potentials should converge for a far enough distance.

  !  Go over all atom types

  do ity = 1, clust%type_num
     !  Go over all atoms, superimpose the contribution to the local
     !  ionic potential from each atom over each point of the grid.
     do iat = 1, clust%natmi(ity)
        ja = ja + 1

        !  Choose Gaussian distributions with dispersion proportional to
        !  the pseudoatom cut-off radius.

        sigma = sigmafac * p_pot%rcore(ity)
        charge = two * p_pot%zion(ity)
        netchrg = netchrg + p_pot%zion(ity)
        fdens = p_pot%zion(ity) / sigma / sigma / sigma / sqrt(pi) / pi

        !  The potential of Gaussian distributions cancels the long-range
        !  local potential very fast. We do not need to sum over more than
        !  2 or 3 replicas in each direction (u,v) around the periodic cell.
        !  To be on the safer side, we shall take up to 6 sigmas around.

        nreplicau = nint( six * sigma / pbc%box_size(1) ) + 1
        nreplicav = nint( six * sigma / pbc%box_size(2) ) + 1
 
        do icell = -nreplicau, nreplicau
           do jcell = -nreplicav, nreplicav

              do jgrid = 1, parallel%mydim
                 igrid = jgrid + parallel%irows(parallel%group_iam) - 1

                 rr(1) = (grid%shift(1) + grid%kx(igrid))*grid%step(1) + &
                      real(icell,dp)*pbc%box_size(1)
                 rr(2) = (grid%shift(2) + grid%ky(igrid))*grid%step(2) + &
                      real(jcell,dp)*pbc%box_size(2)
                 rr(3) = (grid%shift(3) + grid%kz(igrid))*grid%step(3)

                 ! Convert into cartesian coordinates
                 call matvec3('N',pbc%avec_norm,rr,rr1)

                 ! Calculate the vector from the grid point to the atom
                 rr1(1) = rr1(1) - clust%xatm(ja)
                 rr1(2) = rr1(2) - clust%yatm(ja)
                 rr1(3) = rr1(3) - clust%zatm(ja)

                 dist = sqrt(dot_product(rr1,rr1))

                 !
                 !  If dist is larger than the largest value of r given in the
                 !  pseudopotential file, rs(i,ity), use -zion/r as the ionic
                 !  potential, vtemp. The factor of two in zion is conversion from
                 !  Hartree to Rydberg.
                 if (dist >= p_pot%rs(p_pot%ns(ity)-1,ity)) then

                    vtemp = -two*p_pot%zion(ity)/dist

                 else
                    !
                    !  If dist is smaller than the largest value of r, find the
                    !  index of dist in the pseudopotential by inverting the
                    !  logarithmic pseudopotential grid, then interpolate the
                    !  LOCAL pseudopotential to determine the potential and charge
                    !  density.
                    !
                    !  rindex = 1/B(ity)*log((C(ity)+dist)/A(ity))
                    !  A = p_pot%par_a
                    !  B = p_pot%par_b
                    !  C = p_pot%par_c

                    rindex = one/p_pot%par_b(ity)* &
                         log((p_pot%par_c(ity)+dist)/p_pot%par_a(ity))
                    jok  = idint(rindex) + 1
                    delr = (dist-p_pot%rs(jok,ity))/ &
                         (p_pot%rs(jok+1,ity)-p_pot%rs(jok,ity))

                    !  Place pseudopotential on the 3d uniform grid by
                    !  interpolating the 1d pseudopotential given on the
                    !  logarithmic grid.
                    vtemp = p_pot%vion(jok,ity)*p_pot%rs(jok,ity)/dist + &
                         delr*(p_pot%vion(jok+1,ity)*p_pot%rs(jok+1,ity)- &
                         p_pot%vion(jok,ity)*p_pot%rs(jok,ity))/dist
                    if (jok == 1) vtemp = p_pot%vion(2,ity)
                 endif

                 vion_local(jgrid) = vion_local(jgrid) + vtemp
                 distsigma = dist/sigma

                 !  Add the potential of a Gaussian distribution.
                 !  if dist == 0, use l'hospital law to calculate the limit 
                 !  of (one - erfc(distsigma)) * charge / dist when distsigma --> 0

                 if (dist == 0) then
                    vion_local(jgrid) = vion_local(jgrid) + (erfzero * charge) /sigma
                 else
                    vion_local(jgrid) = vion_local(jgrid) + &
                      ( erf(distsigma) ) / dist * charge
                 endif

                 !  Calculate the Gaussian charge at this grid point
                 distsigma = distsigma*distsigma
                 rho_gauss(jgrid) = rho_gauss(jgrid) + fdens * exp(-distsigma)

              enddo               ! jgrid = 1, parallel%mydim
           enddo                  ! jcell = -nreplicav, nreplicav
        enddo                     ! icell = -nreplicau, nreplicau
     enddo                        ! iat = 1, clust%natmi(ity)
  enddo                           ! ity = 1, clust%type_num



  !  Check charge neutrality in the slab.

  charge = sum(rho_gauss(1:parallel%mydim)) * grid%hcub * real(rsymm%ntrans,dp)
  psumvec = charge
  call psum(psumvec,1,parallel%group_size,parallel%group_comm)
  charge = psumvec(1)

  if (parallel%iammaster) write(7,'(/,a,/,a,g11.3,a,g11.3,/)') &
       ' Charge neutrality check.',' net charge in unit cell [e] = ', &
       netchrg - charge,'   out of ',netchrg

  if ( abs(netchrg - charge) > tol_c .and. parallel%iammaster) &
       write(7,*) ' WARNING: net charge is not zero. The ', &
       ' compensated ionic charge is not neutral: ',netchrg, charge


  ! Remove the effect of the fake density added. 
  ! To calculate the potential the fake gaussian charges created,
  ! use the hartset_slab subroutine to get the boundary conditions,
  ! and solve Poisson's equation with those conditions to get the fake potential.

  allocate(vect(parallel%nwedge+1),stat=alcstat)
  call alccheck('vect',parallel%nwedge+1,alcstat)
  allocate(brho(parallel%ldn),stat=alcstat)
  call alccheck('brho',parallel%ldn,alcstat)
  brho(:) = zero

  call hartset_slab(grid,pbc,parallel,grid%norder,lpole, &
       grid%coe2,rho_gauss,brho)

  call hpotcg(parallel,brho,rho_gauss,grid%norder,grid%coe2, &
       grid%lap_dir_num,vect)
  deallocate(brho, vect)


  ! The difference between the two potentials - the local potential and the fake potential,
  ! is the Vion potential.

  do jgrid = 1, parallel%mydim
     pot%vion(jgrid) = vion_local(jgrid) - rho_gauss(jgrid)
  enddo

  ! vtemp=sum(pot%vion)/parallel%mydim ! AMIR - should be changed for MPI

  ! pot%vion=pot%vion-vtemp


  !=======================BEGIN CHARGED-SHEET=====================

  if (clust%has_charged_sheet) then
     ! Initialize ics - counter of all charged sheets in the system.
     ics = 0

     ! Go over all charged sheet types.
     do ity = 1, clust%number_charged_sheet_type
        ! Go over sheets within each type, superimpose the
        ! contribution to the local ionic potential from each sheet
        ! over each point of the grid
        do i = 1, clust%number_charged_sheets_per_type(ity)
           ics = ics + 1
           do jgrid = 1,parallel%mydim
              igrid = jgrid + parallel%irows(parallel%group_iam) - 1

              dzz = (grid%shift(3) + grid%kz(igrid))*grid%step(3) &
                   - clust%z_charged_sheets(ics)
              dzabs = abs(dzz)

              ! The potential is simply 2*pi*q_sheet*|dzz|, where q_sheet is the
              ! per-area charge of the sheet. The extra factor of two is
              ! conversion from hartrees to rydbergs.
              qsarea = clust%sheet_charge_of_type(ity)/pbc%a_surface_cell
              vtemp = two*pi*two*qsarea*dzabs
              pot%vion(jgrid) = pot%vion(jgrid) + vtemp

           enddo
        enddo
     enddo

  endif
  !===================END CHARGED-SHEET=========================


  deallocate(vion_local, rho_gauss)

end subroutine ionpot_slab
!===============================================================


