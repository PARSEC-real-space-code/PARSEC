!===============================================================
!
!  Copyright (C) 2005 Finite Difference Research Group
!  This file is part of parsec, http://www.ices.utexas.edu/parsec/
!
!  Version of ionpot for one-dimensional (wire) systems. Sets up the
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
!  author: Murilo Tiago, UTexas, February 2007
!
!---------------------------------------------------------------
subroutine ionpot_wire(clust,grid,pot,p_pot,pbc,rsymm,parallel,lpole)

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
  !  number of cell images used to calculate vion_local
  integer nreplica
  !  temporary storage variables
  real(dp) :: dist, rindex, delr, vtemp, dvtemp, rr(3), rr1(3)
  !  temporary holders of Gaussian dispersion, charge, rr/sigma, density
  real(dp) :: sigma, charge, distsigma, fdens
  !  total electric charge in the unit cell
  real(dp) :: netchrg
  !  counters
  integer ja, ity, iat, igrid, jgrid, jok, icell
  !  allocation check
  integer alcstat
  !  local pseudopotential + potential of the Gaussian charge
  real(dp), allocatable :: vion_local(:)
  !  total Gaussian charge
  real(dp), allocatable :: rho_gauss(:)
  !  work arrays for the Poisson solver
  real(dp), allocatable :: vect(:),brho(:)
  !  intrinsic function
  real(dp), intrinsic :: erf
  !  tolerance in compensated ionic charge
  real(dp), parameter :: tol_c = 1.d-2

  ! work variable
  real(dp), dimension(1) :: chargevec(1)

  !---------------------------------------------------------------
  !
  !  Start by calculating vion_local and setting the total Gaussian charge.
  !
  allocate(vion_local(parallel%mydim),stat=alcstat)
  call alccheck('vion_local',parallel%mydim, alcstat)
  vion_local(:) = zero
  allocate(rho_gauss(parallel%ldn),stat=alcstat)
  call alccheck('rho_gauss',parallel%ldn, alcstat)
  rho_gauss(:) = zero

  !  Initialize ja - counter of all atoms in the system.
  ja = 0
  !  Initialize cumulator for the total electric charge
  netchrg = zero
  !  Go over all atom types.

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
        !  2 or 3 replicas around the periodic cell.
        nreplica = nint( six * sigma / pbc%box_size(1) ) + 1
        
        do icell = -nreplica, nreplica

           do jgrid = 1, parallel%mydim
              igrid = jgrid + parallel%irows(parallel%group_iam) - 1

              rr(1) = (grid%shift(1) + grid%kx(igrid))*grid%step(1) - &
                   clust%xatm(ja) + real(icell,dp)*pbc%box_size(1)
              rr(2) = (grid%shift(2) + grid%ky(igrid))*grid%step(2) - &
                   clust%yatm(ja)
              rr(3) = (grid%shift(3) + grid%kz(igrid))*grid%step(3) - &
                   clust%zatm(ja)
              call matvec3('N',pbc%avec_norm,rr,rr1)
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
                 !  logarithmic grid. If "old style", interpolate the
                 !  pseudopotential directly.
                 vtemp = p_pot%vion(jok,ity)*p_pot%rs(jok,ity)/dist + &
                      delr*(p_pot%vion(jok+1,ity)*p_pot%rs(jok+1,ity)- &
                      p_pot%vion(jok,ity)*p_pot%rs(jok,ity))/dist
                 if (p_pot%uspline(ity)) then
                    call splint(p_pot%rs(-p_pot%norder:p_pot%ns(ity),ity), &
                                p_pot%vion(-p_pot%norder:p_pot%ns(ity),ity), &
                                p_pot%d2vion(-p_pot%norder:p_pot%ns(ity),ity), &
                                p_pot%ns(ity)+1+p_pot%norder,dist,vtemp,dvtemp)
                 endif
                 if (jok == 1) vtemp = p_pot%vion(2,ity)
              endif

              vion_local(jgrid) = vion_local(jgrid) + vtemp
              distsigma = dist/sigma

              !  Add the potential of a Gaussian distribution.
              if( dist .ne. 0 ) then
                  vion_local(jgrid) = vion_local(jgrid) + &
                       ( erf(distsigma) ) / dist * charge

              else
              ! potential limit z-->0 according to ionpot_slab / basic math
                  vion_local(jgrid) = vion_local(jgrid) + &
                       (two/sqrt(pi) * charge) / sigma
              end if
              !  Calculate the Gaussian charge at this grid point
              distsigma = distsigma*distsigma
              rho_gauss(jgrid) = rho_gauss(jgrid) + fdens * exp(-distsigma)
           enddo               ! jgrid = 1, parallel%mydim
        enddo                  ! icell = -nreplica, nreplica
     enddo                  ! iat = 1, clust%natmi(ity)
  enddo                     ! ity = 1, clust%type_num

  !  Check charge neutrality in the wire.
  charge = sum(rho_gauss(1:parallel%mydim)) * grid%hcub * real(rsymm%ntrans,dp)
  chargevec = charge
  call psum(chargevec,1,parallel%group_size,parallel%group_comm)
  charge = chargevec(1)
  if (parallel%iammaster) write(7,'(/,a,/,a,g11.3,a,g11.3,/)') &
       ' Charge neutrality check.',' net charge in unit cell [e] = ', &
       netchrg - charge,'   out of ',netchrg
  if ( abs(netchrg - charge) > tol_c .and. parallel%iammaster) &
       write(7,*) ' WARNING: net charge is not zero. The ', &
       ' compensated ionic charge is not neutral: ',netchrg, charge

  !  Solve Poisson's equation for the Gaussian charge distribution.
  allocate(vect(parallel%nwedge+1),stat=alcstat)
  call alccheck('vect',parallel%nwedge+1,alcstat)
  allocate(brho(parallel%ldn),stat=alcstat)
  call alccheck('brho',parallel%ldn,alcstat)
  brho(:) = zero

  call hartset_wire(grid,pbc,rsymm,parallel,grid%norder,lpole, &
       grid%coe2,rho_gauss,brho)

  call hpotcg(parallel,brho,rho_gauss,grid%norder,grid%coe2, &
       grid%lap_dir_num,vect)
  deallocate(brho, vect)

  !  Add vion_local to the potential of the Gaussian distribution.
  do jgrid = 1, parallel%mydim
     pot%vion(jgrid) = vion_local(jgrid) - rho_gauss(jgrid)
  enddo

  deallocate(vion_local, rho_gauss)

end subroutine ionpot_wire
!===============================================================
