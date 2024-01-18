!===============================================================
!
!  Copyright (C) 2005 Finite Difference Research Group
!  This file is part of parsec, http://www.ices.utexas.edu/parsec/
!
!  This subroutine calculates the ionic force and total nuclear energy
!  due to the local pseudopotential. Boundary conditions are wire (1-D).
!
!  As in the calculation of the ionic potential, there are two contributions:
!  one arising from the sum of local pseudopotential and Gaussian potential,
!  and another which is the negative of the Gaussian potential.
!
!  author: Murilo Tiago, UTexas, February 2007
!
!---------------------------------------------------------------
subroutine forcloc_slab(clust,grid,pbc,pot,p_pot,rsymm,parallel,rho,ipr)

  use constants
  use cluster_module
  use grid_module
  use potential_module
  use pseudo_potential_module
  use symmetry_module
  use parallel_data_module
  use pbc_module
  implicit none
  !
  !  Input/Output variables:
  !
  !  the cluster
  type (cluster), intent(inout) :: clust
  !  grid related data
  type(grid_data), intent(in) :: grid
  !  periodic boundary conditions data
  type (pbc_data), intent(in) :: pbc
  !  potential related data
  type (potential), intent(in) :: pot
  !  pseudo_potential related data
  type (pseudo_potential), intent(in) :: p_pot
  !  symmetry operations in reduced group:
  type (symmetry), intent(in) :: rsymm
  !  parallel computation related data
  type (parallel_data), intent(in) :: parallel

  !  charge density on the 3-d grid (total charge only, 
  !  even if spin-polarized)
  real(dp), intent(in) :: rho(parallel%mydim)
  !  print flag
  integer, intent(in) :: ipr
  !
  !  Work variables:
  !
  !  hcub = h**3
  real(dp) :: hcub

  !  temporary variables defined at the beginning of the loop over ity
  integer jcore, npoint
  real(dp) ::  zcore, rmin, rlarge, ainv, binv, c_a
  !  temporary holders of pseudopotential, atomic charge, 
  !  core charge derivatives
  real(dp) :: dv0, dv, rhoc, drhoc, rhoat, drhoat
  !  number of cell images used to calculate vion_local
  integer nreplicax, nreplicay
  !  temporary holders of Gaussian dispersion, charge, rr/sigma, density
  real(dp) :: sigma, distsigma, fdens
  !  distance between grid point and atom, its square, and
  !  projections
  real(dp) :: r1, r1sq, rr1(3), invr1
  !  index of radial grid, distance to previous grid point
  integer ind
  real(dp) :: delinv

  !  force variables:
  !  Term explicitly involving the local pseudopotential, and
  !  components.
  real(dp) :: fl, f_l(3)
  !  Term involving the Gaussian distribution, and
  !  components.
  real(dp) :: f_g(3)
  !  term resulting from core-correction, and components
  real(dp) :: fcore, f_core(3)
  !  term resulting from non-consistency correction, and components
  real(dp) :: fc, f_c(3)
  !  counters
  integer ii, ja, ity, iat, igrid, jgrid, itrans, icellx, icelly
  real(dp) :: rw(3), rrw(3), rrw1(3)
  !  allocation check
  integer alcstat

  !  work arrays:
  !  Difference of Vhxc from current iteration and the previous one
  !  (average of up and down components, if spin-polarized computation)
  real(dp), allocatable :: dvhxc(:)
  !  Average of spin-up and spin-down vxc. equal to just vxc in the 
  !  absence of spin-polarization; needed only if core-correction is
  !  used.
  real(dp), allocatable :: vxcav(:)
  !  Gradient of the Hartree potential.
  real(dp), allocatable :: grad_vh(:,:)
  !  Work array.
  real(dp), allocatable :: wvec(:)
  !  Characters of the gradient operator
  integer :: chi_grad(rsymm%ntrans,3)
  !  constant
  real(dp) :: twosqrtpi
  !
  ! External functions:
  !
  real(dp), external :: erfc

  !---------------------------------------------------------------

  hcub = grid%hcub

  allocate(dvhxc(parallel%mydim),stat=alcstat)
  call alccheck('dvhxc',parallel%mydim,alcstat)
  allocate(vxcav(parallel%mydim),stat=alcstat)
  call alccheck('vxcav',parallel%mydim,alcstat)

  !  Prepare difference of Hartree-exchange-correlation with respect
  !  to previous iteration (for accelerating force convergence).
  do ii = 1, parallel%mydim
     dvhxc(ii) = (pot%vhart(ii) + pot%vxc(ii,1)) - pot%vhxcold(ii,1)
  enddo
  !  If spin-polarized, average the spin-up and spin-down dvhxc.
  if (pot%nspin == 2) then
     do ii = 1, parallel%mydim
        dvhxc(ii)=dvhxc(ii)+(pot%vhart(ii)+pot%vxc(ii,2)) - &
             pot%vhxcold(ii,2)
        dvhxc(ii)=half*dvhxc(ii)
     enddo
  endif

  !  Prepare vxc that is average of up and down vxc.
  if (pot%nspin == 1) then
     do ii = 1, parallel%mydim
        vxcav(ii) = pot%vxc(ii,1)
     enddo
  else
     do ii = 1, parallel%mydim
        vxcav(ii) = half*(pot%vxc(ii,1) + pot%vxc(ii,2))
     enddo
  endif

  !  Define the characters of the gradient operator
  call grad_char(grid,rsymm,parallel,chi_grad)

  !  Calculate the gradient of the Hartree potential. It will be used to
  !  compute the force term from the Gaussian charge.
  allocate(grad_vh(3,parallel%mydim),stat=alcstat)
  call alccheck('grad_vh',3*parallel%mydim,alcstat)
  allocate(wvec(parallel%nwedge+1),stat=alcstat)
  call alccheck('wvec',parallel%nwedge+1,alcstat)
  call gradmvs(grid,parallel,pot%vhart,grad_vh,grid%coe1,grid%norder,wvec)
  deallocate(wvec)

  !  Header, if printing.
  if (parallel%iammaster .and. ipr >= 2) then
     write(7,*) 'Format is: atom #, followed by x,y,z components of'
     write(7,*) 'ion, self-consistency, core-correction forces, and'
     write(7,*) 'Gaussian forces.'
  endif

  !  Initialize atom counter, ja.
  ja = 0
  !  For all atom types...
  do ity = 1, clust%type_num
     !  Define temporary variables.
     jcore = p_pot%icore(ity)

     npoint = p_pot%ns(ity)
     !  Factor "two" comes from unit change (hartrees -> rydbergs).
     zcore  = -two*p_pot%zion(ity)
     rmin = p_pot%rs(2,ity) 
     rlarge = p_pot%rs(npoint-6,ity)

     !  Choose Gaussian distributions with dispersion proportional to
     !  the pseudoatom cut-off radius.
     sigma = sigmafac * p_pot%rcore(ity)
     fdens = p_pot%zion(ity) / sigma / sigma / sigma / sqrt(pi) / pi
     twosqrtpi = two / sqrt(pi) / sigma
     !  The potential of Gaussian distributions cancels the long-range
     !  local potential very fast. We do not need to sum over more than
     !  2 or 3 replicas around the periodic cell.
     nreplicax = nint( six * sigma / pbc%box_size(1) ) + 1
     nreplicay = nint( six * sigma / pbc%box_size(2) ) + 1
     ainv = one/p_pot%par_a(ity)
     binv = one/p_pot%par_b(ity)
     c_a = p_pot%par_c(ity) / p_pot%par_a(ity)
     !  For all atoms within each type...
     do iat   = 1, clust%natmi(ity)
        !  Update atom counter.
        ja = ja + 1
        !  Initialize all force accumulators (and core charge derivative)
        !  to zero.
        f_l(:) = zero
        f_g(:) = zero
        f_c(:) = zero
        f_core(:) = zero
        do icellx = -nreplicax, nreplicax
         do icelly = -nreplicay, nreplicay

           !  For all grid points...
           do jgrid = 1, parallel%mydim
              igrid = jgrid + parallel%irows(parallel%group_iam) - 1
              rw(1) = (grid%shift(1) + grid%kx(igrid)) * grid%step(1) &
                   + real(icellx,dp)*pbc%box_size(1)
              rw(2) = (grid%shift(2) + grid%ky(igrid)) * grid%step(2)  &
                   + real(icelly,dp)*pbc%box_size(2)
              rw(3) = (grid%shift(3) + grid%kz(igrid)) * grid%step(3)
              call matvec3('N',pbc%avec_norm,rw,rrw1)
              do itrans = 1, rsymm%ntrans
                 call symop(rsymm,itrans,rrw1,rrw)

                 !  Calculate distance from grid point to atom.
                 rr1(1) = rrw(1) - clust%xatm(ja)
                 rr1(2) = rrw(2) - clust%yatm(ja)
                 rr1(3) = rrw(3) - clust%zatm(ja)

                 r1sq = dot_product(rr1,rr1)
                 r1 = sqrt(r1sq)

                 !  If distance to atom negligible, all forces are zero all 
                 !  charge and potential derivatives vanish, so skip.
                 if(r1 < rmin) cycle
                 invr1 = one/r1
                 !
                 !  If distance to atom larger than the largest value in 
                 !  pseudopotential file, zion/r is used as the potential. 
                 if (r1 >= rlarge) then
                    dv = -zcore/r1sq
                    drhoat = zero 
                 else
                    !  Define index variables for interpolation between radial
                    !  grid points.
                    ind = idint(binv*log(c_a+ainv*r1)) + 1
                    delinv=one/(p_pot%rs(ind+1,ity)-p_pot%rs(ind,ity))

                    dv = ( (p_pot%rs(ind+1,ity)-r1)*p_pot%dvion(ind,ity)* &
                         p_pot%rs(ind,ity)*p_pot%rs(ind,ity) + &
                         (r1-p_pot%rs(ind,ity))*p_pot%dvion(ind+1,ity) &
                         *p_pot%rs(ind+1,ity)*p_pot%rs(ind+1,ity))* &
                         delinv*invr1*invr1
                    drhoat = ((p_pot%rs(ind+1,ity)-r1) *  &
                         p_pot%drhodr(ind,ity) + (r1-p_pot%rs(ind,ity))* &
                         p_pot%drhodr(ind+1,ity)) * delinv
                    if (p_pot%uspline(ity)) then
                       call splint(p_pot%rs(-p_pot%norder:npoint,ity), &
                                   p_pot%vion(-p_pot%norder:npoint,ity), &
                                   p_pot%d2vion(-p_pot%norder:npoint,ity), &
                                   npoint+1+p_pot%norder,r1,dv0,dv)
                       call splint(p_pot%rs(-p_pot%norder:npoint,ity), &
                                   p_pot%drhodr(-p_pot%norder:npoint,ity), &
                                   p_pot%d2rhodr(-p_pot%norder:npoint,ity), &
                                   npoint+1+p_pot%norder,r1,rhoat,drhoat)
                    endif
                 endif

                 !
                 !  Add contribution from the Gaussian charge.
                 distsigma = r1/sigma
                 dv = dv + zcore * ( one - erfc(distsigma) ) / r1sq - &
                      zcore * twosqrtpi * exp(-distsigma*distsigma) / r1
                 !
                 !  Interpolate core charge...
                 if (jcore > 0) then
                    if (r1 >= rlarge) then
                       drhoc = zero
                    else
                       ind = idint(binv*log(c_a+ainv*r1)) + 1
                       delinv = one/(p_pot%rs(ind+1,ity)-p_pot%rs(ind,ity))

                       drhoc = ((p_pot%rs(ind+1,ity)-r1)* &
                            p_pot%ddenc(ind,ity) + (r1-p_pot%rs(ind,ity))* &
                            p_pot%ddenc(ind+1,ity)) * delinv
                       if (p_pot%uspline(ity)) then
                          call splint(p_pot%rs(-p_pot%norder:npoint,ity), &
                                      p_pot%denc(-p_pot%norder:npoint,ity), &
                                      p_pot%d2denc(-p_pot%norder:npoint,ity), &
                                      npoint+1+p_pot%norder,r1,rhoc,drhoc)
                       endif
                    endif

                 endif
                 !
                 !  Update force components according to the formulas given
                 !  in the beginning.
                 !
                 fl = dv*rho(jgrid)/r1
                 fc = drhoat*dvhxc(jgrid)/r1
                 f_l = f_l + fl*rr1
                 f_c = f_c + fc*rr1
                 if (jcore > 0) then
                    fcore = drhoc*vxcav(jgrid)/r1
                    f_core = f_core + fcore*rr1
                 endif

                 !  Add contribution from the negative of Gaussian charge.
                 !  We get this from the Hartree potential.
                 distsigma = distsigma*distsigma
                 do ii = 1, 3
                    f_g(ii) = f_g(ii) + grad_vh(ii,jgrid) * &
                         chi_grad(itrans,ii) * fdens * exp(-distsigma)
                 enddo

              enddo            ! itrans = 1, rsymm%ntrans
           enddo               ! jgrid = 1, parallel%mydim
        enddo                  ! icelly = -nreplicay, nreplicay
       enddo                   ! icellx = -nreplicax, nreplicax

        call psum(f_l,3,parallel%group_size,parallel%group_comm)
        call psum(f_g,3,parallel%group_size,parallel%group_comm)
        call psum(f_c,3,parallel%group_size,parallel%group_comm)
        call psum(f_core,3,parallel%group_size,parallel%group_comm)

        !  Report forces explicitly, if asked (only for debugging).
        if (parallel%iammaster .and. ipr >= 2) then
           write(7,20) ja,f_l*hcub,f_c*hcub,f_core*hcub,f_g*hcub
        endif

        !  Update total force arrays (adding the three components to the
        !  ion-ion forces).
        clust%force(:,ja) = clust%force(:,ja) + (f_l + f_c + f_g)*hcub
        if (jcore > 0) &
           clust%force(:,ja) = clust%force(:,ja) + f_core*hcub
     enddo                  ! jat = 1, clust%natmi(ity)
  enddo                     ! ity = 1, clust%type_num

20 format('at',i4,1x,3(f11.6,1x),1x,3(f11.6,1x),/,7x,3(f11.6,1x),1x,3(f11.6,1x))

  if (parallel%iammaster .and. ipr >= 2) write(7,*)

  deallocate(dvhxc, vxcav)
  deallocate(grad_vh)

end subroutine forcloc_slab

!===============================================================
