!===============================================================
!
! Copyright (C) 2005 Finite Difference Research Group
! This file is part of parsec, http://www.ices.utexas.edu/parsec/
!
! This subroutine sets up the superposition of the LOCAL ionic
! potentials on the grid. This part of the potential changes only
! when the atoms move. It is fixed throughout the iterations to
! self-consistency for a given atomic configuration.
!
!---------------------------------------------------------------
subroutine ionpot(clust,grid,p_pot,parallel,vion_d,oldinpformat)

  use constants
  use cluster_module
  use grid_module
  use pseudo_potential_module
  use parallel_data_module
  implicit none
  !
  ! Input/Output variables:
  !
  ! the cluster
  type (cluster), intent(in) :: clust
  ! grid related data
  type (grid_data), intent(in) :: grid
  ! pseudo_potential related data
  type (pseudo_potential), intent(in) :: p_pot
  ! parallel computation related data
  type (parallel_data), intent(in) :: parallel
  ! local part of ionic potential
  real(dp), intent(out) :: vion_d(parallel%mydim)
  ! type of interpolation flag
  logical, intent(in) :: oldinpformat
  !
  ! Work variables:
  !
  ! temporary storage variables
  real(dp) :: dist, rindex, delr, vtemp, dvtemp, dxx, dyy, dzz
  integer jok
  ! counters
  integer ja, ity, iat, igrid, jgrid, ipt

  !---------------------------------------------------------------

  vion_d(:) = zero

  ! initialize ja - counter of all atoms in the system
  ja = 0
  ! go over all atom types

  do ity = 1, clust%type_num
     ! Go over all atoms, superimpose the contribution to the local
     ! ionic potential from each atom over each point of the grid.
     do iat = 1, clust%natmi(ity)
        ja = ja + 1
        do jgrid = 1, parallel%mydim
           igrid = jgrid + parallel%irows(parallel%group_iam) - 1

           dxx = (grid%shift(1) + grid%kx(igrid))*grid%step(1) &
                - clust%xatm(ja)
           dyy = (grid%shift(2) + grid%ky(igrid))*grid%step(2) &
                - clust%yatm(ja)
           dzz = (grid%shift(3) + grid%kz(igrid))*grid%step(3) &
                - clust%zatm(ja)

           dist = sqrt(dxx*dxx+dyy*dyy+dzz*dzz)
           !
           ! If dist is larger than the largest value of r given in the
           ! pseudopotential file, rs(i,ity), use -zion/r as the ionic
           ! potential, vtemp. The factor of two in zion is conversion from
           ! Hartree to Rydberg.
           if (dist >= p_pot%rs(p_pot%ns(ity)-1,ity)) then
              vtemp = -two*p_pot%zion(ity)/dist
           else
              !
              ! if dist is smaller than the largest value of r, find the index
              ! of dist in the pseudopotential by inverting the logarithmic
              ! pseudopotential grid, then interpolate the LOCAL
              ! pseudopotential to determine the potential and charge 
              ! density
              !
              ! rindex = 1/B(ity)*log((C(ity)+dist)/A(ity))
              ! A = p_pot%par_a
              ! B = p_pot%par_b
              ! C = p_pot%par_c
              !
              ! AJB: however, if par_c = 0 and dist -> 0 you are in danger
              ! of getting a negative or -Inf index.

              rindex = one/p_pot%par_b(ity)* &
                   log((p_pot%par_c(ity)+dist)/p_pot%par_a(ity))
              jok  = idint(rindex) + 1
              ! AJB: Hack. That is not what we would like eventually
              if (jok <1 ) then
               write(9,*) "ionpbc, warning: tough to interpolate this pseudo"
               write(9,*) "dist, rindex,jok:", dist, rindex,jok
               write(9,*) "setting jok to be 1"
               jok = 1
              !write(9,*) jok,p_pot%rs(jok,ity)
              endif
              delr = (dist-p_pot%rs(jok,ity))/ &
                   (p_pot%rs(jok+1,ity)-p_pot%rs(jok,ity))
              ! place pseudopotential on the 3d uniform grid by interpolating
              ! the 1d pseudopotential given on the logarithmic grid
              ! If "old style", interpolate the pseudopotential directly
              if (oldinpformat) then
                 vtemp = p_pot%vion(jok,ity) + delr* &
                      (p_pot%vion(jok+1,ity) - &
                      p_pot%vion(jok,ity))
                 ! If "new style", do the interpolation on r*pseudopotential -
                 ! improves the interpolation accuracy because it is exact
                 ! beyond pseudopotential cutoff radius.
              else
                 vtemp = p_pot%vion(jok,ity)*p_pot%rs(jok,ity)/dist + &
                      delr*(p_pot%vion(jok+1,ity)*p_pot%rs(jok+1,ity) &
                      -p_pot%vion(jok,ity)*p_pot%rs(jok,ity))/dist
                 if (jok == 1) vtemp = p_pot%vion(2,ity)
              endif
              if (p_pot%uspline(ity)) then
                 call splint(p_pot%rs(-p_pot%norder:p_pot%ns(ity),ity), &
                             p_pot%vion(-p_pot%norder:p_pot%ns(ity),ity), &
                             p_pot%d2vion(-p_pot%norder:p_pot%ns(ity),ity), &
                             p_pot%ns(ity)+1+p_pot%norder,dist,vtemp,dvtemp)
              endif
           endif

           vion_d(jgrid) = vion_d(jgrid) + vtemp  

        enddo               ! jgrid = 1, parallel%mydim
     enddo                  ! iat = 1, clust%natmi(ity)
  enddo                     ! ity = 1, clust%type_num

  if (clust%has_ptchrg) then
     ! Initialize ipt - counter of all point charges in the system.
     ipt = 0

     ! Go over all point charge types.

     do ity = 1, clust%npttyp
        ! Go over all point charges within each type, superimpose the
        ! contribution to the local ionic potential from each point
        ! charge over each point of the grid
        do iat = 1, clust%nptt(ity)
           ipt = ipt + 1
           do jgrid = 1,parallel%mydim
              igrid = jgrid + parallel%irows(parallel%group_iam) - 1

              dxx = (grid%shift(1) + grid%kx(igrid))*grid%step(1) &
                   - clust%xpt(ipt)
              dyy = (grid%shift(2) + grid%ky(igrid))*grid%step(2) &
                   - clust%ypt(ipt)
              dzz = (grid%shift(3) + grid%kz(igrid))*grid%step(3) &
                   - clust%zpt(ipt)

              dist = sqrt(dxx*dxx+dyy*dyy+dzz*dzz)

              ! The potential is simply -qpt/r, where qpt is the point
              ! charge. The factor of two is conversion from hartrees
              ! to rydbergs.
              ! If PC is too close to grid point -> small dist -> numerical
              ! instability. Threshold set to one -> further testing needed

              if(dist > one) then
                 vtemp = -two*clust%qpt(ity)/dist

                 vion_d(jgrid) = vion_d(jgrid) + vtemp

              endif         ! if(dist > 1) then

           enddo            ! jgrid = 1,parallel%mydim
        enddo               ! iat = 1, clust%nptt(ity)
     enddo                  ! ity = 1, clust%npttyp

  endif                  ! ptflag == .true.

end subroutine ionpot
!===============================================================
