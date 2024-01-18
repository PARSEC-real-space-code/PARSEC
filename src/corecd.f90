!===============================================================
!
! Copyright (C) 2005 Finite Difference Research Group
! This file is part of parsec, http://www.ices.utexas.edu/parsec/
!
! This subroutine sets up the superposition of core-correction
! charge densities, if any, associated with all atoms, on the
! grid. This subroutine handles both the case of a periodic
! system and a confined system.
!
! Core correction is described in: S. G. Louie, S. Froyen, and M.
! L. Cohen, Phys. Rev. B26, 1738 (1982).
!
! NOTE: The core charge density is NOT divided by 4*pi*r^2 to
! convert to a volume density. This means that we expect the
! input from the pseudopotential code to be given as a VOLUME
! density even when still given on the radial pseudopotential
! grid (compare with the computation of rho_r in initchrg.f)
!
!---------------------------------------------------------------
subroutine corecd(clust,grid,p_pot,pbc,parallel,rhoc,ierr)

  use constants
  use cluster_module
  use grid_module
  use pseudo_potential_module
  use pbc_module
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
  ! periodic structure related data
  type (pbc_data), intent(in) :: pbc
  ! parallel computation related data
  type (parallel_data), intent(in) :: parallel
  ! superposition of the core-correction charge density,
  ! given on the 3-d grid points
  real(dp), intent(out) :: rhoc(parallel%mydim)
  ! error flag, 380 < ierr < 391
  integer, intent(out) :: ierr
  !
  ! Work variables:
  !
  ! temporary storage variables
  real(dp) :: dist,rindex,delr,corch,dcorch,dxx,dyy,dzz,xa,ya,za
  ! coordinates of grid point with lowest charge density
  real(dp) :: rmin(parallel%group_size,3)
  real(dp), dimension(parallel%group_size) :: dmax,dmin
  real(dp) :: uvec(3), xvec(3)
  integer jmin
  ! number of replicas of the periodic cell to be used in the
  ! construction of non-local spheres around each atom
  ! nreplica = 0 if .not. is_pbc
  integer nreplica(3)
  ! counters
  integer ja,itype,i,iat,iw,jw,icellx,icelly,icellz,jok
  ! minimal numerical tolerance for zero in this subroutine
  real(dp), parameter :: eps = 1.d-1

  !---------------------------------------------------------------
  !
  ! If there is not core-correction charge, skip all this.
  !
  rhoc(:) = zero
  if (maxval(p_pot%icore) == 0) return

  ! Define parameters for periodic boundary conditions.
  nreplica = 0
  if (pbc%per > 0) nreplica(1:pbc%per) = 1
  ! initialize ja - counter of all atoms in the system
  ja = 0
  ! go over all atom types
  do itype = 1, clust%type_num
     ! If core-correction flag is on:
     if (p_pot%icore(itype) == 1) then
        ! go over all atoms and superimpose the contribution 
        ! from each atom over each point of the grid
        do iat = 1, clust%natmi(itype)
           ja = ja + 1
           xa = clust%xatm(ja)
           ya = clust%yatm(ja)
           za = clust%zatm(ja)
           do jw = 1, parallel%mydim
              iw = jw + parallel%irows(parallel%group_iam) - 1

              do icellx = -nreplica(1),nreplica(1)
                 do icelly = -nreplica(2),nreplica(2)
                    do icellz = -nreplica(3),nreplica(3)
                       uvec(1) = (grid%shift(1) + grid%kx(iw)) &
                            *grid%step(1) + real(icellx,dp)*pbc%box_size(1)
                       uvec(2) = (grid%shift(2) + grid%ky(iw)) &
                            *grid%step(2) + real(icelly,dp)*pbc%box_size(2)
                       uvec(3) = (grid%shift(3) + grid%kz(iw)) &
                            *grid%step(3) + real(icellz,dp)*pbc%box_size(3)

                       if(pbc%is_on) then
                          call matvec3('N',pbc%avec_norm,uvec,xvec)
                       else
                          xvec=uvec
                       endif

                       dxx=xvec(1)-xa
                       dyy=xvec(2)-ya
                       dzz=xvec(3)-za

                       dist = sqrt(dxx*dxx+dyy*dyy+dzz*dzz)
                       ! If dist is larger than the largest value of r given
                       ! in the pseudopotential file, rs(i,itype), there is
                       ! no core charge correction.

                       if (dist >= p_pot%rs(p_pot%ns(itype)-1,itype)) cycle
                       !
                       ! if dist is smaller than the largest value of r, find
                       ! the index of dist in the pseudopotential by inverting
                       ! the logarithmic pseudopotential grid, then
                       ! interpolate the radial core-correction charge, denc,
                       ! to determine the charge contribution.

                       ! rindex = 1/B(itype)*log((C(itype)+dist)/A(itype))
                       ! A = p_pot%par_a
                       ! B = p_pot%par_b
                       ! C = p_pot%par_c
                       rindex = one/p_pot%par_b(itype)* &
                            log((p_pot%par_c(itype)+dist)/p_pot%par_a(itype))
                       jok  = idint(rindex) + 1
                       delr=(dist-p_pot%rs(jok,itype))/ &
                            (p_pot%rs(jok+1,itype)-p_pot%rs(jok,itype))
                       corch = p_pot%denc(jok,itype) + delr &
                            *(p_pot%denc(jok+1,itype)-p_pot%denc(jok,itype))
                       if (p_pot%uspline(itype)) then
                          call splint(p_pot%rs(-p_pot%norder:p_pot%ns(itype),itype), &
                                      p_pot%denc(-p_pot%norder:p_pot%ns(itype),itype), &
                                      p_pot%d2denc(-p_pot%norder:p_pot%ns(itype),itype), &
                                      p_pot%ns(itype)+1+p_pot%norder,dist,corch,dcorch)
                       endif
                       rhoc(jw) = rhoc(jw) + corch
                    enddo
                 enddo
              enddo

           enddo            ! jw = 1, parallel%mydim
        enddo               ! iat = 1, clust%natmi(itype)
        ! If no core correction for this type just update the atom counter
     else
        ja = ja + clust%natmi(itype)
     endif
  enddo                     ! itype = 1, clust%type_num
  !
  ! Find minimum and maximum values of core-correction charge
  ! density, and the grid index where these values occur.
  ! Must initialize all arrays to zero since we are doing global
  ! reduction.
  !
  rmin(:,:) = zero
  dmax(:) = zero
  dmin(:) = zero
  dmax(parallel%group_iam+1) = maxval (rhoc)
  dmin(parallel%group_iam+1) = minval (rhoc)
  jmin = maxval (minloc (rhoc))
  rmin(parallel%group_iam+1,1) = &
       (grid%shift(1) + grid%kx(jmin))*grid%step(1)
  rmin(parallel%group_iam+1,2) = &
       (grid%shift(2) + grid%ky(jmin))*grid%step(2)
  rmin(parallel%group_iam+1,3) = &
       (grid%shift(3) + grid%kz(jmin))*grid%step(3)

  call psum(rmin,parallel%group_size*3, &
       parallel%group_size,parallel%group_comm)
  call psum(dmax,parallel%group_size, &
       parallel%group_size,parallel%group_comm)
  call psum(dmin,parallel%group_size, &
       parallel%group_size,parallel%group_comm)

  if (parallel%iammaster) then
     write(7,*) 'Setup of inital charge and Hartree potential:'
     write(7,*) '---------------------------------------------'
     !
     ! Report smallest and largest values of core-correction charge
     ! density.
     !
     write(7,*) 'extremal values of core charge density [e/bohr^3] :'
     write(7,'(g12.5,3x,g12.5)') maxval(dmax),minval(dmin)
11   format(1x,'x,y,z =',1x,f6.3,1x,f6.3,1x,f6.3)
  endif
  !
  ! Look for negative densities. If the negative density is within
  ! numerical accuracy of zero, set it to zero. If not, declare an
  ! error and quit.
  !
  if (minval(dmin) < zero) then
     if (abs(minval(dmin)) < eps) then
        do i = 1, parallel%mydim
           if (rhoc(i) < zero) rhoc(i) = zero
        enddo
     else
        write(9,*)
        write(9,*) 'ERROR: negative core-correction charge density'
        write(9,*) 'Charge is most negative at:'
        write(9,11) rmin(parallel%group_iam+1,:)
        write(9,*) 'STOP in corecd.'
        call myflush(9)
        ierr = 381
        return
     endif
  endif

end subroutine corecd
!===============================================================
