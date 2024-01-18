!===============================================================
!
! Copyright (C) 2005 Finite Difference Research Group
! This file is part of parsec, http://www.ices.utexas.edu/parsec/
!
! This subroutine sets up the superposition of the atomic charge
! densities on the grid, to be used as an initial guess for the
! self -consistent loop, in case of a fresh calculation.
!
! This subroutine handles both boundary conditions: periodic and
! confined.
!
!---------------------------------------------------------------
subroutine initchrg(clust,elec_st,grid,p_pot,pbc,parallel,rho)

  use constants
  use cluster_module
  use electronic_struct_module
  use grid_module
  use pseudo_potential_module
  use pbc_module
  use parallel_data_module
  implicit none 
  !
  ! Input variables:
  !
  ! the cluster
  type(cluster), intent(in) :: clust
  ! electronic structure
  type (electronic_struct), intent(inout) :: elec_st
  ! grid related data
  type (grid_data), intent(in) :: grid
  ! pseudo_potential related data
  type (pseudo_potential), intent(in) :: p_pot
  ! periodic structure related data
  type (pbc_data), intent(in) :: pbc
  ! parallel computation related data
  type (parallel_data), intent(in) :: parallel
  ! superposition of atomic charges, given on the 3-d grid points
  real(dp), intent(out) :: rho(parallel%mydim,2*elec_st%nspin-1) 
  !
  ! Work variables:
  !
  ! total number of electrons calculated from integrating 
  ! over the charge density
  real(dp) :: tele, teleinv
  real(dp), dimension(1) :: televec
  ! temporary storage variables
  real(dp) :: rindex,delr,acha,dxx,dyy,dzz,xa,ya,za
  ! polarization densities: p_fac(1) = p_fac(2) + p_fac(3)
  real(dp) :: p_fac(3)
  integer jok
  ! counters
  integer ja,itype,isp,jj,iat,icellx,icelly,icellz,iw,jw
  ! number of replicas of the periodic cell to be used in the
  ! construction of non-local spheres around each atom
  ! nreplica = 0 if .not. pbc%is_on
  integer nreplica(3)
  ! distance from a grid point to atoms
  real(dp) :: dist, nclm
  real(dp) :: uvec(3), xvec(3)
  real(dp) :: tmpmat(3,3), tmpnorm

  !---------------------------------------------------------------
  !
  ! Define parameters for periodic boundary conditions.
  nreplica = 0
  if (pbc%per > 0) nreplica(1:pbc%per) = 1
  tmpmat=pbc%latt_vec

  do jj = 1, 3
     tmpnorm = sqrt(sum(tmpmat(:,jj)**2))
     tmpmat(:,jj) = tmpmat(:,jj)/tmpnorm
  end do
  if (elec_st%ncl) elec_st%spin3d = zero
  ! initialize charge density
  rho(:,:) = zero

  ! initialize ja - counter of all atoms in the system
  ja = 0
  ! go over all atom types
  do itype = 1, clust%type_num
     ! If spin-polarized, initialize the spin-up and spin-down polarization
     ! densities according to input spin polarization.
     p_fac(1) = one
     p_fac(2) = one + p_pot%spol(itype)
     p_fac(3) = one - p_pot%spol(itype)
     ! go over all atoms superimpose the contribution of atomic charge
     ! and ionic potential from each atom over each point of the grid
     do iat = 1, clust%natmi(itype)
        ja = ja + 1
        xa = clust%xatm(ja)
        ya = clust%yatm(ja)
        za = clust%zatm(ja)

        if(elec_st%ncl) then
           nclm = sqrt(dot_product(elec_st%init_mag(:,ja), &
             elec_st%init_mag(:,ja)))
           nclm = nclm / grid%hcub
        endif

!$OMP PARALLEL DO &
!$OMP& SCHEDULE(RUNTIME)  &
!$OMP& DEFAULT(SHARED) &
!$OMP& PRIVATE(jw,iw,icellx,icelly,icellz,uvec,dxx,dyy,dzz,xvec,dist,jok,rindex,delr,acha,isp)
        do jw = 1, parallel%mydim
           iw = jw + parallel%irows(parallel%group_iam) - 1

           ! For an improved guess for the initial charge density, add in
           ! contributions from the nearest REPLICAS of the atoms in the
           ! original supercell, if periodic boundary conditions are used.
           do icellx = -nreplica(1),nreplica(1)
              do icelly = -nreplica(2),nreplica(2)
                 do icellz = -nreplica(3),nreplica(3)
                    uvec(1) = (grid%shift(1) + grid%kx(iw)) &
                         *grid%step(1) + real(icellx,dp)*pbc%box_size(1)
                    uvec(2) = (grid%shift(2) + grid%ky(iw)) &
                         *grid%step(2) + real(icelly,dp)*pbc%box_size(2)
                    uvec(3) = (grid%shift(3) + grid%kz(iw)) &
                         *grid%step(3) + real(icellz,dp)*pbc%box_size(3)
                    call matvec3('N',pbc%avec_norm,uvec,xvec)

                    dxx = xvec(1) - xa
                    dyy = xvec(2) - ya
                    dzz = xvec(3) - za

                    dist = sqrt(dxx*dxx+dyy*dyy+dzz*dzz)

                    ! If dist is larger than the largest value of r given in
                    ! the pseudopotential file, rs(i,itype), drop this point

                    if (dist >= p_pot%rs(p_pot%ns(itype)-1 &
                         ,itype)) cycle
                    ! If dist is smaller than the largest value of r, find
                    ! the index of dist in the pseudopotential by inverting
                    ! the logarithmic pseudopotential grid, then interpolate
                    ! the atomic charge density to determine the charge
                    ! density on the real-space grid.
                    !
                    ! rindex = 1/B(itype)*log((C(itype)+dist)/A(itype))
                    ! A = p_pot%par_a
                    ! B = p_pot%par_b
                    ! C = p_pot%par_c
                    if (dist < 1D-6) then
                       jok = 1
                    else
                       rindex=one/p_pot%par_b(itype)* &
                         log((p_pot%par_c(itype)+dist)/p_pot%par_a(itype))
                       jok  = int(rindex) + 1
                    endif
                    delr=(dist-p_pot%rs(jok,itype))/ &
                         (p_pot%rs(jok+1,itype)-p_pot%rs(jok,itype))
                    acha = p_pot%rho_r(jok,itype) + delr* &
                         (p_pot%rho_r(jok+1,itype) - p_pot%rho_r(jok,itype))

                    do isp = 1, elec_st%nspin*2 - 1
                       rho(jw,isp) = rho(jw,isp) + acha*p_fac(isp)
                    enddo

                    if(elec_st%ncl) elec_st%spin3d(jw,:) = elec_st%spin3d(jw,:) + &
                         (p_fac(2)-p_fac(3))*acha*elec_st%init_mag(:,ja)/nclm
                 enddo
              enddo
           enddo
        enddo               ! jw = 1, parallel%mydim
!$OMP END PARALLEL DO

     enddo                  ! iat = 1, clust%natmi(itype)
  enddo                     ! itype = 1, clust%type_num

  tele = zero
  teleinv = zero
!$OMP PARALLEL &
!$OMP& DEFAULT(SHARED) 

  do isp = 1, 2*elec_st%nspin - 1
!$OMP DO &
!$OMP& SCHEDULE(STATIC)  &
!$OMP& PRIVATE(jw,iw) 
     do jw =  1, parallel%mydim
        if (rho(jw,isp) < zero) then
           write(9,*) ' WARNING in initchrg : negative ' &
                ,'electron density ',rho(jw,isp),' at '
           write(9,*) grid%kx(iw),grid%ky(iw) ,grid%kz(iw) &
                ,' ispin = ',isp-1
           call myflush(9)
        endif
     enddo
!$OMP END DO 
  enddo
  !
  ! Renormalize the total electron count to the correct number:
  ! tele = total electrons = integral of rho over all space
  !      = h**3 * sum(rho)
  !
!$OMP DO &
!$OMP& SCHEDULE(STATIC)  &
!$OMP& PRIVATE(jw) &
!$OMP& REDUCTION(+:tele) 
     do jw =  1, parallel%mydim
         tele = tele + rho(jw,1)
     end do
!$OMP END DO

!$OMP MASTER
  tele=  tele * real(elec_st%nrep,dp)
  tele = tele * grid%hcub
  televec = tele

  call psum(televec,1,parallel%group_size,parallel%group_comm)
  tele = televec(1)

  !rho(:,1) = elec_st%xele/tele*rho(:,1)
  if (tele>zero) then
      teleinv = 1/tele * elec_st%xele
  else
      write(9,*) "WARNING, no charge?"
      teleinv = zero
  endif
!$OMP END MASTER
!$OMP BARRIER

!$OMP DO &
!$OMP& SCHEDULE(STATIC)  &
!$OMP& PRIVATE(jw) 
     do jw =  1, parallel%mydim
        rho(jw,1)=rho(jw,1)*teleinv
    end do
!$OMP END DO

!$OMP MASTER
  if (parallel%iammaster) then
     write(7,10) tele
10   format(' Tot. electron charge from atomic density [e] = ',g16.5)
     write(7,*) 'Renormalizing charge density'
  endif
!$OMP END MASTER

  if (elec_st%nspin == 2) then
!$OMP DO &
!$OMP& SCHEDULE(STATIC)  &
!$OMP& PRIVATE(jw) 
     do jw =  1, parallel%mydim
     rho(jw,2) = teleinv*rho(jw,2)
     rho(jw,3) = teleinv*rho(jw,3)
     end do
!$OMP END DO
!$OMP MASTER
! also parallelize this one day
     p_fac = sum(rho,dim=1)

     p_fac = p_fac*real(elec_st%nrep,dp)*grid%hcub

     call psum(p_fac,3,parallel%group_size,parallel%group_comm)
     if (parallel%iammaster .and. (.not. elec_st%ncl)) then
        write(7,20) p_fac(2) - p_fac(3)
20      format(' Initial magnetic moment [Bohr magnetons] = ',f12.5)
     endif
!$OMP END MASTER
  endif

!$OMP END PARALLEL
end subroutine initchrg
!===============================================================
