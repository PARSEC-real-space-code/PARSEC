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
!-----------------------------------------------------------
! This subroutine was changed by Ido Azuri, applying it
! for the calculation of the Hirshfeld function for the
! Tkatchenko-Scheffler vdW correction,
! at the Weizmann Institue of Science, at Leeor Kronik group
!-----------------------------------------------------------
subroutine initchrg_f(clust,elec_st,grid,p_pot,pbc,parallel,rhof,dist_atm)

  use constants
  use cluster_module
  use electronic_struct_module
  use grid_module
  use pseudo_potential_module
  use pbc_module
  use parallel_data_module
  implicit none 
  !
  ! Input output variables:
  !
  ! the cluster
  type(cluster), intent(in) :: clust
  ! electronic structure
  type (electronic_struct), intent(in) :: elec_st
  ! grid related data
  type (grid_data), intent(in) :: grid
  ! pseudo_potential related data
  type (pseudo_potential), intent(in) :: p_pot
  ! periodic structure related data
  type (pbc_data), intent(in) :: pbc
  ! parallel computation related data
  type (parallel_data), intent(in) :: parallel
  !  atomic charge for each atom, given on the 3-d grid points
  real(dp), intent(out) :: rhof(parallel%mydim,2*elec_st%nspin-1,clust%atom_num) 
  ! distance of atoms from grid points
  real(dp), intent(out) :: dist_atm(parallel%mydim,clust%atom_num)
 
  ! Work variables:
  !
  integer :: num_cells_super 
  ! total number of valence electrons for each atom calculated from integrating 
  ! over the charge density of each atom.
  real(dp) ::  tele(clust%atom_num),tele_c(clust%atom_num)
  real(dp) ::  tele_s( &
          (elec_st%hirsh_3d_cell(1,1)*2+1)* &
          (elec_st%hirsh_3d_cell(1,2)*2+1)* &
          (elec_st%hirsh_3d_cell(1,3)*2+1),clust%atom_num)
  ! total number of valence electrons for each atom from p_pot%zion
  real(dp) :: xele_atm(clust%atom_num)
  ! temporary storage variables
  real(dp) :: rindex,delr,acha,dxx,dyy,dzz,xa,ya,za,xx,yy,zz
  real(dp) :: rhof_cell(parallel%mydim, &
          (elec_st%hirsh_3d_cell(1,1)*2+1)* &
          (elec_st%hirsh_3d_cell(1,2)*2+1)* &
          (elec_st%hirsh_3d_cell(1,3)*2+1),clust%atom_num)
  real(dp) :: dist_atm_cell(parallel%mydim, &
          (elec_st%hirsh_3d_cell(1,1)*2+1)* &
          (elec_st%hirsh_3d_cell(1,2)*2+1)* &
          (elec_st%hirsh_3d_cell(1,3)*2+1),clust%atom_num)
  real(dp) :: hirsh_denum(parallel%mydim,1)
  real(dp) :: hirsh_func(parallel%mydim, &
          (elec_st%hirsh_3d_cell(1,1)*2+1)* &
          (elec_st%hirsh_3d_cell(1,2)*2+1)* &
          (elec_st%hirsh_3d_cell(1,3)*2+1),clust%atom_num)
  real(dp) :: hirsh_denum_per(parallel%mydim,1)
  real(dp) :: v_free_per(1,clust%atom_num)
  ! polarization densities: p_fac(1) = p_fac(2) + p_fac(3)
  real(dp) :: p_fac(3)
  integer jok
  ! counters
  integer ja,itype,isp,jj,iat,icellx,icelly,icellz,iw,jw,i,j,k,i_per,i_per1,i_per2,i_per3,ja_per
  ! number of replicas of the periodic cell to be used in the
  ! construction of non-local spheres around each atom
  ! nreplica = 0 if .not. pbc%is_on
  integer nreplica(3)
  ! distance from a grid point to atoms
  real(dp) :: dist
  real(dp) :: uvec(3), xvec(3)
  real(dp) :: tmpmat(3,3), tmpnorm
  !---------------------------------------------------------------
  !
  ! Define parameters for periodic boundary conditions.
  num_cells_super = (elec_st%hirsh_3d_cell(1,1)*2+1)*(elec_st%hirsh_3d_cell(1,2)*2+1)*(elec_st%hirsh_3d_cell(1,3)*2+1)
  nreplica = 0
  if (pbc%per > 0) then
        nreplica(1) = elec_st%hirsh_3d_cell(1,1)
        nreplica(2) = elec_st%hirsh_3d_cell(1,2)
        nreplica(3) = elec_st%hirsh_3d_cell(1,3)
  end if
  tmpmat=pbc%latt_vec

  do jj = 1, 3
     tmpnorm = sqrt(sum(tmpmat(:,jj)**2))
     tmpmat(:,jj) = tmpmat(:,jj)/tmpnorm
  end do
  

  ! initialize charge density
  rhof_cell(:,:,:)=zero
  dist_atm_cell(:,:,:)=zero
  hirsh_denum(:,1)=zero
  hirsh_denum_per(:,1)=zero
  hirsh_func(:,:,:)=zero
  rhof(:,:,:) = zero
  dist_atm(:,:)=zero
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
        do jw = 1, parallel%mydim
           iw = jw + parallel%irows(parallel%group_iam) - 1
                
           ! For an improved guess for the initial charge density, add in
           ! contributions from the nearest REPLICAS of the atoms in the
           ! original supercell, if periodic boundary conditions are used.
           i_per=0
           do icellx = -nreplica(1),nreplica(1)
              do icelly = -nreplica(2),nreplica(2)
                 do icellz = -nreplica(3),nreplica(3)
                       i_per=i_per+1
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
                    
                    


                    rindex=one/p_pot%par_b(itype)* &
                         log((p_pot%par_c(itype)+dist)/p_pot%par_a(itype))
                    jok  = int(rindex) + 1
                    delr=(dist-p_pot%rs(jok,itype))/ &
                         (p_pot%rs(jok+1,itype)-p_pot%rs(jok,itype))
                    acha = p_pot%rho_r(jok,itype) + delr* &
                         (p_pot%rho_r(jok+1,itype) - p_pot%rho_r(jok,itype))

   
                               dist_atm_cell(jw,i_per,ja)=dist
                    do isp = 1, elec_st%nspin*2 - 1                    
                       rhof_cell(jw,i_per,ja) = rhof_cell(jw,i_per,ja) + acha*p_fac(isp)
                    enddo
                         if ((icellx==0 .and. icelly==0) .and. icellz==0) then
                               dist_atm(jw,ja)=dist
                         

                    do isp = 1, elec_st%nspin*2 - 1                    
                       rhof(jw,isp,ja) = rhof(jw,isp,ja) + acha*p_fac(isp)
                    enddo
                         end if
                 enddo
              enddo
               
           enddo

        enddo               ! jw = 1, parallel%mydim       
     enddo                  ! iat = 1, clust%natmi(itype)
  enddo                     ! itype = 1, clust%type_num
  elec_st%dist_periodic(:,:,:)=dist_atm_cell
  !
  ! Renormalize the electron count to the correct number:
  ! tele = valence electrons for each atom = integral of rhof over all space
  !      = h**3 * sum(rho)
  !

!  open(2,file='init_f.dat')

  do ja=1,clust%atom_num
  tele(ja) = sum(rhof(:,1,ja))*real(elec_st%nrep,dp)
  end do
  tele = tele*grid%hcub
 
  call psum(tele,clust%atom_num,parallel%group_size,parallel%group_comm)

  do ja=1,clust%atom_num
  do i_per=1,num_cells_super
  tele_s(i_per,ja)=sum(rhof_cell(:,i_per,ja))
  end do
  end do

  do ja=1,clust%atom_num
  tele_c(ja) = sum(tele_s(:,ja))*real(elec_st%nrep,dp)
  end do
  tele_c = tele_c*grid%hcub
 
  call psum(tele_c,clust%atom_num,parallel%group_size,parallel%group_comm)


  ja=0
  do itype=1,clust%type_num
     do iat =1,clust%natmi(itype)
        ja=ja+1
        xele_atm(ja)=p_pot%zion(itype)
     end do
  end do

  do ja=1,clust%atom_num
    rhof(:,1,ja) = xele_atm(ja)/tele(ja)*rhof(:,1,ja)
  end do

  
  do ja=1,clust%atom_num
    rhof_cell(:,:,ja) = xele_atm(ja)/tele_c(ja)*rhof_cell(:,:,ja)
  end do
 

  v_free_per(1,:)=zero

  do ja=1,clust%atom_num
        do i_per=1,num_cells_super
              v_free_per(1,ja)=v_free_per(1,ja)+sum(dist_atm_cell(:,i_per,ja)**3D0*rhof_cell(:,i_per,ja)*grid%hcub)
        end do
  end do


  call psum(v_free_per,clust%atom_num,parallel%group_size,parallel%group_comm)


  elec_st%v_free_periodic(1,:)=v_free_per(1,:)
   

  do i_per=1,num_cells_super
        do ja=1,clust%atom_num
              hirsh_denum_per(:,1)=hirsh_denum_per(:,1)+rhof_cell(:,i_per,ja)
        end do
  end do


  do ja=1,clust%atom_num
        do i_per=1,num_cells_super
              hirsh_func(:,i_per,ja)=rhof_cell(:,i_per,ja)/hirsh_denum_per(:,1)
        end do
  end do


  call psum(hirsh_func,clust%atom_num,parallel%group_size,parallel%group_comm)


  elec_st%hirshfeld_periodic(:,:,:)=hirsh_func(:,:,:)



  if (elec_st%nspin == 2) then
     do ja=1,clust%atom_num
       rhof(:,2,ja) = xele_atm(ja)/tele(ja)*rhof(:,2,ja)
       rhof(:,3,ja) = xele_atm(ja)/tele(ja)*rhof(:,3,ja)
     end do
  endif
 
      
end subroutine initchrg_f
!===============================================================
