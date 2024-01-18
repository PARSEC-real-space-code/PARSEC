!===============================================================
!
!  Copyright (C) 2005 Finite Difference Research Group
!  This file is part of parsec, http://www.ices.utexas.edu/parsec/
!
!  Define a series of quantities over the processors.
!  Coefficients for the numerical derivatives are also defined.
!
!---------------------------------------------------------------
subroutine init_var(clust,band_st,elec_st,mixer,solver,symm,rsymm,parallel, &
     nloc_p_pot,grid,p_pot,pbc,ipr,export_griddata_flag,nscf,outputgw,enable_data_out,ierr)

  use constants
  use cluster_module
  use electronic_struct_module
  use mixer_module
  use eigen_solver_module
  use parallel_data_module
  use non_local_psp_module
  use grid_module
  use pseudo_potential_module
  use pbc_module
  use symmetry_module
  use bandstruc_module
  use nscf_module
#ifdef MPI
  !  include mpi definitions
  use mpi
#endif
  implicit none

  !
  !  Input/Output variables:
  !
  !  the cluster
  type (cluster), intent(inout) :: clust
  !  electronic structure
  type (electronic_struct), intent(inout) :: elec_st
  !  mixer related data
 ! bandstructure structure
  type (bandstruc), intent(inout) :: band_st
  type (mixer_data), intent(inout) :: mixer
  !  eigen solver data
  type (eigen_solver), intent(inout) :: solver
  !  parallel data
  type (parallel_data), intent(in) :: parallel
  !  non local pseudopotential related data
  type (nonloc_pseudo_potential), intent(inout) :: nloc_p_pot 
  !  grid related data
  type (grid_data), intent(inout) :: grid
  !  pseudopotential related data
  type (pseudo_potential), intent(inout) :: p_pot 
  !  periodic boundary conditions data
  type (pbc_data), intent(inout) :: pbc
  !  symmetry operations in the reduced group
  type (symmetry), intent(inout) :: symm,rsymm
  !  printing flag
  integer, intent(inout) :: ipr
  !  print/export flag
  integer, intent(inout) :: export_griddata_flag(MAX_EXPORT_OPTS)
  ! nscf data structure
  type(nscf_data), intent(inout) :: nscf
  !  logical flag for whether gwoutput is written
  logical, intent(inout) :: outputgw
  !  logical flag for whether raw data is written
  logical, intent(inout) :: enable_data_out
  !  error flag
  integer, intent(out) :: ierr
  !
  !  Work variables:
  !
  !  communicator
  integer comm
  !  rank of the master in the above communicator
  integer masterid
  !  exit code for mpi calls
  integer mpinfo
  !  counters
  integer ii, jj, kk, ity, irp, nmax, kpnum, kplp, isp
  real(dp) kpnorm, ktmp(3)

  !---------------------------------------------------------------

  !
  !  Master updates some variables.
  !
  if (parallel%iammaster) then
     !  Change ekbi into h^3/ekbi, to be used from now on.
     do ity  = 1, clust%type_num
        do ii = 1, p_pot%nlocp(ity)
           if (ii /= p_pot%loc(ity)) then
               p_pot%ekbi(ii,ity) =  grid%hcub / p_pot%ekbi(ii,ity)
           endif
        enddo
     enddo
     !  For confined system, define grid%step; grid%step(1) will be
     !  used in place of grid%stepin.
     if (.not. pbc%is_on) then
        grid%step = grid%stepin
        grid%n1 = idint(two*grid%rmax/grid%stepin) + 2
        grid%n2 = grid%n1
        grid%n3 = grid%n1
        pbc%latt_vec = zero
        do ii = 1, 3
           pbc%latt_vec(ii,ii) = one
        enddo
        pbc%avec_norm = pbc%latt_vec
        pbc%box_size = zero
     endif
     !  Transport important information from rsymm to elec_st.
     elec_st%nrep = rsymm%ntrans
     allocate(elec_st%chi(elec_st%nrep,elec_st%nrep))
     elec_st%chi = rsymm%chi
  endif

#ifdef MPI
  masterid = parallel%masterid
  comm = parallel%comm
  !
  !  Broadcast info to all PEs.
  !
#ifdef AJB_DEBUG
     write(9,*) ' beggining init_var BCASTs'
#endif
  call MPI_BCAST(pbc%is_on, 1,MPI_LOGICAL,masterid,comm,mpinfo)
  call MPI_BCAST(pbc%per, 1,MPI_INTEGER,masterid,comm,mpinfo)
  call MPI_BCAST(pbc%box_size, 3,MPI_DOUBLE_PRECISION,masterid,comm,mpinfo)
  call MPI_BCAST(pbc%latt_vec, 9,MPI_DOUBLE_PRECISION,masterid,comm,mpinfo)
  call MPI_BCAST(pbc%bvec, 9,MPI_DOUBLE_PRECISION,masterid,comm,mpinfo)
  call MPI_BCAST(pbc%adot,  9,MPI_DOUBLE_PRECISION,masterid,comm,mpinfo)
  call MPI_BCAST(pbc%bdot,  9,MPI_DOUBLE_PRECISION,masterid,comm,mpinfo)
  call MPI_BCAST(pbc%avec_norm,  9,MPI_DOUBLE_PRECISION,masterid,comm,mpinfo)
  call MPI_BCAST(grid%grad_bvec_norm,  9,MPI_DOUBLE_PRECISION,masterid,comm,mpinfo)
  call MPI_BCAST(pbc%vcell,  1,MPI_DOUBLE_PRECISION,masterid,comm,mpinfo)
  call MPI_BCAST(pbc%a_surface_cell,  1,MPI_DOUBLE_PRECISION,masterid,comm,mpinfo)
  call MPI_BCAST(rsymm%ntrans,  1,MPI_INTEGER,masterid,comm,mpinfo)
  call MPI_BCAST(rsymm%alatt,  9,MPI_DOUBLE_PRECISION,masterid,comm,mpinfo)
  call MPI_BCAST(rsymm%invlat,  9,MPI_DOUBLE_PRECISION,masterid,comm,mpinfo)
  call MPI_BCAST(symm%ntrans,  1,MPI_INTEGER,masterid,comm,mpinfo)
  call MPI_BCAST(symm%alatt,  9,MPI_DOUBLE_PRECISION,masterid,comm,mpinfo)
  call MPI_BCAST(symm%invlat,  9,MPI_DOUBLE_PRECISION,masterid,comm,mpinfo)
  call MPI_BCAST(clust%atom_num,  1,MPI_INTEGER,masterid,comm,mpinfo)
  call MPI_BCAST(clust%has_ptchrg,  1,MPI_LOGICAL,masterid,comm,mpinfo)
  call MPI_BCAST(clust%has_charged_sheet,  1,MPI_LOGICAL,masterid,comm,mpinfo)
  call MPI_BCAST(grid%rmax,  1,MPI_DOUBLE_PRECISION,masterid,comm,mpinfo)
  call MPI_BCAST(grid%norder,  1,MPI_INTEGER,masterid,comm,mpinfo)

  call MPI_BCAST(grid%ndouble,  1,MPI_INTEGER,masterid,comm,mpinfo)
  call MPI_BCAST(grid%domain_shape,  1,MPI_INTEGER,masterid,comm,mpinfo)
  call MPI_BCAST(grid%d_shape_param,  3,MPI_DOUBLE_PRECISION,masterid,comm,mpinfo)
  call MPI_BCAST(grid%i_shape_param,  1,MPI_INTEGER,masterid,comm,mpinfo)
  call MPI_BCAST(grid%hartree_neibs_flag,  1,MPI_LOGICAL,masterid,comm,mpinfo)
!  call MPI_BCAST(p_pot%rws,clust%type_num,MPI_DOUBLE_PRECISION, 
!       parallel%masterid,parallel%comm,mpinfo)
  call MPI_BCAST(elec_st%scf_so,  1,MPI_LOGICAL,masterid,comm,mpinfo)
  call MPI_BCAST(elec_st%cplx,  1,MPI_LOGICAL,masterid,comm,mpinfo)
  call MPI_BCAST(elec_st%nsave,  1,MPI_INTEGER,masterid,comm,mpinfo)
  call MPI_BCAST(elec_st%nstate,  1,MPI_INTEGER,masterid,comm,mpinfo)

  call MPI_BCAST(elec_st%use_plain_sre,  1,MPI_LOGICAL,masterid,comm,mpinfo)
  call MPI_BCAST(elec_st%explicit_hartree_pbc,  1,MPI_LOGICAL,masterid,comm,mpinfo)

  call MPI_BCAST(elec_st%nspin,  1,MPI_INTEGER,masterid,comm,mpinfo)
  call MPI_BCAST(elec_st%ncl,  1,MPI_LOGICAL,masterid,comm,mpinfo)
  call MPI_BCAST(solver%ncl,  1,MPI_LOGICAL,masterid,comm,mpinfo)
  call MPI_BCAST(elec_st%nrep,  1,MPI_INTEGER,masterid,comm,mpinfo)
  call MPI_BCAST(elec_st%dioniz,  1,MPI_DOUBLE_PRECISION,masterid,comm,mpinfo)
  call MPI_BCAST(elec_st%icorr,  2,MPI_CHARACTER,masterid,comm,mpinfo)

  call MPI_BCAST(band_st%bands_on,  1,MPI_LOGICAL,masterid,comm,mpinfo)
if (band_st%bands_on) then 
  call MPI_BCAST(band_st%npoints,  1,MPI_INTEGER,masterid,comm,mpinfo)
  call MPI_BCAST(band_st%nlines,  1,MPI_INTEGER,masterid,comm,mpinfo)
  if (.not. parallel%iammaster)  then
     allocate(band_st%blines(band_st%nlines))
  endif
  do ii = 1, band_st%nlines
      call MPI_BCAST(band_st%blines(ii)%start,  3,MPI_DOUBLE_PRECISION,masterid,comm,mpinfo)
      call MPI_BCAST(band_st%blines(ii)%end,  3,MPI_DOUBLE_PRECISION,masterid,comm,mpinfo)
  enddo
endif

  call MPI_BCAST(nscf%nscf_on,  1,MPI_LOGICAL,masterid,comm,mpinfo)
if (nscf%nscf_on) then 
  call MPI_BCAST(nscf%nkpt,  1,MPI_INTEGER,masterid,comm,mpinfo)
  call MPI_BCAST(nscf%nstate,  1,MPI_INTEGER,masterid,comm,mpinfo)
  call MPI_BCAST(nscf%kgrid, 3,MPI_INTEGER  ,masterid,comm,mpinfo)
  call MPI_BCAST(nscf%kshift, 3,MPI_DOUBLE_PRECISION  ,masterid,comm,mpinfo)
  if (.not. parallel%iammaster) then
     allocate(nscf%kpts(3,nscf%nkpt))
     allocate(nscf%kpwt(nscf%nkpt))
  endif
  call MPI_BCAST(nscf%kpts, 3*nscf%nkpt,  MPI_DOUBLE_PRECISION,masterid,comm,mpinfo)
  call MPI_BCAST(nscf%kpwt, nscf%nkpt,  MPI_DOUBLE_PRECISION,masterid,comm,mpinfo)
endif

  call MPI_BCAST(solver%nadd,  1,MPI_INTEGER,masterid,comm,mpinfo)
  call MPI_BCAST(solver%maxmv,  1,MPI_INTEGER,masterid,comm,mpinfo)
  call MPI_BCAST(solver%lpole,  1,MPI_INTEGER,masterid,comm,mpinfo)
  call MPI_BCAST(solver%full_hartree_flag,  1,MPI_LOGICAL,masterid,comm,mpinfo)
  call MPI_BCAST(mixer%name,  1,MPI_INTEGER,masterid,comm,mpinfo)
  call MPI_BCAST(mixer%memory,  1,MPI_INTEGER,masterid,comm,mpinfo)
  call MPI_BCAST(mixer%restart,  1,MPI_INTEGER,masterid,comm,mpinfo)
  call MPI_BCAST(mixer%block,  1,MPI_INTEGER,masterid,comm,mpinfo)
  call MPI_BCAST(mixer%param,  1,MPI_DOUBLE_PRECISION,masterid,comm,mpinfo)
  call MPI_BCAST(mixer%en_stage,  1,MPI_INTEGER,masterid,comm,mpinfo)
  call MPI_BCAST(mixer%group_size,  1,MPI_INTEGER,masterid,comm,mpinfo)
  call MPI_BCAST(mixer%preferred_type,  1,MPI_INTEGER,masterid,comm,mpinfo)
  call MPI_BCAST(mixer%update_type,  1,MPI_INTEGER,masterid,comm,mpinfo)
  call MPI_BCAST(mixer%restart_factor,  1,MPI_DOUBLE_PRECISION,masterid,comm,mpinfo)
  call MPI_BCAST(mixer%expand_factor,  1,MPI_DOUBLE_PRECISION,masterid,comm,mpinfo)
  call MPI_BCAST(elec_st%nkpt,  1,MPI_INTEGER,masterid,comm,mpinfo)
  call MPI_BCAST(solver%toler,  1,MPI_DOUBLE_PRECISION,masterid,comm,mpinfo)
  call MPI_BCAST(solver%dyntol,  1,MPI_LOGICAL,masterid,comm,mpinfo)
  call MPI_BCAST(solver%polym0,  1,MPI_INTEGER,masterid,comm,mpinfo)
  call MPI_BCAST(solver%ff_maxiter,  1,MPI_INTEGER,masterid,comm,mpinfo)
  call MPI_BCAST(solver%polym,  1,MPI_INTEGER,masterid,comm,mpinfo)
  call MPI_BCAST(solver%polym_t,  1,MPI_INTEGER,masterid,comm,mpinfo)
  call MPI_BCAST(solver%dpm,  1,MPI_INTEGER,masterid,comm,mpinfo)
  call MPI_BCAST(solver%mv_blksize,  1,MPI_INTEGER,masterid,comm,mpinfo)
  call MPI_BCAST(solver%fix_neig,  1,MPI_LOGICAL,masterid,comm,mpinfo)
  call MPI_BCAST(solver%do_subsp,  1,MPI_LOGICAL,masterid,comm,mpinfo)
  call MPI_BCAST(solver%experimental,  1,MPI_LOGICAL,masterid,comm,mpinfo)
  call MPI_BCAST(grid%max_lap_buffers,  1,MPI_LOGICAL,masterid,comm,mpinfo)
  call MPI_BCAST(grid%experimental,  1,MPI_LOGICAL,masterid,comm,mpinfo)
  call MPI_BCAST(grid%hcub,  1,MPI_DOUBLE_PRECISION,masterid,comm,mpinfo)
  call MPI_BCAST(grid%hcub2,  1,MPI_DOUBLE_PRECISION,masterid,comm,mpinfo)
  call MPI_BCAST(grid%h_2,  1,MPI_DOUBLE_PRECISION,masterid,comm,mpinfo)
  call MPI_BCAST(solver%name,6,MPI_CHARACTER,masterid,comm,mpinfo)
  call MPI_BCAST(ipr,1,MPI_INTEGER,masterid,comm,mpinfo)
  call MPI_BCAST(export_griddata_flag,MAX_EXPORT_OPTS,MPI_INTEGER,masterid,comm,mpinfo)

  call MPI_BCAST(outputgw,  1,MPI_LOGICAL,masterid,comm,mpinfo)
  call MPI_BCAST(enable_data_out,  1,MPI_LOGICAL,masterid,comm,mpinfo)

  !  Send the step sizes. NOTE: grid%stepin is not broadcasted.
  call MPI_BCAST(grid%step,  3,MPI_DOUBLE_PRECISION,masterid,comm,mpinfo)
#endif
  !
  !  All PEs (except master) define/allocate variables.
  !
  if (.not. parallel%iammaster) then
     !  Allocate arrays and define structures.
#ifdef AJB_DEBUG
write(9,*) ' Creating elec struct here because I am not master '
#endif
     call create_electronic_struct(elec_st%nspin,elec_st%nstate,  elec_st,elec_st%mxwd)
     if (elec_st%nsave > 0) allocate(elec_st%indxsave(elec_st%nsave))
     allocate(elec_st%chi(elec_st%nrep,elec_st%nrep))
     call create_symmetry (rsymm%ntrans,rsymm)
     call create_symmetry (symm%ntrans,symm)
  endif
#ifdef AJB_DEBUG
write(9,*) ' doing some more BCASTs'
#endif
#ifdef MPI
  call MPI_BCAST(elec_st%do_vdw,1,MPI_LOGICAL,  parallel%masterid,parallel%comm,mpinfo)
  call MPI_BCAST(elec_st%hirsh_3d_cell,3,MPI_INTEGER,  parallel%masterid,parallel%comm,mpinfo)
  call MPI_BCAST(elec_st%ncharge,  1,MPI_DOUBLE_PRECISION,masterid,comm, mpinfo)
  call MPI_BCAST(elec_st%xele,  1,MPI_DOUBLE_PRECISION,masterid,comm, mpinfo)
  call MPI_BCAST(elec_st%tfermi,  1,MPI_DOUBLE_PRECISION,masterid,comm, mpinfo)
  if (elec_st%nsave > 0) then
       call MPI_BCAST(elec_st%indxsave, elec_st%nsave,MPI_INTEGER  ,masterid,comm,mpinfo)
   endif

  call MPI_BCAST(elec_st%chi, elec_st%nrep**2,MPI_INTEGER  ,masterid,comm,mpinfo)
  call MPI_BCAST(rsymm%gmtrx,3*3*rsymm%ntrans  ,MPI_INTEGER,masterid,comm,mpinfo)
  call MPI_BCAST(rsymm%rmtrx,3*3*rsymm%ntrans  ,MPI_INTEGER,masterid,comm,mpinfo)
  call MPI_BCAST(rsymm%trans,3*3*rsymm%ntrans  ,MPI_DOUBLE_PRECISION,masterid,comm,mpinfo)
  call MPI_BCAST(rsymm%tnp,3*rsymm%ntrans  ,MPI_DOUBLE_PRECISION,masterid,comm,mpinfo)
  call MPI_BCAST(rsymm%chi,rsymm%ntrans*rsymm%ntrans  ,MPI_INTEGER,masterid,comm,mpinfo)
  call MPI_BCAST(symm%gmtrx,3*3*symm%ntrans  ,MPI_INTEGER,masterid,comm,mpinfo)
  call MPI_BCAST(symm%rmtrx,3*3*symm%ntrans  ,MPI_INTEGER,masterid,comm,mpinfo)
  call MPI_BCAST(symm%trans,3*3*symm%ntrans  ,MPI_DOUBLE_PRECISION,masterid,comm,mpinfo)
  call MPI_BCAST(symm%tnp,3*symm%ntrans  ,MPI_DOUBLE_PRECISION,masterid,comm,mpinfo)
#endif

  !  If a Monkhorst-Pack grid is input, construct the irreducible
  !  Brillouin zone. Otherwise, distribute k-points.
  if (elec_st%nkpt /= 0) then
#ifdef MPI
     call MPI_BCAST(elec_st%kptmethod, 1,MPI_INTEGER  ,masterid,comm,mpinfo)
#endif
     if (elec_st%kptmethod == MONKHORST_PACK) then
#ifdef MPI
        call MPI_BCAST(elec_st%mpgrid, 3,MPI_INTEGER  ,masterid,comm,mpinfo)
        call MPI_BCAST(elec_st%mpshift, 3,MPI_DOUBLE_PRECISION  ,masterid,comm,mpinfo)
#endif
        call irrbz(elec_st,symm,pbc,parallel%iammaster,ierr)
        if (ierr /= 0) return
     else
        if (.not. parallel%iammaster) then
           allocate(elec_st%kpts(3,elec_st%nkpt))
           allocate(elec_st%kpwt(elec_st%nkpt))
        endif
#ifdef MPI
        call MPI_BCAST(elec_st%kpts, 3*elec_st%nkpt,  MPI_DOUBLE_PRECISION,masterid,comm,mpinfo)
        call MPI_BCAST(elec_st%kpwt, elec_st%nkpt,  MPI_DOUBLE_PRECISION,masterid,comm,mpinfo)
#endif
     endif
  endif

#ifdef AJB_DEBUG
write(9,*) ' not doing BCASTs for a while'
#endif
  if (elec_st%ncl) then
     if(.not. parallel%iammaster) allocate(elec_st%init_mag(3,clust%atom_num))
  endif
  !
  !  Define kpnum = max(1,elec_st%nkpt).
  !
  kpnum = max(elec_st%nkpt,1)

  allocate (elec_st%magmom (2*elec_st%nstate,kpnum))

  kpnorm = zero
  !
  !  Get coefficients for numerical first derivative based on the
  !  order specified by the user. These coefficients are used only
  !  for GGA purposes.
  !
  allocate(grid%coe1(-grid%norder:grid%norder,3))
  call fornberg(1,grid%norder,grid%coe1(-grid%norder,1),ierr)
#ifdef MPI
#ifdef AJB_DEBUG
write(9,*) ' allreducing ! how fun'
#endif
  call MPI_ALLREDUCE(ierr,ii,1,MPI_INTEGER,MPI_MAX,comm,mpinfo)
  ierr = ii
#endif

  if (ierr /= 0) then
     write(9,*) ' ERROR: Wrong parameter in call of Fornberg'
     write(9,*) ' STOP in init_var'
     return
  endif

  !  Renormalize coes with respect to the grid spacing.
  grid%coe1(:,3) = grid%coe1(:,1)/grid%step(3)
  grid%coe1(:,2) = grid%coe1(:,1)/grid%step(2)
  grid%coe1(:,1) = grid%coe1(:,1)/grid%step(1)

  !  Create coefficients for the gradient term of 
  !  kpoint kinetic energy.

  pbc%nkpt = elec_st%nkpt
  nloc_p_pot%nkpt = elec_st%nkpt
  if (elec_st%nkpt /= 0) then

     allocate(pbc%kpts(3,elec_st%nkpt))
     allocate(nloc_p_pot%kpts(3,elec_st%nkpt))
     pbc%kpts(:,:) = elec_st%kpts(:,:)
     nloc_p_pot%kpts(:,:) = elec_st%kpts(:,:)

     allocate(grid%kecoe1(elec_st%nkpt,6,0:grid%norder))
     grid%kecoe1(:,:,:) = zero
     do kplp = 1, elec_st%nkpt
        do jj =0,grid%norder
           call matvec3('T',grid%grad_bvec_norm,elec_st%kpts(1,kplp),ktmp)
           grid%kecoe1(kplp,1:3,jj) = -2*zi*grid%coe1(jj,1:3)*ktmp
        enddo
     enddo

     do kplp = 1, kpnum
        kpnorm = kpnorm + elec_st%kpwt(kplp)
     enddo
     elec_st%kpwt(:)=elec_st%kpwt(:)/kpnorm

  else
     allocate(elec_st%kpwt(1))
     elec_st%kpwt(1) = one
  endif

  !  Calculate coefficients for the laplacian in case of
  !  non-orthogonal grid.

  if(pbc%is_on .and. pbc%per>1) then
     call pbc_grid_coefs(pbc, grid, parallel%iammaster, 1.d-9,0)
  else
     grid%lap_dir_num = 0
     grid%lap_dir = 0
     grid%b_lap(1:3) = one
     grid%b_lap(4:6) = zero
  end if

  write(9,'(a,i2)') ' Laplacian number of directions: ',grid%lap_dir_num
  write(9,*) "Laplacian coefficients data:"
  write(9,'(6(f5.2, 1x))') grid%b_lap
  !
  !  Get coefficients for numerical second derivative based on the 
  !  order specified by the user.
  !
  allocate(grid%coe2(-grid%norder:grid%norder,3+grid%lap_dir_num))
  call fornberg(2,grid%norder,grid%coe2(-grid%norder,1),ierr)
#ifdef MPI
  call MPI_ALLREDUCE(ierr,ii,1,MPI_INTEGER,MPI_MAX,comm,mpinfo)
  ierr = ii
! #ifdef AJB_DEBUG
     write(9,*) ' survived first allreduce'
! #endif
#endif

  if (ierr /= 0) then
     write(9,*) ' ERROR: Wrong parameter in call of Fornberg'
     write(9,*) ' STOP in init_var'
     return
  endif

  !  Renormalize coekes with respect to the grid spacing. First to the
  !  new directions for non-orthogonal grid.
  if(grid%lap_dir_num > 0) then
     write(9,*) ' grid%lap_dir_num > 0'
     do ii = 1,grid%lap_dir_num
        grid%coe2(:,3+ii) = -grid%b_lap(3+ii)*grid%coe2(:,1)/ (grid%lap_dir_step(ii)**2)
     end do
  end if

  grid%coe2(:,3) = -grid%b_lap(3)*grid%coe2(:,1)/(grid%step(3)**2)
  grid%coe2(:,2) = -grid%b_lap(2)*grid%coe2(:,1)/(grid%step(2)**2)
  grid%coe2(:,1) = -grid%b_lap(1)*grid%coe2(:,1)/(grid%step(1)**2)
  !
  !  Get coefficients for double-grid interpolation. If there are
  !  not interpolating point between neighbors in the grid (that is,
  !  grid%ndouble = 1), then the coefficients are 1.
  !
  allocate(grid%fdouble(-grid%ndouble+1:grid%ndouble-1))
  grid%invd = one/real(grid%ndouble,dp)
  do ii = -grid%ndouble + 1, grid%ndouble - 1
     grid%fdouble(ii) = (one - real(abs(ii),dp)*grid%invd)*grid%invd
  enddo
#ifdef MPI
  !
  !  Distribute cluster structure (except clust%name).
  !
#ifdef AJB_DEBUG
write(9,*) ' One last BCAST before create cluster'
#endif
  call MPI_BCAST(clust%type_num, 1,MPI_INTEGER  ,masterid,comm,mpinfo)
  if (.not. parallel%iammaster) then
     call create_cluster(clust%atom_num,clust)
     allocate (clust%natmi (clust%type_num))
#ifdef AJB_DEBUG
     write(9,*) ' created cluster because I am not master'
#endif
  endif
  if(elec_st%ncl)   call MPI_BCAST(elec_st%init_mag,3*clust%atom_num  ,MPI_DOUBLE_PRECISION,masterid,comm,mpinfo)

if (solver%name /= TEST) then
  call MPI_BCAST(clust%natmi, clust%type_num,MPI_INTEGER ,masterid,comm,mpinfo)
  call MPI_BCAST(clust%mvat, clust%atom_num,MPI_INTEGER ,masterid,comm,mpinfo)
  call MPI_BCAST(clust%xatm, clust%atom_num,MPI_DOUBLE_PRECISION ,masterid,comm,mpinfo)
  call MPI_BCAST(clust%yatm, clust%atom_num,MPI_DOUBLE_PRECISION ,masterid,comm,mpinfo)
  call MPI_BCAST(clust%zatm, clust%atom_num,MPI_DOUBLE_PRECISION ,masterid,comm,mpinfo)
  call MPI_BCAST(clust%amass, clust%atom_num,MPI_DOUBLE_PRECISION ,masterid,comm,mpinfo)
  call MPI_BCAST(clust%force, clust%atom_num*3,MPI_DOUBLE_PRECISION ,masterid,comm,mpinfo)
#ifdef AJB_DEBUG
     write(9,*) ' just BCASTed the first batch of clust data'
#endif
endif
  if (clust%has_ptchrg) then
     call MPI_BCAST(clust%npttyp, 1,MPI_INTEGER ,masterid,comm,mpinfo)
     call MPI_BCAST(clust%nptchrg, 1,MPI_INTEGER ,masterid,comm,mpinfo)
     if (.not. parallel%iammaster) then
        allocate (clust%nptt (clust%npttyp))
        allocate (clust%qpt (clust%npttyp))
        allocate (clust%xpt (clust%nptchrg))
        allocate (clust%ypt (clust%nptchrg))
        allocate (clust%zpt (clust%nptchrg))
     endif
     call MPI_BCAST(clust%nptt, clust%npttyp,MPI_INTEGER ,masterid,comm,mpinfo)
     call MPI_BCAST(clust%qpt, clust%npttyp,MPI_DOUBLE_PRECISION ,masterid,comm,mpinfo)
     call MPI_BCAST(clust%xpt, clust%nptchrg,MPI_DOUBLE_PRECISION ,masterid,comm,mpinfo)
     call MPI_BCAST(clust%ypt, clust%nptchrg,MPI_DOUBLE_PRECISION ,masterid,comm,mpinfo)
     call MPI_BCAST(clust%zpt, clust%nptchrg,MPI_DOUBLE_PRECISION ,masterid,comm,mpinfo)
  endif
  if (clust%has_charged_sheet) then
     call MPI_BCAST(clust%number_charged_sheet_type, 1,MPI_INTEGER ,masterid,comm,mpinfo)
     call MPI_BCAST(clust%number_charged_sheets, 1,MPI_INTEGER  ,masterid,comm,mpinfo)
     if (.not. parallel%iammaster) then
        allocate (clust%number_charged_sheets_per_type (clust%number_charged_sheet_type))
        allocate (clust%sheet_charge_of_type (clust%number_charged_sheet_type))
        allocate (clust%z_charged_sheets (clust%number_charged_sheets))
     endif
     call MPI_BCAST(clust%number_charged_sheets_per_type, clust%number_charged_sheet_type,MPI_INTEGER ,masterid,comm,mpinfo)
     call MPI_BCAST(clust%sheet_charge_of_type, clust%number_charged_sheet_type,MPI_DOUBLE_PRECISION ,masterid,comm,mpinfo)
     call MPI_BCAST(clust%z_charged_sheets, clust%number_charged_sheets,MPI_DOUBLE_PRECISION ,masterid,comm,mpinfo)
  endif
#ifdef AJB_DEBUG
write(9,*) ' almost done with init bcasts'
#endif
if (solver%name /= TEST) then
  !
  !  Distribute entire p_pot structure.
  !
  call MPI_BCAST(p_pot%type_num, 1,MPI_INTEGER ,masterid,comm,mpinfo)
  call MPI_BCAST(p_pot%norder, 1,MPI_INTEGER ,masterid,comm,mpinfo)
  call MPI_BCAST(p_pot%is_so, 1,MPI_LOGICAL  ,masterid,comm,mpinfo)
  if (.not. parallel%iammaster) then
       call create_pseudo_potential(p_pot%type_num, p_pot%norder,p_pot)
   endif
  call MPI_BCAST(p_pot%ns, p_pot%type_num,MPI_INTEGER ,masterid,comm,mpinfo)

  if (.not. parallel%iammaster) then
     nmax = maxval(p_pot%ns)
     call pseudo_potential_set_mxpot (nmax,p_pot)
#ifdef AJB_DEBUG
     write(9,*) ' survived set_mxpot'
#endif
  endif
  call MPI_BCAST(p_pot%nlocp, p_pot%type_num,MPI_INTEGER ,masterid,comm,mpinfo)
  call MPI_BCAST(p_pot%loc, p_pot%type_num,MPI_INTEGER ,masterid,comm,mpinfo)
  call MPI_BCAST(p_pot%icore, p_pot%type_num,MPI_INTEGER ,masterid,comm,mpinfo)
  call MPI_BCAST(p_pot%rcore, p_pot%type_num,MPI_DOUBLE_PRECISION ,masterid,comm,mpinfo)
  call MPI_BCAST(p_pot%eleatm, p_pot%type_num*4,MPI_DOUBLE_PRECISION ,masterid,comm,mpinfo)
  call MPI_BCAST(p_pot%spol, p_pot%type_num,MPI_DOUBLE_PRECISION ,masterid,comm,mpinfo)
  call MPI_BCAST(p_pot%zion, p_pot%type_num,MPI_DOUBLE_PRECISION ,masterid,comm,mpinfo)
  call MPI_BCAST(p_pot%ekbi, p_pot%type_num*4,MPI_DOUBLE_PRECISION ,masterid,comm,mpinfo)
  call MPI_BCAST(p_pot%par_a, p_pot%type_num,MPI_DOUBLE_PRECISION ,masterid,comm,mpinfo)
  call MPI_BCAST(p_pot%par_b, p_pot%type_num,MPI_DOUBLE_PRECISION ,masterid,comm,mpinfo)
  call MPI_BCAST(p_pot%par_c, p_pot%type_num,MPI_DOUBLE_PRECISION ,masterid,comm,mpinfo)
  call MPI_BCAST(p_pot%alpha, p_pot%type_num,MPI_DOUBLE_PRECISION ,masterid,comm,mpinfo)
  call MPI_BCAST(p_pot%beta1, p_pot%type_num,MPI_DOUBLE_PRECISION ,masterid,comm,mpinfo)
  call MPI_BCAST(p_pot%acore, p_pot%type_num,MPI_DOUBLE_PRECISION ,masterid,comm,mpinfo)
  call MPI_BCAST(p_pot%uspline, p_pot%type_num,MPI_LOGICAL ,masterid,comm,mpinfo)

#ifdef AJB_DEBUG
     write(9,*) ' midway through p_pot BCASTs'
#endif
  nmax = (p_pot%mxpot + p_pot%norder + 1) * p_pot%type_num
  call MPI_BCAST(p_pot%rs, nmax,MPI_DOUBLE_PRECISION ,masterid,comm,mpinfo)
  call MPI_BCAST(p_pot%denc, nmax,MPI_DOUBLE_PRECISION ,masterid,comm,mpinfo)
  call MPI_BCAST(p_pot%vion, nmax,MPI_DOUBLE_PRECISION ,masterid,comm,mpinfo)
  call MPI_BCAST(p_pot%rho_r, nmax,MPI_DOUBLE_PRECISION ,masterid,comm,mpinfo)
  call MPI_BCAST(p_pot%drhodr, nmax,MPI_DOUBLE_PRECISION ,masterid,comm,mpinfo)
  call MPI_BCAST(p_pot%d2rhodr, nmax,MPI_DOUBLE_PRECISION ,masterid,comm,mpinfo)
  call MPI_BCAST(p_pot%dvion, nmax,MPI_DOUBLE_PRECISION ,masterid,comm,mpinfo)
  call MPI_BCAST(p_pot%d2vion, nmax,MPI_DOUBLE_PRECISION ,masterid,comm,mpinfo)
  call MPI_BCAST(p_pot%ddenc, nmax,MPI_DOUBLE_PRECISION ,masterid,comm,mpinfo)
  call MPI_BCAST(p_pot%d2denc, nmax,MPI_DOUBLE_PRECISION ,masterid,comm,mpinfo)
  call MPI_BCAST(p_pot%vw, nmax*4,MPI_DOUBLE_PRECISION ,masterid,comm,mpinfo)
  call MPI_BCAST(p_pot%dvw, nmax*4,MPI_DOUBLE_PRECISION ,masterid,comm,mpinfo)
  call MPI_BCAST(p_pot%d2vw, nmax*4,MPI_DOUBLE_PRECISION ,masterid,comm,mpinfo)
  call MPI_BCAST(p_pot%wfspd, nmax*4,MPI_DOUBLE_PRECISION ,masterid,comm,mpinfo)
  call MPI_BCAST(p_pot%dvr_ion, 2*nmax,MPI_DOUBLE_PRECISION ,masterid,comm,mpinfo)
  call MPI_BCAST(p_pot%d2vr_ion, 2*nmax,MPI_DOUBLE_PRECISION ,masterid,comm,mpinfo)
  call MPI_BCAST(p_pot%dvr_so, 2*nmax,MPI_DOUBLE_PRECISION ,masterid,comm,mpinfo)
  call MPI_BCAST(p_pot%d2vr_so, 2*nmax,MPI_DOUBLE_PRECISION ,masterid,comm,mpinfo)

  if (p_pot%is_so) then
     call MPI_BCAST(p_pot%so, p_pot%type_num,MPI_LOGICAL ,masterid,comm,mpinfo)
     call MPI_BCAST(p_pot%cc, 2*p_pot%type_num ,MPI_DOUBLE_PRECISION,masterid,comm,mpinfo)

     nmax = (p_pot%mxpot + p_pot%norder + 1) * p_pot%type_num * 2

     call MPI_BCAST(nmax, 1,MPI_INTEGER,masterid,comm,mpinfo)
     call MPI_BCAST(p_pot%vsor, nmax,MPI_DOUBLE_PRECISION ,masterid,comm,mpinfo)
     call MPI_BCAST(p_pot%vionr, nmax,MPI_DOUBLE_PRECISION  ,masterid,comm,mpinfo)
  else
     p_pot%so = .false.
  endif

endif !solver%name /= TEST
#ifdef AJB_DEBUG
     write(9,*) ' survived BCASTS'
#endif
#endif
  !
  !  Create nloc_p_pot structure.
  !
  nmax = 0
  do ity = 1, clust%type_num
     do ii = 1, p_pot%nlocp(ity)
        if (ii /= p_pot%loc(ity)) nmax = nmax + (2*ii - 1)*clust%natmi(ity)
     enddo
  enddo
  call create_nonloc_pseudo_pot (nmax,clust%atom_num,nloc_p_pot)

  allocate(elec_st%eig(elec_st%nrep,kpnum,elec_st%nspin))
  nmax = (solver%winsize + solver%nadd + elec_st%nstate) * elec_st%nrep
  allocate(elec_st%irep(nmax,kpnum,elec_st%nspin))

  elec_st%irep = 1



  !  Create eigenstate structures with empty arrays.
  do isp = 1, elec_st%nspin
     do kplp=1,kpnum
        do irp = 1, elec_st%nrep
           allocate(elec_st%eig(irp,kplp,isp)%en(1))
           allocate(elec_st%eig(irp,kplp,isp)%occ(1))
           elec_st%eig(irp,kplp,isp)%occ = zero
           allocate(elec_st%eig(irp,kplp,isp)%zwf(1,1))
           allocate(elec_st%eig(irp,kplp,isp)%wf(1,1))
           call destroy_eigenstate(elec_st%eig(irp,kplp,isp))
        enddo
     enddo
  enddo
  elec_st%eig(:,:,:)%mm = -1
  elec_st%eig(:,:,:)%nn = -1
  elec_st%eig(:,:,:)%nec = -1

end subroutine init_var
!===============================================================
