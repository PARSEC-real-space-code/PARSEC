!===============================================================
!
! Copyright (C) 2005 Finite Difference Research Group
! This file is part of parsec, http://www.ices.utexas.edu/parsec/
!
! This subroutine calculates various quantities associated with
! the non-local components of the pseudopotential. They are used
! later for determining this non-local contribution to the
! Hamiltonian matrix and to the calculation of non-local forces.
!
! This subroutine handles both the case of a periodic system
! and a confined system.
!
!---------------------------------------------------------------
subroutine nonloc(clust,grid,p_pot,nloc_p_pot,pbc,parallel,ierr)

  use constants
  use cluster_module
  use grid_module
  use pseudo_potential_module
  use non_local_psp_module
  use pbc_module
  use parallel_data_module
#ifdef MPI
  ! include mpi definitions
  use mpi
#endif
  implicit none
  !
  !  Input/Output variables:
  !
  !  the cluster
  type (cluster), intent(in) :: clust
  !  grid related data
  type (grid_data), intent(in) :: grid
  !  pseudo_potential related data
  type (pseudo_potential), intent(in) :: p_pot
  !  non local pseudo_potential related data
  type (nonloc_pseudo_potential), intent(inout) :: nloc_p_pot
  !  periodic boundary conditions data
  type (pbc_data), intent(inout) :: pbc
  !  parallel computation related data
  type (parallel_data), intent(in) :: parallel
  !  error flag, 300 < ierr < 311
  integer, intent(out) :: ierr
  integer, dimension(1):: ierrvec
  !
  !  Work variables:
  !
  !   counters for k-points
  integer kplp, kpnum

  !  temporary storage variables
  integer npoint,npt
  real(dp) :: rr(3), ylm(7), ylmd(3,7), rnorm, an_tmp(7), vylmd_tmp(3,7)
  !  double grid parameters
  integer idgx, idgy, idgz, ndouble
  real(dp), allocatable :: fdouble(:)
  real(dp) :: invd, ar
  !  true if this grid point is inside the non-local sphere
  logical :: lnloc
  ! partial openmp work to make this faster
  integer :: lnloc_count

  !  true if this grid point is inside the Wigner-Seitz sphere
  logical :: wslnloc
  ! partial openmp work to make this faster
  integer :: wslnloc_count

  !  inloc replaces rnloc and ratom
  integer, dimension(:,:,:), allocatable :: inloc

  real(dp) :: fact1,vwr,dvwr,x,y,z,xa,ya,za
  real(dp), dimension(clust%type_num) :: rc, ainv, binv, rrmin
  !  number of replicas of the periodic cell to be used in the
  !  construction of non-local spheres around each atom
  !  nreplica = 0 if .not. pbc%is_on
  integer nreplica(3)
  !  temporary counter of number of non-local points around each
  !  atom
  integer nloc, wsnloc
  !  loop counters 
  integer i,iat,ity,ja,lm,lp,mg,icellx,icelly,icellz,im
  !  transformation of non-local arrays
  integer ioffset,itot,icur,icount,irow,itran,wsirow,wsitran,wsicur,wsitot,wsicount
  !  local aliases of variables contained in structures
  integer maxnloc,natom,nnodes,comm,mpinfo,wsmaxnloc
  !  counter for the number of non-local projectors
  integer ilm

  real(dp) :: rr_tmp
  real(dp) :: uvec(3), xvec(3), xholdvec(3)
  real(dp) :: tmpmat(3,3), tmpnorm
  !  spin-orbit vriables
  real(dp) :: vwr_ion,vwr_so
  complex(dpc) :: dvr_so, dvr_ion
  complex(dpc), dimension(5) :: zylm, an_tmp_ion, an_tmp_so
  complex(dpc), dimension(3,5) :: zylmd, zvylmd_tmpi, zvylmd_tmps
  integer isy,solm,iso,sol

  !---------------------------------------------------------------

  kpnum = max(pbc%nkpt,1)
  ierr = 0

  natom = clust%atom_num
  nnodes = parallel%procs_num
  comm = parallel%comm
  !
  ! Initialize spin-orbit data.
  !
  nloc_p_pot%is_so = p_pot%is_so
  !
  !  Initialize double grid data.
  !
  ndouble = grid%ndouble
  allocate(fdouble(-ndouble+1:ndouble-1))
  invd = one/real(ndouble,dp)
  do idgx = -ndouble + 1, ndouble - 1
     fdouble(idgx) = (one - real(abs(idgx),dp)*invd)*invd
  enddo

  !  Define parameters for periodic boundary conditions.
  nreplica = 0
  if (pbc%per > 0) nreplica(1:pbc%per) = 1
  !
  !  Determine the spatial extent of the non-local contribution on
  !  the radial grid and translate it to the size of the non-local
  !  block around each atom.
  !
  !  Initialize the atom counter, ja, and arrays in nloc_p_pot structure.
  ja = 0
  nloc_p_pot%wsnlatom(:) = 0
  nloc_p_pot%nlatom(:) = 0
  nloc_p_pot%nlmatm(:) = 0
  nloc_p_pot%skbi(:) = zero
  nloc_p_pot%so_num = 0
  !  For each atom type 
  do ity = 1, clust%type_num
     !  Update temporary variables for number of points in
     !  pseudopotential, npoint, (inverse) radial grid parameters, ainv
     !  and binv, smallest distance from the atom recorded in
     !  psedupotential file (rrmin) and the smaller between the core
     !  -cutoff radius and the largest radius still available in the
     !  pseudopotential file, rc. 
     npoint = p_pot%ns(ity)
     ainv(ity) = one / p_pot%par_a(ity)
     binv(ity) = one / p_pot%par_b(ity)
     rrmin(ity) = p_pot%rs(2, ity)
     rc(ity) = p_pot%rcore(ity)
     if (p_pot%rs(npoint-1, ity) < rc(ity)) then
        rc(ity) = p_pot%rs(npoint-1, ity)
     endif
     !  Initialize matrix for geometrical grid transformations
     !  in case of periodic boundary conditions.
     tmpmat = pbc%latt_vec

     do i = 1, 3
        tmpnorm = sqrt(sum(tmpmat(:,i)**2))
        tmpmat(:,i) = tmpmat(:,i)/tmpnorm
     end do
     !
     !  For each atom inside a given type...
     !
     !  Update atom # counter, ja, and temporary variables for storing
     !  the atomic position, xa, ya, za.
     do iat = 1, clust%natmi(ity)
        ja = ja + 1
        !  Setting a boolean vector to indicate if the atom has
        !  spin-orbit correction psp.
        if(nloc_p_pot%is_so) then
           nloc_p_pot%so(ja) = p_pot%so(ity)
           if(p_pot%so(ity)) nloc_p_pot%so_num = nloc_p_pot%so_num + 1
        else
           nloc_p_pot%so(ja) = .false.
        endif
        !  Distribute atoms round-robin around processors.
        if (mod(ja-1,nnodes) /= parallel%iam) cycle
        xa = clust%xatm(ja)
        ya = clust%yatm(ja)
        za = clust%zatm(ja)
        !
        !  Check if the grid point is within the core-cutoff radius around
        !  the atom. If it is not, it has no non-local contribution because
        !  outside the core all pseudopotential components are equal.
        !  If we are dealing with periodic boundary conditions, we
        !  also need to see whether a given grid point is inside the
        !  core-cutoff radius if shifted by the size of the box in the
        !  +x,-x,+y,-y,+z,-z directions (and any combination thereof). In
        !  theory, one would need to check for shifts of twice the size of
        !  the box, three times the size of the box, etc. In practice this
        !  is NOT needed, because it would imply a core cutoff radius
        !  larger than the box - which makes for a very bad calculation
        !  anyway.

        !  Reset counter of non-local points per given atom.
        nloc = 0
        wsnloc = 0
        !  For each grid point in the full grid...
!$OMP PARALLEL DO &
!$OMP& SCHEDULE(RUNTIME)  &
!$OMP& DEFAULT(SHARED) &
!$OMP& PRIVATE(mg,icellx,icelly,icellz,uvec,x,y,z,xvec,rr_tmp,lnloc_count,wslnloc_count) &
!$OMP& REDUCTION(+:nloc,wsnloc)
        do mg = 1, grid%ndim
           !  Allow each point to be shifted by the length of the box in each
           !  direction.
           do icellx = -nreplica(1),nreplica(1)
              do icelly = -nreplica(2),nreplica(2)
                 do icellz = -nreplica(3),nreplica(3)
                    !  Compute distance between point (or its replica if PBC)
                    !  and atom and increase nloc if it is inside the
                    !  non-local range.
                    !wslnloc = .false.
                    !lnloc = .false.
                    wslnloc_count = 0
                    lnloc_count = 0
                    do idgx = -ndouble + 1, ndouble - 1
                       do idgy = -ndouble + 1, ndouble - 1
                          do idgz = -ndouble + 1, ndouble - 1
                            !no need to keep looking if this is a nonloc point
                            if((wslnloc_count == 0 ) .AND. ( lnloc_count == 0)) then
                             uvec(1)=(grid%shift(1) + grid%fx(mg))* &
                                  grid%step(1)+ &
                                  real(icellx,dp)*pbc%box_size(1)
                             uvec(2)=(grid%shift(2) + grid%fy(mg))* &
                                  grid%step(2)+ &
                                  real(icelly,dp)*pbc%box_size(2)
                             uvec(3)=(grid%shift(3) + grid%fz(mg))* &
                                  grid%step(3)+ &
                                  real(icellz,dp)*pbc%box_size(3)

                             uvec(1) = uvec(1) &
                                  + real(idgx,dp)*grid%step(1)*invd
                             uvec(2) = uvec(2) &
                                  + real(idgy,dp)*grid%step(2)*invd
                             uvec(3) = uvec(3) &
                                  + real(idgz,dp)*grid%step(3)*invd

                             call matvec3('N',pbc%avec_norm,uvec,xvec)

                             x = xvec(1) - xa
                             y = xvec(2) - ya
                             z = xvec(3) - za
                             rr_tmp = sqrt(x*x+y*y+z*z)
                             !if (rr_tmp <= rc(ity)) lnloc = .true.
                             if (rr_tmp <= rc(ity)) lnloc_count = lnloc_count+1
                             !if (rr_tmp <= p_pot%rws(ity)) wslnloc = .true.
                             if (rr_tmp <= p_pot%rws(ity)) wslnloc_count = wslnloc_count+1
                            endif
                          enddo
                       enddo
                    enddo
                    if (lnloc_count > 0 ) nloc = nloc + 1
                    if (wslnloc_count > 0  ) wsnloc = wsnloc + 1
                 enddo
              enddo
           enddo
        enddo
!$OMP END PARALLEL DO
        !  Update array of total number of non-local points per atom and
        !  maximum size of non-local block associated with any atom.
  !DEBUG:
  !write(9,*)'atom, nloc,wsnloc ', ja,nloc,wsnloc
  
        nloc_p_pot%nlatom(ja) = nloc
        nloc_p_pot%wsnlatom(ja) = wsnloc
     enddo                  ! iat = 1, clust%natmi(ity)
  enddo                     ! ity = 1, clust%type_num
  call pisum(nloc_p_pot%wsnlatom,natom,nnodes,comm)
  call pisum(nloc_p_pot%nlatom,natom,nnodes,comm)

#ifdef AJB_DEBUG
     if (parallel%iammaster) then
     write(7,*)
     write(7,*) ' now calling nonloc_pseudo_pot_set_maxnloc'
     write(7,*)
     endif
#endif
  call nonloc_pseudo_pot_set_maxnloc (nloc_p_pot)
#ifdef AJB_DEBUG
     if (parallel%iammaster) then
     write(7,*)
     write(7,*) ' finished nonloc_pseudo_pot_set_maxnloc'
     write(7,*)
     endif
#endif
  maxnloc = nloc_p_pot%maxnloc
  wsmaxnloc = nloc_p_pot%wsmaxnloc
  allocate(inloc(4,maxnloc,natom))
  inloc = 0

  ja = 0
  iso = 0
  isy = 0
  !  For each atom type
  do ity = 1, clust%type_num
     if (p_pot%so(ity)) isy = isy + 1
     do iat = 1, clust%natmi(ity)
        ja = ja + 1
        !  Set some identifying arrays for atoms with spin-orbit coupling.
        if(nloc_p_pot%is_so) then
           if (nloc_p_pot%so(ja)) then
              iso = iso + 1
              nloc_p_pot%so_indx(iso) = ja
              nloc_p_pot%cc(:,iso) = p_pot%cc(:,isy)
           endif
        endif

        !  Distribute atoms round-robin around processors.
        if (mod(ja-1,nnodes) /= parallel%iam) cycle
        xa = clust%xatm(ja)
        ya = clust%yatm(ja)
        za = clust%zatm(ja)
        nloc = 0
        wsnloc = 0
        do mg = 1, grid%ndim
           do icellx = -nreplica(1),nreplica(1)
              do icelly = -nreplica(2),nreplica(2)
                 do icellz = -nreplica(3),nreplica(3)
  !  Compute distance between point (or its replica if PBC) and atom.
  !  If the point is inside non-local range:
  !  1. Increase the non-local point counter, nloc, by one.
  !  2. Update the output arrays, indw, inloc, with the grid point, (in
  !     the irreducible wedge) and its location on the cell, respectively.

                    lnloc = .false.
                    wslnloc = .false.
                    do idgx = -ndouble + 1, ndouble - 1
                       do idgy = -ndouble + 1, ndouble - 1
                          do idgz = -ndouble + 1, ndouble - 1
                          !if this point is flagged, no need to keep comparing
                            if ( (.NOT. lnloc) .AND. (.NOT. wslnloc) ) then

                             uvec(1)=(grid%shift(1) + grid%fx(mg))* &
                                  grid%step(1)+ &
                                  real(icellx,dp)*pbc%box_size(1)
                             uvec(2)=(grid%shift(2) + grid%fy(mg))* &
                                  grid%step(2)+ &
                                  real(icelly,dp)*pbc%box_size(2)
                             uvec(3)=(grid%shift(3) + grid%fz(mg))* &
                                  grid%step(3)+ &
                                  real(icellz,dp)*pbc%box_size(3)

                             uvec(1) = uvec(1) &
                                  + real(idgx,dp)*grid%step(1)*invd
                             uvec(2) = uvec(2) &
                                  + real(idgy,dp)*grid%step(2)*invd
                             uvec(3) = uvec(3) &
                                  + real(idgz,dp)*grid%step(3)*invd

                             call matvec3('N',pbc%avec_norm,uvec,xvec)

                             x = xvec(1) - xa
                             y = xvec(2) - ya
                             z = xvec(3) - za
                             
                             if ( (idgx == 0) & 
                            .and. (idgy == 0) &  
                            .and. (idgz == 0) ) then
                                xholdvec(1)=xvec(1)
                                xholdvec(2)=xvec(2)
                                xholdvec(3)=xvec(3)
                             endif
                             
                             rr_tmp = sqrt(x*x+y*y+z*z)
                             if (rr_tmp <= rc(ity)) lnloc = .true.
                             !if (rr_tmp <= rc(ity)) lnloc_count = lnloc_count + 1
                             if (rr_tmp <= p_pot%rws(ity)) wslnloc = .true.
                             !if (rr_tmp <= p_pot%rws(ity)) wslnloc_count = wslnloc_count +1
                            endif
                          enddo
                       enddo
                    enddo
                    if (wslnloc) then
                       wsnloc = wsnloc + 1
                       nloc_p_pot%wsindw(wsnloc,ja) = grid%rindex(mg)
                       nloc_p_pot%wstran(wsnloc,ja) = grid%rtrans(mg)
                    endif

                    if (lnloc) then
                       nloc = nloc + 1
                       nloc_p_pot%indw(nloc,ja) = grid%rindex(mg)
                       nloc_p_pot%tran(nloc,ja) = grid%rtrans(mg)
                       inloc(:,nloc,ja) = (/ mg, icellx, icelly, icellz /)

                       if(pbc%nkpt > 0) then
                          do kplp = 1, kpnum

                           nloc_p_pot%right(nloc,kplp,ja)= &
                   exp( zi * (pbc%kpts(1,kplp)*xholdvec(1)+ &
                              pbc%kpts(2,kplp)*xholdvec(2)+ &
                              pbc%kpts(3,kplp)*xholdvec(3)) )
                           
                           nloc_p_pot%left(nloc,kplp,ja)= &
                   exp(-1*zi*(pbc%kpts(1,kplp)*xholdvec(1)+ &
                              pbc%kpts(2,kplp)*xholdvec(2)+ &
                              pbc%kpts(3,kplp)*xholdvec(3)) )
                          enddo
                       endif

                    endif
                 enddo
              enddo
           enddo
        enddo
        if (nloc /= nloc_p_pot%nlatom(ja)) then
           write(9,*) ' ERROR in nonloc ', &
                nloc,iat,ja,ity,nloc_p_pot%nlatom(ja)
           ierr = 301
        endif

        if (wsnloc /= nloc_p_pot%wsnlatom(ja)) then
           write(9,*) ' ERROR in nonloc ', &
                wsnloc,iat,ja,ity,nloc_p_pot%wsnlatom(ja)
           ierr = 555
        endif
     enddo                  ! iat = 1, clust%natmi(ity)
  enddo                     ! ity = 1, clust%type_num
  ierrvec = ierr
  call pisum(ierrvec,1,nnodes,comm)
  ierr = ierrvec(1)
  if (ierr /= 0) then
     ierr = 301
     if (parallel%iammaster) write(9,*) ' ERROR in nonloc '
     return
  endif
  !
  !  Report non-local block sizes to file.
  !
  write(9,*)
  write(9,*) 'Non-local pseudopotential messages:'
  write(9,*) '-----------------------------------'
  write(9,*) 'Max # of nonlocal points for one atom = ',nloc_p_pot%maxnloc
  write(9,*)
  write(9,*) 'Sizes of all non-local blocks '
  write(9,'(12(1x,i5))') (nloc_p_pot%nlatom(i),i=1,clust%atom_num)
  if (nloc_p_pot%is_so) write(9,*) &
       'Tot. number of atoms with spin-orbit term =',nloc_p_pot%so_num
  write(9,*)
  if (parallel%iammaster) then
     write(7,*)
     write(7,*) 'Non-local pseudopotential messages:'
     write(7,*) '-----------------------------------'
     write(7,*) 'Max # of nonlocal points for one atom = ',nloc_p_pot%maxnloc
     write(7,*)
     write(7,*) 'Sizes of all non-local blocks '
     write(7,'(12(1x,i5))') (nloc_p_pot%nlatom(i),i=1,clust%atom_num)
     if (nloc_p_pot%is_so) write(7,*) &
          'Tot. number of atoms with spin-orbit term =',nloc_p_pot%so_num
     write(7,*)
  endif
  !
  !  Calculate the non-local potential vector for each angular
  !  component. The result is stored in anloc.
  !  The non-local pseudopotential to the hamiltonian is
  !  subsequently calculated as the outer product of these vectors,
  !  summed over all angular components l,m.
  !
  ja = 0
  !  Initialize the counter of non-local projectors.
  ilm = 0
  iso = 0
  isy = 0

  if (nloc_p_pot%is_so) then
     nloc_p_pot%v_so(:,:,:) = zzero
     nloc_p_pot%dv_so(:,:,:,:) = zzero
     nloc_p_pot%v_ion(:,:,:) = zzero
     nloc_p_pot%dv_ion(:,:,:,:) = zzero
  endif
  !
  !  For each atom type, for each atom within each type
  do ity = 1, clust%type_num
     npoint = p_pot%ns(ity)
     if (p_pot%so(ity)) isy = isy + 1
     do iat = 1, clust%natmi(ity)
        !  Again update temporary storage of atom # and position
        ja = ja + 1

        if (nloc_p_pot%so(ja)) iso = iso + 1

        !  Distribute atoms round-robin around processors
        if (mod(ja-1,nnodes) /= parallel%iam) then
           do lp = 1, p_pot%nlocp(ity)
              if (lp /= p_pot%loc(ity)) ilm = ilm + 2*lp - 1
           enddo
           cycle
        endif

        nloc = nloc_p_pot%nlatom(ja)
        !  Initialize angular momentum (l,m) count
        lm = 0
        !
        !  Loop over each angular momentum component, l (and within
        !  the loop calculate for different m's too).
        ! OMP DO REDUCTION(+:lm) from here...
        do lp = 1, p_pot%nlocp(ity)
           !  Skip the local component of the pseudopotential!
           if (lp == p_pot%loc(ity) .and. lp == 1) cycle
           !  Calculate the multiplicand - (non-local pseudopotential)*
           !  (radial pseudowavefunction)*(angular wave function), put the
           !  results in anloc. That is, calculate explicitly |Vl-Vloc|phi_lm>.

           do i = 1, nloc
              an_tmp(:) = zero
              vylmd_tmp(:,:) = zero

              an_tmp_ion(:) = zzero
              zvylmd_tmpi(:,:) = zzero

              an_tmp_so(:) = zzero
              zvylmd_tmps(:,:) = zzero

              mg = inloc(1,i,ja)
              icellx = inloc(2,i,ja)
              icelly = inloc(3,i,ja)
              icellz = inloc(4,i,ja)
              do idgx = -ndouble + 1, ndouble - 1
                 do idgy = -ndouble + 1, ndouble - 1
                    do idgz = -ndouble + 1, ndouble - 1
                       uvec(1)=(grid%shift(1) + grid%fx(mg))* &
                            grid%step(1) + real(icellx,dp)*pbc%box_size(1)
                       uvec(2)=(grid%shift(2) + grid%fy(mg))* &
                            grid%step(2) + real(icelly,dp)*pbc%box_size(2)
                       uvec(3)=(grid%shift(3) + grid%fz(mg))* &
                            grid%step(3) + real(icellz,dp)*pbc%box_size(3)
                       uvec(1) = uvec(1) + real(idgx,dp)*grid%step(1)*invd
                       uvec(2) = uvec(2) + real(idgy,dp)*grid%step(2)*invd
                       uvec(3) = uvec(3) + real(idgz,dp)*grid%step(3)*invd

                       call matvec3('N',pbc%avec_norm,uvec,xvec)
                       rr(1) = xvec(1) - clust%xatm(ja)
                       rr(2) = xvec(2) - clust%yatm(ja)
                       rr(3) = xvec(3) - clust%zatm(ja)

                       rnorm = sqrt(dot_product(rr,rr))
                       !  In order to avoid dividing by zero:
                       if (rnorm < rrmin(ity))  then
                          rnorm = rrmin(ity)
                       endif
                       if (rnorm > rc(ity)) cycle

                       if(lp /= p_pot%loc(ity)) then
                          npt = idint(binv(ity) * log(one+ainv(ity)*rnorm))+1
                          fact1 = (rnorm - p_pot%rs(npt,ity)) /  &
                               (p_pot%rs(npt+1,ity) - p_pot%rs(npt,ity))

                          vwr  = p_pot%vw(npt, ity, lp) + fact1* &
                               (p_pot%vw(npt+1,ity,lp)-p_pot%vw(npt,ity,lp))
                          dvwr  = p_pot%dvw(npt, ity, lp) + fact1* &
                               (p_pot%dvw(npt+1,ity,lp)-p_pot%dvw(npt,ity,lp))
                          if (p_pot%uspline(ity)) then
                             call splint(p_pot%rs(-p_pot%norder:npoint,ity), &
                                         p_pot%vw(-p_pot%norder:npoint,ity,lp), &
                                         p_pot%d2vw(-p_pot%norder:npoint,ity,lp), &
                                         npoint+1+p_pot%norder,rnorm,vwr,dvwr)
                          endif

                          vwr = vwr * fdouble(idgx)*fdouble(idgy)*fdouble(idgz)
                          dvwr = dvwr * fdouble(idgx)*fdouble(idgy)*fdouble(idgz)

                          call ylm_real(lp,rr,rnorm,ylm,ylmd)
                          do im = 1, 2*lp - 1
                             an_tmp(im) = an_tmp(im) + ylm(im)*vwr
                             vylmd_tmp(:,im) = vylmd_tmp(:,im) + &
                                  ylmd(:,im)*vwr + dvwr*ylm(im)*rr/rnorm
                          enddo
                       endif

                       !   Do the same for the spin-orbit potentials, this
                       !   time with complex arrays here s orbitals have no
                       !   contribution to spin-orbit and local component has
                       !   a contribution.
                       if(nloc_p_pot%is_so .and. lp > 1) then 
                          if (nloc_p_pot%so(ja)) then
                             sol = lp - 1
                             vwr_ion = p_pot%vionr(npt,isy,sol)+ &
                                  fact1*(p_pot%vionr(npt+1,isy,sol) &
                                  -p_pot%vionr(npt,isy,sol))

                             vwr_so = p_pot%vsor(npt,isy,sol)+ &
                                  fact1*(p_pot%vsor(npt+1,isy,sol) &
                                  -p_pot%vsor(npt,isy,sol))

                             dvr_so  = p_pot%dvr_so(npt,isy,sol) + fact1* &
                               (p_pot%dvr_so(npt+1,isy,sol)-p_pot%dvr_so(npt,isy,sol))

                             dvr_ion  = p_pot%dvr_ion(npt,isy,sol) + fact1* &
                               (p_pot%dvr_ion(npt+1,isy,sol)-p_pot%dvr_ion(npt,isy,sol))

                             if (p_pot%uspline(ity)) then
                                call splint(p_pot%rs(-p_pot%norder:npoint,ity), &
                                            p_pot%vionr(-p_pot%norder:npoint,isy,sol), &
                                            p_pot%d2vr_ion(-p_pot%norder:npoint,isy,sol), &
                                            npoint+1+p_pot%norder,rnorm,vwr_ion,dvr_ion)
                                call splint(p_pot%rs(-p_pot%norder:npoint,ity), &
                                            p_pot%vsor(-p_pot%norder:npoint,isy,sol), &
                                            p_pot%d2vr_so(-p_pot%norder:npoint,isy,sol), &
                                            npoint+1+p_pot%norder,rnorm,vwr_so,dvr_so)
                             endif

                             vwr_ion = vwr_ion &
                                  *fdouble(idgx)*fdouble(idgy)*fdouble(idgz)

                             vwr_so = vwr_so &
                                  *fdouble(idgx)*fdouble(idgy)*fdouble(idgz)

                             dvr_so = dvr_so * fdouble(idgx)*fdouble(idgy)*fdouble(idgz)


                             dvr_ion = dvr_ion * fdouble(idgx)*fdouble(idgy)*fdouble(idgz)


                             call ylm_cplx(sol,rr,rnorm,zylm,zylmd)
                             rr = rr / rnorm

                             do im = 1, 2*sol + 1
                                an_tmp_ion(im) = an_tmp_ion(im) + zylm(im) * &
                                     vwr_ion

                                an_tmp_so(im) = an_tmp_so(im) + zylm(im) * &
                                     vwr_so

                                zvylmd_tmps(:,im) = zvylmd_tmps(:,im) + &
                                  zylmd(:,im)*vwr_so + dvr_so*zylm(im)*rr

                                zvylmd_tmpi(:,im) = zvylmd_tmpi(:,im) + &
                                  zylmd(:,im)*vwr_ion + dvr_ion*zylm(im)*rr

                             enddo
                          endif !if(p_pot%spinorbit .and. ... .ne....) then
                       endif !if (nloc_p_pot%so(ja))then

                    enddo !idgx = -ndouble + 1, ndouble - 1
                 enddo !idgy = -ndouble + 1, ndouble - 1
              enddo !idgz = -ndouble + 1, ndouble - 1

              if(nloc_p_pot%is_so .and. lp > 1) then 
                 if (nloc_p_pot%so(ja)) then
                   sol = lp - 1
                   do im = 1, 2*sol + 1
                      solm = 3*(sol - 1) + im
                      nloc_p_pot%v_so(i,solm,iso) = an_tmp_so(im)
                      nloc_p_pot%v_ion(i,solm,iso) = an_tmp_ion(im)
                      nloc_p_pot%dv_so(:,i,solm,iso) = zvylmd_tmps(:,im)
                      nloc_p_pot%dv_ion(:,i,solm,iso) = zvylmd_tmpi(:,im)
                   enddo
                endif
              endif

              if (lp /= p_pot%loc(ity)) then
                 do im = 1, 2*lp - 1
                    nloc_p_pot%anloc(i, ilm+im) = an_tmp(im)
                    nloc_p_pot%vylmd(:, i, ilm+im) = vylmd_tmp(:,im)
                 enddo
              endif

           enddo  ! do i = 1, nloc
           if (lp == p_pot%loc(ity)) cycle
           do im = 1, 2*lp - 1
              nloc_p_pot%skbi(ilm+im) = p_pot%ekbi(lp,ity)
           enddo

           lm = lm + 2*lp - 1
           ilm = ilm + 2*lp - 1

        enddo               ! lp = 1, p_pot%nlocp(ity)
        !  Update array with number of lm components for each atom.
        nloc_p_pot%nlmatm(ja) = lm

     enddo                  ! iat = 1, clust%natmi(ity)
  enddo                     ! ity = 1, clust%type_num

  deallocate(inloc)
  deallocate(fdouble)
  !
  !  Share information. Follow the round-robin distribution of atoms.
  !
  call pisum(nloc_p_pot%indw,maxnloc*nloc_p_pot%atom_num,nnodes,comm)
  call pisum(nloc_p_pot%tran,maxnloc*nloc_p_pot%atom_num,nnodes,comm)
  if(wsmaxnloc>0) then
  call pisum(nloc_p_pot%wsindw,wsmaxnloc*nloc_p_pot%atom_num,nnodes,comm)
  call pisum(nloc_p_pot%wstran,wsmaxnloc*nloc_p_pot%atom_num,nnodes,comm)
  endif
  call pisum(nloc_p_pot%nlmatm,nloc_p_pot%atom_num,nnodes,comm)
  call psum(nloc_p_pot%anloc,maxnloc*nloc_p_pot%nlm,nnodes,comm)
  call psum(nloc_p_pot%vylmd,3*maxnloc*nloc_p_pot%nlm,nnodes,comm)
  call psum(nloc_p_pot%skbi,nloc_p_pot%nlm,nnodes,comm)
  if (pbc%nkpt > 0) then
     call zpsum(nloc_p_pot%right,maxnloc*pbc%nkpt*nloc_p_pot%atom_num, &
          nnodes,comm)
     call zpsum(nloc_p_pot%left,maxnloc*pbc%nkpt*nloc_p_pot%atom_num, &
          nnodes,comm)
  endif

#ifdef MPI
  do iso = 1, nloc_p_pot%so_num
     ja = nloc_p_pot%so_indx(iso)
     i = mod(ja-1,nnodes)
     call MPI_BCAST(nloc_p_pot%v_so(1,1,iso),maxnloc*8,   MPI_DOUBLE_COMPLEX,i,comm,mpinfo)
     call MPI_BCAST(nloc_p_pot%v_ion(1,1,iso),maxnloc*8,  MPI_DOUBLE_COMPLEX,i,comm,mpinfo)
     call MPI_BCAST(nloc_p_pot%cc(1,iso),2,  MPI_DOUBLE_PRECISION,i,comm,mpinfo)
     call MPI_BCAST(nloc_p_pot%dv_so(1,1,1,iso),maxnloc*24,  MPI_DOUBLE_COMPLEX,i,comm,mpinfo)
     call MPI_BCAST(nloc_p_pot%dv_ion(1,1,1,iso),maxnloc*24, MPI_DOUBLE_COMPLEX,i,comm,mpinfo)
  enddo
#endif
#ifdef OLDDEBUG
  if (nloc_p_pot%so_num > 0) nloc_p_pot%cc = nloc_p_pot%cc*1.d1
#endif
  !
  !  -------- Transformation of the non local arrays -------------
  !  Transform the global "non-local" arrays into arrays containing
  !  information about the local rows of the processor.i.e.,
  !  For all atoms on my processor:
  !  a) all rows of Anloc that are not my rows are removed.
  !  b) my rows are shifted in adjacent positions.
  !  c) indw is changed to include only the above local information
  !     The offset is subtracted to yield the local row locations
  !  d) nlatom contains the number of my rows in Anloc
  !  e) nlmatm is the same since it does not depend on rows:
  !
  write(9,*) 'Transform atoms', nloc_p_pot%atom_num
  ioffset = parallel%irows(parallel%group_iam) - 1
  icount = 0
  itot = 0
  iso = 0
  ilm = 0
  do ja = 1, nloc_p_pot%atom_num
     if(nloc_p_pot%so(ja)) iso = iso + 1     
     !  ..for all atoms
     icur = 0
     nloc = nloc_p_pot%nlatom(ja)
     do i = 1, nloc
        !  ..for each row in Anloc
        irow = nloc_p_pot%indw(i,ja)
        itran = nloc_p_pot%tran(i,ja)

        if (parallel%irows(parallel%group_iam) <= irow.and.irow < &
             parallel%irows(parallel%group_iam+1)) then
           !  ..if irow on my processor shift it up
           icur = icur + 1
           do lm = 1, nloc_p_pot%nlmatm(ja)
              nloc_p_pot%anloc(icur,ilm+lm) = nloc_p_pot%anloc(i,ilm+lm)
              nloc_p_pot%vylmd(:,icur,ilm+lm) = nloc_p_pot%vylmd(:,i,ilm+lm)
           enddo
           do lm = 1, 8
              if(nloc_p_pot%so(ja)) then
                 nloc_p_pot%v_so(icur,lm,iso) = nloc_p_pot%v_so(i,lm,iso)
                 nloc_p_pot%v_ion(icur,lm,iso) = nloc_p_pot%v_ion(i,lm,iso)
                 nloc_p_pot%dv_so(:,icur,lm,iso) = nloc_p_pot%dv_so(:,i,lm,iso)
                 nloc_p_pot%dv_ion(:,icur,lm,iso) = nloc_p_pot%dv_ion(:,i,lm,iso)
              endif
           enddo
           if (pbc%nkpt > 0) then
              nloc_p_pot%right(icur,:,ja) = nloc_p_pot%right(i,:,ja)
              nloc_p_pot%left(icur,:,ja) = nloc_p_pot%left(i,:,ja)
           endif
           nloc_p_pot%indw(icur,ja) = irow - ioffset
           nloc_p_pot%tran(icur,ja) = itran
        endif

     enddo

     ilm = ilm + nloc_p_pot%nlmatm(ja)
     !  ..new number of rows in Anloc for this atom.
     nloc_p_pot%nlatom(ja) = icur
     icount = icount + icur
     itot   = itot + nloc
  enddo
  write(9,*) 'My number of non-local rows: ',icount,' of total ',itot

  !
  ! do the same to Wigner-Seitz spheres ...
  wsicount = 0
  wsitot = 0
  do ja = 1, nloc_p_pot%atom_num
     wsicur = 0
     wsnloc = nloc_p_pot%wsnlatom(ja)
     do i = 1, wsnloc 
        !  ..for each row 
        wsirow = nloc_p_pot%wsindw(i,ja)
        wsitran = nloc_p_pot%wstran(i,ja)
        if (parallel%irows(parallel%group_iam) <= wsirow .and. wsirow < &
             parallel%irows(parallel%group_iam+1)) then
           wsicur = wsicur + 1
           nloc_p_pot%wsindw(wsicur,ja) = wsirow - ioffset
           nloc_p_pot%wstran(wsicur,ja) = wsitran
        endif
     enddo
     nloc_p_pot%wsnlatom(ja) = wsicur
     wsicount = wsicount + wsicur
     wsitot   = wsitot + wsnloc
  enddo
  write(9,*) 'My number of rows in WS sphere: ',wsicount,' of total ',wsitot
end subroutine nonloc
!===============================================================
subroutine ylm_real(lp,rr,rnorm,ylm,ylmd)
!
!  For an input point with cartesian coordinates rr(1:3),
!  calculate the real combination of spherical harmonics
!  ylm(1:2*lp-1) and their gradients, ylmd(1:3, 1:2*lp-1).
!
!---------------------------------------------------------------
  use constants
  implicit none 
  !
  !  Input variables:
  !
  integer, intent(in) :: lp ! l + 1, orbital quantum number
  real(dp), intent(in) :: rr(3),rnorm
  real(dp), intent(out) :: ylm(7),ylmd(3,7)
  !
  !  Work variables:
  !
  !  coordinate factors
  real(dp) :: xx, yy, zz, x2, y2, z2, xy, xz, yz, xyz, rinv, thrinv
  !
  !  Spherical harmonics coefficients - for s, p, and d
  !  formulae, see, e.g., pp. 11-12 in W. A. Harrison, "Electronic
  !  Structure and the Properties of Solids: The Physics of the
  !  Chemical Bond". For f formulae, see, e.g., 
  !  S. F. A. Kettle, "Physical Inorganic Chemistry", p. 244.
  !  Notice that they *have* the factor sqrt( (2*l+1)/4*pi ).
  !
  !  1/(2*sqrt(pi)) for s
  real(dp), parameter :: c0  = 0.282094791773878143474039725779d0
  !  sqrt(3/(4*pi)), for px, py, pz 
  real(dp), parameter :: c1  = 0.488602511902919921586384622837d0
  !  sqrt(15/(4*pi)), dxy, dzy, dzx 
  real(dp), parameter :: c21 = 1.092548430592079070543385705802d0
  !  sqrt(5/(16*pi)), for dz^2 
  real(dp), parameter :: c22 = 0.315391565252520006030893690294d0
  !  sqrt(7/(16*pi)), for f_z(5z^2-3r^2)
  real(dp), parameter :: c31 = 0.373176332590115391414395913197d0
  !  sqrt(21/(32*pi)), for f_x(5z^2-r^2),f_y(5z^2-r^2)
  real(dp), parameter :: c32 = 0.457045799464465736158020696916d0
  !  sqrt(105/(16*pi)), for f_z(x^2-y^2)
  real(dp), parameter :: c33 = 1.445305721320277027694690077198d0
  !  sqrt(105/(4*pi)), for f_xyz
  real(dp), parameter :: c34 = 2.890611442640554055389380154398d0
  !  sqrt(35/(32*pi)), for f_x(x^2-3y^2),f_y(y^2-3x^2)
  real(dp), parameter :: c35 = 0.590043589926643510345610277541d0

  !---------------------------------------------------------------

  ylm(:) = zero
  ylmd(:,:) = zero
  rinv = one/rnorm
  xx = rr(1) * rinv
  yy = rr(2) * rinv
  zz = rr(3) * rinv

  select case(lp)
  case (1)
     !
     !  l = 0 (lp = 1) , s orbital
     !
     ylm(1) = c0
  case (2)
     !
     !  l = 1 (lp = 2), p orbitals
     !
     !  --- px, py, pz
     ylm(1) = c1*xx
     ylm(2) = c1*yy
     ylm(3) = c1*zz
     !  --- derivatives of px,py,pz
     xy = -rinv*xx*yy
     xz = -rinv*xx*zz
     yz = -rinv*yy*zz

     ylmd(1,1) = c1*rinv*(one-xx*xx) 
     ylmd(2,1) = c1*xy
     ylmd(3,1) = c1*xz

     ylmd(1,2) = c1*xy
     ylmd(2,2) = c1*rinv*(one-yy*yy)
     ylmd(3,2) = c1*yz

     ylmd(1,3) = c1*xz
     ylmd(2,3) = c1*yz
     ylmd(3,3) = c1*rinv*(one-zz*zz)
  case (3)
     !
     !  l = 2 (lp = 3) , d orbitals
     !
     !  --- dxy, dyz, dxz, dz2, dx2-y2 
     ylm(1) = c21*xx*yy
     ylm(2) = c21*zz*yy
     ylm(3) = c21*xx*zz
     ylm(4) = c22*(three*zz*zz-one)
     ylm(5) = c21*half*(xx*xx-yy*yy)
     !  --- derivatives of dxy, dyz, dxz, dz2, dx2-y2 
     x2 = xx*xx
     y2 = yy*yy
     z2 = zz*zz
     xyz = -two*rinv*xx*yy*zz

     ylmd(1,1) = c21*rinv*yy*(one-two*x2) 
     ylmd(2,1) = c21*rinv*xx*(one-two*y2) 
     ylmd(3,1) = c21*xyz

     ylmd(1,2) = c21*xyz
     ylmd(2,2) = c21*rinv*zz*(one-two*y2) 
     ylmd(3,2) = c21*rinv*yy*(one-two*z2) 

     ylmd(1,3) = c21*rinv*zz*(one-two*x2) 
     ylmd(2,3) = c21*xyz
     ylmd(3,3) = c21*rinv*xx*(one-two*z2) 

     ylmd(1,4) = -c22*six*rinv*xx*z2 
     ylmd(2,4) = -c22*six*rinv*yy*z2 
     ylmd(3,4) =  c22*six*rinv*zz*(one-z2) 

     ylmd(1,5) =  c21*rinv*xx*(one-x2+y2) 
     ylmd(2,5) = -c21*rinv*yy*(one+x2-y2) 
     ylmd(3,5) = -c21*rinv*zz*(x2-y2) 
  case (4)
     !
     !  l = 3 (lp = 4) , f orbitals
     !
     !  --- fz(5z^2-3r^2), fx(5z^2-r^2), fy(5z^2-r^2),
     !      fz(x^2-y^2), fxyz, fx(x^2-3y^2), fy(y^2-3x^2)
     ylm(1) = c31*zz*(five*zz*zz-three)
     ylm(2) = c32*xx*(five*zz*zz-one)
     ylm(3) = c32*yy*(five*zz*zz-one)
     ylm(4) = c33*zz*(xx*xx-yy*yy)
     ylm(5) = c34*xx*yy*zz
     ylm(6) = c35*xx*(xx*xx-three*yy*yy)
     ylm(7) = c35*yy*(yy*yy-three*xx*xx)
     !  --- derivatives of the seven f orbitals
     x2 = xx*xx
     y2 = yy*yy
     z2 = zz*zz
     xy = rinv*xx*yy
     xz = rinv*xx*zz
     yz = rinv*yy*zz
     thrinv = three*rinv

     ylmd(1,1) = c31*three*xz*(one-five*z2)
     ylmd(2,1) = c31*three*yz*(one-five*z2)
     ylmd(3,1) = c31*thrinv*(z2-one)*(one-five*z2) 

     ylmd(1,2) = c32*rinv*(five*z2-one+x2*(one-15.d0*z2))
     ylmd(2,2) = c32*xy*(one-15.d0*z2)
     ylmd(3,2) = c32*xz*(11.d0-15.d0*z2)

     ylmd(1,3) = c32*xy*(one-15.d0*z2)
     ylmd(2,3) = c32*rinv*(five*z2-one+y2*(one-15.d0*z2))
     ylmd(3,3) = c32*yz*(11.d0-15.d0*z2)

     ylmd(1,4) =  c33*xz*(two-three*(x2-y2))
     ylmd(2,4) = -c33*yz*(two+three*(x2-y2))
     ylmd(3,4) =  c33*rinv*(one-three*z2)*(x2-y2)

     ylmd(1,5) = c34*yz*(one-three*x2)
     ylmd(2,5) = c34*xz*(one-three*y2)
     ylmd(3,5) = c34*xy*(one-three*z2)

     ylmd(1,6) = c35*thrinv*(x2-y2-x2*(x2-three*y2))
     ylmd(2,6) = c35*three*xy*(three*y2-x2-two)
     ylmd(3,6) = c35*three*xz*(three*y2-x2)

     ylmd(1,7) = c35*three*xy*(three*x2-y2-two)
     ylmd(2,7) = c35*thrinv*(y2-x2-y2*(y2-three*x2))
     ylmd(3,7) = c35*three*yz*(three*x2-y2)
  end select

end subroutine ylm_real
!===============================================================
subroutine ylm_cplx(lp,rr,rnorm,zylm,zylmd)
!
!  For an input point with cartesian coordinates rr(1:3),
!  calculate the real combination of spherical harmonics
!  ylm(1:2*lp-1) and their gradients, ylmd(1:3, 1:2*lp-1).
!
!---------------------------------------------------------------
  use constants
  implicit none
  !
  !  Input variables:
  !
  integer, intent(in) :: lp ! l + 1, orbital quantum number
  real(dp), intent(in) :: rr(3)
  real(dp), intent(in) :: rnorm
  complex(dpc), intent(out) :: zylm(5)
  complex(dpc), intent(out) :: zylmd(3,5)
  !
  !  Work variables:
  !
  !  coordinate factors
  real(dp) :: xx, yy, zz, x2, y2, z2, xy, xz, yz, xyz, rinv
  !
  !  Spherical harmonics coefficients - for s, p, and d
  !  formulae, see, e.g., pp. 11-12 in W. A. Harrison, "Electronic
  !  Structure and the Properties of Solids: The Physics of the
  !  Chemical Bond". For f formulae, see, e.g.,
  !  S. F. A. Kettle, "Physical Inorganic Chemistry", p. 244.
  !  Notice that they *have* the factor sqrt( (2*l+1)/4*pi ).
  !

  !  sqrt(3/(4*pi)), for Y10
  real(dp), parameter :: c10  = 0.488602511902919921586384622838d0
  !  sqrt(3/(8*pi)), for Y11
  real(dp), parameter :: c11 =  0.345494149471335479265244646031d0
  !  sqrt(5/(16*pi)), for Y20
  real(dp), parameter :: c20 =  0.315391565252520006030893690294d0
  !  sqrt(15/(8*pi)), for Y21
  real(dp), parameter :: c21 =  0.772548404046379160684385470622d0
  !  sqrt(15/(32*pi)), for Y22
  real(dp), parameter :: c22 =  0.386274202023189580342192735311d0

  !---------------------------------------------------------------

  zylm(:) = zzero
  zylmd(:,:) = zzero
  rinv = one/rnorm
  xx = rr(1) * rinv
  yy = rr(2) * rinv
  zz = rr(3) * rinv
  xy = xx*yy
  xz = xx*zz
  yz = yy*zz
  x2 = xx*xx
  y2 = yy*yy
  z2 = zz*zz
  xyz = two*xx*yy*zz
  select case(lp)
  case (1)
     !
     !  l = 1, p orbitals
     !
     !  --- Y1m
     zylm(1) = c11*cmplx(xx,-yy)  ! Y(1,-1)
     zylm(2) = c10*cmplx(zz,zero) ! Y(1, 0)
     zylm(3) = c11*cmplx(-xx,-yy) ! Y(1, 1)
     !  --- derivatives of Y1m

     zylmd(1,1) = c11*rinv*cmplx(one-xx*xx,-xy)
     zylmd(2,1) = c11*rinv*cmplx(xy,yy*yy-one)
     zylmd(3,1) = c11*rinv*cmplx(xz,-yz)

     zylmd(1,2) = c10*rinv*cmplx(xz,zero)
     zylmd(2,2) = c10*rinv*cmplx(yz,zero)
     zylmd(3,2) = c10*rinv*cmplx(one-zz*zz,zero)

     zylmd(1,3) = -conjg(zylmd(1,1))
     zylmd(2,3) = -conjg(zylmd(1,2))
     zylmd(3,3) = -conjg(zylmd(1,3))
  case (2)
     !
     !  l = 2 , d orbitals
     !
     !  --- Y2m

     zylm(1) = c22*cmplx(x2-y2,-two*xy)         ! Y(2,-2)
     zylm(2) = c21*cmplx(xz,-yz)                ! Y(2,-1)
     zylm(3) = c20*cmplx(three*z2-one,zero)     ! Y(2. 0)
     zylm(4) = c21*cmplx(-xz,-yz)               ! Y(2, 1)
     zylm(5) = c22*cmplx(x2-y2,two*xy)          ! Y(2, 2)

     !  --- derivatives of Y2m

     zylmd(1,1) = c22*rinv*two*cmplx(xx*(two*y2+z2),yy*(x2-y2-z2))
     zylmd(2,1) = c22*rinv*two*cmplx(-yy*(two*x2+z2),-xx*(x2-y2+z2))
     zylmd(3,1) = c22*rinv*two*cmplx(-zz*(x2-y2),xyz)

     zylmd(1,2) = c21*rinv*cmplx(zz*(x2-y2-z2),xyz)
     zylmd(2,2) = c21*rinv*cmplx(xyz,-zz*(x2-y2+z2))
     zylmd(3,2) = c21*rinv*cmplx(-xx*(x2+y2-z2),-yy*(x2+y2-z2))

     zylmd(1,3) = c20*rinv*cmplx(-six*xx*z2,zero)
     zylmd(2,3) = c20*rinv*cmplx(-six*yy*z2,zero)
     zylmd(3,3) = c20*rinv*cmplx(six*zz*(x2+y2),zero)

     zylmd(1,4) = -conjg(zylmd(1,2))
     zylmd(2,4) = -conjg(zylmd(2,2))
     zylmd(3,4) = -conjg(zylmd(3,2))

     zylmd(1,5) = conjg(zylmd(1,1))
     zylmd(2,5) = conjg(zylmd(2,1))
     zylmd(3,5) = conjg(zylmd(3,1))

  end select

end subroutine ylm_cplx
!===============================================================
