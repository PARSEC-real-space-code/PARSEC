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
subroutine upot(clust,grid,p_pot,u_pot,pbc,parallel,nspin,cplx,ierr)

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
  ! Input/Output variables:
  !
  !  the cluster
  type (cluster), intent(in) :: clust
  !  grid related data
  type (grid_data), intent(in) :: grid
  !  pseudo_potential related data
  type (pseudo_potential), intent(in) :: p_pot
  !  non local pseudo_potential related data
  type (nonloc_pseudo_potential), intent(inout) :: u_pot
  !  periodic boundary conditions data
  type (pbc_data), intent(inout) :: pbc
  !  parallel computation related data
  type (parallel_data), intent(in) :: parallel
  !  number of spins
  integer, intent(in) :: nspin
  !  real/complex flag
  logical, intent(in) :: cplx
  !  error flag, 310 < ierr < 321
  integer, intent(out) :: ierr
  integer, dimension(1) :: ierrvec
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
  real(dp) :: invd
  !  true if this grid point is inside the non-local sphere
  logical :: lnloc

  !  inloc replaces rnloc and ratom
  integer, dimension(:,:,:), allocatable :: inloc

  real(dp) :: fact1,vwr,dvwr,x,y,z,xa,ya,za
  real(dp), dimension(clust%type_num) :: rc, ainv, binv, rrmin
  !  number of replicas of the periodic cell to be used in the
  !  construction of non-local spheres around each atom
  !  nreplica = 0 if .not. pbc%is_on
  integer nreplica(3)
  !  temporary counter of number of non-local points around each atom
  integer nloc
  !  loop counters 
  integer i,iat,ity,ja,jja,lm,lp,mg,icellx,icelly,icellz,im
  !  transformation of non-local arrays
  integer ioffset,itot,icur,icount,irow,itran
  !  local aliases of variables contained in structures
  integer maxnloc,natom,nnodes,comm,mpinfo
  !  counter for the number of non-local projectors
  integer ilm
  !  flag for the counting of non-local projectors
  logical ldef

  real(dp) :: rr_tmp
  real(dp) :: uvec(3), xvec(3), xholdvec(3)
  real(dp) :: tmpmat(3,3), tmpnorm

  !---------------------------------------------------------------

  ierr = 0

  nnodes = parallel%procs_num
  comm = parallel%comm

#ifdef MPI
  call MPI_BCAST(p_pot%uu, p_pot%type_num*4, &
       MPI_DOUBLE_PRECISION,parallel%masterid,comm,mpinfo)
  call MPI_BCAST(p_pot%jj, p_pot%type_num*4, &
       MPI_DOUBLE_PRECISION,parallel%masterid,comm,mpinfo)
#endif
  !
  ! Create an empty non_local_psp module for u_pot
  i = 0
  do ity = 1, p_pot%type_num
     do lp = 1, p_pot%nlocp(ity)
        if ( p_pot%uu(lp,ity) /= zero .or. p_pot%jj(lp,ity) /= zero ) &
             i = i + (2*lp - 1)*clust%natmi(ity)
     enddo
  enddo
  u_pot%maxnloc = 0

  if (i == 0) return

  natom = 0
  do ity = 1, p_pot%type_num
     ldef = .false.
     do lp = 1, p_pot%nlocp(ity)
        if ( p_pot%uu(lp,ity) /= zero .or. p_pot%jj(lp,ity) /= zero ) &
             ldef = .true.
     enddo
     if (ldef) natom = natom + clust%natmi(ity)
  enddo
  
  call create_nonloc_pseudo_pot (i,natom,u_pot)

  if (associated(u_pot%denmat)) deallocate(u_pot%denmat)
  if (associated(u_pot%zdenmat)) deallocate(u_pot%zdenmat)

  ilm = 0
  do lp = 1, 5
     ldef = .false.
     do ity = 1, p_pot%type_num
        if (lp > p_pot%nlocp(ity)) cycle
        if ( p_pot%uu(lp,ity) /= zero .or. p_pot%jj(lp,ity) /= zero ) &
             ldef = .true.
     enddo
     if (ldef) then
        ilm = ilm + 2*lp - 1
     else
        cycle
     endif
!
!   For the moment, the U potential is implemented for d orbitals only.
!
     if (lp == 3) then
        call initialize_p_matrices(u_pot)
     else
        if (parallel%iammaster) then
           write(7,*) 'ERROR: U potential defined for orbital l = ',lp-1
           write(7,*) 'For the moment, U can be defined only for p orbital', &
                ' ( l = 2 ).'
           write(7,*) 'Stop.'
        endif
        ierr = 311
        return
     endif
  enddo

  if (cplx) then
     allocate (u_pot%zdenmat (ilm,ilm,natom,nspin))
     u_pot%zdenmat(:,:,:,:) = zzero
  else
     allocate (u_pot%denmat (ilm,ilm,natom,nspin))
     u_pot%denmat(:,:,:,:) = zero
  endif
  
  kpnum = max(pbc%nkpt,1)
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
  !  Initialize the atom counter, ja, and arrays in u_pot structure.
  ja = 0
  jja = 0
  u_pot%nlatom(:) = 0
  u_pot%nlmatm(:) = 0
  u_pot%skbi(:) = zero

  !  for each atom type 
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
     !  Choose core-cutoff radius such that 90% of the wavefunction
     !  distribution is within that radius.
     rc(ity) = zero
     do lp = 1, p_pot%nlocp(ity)
        if ( p_pot%uu(lp,ity) == zero .and. p_pot%jj(lp,ity) == zero ) cycle
        fact1 = p_pot%rs(2,ity) - p_pot%rs(1,ity)
        tmpnorm = ( p_pot%wfspd(1,ity,lp)*p_pot%rs(1,ity) )**2
        do i = 2, p_pot%ns(ity)
           fact1 = half*(p_pot%rs(i+1,ity)-p_pot%rs(i-1,ity))
           tmpnorm = tmpnorm + &
                fact1*( p_pot%wfspd(i,ity,lp)*p_pot%rs(i,ity) )**2
           if (tmpnorm > 0.95d0) exit
        enddo
        rr_tmp = p_pot%rs(i,ity)
        if ( rc(ity) < rr_tmp ) rc(ity) = rr_tmp
     enddo

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
        !if no on-site Coulomb interaction for an atom, skip the loop
        if ( rc(ity) == zero ) cycle
        jja = jja + 1
        !  distribute atoms round-robin around processors
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

        !  reset counter of non-local points per given atom
        nloc = 0
        !  for each grid point in the full grid...
        do mg = 1, grid%ndim
           !  allow each point to be shifted by the length of the box in each
           !  direction
           do icellx = -nreplica(1),nreplica(1)
              do icelly = -nreplica(2),nreplica(2)
                 do icellz = -nreplica(3),nreplica(3)
                    !  Compute distance between point (or its replica if PBC)
                    !  and atom and increase nloc if it is inside the
                    !  non-local range.

                    lnloc = .false.
                    do idgx = -ndouble + 1, ndouble - 1
                       do idgy = -ndouble + 1, ndouble - 1
                          do idgz = -ndouble + 1, ndouble - 1

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

                             if (rr_tmp <= rc(ity)) lnloc = .true.
                          enddo
                       enddo
                    enddo
                    if (lnloc) nloc = nloc + 1

                 enddo
              enddo
           enddo
        enddo
        !  Update array of total number of non-local points per atom and
        !  maximum size of non-local block associated with any atom.
        u_pot%nlatom(jja) = nloc
     enddo                  ! iat = 1, clust%natmi(ity)
  enddo                     ! ity = 1, clust%type_num
  call pisum(u_pot%nlatom,natom,nnodes,comm)
  call nonloc_pseudo_pot_set_maxnloc (u_pot)
  maxnloc = u_pot%maxnloc
  if (maxnloc > 0) then
     allocate(inloc(4,maxnloc,natom))
     inloc = 0
  endif

  !  Use u_pot%so_indx to store the list of atoms which have U potential.
  allocate(u_pot%so_indx(clust%atom_num))
  u_pot%so_indx = 0

  ja = 0
  jja = 0
  !  for each atom type
  do ity = 1, clust%type_num
     do iat = 1, clust%natmi(ity)
        ja = ja + 1
        !if no on-site Coulomb interaction for an atom, skip the loop
        if ( rc(ity) == zero ) cycle
        jja = jja + 1
        u_pot%so_indx(ja) = jja
        !  distribute atoms round-robin around processors
        if (mod(ja-1,nnodes) /= parallel%iam) cycle
        xa = clust%xatm(ja)
        ya = clust%yatm(ja)
        za = clust%zatm(ja)
        nloc = 0                  
        do mg = 1, grid%ndim
           do icellx = -nreplica(1),nreplica(1)
              do icelly = -nreplica(2),nreplica(2)
                 do icellz = -nreplica(3),nreplica(3)
                    !  Compute distance between point (or its replica if PBC)
                    !  and atom.
                    !  If the point is inside non-local range:
                    !  1. Increase the non-local point counter, nloc, by one.
                    !  Update the four output arrays, indw, ratom,
                    !  with the grid point, its distance from the atom, the
                    !  closest radial grid point, and the fractional distance
                    !  to the next grid point (for extrapolation),
                    !  respectively.
                    !  2. Store the cartesian coordinates of this point in!
                    !  rnloc, without bringing it back to the initial cell
                    !  (if PBC).

                    lnloc = .false.
                    do idgx = -ndouble + 1, ndouble - 1
                       do idgy = -ndouble + 1, ndouble - 1
                          do idgz = -ndouble + 1, ndouble - 1

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
                          enddo
                       enddo
                    enddo
                    if (lnloc) then
                       nloc = nloc + 1
                       u_pot%indw(nloc,jja) = grid%rindex(mg)
                       u_pot%tran(nloc,jja) = grid%rtrans(mg)
                       inloc(:,nloc,jja) = (/ mg, icellx, icelly, icellz /)

                       if(pbc%nkpt > 0) then
                          do kplp = 1, kpnum

                           u_pot%right(nloc,kplp,jja)= &
                   exp( zi * (pbc%kpts(1,kplp)*xholdvec(1)+ &
                              pbc%kpts(2,kplp)*xholdvec(2)+ &
                              pbc%kpts(3,kplp)*xholdvec(3)) )
                           
                           u_pot%left(nloc,kplp,jja)= &
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
        if (nloc /= u_pot%nlatom(jja)) then
           write(9,*) ' ERROR in upot ', &
                nloc,iat,jja,ity,u_pot%nlatom(jja)
           ierr = 312
        endif
     enddo                  ! iat = 1, clust%natmi(ity)
  enddo                     ! ity = 1, clust%type_num
  ierrvec = ierr
  call pisum(ierrvec,1,nnodes,comm)
  ierr = ierrvec(1)
  if (ierr /= 0) then
     ierr = 312
     if (parallel%iammaster) write(9,*) ' ERROR in upot '
     return
  endif
  !
  !  Report non-local block sizes to file.
  !
  if (parallel%iammaster .and. u_pot%maxnloc > 0) then
     write(7,*)
     write(7,*) 'On-site Coulomb interaction (LDA+U) messages:'
     write(7,*) '---------------------------------------------'
     write(7,*)
     write(7,*) ' Cut-off radii [bohr] '
     do ity = 1, clust%type_num
        if ( maxval(abs(p_pot%uu(:,ity))) == zero .and. &
             maxval(abs(p_pot%jj(:,ity))) == zero ) cycle
        write(7,'(a2,f12.4)') clust%name(ity), rc(ity)
     enddo
     write(7,*)
     write(7,*) ' Potential constant (eV) '
     do ity = 1, clust%type_num
        do lp = 1, p_pot%nlocp(ity)
           if ( p_pot%uu(lp,ity) == zero .and. p_pot%jj(lp,ity) == zero ) cycle
           write(7,'(a2,a,f12.4,a,f12.4,a,i2,a)') clust%name(ity), &
                ' U = ', p_pot%uu(lp,ity), ' Ry       J = ',  &
                p_pot%jj(lp,ity), ' Ry  ( l = ',lp-1, ' )'
        enddo
     enddo
     write(7,*)
     write(7,*) 'Max # of non-local points for one atom = ',u_pot%maxnloc
     write(7,*)
     write(7,*) 'Sizes of all non-local blocks '
     write(7,'(12(1x,i5))') (u_pot%nlatom(i),i=1,u_pot%atom_num)
     write(7,*)
  endif
  !
  !  Calculate the non-local potential vector for each angular component.
  !  The result is stored in anloc. The non-local pseudopotential to the
  !  hamiltonian is subsequently calculated as the outer product of these
  !  vectors, summed over all angular components l,m.
  !
  ja = 0
  jja = 0
  !  Initialize the counter of non-local projectors.
  ilm = 0
  !
  !  for each atom type, for each atom within each type
  do ity = 1, clust%type_num
     do iat = 1, clust%natmi(ity)
        !  again update temporary storage of atom # and position
        ja = ja + 1
        !if no on-site Coulomb interaction for an atom, skip the loop
        if (u_pot%so_indx(ja) == 0) cycle
        jja = jja + 1

        !  distribute atoms round-robin around processors
        if (mod(ja-1,nnodes) /= parallel%iam) then
           do lp = 1, p_pot%nlocp(ity)
              if (p_pot%uu(lp,ity) /= zero .or. p_pot%jj(lp,ity) /= zero) &
                   ilm = ilm + 2*lp - 1
           enddo
           cycle
        endif

        nloc = u_pot%nlatom(jja)
        !  initialize angular momentum (l,m) count
        lm = 0
        !
        !  loop over each angular momentum component, l (and within
        !  the loop calculate for different m's too).
        do lp = 1, p_pot%nlocp(ity)
           !  skip if there is no U or J potential
           if ( p_pot%uu(lp,ity) == zero .and. p_pot%jj(lp,ity) == zero ) cycle

           !  Calculate the multiplicand - (non-local pseudopotential)*
           !  (radial pseudowavefunction)*(angular wave function), put the
           !  results in anloc. That is, calculate explicitly |Vl-Vloc|phi_lm>.

           do i = 1, nloc
              an_tmp(:) = zero
              vylmd_tmp(:,:) = zero

              mg = inloc(1,i,jja)
              icellx = inloc(2,i,jja)
              icelly = inloc(3,i,jja)
              icellz = inloc(4,i,jja)
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
                       !  in order to avoid dividing by zero:
                       if (rnorm < rrmin(ity))  then
                          rnorm = rrmin(ity)
                       endif
                       if (rnorm > rc(ity)) cycle

                       npt = idint(binv(ity) * log(one+ainv(ity)*rnorm))+1
                       fact1 = (rnorm - p_pot%rs(npt,ity)) /  &
                            (p_pot%rs(npt+1,ity) - p_pot%rs(npt,ity))

                       vwr  = p_pot%wfspd(npt, ity, lp) + fact1* &
                            (p_pot%wfspd(npt+1,ity,lp)-p_pot%wfspd(npt,ity,lp))
                       dvwr  = p_pot%dwfspd(npt, ity, lp) + fact1* &
                            (p_pot%dwfspd(npt+1,ity,lp)-p_pot%dwfspd(npt,ity,lp))

                       vwr = vwr * fdouble(idgx)*fdouble(idgy)*fdouble(idgz)
                       dvwr = dvwr * fdouble(idgx)*fdouble(idgy)*fdouble(idgz)

                       call ylm_real(lp,rr,rnorm,ylm,ylmd)
                       rr = rr / rnorm
                       do im = 1, 2*lp - 1
                          an_tmp(im) = an_tmp(im) + ylm(im)*vwr
                          vylmd_tmp(:,im) = vylmd_tmp(:,im) + &
                               ylmd(:,im)*vwr + dvwr*ylm(im)*rr
                       enddo

                    enddo
                 enddo
              enddo

              mg = inloc(1,i,jja)
              icellx = inloc(2,i,jja)
              icelly = inloc(3,i,jja)
              icellz = inloc(4,i,jja)
              uvec(1)=(grid%shift(1) + grid%fx(mg))*  &
                   grid%step(1) + real(icellx,dp)*pbc%box_size(1)
              uvec(2)=(grid%shift(2) + grid%fy(mg))*  &
                   grid%step(2) + real(icelly,dp)*pbc%box_size(2)
              uvec(3)=(grid%shift(3) + grid%fz(mg))*  &
                   grid%step(3) + real(icellz,dp)*pbc%box_size(3)

              call matvec3('N',pbc%avec_norm,uvec,xvec)
              rr(1) = xvec(1) - clust%xatm(ja)
              rr(2) = xvec(2) - clust%yatm(ja)
              rr(3) = xvec(3) - clust%zatm(ja)
              rnorm = sqrt(dot_product(rr,rr))

              do im = 1, 2*lp - 1
                 u_pot%anloc(i, ilm+im) = an_tmp(im)
                 u_pot%vylmd(:, i, ilm+im) = vylmd_tmp(:,im)
              enddo

           enddo

           !  store p_pot%uu/jj(ity) in u_pot%uu/jj(jja)
           do im = 1, 2*lp - 1
              u_pot%uu(ilm+im) = p_pot%uu(lp,ity)
              u_pot%jj(ilm+im) = p_pot%jj(lp,ity)
           enddo

           lm = lm + 2*lp - 1
           ilm = ilm + 2*lp - 1

        enddo               ! lp = 1, p_pot%nlocp(ity)
        !  update array with number of lm components for each atom
        u_pot%nlmatm(jja) = lm
     enddo                  ! iat = 1, clust%natmi(ity)
  enddo                     ! ity = 1, clust%type_num

  deallocate(fdouble)
  !
  ! If on-site U potential is not present, bail out. There is nothing else
  ! to do.
  !
  if (u_pot%maxnloc == 0) return

  deallocate(inloc)
  !
  ! Must rescale projectors because anloc and Kohn-Sham eigenvectors have
  ! different normalizations (nloc_p_pot%anloc does not need to be rescaled
  ! because of nloc_p_pot%skbi).
  !
  rr_tmp = sqrt( grid%hcub )
  u_pot%anloc = u_pot%anloc * rr_tmp
  u_pot%vylmd = u_pot%vylmd * rr_tmp

  !
  !  Share information. Follow the round-robin distribution of atoms.
  !
  call pisum(u_pot%indw,maxnloc*u_pot%atom_num,nnodes,comm)
  call pisum(u_pot%tran,maxnloc*u_pot%atom_num,nnodes,comm)
  call pisum(u_pot%nlmatm,u_pot%atom_num,nnodes,comm)
  call psum(u_pot%anloc,maxnloc*u_pot%nlm,nnodes,comm)
  call psum(u_pot%vylmd,3*maxnloc*u_pot%nlm,nnodes,comm)
  call psum(u_pot%uu,u_pot%nlm,nnodes,comm)
  call psum(u_pot%jj,u_pot%nlm,nnodes,comm)
  if (pbc%nkpt > 0) then
     call zpsum(u_pot%right,maxnloc*pbc%nkpt*u_pot%atom_num, &
          nnodes,comm)
     call zpsum(u_pot%left,maxnloc*pbc%nkpt*u_pot%atom_num, &
          nnodes,comm)
  endif
  !
  !  -------- Transformation of the non local arrays -------------
  !  Transform the global "non-local" arrays into arrays containing
  !  information about the local rows of the processor.i.e.,
  !  For all atoms on my processor:
  !  a) all rows of Anloc that are not my rows are removed.
  !  b) my rows are shifted in adjacent positions.
  !  c) indw is changed to include only the above local information
  !  The offset is subtracted to yield the local row locations
  !  d) nlatom contains the number of my rows in Anloc
  !  e) nlmatm is the same since it does not depend on rows:
  !
  write(9,*) 'Transform atoms', u_pot%atom_num
  ioffset = parallel%irows(parallel%group_iam) - 1

  icount = 0
  itot = 0
  ilm = 0
  do jja = 1, u_pot%atom_num
     !  ..for all atoms
     icur = 0
     nloc = u_pot%nlatom(jja)
     if (nloc == 0) cycle
     do i = 1, nloc
        !  ..for each row in Anloc
        irow = u_pot%indw(i,jja)
        itran = u_pot%tran(i,jja)
        if ( parallel%irows(parallel%group_iam) <= irow .and. &
             irow < parallel%irows(parallel%group_iam+1) ) then
           !  ..if irow on my processor shift it up
           icur = icur + 1
           do lm = 1, u_pot%nlmatm(jja)
              u_pot%anloc(icur,ilm+lm) = u_pot%anloc(i,ilm+lm)
              u_pot%vylmd(:,icur,ilm+lm) = u_pot%vylmd(:,i,ilm+lm)
           enddo
           if (pbc%nkpt > 0) then
              u_pot%right(icur,:,jja) = u_pot%right(i,:,jja)
              u_pot%left(icur,:,jja) = u_pot%left(i,:,jja)
           endif
           u_pot%indw(icur,jja) = irow - ioffset
           u_pot%tran(icur,jja) = itran
        endif
     enddo
     ilm = ilm + u_pot%nlmatm(jja)
     !  ..new number of rows in Anloc for this atom.
     u_pot%nlatom(jja) = icur
     icount = icount + icur
     itot   = itot + nloc
  enddo
  write(9,*) 'My number of non-local rows: ',icount,' of total ',itot

end subroutine upot
!===============================================================
!
!  Performs the operation Sum_{m,n} A_{m,n} * B _{m,n}, with
!  m,n = 1,...,5 and complex matrices A and B.
!
!---------------------------------------------------------------
function contraction(aa,bb)

  use constants
  implicit none
  !
  !  Input/Output variables:
  !
  complex(dpc), dimension(5,5), intent(in) :: aa, bb
  complex(dpc) :: contraction
  !
  !  Work variables:
  !
  integer :: ii, jj

  !---------------------------------------------------------------

  contraction = zzero
  do ii = 1, 5
     do jj = 1, 5
        contraction = contraction + aa(ii,jj) * bb(ii,jj)
     enddo
  enddo

end function contraction
!===============================================================
!
! Initializes matrices p3, p4 for the U potential. Definitions
! are consistent with: Liechtenstein et al., 64, R5467 (1995).
!
!---------------------------------------------------------------
subroutine initialize_p_matrices(u_pot)

  use constants
  use non_local_psp_module
  implicit none
  !
  !  Input/Output variables:
  !
  !  non local pseudo_potential related data
  type (nonloc_pseudo_potential), intent(inout) :: u_pot

  !---------------------------------------------------------------

  if (.not. associated(u_pot%p3)) allocate(u_pot%p3(5,5,5,5))
  u_pot%p3 = zzero
  if (.not. associated(u_pot%p4)) allocate(u_pot%p4(5,5,5,5))
  u_pot%p4 = zzero

  u_pot%p3(1,1,1,1) = cmplx(   0.500000000000000,   0.000000000000000)
  u_pot%p3(1,1,2,2) = cmplx(   0.141156462585034,   0.087912087912088)
  u_pot%p3(1,1,2,3) = cmplx(   0.172553636839351,  -0.087912087912088)
  u_pot%p3(1,1,3,2) = cmplx(   0.172553636839351,  -0.087912087912088)
  u_pot%p3(1,1,3,3) = cmplx(   0.141156462585034,   0.087912087912088)
  u_pot%p3(1,1,4,2) = cmplx(  -0.186489700532716,   0.000000000000000)
  u_pot%p3(1,1,4,3) = cmplx(  -0.186489700532716,   0.000000000000000)
  u_pot%p3(1,1,4,4) = cmplx(   0.049973835688122,   0.000000000000000)
  u_pot%p3(1,1,5,5) = cmplx(   0.821559392987964,   0.000000000000000)
  u_pot%p3(1,2,1,4) = cmplx(   0.186489700532716,   0.000000000000000)
  u_pot%p3(1,2,2,1) = cmplx(   0.358843537414966,  -0.087912087912088)
  u_pot%p3(1,2,2,3) = cmplx(   0.000000000000000,  -0.107669878803656)
  u_pot%p3(1,2,2,4) = cmplx(  -0.076134101431599,  -0.222737564604946)
  u_pot%p3(1,2,2,5) = cmplx(   0.091182626896913,  -0.087912087912088)
  u_pot%p3(1,2,3,1) = cmplx(  -0.172553636839351,   0.087912087912088)
  u_pot%p3(1,2,3,3) = cmplx(   0.000000000000000,  -0.107669878803656)
  u_pot%p3(1,2,3,4) = cmplx(  -0.228402304294797,  -0.222737564604946)
  u_pot%p3(1,2,3,5) = cmplx(  -0.084641548927263,   0.087912087912088)
  u_pot%p3(1,2,4,1) = cmplx(   0.186489700532716,   0.000000000000000)
  u_pot%p3(1,2,4,3) = cmplx(  -0.011329479379702,   0.000000000000000)
  u_pot%p3(1,2,4,5) = cmplx(   0.186489700532716,   0.000000000000000)
  u_pot%p3(1,2,5,4) = cmplx(  -0.177239219355498,   0.000000000000000)
  u_pot%p3(1,3,1,4) = cmplx(   0.186489700532716,   0.000000000000000)
  u_pot%p3(1,3,2,1) = cmplx(  -0.172553636839351,   0.087912087912088)
  u_pot%p3(1,3,2,2) = cmplx(   0.000000000000000,   0.107669878803656)
  u_pot%p3(1,3,2,4) = cmplx(   0.228402304294797,   0.222737564604946)
  u_pot%p3(1,3,2,5) = cmplx(  -0.084641548927263,   0.087912087912088)
  u_pot%p3(1,3,3,1) = cmplx(   0.358843537414966,  -0.087912087912088)
  u_pot%p3(1,3,3,2) = cmplx(   0.000000000000000,   0.107669878803656)
  u_pot%p3(1,3,3,4) = cmplx(   0.076134101431599,   0.222737564604946)
  u_pot%p3(1,3,3,5) = cmplx(   0.091182626896913,  -0.087912087912088)
  u_pot%p3(1,3,4,1) = cmplx(   0.186489700532716,   0.000000000000000)
  u_pot%p3(1,3,4,2) = cmplx(   0.011329479379702,   0.000000000000000)
  u_pot%p3(1,3,4,5) = cmplx(   0.186489700532716,   0.000000000000000)
  u_pot%p3(1,3,5,4) = cmplx(  -0.177239219355498,   0.000000000000000)
  u_pot%p3(1,4,1,2) = cmplx(  -0.186489700532716,   0.000000000000000)
  u_pot%p3(1,4,1,3) = cmplx(  -0.186489700532716,   0.000000000000000)
  u_pot%p3(1,4,2,2) = cmplx(   0.076134101431599,   0.222737564604946)
  u_pot%p3(1,4,2,3) = cmplx(  -0.228402304294797,  -0.222737564604946)
  u_pot%p3(1,4,3,2) = cmplx(   0.228402304294797,   0.222737564604946)
  u_pot%p3(1,4,3,3) = cmplx(  -0.076134101431599,  -0.222737564604946)
  u_pot%p3(1,4,4,1) = cmplx(   0.450026164311878,   0.000000000000000)
  u_pot%p3(1,4,4,5) = cmplx(   0.257195185766614,   0.000000000000000)
  u_pot%p3(1,4,5,2) = cmplx(   0.177239219355498,   0.000000000000000)
  u_pot%p3(1,4,5,3) = cmplx(   0.177239219355498,   0.000000000000000)
  u_pot%p3(1,5,2,2) = cmplx(  -0.091182626896913,   0.087912087912088)
  u_pot%p3(1,5,2,3) = cmplx(   0.084641548927263,  -0.087912087912088)
  u_pot%p3(1,5,3,2) = cmplx(   0.084641548927263,  -0.087912087912088)
  u_pot%p3(1,5,3,3) = cmplx(  -0.091182626896913,   0.087912087912088)
  u_pot%p3(1,5,4,2) = cmplx(  -0.186489700532716,   0.000000000000000)
  u_pot%p3(1,5,4,3) = cmplx(  -0.186489700532716,   0.000000000000000)
  u_pot%p3(1,5,4,4) = cmplx(  -0.257195185766614,   0.000000000000000)
  u_pot%p3(1,5,5,1) = cmplx(  -0.321559392987964,   0.000000000000000)
  u_pot%p3(2,1,1,2) = cmplx(   0.490711669283098,  -0.087912087912088)
  u_pot%p3(2,1,1,3) = cmplx(  -0.040685504971219,   0.087912087912088)
  u_pot%p3(2,1,2,2) = cmplx(   0.000000000000000,  -0.053834939401828)
  u_pot%p3(2,1,2,3) = cmplx(   0.000000000000000,   0.053834939401828)
  u_pot%p3(2,1,2,4) = cmplx(   0.076134101431599,   0.000000000000000)
  u_pot%p3(2,1,2,5) = cmplx(   0.000000000000000,   0.175824175824176)
  u_pot%p3(2,1,3,2) = cmplx(   0.000000000000000,  -0.053834939401828)
  u_pot%p3(2,1,3,3) = cmplx(   0.000000000000000,   0.053834939401828)
  u_pot%p3(2,1,3,4) = cmplx(   0.217072824915095,  -0.304536405726396)
  u_pot%p3(2,1,3,5) = cmplx(   0.000000000000000,  -0.175824175824176)
  u_pot%p3(2,1,4,2) = cmplx(   0.000000000000000,   0.146603463173347)
  u_pot%p3(2,1,4,3) = cmplx(   0.000000000000000,  -0.005664739689851)
  u_pot%p3(2,1,4,4) = cmplx(   0.107669878803656,   0.000000000000000)
  u_pot%p3(2,1,4,5) = cmplx(  -0.363728919888213,   0.000000000000000)
  u_pot%p3(2,1,5,2) = cmplx(  -0.216509680795395,   0.087912087912088)
  u_pot%p3(2,1,5,3) = cmplx(  -0.040685504971219,  -0.087912087912088)
  u_pot%p3(2,2,1,1) = cmplx(   0.009288330716902,   0.087912087912088)
  u_pot%p3(2,2,1,3) = cmplx(   0.000000000000000,   0.107669878803656)
  u_pot%p3(2,2,1,4) = cmplx(   0.076134101431599,  -0.005664739689851)
  u_pot%p3(2,2,1,5) = cmplx(   0.216509680795395,   0.087912087912088)
  u_pot%p3(2,2,2,1) = cmplx(   0.000000000000000,   0.053834939401828)
  u_pot%p3(2,2,2,2) = cmplx(   0.500000000000000,   0.000000000000000)
  u_pot%p3(2,2,2,3) = cmplx(  -0.043956043956044,   0.000000000000000)
  u_pot%p3(2,2,2,4) = cmplx(   0.000000000000000,  -0.031081616755453)
  u_pot%p3(2,2,2,5) = cmplx(   0.000000000000000,  -0.053834939401828)
  u_pot%p3(2,2,3,1) = cmplx(   0.000000000000000,   0.053834939401828)
  u_pot%p3(2,2,3,3) = cmplx(   0.608320251177394,   0.000000000000000)
  u_pot%p3(2,2,3,4) = cmplx(   0.000000000000000,   0.093244850266358)
  u_pot%p3(2,2,3,5) = cmplx(   0.000000000000000,  -0.053834939401828)
  u_pot%p3(2,2,4,1) = cmplx(   0.000000000000000,  -0.146603463173347)
  u_pot%p3(2,2,4,4) = cmplx(   0.614861329147043,   0.000000000000000)
  u_pot%p3(2,2,4,5) = cmplx(   0.000000000000000,  -0.146603463173347)
  u_pot%p3(2,2,5,1) = cmplx(   0.216509680795395,  -0.087912087912088)
  u_pot%p3(2,2,5,3) = cmplx(   0.000000000000000,  -0.107669878803656)
  u_pot%p3(2,2,5,4) = cmplx(  -0.076134101431599,  -0.005664739689851)
  u_pot%p3(2,2,5,5) = cmplx(   0.009288330716902,  -0.087912087912088)
  u_pot%p3(2,3,1,1) = cmplx(   0.040685504971219,  -0.087912087912088)
  u_pot%p3(2,3,1,2) = cmplx(   0.000000000000000,  -0.107669878803656)
  u_pot%p3(2,3,1,4) = cmplx(  -0.228402304294797,   0.005664739689851)
  u_pot%p3(2,3,1,5) = cmplx(   0.040685504971219,  -0.087912087912088)
  u_pot%p3(2,3,2,1) = cmplx(   0.000000000000000,  -0.053834939401828)
  u_pot%p3(2,3,2,2) = cmplx(   0.043956043956044,   0.000000000000000)
  u_pot%p3(2,3,2,4) = cmplx(   0.000000000000000,  -0.031081616755453)
  u_pot%p3(2,3,2,5) = cmplx(   0.000000000000000,   0.053834939401828)
  u_pot%p3(2,3,3,1) = cmplx(   0.000000000000000,  -0.053834939401828)
  u_pot%p3(2,3,3,2) = cmplx(  -0.108320251177394,   0.000000000000000)
  u_pot%p3(2,3,3,4) = cmplx(   0.000000000000000,  -0.031081616755453)
  u_pot%p3(2,3,3,5) = cmplx(   0.000000000000000,   0.053834939401828)
  u_pot%p3(2,3,4,1) = cmplx(   0.000000000000000,   0.005664739689851)
  u_pot%p3(2,3,4,4) = cmplx(  -0.037414965986395,   0.000000000000000)
  u_pot%p3(2,3,4,5) = cmplx(   0.000000000000000,   0.005664739689851)
  u_pot%p3(2,3,5,1) = cmplx(   0.040685504971219,   0.087912087912088)
  u_pot%p3(2,3,5,2) = cmplx(   0.000000000000000,   0.107669878803656)
  u_pot%p3(2,3,5,4) = cmplx(   0.228402304294797,   0.005664739689851)
  u_pot%p3(2,3,5,5) = cmplx(   0.040685504971219,   0.087912087912088)
  u_pot%p3(2,4,1,2) = cmplx(  -0.076134101431599,   0.005664739689851)
  u_pot%p3(2,4,1,3) = cmplx(   0.228402304294797,  -0.005664739689851)
  u_pot%p3(2,4,2,1) = cmplx(  -0.076134101431599,   0.000000000000000)
  u_pot%p3(2,4,2,2) = cmplx(   0.000000000000000,   0.031081616755453)
  u_pot%p3(2,4,2,3) = cmplx(   0.000000000000000,   0.031081616755453)
  u_pot%p3(2,4,2,5) = cmplx(   0.076134101431599,   0.000000000000000)
  u_pot%p3(2,4,3,1) = cmplx(  -0.217072824915095,   0.304536405726396)
  u_pot%p3(2,4,3,2) = cmplx(   0.000000000000000,  -0.093244850266358)
  u_pot%p3(2,4,3,3) = cmplx(   0.000000000000000,   0.031081616755453)
  u_pot%p3(2,4,3,5) = cmplx(   0.217072824915095,   0.304536405726396)
  u_pot%p3(2,4,4,1) = cmplx(  -0.107669878803656,   0.000000000000000)
  u_pot%p3(2,4,4,2) = cmplx(  -0.114861329147043,   0.000000000000000)
  u_pot%p3(2,4,4,3) = cmplx(   0.037414965986395,   0.000000000000000)
  u_pot%p3(2,4,4,5) = cmplx(  -0.107669878803656,   0.000000000000000)
  u_pot%p3(2,4,5,2) = cmplx(   0.076134101431599,   0.005664739689851)
  u_pot%p3(2,4,5,3) = cmplx(  -0.228402304294797,  -0.005664739689851)
  u_pot%p3(2,5,1,2) = cmplx(  -0.216509680795395,  -0.087912087912088)
  u_pot%p3(2,5,1,3) = cmplx(  -0.040685504971219,   0.087912087912088)
  u_pot%p3(2,5,2,1) = cmplx(   0.000000000000000,  -0.175824175824176)
  u_pot%p3(2,5,2,2) = cmplx(   0.000000000000000,   0.053834939401828)
  u_pot%p3(2,5,2,3) = cmplx(   0.000000000000000,  -0.053834939401828)
  u_pot%p3(2,5,2,4) = cmplx(  -0.076134101431599,   0.000000000000000)
  u_pot%p3(2,5,3,1) = cmplx(   0.000000000000000,   0.175824175824176)
  u_pot%p3(2,5,3,2) = cmplx(   0.000000000000000,   0.053834939401828)
  u_pot%p3(2,5,3,3) = cmplx(   0.000000000000000,  -0.053834939401828)
  u_pot%p3(2,5,3,4) = cmplx(  -0.217072824915095,  -0.304536405726396)
  u_pot%p3(2,5,4,1) = cmplx(   0.363728919888213,   0.000000000000000)
  u_pot%p3(2,5,4,2) = cmplx(   0.000000000000000,   0.146603463173347)
  u_pot%p3(2,5,4,3) = cmplx(   0.000000000000000,  -0.005664739689851)
  u_pot%p3(2,5,4,4) = cmplx(   0.107669878803656,   0.000000000000000)
  u_pot%p3(2,5,5,2) = cmplx(   0.490711669283098,   0.087912087912088)
  u_pot%p3(2,5,5,3) = cmplx(  -0.040685504971219,  -0.087912087912088)
  u_pot%p3(3,1,1,2) = cmplx(  -0.040685504971219,   0.087912087912088)
  u_pot%p3(3,1,1,3) = cmplx(   0.490711669283098,  -0.087912087912088)
  u_pot%p3(3,1,2,2) = cmplx(   0.000000000000000,  -0.053834939401828)
  u_pot%p3(3,1,2,3) = cmplx(   0.000000000000000,   0.053834939401828)
  u_pot%p3(3,1,2,4) = cmplx(  -0.217072824915095,   0.304536405726396)
  u_pot%p3(3,1,2,5) = cmplx(   0.000000000000000,  -0.175824175824176)
  u_pot%p3(3,1,3,2) = cmplx(   0.000000000000000,  -0.053834939401828)
  u_pot%p3(3,1,3,3) = cmplx(   0.000000000000000,   0.053834939401828)
  u_pot%p3(3,1,3,4) = cmplx(  -0.076134101431599,   0.000000000000000)
  u_pot%p3(3,1,3,5) = cmplx(   0.000000000000000,   0.175824175824176)
  u_pot%p3(3,1,4,2) = cmplx(   0.000000000000000,   0.005664739689851)
  u_pot%p3(3,1,4,3) = cmplx(   0.000000000000000,  -0.146603463173347)
  u_pot%p3(3,1,4,4) = cmplx(  -0.107669878803656,   0.000000000000000)
  u_pot%p3(3,1,4,5) = cmplx(  -0.363728919888213,   0.000000000000000)
  u_pot%p3(3,1,5,2) = cmplx(  -0.040685504971219,  -0.087912087912088)
  u_pot%p3(3,1,5,3) = cmplx(  -0.216509680795395,   0.087912087912088)
  u_pot%p3(3,2,1,1) = cmplx(   0.040685504971219,  -0.087912087912088)
  u_pot%p3(3,2,1,3) = cmplx(   0.000000000000000,   0.107669878803656)
  u_pot%p3(3,2,1,4) = cmplx(   0.228402304294797,  -0.005664739689851)
  u_pot%p3(3,2,1,5) = cmplx(   0.040685504971219,  -0.087912087912088)
  u_pot%p3(3,2,2,1) = cmplx(   0.000000000000000,   0.053834939401828)
  u_pot%p3(3,2,2,3) = cmplx(  -0.108320251177394,   0.000000000000000)
  u_pot%p3(3,2,2,4) = cmplx(   0.000000000000000,  -0.031081616755453)
  u_pot%p3(3,2,2,5) = cmplx(   0.000000000000000,  -0.053834939401828)
  u_pot%p3(3,2,3,1) = cmplx(   0.000000000000000,   0.053834939401828)
  u_pot%p3(3,2,3,3) = cmplx(   0.043956043956044,   0.000000000000000)
  u_pot%p3(3,2,3,4) = cmplx(   0.000000000000000,  -0.031081616755453)
  u_pot%p3(3,2,3,5) = cmplx(   0.000000000000000,  -0.053834939401828)
  u_pot%p3(3,2,4,1) = cmplx(   0.000000000000000,  -0.005664739689851)
  u_pot%p3(3,2,4,4) = cmplx(  -0.037414965986395,   0.000000000000000)
  u_pot%p3(3,2,4,5) = cmplx(   0.000000000000000,  -0.005664739689851)
  u_pot%p3(3,2,5,1) = cmplx(   0.040685504971219,   0.087912087912088)
  u_pot%p3(3,2,5,3) = cmplx(   0.000000000000000,  -0.107669878803656)
  u_pot%p3(3,2,5,4) = cmplx(  -0.228402304294797,  -0.005664739689851)
  u_pot%p3(3,2,5,5) = cmplx(   0.040685504971219,   0.087912087912088)
  u_pot%p3(3,3,1,1) = cmplx(   0.009288330716902,   0.087912087912088)
  u_pot%p3(3,3,1,2) = cmplx(   0.000000000000000,  -0.107669878803656)
  u_pot%p3(3,3,1,4) = cmplx(  -0.076134101431599,   0.005664739689851)
  u_pot%p3(3,3,1,5) = cmplx(   0.216509680795395,   0.087912087912088)
  u_pot%p3(3,3,2,1) = cmplx(   0.000000000000000,  -0.053834939401828)
  u_pot%p3(3,3,2,2) = cmplx(   0.608320251177394,   0.000000000000000)
  u_pot%p3(3,3,2,4) = cmplx(   0.000000000000000,   0.093244850266358)
  u_pot%p3(3,3,2,5) = cmplx(   0.000000000000000,   0.053834939401828)
  u_pot%p3(3,3,3,1) = cmplx(   0.000000000000000,  -0.053834939401828)
  u_pot%p3(3,3,3,2) = cmplx(  -0.043956043956044,   0.000000000000000)
  u_pot%p3(3,3,3,3) = cmplx(   0.500000000000000,   0.000000000000000)
  u_pot%p3(3,3,3,4) = cmplx(   0.000000000000000,  -0.031081616755453)
  u_pot%p3(3,3,3,5) = cmplx(   0.000000000000000,   0.053834939401828)
  u_pot%p3(3,3,4,1) = cmplx(   0.000000000000000,   0.146603463173347)
  u_pot%p3(3,3,4,4) = cmplx(   0.614861329147043,   0.000000000000000)
  u_pot%p3(3,3,4,5) = cmplx(   0.000000000000000,   0.146603463173347)
  u_pot%p3(3,3,5,1) = cmplx(   0.216509680795395,  -0.087912087912088)
  u_pot%p3(3,3,5,2) = cmplx(   0.000000000000000,   0.107669878803656)
  u_pot%p3(3,3,5,4) = cmplx(   0.076134101431599,   0.005664739689851)
  u_pot%p3(3,3,5,5) = cmplx(   0.009288330716902,  -0.087912087912088)
  u_pot%p3(3,4,1,2) = cmplx(  -0.228402304294797,   0.005664739689851)
  u_pot%p3(3,4,1,3) = cmplx(   0.076134101431599,  -0.005664739689851)
  u_pot%p3(3,4,2,1) = cmplx(   0.217072824915095,  -0.304536405726396)
  u_pot%p3(3,4,2,2) = cmplx(   0.000000000000000,   0.031081616755453)
  u_pot%p3(3,4,2,3) = cmplx(   0.000000000000000,  -0.093244850266358)
  u_pot%p3(3,4,2,5) = cmplx(  -0.217072824915095,  -0.304536405726396)
  u_pot%p3(3,4,3,1) = cmplx(   0.076134101431599,   0.000000000000000)
  u_pot%p3(3,4,3,2) = cmplx(   0.000000000000000,   0.031081616755453)
  u_pot%p3(3,4,3,3) = cmplx(   0.000000000000000,   0.031081616755453)
  u_pot%p3(3,4,3,5) = cmplx(  -0.076134101431599,   0.000000000000000)
  u_pot%p3(3,4,4,1) = cmplx(   0.107669878803656,   0.000000000000000)
  u_pot%p3(3,4,4,2) = cmplx(   0.037414965986395,   0.000000000000000)
  u_pot%p3(3,4,4,3) = cmplx(  -0.114861329147043,   0.000000000000000)
  u_pot%p3(3,4,4,5) = cmplx(   0.107669878803656,   0.000000000000000)
  u_pot%p3(3,4,5,2) = cmplx(   0.228402304294797,   0.005664739689851)
  u_pot%p3(3,4,5,3) = cmplx(  -0.076134101431599,  -0.005664739689851)
  u_pot%p3(3,5,1,2) = cmplx(  -0.040685504971219,   0.087912087912088)
  u_pot%p3(3,5,1,3) = cmplx(  -0.216509680795395,  -0.087912087912088)
  u_pot%p3(3,5,2,1) = cmplx(   0.000000000000000,   0.175824175824176)
  u_pot%p3(3,5,2,2) = cmplx(   0.000000000000000,   0.053834939401828)
  u_pot%p3(3,5,2,3) = cmplx(   0.000000000000000,  -0.053834939401828)
  u_pot%p3(3,5,2,4) = cmplx(   0.217072824915095,   0.304536405726396)
  u_pot%p3(3,5,3,1) = cmplx(   0.000000000000000,  -0.175824175824176)
  u_pot%p3(3,5,3,2) = cmplx(   0.000000000000000,   0.053834939401828)
  u_pot%p3(3,5,3,3) = cmplx(   0.000000000000000,  -0.053834939401828)
  u_pot%p3(3,5,3,4) = cmplx(   0.076134101431599,   0.000000000000000)
  u_pot%p3(3,5,4,1) = cmplx(   0.363728919888213,   0.000000000000000)
  u_pot%p3(3,5,4,2) = cmplx(   0.000000000000000,   0.005664739689851)
  u_pot%p3(3,5,4,3) = cmplx(   0.000000000000000,  -0.146603463173347)
  u_pot%p3(3,5,4,4) = cmplx(  -0.107669878803656,   0.000000000000000)
  u_pot%p3(3,5,5,2) = cmplx(  -0.040685504971219,  -0.087912087912088)
  u_pot%p3(3,5,5,3) = cmplx(   0.490711669283098,   0.087912087912088)
  u_pot%p3(4,1,1,2) = cmplx(   0.186489700532716,   0.000000000000000)
  u_pot%p3(4,1,1,3) = cmplx(   0.186489700532716,   0.000000000000000)
  u_pot%p3(4,1,1,4) = cmplx(   0.186289900575615,   0.000000000000000)
  u_pot%p3(4,1,2,2) = cmplx(   0.000000000000000,  -0.222737564604946)
  u_pot%p3(4,1,2,3) = cmplx(   0.000000000000000,  -0.070469361741748)
  u_pot%p3(4,1,2,4) = cmplx(  -0.107669878803656,   0.000000000000000)
  u_pot%p3(4,1,2,5) = cmplx(  -0.009250481177218,   0.000000000000000)
  u_pot%p3(4,1,3,2) = cmplx(   0.000000000000000,   0.070469361741748)
  u_pot%p3(4,1,3,3) = cmplx(   0.000000000000000,   0.222737564604946)
  u_pot%p3(4,1,3,4) = cmplx(   0.107669878803656,   0.000000000000000)
  u_pot%p3(4,1,3,5) = cmplx(  -0.009250481177218,   0.000000000000000)
  u_pot%p3(4,1,4,2) = cmplx(  -0.107669878803656,   0.000000000000000)
  u_pot%p3(4,1,4,3) = cmplx(   0.107669878803656,   0.000000000000000)
  u_pot%p3(4,1,5,2) = cmplx(  -0.186489700532716,   0.000000000000000)
  u_pot%p3(4,1,5,3) = cmplx(  -0.186489700532716,   0.000000000000000)
  u_pot%p3(4,1,5,4) = cmplx(  -0.006541077969649,   0.000000000000000)
  u_pot%p3(4,2,1,1) = cmplx(  -0.186489700532716,   0.000000000000000)
  u_pot%p3(4,2,1,3) = cmplx(  -0.445475129209892,   0.000000000000000)
  u_pot%p3(4,2,1,5) = cmplx(  -0.186489700532716,   0.000000000000000)
  u_pot%p3(4,2,2,1) = cmplx(   0.000000000000000,   0.222737564604946)
  u_pot%p3(4,2,2,4) = cmplx(   0.192830978545264,   0.000000000000000)
  u_pot%p3(4,2,2,5) = cmplx(   0.000000000000000,   0.222737564604946)
  u_pot%p3(4,2,3,1) = cmplx(   0.000000000000000,  -0.070469361741748)
  u_pot%p3(4,2,3,4) = cmplx(   0.257195185766614,   0.000000000000000)
  u_pot%p3(4,2,3,5) = cmplx(   0.000000000000000,  -0.070469361741748)
  u_pot%p3(4,2,4,1) = cmplx(   0.107669878803656,   0.000000000000000)
  u_pot%p3(4,2,4,5) = cmplx(   0.107669878803656,   0.000000000000000)
  u_pot%p3(4,2,5,1) = cmplx(   0.186489700532716,   0.000000000000000)
  u_pot%p3(4,2,5,3) = cmplx(   0.445475129209892,   0.000000000000000)
  u_pot%p3(4,2,5,5) = cmplx(   0.186489700532716,   0.000000000000000)
  u_pot%p3(4,3,1,1) = cmplx(  -0.186489700532716,   0.000000000000000)
  u_pot%p3(4,3,1,2) = cmplx(   0.445475129209892,   0.000000000000000)
  u_pot%p3(4,3,1,5) = cmplx(  -0.186489700532716,   0.000000000000000)
  u_pot%p3(4,3,2,1) = cmplx(   0.000000000000000,   0.070469361741748)
  u_pot%p3(4,3,2,4) = cmplx(   0.257195185766614,   0.000000000000000)
  u_pot%p3(4,3,2,5) = cmplx(   0.000000000000000,   0.070469361741748)
  u_pot%p3(4,3,3,1) = cmplx(   0.000000000000000,  -0.222737564604946)
  u_pot%p3(4,3,3,4) = cmplx(   0.192830978545264,   0.000000000000000)
  u_pot%p3(4,3,3,5) = cmplx(   0.000000000000000,  -0.222737564604946)
  u_pot%p3(4,3,4,1) = cmplx(  -0.107669878803656,   0.000000000000000)
  u_pot%p3(4,3,4,5) = cmplx(  -0.107669878803656,   0.000000000000000)
  u_pot%p3(4,3,5,1) = cmplx(   0.186489700532716,   0.000000000000000)
  u_pot%p3(4,3,5,2) = cmplx(  -0.445475129209892,   0.000000000000000)
  u_pot%p3(4,3,5,5) = cmplx(   0.186489700532716,   0.000000000000000)
  u_pot%p3(4,4,1,1) = cmplx(   0.313710099424385,   0.000000000000000)
  u_pot%p3(4,4,1,5) = cmplx(   0.006541077969649,   0.000000000000000)
  u_pot%p3(4,4,2,1) = cmplx(   0.107669878803656,   0.000000000000000)
  u_pot%p3(4,4,2,2) = cmplx(   0.307169021454736,   0.000000000000000)
  u_pot%p3(4,4,2,3) = cmplx(  -0.257195185766614,   0.000000000000000)
  u_pot%p3(4,4,2,5) = cmplx(   0.107669878803656,   0.000000000000000)
  u_pot%p3(4,4,3,1) = cmplx(  -0.107669878803656,   0.000000000000000)
  u_pot%p3(4,4,3,2) = cmplx(  -0.257195185766614,   0.000000000000000)
  u_pot%p3(4,4,3,3) = cmplx(   0.307169021454736,   0.000000000000000)
  u_pot%p3(4,4,3,5) = cmplx(  -0.107669878803656,   0.000000000000000)
  u_pot%p3(4,4,4,4) = cmplx(   0.500000000000000,   0.000000000000000)
  u_pot%p3(4,4,5,1) = cmplx(   0.006541077969649,   0.000000000000000)
  u_pot%p3(4,4,5,5) = cmplx(   0.313710099424385,   0.000000000000000)
  u_pot%p3(4,5,1,2) = cmplx(   0.186489700532716,   0.000000000000000)
  u_pot%p3(4,5,1,3) = cmplx(   0.186489700532716,   0.000000000000000)
  u_pot%p3(4,5,1,4) = cmplx(  -0.006541077969649,   0.000000000000000)
  u_pot%p3(4,5,2,1) = cmplx(   0.009250481177218,   0.000000000000000)
  u_pot%p3(4,5,2,2) = cmplx(   0.000000000000000,  -0.222737564604946)
  u_pot%p3(4,5,2,3) = cmplx(   0.000000000000000,  -0.070469361741748)
  u_pot%p3(4,5,2,4) = cmplx(  -0.107669878803656,   0.000000000000000)
  u_pot%p3(4,5,3,1) = cmplx(   0.009250481177218,   0.000000000000000)
  u_pot%p3(4,5,3,2) = cmplx(   0.000000000000000,   0.070469361741748)
  u_pot%p3(4,5,3,3) = cmplx(   0.000000000000000,   0.222737564604946)
  u_pot%p3(4,5,3,4) = cmplx(   0.107669878803656,   0.000000000000000)
  u_pot%p3(4,5,4,2) = cmplx(  -0.107669878803656,   0.000000000000000)
  u_pot%p3(4,5,4,3) = cmplx(   0.107669878803656,   0.000000000000000)
  u_pot%p3(4,5,5,2) = cmplx(  -0.186489700532716,   0.000000000000000)
  u_pot%p3(4,5,5,3) = cmplx(  -0.186489700532716,   0.000000000000000)
  u_pot%p3(4,5,5,4) = cmplx(   0.186289900575615,   0.000000000000000)
  u_pot%p3(5,1,1,5) = cmplx(  -0.321559392987964,   0.000000000000000)
  u_pot%p3(5,1,2,2) = cmplx(  -0.091182626896913,  -0.087912087912088)
  u_pot%p3(5,1,2,3) = cmplx(   0.084641548927263,   0.087912087912088)
  u_pot%p3(5,1,3,2) = cmplx(   0.084641548927263,   0.087912087912088)
  u_pot%p3(5,1,3,3) = cmplx(  -0.091182626896913,  -0.087912087912088)
  u_pot%p3(5,1,4,2) = cmplx(   0.186489700532716,   0.000000000000000)
  u_pot%p3(5,1,4,3) = cmplx(   0.186489700532716,   0.000000000000000)
  u_pot%p3(5,1,4,4) = cmplx(  -0.257195185766614,   0.000000000000000)
  u_pot%p3(5,2,1,4) = cmplx(   0.177239219355498,   0.000000000000000)
  u_pot%p3(5,2,2,1) = cmplx(   0.091182626896913,   0.087912087912088)
  u_pot%p3(5,2,2,3) = cmplx(   0.000000000000000,   0.107669878803656)
  u_pot%p3(5,2,2,4) = cmplx(   0.076134101431599,  -0.222737564604946)
  u_pot%p3(5,2,2,5) = cmplx(   0.358843537414966,   0.087912087912088)
  u_pot%p3(5,2,3,1) = cmplx(  -0.084641548927263,  -0.087912087912088)
  u_pot%p3(5,2,3,3) = cmplx(   0.000000000000000,   0.107669878803656)
  u_pot%p3(5,2,3,4) = cmplx(   0.228402304294797,  -0.222737564604946)
  u_pot%p3(5,2,3,5) = cmplx(  -0.172553636839351,  -0.087912087912088)
  u_pot%p3(5,2,4,1) = cmplx(  -0.186489700532716,   0.000000000000000)
  u_pot%p3(5,2,4,3) = cmplx(   0.011329479379702,   0.000000000000000)
  u_pot%p3(5,2,4,5) = cmplx(  -0.186489700532716,   0.000000000000000)
  u_pot%p3(5,2,5,4) = cmplx(  -0.186489700532716,   0.000000000000000)
  u_pot%p3(5,3,1,4) = cmplx(   0.177239219355498,   0.000000000000000)
  u_pot%p3(5,3,2,1) = cmplx(  -0.084641548927263,  -0.087912087912088)
  u_pot%p3(5,3,2,2) = cmplx(   0.000000000000000,  -0.107669878803656)
  u_pot%p3(5,3,2,4) = cmplx(  -0.228402304294797,   0.222737564604946)
  u_pot%p3(5,3,2,5) = cmplx(  -0.172553636839351,  -0.087912087912088)
  u_pot%p3(5,3,3,1) = cmplx(   0.091182626896913,   0.087912087912088)
  u_pot%p3(5,3,3,2) = cmplx(   0.000000000000000,  -0.107669878803656)
  u_pot%p3(5,3,3,4) = cmplx(  -0.076134101431599,   0.222737564604946)
  u_pot%p3(5,3,3,5) = cmplx(   0.358843537414966,   0.087912087912088)
  u_pot%p3(5,3,4,1) = cmplx(  -0.186489700532716,   0.000000000000000)
  u_pot%p3(5,3,4,2) = cmplx(  -0.011329479379702,   0.000000000000000)
  u_pot%p3(5,3,4,5) = cmplx(  -0.186489700532716,   0.000000000000000)
  u_pot%p3(5,3,5,4) = cmplx(  -0.186489700532716,   0.000000000000000)
  u_pot%p3(5,4,1,2) = cmplx(  -0.177239219355498,   0.000000000000000)
  u_pot%p3(5,4,1,3) = cmplx(  -0.177239219355498,   0.000000000000000)
  u_pot%p3(5,4,2,2) = cmplx(  -0.076134101431599,   0.222737564604946)
  u_pot%p3(5,4,2,3) = cmplx(   0.228402304294797,  -0.222737564604946)
  u_pot%p3(5,4,3,2) = cmplx(  -0.228402304294797,   0.222737564604946)
  u_pot%p3(5,4,3,3) = cmplx(   0.076134101431599,  -0.222737564604946)
  u_pot%p3(5,4,4,1) = cmplx(   0.257195185766614,   0.000000000000000)
  u_pot%p3(5,4,4,5) = cmplx(   0.450026164311878,   0.000000000000000)
  u_pot%p3(5,4,5,2) = cmplx(   0.186489700532716,   0.000000000000000)
  u_pot%p3(5,4,5,3) = cmplx(   0.186489700532716,   0.000000000000000)
  u_pot%p3(5,5,1,1) = cmplx(   0.821559392987964,   0.000000000000000)
  u_pot%p3(5,5,2,2) = cmplx(   0.141156462585034,  -0.087912087912088)
  u_pot%p3(5,5,2,3) = cmplx(   0.172553636839351,   0.087912087912088)
  u_pot%p3(5,5,3,2) = cmplx(   0.172553636839351,   0.087912087912088)
  u_pot%p3(5,5,3,3) = cmplx(   0.141156462585034,  -0.087912087912088)
  u_pot%p3(5,5,4,2) = cmplx(   0.186489700532716,   0.000000000000000)
  u_pot%p3(5,5,4,3) = cmplx(   0.186489700532716,   0.000000000000000)
  u_pot%p3(5,5,4,4) = cmplx(   0.049973835688122,   0.000000000000000)
  u_pot%p3(5,5,5,5) = cmplx(   0.500000000000000,   0.000000000000000)
  u_pot%p4(1,1,1,1) = cmplx(   0.367346938775510,   0.000000000000000)
  u_pot%p4(1,1,2,2) = cmplx(  -0.263300191871620,   0.000000000000000)
  u_pot%p4(1,1,2,3) = cmplx(   0.085731728588871,   0.000000000000000)
  u_pot%p4(1,1,3,2) = cmplx(   0.085731728588871,   0.000000000000000)
  u_pot%p4(1,1,3,3) = cmplx(  -0.263300191871620,   0.000000000000000)
  u_pot%p4(1,1,4,4) = cmplx(  -0.177568463282749,   0.000000000000000)
  u_pot%p4(1,1,5,5) = cmplx(   0.336821908250480,   0.000000000000000)
  u_pot%p4(1,2,1,2) = cmplx(  -0.219780219780220,  -0.091182626896913)
  u_pot%p4(1,2,1,3) = cmplx(  -0.043956043956044,   0.091182626896913)
  u_pot%p4(1,2,1,4) = cmplx(  -0.006166987451479,   0.000000000000000)
  u_pot%p4(1,2,2,1) = cmplx(   0.095543345543346,  -0.087912087912088)
  u_pot%p4(1,2,2,2) = cmplx(   0.000000000000000,   0.052499747452378)
  u_pot%p4(1,2,2,3) = cmplx(   0.000000000000000,  -0.055170131351278)
  u_pot%p4(1,2,2,4) = cmplx(  -0.076134101431599,  -0.148491709736631)
  u_pot%p4(1,2,2,5) = cmplx(   0.091182626896913,  -0.087912087912088)
  u_pot%p4(1,2,3,1) = cmplx(  -0.086821908250480,   0.087912087912088)
  u_pot%p4(1,2,3,2) = cmplx(   0.000000000000000,   0.055170131351278)
  u_pot%p4(1,2,3,3) = cmplx(   0.000000000000000,  -0.052499747452378)
  u_pot%p4(1,2,3,4) = cmplx(  -0.228402304294797,  -0.148491709736631)
  u_pot%p4(1,2,3,5) = cmplx(  -0.084641548927263,   0.087912087912088)
  u_pot%p4(1,2,4,1) = cmplx(   0.186489700532716,   0.000000000000000)
  u_pot%p4(1,2,4,2) = cmplx(  -0.001888246563284,   0.000000000000000)
  u_pot%p4(1,2,4,3) = cmplx(  -0.005664739689851,   0.000000000000000)
  u_pot%p4(1,2,4,5) = cmplx(   0.186489700532716,   0.000000000000000)
  u_pot%p4(1,2,5,2) = cmplx(   0.219780219780220,   0.087912087912088)
  u_pot%p4(1,2,5,3) = cmplx(   0.043956043956044,  -0.087912087912088)
  u_pot%p4(1,2,5,4) = cmplx(   0.004625240588609,   0.000000000000000)
  u_pot%p4(1,3,1,2) = cmplx(  -0.043956043956044,   0.091182626896913)
  u_pot%p4(1,3,1,3) = cmplx(  -0.219780219780220,  -0.091182626896913)
  u_pot%p4(1,3,1,4) = cmplx(  -0.006166987451479,   0.000000000000000)
  u_pot%p4(1,3,2,1) = cmplx(  -0.086821908250480,   0.087912087912088)
  u_pot%p4(1,3,2,2) = cmplx(   0.000000000000000,   0.052499747452378)
  u_pot%p4(1,3,2,3) = cmplx(   0.000000000000000,  -0.055170131351278)
  u_pot%p4(1,3,2,4) = cmplx(   0.228402304294797,   0.148491709736631)
  u_pot%p4(1,3,2,5) = cmplx(  -0.084641548927263,   0.087912087912088)
  u_pot%p4(1,3,3,1) = cmplx(   0.095543345543346,  -0.087912087912088)
  u_pot%p4(1,3,3,2) = cmplx(   0.000000000000000,   0.055170131351278)
  u_pot%p4(1,3,3,3) = cmplx(   0.000000000000000,  -0.052499747452378)
  u_pot%p4(1,3,3,4) = cmplx(   0.076134101431599,   0.148491709736631)
  u_pot%p4(1,3,3,5) = cmplx(   0.091182626896913,  -0.087912087912088)
  u_pot%p4(1,3,4,1) = cmplx(   0.186489700532716,   0.000000000000000)
  u_pot%p4(1,3,4,2) = cmplx(   0.005664739689851,   0.000000000000000)
  u_pot%p4(1,3,4,3) = cmplx(   0.001888246563284,   0.000000000000000)
  u_pot%p4(1,3,4,5) = cmplx(   0.186489700532716,   0.000000000000000)
  u_pot%p4(1,3,5,2) = cmplx(   0.043956043956044,  -0.087912087912088)
  u_pot%p4(1,3,5,3) = cmplx(   0.219780219780220,   0.087912087912088)
  u_pot%p4(1,3,5,4) = cmplx(   0.004625240588609,   0.000000000000000)
  u_pot%p4(1,4,1,2) = cmplx(  -0.192656687984194,   0.000000000000000)
  u_pot%p4(1,4,1,3) = cmplx(  -0.192656687984194,   0.000000000000000)
  u_pot%p4(1,4,2,2) = cmplx(   0.000000000000000,   0.074245854868315)
  u_pot%p4(1,4,2,3) = cmplx(   0.000000000000000,  -0.074245854868315)
  u_pot%p4(1,4,2,4) = cmplx(   0.104999494904756,   0.000000000000000)
  u_pot%p4(1,4,3,2) = cmplx(   0.000000000000000,   0.074245854868315)
  u_pot%p4(1,4,3,3) = cmplx(   0.000000000000000,  -0.074245854868315)
  u_pot%p4(1,4,3,4) = cmplx(  -0.104999494904756,   0.000000000000000)
  u_pot%p4(1,4,4,1) = cmplx(   0.272457701029130,   0.000000000000000)
  u_pot%p4(1,4,4,5) = cmplx(   0.257195185766614,   0.000000000000000)
  u_pot%p4(1,4,5,2) = cmplx(   0.181864459944107,   0.000000000000000)
  u_pot%p4(1,4,5,3) = cmplx(   0.181864459944107,   0.000000000000000)
  u_pot%p4(1,5,1,5) = cmplx(  -0.015262515262515,   0.000000000000000)
  u_pot%p4(1,5,5,1) = cmplx(   0.015262515262515,   0.000000000000000)
  u_pot%p4(2,1,1,2) = cmplx(   0.227411477411477,  -0.087912087912088)
  u_pot%p4(2,1,1,3) = cmplx(   0.045046223617652,   0.087912087912088)
  u_pot%p4(2,1,2,1) = cmplx(   0.087912087912088,   0.091182626896913)
  u_pot%p4(2,1,2,2) = cmplx(   0.000000000000000,  -0.052499747452378)
  u_pot%p4(2,1,2,3) = cmplx(   0.000000000000000,   0.052499747452378)
  u_pot%p4(2,1,2,4) = cmplx(   0.074245854868315,  -0.152268202863198)
  u_pot%p4(2,1,2,5) = cmplx(   0.087912087912088,   0.087912087912088)
  u_pot%p4(2,1,3,1) = cmplx(  -0.087912087912088,  -0.091182626896913)
  u_pot%p4(2,1,3,2) = cmplx(   0.000000000000000,  -0.052499747452378)
  u_pot%p4(2,1,3,3) = cmplx(   0.000000000000000,   0.052499747452378)
  u_pot%p4(2,1,3,4) = cmplx(   0.222737564604946,  -0.152268202863198)
  u_pot%p4(2,1,3,5) = cmplx(  -0.087912087912088,  -0.087912087912088)
  u_pot%p4(2,1,4,1) = cmplx(  -0.192656687984194,   0.000000000000000)
  u_pot%p4(2,1,4,2) = cmplx(   0.000000000000000,  -0.003776493126567)
  u_pot%p4(2,1,4,3) = cmplx(   0.000000000000000,  -0.003776493126567)
  u_pot%p4(2,1,4,5) = cmplx(  -0.181864459944107,   0.000000000000000)
  u_pot%p4(2,1,5,2) = cmplx(  -0.216509680795395,   0.087912087912088)
  u_pot%p4(2,1,5,3) = cmplx(  -0.040685504971219,  -0.087912087912088)
  u_pot%p4(2,2,1,1) = cmplx(  -0.263300191871620,   0.000000000000000)
  u_pot%p4(2,2,1,2) = cmplx(   0.000000000000000,   0.106334686854206)
  u_pot%p4(2,2,1,3) = cmplx(   0.000000000000000,   0.106334686854206)
  u_pot%p4(2,2,1,4) = cmplx(   0.000000000000000,  -0.001888246563284)
  u_pot%p4(2,2,2,1) = cmplx(   0.000000000000000,   0.001335191949450)
  u_pot%p4(2,2,2,2) = cmplx(   0.248735391592534,   0.000000000000000)
  u_pot%p4(2,2,2,3) = cmplx(  -0.112070469213326,   0.000000000000000)
  u_pot%p4(2,2,2,4) = cmplx(   0.000000000000000,  -0.066788474099514)
  u_pot%p4(2,2,2,5) = cmplx(   0.000000000000000,  -0.001335191949450)
  u_pot%p4(2,2,3,1) = cmplx(   0.000000000000000,   0.001335191949450)
  u_pot%p4(2,2,3,2) = cmplx(  -0.024158381301238,   0.000000000000000)
  u_pot%p4(2,2,3,3) = cmplx(   0.147741147741148,   0.000000000000000)
  u_pot%p4(2,2,3,4) = cmplx(   0.000000000000000,   0.066788474099514)
  u_pot%p4(2,2,3,5) = cmplx(   0.000000000000000,  -0.001335191949450)
  u_pot%p4(2,2,4,1) = cmplx(   0.000000000000000,  -0.150379956299914)
  u_pot%p4(2,2,4,2) = cmplx(   0.000000000000000,   0.004625240588609)
  u_pot%p4(2,2,4,3) = cmplx(   0.000000000000000,  -0.004625240588609)
  u_pot%p4(2,2,4,4) = cmplx(   0.130123844409559,   0.000000000000000)
  u_pot%p4(2,2,4,5) = cmplx(   0.000000000000000,  -0.150379956299914)
  u_pot%p4(2,2,5,2) = cmplx(   0.000000000000000,  -0.106334686854206)
  u_pot%p4(2,2,5,3) = cmplx(   0.000000000000000,  -0.106334686854206)
  u_pot%p4(2,2,5,4) = cmplx(   0.000000000000000,  -0.001888246563284)
  u_pot%p4(2,2,5,5) = cmplx(  -0.263300191871620,   0.000000000000000)
  u_pot%p4(2,3,1,1) = cmplx(   0.085731728588871,   0.000000000000000)
  u_pot%p4(2,3,1,2) = cmplx(   0.000000000000000,  -0.001335191949450)
  u_pot%p4(2,3,1,3) = cmplx(   0.000000000000000,  -0.001335191949450)
  u_pot%p4(2,3,1,4) = cmplx(   0.000000000000000,   0.001888246563284)
  u_pot%p4(2,3,2,1) = cmplx(   0.000000000000000,  -0.001335191949450)
  u_pot%p4(2,3,2,2) = cmplx(  -0.068114425257282,   0.000000000000000)
  u_pot%p4(2,3,2,3) = cmplx(   0.026338740624455,   0.000000000000000)
  u_pot%p4(2,3,2,4) = cmplx(   0.000000000000000,  -0.004625240588609)
  u_pot%p4(2,3,2,5) = cmplx(   0.000000000000000,   0.001335191949450)
  u_pot%p4(2,3,3,1) = cmplx(   0.000000000000000,  -0.001335191949450)
  u_pot%p4(2,3,3,2) = cmplx(   0.039420896563754,   0.000000000000000)
  u_pot%p4(2,3,3,3) = cmplx(  -0.068114425257282,   0.000000000000000)
  u_pot%p4(2,3,3,4) = cmplx(   0.000000000000000,   0.004625240588609)
  u_pot%p4(2,3,3,5) = cmplx(   0.000000000000000,   0.001335191949450)
  u_pot%p4(2,3,4,1) = cmplx(   0.000000000000000,   0.001888246563284)
  u_pot%p4(2,3,4,2) = cmplx(   0.000000000000000,  -0.004625240588609)
  u_pot%p4(2,3,4,3) = cmplx(   0.000000000000000,   0.004625240588609)
  u_pot%p4(2,3,4,4) = cmplx(  -0.035234606663178,   0.000000000000000)
  u_pot%p4(2,3,4,5) = cmplx(   0.000000000000000,   0.001888246563284)
  u_pot%p4(2,3,5,2) = cmplx(   0.000000000000000,   0.001335191949450)
  u_pot%p4(2,3,5,3) = cmplx(   0.000000000000000,   0.001335191949450)
  u_pot%p4(2,3,5,4) = cmplx(   0.000000000000000,   0.001888246563284)
  u_pot%p4(2,3,5,5) = cmplx(   0.085731728588871,   0.000000000000000)
  u_pot%p4(2,4,1,2) = cmplx(  -0.076134101431599,   0.003776493126567)
  u_pot%p4(2,4,1,3) = cmplx(   0.228402304294797,  -0.003776493126567)
  u_pot%p4(2,4,1,4) = cmplx(  -0.002670383898900,   0.000000000000000)
  u_pot%p4(2,4,2,1) = cmplx(  -0.001888246563284,  -0.152268202863198)
  u_pot%p4(2,4,2,2) = cmplx(   0.000000000000000,  -0.035706857344062)
  u_pot%p4(2,4,2,3) = cmplx(   0.000000000000000,   0.026456376166844)
  u_pot%p4(2,4,2,4) = cmplx(  -0.307692307692308,   0.000000000000000)
  u_pot%p4(2,4,2,5) = cmplx(   0.001888246563284,  -0.152268202863198)
  u_pot%p4(2,4,3,1) = cmplx(   0.005664739689851,   0.152268202863198)
  u_pot%p4(2,4,3,2) = cmplx(   0.000000000000000,  -0.026456376166844)
  u_pot%p4(2,4,3,3) = cmplx(   0.000000000000000,   0.035706857344062)
  u_pot%p4(2,4,3,4) = cmplx(  -0.219780219780220,   0.000000000000000)
  u_pot%p4(2,4,3,5) = cmplx(  -0.005664739689851,   0.152268202863198)
  u_pot%p4(2,4,4,1) = cmplx(  -0.107669878803656,   0.000000000000000)
  u_pot%p4(2,4,4,2) = cmplx(   0.015262515262515,   0.000000000000000)
  u_pot%p4(2,4,4,3) = cmplx(   0.002180359323216,   0.000000000000000)
  u_pot%p4(2,4,4,5) = cmplx(  -0.107669878803656,   0.000000000000000)
  u_pot%p4(2,4,5,2) = cmplx(   0.076134101431599,   0.003776493126567)
  u_pot%p4(2,4,5,3) = cmplx(  -0.228402304294797,  -0.003776493126567)
  u_pot%p4(2,4,5,4) = cmplx(  -0.002670383898900,   0.000000000000000)
  u_pot%p4(2,5,1,2) = cmplx(  -0.216509680795395,  -0.087912087912088)
  u_pot%p4(2,5,1,3) = cmplx(  -0.040685504971219,   0.087912087912088)
  u_pot%p4(2,5,2,1) = cmplx(   0.087912087912088,  -0.087912087912088)
  u_pot%p4(2,5,2,2) = cmplx(   0.000000000000000,   0.052499747452378)
  u_pot%p4(2,5,2,3) = cmplx(   0.000000000000000,  -0.052499747452378)
  u_pot%p4(2,5,2,4) = cmplx(  -0.074245854868315,  -0.152268202863198)
  u_pot%p4(2,5,2,5) = cmplx(   0.087912087912088,  -0.091182626896913)
  u_pot%p4(2,5,3,1) = cmplx(  -0.087912087912088,   0.087912087912088)
  u_pot%p4(2,5,3,2) = cmplx(   0.000000000000000,   0.052499747452378)
  u_pot%p4(2,5,3,3) = cmplx(   0.000000000000000,  -0.052499747452378)
  u_pot%p4(2,5,3,4) = cmplx(  -0.222737564604946,  -0.152268202863198)
  u_pot%p4(2,5,3,5) = cmplx(  -0.087912087912088,   0.091182626896913)
  u_pot%p4(2,5,4,1) = cmplx(   0.181864459944107,   0.000000000000000)
  u_pot%p4(2,5,4,2) = cmplx(   0.000000000000000,  -0.003776493126567)
  u_pot%p4(2,5,4,3) = cmplx(   0.000000000000000,  -0.003776493126567)
  u_pot%p4(2,5,4,5) = cmplx(   0.192656687984194,   0.000000000000000)
  u_pot%p4(2,5,5,2) = cmplx(   0.227411477411477,   0.087912087912088)
  u_pot%p4(2,5,5,3) = cmplx(   0.045046223617652,  -0.087912087912088)
  u_pot%p4(3,1,1,2) = cmplx(   0.045046223617652,   0.087912087912088)
  u_pot%p4(3,1,1,3) = cmplx(   0.227411477411477,  -0.087912087912088)
  u_pot%p4(3,1,2,1) = cmplx(  -0.087912087912088,  -0.091182626896913)
  u_pot%p4(3,1,2,2) = cmplx(   0.000000000000000,  -0.052499747452378)
  u_pot%p4(3,1,2,3) = cmplx(   0.000000000000000,   0.052499747452378)
  u_pot%p4(3,1,2,4) = cmplx(  -0.222737564604946,   0.152268202863198)
  u_pot%p4(3,1,2,5) = cmplx(  -0.087912087912088,  -0.087912087912088)
  u_pot%p4(3,1,3,1) = cmplx(   0.087912087912088,   0.091182626896913)
  u_pot%p4(3,1,3,2) = cmplx(   0.000000000000000,  -0.052499747452378)
  u_pot%p4(3,1,3,3) = cmplx(   0.000000000000000,   0.052499747452378)
  u_pot%p4(3,1,3,4) = cmplx(  -0.074245854868315,   0.152268202863198)
  u_pot%p4(3,1,3,5) = cmplx(   0.087912087912088,   0.087912087912088)
  u_pot%p4(3,1,4,1) = cmplx(  -0.192656687984194,   0.000000000000000)
  u_pot%p4(3,1,4,2) = cmplx(   0.000000000000000,   0.003776493126567)
  u_pot%p4(3,1,4,3) = cmplx(   0.000000000000000,   0.003776493126567)
  u_pot%p4(3,1,4,5) = cmplx(  -0.181864459944107,   0.000000000000000)
  u_pot%p4(3,1,5,2) = cmplx(  -0.040685504971219,  -0.087912087912088)
  u_pot%p4(3,1,5,3) = cmplx(  -0.216509680795395,   0.087912087912088)
  u_pot%p4(3,2,1,1) = cmplx(   0.085731728588871,   0.000000000000000)
  u_pot%p4(3,2,1,2) = cmplx(   0.000000000000000,   0.001335191949450)
  u_pot%p4(3,2,1,3) = cmplx(   0.000000000000000,   0.001335191949450)
  u_pot%p4(3,2,1,4) = cmplx(   0.000000000000000,  -0.001888246563284)
  u_pot%p4(3,2,2,1) = cmplx(   0.000000000000000,   0.001335191949450)
  u_pot%p4(3,2,2,2) = cmplx(  -0.068114425257282,   0.000000000000000)
  u_pot%p4(3,2,2,3) = cmplx(   0.039420896563754,   0.000000000000000)
  u_pot%p4(3,2,2,4) = cmplx(   0.000000000000000,   0.004625240588609)
  u_pot%p4(3,2,2,5) = cmplx(   0.000000000000000,  -0.001335191949450)
  u_pot%p4(3,2,3,1) = cmplx(   0.000000000000000,   0.001335191949450)
  u_pot%p4(3,2,3,2) = cmplx(   0.026338740624455,   0.000000000000000)
  u_pot%p4(3,2,3,3) = cmplx(  -0.068114425257282,   0.000000000000000)
  u_pot%p4(3,2,3,4) = cmplx(   0.000000000000000,  -0.004625240588609)
  u_pot%p4(3,2,3,5) = cmplx(   0.000000000000000,  -0.001335191949450)
  u_pot%p4(3,2,4,1) = cmplx(   0.000000000000000,  -0.001888246563284)
  u_pot%p4(3,2,4,2) = cmplx(   0.000000000000000,   0.004625240588609)
  u_pot%p4(3,2,4,3) = cmplx(   0.000000000000000,  -0.004625240588609)
  u_pot%p4(3,2,4,4) = cmplx(  -0.035234606663178,   0.000000000000000)
  u_pot%p4(3,2,4,5) = cmplx(   0.000000000000000,  -0.001888246563284)
  u_pot%p4(3,2,5,2) = cmplx(   0.000000000000000,  -0.001335191949450)
  u_pot%p4(3,2,5,3) = cmplx(   0.000000000000000,  -0.001335191949450)
  u_pot%p4(3,2,5,4) = cmplx(   0.000000000000000,  -0.001888246563284)
  u_pot%p4(3,2,5,5) = cmplx(   0.085731728588871,   0.000000000000000)
  u_pot%p4(3,3,1,1) = cmplx(  -0.263300191871620,   0.000000000000000)
  u_pot%p4(3,3,1,2) = cmplx(   0.000000000000000,  -0.106334686854206)
  u_pot%p4(3,3,1,3) = cmplx(   0.000000000000000,  -0.106334686854206)
  u_pot%p4(3,3,1,4) = cmplx(   0.000000000000000,   0.001888246563284)
  u_pot%p4(3,3,2,1) = cmplx(   0.000000000000000,  -0.001335191949450)
  u_pot%p4(3,3,2,2) = cmplx(   0.147741147741148,   0.000000000000000)
  u_pot%p4(3,3,2,3) = cmplx(  -0.024158381301238,   0.000000000000000)
  u_pot%p4(3,3,2,4) = cmplx(   0.000000000000000,   0.066788474099514)
  u_pot%p4(3,3,2,5) = cmplx(   0.000000000000000,   0.001335191949450)
  u_pot%p4(3,3,3,1) = cmplx(   0.000000000000000,  -0.001335191949450)
  u_pot%p4(3,3,3,2) = cmplx(  -0.112070469213326,   0.000000000000000)
  u_pot%p4(3,3,3,3) = cmplx(   0.248735391592534,   0.000000000000000)
  u_pot%p4(3,3,3,4) = cmplx(   0.000000000000000,  -0.066788474099514)
  u_pot%p4(3,3,3,5) = cmplx(   0.000000000000000,   0.001335191949450)
  u_pot%p4(3,3,4,1) = cmplx(   0.000000000000000,   0.150379956299914)
  u_pot%p4(3,3,4,2) = cmplx(   0.000000000000000,  -0.004625240588609)
  u_pot%p4(3,3,4,3) = cmplx(   0.000000000000000,   0.004625240588609)
  u_pot%p4(3,3,4,4) = cmplx(   0.130123844409559,   0.000000000000000)
  u_pot%p4(3,3,4,5) = cmplx(   0.000000000000000,   0.150379956299914)
  u_pot%p4(3,3,5,2) = cmplx(   0.000000000000000,   0.106334686854206)
  u_pot%p4(3,3,5,3) = cmplx(   0.000000000000000,   0.106334686854206)
  u_pot%p4(3,3,5,4) = cmplx(   0.000000000000000,   0.001888246563284)
  u_pot%p4(3,3,5,5) = cmplx(  -0.263300191871620,   0.000000000000000)
  u_pot%p4(3,4,1,2) = cmplx(  -0.228402304294797,   0.003776493126567)
  u_pot%p4(3,4,1,3) = cmplx(   0.076134101431599,  -0.003776493126567)
  u_pot%p4(3,4,1,4) = cmplx(   0.002670383898900,   0.000000000000000)
  u_pot%p4(3,4,2,1) = cmplx(  -0.005664739689851,  -0.152268202863198)
  u_pot%p4(3,4,2,2) = cmplx(   0.000000000000000,   0.035706857344062)
  u_pot%p4(3,4,2,3) = cmplx(   0.000000000000000,  -0.026456376166844)
  u_pot%p4(3,4,2,4) = cmplx(  -0.219780219780220,   0.000000000000000)
  u_pot%p4(3,4,2,5) = cmplx(   0.005664739689851,  -0.152268202863198)
  u_pot%p4(3,4,3,1) = cmplx(   0.001888246563284,   0.152268202863198)
  u_pot%p4(3,4,3,2) = cmplx(   0.000000000000000,   0.026456376166844)
  u_pot%p4(3,4,3,3) = cmplx(   0.000000000000000,  -0.035706857344062)
  u_pot%p4(3,4,3,4) = cmplx(  -0.307692307692308,   0.000000000000000)
  u_pot%p4(3,4,3,5) = cmplx(  -0.001888246563284,   0.152268202863198)
  u_pot%p4(3,4,4,1) = cmplx(   0.107669878803656,   0.000000000000000)
  u_pot%p4(3,4,4,2) = cmplx(   0.002180359323216,   0.000000000000000)
  u_pot%p4(3,4,4,3) = cmplx(   0.015262515262515,   0.000000000000000)
  u_pot%p4(3,4,4,5) = cmplx(   0.107669878803656,   0.000000000000000)
  u_pot%p4(3,4,5,2) = cmplx(   0.228402304294797,   0.003776493126567)
  u_pot%p4(3,4,5,3) = cmplx(  -0.076134101431599,  -0.003776493126567)
  u_pot%p4(3,4,5,4) = cmplx(   0.002670383898900,   0.000000000000000)
  u_pot%p4(3,5,1,2) = cmplx(  -0.040685504971219,   0.087912087912088)
  u_pot%p4(3,5,1,3) = cmplx(  -0.216509680795395,  -0.087912087912088)
  u_pot%p4(3,5,2,1) = cmplx(  -0.087912087912088,   0.087912087912088)
  u_pot%p4(3,5,2,2) = cmplx(   0.000000000000000,   0.052499747452378)
  u_pot%p4(3,5,2,3) = cmplx(   0.000000000000000,  -0.052499747452378)
  u_pot%p4(3,5,2,4) = cmplx(   0.222737564604946,   0.152268202863198)
  u_pot%p4(3,5,2,5) = cmplx(  -0.087912087912088,   0.091182626896913)
  u_pot%p4(3,5,3,1) = cmplx(   0.087912087912088,  -0.087912087912088)
  u_pot%p4(3,5,3,2) = cmplx(   0.000000000000000,   0.052499747452378)
  u_pot%p4(3,5,3,3) = cmplx(   0.000000000000000,  -0.052499747452378)
  u_pot%p4(3,5,3,4) = cmplx(   0.074245854868315,   0.152268202863198)
  u_pot%p4(3,5,3,5) = cmplx(   0.087912087912088,  -0.091182626896913)
  u_pot%p4(3,5,4,1) = cmplx(   0.181864459944107,   0.000000000000000)
  u_pot%p4(3,5,4,2) = cmplx(   0.000000000000000,   0.003776493126567)
  u_pot%p4(3,5,4,3) = cmplx(   0.000000000000000,   0.003776493126567)
  u_pot%p4(3,5,4,5) = cmplx(   0.192656687984194,   0.000000000000000)
  u_pot%p4(3,5,5,2) = cmplx(   0.045046223617652,  -0.087912087912088)
  u_pot%p4(3,5,5,3) = cmplx(   0.227411477411477,   0.087912087912088)
  u_pot%p4(4,1,1,2) = cmplx(   0.186489700532716,   0.000000000000000)
  u_pot%p4(4,1,1,3) = cmplx(   0.186489700532716,   0.000000000000000)
  u_pot%p4(4,1,1,4) = cmplx(   0.008721437292866,   0.000000000000000)
  u_pot%p4(4,1,2,1) = cmplx(  -0.006166987451479,   0.000000000000000)
  u_pot%p4(4,1,2,2) = cmplx(   0.000000000000000,  -0.074245854868315)
  u_pot%p4(4,1,2,3) = cmplx(   0.000000000000000,   0.078022347994883)
  u_pot%p4(4,1,2,4) = cmplx(  -0.107669878803656,   0.000000000000000)
  u_pot%p4(4,1,2,5) = cmplx(  -0.004625240588609,   0.000000000000000)
  u_pot%p4(4,1,3,1) = cmplx(  -0.006166987451479,   0.000000000000000)
  u_pot%p4(4,1,3,2) = cmplx(   0.000000000000000,  -0.078022347994883)
  u_pot%p4(4,1,3,3) = cmplx(   0.000000000000000,   0.074245854868315)
  u_pot%p4(4,1,3,4) = cmplx(   0.107669878803656,   0.000000000000000)
  u_pot%p4(4,1,3,5) = cmplx(  -0.004625240588609,   0.000000000000000)
  u_pot%p4(4,1,4,1) = cmplx(  -0.263736263736264,   0.000000000000000)
  u_pot%p4(4,1,4,2) = cmplx(  -0.002670383898900,   0.000000000000000)
  u_pot%p4(4,1,4,3) = cmplx(   0.002670383898900,   0.000000000000000)
  u_pot%p4(4,1,4,5) = cmplx(  -0.263736263736264,   0.000000000000000)
  u_pot%p4(4,1,5,2) = cmplx(  -0.186489700532716,   0.000000000000000)
  u_pot%p4(4,1,5,3) = cmplx(  -0.186489700532716,   0.000000000000000)
  u_pot%p4(4,1,5,4) = cmplx(  -0.006541077969649,   0.000000000000000)
  u_pot%p4(4,2,1,2) = cmplx(   0.074245854868315,   0.000000000000000)
  u_pot%p4(4,2,1,3) = cmplx(  -0.222737564604946,   0.000000000000000)
  u_pot%p4(4,2,2,1) = cmplx(   0.000000000000000,   0.148491709736631)
  u_pot%p4(4,2,2,2) = cmplx(   0.000000000000000,   0.035706857344062)
  u_pot%p4(4,2,2,3) = cmplx(   0.000000000000000,  -0.035706857344062)
  u_pot%p4(4,2,2,4) = cmplx(   0.322954822954823,   0.000000000000000)
  u_pot%p4(4,2,2,5) = cmplx(   0.000000000000000,   0.148491709736631)
  u_pot%p4(4,2,3,1) = cmplx(   0.000000000000000,  -0.148491709736631)
  u_pot%p4(4,2,3,2) = cmplx(   0.000000000000000,   0.035706857344062)
  u_pot%p4(4,2,3,3) = cmplx(   0.000000000000000,  -0.035706857344062)
  u_pot%p4(4,2,3,4) = cmplx(   0.221960579103436,   0.000000000000000)
  u_pot%p4(4,2,3,5) = cmplx(   0.000000000000000,  -0.148491709736631)
  u_pot%p4(4,2,4,1) = cmplx(   0.104999494904756,   0.000000000000000)
  u_pot%p4(4,2,4,5) = cmplx(   0.104999494904756,   0.000000000000000)
  u_pot%p4(4,2,5,2) = cmplx(  -0.074245854868315,   0.000000000000000)
  u_pot%p4(4,2,5,3) = cmplx(   0.222737564604946,   0.000000000000000)
  u_pot%p4(4,3,1,2) = cmplx(   0.222737564604946,   0.000000000000000)
  u_pot%p4(4,3,1,3) = cmplx(  -0.074245854868315,   0.000000000000000)
  u_pot%p4(4,3,2,1) = cmplx(   0.000000000000000,   0.148491709736631)
  u_pot%p4(4,3,2,2) = cmplx(   0.000000000000000,  -0.035706857344062)
  u_pot%p4(4,3,2,3) = cmplx(   0.000000000000000,   0.035706857344062)
  u_pot%p4(4,3,2,4) = cmplx(   0.221960579103436,   0.000000000000000)
  u_pot%p4(4,3,2,5) = cmplx(   0.000000000000000,   0.148491709736631)
  u_pot%p4(4,3,3,1) = cmplx(   0.000000000000000,  -0.148491709736631)
  u_pot%p4(4,3,3,2) = cmplx(   0.000000000000000,  -0.035706857344062)
  u_pot%p4(4,3,3,3) = cmplx(   0.000000000000000,   0.035706857344062)
  u_pot%p4(4,3,3,4) = cmplx(   0.322954822954823,   0.000000000000000)
  u_pot%p4(4,3,3,5) = cmplx(   0.000000000000000,  -0.148491709736631)
  u_pot%p4(4,3,4,1) = cmplx(  -0.104999494904756,   0.000000000000000)
  u_pot%p4(4,3,4,5) = cmplx(  -0.104999494904756,   0.000000000000000)
  u_pot%p4(4,3,5,2) = cmplx(  -0.222737564604946,   0.000000000000000)
  u_pot%p4(4,3,5,3) = cmplx(   0.074245854868315,   0.000000000000000)
  u_pot%p4(4,4,1,1) = cmplx(  -0.177568463282749,   0.000000000000000)
  u_pot%p4(4,4,2,2) = cmplx(   0.130123844409559,   0.000000000000000)
  u_pot%p4(4,4,2,3) = cmplx(  -0.035234606663178,   0.000000000000000)
  u_pot%p4(4,4,3,2) = cmplx(  -0.035234606663178,   0.000000000000000)
  u_pot%p4(4,4,3,3) = cmplx(   0.130123844409559,   0.000000000000000)
  u_pot%p4(4,4,4,4) = cmplx(   0.094889237746381,   0.000000000000000)
  u_pot%p4(4,4,5,5) = cmplx(  -0.177568463282749,   0.000000000000000)
  u_pot%p4(4,5,1,2) = cmplx(   0.186489700532716,   0.000000000000000)
  u_pot%p4(4,5,1,3) = cmplx(   0.186489700532716,   0.000000000000000)
  u_pot%p4(4,5,1,4) = cmplx(  -0.006541077969649,   0.000000000000000)
  u_pot%p4(4,5,2,1) = cmplx(   0.004625240588609,   0.000000000000000)
  u_pot%p4(4,5,2,2) = cmplx(   0.000000000000000,  -0.074245854868315)
  u_pot%p4(4,5,2,3) = cmplx(   0.000000000000000,   0.078022347994883)
  u_pot%p4(4,5,2,4) = cmplx(  -0.107669878803656,   0.000000000000000)
  u_pot%p4(4,5,2,5) = cmplx(   0.006166987451479,   0.000000000000000)
  u_pot%p4(4,5,3,1) = cmplx(   0.004625240588609,   0.000000000000000)
  u_pot%p4(4,5,3,2) = cmplx(   0.000000000000000,  -0.078022347994883)
  u_pot%p4(4,5,3,3) = cmplx(   0.000000000000000,   0.074245854868315)
  u_pot%p4(4,5,3,4) = cmplx(   0.107669878803656,   0.000000000000000)
  u_pot%p4(4,5,3,5) = cmplx(   0.006166987451479,   0.000000000000000)
  u_pot%p4(4,5,4,1) = cmplx(  -0.263736263736264,   0.000000000000000)
  u_pot%p4(4,5,4,2) = cmplx(  -0.002670383898900,   0.000000000000000)
  u_pot%p4(4,5,4,3) = cmplx(   0.002670383898900,   0.000000000000000)
  u_pot%p4(4,5,4,5) = cmplx(  -0.263736263736264,   0.000000000000000)
  u_pot%p4(4,5,5,2) = cmplx(  -0.186489700532716,   0.000000000000000)
  u_pot%p4(4,5,5,3) = cmplx(  -0.186489700532716,   0.000000000000000)
  u_pot%p4(4,5,5,4) = cmplx(   0.008721437292866,   0.000000000000000)
  u_pot%p4(5,1,1,5) = cmplx(   0.015262515262515,   0.000000000000000)
  u_pot%p4(5,1,5,1) = cmplx(  -0.015262515262515,   0.000000000000000)
  u_pot%p4(5,2,1,2) = cmplx(   0.219780219780220,  -0.087912087912088)
  u_pot%p4(5,2,1,3) = cmplx(   0.043956043956044,   0.087912087912088)
  u_pot%p4(5,2,1,4) = cmplx(  -0.004625240588609,   0.000000000000000)
  u_pot%p4(5,2,2,1) = cmplx(   0.091182626896913,   0.087912087912088)
  u_pot%p4(5,2,2,2) = cmplx(   0.000000000000000,  -0.052499747452378)
  u_pot%p4(5,2,2,3) = cmplx(   0.000000000000000,   0.055170131351278)
  u_pot%p4(5,2,2,4) = cmplx(   0.076134101431599,  -0.148491709736631)
  u_pot%p4(5,2,2,5) = cmplx(   0.095543345543346,   0.087912087912088)
  u_pot%p4(5,2,3,1) = cmplx(  -0.084641548927263,  -0.087912087912088)
  u_pot%p4(5,2,3,2) = cmplx(   0.000000000000000,  -0.055170131351278)
  u_pot%p4(5,2,3,3) = cmplx(   0.000000000000000,   0.052499747452378)
  u_pot%p4(5,2,3,4) = cmplx(   0.228402304294797,  -0.148491709736631)
  u_pot%p4(5,2,3,5) = cmplx(  -0.086821908250480,  -0.087912087912088)
  u_pot%p4(5,2,4,1) = cmplx(  -0.186489700532716,   0.000000000000000)
  u_pot%p4(5,2,4,2) = cmplx(   0.001888246563284,   0.000000000000000)
  u_pot%p4(5,2,4,3) = cmplx(   0.005664739689851,   0.000000000000000)
  u_pot%p4(5,2,4,5) = cmplx(  -0.186489700532716,   0.000000000000000)
  u_pot%p4(5,2,5,2) = cmplx(  -0.219780219780220,   0.091182626896913)
  u_pot%p4(5,2,5,3) = cmplx(  -0.043956043956044,  -0.091182626896913)
  u_pot%p4(5,2,5,4) = cmplx(   0.006166987451479,   0.000000000000000)
  u_pot%p4(5,3,1,2) = cmplx(   0.043956043956044,   0.087912087912088)
  u_pot%p4(5,3,1,3) = cmplx(   0.219780219780220,  -0.087912087912088)
  u_pot%p4(5,3,1,4) = cmplx(  -0.004625240588609,   0.000000000000000)
  u_pot%p4(5,3,2,1) = cmplx(  -0.084641548927263,  -0.087912087912088)
  u_pot%p4(5,3,2,2) = cmplx(   0.000000000000000,  -0.052499747452378)
  u_pot%p4(5,3,2,3) = cmplx(   0.000000000000000,   0.055170131351278)
  u_pot%p4(5,3,2,4) = cmplx(  -0.228402304294797,   0.148491709736631)
  u_pot%p4(5,3,2,5) = cmplx(  -0.086821908250480,  -0.087912087912088)
  u_pot%p4(5,3,3,1) = cmplx(   0.091182626896913,   0.087912087912088)
  u_pot%p4(5,3,3,2) = cmplx(   0.000000000000000,  -0.055170131351278)
  u_pot%p4(5,3,3,3) = cmplx(   0.000000000000000,   0.052499747452378)
  u_pot%p4(5,3,3,4) = cmplx(  -0.076134101431599,   0.148491709736631)
  u_pot%p4(5,3,3,5) = cmplx(   0.095543345543346,   0.087912087912088)
  u_pot%p4(5,3,4,1) = cmplx(  -0.186489700532716,   0.000000000000000)
  u_pot%p4(5,3,4,2) = cmplx(  -0.005664739689851,   0.000000000000000)
  u_pot%p4(5,3,4,3) = cmplx(  -0.001888246563284,   0.000000000000000)
  u_pot%p4(5,3,4,5) = cmplx(  -0.186489700532716,   0.000000000000000)
  u_pot%p4(5,3,5,2) = cmplx(  -0.043956043956044,  -0.091182626896913)
  u_pot%p4(5,3,5,3) = cmplx(  -0.219780219780220,   0.091182626896913)
  u_pot%p4(5,3,5,4) = cmplx(   0.006166987451479,   0.000000000000000)
  u_pot%p4(5,4,1,2) = cmplx(  -0.181864459944107,   0.000000000000000)
  u_pot%p4(5,4,1,3) = cmplx(  -0.181864459944107,   0.000000000000000)
  u_pot%p4(5,4,2,2) = cmplx(   0.000000000000000,   0.074245854868315)
  u_pot%p4(5,4,2,3) = cmplx(   0.000000000000000,  -0.074245854868315)
  u_pot%p4(5,4,2,4) = cmplx(   0.104999494904756,   0.000000000000000)
  u_pot%p4(5,4,3,2) = cmplx(   0.000000000000000,   0.074245854868315)
  u_pot%p4(5,4,3,3) = cmplx(   0.000000000000000,  -0.074245854868315)
  u_pot%p4(5,4,3,4) = cmplx(  -0.104999494904756,   0.000000000000000)
  u_pot%p4(5,4,4,1) = cmplx(   0.257195185766614,   0.000000000000000)
  u_pot%p4(5,4,4,5) = cmplx(   0.272457701029130,   0.000000000000000)
  u_pot%p4(5,4,5,2) = cmplx(   0.192656687984194,   0.000000000000000)
  u_pot%p4(5,4,5,3) = cmplx(   0.192656687984194,   0.000000000000000)
  u_pot%p4(5,5,1,1) = cmplx(   0.336821908250480,   0.000000000000000)
  u_pot%p4(5,5,2,2) = cmplx(  -0.263300191871620,   0.000000000000000)
  u_pot%p4(5,5,2,3) = cmplx(   0.085731728588871,   0.000000000000000)
  u_pot%p4(5,5,3,2) = cmplx(   0.085731728588871,   0.000000000000000)
  u_pot%p4(5,5,3,3) = cmplx(  -0.263300191871620,   0.000000000000000)
  u_pot%p4(5,5,4,4) = cmplx(  -0.177568463282749,   0.000000000000000)
  u_pot%p4(5,5,5,5) = cmplx(   0.367346938775510,   0.000000000000000)

end subroutine initialize_p_matrices
!===============================================================
