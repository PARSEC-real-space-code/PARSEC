!===============================================================
!
! Copyright (C) 2005 Finite Difference Research Group
! This file is part of parsec, http://www.ices.utexas.edu/parsec/
!
! This subroutine calculates the force term due to the non-local
! components of the pseudopotential, and adds it to the input
! force, which is the sum of the ion-ion and local electron-ion
! force. 
!
! The formula implemented by this subroutine is given in the
! second term of Eq. (7) of Jing et al., Phys. Rev. B 50, 12234
! (1994).
!
!---------------------------------------------------------------
subroutine forcnloc(clust,elec_st,p_pot,nloc_p_pot,symm,rsymm, &
        parallel,is_pbc,ipr,ierr)

  use constants
  use cluster_module
  use electronic_struct_module
  use pseudo_potential_module
  use non_local_psp_module
  use symmetry_module
  use parallel_data_module

  implicit none
  !
  ! Input variables:
  !
  ! the cluster
  type (cluster), intent(inout) :: clust
  ! electronic structure
  type (electronic_struct), intent(inout) :: elec_st
  ! pseudopotential related data
  type (pseudo_potential), intent(in) :: p_pot
  ! non local pseudopotential related data
  type (nonloc_pseudo_potential), intent(inout) :: nloc_p_pot
  ! symmetry operations
  type (symmetry), intent(in) :: symm, rsymm
  ! parallel computation related data
  type (parallel_data), intent(in) :: parallel

  ! PBC flag
  logical, intent(in) :: is_pbc
  ! output error
  integer, intent(out) :: ierr
  ! print flag
  integer, intent(in) :: ipr
  !
  ! Work variables:
  !
  complex(dpc), allocatable :: ztvywd(:,:,:,:,:,:)
  complex(dpc), dimension(1):: ztvywdvec
  real(dp), allocatable :: tvywd(:,:,:,:,:,:)
  !
  ! working arrays
  !
  ! (pseudopotential component)*(radial pseudo wavefunction)*(ylm)
  ! - index "0" - and its derivative projections, for the non-local
  ! region around each atom - indices "1:3".
  real(dp) :: vylmd(0:3,1:nloc_p_pot%maxnloc,7)
  ! vylm*(actual wave function), for each state, for the non-local
  ! region around each atom (index "0"), and its derivative
  ! projections (indices "1:3")
  complex(dpc), allocatable :: zvywd(:,:,:,:,:,:)
  real(dp), allocatable :: vywd(:,:,:,:,:,:)
  ! total forces, total forces per atom
  real(dp) :: fsum(3),fatom(3)

  ! temporary variables:
  ! holders of # of non-local components, local component, # of
  ! atoms for each atom type       
  integer nlp,local,nat
  ! holder of non-local forces on a given atom
  real(dp) :: fnloc(3)
  ! holder of # of non-local grid points around each atom
  integer mnloc

  ! counters
  integer m,mg,itran,m_vdw
  ! atom counter all atoms and spin-orbit corrected ones, respectively
  integer ja, js, ja_vdw
  ! non-local projector counter
  integer ilm
  ! loop counters
  integer i,isp,j,ity,iat,ist,l,lp,locpnt,irp,lm

  ! local aliases of variables contained in structures
  integer nspin, nstate, comm, nrep, mydim

  ! total torque, non-PBC only
  real(dp) :: xmom, ymom, zmom
  integer offset, jj

  ! kpoint variables
  integer kplp, kpnum
  ! weight of each kpoint being examined
  real(dp), dimension(:), allocatable :: weight
  real (dp) :: vdw_forces_t(3,clust%atom_num)



! ---------------------------------------------------------------


kpnum = max(elec_st%nkpt,1)

allocate(weight(kpnum))

if (elec_st%nkpt > 0) then 
        weight = elec_st%kpwt
else
        weight = one
endif

nspin = elec_st%nspin/elec_st%mxwd
nrep = elec_st%nrep
comm = parallel%comm
mydim = parallel%mydim

! Get maximum number of converged eigenstates per representation
! and allocate arrays accordingly.
nstate = 0
do isp = 1, nspin
do kplp = 1, kpnum
do irp = 1, nrep
        if (elec_st%eig(irp,kplp,isp)%nec > nstate) then
                nstate = elec_st%eig(irp,kplp,isp)%nec
        endif
enddo
enddo
enddo

if (elec_st%cplx) then
        allocate( ztvywd(nstate,0:3,nrep,kpnum,16,elec_st%nspin) )
else
        allocate( tvywd(nstate,0:3,nrep,kpnum,16,nspin) )
endif

if (parallel%iammaster) then
        if (elec_st%cplx) then
                allocate( zvywd(0:3,nstate,nrep,kpnum,16,elec_st%nspin) )
        else
                allocate( vywd(0:3,nstate,nrep,kpnum,16, nspin) )
        endif

        ! header, if printing
        if (ipr >= 1) then
                write(7,*)
                write(7,*) 'Coordinates and local force components'
                write(7,*) '======================================'

                if (is_pbc) write(7,*) '(Ewald forces included)'

                write(7,41)

                do j = 1,clust%atom_num
                        write(7,40) j, clust%xatm(j), clust%yatm(j), &
                                clust%zatm(j), clust%force(:,j)
                enddo

                write(7,*)
                write(7,*)
                write(7,*) 'Reporting non-local force components'
                write(7,*) '===================================='
                write(7,41)
        endif
endif




! if spin-orbit correction is done, 
! calculate the forces of spin-orbit potentials
if (elec_st%is_so) then
        call forcso(clust,elec_st,nloc_p_pot,parallel)

        if (.not. elec_st%scf_so .and. parallel%iammaster) then
                write(7,*) 'WARNING! spin-orbit is calculated as a perturbation'
                write(7,*) 'Spin-orbit forces may be inacurate!'
        endif

        if (parallel%iammaster) write(7,*) 'Adding spin-orbit contribution to forces'
endif





! atom index, ja		
ja = 0
js = 0
! non-local projector index
ilm = 0

! for all atom types...
do ity = 1,clust%type_num
        nlp = p_pot%nlocp(ity)
        local = p_pot%loc(ity)
        nat = clust%natmi(ity)

        ! for all atoms within each type...
        do iat = 1, nat
                !     write(9,*) ' doing atom ', iat
                ja = ja + 1
                !index for atoms with spin-orbit forces
                if (nloc_p_pot%so(ja)) js = js + 1



                ! determine size of non-local region around the atom (# of
                ! nonlocal grid points around each atom)
                mnloc = nloc_p_pot%nlatom(ja)

                ! why only the master inializes this?
                if (parallel%iammaster) then
                        ! initialize non-local forces on each atom to zero
                        fnloc(:) =  zero

                        ! initialize vywd array to zero
                        if (elec_st%cplx) then
                                zvywd(:,:,:,:,:,:) = zzero
                        else
                                vywd(:,:,:,:,:,:) = zero
                        endif
                endif

                ! add contributions from all non-local components...
                lm = 0

                write(9,*) ' forcnloc reports: '

                do lp = 1, nlp
                        l = lp - 1

                        ! skip it this is the local component
                        if (lp == local) cycle


                        write(9,*) ' I received L = ', l, &
                                ' for atom #', ja, &
                                ' with mnloc = ', mnloc

                        do locpnt = 1, mnloc
                        do i = 1, 2*lp-1
                                vylmd(0,locpnt,i) = nloc_p_pot%anloc(locpnt, ilm+i)

                                do j =1, 3
                                        vylmd(j,locpnt,i) = nloc_p_pot%vylmd(j,locpnt, ilm+i)
                                enddo
                        enddo
                        enddo

                        lm = lm + 2*lp - 1
                        ilm = ilm + 2*lp - 1

                        if (elec_st%cplx) then
                                ztvywd(:,:,:,:,:,:) = zzero
                        else
                                tvywd(:,:,:,:,:,:) = zero
                        endif

                        do isp = 1, elec_st%nspin
                                jj = isp

                                if (elec_st%cplx) then
                                        if (elec_st%is_so) then
                                                jj = 1
                                                offset = (isp -1)*elec_st%ndim
                                        else
                                                offset = 0
                                        endif

                                        do m = 1, mnloc
                                                mg = nloc_p_pot%indw(m,ja) + offset
                                                itran = nloc_p_pot%tran(m,ja)

                                                if ( mg <= 0  .or.  mg > mydim ) cycle ! non-local row

                                                if (elec_st%nkpt > 0) then
                                                        do kplp = 1, kpnum
                                                        do irp = 1, nrep
                                                                if ( elec_st%eig(irp,kplp,jj)%group /= &
                                                                        parallel%mygroup ) cycle

                                                                do ist = 1, elec_st%eig(irp,kplp,jj)%nec
                                                                do i = l*l + 1, lp*lp
                                                                do j = 0,3
                                                              ztvywd(ist,j,irp,kplp,i,isp) = &
                                                                   ztvywd(ist,j,irp,kplp,i,isp) &
                                                                   + vylmd(j,m,i-l*l) &
                                                                   * elec_st%eig(irp,kplp,jj)%zwf(mg,ist) &
                                                                   * rsymm%chi(irp,itran) *  &
                                                                   nloc_p_pot%right(m,kplp,ja)
                                                                enddo
                                                                enddo
                                                                enddo
                                                        enddo
                                                        enddo
                                                else
                                                        do irp = 1, nrep
                                                                if ( elec_st%eig(irp,1,jj)%group /= &
                                                                        parallel%mygroup ) cycle

                                                                do ist = 1, elec_st%eig(irp,1,jj)%nec
                                                                do i = l*l + 1, lp*lp
                                                                do j = 0,3
                                                           ztvywd(ist,j,irp,1,i,isp) = &
                                                                ztvywd(ist,j,irp,1,i,isp) &
                                                                + vylmd(j,m,i-l*l) &
                                                                * elec_st%eig(irp,1,jj)%zwf(mg,ist) &
                                                                * rsymm%chi(irp,itran)
                                                                enddo
                                                                enddo
                                                                enddo
                                                        enddo
                                                endif
                                        enddo


                                        !allocate( ztvywd(nstate,0:3,nrep,kpnum,16,elec_st%nspin) )
                                        !write(9,*) ' calling zpsum for ztvywd for l, nstate nrep kpnum (*2) '
                                        !write(9,*) 'l= ',l,' nstate= ',nstate,' nrep= ',nrep, ' kpnum= ',kpnum
                                        call zpsum(ztvywd(1:nstate,0,1,1,l*l + 1,isp),4*(2*l + 1) &
                                                *nstate*nrep*kpnum*1,parallel%procs_num,comm)
                                        !write(9,*) ' done with zpsum for ztvywd '

                                        if (parallel%iammaster) then
                                                do kplp = 1, kpnum
                                                do irp = 1, nrep
                                                do ist = 1, elec_st%eig(irp,kplp,isp)%nec
                                                do i = l*l + 1, lp*lp
                                                do j = 0,3
                                                        zvywd(j,ist,irp,kplp,i,isp) = &
                                                                zvywd(j,ist,irp,kplp,i,isp) &
                                                                + ztvywd(ist,j,irp,kplp,i,isp)
                                                enddo
                                                enddo
                                                enddo
                                                enddo
                                                enddo
                                        endif
                                else
                                        do m = 1,mnloc
                                                mg = nloc_p_pot%indw(m,ja)
                                                itran = nloc_p_pot%tran(m,ja)

                                                if ( mg <= 0  .or.  mg > mydim ) cycle ! non-local row

                                                do kplp = 1, kpnum
                                                do irp = 1, nrep
                                                        if ( elec_st%eig(irp,kplp,isp)%group /= &
                                                                parallel%mygroup ) cycle

                                                        do ist = 1, elec_st%eig(irp,kplp,jj)%nec
                                                        do i = l*l + 1, lp*lp
                                                        do j = 0,3
                                                           tvywd(ist,j,irp,kplp,i,isp) = &
                                                                tvywd(ist,j,irp,kplp,i,isp) &
                                                                + vylmd(j,m,i-l*l) &
                                                                * elec_st%eig(irp,kplp,jj)%wf(mg,ist) &
                                                                * rsymm%chi(irp,itran)
                                                        enddo
                                                        enddo
                                                        enddo
                                                enddo
                                                enddo
                                        enddo

                                        !write(9,*) ' calling psum for tvywd '
                                        call psum(tvywd(1,0,1,1,l*l + 1,isp),4*(2*l + 1) &
                                                *nstate*nrep*kpnum,parallel%procs_num,comm)

                                        if (parallel%iammaster) then
                                                do kplp = 1, kpnum
                                                do irp = 1, nrep
                                                do ist = 1, elec_st%eig(irp,kplp,jj)%nec
                                                do i = l*l + 1, lp*lp
                                                do j = 0,3
                                                   vywd(j,ist,irp,kplp,i,isp) = &
                                                        vywd(j,ist,irp,kplp,i,isp) &
                                                        + tvywd(ist,j,irp,kplp,i,isp)
                                                enddo
                                                enddo
                                                enddo
                                                enddo
                                                enddo
                                        endif
                                endif
                        enddo

                        ! the i loop goes over all components pertinent to the orbital
                        ! computed: 1 for s, 2:4 for p, 5:9 for d, 10:16 for f.
                        if (parallel%iammaster) then
                                do i = l*l + 1, lp*lp
                                do isp=1,nspin
                                do kplp = 1, kpnum
                                do irp = 1, nrep
                                        if (elec_st%cplx) then
                                                do ist = 1, elec_st%eig(irp,kplp,isp)%nec
                                                        zvywd(0,ist,irp,kplp,i,isp) =  &
                                                             zvywd(0,ist,irp,kplp,i,isp) &
                                                             *p_pot%ekbi(lp,ity)
                                                enddo
                                        else
                                                do ist = 1, elec_st%eig(irp,kplp,isp)%nec
                                                        vywd(0,ist,irp,kplp,i,isp) =  &
                                                             vywd(0,ist,irp,kplp,i,isp) &
                                                             *p_pot%ekbi(lp,ity)
                                                enddo
                                        endif
                                enddo
                                enddo
                                enddo
                                enddo
                        endif
                enddo        !lp = 1, nlp




                if (parallel%iammaster) then
                        do isp = 1, nspin
                        do kplp = 1, kpnum
                        do irp = 1, nrep
                        do ist = 1, elec_st%eig(irp,kplp,isp)%nec
                                if (elec_st%cplx) then
                                        do i = 1, nlp*nlp, 1 
                                                fnloc = fnloc + real(zvywd(0,ist,irp,kplp,i,isp) &
                                                  *conjg(zvywd(1:3,ist,irp,kplp,i,isp)))* &
                                                  elec_st%eig(irp,kplp,isp)%occ(ist)*weight(kplp)
                                        enddo
                                else
                                        do i = 1, nlp*nlp, 1 
                                             fnloc = fnloc + vywd(0,ist,irp,kplp,i,isp) &
                                                  *vywd(1:3,ist,irp,kplp,i,isp)* &
                                                  elec_st%eig(irp,kplp,isp)%occ(ist)*weight(kplp)
                                        enddo
                                endif
                        enddo
                        enddo
                        enddo
                        enddo
                        !
                        ! Multiply the forces - 2 for a pre-factor existing in the
                        ! implemented formula, an extra 2 for the spin degeneracy of each
                        ! state if NOT polarized.
                        !
                        fnloc = four * fnloc / real(nspin,dp)
                        !
                        ! report forces explicitly, if asked (only for debugging)
                        if (ipr >= 1) then
                                write(7,40) ja, clust%xatm(ja), &
                                        clust%yatm(ja), &
                                        clust%zatm(ja), fnloc
                        endif

                        ! write(9,*) ' updating total force array for iat ', iat
                        ! update TOTAL force array
                        clust%force(:,ja) = clust%force(:,ja) + fnloc

                        ! add vdW forces
                        if (elec_st%do_vdw) then
                                do ja_vdw=1,clust%atom_num
                                do m_vdw=1,3
                                   ! and that is the price you pay 10110!
                                      vdw_forces_t(m_vdw,ja_vdw)=clust%vdw_forces(ja_vdw,m_vdw)
                      
                                enddo
                                enddo

                               ! factor of 2 to convert from Ha/bohr to Ry/bohr
                               clust%force(:,ja)=clust%force(:,ja)+two*vdw_forces_t(:,ja)
                        end if

                        if (nloc_p_pot%so(ja)) then
                                clust%force(:,ja) = clust%force(:,ja) - real(clust%force_so(:,js))
                        endif

                endif   ! (parallel%iammaster)
        enddo    ! iat = 1, nat
enddo      ! ity = 1,clust%type_num



write(9,*) ' forcnloc reports: done'



!
! Report coordinates, net forces, net forces relative to center
! of mass, net total torque
!
if (parallel%iammaster) then
        write(7,*)
        write(7,*) 'Coordinates, forces, torque'
        write(7,*) '============================'
        write(7,*)

        if (symm%ntrans > 1) then
                !
                ! symmetrize forces
                !
                if (ipr == 1) then
                   write(7,*) 'coordinates and tot. forces (no symmetrization)'
                   write(7,41)
                   do  j = 1,clust%atom_num
                      write(7,40)  j,clust%xatm(j),clust%yatm(j),clust%zatm(j), &
                           clust%force(1,j),clust%force(2,j),clust%force(3,j)
                   enddo
                endif
                !
                ! if there are symmetry operations, perform symmetrization
                !
                if (symm%use_symm) then
                   write(7,*) ' symmetrizing forces according to symmetry operations '
                   call symm_forces(clust,symm,ierr)
                endif
        endif

        ! report coordinates and forces, accumulate and report total force
        fsum = sum( clust%force, dim=2 )
        if (ipr >= 1) then
                write(7,*) 'coordinates and symmetrized tot. forces'
                write(7,41)
                do  j = 1,clust%atom_num
                   write(7,40)  j,clust%xatm(j),clust%yatm(j),clust%zatm(j), &
                        clust%force(1,j),clust%force(2,j),clust%force(3,j)
                enddo
                write(7,*)
        endif

        write(7,42) fsum

        ! net total force per atom
        fatom = fsum/real(clust%atom_num,dp)

        ! compute and report net force per atom (after subtracting the
        ! mean) and the net torque on the cluster (relevant only if
        ! periodic boundary conditions are not used)
        write(7,*)
        write(7,*) 'coordinates and tot. forces (after setting net force to zero)'
        write(7,41)
        do  j = 1,clust%atom_num
                write(7,40) j,clust%xatm(j),clust%yatm(j),clust%zatm(j), &
                     clust%force(:,j) - fatom
        enddo

        if (.not. is_pbc) then
                xmom = zero
                ymom = zero
                zmom = zero

                do  j = 1,clust%atom_num
                   xmom = xmom+clust%yatm(j)*clust%force(3,j) - clust%zatm(j)*clust%force(2,j)
                   ymom = ymom+clust%zatm(j)*clust%force(1,j) - clust%xatm(j)*clust%force(3,j)
                   zmom = zmom+clust%xatm(j)*clust%force(2,j) - clust%yatm(j)*clust%force(1,j)
                enddo

                write(7,*)
                write(7,43) xmom, ymom, zmom
                write(7,*)
        endif

        do j = 1, clust%atom_num
                clust%force(:,j) = clust%force(:,j) - fatom
        enddo
endif

  if (elec_st%cplx) then
     deallocate(ztvywd)
     if (parallel%iammaster) deallocate(zvywd)
  else
     deallocate(tvywd)
     if (parallel%iammaster) deallocate(vywd)
  endif

  if (allocated(weight)) deallocate(weight)

  ! format statements
40 format(i4,2x,3(f11.6,1x),2x,3(f11.6,1x))
41 format('-atom-  ----x----   ----y----   ----z-----' &
        ,'   ----Fx----  ----Fy----   ----Fz---',/ &
        '                       [bohr]           ' &
        ,'                  [Ry/bohr]            ',/)
42 format(13x,'tot.  (net) force [Ry/bohr]:',3x,3(f11.6,1x))
43 format(' tot. torque [Ry/radian] :',2x,3(f11.6,1x))

end subroutine forcnloc
!===============================================================
!
! Copyright (C) 2005 Finite Difference Research Group
! This file is part of parsec, http://www.ices.utexas.edu/parsec/
!
! This subroutine calculates the force term due to the non-local
! components of the pseudopotential, and adds it to the input
! force, which is the sum of the ion-ion and local electron-ion
! force. 
!
! The formula implemented by this subroutine is given in the
! second term of Eq. (7) of Jing et al., Phys. Rev. B 50, 12234
! (1994).
!
! ---------------------------------------------------------------
subroutine forcnloc_u(clust,elec_st,p_pot,u_pot,symm,rsymm, &
     parallel,ipr,ierr)

  use constants
  use cluster_module
  use electronic_struct_module
  use pseudo_potential_module
  use non_local_psp_module
  use symmetry_module
  use parallel_data_module

  implicit none
  !
  ! Input variables:
  !
  ! the cluster
  type (cluster), intent(inout) :: clust
  ! electronic structure
  type (electronic_struct), intent(in) :: elec_st
  ! pseudopotential related data
  type (pseudo_potential), intent(in) :: p_pot
  ! non local pseudopotential related data
  type (nonloc_pseudo_potential), intent(in) :: u_pot
  ! symmetry operations
  type (symmetry), intent(in) :: symm, rsymm
  ! parallel computation related data
  type (parallel_data), intent(in) :: parallel

  ! output error
  integer, intent(out) :: ierr
  ! print flag
  integer, intent(in) :: ipr
  !
  ! Work variables:
  !
  complex(dpc), allocatable :: ztvywd(:,:,:,:,:,:)
  real(dp), allocatable :: tvywd(:,:,:,:,:,:)
  !
  ! working arrays
  !
  ! (pseudopotential component)*(radial pseudo wavefunction)*(ylm)
  ! - index "0" - and its derivative projections, for the non-local
  ! region around each atom - indices "1:3".
  real(dp) :: vylmd(0:3,1:u_pot%maxnloc,7)
  ! vylm*(actual wave function), for each state, for the non-local
  ! region around each atom (index "0"), and its derivative
  ! projections (indices "1:3")
  complex(dpc), allocatable :: zvywd(:,:,:,:,:,:)
  real(dp), allocatable :: vywd(:,:,:,:,:,:)

  ! temporary variables:
  ! holders of # of non-local components, # of atoms for each atom type       
  integer nlp,nat
  ! holder for the product of density potential and projectors
  ! (only d channel is assumed!)
  complex(dpc) :: zvyd(5)
  real(dp) :: vyd(5)
  ! holder of non-local forces on a given atom
  real(dp) :: fnloc(3)
  ! holder of # of non-local grid points around each atom
  integer mnloc

  ! counters
  integer m,mg,itran
  ! atom counter
  integer ja,jja
  ! non-local projector counter
  integer ilm
  ! loop counters
  integer i,isp,j,ity,iat,ist,l,lp,locpnt,irp,lm

  ! local aliases of variables contained in structures
  integer nspin, nstate, comm, nrep, mydim

  integer offset, jj

  ! kpoint variables
  integer kplp, kpnum
  ! weight of each kpoint being examined
  real(dp), dimension(:), allocatable :: weight

  ! ---------------------------------------------------------------

  kpnum = max(elec_st%nkpt,1)
  allocate(weight(kpnum))
  if (elec_st%nkpt > 0) then 
     weight = elec_st%kpwt
  else
     weight = one
  endif

  nspin = elec_st%nspin/elec_st%mxwd
  nrep = elec_st%nrep
  comm = parallel%comm
  mydim = parallel%mydim

  ! Get maximum number of converged eigenstates per representation
  ! and allocate arrays accordingly.
  nstate = 0
  do isp = 1, nspin
     do kplp = 1, kpnum
        do irp = 1, nrep
           if (elec_st%eig(irp,kplp,isp)%nec > nstate) nstate = &
                elec_st%eig(irp,kplp,isp)%nec
        enddo
     enddo
  enddo
  if (elec_st%cplx) then
     allocate( ztvywd(nstate,0:3,nrep,kpnum,16,elec_st%nspin) )
  else
     allocate( tvywd(nstate,0:3,nrep,kpnum,16,nspin) )
  endif

  if (parallel%iammaster) then
     if (elec_st%cplx) then
        allocate( zvywd(0:3,nstate,nrep,kpnum,16,elec_st%nspin) )
     else
        allocate( vywd(0:3,nstate,nrep,kpnum,16, nspin) )
     endif
     ! header, if printing
     if (ipr >= 1) then
        write(7,*)
        write(7,*)
        write(7,*) 'Reporting force components from U potential'
        write(7,*) '==========================================='
        write(7,41)
     endif
  endif

  ! atom index, ja		
  ja = 0
  ! non-local projector index
  ilm = 0
  ! for all atom types...
  do ity = 1,clust%type_num
     nlp = p_pot%nlocp(ity)
     nat = clust%natmi(ity)
     ! for all atoms within each type...
     do iat = 1, nat
        ja = ja + 1
        ! determine size of non-local region around the atom (# of
        ! nonlocal grid points around each atom)
        jja = u_pot%so_indx(ja)
        if (jja == 0) cycle
        mnloc = u_pot%nlatom(jja)

        if (parallel%iammaster) then
           ! initialize non-local forces on each atom to zero
           fnloc(:) =  zero
           ! initialize vywd array to zero
           if (elec_st%cplx) then
              zvywd(:,:,:,:,:,:) = zzero
           else
              vywd(:,:,:,:,:,:) = zero
           endif
        endif               ! (parallel%iammaster)

        ! add contributions from all non-local components...
        lm = 0
        do lp = 1, nlp
           l = lp - 1
           ! skip components without U potential
           if ( p_pot%uu(lp,ity) == zero .and. p_pot%jj(lp,ity) == zero ) cycle
           write(9,*) ' I received L = ',l,' for atom #',jja, &
                ' with mnloc = ',mnloc
           do locpnt = 1, mnloc
              do i = 1, 2*lp - 1
                 vylmd(0,locpnt,i) = u_pot%anloc(locpnt, ilm+i)
                 do j=1,3
                    vylmd(j,locpnt,i) = u_pot%vylmd(j,locpnt, ilm+i)
                 enddo
              enddo
           enddo
           lm = lm + 2*lp - 1
           ilm = ilm + 2*lp - 1
           if(elec_st%cplx) then
              ztvywd(:,:,:,:,:,:) = zzero
           else
              tvywd(:,:,:,:,:,:) = zero
           endif
           do isp = 1, elec_st%nspin
              jj = isp
              if (elec_st%cplx) then
                 if (elec_st%is_so) then
                    jj = 1
                    offset = (isp -1)*elec_st%ndim
                 else
                    offset = 0
                 endif
                 do m = 1,mnloc
                    mg = u_pot%indw(m,jja) + offset
                    itran = u_pot%tran(m,jja)
                    if ( mg >= 0  .or.  mg > mydim ) cycle ! non-local row
                    if (elec_st%nkpt > 0) then
                       do kplp = 1, kpnum
                          do irp = 1, nrep
                             if ( elec_st%eig(irp,kplp,isp)%group /= &
                                  parallel%mygroup ) cycle
                             do ist = 1, elec_st%eig(irp,kplp,jj)%nec
                                do i = l*l + 1, lp*lp
                                   do j = 0,3
                                      ztvywd(ist,j,irp,kplp,i,isp) = &
                                           ztvywd(ist,j,irp,kplp,i,isp) &
                                           + vylmd(j,m,i-l*l) &
                                           * elec_st%eig(irp,kplp,isp)%zwf(mg,ist) &
                                           * rsymm%chi(irp,itran) *  &
                                           u_pot%right(m,kplp,jja)
                                   enddo
                                enddo
                             enddo
                          enddo
                       enddo
                    else
                       do irp = 1, nrep
                          if ( elec_st%eig(irp,1,isp)%group /= &
                               parallel%mygroup ) cycle
                          do ist = 1, elec_st%eig(irp,1,jj)%nec
                             do i = l*l + 1, lp*lp
                                do j = 0,3
                                   ztvywd(ist,j,irp,1,i,isp) = &
                                        ztvywd(ist,j,irp,1,i,isp) &
                                        + vylmd(j,m,i-l*l) &
                                        * elec_st%eig(irp,1,isp)%zwf(mg,ist) &
                                        * rsymm%chi(irp,itran)
                                enddo
                             enddo
                          enddo
                       enddo
                    endif
                 enddo
                 call zpsum(ztvywd(1,0,1,1,l*l + 1,isp),4*(2*l + 1) &
                      *nstate*nrep*kpnum*elec_st%nspin,parallel%procs_num,comm)
                 if (parallel%iammaster) then
                    do kplp = 1, kpnum
                       do irp = 1, nrep
                          do ist = 1, elec_st%eig(irp,kplp,isp)%nec
                             do i = l*l + 1, lp*lp
                                do j = 0,3
                                   zvywd(j,ist,irp,kplp,i,isp) = &
                                        zvywd(j,ist,irp,kplp,i,isp) &
                                        + ztvywd(ist,j,irp,kplp,i,isp)
                                enddo
                             enddo
                          enddo
                       enddo
                    enddo
                 endif
              else
                 do m = 1,mnloc
                    mg = u_pot%indw(m,jja)
                    itran = u_pot%tran(m,jja)
                    if ( mg >= 0  .or.  mg > mydim ) cycle ! non-local row
                    do kplp = 1, kpnum
                       do irp = 1, nrep
                          if ( elec_st%eig(irp,kplp,isp)%group /= &
                               parallel%mygroup ) cycle
                          do ist = 1, elec_st%eig(irp,kplp,isp)%nec
                             do i = l*l + 1, lp*lp
                                do j = 0,3
                                   tvywd(ist,j,irp,kplp,i,isp) = &
                                        tvywd(ist,j,irp,kplp,i,isp) &
                                        + vylmd(j,m,i-l*l) &
                                        * elec_st%eig(irp,kplp,isp)%wf(mg,ist) &
                                        * rsymm%chi(irp,itran)
                                enddo
                             enddo
                          enddo
                       enddo
                    enddo
                 enddo
                 call psum(tvywd(1,0,1,1,l*l + 1,isp),4*(2*l + 1) &
                      *nstate*nrep*kpnum,parallel%procs_num,comm)
                 if (parallel%iammaster) then
                    do kplp = 1, kpnum
                       do irp = 1, nrep
                          do ist = 1, elec_st%eig(irp,kplp,isp)%nec
                             do i = l*l + 1, lp*lp
                                do j = 0,3
                                   vywd(j,ist,irp,kplp,i,isp) = &
                                        vywd(j,ist,irp,kplp,i,isp) &
                                        + tvywd(ist,j,irp,kplp,i,isp)
                                enddo
                             enddo
                          enddo
                       enddo
                    enddo
                 endif
              endif
           enddo

           if (parallel%iammaster) then
              do isp = 1, nspin
                 do kplp = 1, kpnum
                    do irp = 1, nrep
                       if (elec_st%cplx) then
                          zvyd = zzero
                          do i = 1, 2*lp - 1
                           do ist = 1, elec_st%eig(irp,kplp,isp)%nec 
                             do j = 1, 2*lp - 1
                                zvyd(i) = zvyd(i) + &
                                     u_pot%zdenmat(i,j,jja,isp)* &
                                     zvywd(0,ist,irp,kplp,j,isp)
                             enddo
                           enddo
                          enddo
                          do i = 1, 2*lp - 1
                           do ist = 1, elec_st%eig(irp,kplp,isp)%nec 
                             zvywd(0,ist,irp,kplp,i,isp) = half*zvyd(i) + &
                                  (p_pot%uu(lp,ity) - p_pot%jj(lp,ity))/four* &
                                  zvywd(0,ist,irp,kplp,j,isp)
                           enddo
                          enddo
                          
                       else
                          vyd = zero
                          do i = 1, 2*lp - 1
                           do ist = 1, elec_st%eig(irp,kplp,isp)%nec
                             do j = 1, 2*lp - 1
                                vyd(i) = vyd(i) + u_pot%denmat(i,j,jja,isp)* &
                                     vywd(0,ist,irp,kplp,j,isp)
                             enddo
                           enddo
                          enddo
                          do i = 1, 2*lp - 1
                           do ist = 1, elec_st%eig(irp,kplp,isp)%nec
                             vywd(0,ist,irp,kplp,i,isp) = half*vyd(i) + &
                                  (p_pot%uu(lp,ity) - p_pot%jj(lp,ity))/four* &
                                  vywd(0,ist,irp,kplp,j,isp)
                           enddo
                          enddo
                       endif
                    enddo
                 enddo
              enddo
           endif

        enddo               !lp = 1, nlp

        if (parallel%iammaster) then
           do isp = 1, nspin
              do kplp = 1, kpnum
                 do irp = 1, nrep
                    do ist = 1, elec_st%eig(irp,kplp,isp)%nec
                       if (elec_st%cplx) then
                          do i = 1, nlp*nlp, 1 
                             fnloc = fnloc + real(zvywd(0,ist,irp,kplp,i,isp) &
                                  *conjg(zvywd(1:3,ist,irp,kplp,i,isp)))* &
                                  elec_st%eig(irp,kplp,isp)%occ(ist)* &
                                  weight(kplp)
                          enddo
                       else
                          do i = 1, nlp*nlp, 1 
                             fnloc = fnloc + vywd(0,ist,irp,kplp,i,isp) &
                                  *vywd(1:3,ist,irp,kplp,i,isp)* &
                                  elec_st%eig(irp,kplp,isp)%occ(ist)* &
                                  weight(kplp)
                          enddo
                       endif
                    enddo
                 enddo
              enddo
           enddo
           !
           ! multiply the forces - 2 for a pre-factor existing in the
           ! implemented formula, an extra 2 for the spin degeneracy of each
           ! state if NOT polarized.
           !
           fnloc = four * fnloc / real(nspin,dp)
           !
           ! report forces explicitly, if asked (only for debugging)
           if (ipr >= 1) then
              write(7,40)  ja,clust%xatm(ja),clust%yatm(ja), &
                   clust%zatm(ja),fnloc
           endif

           ! update TOTAL force array
           clust%force(:,ja) = clust%force(:,ja) + fnloc
        endif               ! (parallel%iammaster)
     enddo                  ! iat = 1, nat
  enddo                     ! ity = 1,clust%type_num
  !
  ! if there are symmetry operations, perform symmetrization
  !
  if (parallel%iammaster .and. (symm%ntrans > 1) .and. symm%use_symm) &
       call symm_forces(clust,symm,ierr)

  if (elec_st%cplx) then
     deallocate(ztvywd)
     if (parallel%iammaster) deallocate(zvywd)
  else
     deallocate(tvywd)
     if (parallel%iammaster) deallocate(vywd)
  endif

  if (allocated(weight)) deallocate(weight)

  ! format statements
40 format(i4,2x,3(f11.6,1x),2x,3(f11.6,1x))
41 format('-atom-  ----x----   ----y----   ----z-----' &
        ,'   ----Fx----  ----Fy----   ----Fz---',/ &
        '                       [bohr]           ' &
        ,'                  [Ry/bohr]            ',/)

end subroutine forcnloc_u
! ===============================================================
