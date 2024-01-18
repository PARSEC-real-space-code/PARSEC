!===============================================================
!
! Copyright (C) 2005 Finite Difference Research Group
! This file is part of parsec, http://www.ices.utexas.edu/parsec/
!
!---------------------------------------------------------------
subroutine pls(elec_st,solver,nloc_p_pot,parallel)
  use constants
  use electronic_struct_module
  use non_local_psp_module
  use parallel_data_module
  use eigen_solver_module
#ifdef MPI
  use mpi
#endif
  implicit none
  !
  ! Input/Output variables:
  !
  ! electronic structure
  type (electronic_struct), intent(inout) :: elec_st
  ! solver related data
  type (eigen_solver), intent(in) :: solver
  ! non local pseudopotential related data
  type (nonloc_pseudo_potential), intent(in) :: nloc_p_pot 
  ! parallel computation related data
  type (parallel_data), intent(in) :: parallel
  !
  ! Work variables:
  !
  ! distributed spinor (two-component wavefunctions)
  complex(dpc),dimension(:,:,:), allocatable :: p_spin
  ! temporary array for half-spinors
  complex(dpc) :: q_tmp(parallel%mydim)
  ! Kohn-Sham hamiltonian in the space of spinors
  complex(dpc), allocatable :: h_spin(:,:)

  real(dp) :: eval(2*elec_st%nstate)
  complex(dpc) :: qp(4)
  ! local scalars
  integer i,j,k, ii, jj, nn, nmax, info
  integer mdim, ldn, npe, tn, tmp1, tmp2
  integer comm, kpnum, kplp, mpinfo
  real(dp) magmom
  integer, allocatable :: indx_tmp(:)

  !---------------------------------------------------------------
  npe = parallel%group_size
  comm = parallel%group_comm
#ifdef DEBUG_1
  nloc_p_pot%cc = nloc_p_pot%cc*1.d5
#endif
#ifdef DEBUG_1
  write(9,*) ' pls A ',nloc_p_pot%so_num
  write(9,*) ' pls B ',nloc_p_pot%so_indx
  write(9,*) ' pls C ',nloc_p_pot%cc
#else
  write(9,*) "Inside pls"
#endif
  ! number of k-points
  kpnum = max(elec_st%nkpt,1)
  ldn = parallel%ldn
  mdim = parallel%mydim 
  nn = elec_st%nstate
  nmax = 2*nn
  allocate(h_spin(nmax,nmax))
  allocate(p_spin(mdim,elec_st%nstate,2))
  jj = 0
  do kplp = 1, kpnum
     if (elec_st%eig(1,kplp,1)%mm > jj) jj = elec_st%eig(1,kplp,1)%mm
     if (elec_st%eig(1,kplp,2)%mm > jj) jj = elec_st%eig(1,kplp,2)%mm
  enddo
  tn = solver%winsize + 2*jj
  do kplp = 1, kpnum
     !
     ! Only processors belonging to the group in charge of this
     ! representation should work.
     !
     if ( elec_st%eig(1,kplp,1)%group /= parallel%mygroup ) then
        call eig_e_adj_size(elec_st%eig(1,kplp,1),tn)
        cycle
     endif

     p_spin(:,:,:) = zzero
     h_spin(:,:) = zzero
     eval = zero
     !  first put the eigenvalues on the diagonal
     do i = 1, nn
        h_spin(i,i) = cmplx(elec_st%eig(1,kplp,1)%en(i),zero,dpc)
        k = i + nn
        h_spin(k,k) = cmplx(elec_st%eig(1,kplp,2)%en(i),zero,dpc)
     enddo

     if(elec_st%cplx) then
        do i = 1, 2            !spin
           do j = 1, nn        !number of states
              do k = 1, mdim   !elements in each wf, the last is zero
                 p_spin(k,j,i) = elec_st%eig(1,kplp,i)%zwf(k,j)
              enddo
           enddo
        enddo
     else
        do i = 1, 2            !spin
           do j = 1, nn        !number of states
              do k = 1, mdim   !elements in each wf, the last is zero
                 p_spin(k,j,i) = cmplx(elec_st%eig(1,kplp,i)%wf(k,j),zero,dpc)
              enddo
           enddo
        enddo
     endif
     do i = 1, nn
        do j = 1, nn
           qp(:) = zzero
#ifdef DEBUG_1
           write(9,*) ' LzSz ',i,j
#endif
           call lzsz(nloc_p_pot,parallel,p_spin(:,j,1),q_tmp(:),mdim,1,kplp)
           qp(1) = dot_product(p_spin(:,i,1),q_tmp(:))
           call zpsum(qp(1),1,npe,comm)
           h_spin(i,j)=h_spin(i,j) + conjg(qp(1))

           call lsxy(nloc_p_pot,parallel,p_spin(:,j,2),q_tmp,mdim,-1,kplp)
           qp(2) = dot_product(p_spin(:,i,1),q_tmp(:))
           call zpsum(qp(2),1,npe,comm)
           h_spin(i,j+nn) = h_spin(i,j+nn) + conjg(qp(2))

           ii = i + nn
           jj = j + nn

           call lzsz(nloc_p_pot,parallel,p_spin(:,j,2),q_tmp,mdim,-1,kplp)
           qp(3) = dot_product(p_spin(:,i,2),q_tmp(:))
           call zpsum(qp(3),1,npe,comm)
           h_spin(ii,jj) = h_spin(ii,jj) + conjg(qp(3))

           call lsxy(nloc_p_pot,parallel,p_spin(:,j,1),q_tmp,mdim,1,kplp)
           qp(4) = dot_product(p_spin(:,i,2),q_tmp(:))
           call zpsum(qp(4),1,npe,comm)
           h_spin(ii,j) = h_spin(ii,j) + conjg(qp(4))
#ifdef DEBUG_1
           if (parallel%iammaster) then
              write(7,*) ' pls ',i,j,sum(qp),qp
              call myflush(9)
           endif
#endif

        enddo
     enddo
#ifdef MPI
     call MPI_BARRIER(comm,mpinfo)
#endif
#ifdef DEBUG_1
     if (parallel%iammaster) then
        write(7,*) ' hmtrx_so = ',sum(h)
     endif
#endif

     call my_zheev('U', nmax, h_spin, nmax, eval, info)
#ifdef DEBUG_1
     if (parallel%iammaster) then
        write(7,*) ' hmtrx_eig = ',eval
     endif
#endif

     tmp1 = elec_st%eig(1,kplp,1)%nec + elec_st%eig(1,kplp,2)%nec
     tmp2 = elec_st%eig(1,kplp,1)%mm + elec_st%eig(1,kplp,2)%mm
     call destroy_eigenstate(elec_st%eig(1,kplp,2))

     call destroy_eigenstate(elec_st%eig(1,kplp,1))

     call create_eigenstate(elec_st%eig(1,kplp,1),2*ldn,tn,.true.)
 
     elec_st%eig(1,kplp,1)%en(1:nmax)= eval(1:nmax)

     elec_st%eig(1,kplp,1)%zwf = zzero

     do i = 1, nmax
        do j = 1, nn
           jj = j + nn
           elec_st%eig(1,kplp,1)%zwf(1:mdim,i) = elec_st%eig(1,kplp,1) &
                %zwf(1:mdim,i)+h_spin(j,i)*p_spin(1:mdim,j,1)

           elec_st%eig(1,kplp,1)%zwf(1+mdim:2*mdim,i) = elec_st%eig &
                (1,kplp,1)%zwf(1+mdim:2*mdim,i)+h_spin(jj,i)*p_spin(1:mdim,j,2)
        enddo
     enddo

     elec_st%eig(1,kplp,1)%nec = tmp1
     elec_st%eig(1,kplp,1)%mm = tmp2
     elec_st%eig(1,kplp,1)%occ = zero
     elec_st%eig(1,kplp,1)%occ(1:int(elec_st%xele)) = one
     elec_st%eig(1,kplp,1)%nn = nmax

     do i = 1, nmax
        magmom = zero
        do j = 1, mdim
           magmom = magmom + &
                abs(elec_st%eig(1,kplp,1)%zwf(j,i))**2 - &
                abs(elec_st%eig(1,kplp,1)%zwf(j+mdim,i))**2
        enddo
        call psum(magmom,1,npe,comm)
        elec_st%magmom(i,kplp) = magmom
     enddo

  enddo                     ! do kplp = 1, kpnum
  !
  ! Update eigenstate structures across MPI groups.
  !
  do kplp = 1, kpnum
     i = elec_st%eig(1,kplp,1)%group * parallel%group_size
#ifdef MPI
     call MPI_BCAST(elec_st%eig(1,kplp,1)%nec, &
          1,MPI_INTEGER,i,parallel%comm,mpinfo)
     call MPI_BCAST(elec_st%eig(1,kplp,1)%mm, &
          1,MPI_INTEGER,i,parallel%comm,mpinfo)
     call MPI_BCAST(elec_st%eig(1,kplp,1)%nn, &
          1,MPI_INTEGER,i,parallel%comm,mpinfo)
     call MPI_BCAST(elec_st%eig(1,kplp,1)%en, &
          tn,MPI_DOUBLE_PRECISION,i,parallel%comm,mpinfo)
#endif
  enddo

  if(.not. elec_st%cplx)elec_st%cplx=.true.
  elec_st%ifmax(1)  = elec_st%ifmax(1) + elec_st%ifmax(2)
  elec_st%ifmax(2) = 0
  elec_st%ntotal(1) = elec_st%ntotal(2) + elec_st%ntotal(1)
  elec_st%ntotal(2) = 0
  elec_st%nstate = 2 * elec_st%nstate
  jj = (solver%winsize + solver%nadd + nmax)*elec_st%nrep
  deallocate(elec_st%irep)
  allocate(elec_st%irep(jj,kpnum,2))

  elec_st%irep(:,:,1) = 1

  if (elec_st%nsave > 0) then
     nn = elec_st%nsave
     allocate(indx_tmp(nn))
     indx_tmp = elec_st%indxsave
     deallocate(elec_st%indxsave)
     elec_st%nsave = elec_st%nsave * 2
     allocate(elec_st%indxsave(elec_st%nsave))
     do jj = 1, nn
        elec_st%indxsave(2*jj - 1) = indx_tmp(jj)
        elec_st%indxsave(2*jj) = indx_tmp(jj)
     enddo
     deallocate(indx_tmp)
  endif

  deallocate(h_spin)
  deallocate(p_spin)
#ifdef MPI
  call MPI_BARRIER(parallel%comm,mpinfo)
#endif
end subroutine pls
!===============================================================

!---------------------------------------------------------------
subroutine lzsz(nloc_p_pot,parallel,p_spin,q_tmp,ndim,m,kplp)
  use constants
  use non_local_psp_module
  use parallel_data_module
  use pseudo_potential_module
#ifdef MPI
  ! include mpi definitions
  use mpi
#endif
  implicit none
  !
  ! Input/Output variables:
  !
  ! non local pseudopotential related data
  type (nonloc_pseudo_potential), intent(in) :: nloc_p_pot
  ! parallel computation related data
  type (parallel_data), intent(in) :: parallel
  ! input vectors in IW, distributed accross PEs
  integer, intent(in) :: ndim
  complex(dpc), intent(in) :: p_spin(ndim)
  ! output vectors in IW, distributed accross PEs
  complex(dpc), intent(out) :: q_tmp(ndim)
  ! sign of spin direction (+/-)
  integer, intent(in) :: m
  ! index of k-point
  integer, intent(in) :: kplp
  !
  ! Work variables:
  !
  complex(dpc) :: dotso(8,nloc_p_pot%so_num)
  complex(dpc) :: dotio(8,nloc_p_pot%so_num)
  real(dp) :: lz(8)
  complex(dpc) tmp
  !  real(dp) rcn1, rcn2
  integer i, ja, lm, lj, ml, ii, js, mpinfo
  !---------------------------------------------------------------
  dotio = zzero
  dotso = zzero
  lz = zero
  q_tmp = zzero
  i = 0
  do lj = 1, 2
     do ml = -lj, lj
        i = i + 1
        lz(i) = real(m,dp) * half * real(ml,dp)
     enddo
  enddo

  ! if NO k-point sampling ...
  if (nloc_p_pot%nkpt == 0) then
     do js = 1, nloc_p_pot%so_num
        ja = nloc_p_pot%so_indx(js)
        lj = 1
        do lm = 1, 8
           if (lm > 3) lj = 2
           do i = 1, nloc_p_pot%nlatom(ja)
              dotio(lm,js) = dotio(lm,js) + &
                   conjg(nloc_p_pot%v_ion(i,lm,js)) * &
                   p_spin(nloc_p_pot%indw(i,ja))

              dotso(lm,js) = dotso(lm,js) + &
                   conjg(nloc_p_pot%v_SO(i,lm,js)) * &
                   p_spin(nloc_p_pot%indw(i,ja))
           enddo
           dotso(lm,js) = dotso(lm,js)*nloc_p_pot%cc(lj,js)
           dotio(lm,js) = dotio(lm,js)*nloc_p_pot%cc(lj,js)
        enddo
     enddo
     ! in case of k-point sampling ..
  else
     do js = 1, nloc_p_pot%so_num
        ja = nloc_p_pot%so_indx(js)
        lj = 1
        do lm = 1, 8
           if (lm > 3) lj = 2
           do i = 1, nloc_p_pot%nlatom(ja)
              dotio(lm,js) = dotio(lm,js) + &
                   conjg(nloc_p_pot%v_ion(i,lm,js)) * &
                   p_spin(nloc_p_pot%indw(i,ja)) * &
                   nloc_p_pot%right(i,kplp,ja)

              dotso(lm,js) = dotso(lm,js) + &
                   conjg(nloc_p_pot%v_SO(i,lm,js)) * &
                   p_spin(nloc_p_pot%indw(i,ja)) * &
                   nloc_p_pot%right(i,kplp,ja)
           enddo
           dotso(lm,js) = dotso(lm,js)*nloc_p_pot%cc(lj,js)
           dotio(lm,js) = dotio(lm,js)*nloc_p_pot%cc(lj,js)
        enddo
     enddo
  endif ! if (nloc_p_pot%nkpt == 0)
#ifdef DEBUG_1
  write(9,*) ' dotso ',dotso
  write(9,*) ' dotio ',dotio
#endif
  call zpsum(dotso,8*nloc_p_pot%so_num, &
       parallel%group_size,parallel%group_comm)

  call zpsum(dotio,8*nloc_p_pot%so_num, &
       parallel%group_size,parallel%group_comm)

  do js = 1, nloc_p_pot%so_num
     ja=nloc_p_pot%so_indx(js)
     lj = 1
     do lm = 1, 8
        if( lm > 3) lj=2
        do i = 1, nloc_p_pot%nlatom(ja)
           ii = nloc_p_pot%indw(i,ja)
           if (nloc_p_pot%tran(i,ja) == 1) then
              tmp = zzero
              if(lz(lm) /= zero ) then
                 ! first terms that mixed projectors {|Vso><Vion|+h.c.}LS|P>
                 tmp = tmp + lz(lm) * dotso(lm,js) * &
                      nloc_p_pot%v_ion(i,lm,js)
                 tmp = tmp + lz(lm) * dotio(lm,js) *  &
                      nloc_p_pot%v_so(i,lm,js)
                 ! term of the projectors -0.5*|Vso><Vso|LS|P>
                 tmp = tmp - lz(lm) * half * dotso(lm,js)* &
                      nloc_p_pot%v_so(i,lm,js)
              endif         ! if(lz(lm) /= zero)

              ! finally, the constant term 0.25*L(L+1)|Vso><Vso|P>
              tmp = tmp + half * half * dotso(lm,js) * real(lj,dp)* &
                   (one+real(lj,dp)) * nloc_p_pot%v_so(i,lm,js)

              ! in case of k-point sampling ..
              if(nloc_p_pot%nkpt > 0) tmp = tmp * &
                   nloc_p_pot%left(i,kplp,ja)

              q_tmp(ii)=q_tmp(ii) + tmp

           endif            ! if (nloc_p_pot%tran(i,ja) == 1)
        enddo               ! do i = 1, nloc_p_pot%nlatom(ja)
     enddo                  ! do lm=1,8
  enddo                     ! do js = 1, nloc_p_pot%so_num
#ifdef MPI
  call MPI_BARRIER(parallel%group_comm,mpinfo)
#endif
end subroutine lzsz
!===================================================================

!-------------------------------------------------------------------
subroutine lsxy(nloc_p_pot,parallel,p_spin,q_tmp,ndim,m,kplp)
  use constants
  use non_local_psp_module
  use parallel_data_module
  use pseudo_potential_module
#ifdef MPI
  ! include mpi definitions
  use mpi
#endif
  implicit none
  !
  ! Input/Output variables:
  !
  ! non local pseudopotential related data
  type (nonloc_pseudo_potential), intent(in) :: nloc_p_pot
  ! parallel computation related data
  type (parallel_data), intent(in) :: parallel
  integer, intent(in) :: ndim
  ! input vectors in IW, distributed accross PEs
  complex(dpc), intent(in) :: p_spin(ndim)
  ! output vectors in IW, distributed accross PEs
  complex(dpc), intent(out) :: q_tmp(ndim)
  integer, intent(in) :: m !for "up" m=1, for "down" m=-1
  ! index of k-point
  integer, intent(in) :: kplp
  !
  ! Work variables:
  !
  complex(dpc) :: dotso(8,nloc_p_pot%so_num)
  complex(dpc) :: dotio(8,nloc_p_pot%so_num)
  real(dp) :: lpm(8)
  complex(dpc) tmp
  real(dp) :: numr
  integer i, ja, lm, lj, ml, ii, js, num, mpinfo
  !---------------------------------------------------------------
  dotso = zzero
  dotio = zzero
  lpm = zero
  q_tmp = zzero
  i = 0
  do lj = 1, 2
     do ml = -lj, lj
        i = i + 1
        num = lj*(lj+1)-ml*(m+ml)
        numr = real(num,dp)
        lpm(i) = half * sqrt(numr)
     enddo
  enddo
  ! if NO k-point sampling ...
  if (nloc_p_pot%nkpt == 0) then
     do js = 1, nloc_p_pot%so_num
        ja = nloc_p_pot%so_indx(js)
        lj = 1
        do lm = 1, 8
           if (lm > 3) lj = 2
           if (lpm(lm) /= zero) then
              do i = 1, nloc_p_pot%nlatom(ja)

                 dotio(lm,js) = dotio(lm,js) + &
                      conjg(nloc_p_pot%v_ion(i,lm,js)) * &
                      p_spin(nloc_p_pot%indw(i,ja)) 

                 dotso(lm,js) = dotso(lm,js) + &
                      conjg(nloc_p_pot%v_SO(i,lm,js)) * &
                      p_spin(nloc_p_pot%indw(i,ja)) 

              enddo
              dotso(lm,js) = dotso(lm,js)*nloc_p_pot%cc(lj,js)
              dotio(lm,js) = dotio(lm,js)*nloc_p_pot%cc(lj,js)
           endif
        enddo
     enddo
     ! in case of k-point sampling
  else
     do js = 1, nloc_p_pot%so_num
        ja = nloc_p_pot%so_indx(js)
        lj = 1
        do lm = 1, 8
           if (lm > 3) lj = 2
           if (lpm(lm) /= zero) then
              do i = 1, nloc_p_pot%nlatom(ja)

                 dotio(lm,js) = dotio(lm,js) + &
                      conjg(nloc_p_pot%v_ion(i,lm,js)) * &
                      p_spin(nloc_p_pot%indw(i,ja)) * &
                      nloc_p_pot%right(i,kplp,ja)

                 dotso(lm,js) = dotso(lm,js) + &
                      conjg(nloc_p_pot%v_SO(i,lm,js)) * &
                      p_spin(nloc_p_pot%indw(i,ja))* &
                      nloc_p_pot%right(i,kplp,ja)
              enddo
              dotso(lm,js) = dotso(lm,js)*nloc_p_pot%cc(lj,js)
              dotio(lm,js) = dotio(lm,js)*nloc_p_pot%cc(lj,js)
           endif
        enddo
     enddo
  endif ! if (nloc_p_pot%nkpt == 0)

#ifdef DEBUG_1
  write(9,*) ' dotso ',dotso
  write(9,*) ' dotio ',dotio
#endif
  call zpsum(dotso,8*nloc_p_pot%so_num, &
       parallel%group_size,parallel%group_comm)

  call zpsum(dotio,8*nloc_p_pot%so_num, &
       parallel%group_size,parallel%group_comm)

  do js = 1, nloc_p_pot%so_num
     ja = nloc_p_pot%so_indx(js)
     lj = 1
     do lm = 1, 8
        if (lpm(lm) /= zero ) then
           if ( lm > 3) lj=2
           do i = 1, nloc_p_pot%nlatom(ja)
              ii = nloc_p_pot%indw(i,ja)
              if (nloc_p_pot%tran(i,ja) == 1) then
                 tmp = zzero
                 ! first term (mixed projectors) {|Vso><Vion|+hc}LS|P>
                 tmp = tmp + lpm(lm) * dotso(lm,js) * &
                      nloc_p_pot%v_ion(i,m+lm,js)
                 tmp = tmp + lpm(lm) * dotio(lm,js) *  &
                      nloc_p_pot%v_so(i,m+lm,js)
                 ! second term -0.5|Vso><Vso|LS|P>
                 tmp = tmp - lpm(lm) * half * dotso(lm,js) * &
                      nloc_p_pot%v_so(i,m+lm,js)
                 ! in case of k-point sampling ..
                 if (nloc_p_pot%nkpt > 0) &
                      tmp = tmp * nloc_p_pot%left(i,kplp,ja)

                 q_tmp(ii)=q_tmp(ii) + tmp

              endif         ! if (nloc_p_pot%tran(i,ja) == 1)
           enddo            ! do i = 1, nloc_p_pot%nlatom(ja)
        endif               ! if(lpm(lm)/=zero)
     enddo                  ! do lm=1,8
  enddo                     ! do js = 1, nloc_p_pot%so_num
#ifdef MPI
  call MPI_BARRIER(parallel%group_comm,mpinfo)
#endif
end subroutine lsxy
!===============================================================

