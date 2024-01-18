!===============================================================
!
! Copyright (C) 2005 Finite Difference Research Group
! This file is part of parsec, http://www.ices.utexas.edu/parsec/
!
!---------------------------------------------------------------
subroutine forcso(clust,elec_st,nloc_p_pot,parallel)
  use constants
  use cluster_module
  use electronic_struct_module
  use non_local_psp_module
  use parallel_data_module
#ifdef MPI
  use mpi
#endif
  implicit none
  !
  ! Input/Output variables:
  !
  ! the cluster
  type (cluster), intent(inout) :: clust
  ! electronic structure
  type (electronic_struct), intent(inout) :: elec_st
  ! non local pseudopotential related data
  type (nonloc_pseudo_potential), intent(inout) :: nloc_p_pot 
  ! parallel computation related data
  type (parallel_data), intent(in) :: parallel
  !
  ! Work variables:
  !
  ! distributed spinor (two-component wavefunctions)
  complex(dpc),dimension(:,:), allocatable :: p_spin
  ! temporary array for half-spinors
  complex(dpc) :: f_tmp
  ! weight of each kpoint being examined
  real(dp), dimension(:), allocatable :: weight
  !
  ! local scalars
  integer :: i, j, ja, nn
  integer :: mdim, ldn, npe
  integer :: comm, kpnum, kplp, mpinfo
  !---------------------------------------------------------------
  npe = parallel%group_size
  comm = parallel%group_comm
  ldn = parallel%ldn
  mdim = parallel%mydim 
  nn=elec_st%nstate

  kpnum = max(elec_st%nkpt,1)
  allocate(weight(kpnum))
  if (elec_st%nkpt > 0) then
     weight = elec_st%kpwt
  else
     weight = one
  endif

  clust%force_so(:,:) = zzero
  allocate(p_spin(mdim,2))
  do ja = 1, nloc_p_pot%so_num
     do j=1,3
        do kplp = 1, kpnum
           do i=1,nn 
  
              if ( elec_st%eig(1,kplp,1)%group /= parallel%mygroup ) cycle
   
              p_spin(:,1) = elec_st%eig(1,kplp,1)%zwf(1:mdim,i)
              p_spin(:,2) = elec_st%eig(1,kplp,1)%zwf(1+mdim:2*mdim,i)

              call flzsz(nloc_p_pot,parallel,p_spin(:,1),f_tmp,mdim,1,ja,j,kplp)
              clust%force_so(j,ja) = clust%force_so(j,ja) + f_tmp &
              *elec_st%eig(1,kplp,1)%occ(i)*weight(kplp)
              

              call flsxy(nloc_p_pot,parallel,p_spin(:,:),f_tmp,mdim,-1,ja,j,kplp)
              clust%force_so(j,ja) = clust%force_so(j,ja) + f_tmp &
              *elec_st%eig(1,kplp,1)%occ(i)*weight(kplp)
              
              
              call flzsz(nloc_p_pot,parallel,p_spin(:,2),f_tmp,mdim,-1,ja,j,kplp)
              clust%force_so(j,ja) = clust%force_so(j,ja) + f_tmp  &
              *elec_st%eig(1,kplp,1)%occ(i)*weight(kplp)

              
              call flsxy(nloc_p_pot,parallel,p_spin(:,:),f_tmp,mdim,1,ja,j,kplp)
              clust%force_so(j,ja) = clust%force_so(j,ja) + f_tmp &
              *elec_st%eig(1,kplp,1)%occ(i)*weight(kplp)

              
           enddo    !do i=1,nn
        enddo       !do kplp = 1, kpnum
     enddo          !do j=1,3

!#ifdef DEBUG
!  write(7,*) 'force 1',ja,clust%force_so(1,ja)
!  write(7,*) 'force 2',ja,clust%force_so(2,ja)
!  write(7,*) 'force 3',ja,clust%force_so(3,ja)
!#endif

  enddo             !do ja = 1, nloc_p_pot%so_num

  deallocate(p_spin)
  deallocate(weight)

#ifdef MPI
  call MPI_BARRIER(parallel%comm,mpinfo)
#endif
end subroutine forcso
!===============================================================

!---------------------------------------------------------------
subroutine flzsz(nloc_p_pot,parallel,p_spin,force,ndim,m,js,j,kplp)
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
  complex(dpc), intent(out) :: force
  ! sign of spin direction (+/-)
  integer, intent(in) :: m
  ! index of of atom 
  integer, intent(in) :: js
  ! index of derivative direction
  integer, intent(in) :: j
  ! index of k-point
  integer, intent(in) :: kplp
  !
  ! Work variables:
  !
  complex(dpc) :: ddotso(8)
  complex(dpc) :: ddotio(8)
  complex(dpc) :: dotso(8)
  complex(dpc) :: dotio(8)
  real(dp) :: lz(8)
  complex(dpc) tmp
  !  real(dp) rcn1, rcn2
  integer i, ja, lm, lj, ml, mpinfo
  !---------------------------------------------------------------
  force = zzero
  dotio = zzero
  dotso = zzero
  ddotso = zzero
  ddotio = zzero
  ja=nloc_p_pot%so_indx(js)
  lz = zero
  i = 0
  do lj=1,2
     do ml = -lj, lj
        i = i + 1
        lz(i) = real(m,dp) * half * real(ml,dp)
     enddo
  enddo

  ! if NO k-point sampling ...
  if (nloc_p_pot%nkpt == 0) then
     lj = 1
     do lm = 1, 8
        if (lm > 3) lj = 2
        do i = 1, nloc_p_pot%nlatom(ja)
           dotio(lm) = dotio(lm) + &
                conjg(nloc_p_pot%v_ion(i,lm,js)) * &
                p_spin(nloc_p_pot%indw(i,ja))

           dotso(lm) = dotso(lm) + &
                conjg(nloc_p_pot%v_SO(i,lm,js)) * &
                p_spin(nloc_p_pot%indw(i,ja))

           ddotio(lm) = ddotio(lm) + &
                conjg(nloc_p_pot%dv_ion(j,i,lm,js)) * &
                p_spin(nloc_p_pot%indw(i,ja)) 

           ddotso(lm) = ddotso(lm) + &
                conjg(nloc_p_pot%dv_SO(j,i,lm,js)) * &
                p_spin(nloc_p_pot%indw(i,ja))
        enddo
        ddotso(lm) = ddotso(lm)*nloc_p_pot%cc(lj,js)
        ddotio(lm) = ddotio(lm)*nloc_p_pot%cc(lj,js)

        dotso(lm) = dotso(lm)*nloc_p_pot%cc(lj,js)
        dotio(lm) = dotio(lm)*nloc_p_pot%cc(lj,js)
     enddo
     ! in case of k-point sampling ..
  else
     lj = 1
     do lm = 1, 8
        if (lm > 3) lj = 2
        do i = 1, nloc_p_pot%nlatom(ja)
           dotio(lm) = dotio(lm) + &
                conjg(nloc_p_pot%v_ion(i,lm,js)) * &
                p_spin(nloc_p_pot%indw(i,ja)) * &
                nloc_p_pot%right(i,kplp,ja)
           
           dotso(lm) = dotso(lm) + &
                conjg(nloc_p_pot%v_SO(i,lm,js)) * &
                p_spin(nloc_p_pot%indw(i,ja)) * &
                nloc_p_pot%right(i,kplp,ja)

           ddotio(lm) = ddotio(lm) + &
                conjg(nloc_p_pot%dv_ion(j,i,lm,js)) * &
                p_spin(nloc_p_pot%indw(i,ja)) * &
                nloc_p_pot%right(i,kplp,ja)

           ddotso(lm) = ddotso(lm) + &
                conjg(nloc_p_pot%dv_SO(j,i,lm,js)) * &
                p_spin(nloc_p_pot%indw(i,ja)) * &
                nloc_p_pot%right(i,kplp,ja)

        enddo
        ddotso(lm) = ddotso(lm)*nloc_p_pot%cc(lj,js)
        ddotio(lm) = ddotio(lm)*nloc_p_pot%cc(lj,js)

        dotso(lm) = dotso(lm)*nloc_p_pot%cc(lj,js)
        dotio(lm) = dotio(lm)*nloc_p_pot%cc(lj,js)
     enddo
  endif ! if (nloc_p_pot%nkpt == 0)

  call zpsum(ddotso,8,parallel%group_size,parallel%group_comm)

  call zpsum(ddotio,8,parallel%group_size,parallel%group_comm)

  call zpsum(dotso,8,parallel%group_size,parallel%group_comm)

  call zpsum(dotio,8,parallel%group_size,parallel%group_comm)

  lj = 1
  do lm = 1, 8
     if( lm > 3) lj=2
     tmp = zzero
     if(lz(lm) /= zero ) then
        ! first terms that mixed projectors {|Vso><Vion|+h.c.}LS|P>
        tmp = tmp + lz(lm) * dotso(lm) * conjg(ddotio(lm))
        tmp = tmp + lz(lm) * dotio(lm) * conjg(ddotso(lm))
        tmp = tmp + lz(lm) * ddotso(lm) * conjg(dotio(lm))
        tmp = tmp + lz(lm) * ddotio(lm) * conjg(dotso(lm))
        ! term of the projectors -0.5*|Vso><Vso|LS|P>
        tmp = tmp - lz(lm) * half * ddotso(lm) * conjg(dotso(lm))
        tmp = tmp - lz(lm) * half * dotso(lm) * conjg(ddotso(lm))
     endif         ! if(lz(lm) /= zero)
     
     ! finally, the constant term 0.25*L(L+1)|Vso><Vso|P>
     tmp = tmp + half * half * dotso(lm) * real(lj,dp)* &
          (one+real(lj,dp)) * conjg(ddotso(lm))
     tmp = tmp + half * half * ddotso(lm) * real(lj,dp)* &
          (one+real(lj,dp)) * conjg(dotso(lm))
     ! in case of k-point sampling ..
!     if(nloc_p_pot%nkpt > 0) tmp = tmp * &
!          nloc_p_pot%left(i,kplp,ja)
     
     force = force + tmp
     
  enddo                  ! do lm=1,8

#ifdef MPI
  call MPI_BARRIER(parallel%group_comm,mpinfo)
#endif
end subroutine flzsz
!===================================================================

!-------------------------------------------------------------------
subroutine flsxy(nloc_p_pot,parallel,p_spin,force,ndim,m,js,j,kplp)
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
  complex(dpc), intent(in) :: p_spin(ndim,2)
  ! output force
  complex(dpc), intent(out) :: force
  integer, intent(in) :: m !for "up" m=1, for "down" m=-1
  ! index of derivative direction
  integer, intent(in) :: js
  ! index of derivative direction
  integer, intent(in) :: j
  ! index of k-point
  integer, intent(in) :: kplp
  !
  ! Work variables:
  !
  integer isp(2)
  complex(dpc) :: ddotso(8,2)
  complex(dpc) :: ddotio(8,2)
  complex(dpc) :: dotso(8,2)
  complex(dpc) :: dotio(8,2)
  real(dp) :: lpm(8)
  complex(dpc) tmp
  real(dp) :: numr
  integer i, ja, lm, lj, ml, ii, num, mpinfo
  !---------------------------------------------------------------
  force = zzero
  ddotso = zzero
  ddotio = zzero
  dotso = zzero
  dotio = zzero
  ja=nloc_p_pot%so_indx(js)
  lpm = zero
  if (m == 1) isp(1) = 1
  if (m == -1) isp(1) = 2
  isp(2) = 2/isp(1)
  i = 0
  do lj=1,2
     do ml = -lj, lj
        i = i + 1
        num = lj*(lj+1)-ml*(m+ml)
        numr = real(num,dp)
        lpm(i) = half * sqrt(numr)
     enddo
  enddo
  ! if NO k-point sampling ...
  if (nloc_p_pot%nkpt == 0) then
     lj = 1
     do ii = 1, 2
        do lm = 1, 8
           if (lm > 3) lj = 2
           do i = 1, nloc_p_pot%nlatom(ja)
              dotio(lm,ii) = dotio(lm,ii) + &
                conjg(nloc_p_pot%v_ion(i,lm,js)) * &
                p_spin(nloc_p_pot%indw(i,ja),isp(ii))

              dotso(lm,ii) = dotso(lm,ii) + &
                   conjg(nloc_p_pot%v_SO(i,lm,js)) * &
                   p_spin(nloc_p_pot%indw(i,ja),isp(ii))
              
              ddotio(lm,ii) = ddotio(lm,ii) + &
                   conjg(nloc_p_pot%dv_ion(j,i,lm,js)) * &
                   p_spin(nloc_p_pot%indw(i,ja),isp(ii))

              ddotso(lm,ii) = ddotso(lm,ii) + &
                   conjg(nloc_p_pot%dv_SO(j,i,lm,js)) * &
                   p_spin(nloc_p_pot%indw(i,ja),isp(ii))
           enddo
           ddotso(lm,ii) = ddotso(lm,ii)*nloc_p_pot%cc(lj,js)
           ddotio(lm,ii) = ddotio(lm,ii)*nloc_p_pot%cc(lj,js)
           
           dotso(lm,ii) = dotso(lm,ii)*nloc_p_pot%cc(lj,js)
           dotio(lm,ii) = dotio(lm,ii)*nloc_p_pot%cc(lj,js)
        enddo    ! do lm = 1, 8
     enddo       !  do ii = 1, 2
     ! in case of k-point sampling ..
  else
     lj = 1
     do ii = 1, 2
        do lm = 1, 8
           if (lm > 3) lj = 2
           do i = 1, nloc_p_pot%nlatom(ja)
              dotio(lm,ii) = dotio(lm,ii) + &
                   conjg(nloc_p_pot%v_ion(i,lm,js)) * &
                   p_spin(nloc_p_pot%indw(i,ja),isp(ii)) * &
                   nloc_p_pot%right(i,kplp,ja)
              
              dotso(lm,ii) = dotso(lm,ii) + &
                   conjg(nloc_p_pot%v_SO(i,lm,js)) * &
                   p_spin(nloc_p_pot%indw(i,ja),isp(ii)) * &
                   nloc_p_pot%right(i,kplp,ja)
              
              ddotio(lm,ii) = ddotio(lm,ii) + &
                   conjg(nloc_p_pot%dv_ion(j,i,lm,js)) * &
                   p_spin(nloc_p_pot%indw(i,ja),isp(ii)) * &
                   nloc_p_pot%right(i,kplp,ja)
              
              ddotso(lm,ii) = ddotso(lm,ii) + &
                   conjg(nloc_p_pot%dv_SO(j,i,lm,js)) * &
                   p_spin(nloc_p_pot%indw(i,ja),isp(ii)) * &
                   nloc_p_pot%right(i,kplp,ja)
              
           enddo
           ddotso(lm,ii) = ddotso(lm,ii)*nloc_p_pot%cc(lj,js)
           ddotio(lm,ii) = ddotio(lm,ii)*nloc_p_pot%cc(lj,js)
           
           dotso(lm,ii) = dotso(lm,ii)*nloc_p_pot%cc(lj,js)
           dotio(lm,ii) = dotio(lm,ii)*nloc_p_pot%cc(lj,js)
        enddo
     enddo
  endif ! if (nloc_p_pot%nkpt == 0)
  
  call zpsum(ddotso,16,parallel%group_size,parallel%group_comm)

  call zpsum(ddotio,16,parallel%group_size,parallel%group_comm)
  
  call zpsum(dotso,16,parallel%group_size,parallel%group_comm)
  
  call zpsum(dotio,16,parallel%group_size,parallel%group_comm)
  
  
  lj = 1
  do lm = 1, 8
     if(lpm(lm) /= zero ) then
        if( lm > 3) lj=2
        tmp = zzero
        ! first term (mixed projectors) {|Vso><Vion|+hc}LS|P>
        tmp = tmp + lpm(lm) * dotso(lm,1) * conjg(ddotio(lm,2))
        tmp = tmp + lpm(lm) * ddotso(lm,1) * conjg(dotio(lm,2))
        tmp = tmp + lpm(lm) * dotio(lm,1) * conjg(ddotso(lm,2))
        tmp = tmp + lpm(lm) * ddotio(lm,1) * conjg(dotso(lm,2))
        ! second term -0.5|Vso><Vso|LS|P>
        tmp = tmp - lpm(lm) * half * dotso(lm,1) * conjg(ddotso(lm,2))
        tmp = tmp - lpm(lm) * half * ddotso(lm,1) * conjg(dotso(lm,2))
        ! in case of k-point sampling ..
!        if(nloc_p_pot%nkpt > 0) &
!             tmp = tmp * nloc_p_pot%left(i,kplp,ja)
        
        force = force + tmp

     endif               ! if(lpm(lm)/=zero)
  enddo                  ! do lm=1,8
#ifdef MPI
  call MPI_BARRIER(parallel%group_comm,mpinfo)
#endif
end subroutine flsxy
!===============================================================

