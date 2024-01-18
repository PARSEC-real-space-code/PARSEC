!===============================================================
!
! Copyright (C) 2005 Finite Difference Research Group
! This file is part of parsec, http://www.ices.utexas.edu/parsec/
!
! Matrix-vector multiplication for the spin-orbit component of
! pseudopotentials.
!
!---------------------------------------------------------------
subroutine zmatvec_so(kplp,nloc_p_pot,u_pot,solver,parallel,vnew,p,q,m,totldn)

  use constants
  use matvec_module
  use matvec_interface
  use non_local_psp_module
  use eigen_solver_module
  use parallel_data_module
  implicit none
  !
  ! Input/Output variables:
  !
  ! non local pseudopotential related data
  type (nonloc_pseudo_potential), intent(in) :: nloc_p_pot 
  ! U potential related data
  type (nonloc_pseudo_potential), intent(in) :: u_pot
  ! solver related data
  type (eigen_solver), intent(inout) :: solver
  ! parallel computation related data
  type (parallel_data), intent(in) :: parallel
  ! vnew
  real(dp), intent(in) :: vnew(parallel%mydim,2)
  ! number of vectors
  integer, intent(in) :: m
  ! length of vectors
  integer, intent(in) :: totldn
  ! index of k-point
  integer, intent(in) :: kplp
  ! input vectors in IW, distributed accross PEs
  complex(dpc), intent(in) :: p(totldn,m)
  !complex(dpc), intent(in) :: p(2*parallel%ldn,m)
  ! output vectors in IW, distributed accross PEs
  complex(dpc), intent(out) :: q(totldn,m)
  !complex(dpc), intent(out) :: q(2*parallel%ldn,m)
  !
  ! Work variables:
  real(dp), dimension(:), allocatable :: vud
  complex(dpc), dimension(:,:), allocatable :: qup, qdn, pup, pdn
  !
  ! local scalars
  integer i, ldn, ldn2, dim, dim2, k
  !---------------------------------------------------------------
  !ldn = parallel%ldn
  ldn2 = totldn
  ldn = ldn2/2
  allocate(qup(ldn,m))
  allocate(qdn(ldn,m))
  allocate(pup(ldn,m))
  allocate(pdn(ldn,m))
  !allocate(pdn(parallel%ldn,m))
  q(:,:) = zzero
  qup(:,:) = zzero
  qdn(:,:) = zzero
  pup(:,:) = zzero
  pdn(:,:) = zzero

  dim = parallel%mydim
  dim2 = 2*parallel%mydim

  pup(1:dim,1:m) = p(1:dim,1:m)
  !AJB: umm, 1+ldn or 1+dim?
  pdn(1:dim,1:m) = p(1+dim:dim2,1:m)
  if(solver%ncl)then
     allocate(vud(dim))
     vud = half * (vnew(1:dim,1) + vnew(1:dim,2))
  endif
  ! applying the scalar-relativistic Hamiltonian on the "up" part
  !------------------------------------------------------------------
  if (nloc_p_pot%nkpt == 0) then
     if(solver%ncl)then
        call dcopy(parallel%mydim,vud(1),1,solver%adiag,1,ldn)
     else
        call dcopy(parallel%mydim,vnew(1,1),1,solver%adiag,1)
     endif
     solver%adiag(:) = solver%adiag(:) + sum(solver%coe2(0,:))
  else
     if(solver%ncl)then
        call dcopy(parallel%mydim,vud(1),1,solver%adiag,1)
     else
        call dcopy(parallel%mydim,vnew(1,1),1,solver%adiag,1)
     endif
     solver%adiag(:) = solver%adiag(:) + sum(solver%coe2(0,:)) + &
          nloc_p_pot%kpts(1,kplp)*nloc_p_pot%kpts(1,kplp) + &
          nloc_p_pot%kpts(2,kplp)*nloc_p_pot%kpts(2,kplp) + &
          nloc_p_pot%kpts(3,kplp)*nloc_p_pot%kpts(3,kplp)
  endif
  call zmatvec_nloc(kplp, nloc_p_pot, u_pot, solver, parallel, pup, qup, m,ldn)
  if(nloc_p_pot%nkpt == 0) then
     call Zmatvec_ke(solver, parallel, pup, qup, m,ldn)
  else
     call Zmatvec_ke_kpt(kplp, solver, parallel, pup, qup, m,ldn)
  endif
  ! applying the scalar-relativistic Hamiltonian on the "dn" part
  !------------------------------------------------------------------
  if (nloc_p_pot%nkpt == 0) then
     if(solver%ncl)then
        call dcopy(parallel%mydim,vud(1),1,solver%adiag,1)
     else
        call dcopy(parallel%mydim,vnew(1,2),1,solver%adiag,1)
     endif
        solver%adiag(:) = solver%adiag(:) + sum(solver%coe2(0,:))
  else
     if(solver%ncl)then
        call dcopy(parallel%mydim,vud(1),1,solver%adiag,1)
     else
        call dcopy(parallel%mydim,vnew(1,2),1,solver%adiag,1)
     endif
     solver%adiag(:) = solver%adiag(:) + sum(solver%coe2(0,:)) + &
          nloc_p_pot%kpts(1,kplp)*nloc_p_pot%kpts(1,kplp) + &
          nloc_p_pot%kpts(2,kplp)*nloc_p_pot%kpts(2,kplp) + &
          nloc_p_pot%kpts(3,kplp)*nloc_p_pot%kpts(3,kplp)
  endif

  call Zmatvec_nloc(kplp, nloc_p_pot, u_pot, solver, parallel, pdn, qdn, m,ldn)
  if(nloc_p_pot%nkpt == 0) then
     call Zmatvec_ke(solver, parallel, pdn, qdn, m,ldn)
  else
     call Zmatvec_ke_kpt(kplp, solver, parallel, pdn, qdn, m,ldn)
  endif
  q(1:dim,1:m) = qup(1:dim,1:m)
  q(1+dim:dim2,1:m) = qdn(1:dim,1:m)
  !
  ! if non-collinear calculation, apply the B_xc \cdot vec{sigma}
  !--------------------------------------------------------------
  if (solver%ncl) then
     do i = 1, parallel%mydim
        do k = 1,m
           q(i,k) = q(i,k) + solver%bxc(i,1) * pup(i,k)
           q(i,k) = q(i,k) + solver%bxc(i,2) * pdn(i,k)
           
           q(i+dim,k) = q(i+dim,k) - solver%bxc(i,1) * pdn(i,k)
           q(i+dim,k) = q(i+dim,k) + conjg(solver%bxc(i,2)) * pup(i,k)
        enddo
     enddo
     deallocate(vud)
  endif
  ! applying the spin-orbit Hamiltonian on the "up" part
  !---------------------------------------------------------
  if (nloc_p_pot%is_so) then
     qup(:,:) = zzero
     qdn(:,:) = zzero
     do i = 1, m
        call lzsz(nloc_p_pot,parallel,pup(1,i),qup(1,i),dim,1,kplp)
        q(1:dim,i) = q(1:dim,i) + qup(1:dim,i)
        call lsxy(nloc_p_pot,parallel,pup(1,i),qdn(1,i),dim,1,kplp)
        q(1+dim:dim2,i) = q(1+dim:dim2,i) + qdn(1:dim,i)
     enddo
     
     ! applying the spin-orbit Hamiltonian on the "dn" part      
     !-----------------------------------------------------------
     qup(:,:) = zzero
     qdn(:,:) = zzero
     do i = 1, m
        call lzsz(nloc_p_pot,parallel,pdn(1,i),qdn(1,i),dim,-1,kplp)
        q(1+dim:dim2,i) = q(1+dim:dim2,i) + qdn(1:dim,i)
        call lsxy(nloc_p_pot,parallel,pdn(1,i),qup(1,i),dim,-1,kplp)
        q(1:dim,i) = q(1:dim,i) + qup(1:dim,i)
     enddo
  endif
  deallocate(qup, qdn, pup, pdn)
end subroutine zmatvec_so
!===============================================================

