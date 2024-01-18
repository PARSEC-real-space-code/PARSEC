!===============================================================
!
!  Copyright (C) 2005 Finite Difference Research Group
!  This file is part of parsec, http://www.ices.utexas.edu/parsec/
!
!  This subroutine calls subroutines to set up the exchange
!  correlation magnetic field B_xc for non-collinear calculations.
!  When the subroutine is called in the first time, magnetization
!  is taken from initial condition made by user (initcharge.f90).
!  In that case, it only calculates the correct stencile that 
!  operates on the trial functions (solver%bxc=Bxc\cdot\sigma.
!
!  If the subroutine is called not in the first time, it calculates
!  the magnetization with <wf|\vec{\sigma}|wf> (\sigma = Pauly vector
!  of matrices. In that case, the occupation of orbitals is distributed
!  to all PEs and the distribution in newrho.f90p is skipped. 
!  The difference in total electron number of up/dn population is 
!  defined as \int{(\rho_{up}-\rho_{dn})d^3r}=\int{|\vec{m}|d^3r}.
!  This value is put into elec_st%totel. The total magnetization
!  is calculated by integration over magnetization density.
!---------------------------------------------------------------
subroutine bxc(pot,parallel,elec_st,solver,first,nloc_p_pot)

  use constants
  use potential_module
  use parallel_data_module
  use electronic_struct_module
  use eigen_solver_module
  use non_local_psp_module
#ifdef MPI
  use mpi
#endif
  implicit none
  !
  !  Input/Output variables:
  !
  !  potential related data
  type (potential), intent(in) :: pot
  !  parallel computation related data
  type (parallel_data), intent(inout) :: parallel
  !  electronic structure
  type (electronic_struct), intent(inout) :: elec_st
  ! solver related data
  type (eigen_solver), intent(inout) :: solver
  ! non local pseudopotential related data
  type (nonloc_pseudo_potential), intent(in) :: nloc_p_pot
  ! 
  logical, intent(in) :: first
  !
  !  Work variables:
  !
  !  counters
  !
  integer, allocatable :: nlatom(:)
  integer :: i, j, kpnum, kplp, nn, nelems, mpinfo
  real(dp) :: tmp,kpwt,norm
  integer :: p(parallel%mydim)
  integer :: p1(parallel%mydim)
  integer, save :: kk = 0
  real(dp) ::q(parallel%mydim)
  complex(dpc), dimension(2*parallel%mydim) :: wf
  !---------------------------------------------------------------
  kpnum = max(elec_st%nkpt,1)
  wf = zzero
  kpwt = one
  if (first) then
     solver%bxc = zzero
     elec_st%spin3dn = zero
     do i = 1, parallel%mydim
        
        elec_st%spin3dn(i) = sqrt(elec_st%spin3d(i,1)**2 +  &
             elec_st%spin3d(i,2)**2 + elec_st%spin3d(i,3)**2)
        
        if (abs(elec_st%spin3dn(i)) > zero) then
           
           tmp = half * (pot%vxc(i,1)-pot%vxc(i,2)) / elec_st%spin3dn(i)
           
           solver%bxc(i,1) = tmp * cmplx(elec_st%spin3d(i,3),zero)
           
           solver%bxc(i,2) = tmp*cmplx(elec_st%spin3d(i,1),-elec_st%spin3d(i,2))
        endif
     enddo
  else
! Broadcast the occupations calculated in flevel.f90
#ifdef MPI
     do kplp = 1, kpnum
        nelems = elec_st%eig(1,kplp,1)%nn
        if (nelems == 0) cycle
        call MPI_BCAST(elec_st%eig(1,kplp,1)%occ,nelems, &
             MPI_DOUBLE_PRECISION,parallel%masterid, &
             parallel%comm,mpinfo)
     enddo
#endif
     ! Calculate the magnetic moment
     !
     elec_st%spin3d = zero
     do kplp = 1,kpnum
        if ( elec_st%eig(1,kplp,1)%group /= parallel%mygroup ) cycle
        do nn = 1, elec_st%eig(1,kplp,1)%nn
           !
           ! retrive total wave-function
           call zcopy(2*parallel%mydim,elec_st%eig(1,kplp,1)%zwf(1,nn),1,wf,1)
           !
           ! create the 3D spin magnetization and it's norm
           do i = 1,parallel%mydim
              j = i + parallel%mydim
              elec_st%spin3d(i,1) =  elec_st%spin3d(i,1) + &
                   (two*(real(wf(i))*real(wf(j))+aimag(wf(i))*aimag(wf(j)))) * &
                   elec_st%eig(1,kplp,1)%occ(nn)*elec_st%kpwt(kplp)
              
              elec_st%spin3d(i,2) = elec_st%spin3d(i,2) + &
                   (two*(real(wf(i))*aimag(wf(j))-real(wf(j))*aimag(wf(i)))) * &
                   elec_st%eig(1,kplp,1)%occ(nn)*elec_st%kpwt(kplp)
              
              elec_st%spin3d(i,3) = elec_st%spin3d(i,3) + &
                   (abs(wf(i))*abs(wf(i)) - abs(wf(j))*abs(wf(j))) * &
                   elec_st%eig(1,kplp,1)%occ(nn)*elec_st%kpwt(kplp)
           enddo
        enddo
     enddo
     elec_st%spin3dn = zero
     do i = 1,parallel%mydim
        tmp = dot_product(elec_st%spin3d(i,1:3),elec_st%spin3d(i,1:3))
        elec_st%spin3dn(i) = sqrt(tmp)
     enddo
  endif
  
  elec_st%net_magmom = zero
  elec_st%net_magmom = sum(elec_st%spin3dn(:))
  
  call psum(elec_st%net_magmom,1,parallel%group_size,parallel%group_comm)  
  
#ifdef MPI
  ! AJB: uneeded, psum is a barrier.
  !call MPI_BARRIER(parallel%comm,mpinfo)
#endif
  
  call  nclmom(pot,parallel,elec_st,solver,nloc_p_pot)
  
  do i = 1,2
     tmp = (-one)**i
     elec_st%totel(i) = half*(elec_st%xele - tmp*elec_st%net_magmom)
  enddo
  
  write(7,*)
  write(7,52) elec_st%net_magmom
  
52 format(1x,'Integral{|m|} =',1x,f7.2)
  
end subroutine bxc
!===============================================================

!===============================================================
!
!  Copyright (C) 2005 Finite Difference Research Group
!  This file is part of parsec, http://www.ices.utexas.edu/parsec/
!
!  This subroutine calls subroutines to set up the exchange
!  correlation magnetic field B_xc for non-collinear calculations
!  
!---------------------------------------------------------------
subroutine nclmom(pot,parallel,elec_st,solver,nloc_p_pot)

  use constants
  use potential_module
  use parallel_data_module
  use electronic_struct_module
  use eigen_solver_module
  use non_local_psp_module
#ifdef MPI
  use mpi
#endif
  implicit none
  !
  !  Input/Output variables:
  !
  !  potential related data
  type (potential), intent(in) :: pot
  !  parallel computation related data
  type (parallel_data), intent(inout) :: parallel
  !  electronic structure
  type (electronic_struct), intent(inout) :: elec_st
  ! solver related data
  type (eigen_solver), intent(inout) :: solver
  ! non local pseudopotential related data
  type (nonloc_pseudo_potential), intent(in) :: nloc_p_pot
  !
  !  Work variables:
  !
  !  counters
  !
  integer :: i, j, kpnum, kplp, nn, ja, ii, nlatm, w
  real(dp) :: tmp, kpwt
  real(dp) ::atomspin(nloc_p_pot%atom_num,3)
  integer, allocatable :: nlatom(:)
  !---------------------------------------------------------------
  !
  ! calculate the atomic magnetic moment
  ! integrating spin3d in a sphere of WS radius around each atom
  allocate(nlatom(nloc_p_pot%atom_num))
  nlatom(:) = nloc_p_pot%wsnlatom(:)
  atomspin = zero
  call pisum(nlatom,nloc_p_pot%atom_num,parallel%group_size,parallel%group_comm)
  do ja = 1, nloc_p_pot%atom_num
     nlatm = nloc_p_pot%wsnlatom(ja)
     do i = 1, nlatm
        ii = nloc_p_pot%wsindw(i,ja)
        atomspin(ja,1:3) = atomspin(ja,1:3) + elec_st%spin3d(ii,1:3)
     enddo
  enddo
 
  call psum(atomspin,3*nloc_p_pot%atom_num,parallel%group_size,parallel%group_comm)

  write(7,*)
  write(7,42)
  write(7,*)
  write(7,41)
  do ja = 1, nloc_p_pot%atom_num
     write(7,40)ja,atomspin(ja,1),atomspin(ja,2),atomspin(ja,3)
  enddo
  !
  ! compute the total magnetic moment
  elec_st%s3dtot = zero
  do j = 1,3
     elec_st%s3dtot(j) = sum(elec_st%spin3d(1:parallel%mydim,j))
  enddo
  call psum(elec_st%s3dtot,3,parallel%group_size,parallel%group_comm)
  write(7,*)
  write(7,*)
  write(7,87)elec_st%s3dtot
  deallocate(nlatom)

40 format(i4,2x,3(f7.2,1x))
41 format('-atom-  --Mx--  --My--  --Mz--' ,'   [bohr magnetons]       ',/)
42 format('Atomic magnetic moments (Bohr magnetons): ')
88 format(2(i4,3x),3x,f7.2)
87 format(' Tot. Magnetization (Bohr magnetons)', 3(5x,f7.2))
end subroutine nclmom
!===============================================================

