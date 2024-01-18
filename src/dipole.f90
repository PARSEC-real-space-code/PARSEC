!===============================================================
!
! Copyright (C) 2005 Finite Difference Research Group
! This file is part of parsec, http://www.ices.utexas.edu/parsec/
!
! This subroutine calculates the total charge and the total 
! dipole of the system (i.e., the zeroth and first multipole
! terms), based on the charge density.
!
!---------------------------------------------------------------
subroutine dipole(clust,rsymm,elec_st,grid,parallel,zion)

  use constants
  use cluster_module
  use symmetry_module
  use electronic_struct_module
  use grid_module
  use parallel_data_module
  implicit none
  !
  ! Input/Output variables:
  !
  ! the cluster
  type (cluster), intent(in) :: clust
  ! symmetry operations in reduced group:
  type (symmetry), intent(in) :: rsymm
  ! electronic structure
  type (electronic_struct), intent(inout) :: elec_st
  ! grid related data
  type (grid_data), intent(in) :: grid
  ! parallel computation related data
  type (parallel_data), intent(in) :: parallel

  ! positive charge of pseudo ion for each atom type
  real(dp), intent(in) :: zion(1:clust%type_num)
  ! 
  ! Work variables:
  !
  ! charge of cluster
  real(dp) :: charge
  real(dp),dimension(1) :: chargevec
  ! dipole along the three axes
  real(dp) :: dipr(3)
  ! counters
  integer i,j,ity,itrans,ioffset
  real(dp) :: rhoi, rw(3), rrw(3)

  !---------------------------------------------------------------
  !
  ! Sum over all grid points to get electron contribution to the
  ! charge and dipole. Note that the charge density, rho, is
  ! calculated elsewhere in the code as a positive quantity, even
  ! though it is of course negative - hence the minus signs in rhoi.
  !
  charge = zero
  dipr(:) = zero

  ioffset = parallel%irows(parallel%group_iam) - 1
  do j = parallel%irows(parallel%group_iam) &
       ,parallel%irows(parallel%group_iam+1)-1
     rhoi = -elec_st%rho(j-ioffset,1)
     charge = charge + rhoi
     rw(1) = (grid%shift(1) + grid%kx(j)) * grid%step(1)
     rw(2) = (grid%shift(2) + grid%ky(j)) * grid%step(2)
     rw(3) = (grid%shift(3) + grid%kz(j)) * grid%step(3)
     do itrans = 1, rsymm%ntrans
        call symop(rsymm,itrans,rw,rrw)
        dipr = dipr + rrw * rhoi
     enddo
  enddo
  chargevec = charge
  call psum(chargevec,1,parallel%group_size,parallel%group_comm)
  charge = chargevec(1)
  call psum(dipr,3,parallel%group_size,parallel%group_comm)

  ! Multiply by h^3 to approximate volume integration by the
  ! summation; charge has the extra ntrans factor because it was
  ! summed over the wedge above.

  charge = charge*grid%hcub*real(rsymm%ntrans,dp)
  dipr = dipr*grid%hcub

  if (parallel%iammaster) then
     ! Sum over all atoms to get the nuclei contribution to the 
     ! charge and dipole.
     j = 0
     do ity = 1, clust%type_num
        do i  = 1, clust%natmi(ity)
           j = j + 1
           charge = charge + zion(ity)
           dipr(1) = dipr(1) + zion(ity)*clust%xatm(j)
           dipr(2) = dipr(2) + zion(ity)*clust%yatm(j)
           dipr(3) = dipr(3) + zion(ity)*clust%zatm(j)
        enddo
     enddo

     ! Convert dipole units from atomic ones to debyes.
     dipr = dipr*debye 

     elec_st%dipx = dipr(1)
     elec_st%dipy = dipr(2)
     elec_st%dipz = dipr(3)
     elec_st%charge = charge

     ! Report results.
     write(7,50) elec_st%ncharge
     write(7,51) charge 
     write(7,*) 'Dipole components [Debye]:'
     write(7,52) elec_st%dipx
     write(7,53) elec_st%dipy 
     write(7,54) elec_st%dipz 
     call myflush(7)

50   format(' Tot. charge (user input) [e] = ',f12.8)
51   format(' Tot. charge (computed) [e] = ',f12.8)
52   format(' x-component = ',f12.8)
53   format(' y-component = ',f12.8)
54   format(' z-component = ',f12.8)
  endif

end subroutine dipole
!===============================================================
