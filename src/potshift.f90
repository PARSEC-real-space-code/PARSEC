!===============================================================
!
! Copyright (C) 2005 Finite Difference Research Group
! This file is part of parsec, http://www.ices.utexas.edu/parsec/
!
! This subroutine shifts a potential so that its average is zero.
! Needed only with periodic boundary conditions.
!
! author: ?
!
!---------------------------------------------------------------
subroutine potshift(vcell,hcub,nrep,ndim,procs_num,comm,pot,avg)

  use constants
  implicit none
  !
  ! Input/Output variables:
  !
  ! hcub = (grid spacing)^3
  real(dp), intent(in) :: hcub
  ! number of symmetry representations (used as ratio between
  ! volumes of full cell and of irreducible wedge)
  integer, intent(in) :: nrep
  ! size of input potential
  integer, intent(in) :: ndim
  ! communication parameters
  integer, intent(in) :: procs_num,comm
  ! input: potential, output: shifted potential
  real(dp), intent(inout) :: pot(ndim)
  ! volume of unit cell
  real(dp), intent(in) :: vcell
  ! average of potential prior to shifting
  real(dp), intent(out) :: avg
  real(dp), dimension(1):: avgvec
  !---------------------------------------------------------------
  avg = sum(pot) * real(nrep,dp)
  avgvec = avg
  call psum(avgvec,1,procs_num,comm)
  avg = avgvec(1)

  ! Rescale average.
  avg = avg*hcub / vcell

  pot = pot - avg

end subroutine potshift
!===============================================================
