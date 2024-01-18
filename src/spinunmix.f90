!===============================================================
!
! Copyright (C) 2005 Finite Difference Research Group
! This file is part of parsec, http://www.ices.utexas.edu/parsec/
!
! author: ?
!
!---------------------------------------------------------------
subroutine spinunmix(pot,mixer)

  use constants
  use potential_module
  use mixer_module
  implicit none
  !
  ! Input/Output variables:
  !
  type (potential), intent(inout) :: pot
  type (mixer_data), intent(in) :: mixer
  !
  ! Work variables:
  !
  integer i
  !---------------------------------------------------------------
  do i = 1,pot%ndim
     pot%vnew(i,1)=mixer%xout(i)
     pot%vnew(i,2)=mixer%xout(pot%ndim+i)
  enddo

end subroutine spinunmix
!===============================================================
