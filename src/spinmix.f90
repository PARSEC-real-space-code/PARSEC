!===============================================================
!
! Copyright (C) 2005 Finite Difference Research Group
! This file is part of parsec, http://www.ices.utexas.edu/parsec/
!
! author: ?
!
!---------------------------------------------------------------
subroutine spinmix(pot,mixer)

  use constants
  use potential_module
  use mixer_module
  implicit none
  !
  ! Input/Output variables:
  !
  type (potential), intent(in) :: pot
  type (mixer_data), intent(inout) :: mixer
  !
  ! Work variables:
  !
  integer i
  !---------------------------------------------------------------
  do i = 1,pot%ndim
     mixer%xin(i)=pot%vold(i,1)
     mixer%xin(pot%ndim+i)=pot%vold(i,2)
     mixer%xout(i)=pot%vnew(i,1)
     mixer%xout(pot%ndim+i)=pot%vnew(i,2)
  enddo

end subroutine spinmix
!===============================================================
