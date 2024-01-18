!===============================================================
!
! Copyright (C) 2005 Finite Difference Research Group
! This file is part of parsec, http://www.ices.utexas.edu/parsec/
!
! This is a machine-dependent subroutine that should be modified
! when porting to new systems.
!
!---------------------------------------------------------------
subroutine myflush(file_num)

  use constants
  implicit none
  !
  ! Input/Output variables:
  !
  integer file_num 
#ifdef IBM
  external flush
#elif GFORTRAN
  intrinsic flush
#else
  !let's use fortran 2003 flush.
#endif
  !---------------------------------------------------------------
#ifdef IBM
  
  call flush(file_num)
  
#elif GFORTRAN
  call flush(file_num)
#else
  
  flush(file_num)
  
#endif
end subroutine myflush
!===============================================================
