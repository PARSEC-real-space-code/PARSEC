!===============================================================
!
! Copyright (C) 2005 Finite Difference Research Group
! This file is part of parsec, http://www.ices.utexas.edu/parsec/
!
! Author: M. Alemany; interface with FFTW 3.* added by M. Tiago
!
! This file actually has a set of two subroutines, both called fftw:
! cpp option FFTW2 : interface with FFTW 2.1.* (and possibly older versions)
! cpp option FFTW3 : interface with FFTW 3.*
!
! AJB: it seems to me that parsec should be using real-2-complex fftws
!
!---------------------------------------------------------------
subroutine cfftw(chd,n1,n2,n3,mode)

  use constants
  implicit none
  !
  ! Input/Output variables:
  !
  integer, intent(in) :: n1, n2, n3, mode
  complex(dpc), intent(inout) :: chd(n1*n2*n3)
  integer :: k,j,i
  !
  ! Work variables:
  !
  integer dirn, len
  integer (kind=8) plan
  complex(dpc) :: fac
!  character (len=10) :: file_name
! AJB: this checks out:
  integer, parameter :: fftw_forward = -1
  integer, parameter :: fftw_backward = +1

#ifdef USEFFTW2
  !
  ! Interface with FFTW2 library, fftw_lib routines:
  !
  integer, parameter :: fftw_estimate = 0
  integer, parameter :: fftw_measure = 1
  integer, parameter :: fftw_in_place = 8
  integer, parameter :: fftw_use_wisdom = 16
  integer, parameter :: real_to_complex = -1
  integer, parameter :: complex_to_real = 1
#elif USEFFTW3
  !
  ! Interface with FFTW3 library, parameters from fftw3.f
  !
  integer, parameter :: fftw_estimate = 64
#endif
  !---------------------------------------------------------------
  len = n1*n2*n3 

  if (mode > 0) then
     !file_name = 'forwardfft'
     dirn = fftw_forward
  else
     !file_name = 'backwrdfft'
     dirn = fftw_backward
  endif

#ifdef USEFFTW3

  call dfftw_plan_dft_3d(plan, n1, n2, n3,chd,chd,dirn,fftw_estimate)
  call dfftw_execute(plan)
  call dfftw_destroy_plan(plan)

#elif USEFFTW2

!#ifdef SUN
!  call fftw3d_f77_create_plan_(plan,n1,n2,n3,dirn,fftw_estimate + &
!#else
  call fftw3d_f77_create_plan(plan,n1,n2,n3,dirn,fftw_estimate + fftw_in_place)
!#endif

! #ifdef SUN
!   call fftwnd_f77_(plan,1,chd,1,0,0,0,0)
! #else
   call fftwnd_f77(plan,1,chd,1,0,0,0,0)
! #endif
#endif

 !in case we asked for "forward" also normalize:
  if (mode > 0) then
     fac = zone / real(len,dp)
     call zscal(len,fac,chd,1)
  endif

end subroutine cfftw
!===============================================================
