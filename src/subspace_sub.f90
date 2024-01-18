!===================================================================
!
!  Non-recursive implementation of the quick sort method.
!  always sort the array such that the required eigenvalues
!  are located at the beginning of the sorted array.
!
!-------------------------------------------------------------------
subroutine sortreal(which, n, x1, x2, targetval)

  use parsec_global_data
  implicit none
  !
  !  Input/Output variables:
  !
  character(len=2), intent(in) :: which
  integer, intent(in) :: n
  real(dp), intent(inout) :: x1(0:n-1)  ! the array to be sorted
  integer, intent(inout) :: x2(0:n-1)   ! return the index of sorted x1
  real(dp),intent(in) :: targetval
  !
  !  Work variables:
  !
  integer :: i, igap, j
  real(dp) :: temp, xtmp(0:n-1)
  integer :: itemp
  integer :: ierr

  !-----------------------------------------------------------------
  igap = n / 2
  if (igap > 0) call which_xr(which, n, x1, xtmp, targetval,ierr)

  do while (igap > 0) 
              
     do i = igap, n-1
        j = i-igap

        do while (j >= 0) 
           !  always sort the array in non-decreasing order 
           !  according to xtmp, the required eigenvalues are 
           !  located at the beginning of the sorted array x1.
           if (xtmp(j) > xtmp(j+igap)) then
              !  need to reorder xtmp since the next comparisons
              !  are based on the changing x1,
              !  and xtmp need to reflect this change in x1
              temp = xtmp(j)
              xtmp(j) = xtmp(j+igap)
              xtmp(j+igap) = temp

              !  reorder x1
              temp = x1(j)
              x1(j) = x1(j+igap)
              x1(j+igap) = temp

              !  reorder x2
              itemp = x2(j)
              x2(j) = x2(j+igap)
              x2(j+igap) = itemp
           else
              exit
           end if
           j = j-igap
        end do
     end do
     igap = igap / 2
  end do

end subroutine sortreal

!===================================================================
subroutine  which_xr(which, n, x, xtmp, targetval, ierr)

  use parsec_global_data
  implicit none
  !
  !  Input/Output variables:
  !
  character(len=2), intent(in) :: which
  integer, intent(in) :: n
  real(dp), intent(in) :: x(1:n)  
  real(dp), intent(out) :: xtmp(1:n)
  real(dp), intent(in) :: targetval
  integer, intent(out) :: ierr
  !
  !  Work variables:
  !
  integer :: i

  !-----------------------------------------------------------------

  ierr = 0
  if (which == 'SA'  .or.  which == 'SR') then
     do i = 1, n
        xtmp(i) =  x(i)
     end do
  elseif (which == 'LA'  .or.  which == 'LR') then
     do i = 1, n
        xtmp(i) = - x(i)
     end do
  elseif (which == 'SM') then
     do i = 1, n
        xtmp(i) = abs(x(i))
     end do
  elseif (which == 'LM') then
     do i = 1, n
        xtmp(i) = - abs(x(i))
     end do
  elseif (which == 'TA') then
     do i = 1, n
        xtmp(i) =  abs(x(i)-targetval)
     end do
  elseif (which == 'BE') then
     do i = 1, n
        xtmp(i) =  x(i)
     end do
  else
     write(9,*) "Invalid type of which in: subroutine which_xr"
     ierr = 1
  end if

end subroutine which_xr
!=====================================================================
!
! Printout routine.
!
!-------------------------------------------------------------------
subroutine print_array(ldn, m, n, A)

  use constants
  implicit none
  !
  !  Input/Output variables:
  !
  integer, intent(in) :: ldn, m, n
  real(dp), intent(in) :: A(ldn,n)
  !
  !  Work variables:
  !
  integer  :: irow, jcol

  if (ldn < m) then
     write(9,*) "Error in print_array: array ldn < m"
     return
  endif
  if (m /= 1 .and. n /= 1) then
     do irow = 1, m
        write(9,711) (A(irow, jcol), jcol=1,n)
     enddo
  elseif (n == 1) then
     write(9,711) (A(irow,1), irow=1,m)
  else 
     write(9,711) (A(1,jcol), jcol=1,n)
  endif

711 format(5(e16.9, 2x))

end subroutine print_array

