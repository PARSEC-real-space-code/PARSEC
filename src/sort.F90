!===============================================================
!
!     SORTS AN ARRAY BY THE HEAPSORT METHOD
!     W. H. PREUSS ET AL. NUMERICAL RECIPES
!     
!     2007 - Imported from Paratec
!---------------------------------------------------------------!*
subroutine sort(n, arrin, indx)

  use constants
  implicit none
  !
  integer, intent(in) :: n
  real(dp), intent(in) :: arrin(n)
  integer, intent(out) :: indx(n)
  !
  integer :: i, j, l, ir, indxt
  real(dp) :: q

  if (n == 0) return  
  if (n == 1) then  
     indx(1) = 1  
     return  
  end if

  do j = 1, n  
     indx(j) = j  
  end do

  l = n / 2 + 1  
  ir = n

  do
     if (l > 1) then  
        l = l - 1  
        indxt = indx(l)  
        q = arrin(indxt)  
     else  
        indxt = indx(ir)  
        q = arrin(indxt)  
        indx(ir) = indx(1)  
        ir = ir - 1  
        if (ir == 1) then  
           indx(1) = indxt
           return  
        end if
     end if
     i = l  
     j = l + l  
     do while (j <= ir)
        if (j < ir) then  
           if (arrin(indx(j)) < arrin(indx(j + 1))) j = j + 1  
        end if
        if (q < arrin(indx(j))) then  
           indx(i) = indx(j)  
           i = j 
           j = j + j  
        else  
           j = ir + 1  
        end if
     end do
     indx(i) = indxt
  end do

end subroutine sort
!*
!*
subroutine iheapsort(n, iarrin, indx)
  !
  !     SORTS AN ARRAY BY THE HEAPSORT METHOD
  !     W. H. PREUSS ET AL. NUMERICAL RECIPES
  !
  use constants
  implicit none
  !
  integer, intent(in) :: n, iarrin(n)
  integer, intent(out) :: indx(n)
  !
  integer :: i, j, l, indxt, ir, iq

  if (n == 0) return  
  if (n == 1) then  
     indx(1) = 1  
     return  
  end if

  do j = 1, n  
     indx(j) = j  
  end do

  l = n / 2 + 1  
  ir = n

  do
     if (l > 1) then  
        l = l - 1  
        indxt = indx(l)  
        iq = iarrin(indxt)  
     else  
        indxt = indx(ir)  
        iq = iarrin(indxt)  
        indx(ir) = indx(1)  
        ir = ir - 1  
        if (ir == 1) then  
           indx(1) = indxt
           return  
        end if
     end if
     i = l
     j = l + l
     do while (j <= ir)
        if (j < ir) then  
           if (iarrin(indx(j)) < iarrin(indx(j + 1))) j = j + 1
        end if
        if (iq < iarrin(indx(j))) then  
           indx(i) = indx(j)  
           i = j
           j = j + j  
        else  
           j = ir + 1  
        end if
     end do
     indx(i) = indxt
  end do

end subroutine iheapsort
