!=====================================================================
!
!  Cubic spline interpolation
!
!  (C) Copr. 1986-92 Numerical Recipes Software #0).
!
!---------------------------------------------------------------------
subroutine spline(x,y,n,yp1,ypn,y2) 

  use constants
  implicit none
  !
  !  Input/Output variables:
  !
  !  number of tabulated data points
  integer, intent(in) :: n
  !  the grid points and their tabulated function values
  real(dp), intent(in) :: x(n),y(n)
  !  specified first derivatives as the boundary conditions
  real(dp), intent(in) :: yp1,ypn
  !  the second derivative at the grid points
  real(dp), intent(out) :: y2(n) 
  !
  !  Work variables:
  !
  !  variables for tridiagonal algorithm
  real(dp) :: p,qn,sig,un,u(n) 
  !  counters
  integer :: i,k 

  if (yp1.gt.0.99e10) then 
  ! The lower boundary condition is set either to be "natural"
    y2(1)= zero
    u(1) = zero
  else
  ! or else to have a specified first derivative. 
    y2(1) = -half
    u(1) = (three/(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1) 
  endif 

  do i=2,n-1 
  ! This is the decomposition loop of the tridiagonal algorithm.
  ! y2 and u are used for temporary storage of the decomposed factors. 

    sig = (x(i)-x(i-1))/(x(i+1)-x(i-1)) 
    p   = sig*y2(i-1) + two 
    y2(i)= (sig-one)/p 
    u(i) = (six*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1)) &
                /(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p 
  enddo

  if (ypn.gt.0.99e10) then
  ! The upper boundary condition is set either to be "natural"
    qn = zero
    un = zero
  else
  ! or else to have a specified first derivative. 
    qn = half
    un = (three/(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1))) 
  endif 

  y2(n) = (un-qn*u(n-1))/(qn*y2(n-1)+one) 
  do k = n-1,1,-1
  ! This is the backsubstitution loop of the tridiagonal algorithm. 
    y2(k) = y2(k)*y2(k+1)+u(k) 
  enddo
  return 
end subroutine spline
!=====================================================================
!
!  Cubic spline interpolation
!
!  (C) Copr. 1986-92 Numerical Recipes Software #0).
!
!---------------------------------------------------------------------
subroutine splint(xa,ya,y2a,n,x,y,y1) 

  use constants
  implicit none
  !
  !  Input/Output variables:
  !
  !  number of tabulated data points
  integer, intent(in) :: n
  !  the grid points and their tabulated function values 
  !  and second derivative
  real(dp), intent(in) :: xa(n),ya(n),y2a(n)
  !  the position at which the function value will be interpolated
  real(dp), intent(in) :: x
  !  desired function value and its first derviative at x
  real(dp), intent(out) :: y,y1
  !
  !  Work variables:
  !
  !  variables for interpolation
  integer :: k,khi,klo 
  real(dp) a,b,h,dy
     
  klo = 1 
  khi = n

  do while (khi-klo.gt.1)
     k = (khi+klo)/2
     if(xa(k).gt.x)then 
       khi = k 
     else 
       klo = k 
     endif 
  enddo

  ! klo and khi now bracket the input value of x. 
  h = xa(khi)-xa(klo) 
  a = (xa(khi)-x)/h
  dy = ya(khi)-ya(klo)

  ! Cubic spline polynomial is now evaluated. 
  b = (x-xa(klo))/h 
  y = a*ya(klo)+b*ya(khi) + ((a**3-a)*y2a(klo) + &
    (b**3-b)*y2a(khi))*(h**2)/six 
  y1 = dy/h - (three*a**2-one)/six*h*y2a(klo) + &
    (three*b**2-one)/six*h*y2a(khi)

  return 
end subroutine splint
!=====================================================================
