!===============================================================
!
!     Modified by Alberto Garcia (March 1991) to implement
!     the method in Gilat and Bharatiya, PRB 12, 3479 (1975).
!
!     1996 additional documentation. Paul Delaney
!
!     1996 removed compton-profile specific comments and features,
!     and assimilated it into paratec. Resistance is futile...
!     Bernd Pfrommer
!
!		2007 imported to Parsec, from Paratec for DOS computations.
!		Or Cohen
!---------------------------------------------------------------
subroutine spec_tet(vol, e, f, ngrid, ebound, edos, fgrid, neqtet)
  !
  use constants
  implicit none    ! nope.
  !
  !     INPUT:
  !     -----
  !
  integer, intent(in) :: ngrid      ! number of grid points on which energy
                                    ! is stored
  real(dp), intent(in) :: &
       ebound(2), &                 ! start(1) and end(2) energies
                                    ! will be computed.
       f(4), &                      ! the spectral function at the corners
                                    ! of the tetrahedron
       e(4), &                      ! energy levels at the corner of the
                                    ! tetrahedron
       vol                          ! the volume of the tetrahedron in k-space
  !
  !     INPUT/OUTPUT: (NOTE: THE CONTRIBUTION FROM THIS TETRAHEDRON IS ADD
  !     ------------
  !
  integer, intent(inout) :: &
       neqtet                       ! counter for tetrahedra with equal corners
  real(dp), intent(inout) :: &
       edos(ngrid), &               ! the contribution to DOS
       fgrid(ngrid)                 ! the contribution to DOS*F(E)
  !
  !     DESCRIPTION:
  !     -----------
  !
  !     spec_tet calculates the contribution to the spectral function
  !
  !
  !     and the DOS from one tetrahedron. This is one for energies
  !     defined by ebound.
  !
  !
  !     ------------    local variables -----------------------
  !
  real(dp), parameter :: tolerance = 1.0d-9
  real(dp), parameter :: delta = 1.0d-8
  real(dp) :: d21, d31, d32, d41, d42, d43, f21, f31, f32, f41, &
       f42, f43, dn1, dn2, dn4, ei, x, s, fm, fms, denerg, estart
  integer :: i, j, ind(4), istart, iend  
  real(dp) :: es(4), fs(4)  
  real(dp) :: xx, yy, eq_tol  
  logical :: equal  

  equal(xx, yy, eq_tol) = (abs(xx - yy) <= eq_tol)  
  !
  !     --------------------------------------------------------
  !
  !     sort energy levels at corners in increasing order
  call sort(4, e, ind)  
  !
  es(1) = e(ind(1))  
  es(2) = e(ind(2))  
  es(3) = e(ind(3))  
  es(4) = e(ind(4))  
  !
  fs(1) = f(ind(1))  
  fs(2) = f(ind(2))  
  fs(3) = f(ind(3))  
  fs(4) = f(ind(4))
  !
  !     if energies slightly different set them equal.
  !     this is important for numerical stability
  !
  if (equal(es(2), es(1), tolerance)) es(2) = es(1)  
  if (equal(es(3), es(2), tolerance)) es(3) = es(2)  
  if (equal(es(4), es(3), tolerance)) es(4) = es(3)  

  if (es(1) == es(2) .and. es(2) == es(3) .and. es(3) == es(4)) then
     ! found tetrahedron with four corners the same
     neqtet = neqtet + 1  
     return  
  end if

  !     find energy and matrix element differences

  d21 = es(2) - es(1)  
  d31 = es(3) - es(1)  
  d32 = es(3) - es(2)  
  d41 = es(4) - es(1)  
  d42 = es(4) - es(2)  
  d43 = es(4) - es(3)
  
  f21 = fs(2) - fs(1)  
  f31 = fs(3) - fs(1)  
  f32 = fs(3) - fs(2)  
  f41 = fs(4) - fs(1)  
  f42 = fs(4) - fs(2)  
  f43 = fs(4) - fs(3)
  
  !     compute denominators
  !     Hold on computing Gilat's eps's, since we might divide by
  !     zero...

  dn1 = (d41 * d31 * d21) / vol  
  dn2 = (d42 * d32 * d21) / vol  
  dn4 = (d43 * d42 * d41) / vol

  !     Note that these D's have the volume built into them.

  denerg = abs(ebound(2) - ebound(1)) / real(ngrid, dp)  

  estart = ebound(1)  
  !      write(9,*) 'corner values:',es
  !
  !     find start and end value for better speed
  !
  istart = max(1, int((es(1) - ebound(1)) / denerg) - 2)  
  iend = min(ngrid - int((ebound(2) - es(4)) / denerg) + 2, ngrid)

  do i = istart, iend                     !  loop over energy grid
     ei = estart + denerg * real(i, dp)  
     do j = 1, 4  
        if (equal(es(j), ei, delta)) then  
           if (ei < es(j)) then  
              ei = ei + 2.02d0 * delta  
           else  
              ei = ei + 1.01d0 * delta  
           end if
        end if
     end do

     x = ei - es(1)  
     !
     !     distinguish the different energy regions
     !
     !   -------------     E <= E1 : do noth
     if (ei <= es(1)) then  
        s = zero
        fms = zero
        !            write(9,*) 'case1:', i, ei, s
        ! ----------    E1 < E <= E2

     else if (ei <= es(2)) then  
        s = three * x * x / dn1  
        fm = (f21 / d21 + f31 / d31 + f41 / d41) * x * third
        fms = fm * s  
        !            write(9,*) 'case2:', i, ei, s
        ! ----------    E2 < E <= E3

     else if (ei <= es(3)) then  
        if (es(1) /= es(2)) then  
           s = three / dn1 * x * x - three / dn2 * (x - d21)**2  
           fms = (f21 / d21 + f31 / d31 + f41 / d41) * x**3 / dn1 - &
                (f21 / d21 + f32 / d32 + f42 / d42) * (x - d21)**3 / dn2 - &
                three * f21 * (x - d21)**2 / dn2
        else  

           s = three * vol * x / (d41 * d31) * (two - x * (one / d31 + &
                one / d41))
           fms = vol * x / (d41 * d31) * (three * f21 + three * x * &
                (f32 / d32 + f42 / d42) - x**2 / (d41 * d31) * (d41 * &
                (f31 / d31 + f32 / d32 + f41 / d41) + d31 * (f31 / d31 + &
                f41 / d41 + f42 / d42) - f21))
        end if
        !            write(9,*) 'case3:', i, ei, s
        ! -----------  E3 < E <= E4

     else if (ei <= es(4)) then  
        s = three * (x - d41)**2 / dn4  
        fm = (f41 / d41 + f42 / d42 + f43 / d43) * (x - d41) * third + &
             f41
        fms = fm * s  
        !            write(9,*) 'case4:', i, ei, s
        !----   E > E4 : increment integrated DO
     else  
        s = zero
        fms = zero
        !            write(9,*) 'case5:', i, ei, s

     end if

     if (fs(1) * s + fms < zero) then
        write(9,*) 'I ran across a negative value in the calculation'
        write(9,*) 'of the angdos.  You may need to make "tolerance"'
        write(9,*) 'smaller in the file spec_tet.f90 (see below)'
        write(9,*) 'e2 - e1 is only', es(2)-es(1)
     end if

     edos(i) = edos(i) + s  
     fgrid(i) = fgrid(i) + (fs(1) * s + fms)  

  end do

  return  

end subroutine spec_tet
