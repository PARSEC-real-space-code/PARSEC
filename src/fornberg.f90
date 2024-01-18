!===============================================================
!
!  Copyright (C) 2005 Finite Difference Research Group
!  This file is part of parsec, http://www.ices.utexas.edu/parsec/
!
!  Coefficients for the first & second order numerical derivative
!  under the centered finite difference scheme.
!  Bengt Fornberg,  Exxon Res. & Eng. Co., NJ 08801, 'bfornbe@erenj.com'
!  David M. Sloan,  Dept. of Mathematics, U. of Strathclyde,
!  Glasgow G1 1HX, Scotland,  'd.sloan@strath.ac.uk'
!  Acta Numerica 94,  Cambridge Univ. Press (1994)
!
!---------------------------------------------------------------
subroutine fornberg(iddd,norder,coe,ierr)

  use constants
  implicit none
  !
  !  Input/Output variables:
  !
  !  order of expansion of derivative. 
  !  it is the number of neighbors used ON ONE SIDE.
  !  the maximum order implemented is 20.
  integer, intent(in) :: norder

  !  iddd - order of the derivative (iddd = 1 or 2)
  integer, intent(in) :: iddd

  !  coe - coefficients for the derivative
  real(dp), intent(out) :: coe(-norder:norder)

  !  ierr - error flag, 250 < ierr < 261
  integer, intent(out) :: ierr
  !  
  !  Work variables:
  !
  !  counters 
  integer i
  !---------------------------------------------------------------
  !
  !  First order derivative
  !
  if(iddd == 1) then

     select case (norder)
     case (1)
        coe(1) =  0.50000000000000D+00
     case (2)
        coe(1) =  2.d0/3.d0
        coe(2) = -1.d0/12.d0
     case (3)
        coe(1) =  3.d0/4.d0
        coe(2) = -3.d0/20.d0
        coe(3) =  1.d0/60.d0
     case (4)   
        coe(1) =  4.d0/5.d0
        coe(2) = -1.d0/5.d0
        coe(3) =  4.d0/105.d0
        coe(4) = -1.d0/280.d0
     case (5)     
        coe(1) =  0.8333333333D+00
        coe(2) = -0.2380952381D+00
        coe(3) =  0.5952380952D-01
        coe(4) = -0.9920634921D-02
        coe(5) =  0.7936507937D-03
     case (6)    
        coe(1) =  0.8571428571D+00
        coe(2) = -0.2678571429D+00
        coe(3) =  0.7936507937D-01
        coe(4) = -0.1785714286D-01
        coe(5) =  0.2597402597D-02
        coe(6) = -0.1803751804D-03
     case (7)    
        coe(1) =  0.8750000000D+00
        coe(2) = -0.2916666667D+00
        coe(3) =  0.9722222222D-01
        coe(4) = -0.2651515152D-01
        coe(5) =  0.5303030303D-02
        coe(6) = -0.6798756799D-03
        coe(7) =  0.4162504163D-04
     case (8)    
        coe(1) =  0.8888888889D+00
        coe(2) = -0.3111111111D+00
        coe(3) =  0.1131313131D+00
        coe(4) = -0.3535353535D-01
        coe(5) =  0.8702408702D-02
        coe(6) = -0.1554001554D-02
        coe(7) =  0.1776001776D-03
        coe(8) = -0.9712509713D-05
     case (9)      
        coe(1) =  0.9000000000D+00
        coe(2) = -0.3272727273D+00
        coe(3) =  0.1272727273D+00
        coe(4) = -0.4405594406D-01
        coe(5) =  0.1258741259D-01
        coe(6) = -0.2797202797D-02
        coe(7) =  0.4495504496D-03
        coe(8) = -0.4627725216D-04
        coe(9) =  0.2285296403D-05
     case (10)    
        coe(1) =  0.9090909091D+00
        coe(2) = -0.3409090909D+00
        coe(3) =  0.1398601399D+00
        coe(4) = -0.5244755245D-01
        coe(5) =  0.1678321678D-01
        coe(6) = -0.4370629371D-02
        coe(7) =  0.8814714697D-03
        coe(8) = -0.1285479227D-03
        coe(9) =  0.1202787580D-04
        coe(10)= -0.5412544112D-06
     end select

     coe(0) = 0.d0
     do i = 1,norder
        coe(-i) = -coe(i)
     enddo
     !
     !  Second order derivative
     !
  else if (iddd == 2) then

     select case (norder)
     case (1)
        coe(0) = -0.20000000000000D+01
        coe(1) =  0.10000000000000D+01
     case (2) 
        coe(0) = -0.25000000000000D+01
        coe(1) =  0.13333333333333D+01
        coe(2) = -0.83333333333333D-01
     case (3)
        coe(0) = -0.27222222222222D+01
        coe(1) =  0.15000000000000D+01
        coe(2) = -0.15000000000000D+00
        coe(3) =  0.11111111111111D-01
     case (4)
        coe(0) = -0.28472222222222D+01
        coe(1) =  0.16000000000000D+01
        coe(2) = -0.20000000000000D+00
        coe(3) =  0.25396825396825D-01
        coe(4) = -0.17857142857143D-02
     case (5)
        coe(0) = -0.29272222222222D+01
        coe(1) =  0.16666666666667D+01
        coe(2) = -0.23809523809524D+00
        coe(3) =  0.39682539682540D-01
        coe(4) = -0.49603174603175D-02
        coe(5) =  0.31746031746032D-03
     case (6)
        coe(0) = -0.29827777777778D+01
        coe(1) =  0.17142857142857D+01
        coe(2) = -0.26785714285714D+00
        coe(3) =  0.52910052910053D-01
        coe(4) = -0.89285714285714D-02
        coe(5) =  0.10389610389610D-02
        coe(6) = -0.60125060125060D-04
     case (7)
        coe(0) = -0.30235941043084D+01
        coe(1) =  0.17500000000000D+01
        coe(2) = -0.29166666666667D+00
        coe(3) =  0.64814814814815D-01
        coe(4) = -0.13257575757576D-01
        coe(5) =  0.21212121212121D-02
        coe(6) = -0.22662522662523D-03
        coe(7) =  0.11892869035726D-04
     case (8)
        coe(0) = -0.30548441043084D+01
        coe(1) =  0.17777777777778D+01
        coe(2) = -0.31111111111111D+00
        coe(3) =  0.75420875420875D-01
        coe(4) = -0.17676767676768D-01
        coe(5) =  0.34809634809635D-02
        coe(6) = -0.51800051800052D-03
        coe(7) =  0.50742907885765D-04
        coe(8) = -0.24281274281274D-05
     case (9)
        coe(0) = -0.30795354623331D+01
        coe(1) =  0.18000000000000D+01
        coe(2) = -0.32727272727273D+00
        coe(3) =  0.84848484848485D-01
        coe(4) = -0.22027972027972D-01
        coe(5) =  0.50349650349650D-02
        coe(6) = -0.93240093240093D-03
        coe(7) =  0.12844298558584D-03
        coe(8) = -0.11569313039901D-04
        coe(9) =  0.50784364509855D-06
     case (10)
        coe(0) = -0.30995354623331D+01
        coe(1) =  0.18181818181818D+01
        coe(2) = -0.34090909090909D+00
        coe(3) =  0.93240093240093D-01
        coe(4) = -0.26223776223776D-01
        coe(5) =  0.67132867132867D-02
        coe(6) = -0.14568764568765D-02
        coe(7) =  0.25184899134479D-03
        coe(8) = -0.32136980666392D-04
        coe(9) =  0.26728612899924D-05
        coe(10)= -0.10825088224469D-06
     end select

     do i = 1,norder
        coe(-i) = coe(i)
     end do

  else
     write(9,*) ' ERROR: invalid derivative order, iddd = ',iddd
     write(9,*) ' STOP in FORNBERG '
     ierr = 251 
     return
  endif

end subroutine fornberg
!===============================================================
