!===============================================================
!
!  Copyright (C) 2005 Finite Difference Research Group
!  This file is part of parsec, http://www.ices.utexas.edu/parsec/
!
!  FOR CONFINED SYSTEMS:
!  This subroutine computes an effective charge density that is
!  used as the source term in hpotcg for determining the Hartree
!  potential. It is based on a mutipole expansion of the boundary
!  potential in terms of the charge density, providing the
!  appropriate boundary condition. Finite difference terms
!  involving points outside the sphere are computed using the
!  boundary potential and added to the effective charge density.
!
!  For a general discussion of the theory of mutipole expansions,
!  see, e.g., Landau and Lifshitz, the classical theory of fields,
!  pp. 97-99.
!
!  This subroutine uses Eqs. (41.12) and (41.13) of Landau and
!  Lifshitz. However, the spherical harmonics and their properties
!  follow the definition given by Jackson, Classical
!  Electrodynamics, pp. 98-99.
!
!---------------------------------------------------------------
subroutine hartset(grid,rsymm,parallel,norder,lpole,coeke,rho,clm,brho)

  use constants
  use grid_module
  use symmetry_module
  use parallel_data_module
  implicit none
  !
  !  Input/Output variables:
  !
  !  grid related data
  type (grid_data), intent(in) :: grid
  !  symmetry operations in reduced group:
  type (symmetry), intent(in) :: rsymm
  !  parallel computation related data
  type (parallel_data), intent(in) :: parallel

  !  number of neighbors used for derivation (on one side)
  integer, intent(in) :: norder
  !  order of multipole expansion 
  integer, intent(in) :: lpole
  !  coefficients for taking a second derivative
  real(dp), intent(in) :: coeke(-norder:norder,3)
  !  charge density on the 3-d grid
  real(dp), intent(in) :: rho(parallel%ldn)
  !  coefficient array: clm = (l-m)!/(l+m)!
  real(dp), intent(in) :: clm(0:lpole, 0:lpole)
  !  effective charge density = original charge density augmented by
  !  appropriate boundary terms, to be used as the source term in hpotcg
  real(dp), intent(out) :: brho(parallel%ldn)
  !
  !  Work variables:
  !
  !  offset for each processor
  integer ioffset
  !  grid spacing (atomic units)
  real(dp) :: h
  !  temporary storage variables
  !  counters of total and z-axis angular momentum
  integer l,m

  real(dp) :: &     ! for each grid point labeled by j,
       rrho, &               ! rrho = rho(j)
       rw(3),rrw(3), &       ! grid point equiv. to (xx(j),yy(j),zz(j))
       xx2, &                ! xx2 = xx(j)^2
       yy2, &                ! yy2 = yy(j)^2
       zz2, &                ! zz2 = zz(j)^2
       xy, &                 ! xy = sqrt(xx(j)^2 + yy(j)^2)
       r2, &                 ! r2 = xx(j)^2 + yy(j)^2 + zz(j)^2
       x, &                  ! x = cos(theta_j)
       y                     !  y = sin(theta_j)

  real(dp) :: &      ! for (l,m) angular momentum quantum numbers
       r(0:lpole), &         ! r(1) = sqrt(xx(j)^2 + yy(j)^2 + zz(j)^2)
                             ! r(l) = r(1)^l
       tmp1, &               ! tmp1 = rrho*r(l)
       rinv(0:lpole)         ! rinv(1) = 1/r(1)

  !  For each grid point:
  !  temporary storage of the associated Legendre polynom value
  real(dp) :: plm(0:lpole,0:lpole)
  !  The coefficient defined in Eq. (41,13) of Landau and Lifshitz
  complex(dpc) :: qlm(0:lpole, 0:lpole)

  !  Coefficient array: e^{i m \phi_j} = cos(phi_j) + i*sin(phi_j)
  complex(dpc) :: cxy(lpole)

  !  Contribution to the boundary potential (the contribution is the
  !  real part of this quantity).
  complex(dpc) :: cpole

  real(dp) :: cpole_deb, r12, cpole_deb2
  !  Coordinates of relevant neighboring points outside the boundary
  !  sphere.
  real(dp) :: xh, yh, zh, xn, yn, zn, xh1, yh1, zh1
  !  square of grid spacing
  real(dp) :: hsq
  !  2*h^3
  real(dp) :: hcub2
  !  other counters
  integer i,j,jshell,jrow,itrans, i_deb, j_deb, alcstat
  integer i_n,j_n,k_n, mpinfo, j_tmp
  !  constants
  real(dp), parameter :: pi8 = eight*pi

  real(dp), dimension(:), pointer :: boundary_pot

  !---------------------------------------------------------------

  nullify(boundary_pot)
  !  Assume that grid spacing is the same along all directions.
  h = grid%step(1)

  brho(:) = zero
  hsq = h*h
  hcub2 = two*grid%hcub
  ioffset = parallel%irows(parallel%group_iam) - 1

  !  Zero out all of the qlm entries.
  qlm(:,:) = zzero
  plm(:,:) = zero
  rinv(:) = zero
  !
  !  We need to populate the qlm entries by adding in the
  !  contribution from every grid point that we have data for.
  !
  !  qlm(l,m) is given by
  !
  !  \sum_{j=1}^{ndim}e_j r_j^l P_l^m(cos(\theta_j)) e^{i m \phi_j}
  !
  !  where l = 0 to lpole, m = -l to l, e_j is the given by rho(j),
  !  r_j is the position of the grid point calculated by
  !  sqrt(xx(j)^2 + yy(j)^2 + zz(j)^2), P_l^m is the associate
  !  Legendre polynomial and cos(theta_j) = zz(j)/r and
  !  sin(theta_j) = sqrt(xx(j)^2 + yy(j)^2)/r. We use the fact that
  !  e^{i m \phi_j} = cos(phi_j) + i*sin(phi_j) and we can then
  !  write this term (called cxy(m)) as
  !  ((xx(j) + i*yy(j))/sqrt(xx(j)^2 + yy(j)^2))^m.
  !
  !  We also exploit the fact that there is a simple relationship
  !  between -m and +m terms and we do not calculate any m < 0 
  !  terms. It ends up working out that the contribution to brho
  !  from -m elements is exactly equal to the contribution from +m
  !  terms, so later we just double qlm for m > 0, and eventually
  !  consider only the real part of the final result, and we dont
  !  have to do any other calculations for negative m's. 
  !
  !  This do loop adds in contributions to qlm(l,m) from each grid
  !  point.
  !
  do itrans = 1, rsymm%ntrans
     do j = parallel%irows(parallel%group_iam) &
          ,parallel%irows(parallel%group_iam+1)-1
        rrho = rho(j-ioffset)
        rw(1) = (grid%shift(1) + grid%kx(j)) * h
        rw(2) = (grid%shift(2) + grid%ky(j)) * h
        rw(3) = (grid%shift(3) + grid%kz(j)) * h
        call symop(rsymm,itrans,rw,rrw)
        rw = rrw
        r2 = dot_product(rw,rw)

        xx2 = rw(1) * rw(1)
        yy2 = rw(2) * rw(2)
        zz2 = rw(3) * rw(3)
        !
        !  r(1) = sqrt(xx(j)^2 + yy(j)^2 + zz(j)^2) 
        !
        r(1) = sqrt(r2)
        qlm(0,0) = qlm(0,0) + rrho

        if (r(1) > zero) then
           !
           !  r(l) = r^l (higher powers of r)
           !
           do l = 2, lpole
              r(l) = r(l-1)*r(1)
           enddo

           rinv(1) = one/r(1)
           xy = sqrt(xx2+yy2)
           !
           !  x = cos(theta_j), y = sin(theta_j)
           !
           x = rw(3) * rinv(1)
           y = xy*rinv(1)
           !
           !  Evaluate associate Legendre polynomials for x and y with l from
           !  0 to lpole and m from 0 to l. store each in plm(l,m).
           !
           call setplm(x, y, lpole, plm)
           !
           !  First evaluate the contribution to the m=0 terms of qlm(l,m) as
           !  they are simpler since cxy(0) would be 1. Don't bother
           !  calculating cxy(0) or multiplying it into these terms.
           !
           do l = 1, lpole
              qlm(l,0) = qlm(l,0) + rrho*plm(l,0)*r(l)
           enddo

           if (xy > zero) then
              !
              !  cxy(m) = exp(i*m*phi_j) = (cos(phi_j) + i*sin(phi_j))^m
              !  = ((xx(j) + i*yy(j))/(sqrt(xx(j)^2 + yy(j)^2))^m
              !
              !  First calculate cxy(1).
              !
              cxy(1) = cmplx(rw(1),rw(2),dp) / xy
              !
              !  Now calculate cxy(m) for m=2 to lpole.
              !
              do m = 2, lpole
                 cxy(m) = cxy(m-1)*cxy(1)
              enddo
              !
              !  Now we can calculate the contribution to the qlm(l,m) terms
              !  for m non-zero. tmp1 is used to avoid calculating rrho*r(l)
              !  for each m iteration.
              !
              do l = 1, lpole
                 tmp1 = rrho*r(l)
                 do m = 1, l
                    qlm(l,m) = qlm(l,m) + cxy(m)*plm(l,m)*tmp1
                 enddo
              enddo
           endif
        endif
     enddo                  !  %irows
  enddo                     !  itrans = 1, rsymm%ntrans
  !
  !  Global sum for colecting the qlms from each processor.
  call zpsum(qlm,(lpole+1)*(lpole+1),parallel%group_size,parallel%group_comm)
  !
  !  Once all of the contributions to qlm have been calculated we
  !  need to multiply each term by h^3. for the calculation we are
  !  performing it works out that the contributions from -m terms
  !  are equal to +m terms, so we also double the qlm(l,m) terms
  !  where m > 0. also, for the sake of simplicity we multiply in
  !  the clm(l,m) factor to simplify our work later.
  !
  !  First do m = 0 terms.
  !
  do l = 0, lpole
     qlm(l,0) = qlm(l,0)*clm(l,0)*grid%hcub
  enddo
  !
  !  Then do m > 0 terms (factor of 2 accounts for m < 0 terms).
  !
  do l = 1, lpole
     do m = 1, l, 1
        qlm(l,m) = qlm(l,m)*clm(l,m)*hcub2
     enddo
  enddo

  if(grid%hartree_neibs_flag) then
   write(9,*) 'before boundary pot allocation'
   allocate(boundary_pot(grid%neibs_num),stat=alcstat)
   write(9,*) 'after allocation'
   do i_deb=1, grid%neibs_num
    xh = (grid%shift(1) + grid%neibs_fx(i_deb))*h
    yh = (grid%shift(2) + grid%neibs_fy(i_deb))*h
    zh = (grid%shift(3) + grid%neibs_fz(i_deb))*h
    if(mod(i_deb,1000) .eq. 0) write(9,*) 'done neib: ', i_deb 
    cpole_deb=0
    do j_deb = parallel%irows(parallel%group_iam), &
           parallel%irows(parallel%group_iam+1)-1

           xn = (grid%shift(1) + grid%kx(j_deb))*h
           yn = (grid%shift(2) + grid%ky(j_deb))*h
           zn = (grid%shift(3) + grid%kz(j_deb))*h

           r12=sqrt((xn-xh)**2+(yn-yh)**2+(zn-zh)**2)

           cpole_deb=cpole_deb+rho(j_deb-ioffset)/r12
    enddo

      !cpole_deb = cpole_deb * pi8 * grid%hcub / 4.0d0
    cpole_deb = cpole_deb * grid%hcub
    boundary_pot(i_deb)=cpole_deb
  end do
#ifdef MPI
  write(9,*) 'before psum'
  call psum(boundary_pot, grid%neibs_num, parallel%group_size, &
            parallel%group_comm)
  write(9,*) 'after psum'
#endif
 end if


  !
  !  What we really want to calculate here is the boundary
  !  conditions for the determination of the Hartree contribution.
  !  This is based on a multipole expansion. We start with a
  !  boundary rho (brho) that is 8*pi*rho (4*pi with an extra factor
  !  of two for the Hartree->Rydberg conversion) and then find
  !  points which are OUTSIDE our sphere. Since the potential is not
  !  zero outside of the sphere the contribution of these points to
  !  rho must be included in order to produce an accurate Hartree
  !  potential. The contribution of a given point can be found by
  !  psi(l) = 1/R^{l +1} \sum_{m=-l}^{l}
  !           (\frac{(l-m)!}{l+m)!}*(\sum_{j=1}^{ndim}e_j r_j^l 
  !           P_l^m(cos(\theta_j)) e^{i m \phi_j})*
  !           P_l^m(cos(\theta)) e^{-i m \phi}).
  !  according to what we've done above this simplifies to
  !  psi(l) = 1/R^{l+1} qlm(l,m) P_l^m(cos(\theta)) e^{-i m \phi}.
  !  In this calculation R, theta, and phi are the spherical
  !  coordinates of the point outside the sphere which we are
  !  considering. This calculation looks much like the above, 
  !  however the contribution for l=0,lpole and m=0,l are summed to
  !  give the contribution of the particular point to the grid point
  !  which it neighbors. For simplicity we keep most of the
  !  temporary variables the same and this code is laid out pretty
  !  much identical to the above loop, but looking for any neighbors
  !  of grid points which fall outside of the sphere.
  !
  !  Notice that, since all charge densities have the full
  !  symmetry of the system, there is no need to include characters
  !  in the neighbor grid points (i.e., parallel%tneibs is not used).
  !
  !  Start with the boundary rho (brho) equal to 8*pi*rho.
  !

  do i = parallel%irows(parallel%group_iam), &
       parallel%irows(parallel%group_iam+1)-1
     brho(i-ioffset) = rho(i-ioffset)*pi8
  enddo
  !
  !  Go through each neighbour in all 6 directions
  !  (+/-x, +/-y, +/-z) for the order needed.
  !
  do i = parallel%irows(parallel%group_iam), &
       parallel%irows(parallel%group_iam+1)-1
     do jshell = 1,norder
        do jrow = 1, 6
           !
           !  If this neighbor is outside the sphere then its neibs value
           !  will be greater than ndim. we can then find the position of
           !  this point and calculate its contribution to brho.
           !
           if (parallel%neibs(jrow+(jshell-1)*6,i-ioffset) &
                > parallel%nwedge) then
              xh = (grid%shift(1) + grid%kx(i))*h
              yh = (grid%shift(2) + grid%ky(i))*h
              zh = (grid%shift(3) + grid%kz(i))*h
              i_n=grid%kx(i)
              j_n=grid%ky(i)
              k_n=grid%kz(i)
              !write(7,*) ' i , j , k =',i_n,j_n,k_n
              !write(7,*) 'x,y,z = ',xh,yh,zh
              if (jrow == 1) then
                 xh = xh-h*jshell
                 i_n=i_n-jshell
              elseif (jrow == 2) then
                 xh = xh+h*jshell
                 i_n=i_n+jshell
              elseif (jrow == 3) then
                 yh = yh-h*jshell
                 j_n = j_n-jshell
              elseif (jrow == 4) then
                 yh = yh+h*jshell
                 j_n = j_n+jshell
              elseif (jrow == 5) then
                 zh = zh-h*jshell
                 k_n = k_n-jshell
              elseif (jrow == 6) then
                 zh = zh+h*jshell
                 k_n = k_n + jshell
              endif

              !write(7,*) 'jshell=',jshell,',jrow=',jrow
              !write(7,*) ' i , j , k =',i_n,j_n,k_n
              !write(7,*) 'xh,yh,zh = ',xh,yh,zh
            if(grid%hartree_neibs_flag) then
              j_tmp=grid%neibs_index(i_n,j_n,k_n)
              if(j_tmp>0) then
               xh1=(grid%shift(1)+grid%neibs_fx(j_tmp))*h
               yh1=(grid%shift(2)+grid%neibs_fy(j_tmp))*h
               zh1=(grid%shift(3)+grid%neibs_fz(j_tmp))*h
               !write(7,*) 'xh1,yh1,zh1',xh1,yh1,zh1
              endif
            endif

              xx2 = xh*xh
              yy2 = yh*yh
              zz2 = zh*zh
              r2 = xx2 + yy2 + zz2
              r(1) = sqrt(r2)
              !write(7,*) 'point radius is: ',r(1)
              !
              !  rinv(l) = 1/r^(l+1). calculate 1/r and higher powers
              !
              rinv(0) = one/r(1)
              do l = 1,lpole
                 rinv(l) = rinv(l-1)*rinv(0)
              enddo
              xy = sqrt(xx2+yy2)
              !
              !  x = cos(theta), y = sin(theta)
              !
              x = zh*rinv(0)
              y = xy*rinv(0)
              !
              !  Evaluate associate Legendre polynomials for x and y with l
              !  from 0 to lpole and m from 0 to l. store each in plm(l,m).
              !
              call setplm(x, y, lpole, plm)
              
              !  Initialize cpole to zero.
              !
              cpole = zzero
              !  Calculate the contribution of the m=0 terms.
              !
              do l = 0,lpole
                 cpole = cpole + qlm(l,0)*plm(l,0)*rinv(l)
              enddo

              if (xy > zero) then
                 cxy(1) = cmplx(xh,-yh)
                 cxy(1) = cxy(1)/xy
                 do m=2, lpole, 1
                    cxy(m) = cxy(m-1)*cxy(1)
                 enddo

                 !  Calculate the contribution of the m>0 terms.
                 !
                 do l = 1, lpole
                    do m = 1, l, 1
                       cpole=cpole+qlm(l,m)*plm(l,m)*cxy(m)*rinv(l)
                    enddo
                 enddo
              endif

              !  Add the contributions into brho, but scale by two
              !  (hartree->rydberg) and by coeke.
              !

! DEBUG BOUNDARY CONDITIONS - Recalculating cpole according
!                             to the exact integral of rho/r .

            if(grid%hartree_neibs_flag) then
                write(9,*) ' I am doing hartree_neibs things that are marked as debug'

              j_deb=grid%neibs_index(i_n,j_n,k_n)
              if(j_deb>0) then

                cpole_deb=boundary_pot(j_deb)

              else
 
                write(7,*) 'ERROR found non existing neighbor..'
                write(7,*) 'jshell=',jshell,',jrow=',jrow, &
                           ',i,j,k: ',i_n,j_n,k_n
                write(7,*) 'indexg=',grid%indexg(i_n,j_n,k_n)
                write(7,*) 'neib_index=',grid%neibs_index(i_n,j_n,k_n)
                write(7,*) 'ndim=', grid%ndim,',nwedge=',grid%nwedge
                if(grid%indexg(i_n,j_n,k_n) .ne. grid%ndim) then
                  j_deb=grid%indexg(i_n,j_n,k_n)
                  write(7,*) 'grid f indices', grid%fx(j_deb), &
                             grid%fy(j_deb), grid%fz(j_deb)
                  write(7,*) 'grid k indices', grid%kx(j_deb), &
                             grid%ky(j_deb), grid%kz(j_deb)

              xh = (grid%shift(1) + grid%fx(j_deb))*h
              yh = (grid%shift(2) + grid%fy(j_deb))*h
              zh = (grid%shift(3) + grid%fz(j_deb))*h
                  write(7,*) 'x,y,z = ',xh,yh,zh
                  write(7,*) 'indices: ',grid%fx(j_deb), grid%fy(j_deb) &
                              , grid%fz(j_deb)
              xh = (grid%shift(1) + grid%kx(j_deb))*h
              yh = (grid%shift(1) + grid%ky(j_deb))*h
              zh = (grid%shift(1) + grid%kz(j_deb))*h
                endif
                stop 1

              endif


              brho(i-ioffset) = brho(i-ioffset) - &
                   two*coeke(jshell,1)*cpole_deb

!if(mod(i, 10000) .eq. 0) then
!   write(9,*) " *** ", i
!endif

if(abs(one - cpole_deb/real(cpole,dp)) > 1.0d-5) then
              write(9,*) i, cpole_deb/real(cpole,dp)
endif

            else ! grid%hartree_neibs_flag

              brho(i-ioffset) = brho(i-ioffset) - &
                   two*coeke(jshell,1)*real(cpole,dp)
            endif ! hartree_neibs_flag
              !  WARNING: if the grid spacing h is not the same in the x, y, z
              !  directions then coeke(jshell,1) needs to be changed in a way
              !  that takes into account whether the neighboring point is an
              !  x,y, or z neighbor!!
           endif 
        enddo               ! jrow = 1,6
     enddo                  ! jshell = 1,norder
  enddo
  if(associated(boundary_pot)) deallocate(boundary_pot)
end subroutine hartset
!===============================================================
!
!  This subroutine calculates all of the associated Legendre
!  polynomials up to a supplied lpole (maximum lpole is 9). 
!
!---------------------------------------------------------------
subroutine setplm(x, y, lpole, plm)

  use constants
  implicit none
  !
  !  Input/Output variables:
  !
  !  cos(theta),sin(theta), for a given grid point
  real(dp), intent(in) :: x,y
  !  order of multipole expansion
  integer, intent(in) :: lpole
  !  array containing P_{lm}, the associated Legendre polynomials
  real(dp), intent(out) :: plm(0:lpole, 0:lpole)
  !
  !  Work variables:
  !
  !  powers of x, y: xn=x^n, yn=y^n
  real(dp) :: x2, x3, x4, x5, x6, x7, x8, x9
  real(dp) :: y2, y3, y4, y5, y6, y7, y8, y9
  !---------------------------------------------------------------
  !
  !  Recursive Definitions for P_{lm}
  !
  !  P_{M,M}
  !  
  !  p00(x,y) = ONE
  !  p11(x,y) = -y
  !  p22(x,y) = +y*y*3.0
  !  p33(x,y) = -y*y*y*5.0*3.0
  !  p44(x,y) = +y*y*y*y*7.0*5.0*3.0
  !  p55(x,y) = -y*y*y*y*y*9.0*7.0*5.0*3.0
  !  p66(x,y) = +y*y*y*y*y*y*11.0*9.0*7.0*5.0*3.0
  !  p77(x,y) = -y*y*y*y*y*y*y*13.0*11.0*9.0*7.0*5.0*3.0
  !  p88(x,y) = +y*y*y*y*y*y*y*y*15.0*13.0*11.0*9.0*7.0*5.0*3.0
  !  p99(x,y) = -y*y*y*y*y*y*y*y*y*17.0*15.0*13.0*11.0*9.0*7.0*5.0*3.0
  !  
  !  P_{M+1,M}
  !  
  !  p10(x,y) = x*p00(x,y)
  !  p21(x,y) = x*3.0*p11(x,y)
  !  p32(x,y) = x*5.0*p22(x,y)
  !  p43(x,y) = x*7.0*p33(x,y)
  !  p54(x,y) = x*9.0*p44(x,y)
  !  p65(x,y) = x*11.0*p55(x,y)
  !  p76(x,y) = x*13.0*p66(x,y)
  !  p87(x,y) = x*15.0*p77(x,y)
  !  p98(x,y) = x*17.0*p88(x,y)
  !  
  !  P_{L,M} L=M+2, M+3, ... , lpole : M=0
  !  
  !  p20(x,y) = (x*3.0*p10(x,y)-p00(x,y))/2.0
  !  p30(x,y) = (x*5.0*p20(x,y)-2.0*p10(x,y))/3.0
  !  p40(x,y) = (x*7.0*p30(x,y)-3.0*p20(x,y))/4.0
  !  p50(x,y) = (x*9.0*p40(x,y)-4.0*p30(x,y))/5.0
  !  p60(x,y) = (x*11.0*p50(x,y)-5.0*p40(x,y))/6.0
  !  p70(x,y) = (x*13.0*p60(x,y)-6.0*p50(x,y))/7.0
  !  p80(x,y) = (x*15.0*p70(x,y)-7.0*p60(x,y))/8.0
  !  p90(x,y) = (x*17.0*p80(x,y)-8.0*p70(x,y))/9.0
  !  
  !  P_{L,M} L=M+2, M+3, ... , lpole : M=1
  !  
  !  p31(x,y) = (x*5.0*p21(x,y)-3.0*p11(x,y))/2.0
  !  p41(x,y) = (x*7.0*p31(x,y)-4.0*p21(x,y))/3.0
  !  p51(x,y) = (x*9.0*p41(x,y)-5.0*p31(x,y))/4.0
  !  p61(x,y) = (x*11.0*p51(x,y)-6.0*p41(x,y))/5.0
  !  p71(x,y) = (x*13.0*p61(x,y)-7.0*p51(x,y))/6.0
  !  p81(x,y) = (x*15.0*p71(x,y)-8.0*p61(x,y))/7.0
  !  p91(x,y) = (x*17.0*p81(x,y)-9.0*p71(x,y))/8.0
  !  
  !  P_{L,M} L=M+2, M+3, ... , lpole : M=2
  !  
  !  p42(x,y) = (x*7.0*p32(x,y)-5.0*p22(x,y))/2.0
  !  p52(x,y) = (x*9.0*p42(x,y)-6.0*p32(x,y))/3.0
  !  p62(x,y) = (x*11.0*p52(x,y)-7.0*p42(x,y))/4.0
  !  p72(x,y) = (x*13.0*p62(x,y)-8.0*p52(x,y))/5.0
  !  p82(x,y) = (x*15.0*p72(x,y)-9.0*p62(x,y))/6.0
  !  p92(x,y) = (x*17.0*p82(x,y)-10.0*p72(x,y))/7.0
  !  
  !  P_{L,M} L=M+2, M+3, ... , lpole : M=3
  !  
  !  p53(x,y) = (x*9.0*p43(x,y)-7.0*p33(x,y))/2.0
  !  p63(x,y) = (x*11.0*p53(x,y)-8.0*p43(x,y))/3.0
  !  p73(x,y) = (x*13.0*p63(x,y)-9.0*p53(x,y))/4.0
  !  p83(x,y) = (x*15.0*p73(x,y)-10.0*p63(x,y))/5.0
  !  p93(x,y) = (x*17.0*p83(x,y)-11.0*p73(x,y))/6.0
  !  
  !  P_{L,M} L=M+2, M+3, ... , lpole : M=4
  !  
  !  p64(x,y) = (x*11.0*p54(x,y)-9.0*p44(x,y))/2.0
  !  p74(x,y) = (x*13.0*p64(x,y)-10.0*p54(x,y))/3.0
  !  p84(x,y) = (x*15.0*p74(x,y)-11.0*p64(x,y))/4.0
  !  p94(x,y) = (x*17.0*p84(x,y)-12.0*p74(x,y))/5.0
  !  
  !  P_{L,M} L=M+2, M+3, ... , lpole : M=5
  !  
  !  p75(x,y) = (x*13.0*p65(x,y)-11.0*p55(x,y))/2.0
  !  p85(x,y) = (x*15.0*p75(x,y)-12.0*p65(x,y))/3.0
  !  p95(x,y) = (x*17.0*p85(x,y)-13.0*p75(x,y))/4.0
  !  
  !  P_{L,M} L=M+2, M+3, ... , lpole : M=6
  !  
  !  p86(x,y) = (x*15.0*p76(x,y)-13.0*p66(x,y))/2.0
  !  p96(x,y) = (x*17.0*p86(x,y)-14.0*p76(x,y))/3.0
  !  
  !  P_{L,M} L=M+2, M+3, ... , lpole : M=7
  !  
  !  p97(x,y) = (x*17.0*p87(x,y)-15.0*p77(x,y))/2.0
  !
  !  Rather than defining these polynomials recursively as above it
  !  is much more efficient to simplify them all. As it turns out
  !  the numerical coefficients are all powers of two, so there are
  !  no repeating decimals and we avoid having to do computationally
  !  intensive things like divide by 3 or 7. 
  !  If x = cos(theta) and y = sin(theta) we get the following 
  !  simplifications for the associated Legendre polynomials:
  !
  plm(0,0) = 1.0
  if (lpole >= 1) then
     plm(1,0) = x
     plm(1,1) = -y
  endif
  if (lpole >= 2) then
     x2 = x*x
     y2 = y*y
     plm(2,0) = 1.5*x2 - 0.5
     plm(2,1) = -3.0*x*y
     plm(2,2) = 3.0*y2
  endif
  if (lpole >= 3) then
     x3 = x2*x
     y3 = y2*y
     plm(3,0) = 2.5*x3 - 1.5*x
     plm(3,1) = -7.5*x2*y + 1.5*y
     plm(3,2) = 15.0*x*y2
     plm(3,3) = -15.0*y3
  endif
  if (lpole >= 4) then
     x4 = x2*x2
     y4 = y2*y2
     plm(4,0) = 4.375*x4 - 3.75*x2 + 0.375
     plm(4,1) = -17.5*x3*y + 7.5*x*y
     plm(4,2) = 52.5*x2*y2 - 7.5*y2
     plm(4,3) = -105.0*x*y3
     plm(4,4) = 105.0*y4
  endif
  if (lpole >= 5) then
     x5 = x3*x2
     y5 = y3*y2
     plm(5,0) = 7.875*x5 - 8.75*x3 + 1.875*x
     plm(5,1) = -39.375*x4*y + 26.25*x2*y - 1.875*y
     plm(5,2) = 157.5*x3*y2 - 52.5*x*y2
     plm(5,3) = -472.5*x2*y3 + 52.5*y3
     plm(5,4) = 945.0*x*y4
     plm(5,5) = -945.0*y5
  endif
  if (lpole >= 6) then
     x6 = x3*x3
     y6 = y3*y3
     plm(6,0) = 14.4375*x6 - 19.6875*x4 + 6.5625*x2 - 0.3125
     plm(6,1) = -86.625*x5*y + 78.75*x3*y - 13.125*x*y
     plm(6,2) = 433.125*x4*y2 - 236.25*x2*y2 + 13.125*y2
     plm(6,3) = -1732.5*x3*y3 + 472.5*x*y3
     plm(6,4) = 5197.5*x2*y4 - 472.5*y4
     plm(6,5) = -10395.0*x*y5
     plm(6,6) = 10395.0*y6
  endif
  if (lpole >= 7) then
     x7 = x4*x3
     y7 = y4*y3
     plm(7,0) = 26.8125*x7 - 43.3125*x5 + 19.6875*x3 - 2.1875*x
     plm(7,1) = -187.6875*x6*y + 216.5625*x4*y - 59.0625*x2*y + 2.1875*y
     plm(7,2) = 1126.125*x5*y2 - 866.25*x3*y2 + 118.125*x*y2
     plm(7,3) = -5630.625*x4*y3 + 2598.75*x2*y3 - 118.125*y3
     plm(7,4) = 22522.5*x3*y4 - 5197.5*x*y4
     plm(7,5) = -67567.5*x2*y5 + 5197.5*y5
     plm(7,6) = 135135.0*x*y6
     plm(7,7) = -135135.0*y7
  endif
  if (lpole >= 8) then
     x8 = x4*x4
     y8 = y4*y4
     plm(8,0) = 50.2734375*x8 - 93.84375*x6 + 54.140625*x4 - &
          9.84375*x2 + 0.2734375
     plm(8,1) = -402.1875*x7*y + 563.0625*x5*y - 216.5625*x3*y + &
          19.6875*x*y
     plm(8,2) = 2815.3125*x6*y2 - 2815.3125*x4*y2 + 649.6875*x2*y2 - &
          19.6875*y2
     plm(8,3) = -16891.875*x5*y3 + 11261.25*x3*y3 - 1299.375*x*y3
     plm(8,4) = 84459.375*x4*y4 - 33783.75*x2*y4 + 1299.375*y4
     plm(8,5) = -337837.5*x3*y5 + 67567.5*x*y5
     plm(8,6) = 1013512.5*x2*y6 - 67567.5*y6
     plm(8,7) = -2027025.0*x*y7
     plm(8,8) = 2027025.0*y8
  endif
  if (lpole >= 9) then
     y9 = y5*y4
     x9 = x5*x4
     plm(9,0) = 94.9609375*x9 - 201.09375*x7 + 140.765625*x5 - &
          36.09375*x3 + 2.4609375*x
     plm(9,1) = -854.6484375*x8*y + 1407.65625*x6*y - &
          703.828125*x4*y + 108.28125*x2*y - 2.4609375*y
     plm(9,2) = 6837.1875*x7*y2 - 8445.9375*x5*y2 + 2815.3125*x3*y2 &
          - 216.5625*x*y2
     plm(9,3) = -47860.3125*x6*y3 + 42229.6875*x4*y3 &
          - 8445.9375*x2*y3 + 216.5625*y3
     plm(9,4) = 287161.875*x5*y4 - 168918.75*x3*y4 + 16891.875*x*y4
     plm(9,5) = -1435809.375*x4*y5 + 506756.25*x2*y5 - 16891.875*y5
     plm(9,6) = 5743237.5*x3*y6 - 1013512.5*x*y6
     plm(9,7) = -17229712.5*x2*y7 + 1013512.5*y7
     plm(9,8) = 34459425.0*x*y8
     plm(9,9) = -34459425.0*y9
  endif

end subroutine setplm
!===============================================================
