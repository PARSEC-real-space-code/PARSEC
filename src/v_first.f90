!===============================================================
!
! Copyright (C) 2005 Finite Difference Research Group
! This file is part of parsec, http://www.ices.utexas.edu/parsec/
!
! Calculates the local ionic potential in the g-vectors (i.e.
! interpolates vloc and multiplies by the structure factor).
! It also calculates the core charge density (dnc) this is going
! to be used later when calculating the local force on ions.
!
! Adapted from plane-wave programs written by S. Froyen and
! J. L. Martins.
!
!--------------------------------------------------------------
subroutine v_first(pbc,natom,vionz)!{{{

  use constants
  use pbc_module
  implicit none
  !
  ! Input/Output variables:
  !
  ! pbc related data
  type(pbc_data), intent(inout) :: pbc
  ! number of atoms of each type
  integer, intent(in) :: natom(pbc%type_num)
  ! ionic potential
  complex(dpc), dimension(pbc%ng), intent(out) :: vionz
  !
  ! Work variables:
  !
  complex(dpc) :: strfac
  integer ik,ig,iat,j,n,nt,nql,ns
  real(dp) :: vcell,gmax,glmax,delql,fi
  real(dp) :: qj,xn,q2vn,q2vp,q2vm,vqj,dcj
  ! not used at all
  ! real(dp) :: dvqj
  !---------------------------------------------------------------

  ! scaling
  vcell = pbc%vcell
  ! initializations
  ns = pbc%nstar
  vionz(:) = zzero

  do nt=1,pbc%type_num
     ! Check argument limits and compare with gmax.
     delql = pbc%delq(nt)
     nql = pbc%nq(nt)
     gmax = sqrt(two*pbc%ek(ns))
     glmax = delql * real(nql-3,dp)
     if (glmax < gmax) then
         write(7,*) " *** WARNING: in v_first, for nt=",nt
         write(7,10) glmax
         write(7,*) " *** while gmax was ",gmax
     endif
     !
     ! Compute ionic local potential. Ignore first point (G = 0).
     !
     ig = pbc%mstar(1)
     do j=2,ns
        pbc%vql(nt,j) = zero
        !
        ! Quadratic interpolation of q**2 * potential.
        !
        qj = sqrt(two*pbc%ek(j))
        xn = qj / delql + two
        !AJB: should be probperly cast into int.
        n  = xn + half
        if (n < nql) then
           if (n <= 3) n = 4
           xn = xn - real(n,dp)

           q2vn = real((n-2)*(n-2),dp)*pbc%vloc(n,nt)
           q2vp = real((n-1)*(n-1),dp)*pbc%vloc(n+1,nt)
           q2vm = real((n-3)*(n-3),dp)*pbc%vloc(n-1,nt)
           vqj  = q2vn * (one-xn) * (one+xn) + &
                half * (q2vp*(one+xn) - q2vm*(one-xn)) * xn
           vqj  = delql*delql * vqj / (vcell * two*pbc%ek(j))
           ! AJB: DEAD CODE 2015-02-22 Sun 7:00 PM 
           ! dvqj = -two * q2vn * xn + &
           !      q2vp * (half+xn) - q2vm * (half-xn)
           ! dvqj = (delql*dvqj/(two*vcell*qj) - vqj)/(two*pbc%ek(j))

           ! Calculate structure factor by doing direct sum over atoms of
           ! this type and over vectors in this star.
           do ik = 1, pbc%mstar(j)
              strfac = cmplx(zero,zero,dp)
              ig = ig + 1
              do iat = 1, natom(nt)
                 fi = dot_product(real(pbc%kgv(:,ig),dp),pbc%rat(:,iat,nt))
                 strfac = strfac + cmplx(cos(fi),sin(fi))
              enddo
              vionz(ig) = vionz(ig) + vqj*conjg(strfac)
           enddo
           pbc%vql(nt,j) = vqj
        endif
     enddo

     ! compute core charge
     do j=1,ns
        pbc%dnc(nt,j) = zero

        ! interpolate vda
        qj = sqrt(two*pbc%ek(j))
        xn = qj/delql + two
        n = xn + half
        if (n < nql) then
           xn = xn - real(n,dp)
           dcj = pbc%dcor(n,nt) * (one+xn) * (one-xn) &
                + half * (pbc%dcor(n+1,nt)*(one+xn) &
                -pbc%dcor(n-1,nt)*(one-xn)) * xn
           pbc%dnc(nt,j) = dcj
        endif
     enddo
     ! end loop over atomic types
  enddo
10 format('  *** v_first contains information', &
        ' about local potential only up to glmax = ',e10.3)

end subroutine v_first
!===============================================================
!}}}
!===============================================================
!
!
! Calculates the local ionic potential in the g-vectors 
! (i.e.  interpolates vloc and multiplies by the structure factor).
! It also calculates the core charge density (dnc) this is going
! to be used later when calculating the local force on ions.
!
!
!--------------------------------------------------------------
subroutine v_first_omp(pbc,natom,vionz) !{{{

  use constants
  use pbc_module
  implicit none
  !
  ! Input/Output variables:
  !
  ! pbc related data
  type(pbc_data), intent(inout) :: pbc
  ! number of atoms of each type
  integer, intent(in) :: natom(pbc%type_num)
  ! ionic potential
  complex(dpc), dimension(pbc%ng), intent(out) :: vionz

  !
  ! Work variables:
  !
  complex(dpc) :: expfi !temporary for the phase
  ! helper array for the k vector intex, cumulated over stars
  integer :: cum_mstar(pbc%ng)

!  complex(dpc), dimension(pbc%ng) :: vionz_alt
  complex(dpc) :: strfac
  integer, dimension(pbc%ng) :: ig_array
  integer ik,ig,iat,i,j,n,nt,nql,ns
  real(dp) :: vcell,gmax,glmax,delql,fi
  real(dp) :: qj,xn,q2vn,q2vp,q2vm,vqj,dcj,xn_tmp
  ! not used at all
  ! real(dp) :: dvqj
  !---------------------------------------------------------------

  ! scaling
  vcell = pbc%vcell
  ! initializations
  ns = pbc%nstar
  vionz(:) = zzero
  !
  ! create the cum_mstar array:
  cum_mstar(1)=pbc%mstar(1)
  do j=2,ns
!           write(9,*) "v_alt debug, doing star = ",j
  cum_mstar(j)=cum_mstar(j-1)+pbc%mstar(j)
  enddo
!  write(9,*) "v_alt debug, done with cum_mstar creation"
  


  do nt=1,pbc%type_num

     ! Check argument limits and compare with gmax.
     delql = pbc%delq(nt) !a global parameter?
     nql = pbc%nq(nt)
     gmax = sqrt(two*pbc%ek(ns))
     glmax = delql * real(nql-3,dp)
     if (glmax < gmax) then
         write(7,*) " *** WARNING: in v_fist_omp, for nt=",nt
         write(7,10) glmax
         write(7,*) " *** while gmax was ",gmax
     endif
10 format('  *** v_first_omp(test) contains information', &
        ' about local potential only up to glmax = ',f10.5)

     ! compute core charge for star j=1
     j = 1
        pbc%dnc(nt,j) = zero
        qj = sqrt(two*pbc%ek(j))
        xn = qj/delql + two
        n = xn + half
        if (n < nql) then
           xn = xn - real(n,dp)
           dcj = pbc%dcor(n,nt) * (one+xn) * (one-xn) &
                + half * (pbc%dcor(n+1,nt)*(one+xn) &
                -pbc%dcor(n-1,nt)*(one-xn)) * xn
!           write(9,*) "v_first_omp debug, dc_j=1 =  ",dcj
        pbc%dnc(nt,j) = dcj
        endif
     !
     ! Compute ionic local potential. Ignore first point (G = 0).
     ! also deal with core charge
     !
!$OMP PARALLEL DO &
!$OMP& SCHEDULE(DYNAMIC) &
!$OMP& DEFAULT(SHARED) &
!$OMP& PRIVATE(j,qj,xn,xn_tmp,n,dcj,q2vn,q2vp,q2vm,vqj,ik,strfac,iat,fi,expfi) 
! !do not need: !$OMP& REDUCTION(+:vionz) 
     do j=2,ns
!           write(9,*) "v_alt debug, doing star = ",j
        !uncomment after copying over v_first
        pbc%vql(nt,j) = zero
        pbc%dnc(nt,j) = zero !core charge
        !
        ! Quadratic interpolation of q**2 * potential.
        ! (questionable if needed at all)
        !
        qj = sqrt(two*pbc%ek(j))
        xn = qj / delql + two
        n  = xn + half
        if (n < nql) then
           !added core chrge loop here, 
           !dcor(1,ity)=dcor(2,ity)=dcor(3,ity) from fourier_f.f90
           !why is this not interpolated the same way as q2? dunno
           !because of the different schemes the cases for n<=3
           !need some special treatment
           if (n <= 3) then
               xn_tmp = xn - real(n,dp)
               dcj = pbc%dcor(n,nt) * (one+xn_tmp) * (one-xn_tmp) &
                    + half * (pbc%dcor(n+1,nt)*(one+xn_tmp) &
                    -pbc%dcor(n-1,nt)*(one-xn_tmp)) * xn_tmp

                   n = 4 !beacuse of the interpolation?
               xn = xn - real(n,dp)
            else
               xn = xn - real(n,dp)
               dcj = pbc%dcor(n,nt) * (one+xn) * (one-xn) &
                    + half * (pbc%dcor(n+1,nt)*(one+xn) &
                    -pbc%dcor(n-1,nt)*(one-xn)) * xn
            endif
           !vloc(1,ity)=vloc(2,ity)=vloc(3,ity) from fourier_f.f90
           q2vn = real((n-2)*(n-2),dp)*pbc%vloc(n,nt)
           q2vp = real((n-1)*(n-1),dp)*pbc%vloc(n+1,nt)
           q2vm = real((n-3)*(n-3),dp)*pbc%vloc(n-1,nt)
           vqj  = q2vn * (one-xn) * (one+xn) + &
                half * (q2vp*(one+xn) - q2vm*(one-xn)) * xn
           vqj  = delql*delql * vqj / (vcell * two*pbc%ek(j))

           ! AJB: DEAD CODE 2015-02-22 Sun 7:00 PM 
           ! dvqj = -two * q2vn * xn + &
           !      q2vp * (half+xn) - q2vm * (half-xn)
           ! dvqj = (delql*dvqj/(two*vcell*qj) - vqj)/(two*pbc%ek(j))

           pbc%dnc(nt,j) = dcj
           pbc%vql(nt,j) = vqj

           ! Calculate structure factor by doing direct sum:
           do ik = cum_mstar(j-1)+1,cum_mstar(j) !done this way for openmp
           !write(9,*) "v_alt debug, doing ik = ",ik
               strfac = cmplx(zero,zero,dp)
               do iat = 1, natom(nt)
                   fi = dot_product(real(pbc%kgv(:,ik),dp),pbc%rat(:,iat,nt))
                   expfi = cmplx(cos(fi),sin(fi))
                  strfac = strfac + expfi
               enddo
            vionz(ik) = vionz(ik) + vqj*conjg(strfac)
           enddo

        endif
     enddo
!$OMP END PARALLEL DO 

  !   end loop over atomic types
  enddo

end subroutine v_first_omp
!===============================================================
