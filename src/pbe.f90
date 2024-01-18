!===============================================================
!
! Copyright (C) 2005 Finite Difference Research Group
! This file is part of parsec, http://www.ices.utexas.edu/parsec/
! This is a slightly modified version of EASYPBE, a driver for
! the PBE subroutines, using simple inputs.
! author: K. Burke, May 14, 1996.
! Subroutines upgraded to the fortran 95 standard by M. Tiago,
! February 19, 2005.
!
!---------------------------------------------------------------
subroutine pbe(rho,agr,delgr,lap,expbe,vxuppbe,vxdnpbe,ecpbe, &
     vcuppbe,vcdnpbe)

  use constants
  implicit none
  !
  ! Input/Output variables:
  !
  real(dp), intent(inout) :: &
       rho                  ! density
  real(dp), intent(in) :: &
       agr, &               ! |grad rho|
       delgr, &             ! (grad rho).(grad |grad rho|)
       lap                  ! grad^2 rho=Laplacian of rho

  real(dp), intent(out) :: &
       expbe, &             ! exchange energy density, so that
                            ! Ex=int d^3r rho(r) ex(r)
       vxuppbe, &           ! up PBE exchange potential
       vxdnpbe, &           ! down PBE exchange potential
       ecpbe, &             ! PBE correlation energy density
       vcuppbe, &           ! up PBE correlation potential
       vcdnpbe              ! down PBE correlation potential
  !
  ! Work variables:
  !
  real(dp) :: dn,up,dnlap,uplap,agrdn,agrup,dvcdn,dvcup, &
       delgrdn,delgrup,ec,exdnpbe,exuppbe,fk,g,h,rho2,rholap, &
       rs,s,sk,u,uu,t,twoksg,v,vcdn,vcup,vv,ww,zet

  real(dp), parameter :: thrd2 = two*third
  real(dp), parameter :: pi32 = three*pi*pi
  real(dp), parameter :: alpha = 1.91915829267751300662482032624669d0
  !---------------------------------------------------------------

  !
  ! circumvent spin polarization
  up=rho/two
  dn=rho/two
  agrup=agr/two
  agrdn=agr/two
  delgrup=delgr/four
  delgrdn=delgr/four
  uplap=lap/two
  dnlap=lap/two
  !
  ! PBE exchange
  ! use  Ex[up,dn]=0.5*(Ex[2*up]+Ex[2*dn]) (i.e., exact spin-scaling)
  ! do up exchange
  ! fk=local Fermi wavevector for 2*up=(3 pi^2 (2up))^(1/3)
  ! s=dimensionless density gradient=|grad rho|/ (2*fk*rho)_(rho=2*up)
  ! u=delgrad/(rho^2*(2*fk)**3)_(rho=2*up)
  ! v=Laplacian/(rho*(2*fk)**2)_(rho=2*up)
  !
  rho2=two*up
  if(rho2 > 1d-18)then
     fk=(pi32*rho2)**third
     s=two*agrup/(two*fk*rho2)
     u=four*delgrup/(rho2*rho2*(two*fk)**3)
     v=two*uplap/(rho2*four*fk*fk)
     call exchpbe(rho2,s,u,v,exuppbe,vxuppbe)
  else
     exuppbe=zero
     vxuppbe=zero
  endif
  ! repeat for down
  rho2=two*dn
  if(rho2 > 1d-18)then
     fk=(pi32*rho2)**third
     s=two*agrdn/(two*fk*rho2)
     u=four*delgrdn/(rho2*rho2*(two*fk)**3)
     v=two*dnlap/(rho2*four*fk*fk)
     call exchpbe(rho2,s,u,v,exdnpbe,vxdnpbe)
  else
     exdnpbe=zero
     vxdnpbe=zero
  endif
  ! construct total density and contribution to ex
  rho=up+dn
  expbe=(exuppbe*up+exdnpbe*dn)/rho
  !
  ! Now do correlation
  !     zet=(up-dn)/rho
  !     g=phi(zeta)
  !     rs=(3/(4pi*rho))^(1/3)=local Seitz radius=alpha/fk
  !     sk=Ks=Thomas-Fermi screening wavevector=sqrt(4fk/pi)
  !     twoksg=2*Ks*phi
  !     t=correlation dimensionless gradient=|grad rho|/(2*Ks*phi*rho)
  !     uu=delgrad/(rho^2*twoksg^3)
  !     rholap=Laplacian
  !     vv=Laplacian/(rho*twoksg^2)
  !     ww=(|grad up|^2-|grad dn|^2-zet*|grad rho|^2)/(rho*twoksg)^2
  !     ec=LSD correlation energy
  !     vcup=LSD up correlation potential
  !     vcdn=LSD down correlation potential
  !     h=gradient correction to correlation energy
  !     dvcup=gradient correction to up correlation potential
  !     dvcdn=gradient correction to down correlation potential
  !
  if(rho < 1.d-18)return
  zet=(up-dn)/rho
  g=((one+zet)**thrd2+(one-zet)**thrd2)/two
  fk=(pi32*rho)**third
  rs=alpha/fk
  sk=sqrt(four*fk/pi)
  twoksg=two*sk*g
  t=agr/(twoksg*rho)
  uu=delgr/(rho*rho*twoksg**3)
  rholap=uplap+dnlap
  vv=rholap/(rho*twoksg*twoksg)
  ww=(agrup*agrup-agrdn*agrdn-zet*agr*agr)/(rho*rho*twoksg*twoksg)
  call corpbe(rs,zet,t,uu,vv,ww,ec,vcup,vcdn,h,dvcup,dvcdn)
  ecpbe=ec+h
  vcuppbe=vcup+dvcup
  vcdnpbe=vcdn+dvcdn

end subroutine pbe
!===============================================================
!
! official pbe correlation code. K. Burke, May 14, 1996.
!
! References:
! [a] J.P.~Perdew, K.~Burke, and M.~Ernzerhof , PRL 77, 3865 (1996)
! [b] J. P. Perdew, K. Burke, and Y. Wang, {\sl Real-space cutoff
!     construction of a generalized gradient approximation:  The PW91
!     density functional}, submitted to Phys. Rev. B, Feb. 1996.
! [c] J. P. Perdew and Y. Wang, Phys. Rev. B 45, 13244 (1992).
!
!---------------------------------------------------------------
subroutine corpbe(rs,zet,t,uu,vv,ww,ec,vcup,vcdn,h,dvcup,dvcdn)

  use constants
  implicit none
  !
  ! Input/Output variables:
  !
  real(dp), intent(in) :: &
       rs, &   ! Seitz radius=(3/4pi rho)^(1/3)
       zet, &  ! relative spin polarization = (rhoup-rhodn)/rho
       t, &    ! ABS(GRAD rho)/(rho*2.*ks*g)  -- only needed for PBE
       uu, &   ! (GRAD rho)*GRAD(ABS(GRAD rho))/(rho**2 * (2*ks*g)**3)
       vv, &   ! (LAPLACIAN rho)/(rho * (2*ks*g)**2)
       ww      ! (GRAD rho)*(GRAD zet)/(rho * (2*ks*g)**2

  real(dp), intent(out) :: &
       ec, &   ! LSD correlation energy from [a]
       vcup, & ! LSD up correlation potential
       vcdn, & ! LSD dn correlation potential
       h, &    ! nonlocal part of correlation energy per electron
       dvcup, &! nonlocal correction to vcup
       dvcdn   ! nonlocal correction to vcdn
  !
  ! Work variables:
  !
  real(dp) :: alfc,alfm,alfrsm,b,b2,bec,bg,comm,ecrs,eczet, &
       ep,eprs,eu,eurs,f,fac,fact0,fact1,fact2,fact3,fact5,fz,g &
       ,g3,g4,gz,hb,hbt,hrs,hrst,ht,htt,hz,hzt,pon,pref,q4,q5,q8 &
       ,q9,rs2,rs3,rsthrd,rtrs,t2,t4,t6,z4
  !
  ! third*=various multiples of 1/3
  ! numbers for use in LSD energy spin-interpolation formula, [c](9).
  !    gam= 2^(4/3)-2
  !    fzz=f''(0)= 8/(9*gam)
  ! numbers for construction of PBE
  !    gamma=(1-log(2))/pi^2
  !    bet=coefficient in gradient expansion for correlation, [a](4).
  !    eta=small number to stop d phi/ dzeta from blowing up at 
  !        |zeta|=1.
  !
  real(dp), parameter :: thrdm=-third
  real(dp), parameter :: thrd2=two*third
  real(dp), parameter :: sixthm=thrdm/two
  real(dp), parameter :: thrd4=four*third
  real(dp), parameter :: gam=0.5198420997897463295344212145565d0
  real(dp), parameter :: fzz=8.d0/(9.d0*gam)
  real(dp), parameter :: gamma=0.03109069086965489503494086371273d0
  real(dp), parameter :: bet=0.06672455060314922d0
  real(dp), parameter :: delt=bet/gamma
  real(dp), parameter :: eta=1.d-12
  !---------------------------------------------------------------
  !
  ! find LSD energy contributions, using [c](10) and table i[c].
  !     eu=unpolarized LSD correlation energy
  !     eurs=deu/drs
  !     ep=fully polarized LSD correlation energy
  !     eprs=dep/drs
  !     alfm=-spin stiffness, [c](3).
  !     alfrsm=-dalpha/drs
  !     f=spin-scaling factor from [c](9).
  ! construct ec, using [c](8)
  !
  rtrs=sqrt(rs)
  call gcor2(0.0310907d0,0.21370d0,7.5957d0,3.5876d0,1.6382d0, &
       0.49294d0,rtrs,eu,eurs)
  call gcor2(0.01554535d0,0.20548d0,14.1189d0,6.1977d0,3.3662d0, &
       0.62517d0,rtrs,ep,eprs)
  call gcor2(0.0168869d0,0.11125d0,10.357d0,3.6231d0,0.88026d0, &
       0.49671d0,rtrs,alfm,alfrsm)
  alfc = -alfm
  z4 = zet**4
  f=((one+zet)**thrd4+(one-zet)**thrd4-two)/gam
  ec = eu*(one-f*z4)+ep*f*z4-alfm*f*(one-z4)/fzz
  !
  ! LSD potential from [c](a1)
  !     ecrs = dec/drs [c](a2)
  !     eczet=dec/dzeta [c](a3)
  !     fz = df/dzeta [c](a4)
  !
  ecrs = eurs*(one-f*z4)+eprs*f*z4-alfrsm*f*(one-z4)/fzz
  fz = thrd4*((one+zet)**third-(one-zet)**third)/gam
  eczet = four*(zet**3)*f*(ep-eu+alfm/fzz)+fz*(z4*ep-z4*eu &
       -(one-z4)*alfm/fzz)
  comm = ec -rs*ecrs/three-zet*eczet
  vcup = comm + eczet
  vcdn = comm - eczet
  !
  ! PBE correlation energy
  !     g=phi(zeta), given after [a](3)
  !     delt=bet/gamma
  !     b=a of [a](8)
  !
  g=((one+zet)**thrd2+(one-zet)**thrd2)/two
  g3 = g**3
  pon=-ec/(g3*gamma)
  b = delt/(exp(pon)-one)
  b2 = b*b
  t2 = t*t
  t4 = t2*t2
  rs2 = rs*rs
  rs3 = rs2*rs
  q4 = one+b*t2
  q5 = one+b*t2+b2*t4
  h = g3*(bet/delt)*log(one+delt*q4*t2/q5)
  !
  ! energy done. now the potential, using appendix E of [b].
  !
  g4 = g3*g
  t6 = t4*t2
  rsthrd = rs/three
  gz=(((one+zet)*(one+zet)+eta)**sixthm- &
       ((one-zet)*(one-zet)+eta)**sixthm)/three
  fac = delt/b+one
  bg = -three*b2*ec*fac/(bet*g4)
  bec = b2*fac/(bet*g3)
  q8 = q5*q5+delt*q4*q5*t2
  q9 = one+two*b*t2
  hb = -bet*g3*b*t6*(two+b*t2)/q8
  hrs = -rsthrd*hb*bec*ecrs
  fact0 = two*delt-6.d0*b
  fact1 = q5*q9+q4*q9*q9
  hbt = two*bet*g3*t4*((q4*q5*fact0-delt*fact1)/q8)/q8
  hrst = rsthrd*t2*hbt*bec*ecrs
  hz = three*gz*h/g + hb*(bg*gz+bec*eczet)
  ht = two*bet*g3*q9/q8
  hzt = three*gz*ht/g+hbt*(bg*gz+bec*eczet)
  fact2 = q4*q5+b*t2*(q4*q9+q5)
  fact3 = two*b*q5*q9+delt*fact2
  htt = four*bet*g3*t*(two*b/q8-(q9*fact3/q8)/q8)
  comm = h+hrs+hrst+t2*ht/6.d0+7.d0*t2*t*htt/6.d0
  pref = hz-gz*t2*ht/g
  fact5 = gz*(two*ht+t*htt)/g
  comm = comm-pref*zet-uu*htt-vv*ht-ww*(hzt-fact5)
  dvcup = comm + pref
  dvcdn = comm - pref

end subroutine corpbe
!===============================================================
!
! PBE exchange for a spin-unpolarized electronic system
! K Burke's modification of PW91 codes, May 14, 1996
! Modified again by K. Burke, June 29, 1996, with simpler Fx(s)
!
! References:
! [a]J.P.~Perdew, K.~Burke, and M.~Ernzerhof, PRL 77, 3865 (1996);
! [b]J.P. Perdew and Y. Wang, Phys. Rev.  B 33,  8800(1986);
! 40,  3399  (1989) (E).
!
! Formulas:
!      e_x[unif]=ax*rho^(4/3)  [LDA]
!      ax = -0.75*(3/pi)^(1/3)
!      e_x[PBE]=e_x[unif]*FxPBE(s)
!      FxPBE(s)=1+uk-uk/(1+ul*s*s)                 [a](13)
! uk, ul defined after [a](13) 
!
! for u,v, see PW86(24)
!
!---------------------------------------------------------------
subroutine exchpbe(rho,s,u,v,ex,vx)

  use constants
  implicit none
  !
  ! Input/Output variables:
  !
  real(dp), intent(in) :: &
       rho, &  ! density
       s, &    ! ABS(GRAD rho)/(2*KF*rho), where kf=(3 pi^2 rho)^(1/3)
       u, &    ! (GRAD rho)*GRAD(ABS(GRAD rho))/(rho**2 * (2*kf)**3)
       v       ! (LAPLACIAN rho)/(rho*(2*kf)**2)

  real(dp), intent(out) :: &
       ex, &   ! exchange energy per electron
       vx      ! potential
  !
  ! Work variables:
  !
  real(dp) :: exunif,fs,fss,fxpbe,p0,s2

  real(dp), parameter :: thrd4=four/three
  real(dp), parameter :: ax=-0.738558766382022405884230032680836d0
  real(dp), parameter :: um=0.2195149727645171d0
  real(dp), parameter :: uk=0.8040d0
  real(dp), parameter :: ul=um/uk
  !---------------------------------------------------------------
  !
  ! construct lda exchange energy density
  exunif = ax*rho**third
  !
  ! construct PBE enhancement factor
  s2 = s*s
  p0=one+ul*s2
  fxpbe = one+uk-uk/p0
  ex = exunif*fxpbe
  !
  ! energy done. now the potential:
  ! find first and second derivatives of fx w.r.t s.
  !      fs=(1/s)*d fxpbe/ ds
  !      fss=d fs/ds
  !
  fs=two*uk*ul/(p0*p0)
  fss=-four*ul*s*fs/p0
  !
  ! calculate potential from [b](24) 
  vx = exunif*(thrd4*fxpbe-(u-thrd4*s2*s)*fss-v*fs)

end subroutine exchpbe
!===============================================================
!
! Slimmed down version of gcor used in PW91 routines, to
! interpolate LSD correlation energy, as given by (10) of
! J. P. Perdew and Y. Wang, Phys. Rev. B 45, 13244 (1992).
! K. Burke, May 11, 1996.
!
!---------------------------------------------------------------
subroutine gcor2(a,a1,b1,b2,b3,b4,rtrs,gg,ggrs)

  use constants
  implicit none
  !
  ! Input/Output variables:
  !
  real(dp), intent(in) :: a,a1,b1,b2,b3,b4,rtrs

  real(dp), intent(out) :: gg,ggrs
  !
  ! Work variables:
  !
  real(dp) :: q0,q1,q2,q3
  !---------------------------------------------------------------
  q0 = -two*a*(one+a1*rtrs*rtrs)
  q1 = two*a*rtrs*(b1+rtrs*(b2+rtrs*(b3+b4*rtrs)))
  q2 = log(one+one/q1)
  gg = q0*q2
  q3 = a*(b1/rtrs+two*b2+rtrs*(three*b3+four*b4*rtrs))
  ggrs = -two*a*a1*q2-q0*q3/(q1*(one+q1))

end subroutine gcor2
!===============================================================
!
! PW91 spin-unpolarized functional, extracted from spn_pbe
!
!---------------------------------------------------------------
subroutine pw91(rho,agr,delgr,lap,expw91,vxpw91,ecpw91,vcuppw91,vcdnpw91)

  use constants
  implicit none
  !
  ! Input/Output variables:
  !
  real(dp), intent(in) :: &
       rho, &                 ! density
       agr, &                 ! |grad rho|
       delgr, &               ! (grad up).(grad |grad up|)
       lap                    ! grad^2 up=Laplacian of up
  real(dp), intent(out) :: &
       expw91,vxpw91,ecpw91,vcuppw91,vcdnpw91
  ! PW91 quantities
  !
  ! Work variables:
  !
  real(dp) :: alfc,dvcdn,dvcup,ec,ecrs,eczet,fk,g,h, &
       rholap,rs,s,sk,twoksg,t,u,uu,vv,v,vcdn,vcup,ww,zet
  !
  ! constants: pi32=3 pi**2, alpha=(9pi/4)**thrd
  !
  real(dp), parameter :: thrd2 = two*third
  real(dp), parameter :: pi32 = three*pi*pi
  real(dp), parameter :: alpha = 1.91915829267751300662482032624669d0

  !---------------------------------------------------------------
  !
  if(rho < 1.d-18) then
     expw91 = zero
     vxpw91 = zero
     ecpw91 = zero
     vcuppw91 = zero
     vcdnpw91 = zero
  else
     fk=(pi32*rho)**third
     s=agr/(two*fk*rho)
     u=two*delgr/(rho*rho*(two*fk)**3)
     v=lap/(rho*(two*fk)**2)
     call exchpw91(rho,s,u,v,expw91,vxpw91)

     zet=zero
     g=one
     rs=alpha/fk
     sk=sqrt(four*fk/pi)
     twoksg=two*sk*g
     t=agr/(twoksg*rho)
     uu=delgr/(rho*rho*twoksg**3)
     rholap=lap
     vv=rholap/(rho*twoksg**2)
     ww=zero
     call corlsd(rs,zet,ec,vcup,vcdn,ecrs,eczet,alfc)
     call corpw91(rs,zet,g,ec,ecrs,eczet,t,uu,vv,ww,h,dvcup,dvcdn)
     ecpw91=ec+h
     vcuppw91=vcup+dvcup
     vcdnpw91=vcdn+dvcdn
  endif

end subroutine pw91
!===============================================================
