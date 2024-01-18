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
subroutine spn_pbe(up,agrup,delgrup,uplap,dn,agrdn,delgrdn,dnlap, &
     agr,delgr,lcor,lpot, &
     exlsd,vxuplsd,vxdnlsd,eclsd,vcuplsd,vcdnlsd, &
     expw91,vxuppw91,vxdnpw91,ecpw91,vcuppw91,vcdnpw91, &
     expbe,vxuppbe,vxdnpbe,ecpbe,vcuppbe,vcdnpbe)

  use constants
  implicit none
  !
  ! Input/Output variables:
  !
  real(dp), intent(in) :: &
       up, &               ! up density
       agrup, &            ! |grad up|
       delgrup, &          ! (grad up).(grad |grad up|)
       uplap, &            ! grad^2 up=Laplacian of up
       dn,agrdn,delgrdn,dnlap, &
                           ! corresponding for down density
       agr, &              ! |grad rho|
       delgr               ! (grad rho).(grad |grad rho|)
  integer, intent(in) ::  &
       lcor, &             ! flag to do correlation(=0=>don't) 
       lpot                ! flag to do potential(=0=>don't)
  real(dp), intent(out) :: &
       exlsd, &            ! LSD exchange energy density, so that
                           ! ExLSD=int d^3r rho(r) exlsd(r)
       vxuplsd, &          ! up LSD exchange potential
       vxdnlsd, &          ! down LSD exchange potential
       eclsd, &            ! LSD correlation energy density
       vcuplsd, &          ! up LSD correlation potential
       vcdnlsd, &          ! down LSD correlation potential
       expw91,vxuppw91,vxdnpw91,ecpw91,vcuppw91,vcdnpw91, &
                           ! corresponding PW91 quantities
       expbe,vxuppbe,vxdnpbe,ecpbe,vcuppbe,vcdnpbe
                           ! corresponding PBE quantities
  !
  ! Work variables:
  !
  real(dp) :: alfc,dvcdn,dvcup,ec,ecrs,eczet,exdnlsd,exdnpbe, &
       exdnpw91,exuppbe,exuplsd,exuppw91,fk,g,h,rho,rholap,rho2, &
       rs,s,sk,twoksg,t,u,uu,vv,v,vcdn,vcup,ww,zet
  !
  ! constants: pi32=3 pi**2, alpha=(9pi/4)**thrd
  !
  real(dp), parameter :: thrd2 = two*third
  real(dp), parameter :: pi32 = three*pi*pi
  real(dp), parameter ::  alpha = 1.91915829267751300662482032624669d0

  !---------------------------------------------------------------
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
     v=two*uplap/(rho2*(two*fk)**2)

     call spn_exchpbe(rho2,s,u,v,0,lpot,exuplsd,vxuplsd)
     call exchpw91(rho2,s,u,v,exuppw91,vxuppw91)
     call spn_exchpbe(rho2,s,u,v,1,lpot,exuppbe,vxuppbe)
  else
     exuplsd=zero
     vxuplsd=zero
     exuppw91=zero
     vxuppw91=zero
     exuppbe=zero
     vxuppbe=zero
  endif
  ! repeat for down 
  rho2=two*dn
  if(rho2 > 1d-18)then
     fk=(pi32*rho2)**third
     s=two*agrdn/(two*fk*rho2)
     u=four*delgrdn/(rho2*rho2*(two*fk)**3)
     v=two*dnlap/(rho2*(two*fk)**2)

     ! note: 2 following lines are missing in pbe.f (adi)
     call spn_exchpbe(rho2,s,u,v,0,lpot,exdnlsd,vxdnlsd)
     call exchpw91(rho2,s,u,v,exdnpw91,vxdnpw91)
     call spn_exchpbe(rho2,s,u,v,1,lpot,exdnpbe,vxdnpbe)
  else

     ! note: 4 following lines are missing in pbe.f (adi)
     exdnlsd=zero
     vxdnlsd=zero
     exdnpw91=zero
     vxdnpw91=zero
     exdnpbe=zero
     vxdnpbe=zero
  endif
  ! construct total density and contribution to ex 
  ! note: different in pbe.f (adi)
  rho = up + dn
  if(up > 1d-18 .and. dn > 1d-18)then
     exlsd=(exuplsd*up+exdnlsd*dn)/rho
     expw91=(exuppw91*up+exdnpw91*dn)/rho
     expbe=(exuppbe*up+exdnpbe*dn)/rho
  else
     exlsd=zero
     expw91=zero
     expbe=zero
  endif
  if(lcor == 0)return
  !
  ! Now do correlation 
  ! zet=(up-dn)/rho 
  ! g=phi(zeta) 
  ! rs=(3/(4pi*rho))^(1/3)=local Seitz radius=alpha/fk 
  ! sk=Ks=Thomas-Fermi screening wavevector=sqrt(4fk/pi) 
  ! twoksg=2*Ks*phi 
  ! t=correlation dimensionless gradient=|grad rho|/(2*Ks*phi*rho) 
  ! uu=delgrad/(rho^2*twoksg^3) 
  ! rholap=Laplacian 
  ! vv=Laplacian/(rho*twoksg^2) 
  ! ww=(|grad up|^2-|grad dn|^2-zet*|grad rho|^2)/(rho*twoksg)^2 
  ! ec=LSD correlation energy 
  ! vcup=LSD up correlation potential 
  ! vcdn=LSD down correlation potential 
  ! h=gradient correction to correlation energy 
  ! dvcup=gradient correction to up correlation potential 
  ! dvcdn=gradient correction to down correlation potential 
  !
  if (rho < 1.d-18) return
  zet=(up-dn)/rho
  g=((one+zet)**thrd2+(one-zet)**thrd2)/two
  fk=(pi32*rho)**third
  rs=alpha/fk
  sk=sqrt(four*fk/pi)
  twoksg=two*sk*g
  t=agr/(twoksg*rho)
  uu=delgr/(rho*rho*twoksg**3)
  rholap=uplap+dnlap
  vv=rholap/(rho*twoksg**2)
  ww=(agrup**2-agrdn**2-zet*agr**2)/(rho*rho*twoksg**2)
  call spn_corpbe(rs,zet,t,uu,vv,ww,1,lpot,ec,vcup,vcdn,h,dvcup,dvcdn)
  eclsd=ec
  ecpbe=ec+h
  vcuplsd=vcup
  vcdnlsd=vcdn
  ! note: 7 following lines are missing in pbe.f
  vcuppbe=vcup+dvcup
  vcdnpbe=vcdn+dvcdn
  call corlsd(rs,zet,ec,vcup,vcdn,ecrs,eczet,alfc)
  call corpw91(rs,zet,g,ec,ecrs,eczet,t,uu,vv,ww,h,dvcup,dvcdn)
  ecpw91=ec+h
  vcuppw91=vcup+dvcup
  vcdnpw91=vcdn+dvcdn

end subroutine spn_pbe
!===============================================================
!
! PBE EXCHANGE FOR A SPIN-UNPOLARIZED ELECTRONIC SYSTEM 
! K Burke's modification of PW91 codes, May 14, 1996 
! Modified again by K. Burke, June 29, 1996, with simpler Fx(s) 
!
! (for u,v, see PW86(24)) 
!
! References: 
! [a]J.P.~Perdew, K.~Burke, and M.~Ernzerhof, PRL 77, 3865 (1996)
! [b]J.P. Perdew and Y. Wang, Phys. Rev.  B 33,  8800 (1986);
! 40,  3399  (1989) (E). 
!
! Formulas: 
!     e_x[unif]=ax*rho^(4/3)  [LDA] 
!     ax = -0.75*(3/pi)^(1/3) 
!     e_x[PBE]=e_x[unif]*FxPBE(s) 
!     FxPBE(s)=1+uk-uk/(1+ul*s*s)                 [a](13) 
! uk, ul defined after [a](13) 
!
!---------------------------------------------------------------
subroutine spn_exchpbe(rho,s,u,v,lgga,lpot,ex,vx) 

  use constants
  implicit none
  !
  ! Input/Output variables:
  !
  real(dp), intent(in) :: &
       rho, &  ! density
       s, &    ! ABS(GRAD rho)/(2*KF*rho), where kf=(3 pi^2 rho)^(1/3)
       u, &    ! (GRAD rho)*GRAD(ABS(GRAD rho))/(rho**2 * (2*KF)**3)
       v       ! (LAPLACIAN rho)/(rho*(2*KF)**2)
  integer, intent(in) ::  &
       lgga, & ! (=0=>don't put in gradient corrections, just LDA)
       lpot    ! (=0=>don't get potential and don't need U and V)
  real(dp), intent(out) :: &
       ex, &   ! exchange energy per electron
       vx      ! potential
  !
  ! Work variables:
  !
  real(dp) :: uk,p0,s2,fxpbe,fss,fs,exunif,ul

  real(dp), parameter :: thrd4 = four/three
  real(dp), parameter :: ax = -0.738558766382022405884230032680836d0
  real(dp), parameter :: um = 0.2195149727645171d0

  !---------------------------------------------------------------

  ! note: in pbe.f uk is defined as parameter: uk=0.8040d0
  uk = 0.8040d0
  ul=um/uk
  !
  ! construct LDA exchange energy density 
  exunif = AX*rho**third

  ! note: this part is missing in pbe.f (adi)
  if(lgga == 0)then
     ex=exunif
     vx=ex*thrd4
     return
  endif
  !
  ! construct PBE enhancement factor 
  S2 = S*S
  P0=one+ul*S2
  FxPBE = 1d0+uk-uk/P0
  EX = exunif*FxPBE
  ! note: the following line is missing in pbe.f (adi)
  if(lpot == 0)return
  !
  ! energy done. now the potential: 
  ! find first and second derivatives of Fx w.r.t s. 
  ! Fs=(1/s)*d FxPBE/ ds 
  ! Fss=d Fs/ds 
  Fs=two*uk*ul/(P0*P0)
  Fss=-four*ul*S*Fs/P0
  !
  ! calculate potential from [b](24) 
  vx = exunif*(thrd4*fxpbe-(u-thrd4*s2*s)*fss-v*fs)

end subroutine spn_exchpbe
!===============================================================
!
! Official PBE correlation code. K. Burke, May 14, 1996. 
!
! References: 
!   [a] J.P.~Perdew, K.~Burke, and M.~Ernzerhof, PRL 77, 3865 (1996)
!   [b] J. P. Perdew, K. Burke, and Y. Wang, {\sl Real-space cutoff 
!       construction of a generalized gradient approximation:  The PW91 
!       density functional}, submitted to Phys. Rev. B, Feb. 1996. 
!   [c] J. P. Perdew and Y. Wang, Phys. Rev. B {45, 13244 (1992). 
!
!---------------------------------------------------------------
subroutine spn_corpbe(rs,zet,t,uu,vv,ww,lgga,lpot,ec,vcup,vcdn, &
     h,dvcup,dvcdn)
  use constants
  implicit none
  !
  ! Input/Output variables:
  !
  real(dp), intent(in) :: &
       rs, &    ! Seitz radius=(3/4pi rho)^(1/3)
       zet, &   ! relative spin polarization = (rhoup-rhodn)/rho
       t, &     ! ABS(GRAD rho)/(rho*2.*KS*G)  -- only needed for PBE
       uu, &    ! (GRAD rho)*GRAD(ABS(GRAD rho))/(rho**2 * (2*KS*G)**3)
       vv, &    ! (LAPLACIAN rho)/(rho * (2*KS*G)**2)
       ww       ! (GRAD rho)*(GRAD ZET)/(rho * (2*KS*G)**2 
                ! uu,vv,ww, only needed for PBE potential
  integer, intent(in) :: &
       lgga, &  ! flag to do gga (0=>LSD only)
       lpot     ! flag to do potential (0=>energy only)
  real(dp), intent(out) :: &
       ec, &    ! LSD correlation energy from [a]
       vcup, &  ! LSD up correlation potential
       vcdn, &  ! LSD dn correlation potential
       h, &     ! nonlocal part of correlation energy per electron
       dvcup, & ! nonlocal correction to vcup
       dvcdn    ! nonlocal correction to vcdn
  !
  ! Work variables:
  !
  ! thrd*=various multiples of 1/3 
  ! numbers for use in LSD energy spin-interpolation formula, [c](9). 
  !        GAM= 2^(4/3)-2 
  !        FZZ=f''(0)= 8/(9*GAM) 
  ! numbers for construction of PBE 
  !        gamma=(1-log(2))/pi^2 
  !        bet=coefficient in gradient expansion for correlation, [a](4). 
  !        eta=small number to stop d phi/ dzeta from blowing up at 
  !            |zeta|=1. 
  !
  real(dp) :: alfc,alfm,alfrsm,b,b2,bec,bg,comm,ecrs,eczet, &
       ep,eprs,eu,eurs,f,fac,fact0,fact1,fact2,fact3,fact5, &
       fz,g,g3,g4,q5,gz,hb,hbt,hrs,hrst,ht,htt,hz,hzt,pon,pref, &
       q4,q8,q9,rsthrd,rtrs,t2,t4,t6,rs2,rs3,z4

  real(dp), parameter :: thrdm=-third
  real(dp), parameter :: thrd2=two*third
  real(dp), parameter :: sixthm=thrdm/two
  real(dp), parameter :: thrd4=four*third
  real(dp), parameter :: gam = 0.5198420997897463295344212145565d0
  real(dp), parameter :: fzz=8.d0/(9.d0*GAM)
  real(dp), parameter :: gamma = 0.03109069086965489503494086371273d0
  real(dp), parameter :: bet=0.06672455060314922d0
  real(dp), parameter :: delt=bet/gamma
  real(dp), parameter :: eta=1.d-12
  !---------------------------------------------------------------
  !
  ! find LSD energy contributions, using [c](10) and Table I[c]. 
  !   eu=unpolarized LSD correlation energy 
  !   eurs=deu/drs 
  !   ep=fully polarized LSD correlation energy 
  !   eprs=dep/drs 
  !   alfm=-spin stiffness, [c](3). 
  !   alfrsm=-dalpha/drs 
  !   f=spin-scaling factor from [c](9). 
  ! construct ec, using [c](8) 
  !
  rtrs=sqrt(rs)
  call spn_gcor2(0.0310907d0,0.21370d0,7.5957d0,3.5876d0,1.6382d0, &
       0.49294d0,rtrs,eu,eurs)
  call spn_gcor2(0.01554535d0,0.20548d0,14.1189d0,6.1977d0, &
       3.3662d0,0.62517d0,rtrs,ep,eprs)
  call spn_gcor2(0.0168869d0,0.11125d0,10.357d0,3.6231d0, &
       0.88026d0,0.49671d0,rtrs,alfm,alfrsm)
  alfc = -alfm
  z4 = zet**4
  f=((one+zet)**thrd4+(one-zet)**thrd4-two)/gam
  ec = eu*(one-f*z4)+ep*f*z4-alfm*f*(one-z4)/fzz
  !
  ! LSD potential from [c](a1) 
  !      ecrs = dec/drs [c](a2) 
  !      eczet=dec/dzeta [c](a3) 
  !      fz = df/dzeta [c](a4) 
  !
  ecrs = eurs*(one-f*z4)+eprs*f*z4-alfrsm*f*(one-z4)/fzz
  fz = thrd4*((one+zet)**third-(one-zet)**third)/gam
  eczet = four*(zet**3)*f*(ep-eu+alfm/fzz)+fz*(z4*ep-z4*eu &
       -(one-z4)*alfm/fzz)
  comm = ec -rs*ecrs/three-zet*eczet
  vcup = comm + eczet
  vcdn = comm - eczet

  ! note: in pbe.f the following line is missing (adi)
  if(lgga == 0)return
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
  ! note: in pbe.f the following line is missing (adi)
  if(lpot == 0) return
  !
  ! energy done. now the potential, using appendix E of [b]. 
  g4 = g3*g
  t6 = t4*t2
  rsthrd = rs/three
  gz=(((one+zet)**2+eta)**sixthm-((one-zet)**2+eta)**sixthm)/three
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

end subroutine spn_corpbe
!===============================================================
!
! Slimmed down version of gcor used in PW91 routines, to
! interpolate LSD correlation energy, as given by (10) of 
! J. P. Perdew and Y. Wang, Phys. Rev. B 45, 13244 (1992). 
! K. Burke, May 11, 1996. 
!
subroutine spn_gcor2(a,a1,b1,b2,b3,b4,rtrs,gg,ggrs) 

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

  q0 = -two*a*(one+a1*rtrs*rtrs)
  q1 = two*a*rtrs*(b1+rtrs*(b2+rtrs*(b3+b4*rtrs)))
  q2 = log(one+one/q1)
  gg = q0*q2
  q3 = a*(b1/rtrs+two*b2+rtrs*(three*b3+four*b4*rtrs))
  ggrs = -two*a*a1*q2-q0*q3/(q1*(one+q1))

end subroutine spn_gcor2
!===============================================================
!
! gga91 exchange for a spin-unpolarized electronic system 
!
!---------------------------------------------------------------
subroutine exchpw91(d,s,u,v,ex,vx) 

  use constants
  implicit none
  !
  ! Input/Output variables:
  !
  real(dp), intent(in) :: &
       d, &    ! density
       s, &    ! ABS(GRAD d)/(2*kf*d)
       u, &    ! (GRAD d)*GRAD(ABS(GRAD d))/(d**2 * (2*kf)**3)
       v       ! (LAPLACIAN d)/(d*(2*kf)**2)

  real(dp), intent(out) :: &
       ex, &   ! exchange energy per electron
       vx      ! potential
  !
  ! Work variables:
  !
  real(dp) :: f,fac,fs,fss, &
       p0,p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,s2,s3,s4

  real(dp), parameter :: a1=0.19645d0
  real(dp), parameter :: a2=0.27430d0
  real(dp), parameter :: a3=0.15084d0
  real(dp), parameter :: a4=100.d0
  real(dp), parameter :: ax=-0.738558766382022405884230032680836d0
  real(dp), parameter :: a=7.7956d0
  real(dp), parameter :: b1=0.004d0
  real(dp), parameter :: thrd4 = four/three
  !---------------------------------------------------------------
  ! for Becke exchange, set a3=b1=0
  fac = ax*d**third
  s2 = s*s
  s3 = s2*s
  s4 = s2*s2
  p0 = one/sqrt(one+a*a*s2)
  p1 = log(a*s+one/p0)
  p2 = exp(-a4*s2)
  p3 = one/(one+a1*s*p1+b1*s4)
  p4 = one+a1*s*p1+(a2-a3*p2)*s2
  f = p3*p4
  ex = fac*f
  ! local exchange option 
!  ex = fac 
  ! energy done. now the potential: 
  p5 = b1*s2-(a2-a3*p2)
  p6 = a1*s*(p1+a*s*p0)
  p7 = two*(a2-a3*p2)+two*a3*a4*s2*p2-four*b1*s2*f
  fs = p3*(p3*p5*p6+p7)
  p8 = two*s*(b1-a3*a4*p2)
  p9 = a1*p1+a*a1*s*p0*(three-a*a*s2*p0*p0)
  p10 = four*a3*a4*s*p2*(two-a4*s2)-8.d0*b1*s*f-four*b1*s3*fs
  p11 = -p3*p3*(a1*p1+a*a1*s*p0+four*b1*s3)
  fss = p3*p3*(p5*p9+p6*p8)+two*p3*p5*p6*p11+p3*p10+p7*p11
  vx = fac*(thrd4*f-(u-thrd4*s3)*fss-v*fs)
  ! local exchange option: 
!  vx = fac*thrd4 

end subroutine exchpw91
!===============================================================
!
! uniform-gas correlation of Perdew and Wang 1991 
!
!---------------------------------------------------------------
subroutine corlsd(rs,zet,ec,vcup,vcdn,ecrs,eczet,alfc) 

  use constants
  implicit none
  !
  ! Input/Output variables:
  !
  real(dp), intent(in) :: &
       zet, &        ! relative spin polarization
       rs            ! Seitz radius
  real(dp), intent(out) :: &
       ec, &         ! correlation energy per electron
       vcup,vcdn, &  ! up- and down-spin potentials
       ecrs, &       ! derivative of ec w.r.t. rs
       eczet, &      ! derivative of ec w.r.t. zet
       alfc          ! correlation contribution to the spin stiffness
  !
  ! Work variables:
  !
  real(dp) :: alfm,alfrsm,comm,ep,eprs,eu,eurs,f,z4,fz

  real(dp), parameter :: gam=0.5198421d0
  real(dp), parameter :: fzz=1.709921d0
  real(dp), parameter :: thrd4=four/three
  !---------------------------------------------------------------
  f = ((one+zet)**thrd4+(one-zet)**thrd4-two)/gam
  call spn_gcor(0.0310907d0,0.21370d0,7.5957d0,3.5876d0,1.6382d0, &
       0.49294d0,1.00d0,rs,eu,eurs)
  call spn_gcor(0.01554535d0,0.20548d0,14.1189d0,6.1977d0, &
       3.3662d0,0.62517d0,1.00d0,rs,ep,eprs)
  call spn_gcor(0.0168869d0,0.11125d0,10.357d0,3.6231d0,0.88026d0, &
       0.49671d0,1.00d0,rs,alfm,alfrsm)
  ! alfm is minus the spin stiffness alfc 
  alfc = -alfm
  z4 = zet**4
  ec = eu*(one-f*z4)+ep*f*z4-alfm*f*(one-z4)/fzz
  ! energy done. now the potential: 
  ecrs = eurs*(one-f*z4)+eprs*f*z4-alfrsm*f*(one-z4)/fzz
  fz = thrd4*((one+zet)**third-(one-zet)**third)/gam
  eczet = four*(zet**3)*f*(ep-eu+alfm/fzz)+fz*(z4*ep-z4*eu &
       -(one-z4)*alfm/fzz)
  comm = ec -rs*ecrs/three-zet*eczet
  vcup = comm + eczet
  vcdn = comm - eczet
  
end subroutine corlsd
!===============================================================
!
! This subroutine computes the local correlation energy and
! potential for the Perdew-Wang exchange-correlation scheme.
!
!---------------------------------------------------------------
subroutine spn_gcor(a,a1,b1,b2,b3,b4,p,rs,gg,ggrs) 

  use constants
  implicit none
  !
  ! Input/Output variables:
  !
  real(dp), intent(in) :: a,a1,b1,b2,b3,b4,p,rs
  real(dp), intent(out) :: gg,ggrs
  !
  ! Work variables:
  !
  real(dp) :: p1,q0,rsp,q1,q2,q3,rs12,rs32
  !---------------------------------------------------------------
  p1 = p + 1.d0
  q0 = -two*a*(one+a1*rs)
  rs12 = sqrt(rs)
  rs32 = rs12**3
  rsp = rs**p
  q1 = two*a*(b1*rs12+b2*rs+b3*rs32+b4*rs*rsp)
  q2 = log(one+one/q1)
  gg = q0*q2
  q3 = a*(b1/rs12+two*b2+three*b3*rs12+two*b4*p1*rsp)
  ggrs = -two*a*a1*q2-q0*q3/(q1**2+q1)

end subroutine spn_gcor
!===============================================================
!
! PW91 correlation, modified by K. Burke to put all arguments 
! as variables in calling statement, rather than in common block 
! May, 1996. 
!
!---------------------------------------------------------------
subroutine corpw91(rs,zet,g,ec,ecrs,eczet,t,uu,vv,ww,h,dvcup,dvcdn)

  use constants
  implicit none
  !
  ! Input/Output variables:
  !
  real(dp), intent(in) :: &
       rs, &        ! Seitz radius
       zet, &       ! relative spin polarization
       t, &         ! ABS(GRAD d)/(d*2.*ks*g)
       uu, &        ! (GRAD d)*GRAD(ABS(GRAD d))/(d**2 * (2*ks*g)**3)
       vv, &        ! (LAPLACIAN d)/(d * (2*ks*g)**2)
       ww, &        ! (GRAD d)*(GRAD zet)/(d * (2*ks*g)**2
       g,ec,ecrs,eczet

  real(dp), intent(out) :: &
       h, &         ! nonlocal part of correlation energy per el.
       dvcup,dvcdn  ! nonlocal parts of correlation potentials
  !
  ! Work variables:
  !
  real(dp) :: b,b2,bec,bet,bg,cc,ccrs,coeff,comm,delt, &
       fac,fact0,fact1,fact2,fact3,fact4,fact5,g3,g4,gz,h0b,h0bt, &
       h0rs,h0rst,h0t,h0tt,h0z,h0zt,h1,h1rs,h1rst,h1t,h1tt,h1z,h1zt, &
       hrs,hrst,ht,htt,hz,hzt,pon,pref,rs2,rs3,q4,q5,q6,q7,q8,q9, &
       r0,r1,h0,r2,r3,r4,rsthrd,t2,t4,t6

  real(dp), parameter :: xnu=15.75592d0
  real(dp), parameter :: cc0=0.004235d0
  real(dp), parameter :: cx=-0.001667212d0
  real(dp), parameter :: alf=0.09d0
  real(dp), parameter :: c1=0.002568d0
  real(dp), parameter :: c2=0.023266d0
  real(dp), parameter :: c3=7.389d-6
  real(dp), parameter :: c4=8.723d0
  real(dp), parameter :: c5=0.472d0
  real(dp), parameter :: c6=7.389d-2
  real(dp), parameter :: a4=100.d0
  real(dp), parameter :: thrdm=-third
  real(dp), parameter :: thrd2=two/three
  !---------------------------------------------------------------
  bet = xnu*cc0
  delt = two*alf/bet
  g3 = g**3
  g4 = g3*g
  pon = -delt*ec/(g3*bet)
  b = delt/(exp(pon)-one)
  b2 = b*b
  t2 = t*t
  t4 = t2*t2
  t6 = t4*t2
  rs2 = rs*rs
  rs3 = rs2*rs
  q4 = one+b*t2
  q5 = one+b*t2+b2*t4
  q6 = c1+c2*rs+c3*rs2
  q7 = one+c4*rs+c5*rs2+c6*rs3
  cc = -cx + q6/q7
  r0 = 0.663436444d0*rs
  r1 = a4*r0*g4
  coeff = cc-cc0-three*cx/7.d0
  r2 = xnu*coeff*g3
  r3 = exp(-r1*t2)
  h0 = g3*(bet/delt)*log(one+delt*q4*t2/q5)
  h1 = r3*r2*t2
  h = h0+h1
  ! local correlation option: 
!  h = 0.0d0 
  ! energy done. now the potential: 
  ccrs = (c2+2.*c3*rs)/q7 - q6*(c4+2.*c5*rs+3.*c6*rs2)/q7**2
  rsthrd = rs/three
  r4 = rsthrd*ccrs/coeff
  gz = ((one+zet)**thrdm - (one-zet)**thrdm)/three
  fac = delt/b+one
  bg = -three*b2*ec*fac/(bet*g4)
  bec = b2*fac/(bet*g3)
  q8 = q5*q5+delt*q4*q5*t2
  q9 = one+two*b*t2
  h0b = -bet*g3*b*t6*(two+b*t2)/q8
  h0rs = -rsthrd*h0b*bec*ecrs
  fact0 = two*delt-6.d0*b
  fact1 = q5*q9+q4*q9*q9
  h0bt = two*bet*g3*t4*((q4*q5*fact0-delt*fact1)/q8)/q8
  h0rst = rsthrd*t2*h0bt*bec*ecrs
  h0z = three*gz*h0/g + h0b*(bg*gz+bec*eczet)
  h0t = two*bet*g3*q9/q8
  h0zt = three*gz*h0t/g+h0bt*(bg*gz+bec*eczet)
  fact2 = q4*q5+b*t2*(q4*q9+q5)
  fact3 = two*b*q5*q9+delt*fact2
  h0tt = four*bet*g3*t*(two*b/q8-(q9*fact3/q8)/q8)
  h1rs = r3*r2*t2*(-r4+r1*t2/three)
  fact4 = two-r1*t2
  h1rst = r3*r2*t2*(two*r4*(one-r1*t2)-thrd2*r1*t2*fact4)
  h1z = gz*r3*r2*t2*(three-four*r1*t2)/g
  h1t = two*r3*r2*(one-r1*t2)
  h1zt = two*gz*r3*r2*(three-11.d0*r1*t2+four*r1*r1*t4)/g
  h1tt = two*r3*r2*r1*t*(-two+r1*t2)
  hrs = h0rs+h1rs
  hrst = h0rst+h1rst
  ht = h0t+h1t
  htt = h0tt+h1tt
  hz = h0z+h1z
  hzt = h0zt+h1zt
  comm = h+hrs+hrst+t2*ht/6.d0+7.d0*t2*t*htt/6.d0
  pref = hz-gz*t2*ht/g
  fact5 = gz*(two*ht+t*htt)/g
  comm = comm-pref*zet-uu*htt-vv*ht-ww*(hzt-fact5)
  dvcup = comm + pref
  dvcdn = comm - pref
! local correlation option: 
! dvcup = 0.0d0 
! dvcdn = 0.0d0 

end subroutine corpw91
!===============================================================
