!=============================================================================!
! This subroutine computes a non-additive kinetic energy term
! appears in the Frozen Density Embedding (FDE) method.
! See JPC 97, 8050 (1992), for example.
! Part of the AFM project
!                             Jul. 2015, Yuki Sakai (yuki@ices.utexas.edu)
!=============================================================================!
subroutine kinetic_energy(grid,parallel,pot,weight,rhotot,tn,vn,wvec)
  use grid_module
  use parallel_data_module
  use potential_module
  implicit none
  !
  ! Input/Output variables:
  !
  ! grid related data
  type (grid_data), intent(in) :: grid
  ! parallel computation related data
  type (parallel_data), intent(in) :: parallel
  ! potential
  type (potential), intent(in) :: pot

  ! Work variables:
  !
  ! multiplicity of grid points in irreducible wedge (equal to
  ! order of the reduced group)
  real(dp), intent(in) :: weight

  ! total charge density 
  real(dp), intent(in) :: rhotot(parallel%mydim)
  ! work array
  real(dp), intent(out) :: wvec(parallel%nwedge + 1)
  ! non-addotove kinetic energy 
  real(dp), intent(out) :: tn
  ! non-additive kinetic potential
  real(dp), intent(out) :: vn(parallel%mydim)
  !
  !Temp
  real(dp) :: ekin
  real(dp) :: vkin(parallel%mydim)
  real(dp) :: rhot(parallel%mydim)
  ! counters
  integer i

  ! initialize
  vn(:) = zero  
  tn = zero

  !total charge density
  call compute_kin(grid,parallel,pot,weight,rhotot,ekin,vkin,wvec)
  do i = 1, parallel%mydim
     vn(i) = vn(i) + vkin(i)
  end do
  tn = tn + ekin
  
  !sample
  call compute_kin(grid,parallel,pot,weight,pot%rho0,ekin,vkin,wvec)
  tn = tn - ekin

  !tip
  do i = 1, parallel%mydim
     rhot(i) = rhotot(i) - pot%rho0(i)
  end do
  call compute_kin(grid,parallel,pot,weight,rhot,ekin,vkin,wvec)
  do i = 1, parallel%mydim
     vn(i) = vn(i) - vkin(i)
  end do
  tn = tn - ekin

end subroutine kinetic_energy

!===============================================================
!
! This subroutine computes an approximated kinetic energy.
!
!---------------------------------------------------------------
subroutine compute_kin(grid,parallel,pot,weight,rho,ekin,vkin,wv)
  use constants
  use grid_module
  use parallel_data_module
  use potential_module
  implicit none
  ! grid related data
  type (grid_data), intent(in) :: grid
  ! parallel computation related data
  type (parallel_data), intent(in) :: parallel
  ! potential
  type (potential), intent(in) :: pot
  real(dp), intent(in) :: rho(parallel%mydim)
  real(dp), intent(in) :: weight
  real(dp), intent(out) :: ekin
  real(dp), intent(out) :: vkin(parallel%mydim)
  ! GGA variables:
  ! quantities fed to pbe.f:
  ! gradient - grad(n), where n is the charge density
  ! absgrad  - |grad(n)|
  ! gradgrad - components of grad(|grad(n)|
  ! laplac - laplacian(n)
  ! prod  -  grad(n) dot grad(|grad(n)|
  real(dp), dimension(:), allocatable :: absgrad,laplac
  real(dp), dimension(:,:), allocatable :: gradient,gradgrad
  real(dp) :: prod
  ! counters
  integer i
  ! temporary variables
  real(dp), dimension (3) :: temp1,temp2
  real(dp) :: ek, vk
  !Thomas-Fermi coefficient 
  real(dp) :: CTF
  ! work array
  real(dp), intent(out) :: wv(parallel%nwedge + 1)
  ! allocation check
  integer alcstat
  ! temporary parmeters
  real(dp) :: fk, u, v, s, pi32
  integer :: nd
  
  nd = parallel%mydim
  ekin = zero
  vkin(:) = zero

  ! calculate the gradient and the Laplacian of the charge density
  ! at each grid point
  allocate(gradient(3,nd),stat=alcstat)
  call alccheck('gradient',3*nd,alcstat)
  allocate(absgrad(nd),stat=alcstat)
  call alccheck('absgrad',nd,alcstat)
  allocate(gradgrad(3,nd),stat=alcstat)
  call alccheck('gradgrad',3*nd,alcstat)
  allocate(laplac(nd),stat=alcstat)
  call alccheck('laplac',nd,alcstat)

  !Thomas-Fermi Constant
  CTF =  0.6d0 * (three * (pi ** 2)) ** (two / three)
  !
  pi32 = three*pi*pi

  call lapmvs(parallel,rho,laplac,grid%coe2,grid%norder, &
       grid%lap_dir_num,wv)
  laplac = -laplac
  call gradmvs(grid,parallel,rho,gradient,grid%coe1,grid%norder, &
        wv)
  do i = 1,nd
     temp1 = gradient(:,i)
     absgrad(i) = sqrt(dot_product(temp1,temp1))
  enddo

  ! calculate other input quantities and then calculate the 
  ! exchange-correlation

  ! calculate grad(rho)_dot_grad|grad(rho)| and
  ! calculate the exchange-correlation value using pbe
  call gradmvs(grid,parallel,absgrad,gradgrad,grid%coe1, &
       grid%norder,wv)

  do i = 1, nd
     temp1 = gradient(:,i)
     temp2 = gradgrad(:,i)
     prod = dot_product(temp1,temp2)
     fk=(pi32*rho(i))**third
     s=absgrad(i)/(two*fk*rho(i))
     u=two*prod/(rho(i)*rho(i)*(two*fk)**3)
     v=laplac(i)/(rho(i)*(two*fk)**2)
     if(pot%kin_name == 1) then
        call pw91k(rho(i), s, u, v, ek, vk)
     else if(pot%kin_name == 2) then 
        call pbek(rho(i), s, u, v, ek, vk, pot%vum, pot%vuk)
     else if(pot%kin_name == 3) then
        vk=(five/three)*(rho(i)**(two/three))
        ek=rho(i)**(two/three)
     else if(pot%kin_name == 4) then
        call tfwk(rho(i),absgrad(i),laplac(i),ek,vk)
     else if(pot%kin_name ==0) then
        vk = zero
        ek = zero
     else
        write(7,*) "Invalid KE potential"
        stop 612
     end if
     vkin(i) = CTF * vk 
     ekin = ekin + CTF * rho(i) * ek * weight 
  end do
  deallocate(gradient,absgrad,gradgrad,laplac)
end subroutine compute_kin
!===============================================================
!
! pw91k kinetic energy
!
!---------------------------------------------------------------
subroutine pw91k(d,s,u,v,ek,vk)

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
       ek, &   ! exchange energy per electron
       vk      ! potential
  !
  ! Work variables:
  !
  real(dp) :: f,fac,fs,fss, &
       p0,p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,s2,s3,s4

  real(dp), parameter :: a1=0.093907d0
  real(dp), parameter :: a2=0.26608d0
  real(dp), parameter :: a3=0.0809615d0
  real(dp), parameter :: a4=100.d0
!  real(dp), parameter :: ax=-0.7385588d0
  real(dp), parameter :: a=76.32d0
  real(dp), parameter :: b1=0.000057767d0
  real(dp), parameter :: thrd4 = four/three
  !---------------------------------------------------------------
  ! for Becke exchange, set a3=b1=0
  if (d < 1.0d-18) then
     ek = 0
     vk = 0
  else
     fac = d**(two/three)
     s2 = s*s
     s3 = s2*s
     s4 = s2*s2
     p0 = one/sqrt(one+a*a*s2)
     p1 = log(a*s+one/p0)
     p2 = exp(-a4*s2)
     p3 = one/(one+a1*s*p1+b1*s4)
     p4 = one+a1*s*p1+(a2-a3*p2)*s2
     f = p3*p4
     ek = fac*f
     ! local exchange option
     !  ex = fac
     p5 = b1*s2-(a2-a3*p2)
     p6 = a1*s*(p1+a*s*p0)
     p7 = two*(a2-a3*p2)+two*a3*a4*s2*p2-four*b1*s2*f
     fs = p3*(p3*p5*p6+p7)
     p8 = two*s*(b1-a3*a4*p2)
     p9 = a1*p1+a*a1*s*p0*(three-a*a*s2*p0*p0)
     p10 = four*a3*a4*s*p2*(two-a4*s2)-8.d0*b1*s*f-four*b1*s3*fs
     p11 = -p3*p3*(a1*p1+a*a1*s*p0+four*b1*s3)
     fss = p3*p3*(p5*p9+p6*p8)+two*p3*p5*p6*p11+p3*p10+p7*p11
     vk = fac*((five/three)*f-(u-thrd4*s3)*fss-v*fs-fs*s2/three)
     ! local exchange option:
     !  vx = fac*five/three
  end if

end subroutine pw91k
!===============================================================
!PBE-TW functional
!Weizsacker correction
!---------------------------------------------------------------
subroutine pbek(rho,s,u,v,ek,vk,vum,vuk)

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
       ek, &   ! exchange energy per electron
       vk      ! potential
  real(dp), intent(in) :: vum, vuk
  !
  ! Work variables:
  !
  real(dp) :: exunif,fs,fss,fxpbe,p0,s2

  real(dp), parameter :: thrd4=four/three
  real(dp) :: ul
  !variable parameters
  !---------------------------------------------------------------
  !
  ! construct lda kinetic energy density
  exunif =rho**(two/three)
  !
  ! ul
  ul = vum / vuk
  if (rho < 1.0d-18) then
     ek = 0
     vk = 0
  else
     ! construct PBE enhancement factor
     s2 = s*s
     p0=one+ul*s2
     fxpbe = one+vuk-vuk/p0
     ek = exunif*fxpbe
     !
     ! energy done. now the potential:
     ! find first and second derivatives of fx w.r.t s.
     !      fs=(1/s)*d fxpbe/ ds
     !      fss=d fs/ds
     !
     fs=two*vuk*ul/(p0*p0)
     fss=-four*ul*s*fs/p0
     !
     ! calculate potential from [b](24)
     !vx = exunif*((five/three)*fxpbe-(u-thrd4*s2*s)*fss-v*fs)
     vk = exunif*((five/three)*fxpbe-(u-thrd4*s2*s)*fss-v*fs-fs*s2/three)
  end if
end subroutine pbek
!===============================================================
!
! TFW functional
! Does not work currently
!
!---------------------------------------------------------------
subroutine tfwk(rho,grad,lap,ek,vk)
  use constants
  implicit none
  real(dp), intent(in)  :: rho
  real(dp), intent(in)  :: grad
  real(dp), intent(in)  :: lap
  real(dp), intent(out) :: vk
  real(dp), intent(out) :: ek
  real(dp) :: v72
  real(dp) :: CTF

  v72 = one / 72.d0
  CTF =  0.6d0 * (three * (pi ** 2)) ** (two / three)

  if (rho < 1.0d-18) then
     ek = 0
     vk = 0
  else
     ek = rho**(two/three) + (v72*(grad*grad)/(rho*rho)) / CTF
     vk = (five/three)*rho**(two/three) + &
         (v72*(grad*grad/(rho*rho) - two*lap/rho)) /CTF
  end if

end subroutine tfwk
!===============================================================
