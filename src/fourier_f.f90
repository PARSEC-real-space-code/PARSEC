!===============================================================
!
! Copyright (C) 2005 Finite Difference Research Group
! This file is part of parsec, http://www.ices.utexas.edu/parsec/
!
! Performs Fourier filtering in the pseudopotentials (local and
! non-local components) and of core charge density for all atomic
! species according to Briggs et al., PRB 54, 14362 (1996). At
! the end, the input functions (read from pseudopotential files)
! are replaced with filtered functions. Fourier filtering is done
! by performing one-dimensional Fourier transforms with the
! appropriate cut-off parameters. Transforms are calculated using
! an eight point Newton-Cotes integration method. See Abramowitz
! and Stegun Eq. 25.4.17.
!
! If periodic boundary conditions are used, the Fourier transform
! of input (not filtered) pseudopotentials and core charge
! density are saved in structure pbc. They will be used later for
! FFTs.
!
! Adapted from pseudopotential programs written by S. Froyen and
! J. L. Martins.
!
!---------------------------------------------------------------
subroutine fourier_f(clust,p_pot,pbc,step,ipr)

  use constants
  use cluster_module
  use pseudo_potential_module
  use pbc_module
  implicit none
  !
  ! Input/Output variables:
  !
  ! the cluster
  type (cluster), intent(in) :: clust
  ! pseudo_potential structure
  type (pseudo_potential), intent(inout) :: p_pot
  ! periodic boundary conditions structure
  type (pbc_data), intent(inout) :: pbc 
  ! input grid spacing, used to define the range of reciprocal space
  ! AJB: but should we just use a single grid step even if doing PBC?
  real(dp), intent(in) :: step
  ! output flag
  integer, intent(in) :: ipr
  !
  ! Work variables:
  !
  ! counters
  integer ity, ii, j, ll, nr, nql
  real(dp), dimension(p_pot%mxpot) :: &
       vr, &       ! function in real space
       rr, &       ! logarithmic grid in real space
       rab         ! differential of grid, for numerical integration
  ! number of points in 1-D Fourier space mesh (Fourier filtering)
  integer, dimension(p_pot%type_num) :: nq
  ! maxnq = maxval(nq)
  integer :: maxnq
  !
  ! qq - coordinates of mesh points in reciprocal space
  ! fvion - Fourier-like transform of vion
  ! fdenc - Fourier-like transform of denc
  ! frho - Fourier-like transform of rho_r
  ! fvw - Fourier-like transform of vw
  real(dp), dimension (:,:), allocatable :: qq, fvion, fdenc, frho
  real(dp), dimension (:,:,:), allocatable :: fvw

  ! temporary variables
  real(dp), allocatable :: vq(:)
  real(dp) :: vt, xq, alpht, ztot
  logical :: fflag(p_pot%type_num)

  ! spacing of 1-D Fourier space mesh  (Fourier filtering)
  real(dp), parameter :: delq = 0.015d0
  !real(dp), parameter :: delq = 0.010d0
  !
  ! External function:
  real(dp), external :: fint,dfint

  !---------------------------------------------------------------
  !
  ! Fourier filtering is done if:
  ! (1) at least one alpha value is initialized to positive value, or
  ! (2) at least one acore value is initialized to positive value.
  ! (3) there are pbc fourier based stuff here that should be moved
  !
  fflag = .false.
  do ity = 1, p_pot%type_num
     if (p_pot%alpha(ity) > zero) fflag(ity) = .true.
     if (p_pot%acore(ity) > zero) fflag(ity) = .true.
  enddo
  !
  if ( all(.not. fflag ) .and. .not. pbc%is_on ) return

  write(7,*) 'Performing Fourier transform of pseudopotentials:'

  ! start loop over atomic types
  ! range of reciprocal space for Fourier filtering
  do ity = 1, p_pot%type_num
     nq(ity) = 0
     if ( .not. fflag(ity) .and. .not. pbc%is_on ) cycle

     xq = max(one,p_pot%alpha(ity),p_pot%acore(ity)) * pi / step
     nq(ity) = 2*nint( xq/delq ) + 4
  enddo
  maxnq = maxval( nq )
  allocate (qq(maxnq,p_pot%type_num))
  qq(:,:) = zero
  allocate (fdenc(maxnq,p_pot%type_num))
  fdenc(:,:) = zero
  allocate (frho(maxnq,p_pot%type_num))
  frho(:,:) = zero
  allocate (fvion(maxnq,p_pot%type_num))
  fvion(:,:) = zero
  allocate (fvw(maxnq,p_pot%type_num,4))
  fvw(:,:,:) = zero

  do ity = 1, p_pot%type_num
     if ( .not. fflag(ity) .and. .not. pbc%is_on ) cycle

     nr = p_pot%ns(ity)
     nql = nq(ity)

     write(7,*)
     write(7,21) clust%name(ity),nql,delq
     if (p_pot%alpha(ity) > zero) write(7,*) &
          'Actually also performing fourier filtering in pseudopotentials '
     if (p_pot%acore(ity) > zero) write(7,*) &
          'Actually also performing fourier filtering in core charge '
     if (fflag(ity)) write(7,'(3(2x,a,f10.4))') &
          'alpha = ',p_pot%alpha(ity), &
          'beta1 = ',p_pot%beta1(ity),'acore = ',p_pot%acore(ity)
     do j = 2, nql
        qq(j,ity) = delq * real(j-1,dp)
     enddo

     allocate(vq(nql))
     rr(:) = zero
     rr(1:nr) = p_pot%rs(1:nr,ity)
     do ii=1,nr
        rab(ii) = (rr(ii)+p_pot%par_c(ity))*p_pot%par_b(ity)
     enddo
     !
     ! Fourier transform of local potential
     ! fvion is continuous (and non-zero) at q = zero
     vr(:) = zero
     vr(1:nr) = p_pot%vion(1:nr,ity)
     vt = two*p_pot%zion(ity)
     call fourier_1d(nr,rr,rab,vr,0,nql,qq(1,ity),vq,vt)
     do j = 2, nql
        fvion(j,ity) = vq(j)
     enddo
     fvion(1,ity) = -four*p_pot%zion(ity)/pi
     !
     ! Fourier transform of core charge density
     ! fdenc = 0 at q = zero
     if (p_pot%icore(ity) /= 0) then
        vr(:) = zero
        do ii = 2, nr
           vr(ii) = p_pot%denc(ii,ity)
        enddo
        vt = zero
        call fourier_1d(nr,rr,rab,vr,0,nql,qq(1,ity),vq,vt)
        do j = 2, nql
           fdenc(j,ity) = vq(j)
        enddo
        fdenc(1,ity) = zero
     endif
     !
     ! Fourier transform of valence charge density
     ! frho = 0 at q = zero
     vr(:) = zero
     do ii = 2, nr
        vr(ii) = p_pot%rho_r(ii,ity)
     enddo
     vt = zero
     call fourier_1d(nr,rr,rab,vr,0,nql,qq(1,ity),vq,vt)
     do j = 2, nql
        frho(j,ity) = vq(j)
     enddo
     frho(1,ity) = zero
     !
     ! Fourier transform of non-local components
     ! fvw = 0 at q = zero
     do ll = 1, p_pot%nlocp(ity)
        if (ll == p_pot%loc(ity)) cycle
        vr(:) = zero
        vr(1:nr) = p_pot%vw(1:nr,ity,ll)
        vt = zero
        call fourier_1d(nr,rr,rab,vr,ll-1,nql,qq(1,ity),vq,vt)
        do j = 2, nql
           fvw(j,ity,ll) = vq(j)
        enddo
     enddo
     fvw(1,ity,:) = zero

     deallocate(vq)
     if (ipr >= 1) then
        open(20,file='fourier_'//trim(clust%name(ity))//'.dat', &
             form='formatted')
        do j = 1, nql
           write(20,'(10g16.4)') qq(j,ity),fvion(j,ity), &
                fdenc(j,ity),frho(j,ity),(fvw(j,ity,ii),ii=1,4)
        enddo
        close(20)
    endif
     if (ipr >= 0) then
        write(7,*) ' Pseudopotential comparison tests: '
        write(7,22) ' local potential V(0)*r [Ry] = ', &
             sum(fvion(:,ity))*delq &
             *p_pot%rs(2,ity),p_pot%vion(2,ity)*p_pot%rs(2,ity)
        write(7,22) ' V(q=0)/V_coul = ',fvion(1,ity)*pi/four &
             /p_pot%zion(ity),-one
     endif
     if (ipr >= 1) then
        write(7,22) ' alpha_crystal = ',(fvion(2,ity) -  &
             fvion(1,ity))*two*pi*pi/delq**2
        write(7,22) ' integrated core charge ',fdenc(2,ity) &
             *two*pi*pi/qq(2,ity)**2
        write(7,22) ' rho_core(0) = ',sum(fdenc(:,ity)) &
             *delq,p_pot%denc(2,ity)
        write(7,22) ' integrated valence charge ',frho(2,ity) &
             *two*pi*pi/qq(2,ity)**2
        write(7,22) ' rho_val(0) = ',sum(frho(:,ity)) &
             *delq,p_pot%rho_r(2,ity)
        write(7,22) ' non-local potentials, zero value: '
        write(7,23) (fvw(2,ity,ll+1)/qq(2,ity)**2 &
             *exp(-real(ll,dp)*log(qq(2,ity))),ll=0,3)
     endif
     !
     ! Replace old functions with filtered /= and impose the
     ! correct boundary conditions.
     !
     if (p_pot%alpha(ity) > zero) then
        do j = 2, nr
           p_pot%rho_r(j,ity) = fint(rr(j),nq(ity),delq, &
                frho(1,ity),0,p_pot%alpha(ity),p_pot%beta1(ity),step)
           p_pot%drhodr(j,ity) = dfint(rr(j),nq(ity),delq, &
                frho(1,ity),0,p_pot%alpha(ity),p_pot%beta1(ity),step)
           p_pot%vion(j,ity) = fint(rr(j),nq(ity),delq, &
                fvion(1,ity),0,p_pot%alpha(ity),p_pot%beta1(ity),step)
           p_pot%dvion(j,ity) = dfint(rr(j),nq(ity),delq, &
                fvion(1,ity),0,p_pot%alpha(ity),p_pot%beta1(ity),step)
           do ll = 1, p_pot%nlocp(ity)
              p_pot%vw(j,ity,ll) = fint(rr(j),nq(ity),delq,fvw(1,ity,ll), &
                   ll-1,p_pot%alpha(ity),p_pot%beta1(ity),step)
              p_pot%dvw(j,ity,ll) = fint(rr(j),nq(ity),delq,fvw(1,ity,ll), &
                   ll-1,p_pot%alpha(ity),p_pot%beta1(ity),step)
           enddo
        enddo
        do j = -p_pot%norder, 1
           p_pot%rho_r(j,ity) = p_pot%rho_r(2,ity)
           p_pot%drhodr(j,ity) = zero
           p_pot%vion(j,ity) = p_pot%vion(2,ity)
           p_pot%dvion(j,ity) = zero
           do ll = 1, p_pot%nlocp(ity)
              p_pot%vw(j,ity,ll) = zero
              p_pot%dvw(j,ity,ll) = zero
           enddo
        enddo
     endif
     if (p_pot%acore(ity) > zero .and.  p_pot%icore(ity) /= 0) then
        do j = 2, nr
           p_pot%denc(j,ity) = fint(rr(j),nq(ity),delq, &
                fdenc(1,ity),0,p_pot%acore(ity),zero,step)
           p_pot%ddenc(j,ity) = dfint(rr(j),nq(ity),delq, &
                fdenc(1,ity),0,p_pot%acore(ity),zero,step)
        enddo
        do j = -p_pot%norder, 1
           p_pot%denc(j,ity) = p_pot%denc(2,ity)
           p_pot%ddenc(j,ity) = zero
        enddo
     endif
  enddo                     ! ity = 1, p_pot%type_num
  !
  ! Define Fourier transforms of local pseudopotential
  ! components and of core charge density
  !
  if (pbc%is_on) then
     pbc%delq = delq
     do ity = 1, pbc%type_num
        pbc%nq(ity) = int(two*pi/step/delq)+4
     enddo
     maxnq = maxval (pbc%nq)
     !actually alocate vloc and dcor
     call pbc_set_maxdlqp (maxnq, pbc)
     pbc%vloc(:,:) = zero
     pbc%dcor(:,:) = zero

     ! start loop over atomic types
     do ity=1,pbc%type_num
        nql = pbc%nq(ity) - 2

        ! Convert local potential to hartree and save in structure.
        ! Notice that the first few points are padded (for interpolation).
        do ii = 2, nql
           pbc%vloc(ii+1,ity) = fvion(ii,ity) * pi * pi / qq(ii,ity)**2
        enddo
        pbc%vloc(1,ity) = pbc%vloc(3,ity)
        pbc%vloc(2,ity) = pbc%vloc(3,ity)

        ! Core charge density is stored as volume density.
        if (p_pot%icore(ity) /= 0) then
           do ii = 2, nql
              pbc%dcor(ii+1,ity) = fdenc(ii,ity) * pi * pi * two &
                   / qq(ii,ity)**2
           enddo
           pbc%dcor(1,ity) = pbc%dcor(3,ity)
           pbc%dcor(2,ity) = pbc%dcor(3,ity)
        endif
     enddo                     ! ity = 1, pbc%type_num

     ! calculate alpha energy
     alpht = zero
     ztot = zero
     do ity = 1, pbc%type_num
        ztot = ztot + p_pot%zion(ity)*real(clust%natmi(ity),dp)
        alpht = alpht + two*pi*pi*real(clust%natmi(ity),dp)* &
             (fvion(2,ity) - fvion(1,ity))/delq**2
     enddo
     pbc%ealpha = alpht*ztot
  endif

  deallocate(qq,fvion,fdenc,frho,fvw)

21 format(1x,a2,5x,'nql=',i5,5x,'delq=',f6.3)
22 format(1x,a,6(1x,g11.4))
23 format(1x,6(1x,g11.4))

end subroutine fourier_f
!===============================================================
!
! For an input function vr in real space, computed on a logarithmic
! grid rr(1:nr), calculates the function vql defined as:
!
! vql(q) = 2 * q^2 / pi * int_0^infty r^2 vr(r) * J(ll,q * r) dr
!
! where J is the spherical Bessel function of first kind,
! J(0,x) = sin(x)/x.
!
! If vr has a Coulomb component (i.e., vr(r) = -vt/r for large r),
! the divergence is integrated out explicitly. That happens for
! the local component of norm-conserving pseudopotentials.
!
!---------------------------------------------------------------
subroutine fourier_1d(nr,rr,rab,vr,ll,nql,yp,vql,vt)

  use constants
  implicit none
  !
  ! Input/Output variables:
  !
  ! number of grid points
  integer, intent(in) :: nr
  ! coordinates of grid points and differential, dr
  real(dp), intent(in) :: rr(nr), rab(nr)
  ! input function in real space
  real(dp), intent(in) :: vr(nr)
  ! order of Bessel function
  integer, intent(in) :: ll
  ! number of reciprocal space points
  integer, intent(in) :: nql
  ! coordinates of reciprocal space points
  real(dp), intent(in) :: yp(nql)
  ! output function
  real(dp), intent(out) :: vql(nql)
  ! compensating charge
  real(dp), intent(in) :: vt
  !
  ! Work variables:
  !
  ! counters
  integer :: ii, j
  real(dp), dimension(nr+4) :: &
       vtmp, &    ! input function minus compensating function
       y          ! integrand in 1-D Fourier transform
  !
  ! External function:
  real(dp), external :: besselj

  y(:) = zero
  vtmp(:) = zero
  vql(:) = zero
  do ii = 2, nr
     vtmp(ii) = rab(ii)*rr(ii)*(rr(ii)*vr(ii)+vt)
  enddo
  do j = 2, nql
     do ii = 2, nr
        y(ii) = vtmp(ii)*besselj(ll,yp(j)*rr(ii))
     enddo
     do ii = 1, nr, 4
        vql(j) = vql(j) + 7.0*y(ii) + 32.0*y(ii+1) + &
             12.0*y(ii+2) + 32.0*y(ii+3) + 7.0*y(ii+4)
     enddo
     vql(j) = (two*yp(j)*yp(j)*vql(j)/45.0 - vt)*two/pi
  enddo

end subroutine fourier_1d
!===============================================================
!
! For an input function vql in reciprocal space, perform Fourier
! filtering and calculate fint at a given point in real space
! according to:
!
! fint = int_0^infty vql(q) * F_filter(q) * J(ll,q*r) dq
!
! with F_filter(q) = 1                                  q < q0 
!                  = exp( -beta1 * (q/q0 - 1)^2 )       q > q0
! q0 = alpha*pi/h
!
! The filter function was proposed by Briggs et al., PRB 54,
! 14362 (1996). Input function vql is tabulated on a regular grid
! with nql points and spacing delql. First value vql(1) corresponds
! to q = 0.
!
!---------------------------------------------------------------
function fint(r,nql,delql,vql,ll,alpha,beta1,h)

  use constants
  implicit none
  !
  ! Input/Output variables:
  !
  ! radial position at which vr is to be calculated
  real(dp), intent(in) :: r
  ! number of reciprocal space points
  integer, intent(in) :: nql
  ! grid spacing
  real(dp), intent(in) :: delql
  ! input function in reciprocal space
  real(dp), intent(in) :: vql(nql)
  ! order of Bessel function
  integer, intent(in) :: ll
  ! parameters for filter function
  real(dp), intent(in) :: alpha, beta1, h
  ! filtered function
  real(dp) :: fint
  !
  ! Work variables:
  !
  ! counters
  integer :: ii, nfilter
  real(dp) :: xx, ff, q0d, qq0, qq0max, fsum
  ! integrand in Newton-Cottes integration
  real(dp), dimension(nql) :: y
  !
  ! External function:
  real(dp), external :: besselj

  y = zero
  q0d = alpha*pi/h/delql
  nfilter = int( q0d ) + 1
  if (beta1 > zero) then
     qq0max = sqrt(10.d0/beta1) + one
  else
     qq0max = one
  endif

  do ii = 1, nfilter
     xx = r * delql * real(ii-1,dp)
     y(ii) = vql(ii)*besselj(ll,xx)
  enddo

  do ii = nfilter + 1, nql
     qq0 = real(ii-1,dp) / q0d
     if (qq0 > qq0max) exit
     ff = exp( - beta1 * (qq0 - 1) * (qq0 - 1) )
     xx = r * delql * real(ii-1,dp)
     y(ii) = vql(ii)*besselj(ll,xx)*ff
  enddo

  fsum = zero
  do ii = 1, nql-4, 4
     fsum = fsum + 7.0*y(ii) + 32.0*y(ii+1) + &
          12.0*y(ii+2) + 32.0*y(ii+3) + 7.0*y(ii+4)
  enddo
  fint = fsum * delql * two / five / nine

end function fint
!===============================================================
!
! Similar to fint, but calculates the derivative of fint,
! dfint = (d/dr) fint(r).
!
!---------------------------------------------------------------
function dfint(r,nql,delql,vql,ll,alpha,beta1,h)

  use constants
  implicit none
  !
  ! Input/Output variables:
  !
  ! radial position at which vr is to be calculated
  real(dp), intent(in) :: r
  ! number of reciprocal space points
  integer, intent(in) :: nql
  ! grid spacing
  real(dp), intent(in) :: delql
  ! input function in reciprocal space
  real(dp), intent(in) :: vql(nql)
  ! order of Bessel function
  integer, intent(in) :: ll
  ! parameters for filter function
  real(dp), intent(in) :: alpha, beta1, h
  ! derivative of filtered function
  real(dp) :: dfint
  !
  ! Work variables:
  !
  ! counters
  integer :: ii, nfilter
  real(dp) :: xx, ff, q0d, qq0, qq0max, dbessel, fsum
  ! integrand in Newton-Cottes integration
  real(dp), dimension(nql) :: y
  !
  ! External function:
  real(dp), external :: besselj

  y = zero
  q0d = alpha*pi/h/delql
  nfilter = int( q0d ) + 1
  if (beta1 > zero) then
     qq0max = sqrt(10.d0/beta1) + one
  else
     qq0max = one
  endif

  do ii = 2, nfilter
     xx = r * delql * real(ii-1,dp)
     dbessel = real(ll,dp)*besselj(ll,xx)/xx - besselj(ll+1,xx)
     y(ii) = vql(ii)*dbessel*delql * real(ii-1,dp)
  enddo

  do ii = nfilter + 1, nql
     qq0 = real(ii-1,dp) / q0d
     if (qq0 > qq0max) exit
     ff = exp( - beta1 * (qq0 - 1) * (qq0 - 1) )
     xx = r * delql * real(ii-1,dp)
     dbessel = real(ll,dp)*besselj(ll,xx)/xx - besselj(ll+1,xx)
     y(ii) = vql(ii)*dbessel*delql * real(ii-1,dp)*ff
  enddo

  fsum = zero
  do ii = 1, nql-4, 4
     fsum = fsum + 7.0*y(ii) + 32.0*y(ii+1) + &
          12.0*y(ii+2) + 32.0*y(ii+3) + 7.0*y(ii+4)
  enddo
  dfint = fsum * delql * two / five / nine

end function dfint
!===============================================================
!
! Spherical Bessel function of the first kind.
!
!---------------------------------------------------------------
function besselj(l,x)

  use constants
  implicit none
  !
  ! Input/Output variables:
  !
  ! argument
  real(dp), intent(in) :: x
  ! order
  integer, intent(in) :: l
  ! function value
  real(dp) :: besselj
  !
  ! Work variables:
  !
  integer :: j
  real(dp) :: x2,sb0,sb1,by,bym,byp,ux
  real(dp), parameter :: ten = 10.0,fourtn = 14.0,tol=0.001
  !---------------------------------------------------------------
  if(abs(x) > tol) then
     sb0 = sin(x)/x
  else
     x2 = x*x/two
     sb0 = one - (x2/three)*(one - x2/ten)
  endif
  if(l == 0) then
     besselj = sb0
  else
     if(abs(x) > tol) then
        sb1 = (sin(x)/x - cos(x)) / x
     else
        x2 = x*x/two
        sb1 = (x/three)*(one - (x2/five)*(1.0 - x2/fourtn))
     endif
     if(l == 1) then
        besselj = sb1
     elseif(x == zero) then
        besselj = zero
     else
        by = sb1
        bym = sb0
        ux = one / x
        do j=1,l-1
           byp = real(2*j+1,dp)*ux*by - bym
           bym = by
           by = byp
        enddo
        besselj = by
     endif
  endif

end function besselj
!===============================================================
