!===============================================================
!
! Copyright (C) 2005 Finite Difference Research Group
! This file is part of parsec, http://www.ices.utexas.edu/parsec/
!
! This suborutine computes the exchange correlation potential,
! Vxc, and the exchange correlation energy, Exc, both in Rydbergs,
! based on an input charge density given in e/bohr^3, for the
! case of a spin-unpolarized calculation.
!
! NOTE: Extra factor of two in formulas is due to conversion from
! Hartree to Rydberg!
!
! The following exchange-correlation functionals are implemented:
!
! Local density functionals:
! 1. Slater's x-alpha - J. C. Slater, Phys. Rev. 81, 385 (1951).
! 2. Wigner's interpolation - E. Wigner, Phys. Rev. 46, 1002
!           (1934).
! 3. Hedin-Lundqvist -  L. Hedin and B. I. Lundqvist, J. Phys. C,
!           4, 2064 (1971).
! 4. Ceperley-Alder - D. M. Ceperley and B. J. Alder, Phys. Rev.
!           Lett., 45, 566 (1980).
!
! Gradient-corrected functional:
! PBE - Perdew, Burke, Ernzerhof - J. P. Perdew, K. Burke, and 
!           M. Ernzerhof, Phys. Rev. Lett., 77, 3865 (1996).
!
! Asymptotically-corrected functionals:
! 1. Leeuwen-Baerends correction - Leeuwen and Baerends, 
!           Phys. Rev. 49, 2421 (1994).
! 2. Casida-Salahub correction - Casida and Salahub, J. Chem.
!           Phys., 113, 8918 (2000).
!
!---------------------------------------------------------------
subroutine exc_nspn(grid,parallel,nrep,icorr,dioniz,vxc,wvec,exc)

  use constants
  use grid_module
  use parallel_data_module
  implicit none
  !
  ! Input/Output variables:
  !
  ! grid related data
  type (grid_data), intent(in) :: grid
  ! parallel computation related data
  type (parallel_data), intent(in) :: parallel

  ! correlation type
  character (len=2), intent(in) :: icorr
  ! ionization energy difference
  ! (used only if the correlation is of type cs)
  real(dp), intent(in) :: dioniz
  ! order of the reduced group
  integer, intent(in) :: nrep
  ! output exchange-correlation potential
  real(dp), intent(out) :: vxc(parallel%mydim)
  ! work array
  real(dp), intent(out) :: wvec(parallel%nwedge + 1)
  ! total exchange-correlation energy
  real(dp), intent(out) :: exc
  real(dp), dimension(1):: excvec
  !
  ! Work variables:
  !
  ! multiplicity of grid points in irreducible wedge (equal to
  ! order of the reduced group)
  real(dp) :: weight

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

  ! quantitites read from pbe.f:
  ! expbe,vxuppbe,vxdnpbe - exchange enregy and (up and down)
  ! potentials
  ! ecpbe,vcuppbe,vcdnpbe - correlation enregy and (up and down)
  ! potentials
  real(dp) :: expbe,vxuppbe,vxdnpbe,vcuppbe,vcdnpbe,ecpbe
  ! equivalent quantities for pw91
  real(dp) :: expw91,vxpw91,ecpw91,vcuppw91,vcdnpw91

  ! LDA variables:
  ! various recurring parameters in the exchange-correlation
  ! formulas: 
  ! quantities common to all local density expressions:
  ! rs is the local value of the Wigner-Seitz radius
  ! rho temporarily stores the charge at a given grid point
  real(dp) :: a0, twovpia0, p75vpi, alpha, rs, rho
  ! Wigner parameters:
  real(dp), parameter :: aw = -0.875529d0 , bw = 7.8d0
  ! Hedin-Lundqvist parameters:
  real(dp), parameter :: ah = 1.d0/21.d0, ch = 0.045d0
  real(dp) :: x 
  ! Ceperley-Alder parameters
  real(dp), parameter :: g =-0.2846d0, b1 = 1.0529d0
  real(dp), parameter :: b2 = 0.3334d0, c1 = 0.0622d0
  real(dp), parameter :: c2 = 0.096d0, c3 = 0.004d0
  real(dp), parameter :: c4 = 0.0232d0,c5 = 0.0192d0
  real(dp) :: ec, sqrs, eurs
  ! Asymptotic corrections parameters and temporary variables
  real(dp), parameter :: betalb = 0.05d0
  real(dp) :: ac1, ac2 
  real(dp) :: abstemp, corrlb, asinhx

  ! actual number of grid points
  integer ndim
  ! allocation check
  integer alcstat
  ! temporary variables
  real(dp), dimension (3) :: temp1,temp2
  ! counters
  integer i 

  !---------------------------------------------------------------

  ! initialize arrays:

  ndim = parallel%mydim
  weight = real(nrep,dp)

  allocate(gradient(3,ndim),stat=alcstat)
  call alccheck('gradient',3*ndim,alcstat)
  allocate(absgrad(ndim),stat=alcstat)
  call alccheck('absgrad',ndim,alcstat)
  allocate(gradgrad(3,ndim),stat=alcstat)
  call alccheck('gradgrad',3*ndim,alcstat)
  allocate(laplac(ndim),stat=alcstat)
  call alccheck('laplac',ndim,alcstat)

  ! and some more constants:
  a0 = (four/(nine*pi))**third
  twovpia0 = two/(pi*a0)
  p75vpi = 0.75d0/pi
  !
  ! initialize the total exchange-correlation energy to zero
  !
  exc = zero

  if (icorr == 'pw' .or. icorr == 'pb') then

     ! calculate the gradient and the Laplacian of the charge density
     ! at each grid point

     call lapmvs(parallel,vxc,laplac,grid%coe2,grid%norder, &
          grid%lap_dir_num,wvec)
     laplac = -laplac
! maybe -      call daxpy instead
     call gradmvs(grid,parallel,vxc,gradient,grid%coe1,grid%norder, &
          wvec)
     do i = 1,ndim
        temp1 = gradient(:,i)
        absgrad(i) = sqrt(dot_product(temp1,temp1))
     enddo

     ! calculate other input quantities and then calculate the 
     ! exchange-correlation

     ! calculate grad(rho)_dot_grad|grad(rho)| and
     ! calculate the exchange-correlation value using pbe
     exc = zero
     call gradmvs(grid,parallel,absgrad,gradgrad,grid%coe1, &
          grid%norder,wvec)

     if (icorr == 'pb') then
        ! using PBE - a gradient-corrected functional
        do i = 1,ndim
           temp1 = gradient(:,i)
           temp2 = gradgrad(:,i)
           prod = dot_product(temp1,temp2)
           call pbe(vxc(i),absgrad(i),prod,laplac(i), &
                expbe,vxuppbe,vxdnpbe,ecpbe,vcuppbe,vcdnpbe)
           exc=exc+vxc(i)*two*(expbe+ecpbe)*weight
           vxc(i)=vxuppbe+vxdnpbe+vcuppbe+vcdnpbe
        enddo

     elseif (icorr == 'pw') then
        ! using PW91 - a gradient-corrected functional
        do i = 1,ndim
           temp1 = gradient(:,i)
           temp2 = gradgrad(:,i)
           prod = dot_product(temp1,temp2)
           call pw91(vxc(i),absgrad(i),prod,laplac(i), &
                expw91,vxpw91,ecpw91,vcuppw91,vcdnpw91)
           exc=exc+vxc(i)*two*(expw91+ecpw91)*weight
           vxc(i)=two*vxpw91+vcuppw91+vcdnpw91
        enddo
     endif

  else if (icorr == 'xa') then
     ! x-alpha correlation
   
     alpha = one
     ! Note: alpha = 2/3 is pure exchange
     write(7,10) alpha
10   format(' performed x-alpha correlation with alpha = ',f6.3)

     do i=1,ndim
        rho = vxc(i)
        vxc(i) = zero
        if (rho > zero) then
           rs = (p75vpi/rho)**third
           vxc(i) = -1.5d0*twovpia0*alpha/rs
           exc = exc + 0.75d0*rho*vxc(i)*weight
        endif
     enddo

  else if (icorr == 'wi') then
     ! Wigner correlation

     do i=1,ndim
        rho = vxc(i)
        vxc(i) = zero
        if (rho > zero) then
           rs = (p75vpi/rho)**third
           vxc(i) = -twovpia0/rs
           exc = exc + 0.75d0*rho*vxc(i)*weight + rho*aw/(rs + bw)
           vxc(i) = vxc(i) + aw*(four*rs*third+bw)/((rs+bw)*(rs+bw))
        endif
     enddo

  else if (icorr == 'hl') then
     ! Hedin-lundqvist exchange-correlation

     do i=1,ndim
        rho = vxc(i)
        vxc(i) = zero
        if (rho > zero) then
           rs = (p75vpi/rho)**third
           vxc(i) = -twovpia0/rs
           exc = exc + 0.75d0*rho*vxc(i)*weight
           x = rs*ah
           alpha = log(one+one/x)
           vxc(i) = vxc(i) - ch*alpha
           exc = exc - rho*ch*((one+x*x*x)*alpha+x*half-x*x-third)*weight
        endif
     enddo

  else if (icorr == 'ca') then
     ! Ceperly-Alder exchange correlation

     do i=1,ndim
        rho = vxc(i)
        vxc(i) = zero
        if (rho > zero) then
           rs = (p75vpi/rho)**third
           vxc(i) = -twovpia0/rs
           exc = exc + 0.75d0*rho*vxc(i)*weight
           if (rs >= one) then
              sqrs = sqrt(rs)
              ec = g/(one + b1*sqrs + b2*rs)
              vxc(i) = vxc(i) + ec*ec* &
                   (one+3.5d0*b1*sqrs*third+four*b2*rs*third)/g
           else
              alpha = log(rs)
              ec = c1*alpha - c2 + (c3*alpha - c4)*rs
              vxc(i) = vxc(i) + ec - (c1 + (c3*alpha - c5)*rs)*third
           endif
           exc = exc + rho*ec*weight
        endif
     enddo

  elseif (icorr == 'pl') then
     ! Perdew-Wang local exchange correlation
     do i=1,ndim
        rho = vxc(i)
        vxc(i) = zero
        if (rho > zero) then
           rs = (p75vpi/rho)**third
           vxc(i) = -twovpia0/rs
           exc = exc + 0.75d0*rho*vxc(i)*weight
           sqrs=sqrt(rs)
           call gcor2(0.0310907d0,0.21370d0,7.5957d0,3.5876d0,1.6382d0, &
                0.49294d0,sqrs,ec,eurs)
           vxc(i) = vxc(i) + two*( ec -rs*eurs/three )
           exc = exc + rho*ec*weight*two
        endif
     enddo

  else if (icorr == 'lb' .or. icorr == 'cs') then
     ! Asymptotically corrected Ceperly-Alder exchange correlation

     ! define parameters
     ac1=two**third
     ac2=four/three

     ! calculate 2^(1/3) * |grad(rho)|/rho^(4/3)
     call gradmvs(grid,parallel,vxc,gradient,grid%coe1,grid%norder,wvec)
     do i = 1,ndim
        temp1 = gradient(:,i)
        abstemp = sqrt(dot_product(temp1,temp1))
        if (vxc(i) > zero) then
           absgrad(i)=ac1*abstemp/(vxc(i)**ac2)
        endif
     enddo

     ! compute exchange-correlation potential and energy
     do i=1,ndim
        rho = vxc(i)
        vxc(i) = zero
        if (rho > zero) then

           ! compute Ceperely-Alder term
           rs = (p75vpi/rho)**third
           vxc(i) = -twovpia0/rs
           exc = exc + 0.75d0*rho*vxc(i)*weight
           if (rs >= one) then
              sqrs = sqrt(rs)
              ec = g/(one + b1*sqrs + b2*rs)
              vxc(i) = vxc(i) + ec*ec*(one+3.5d0*b1*sqrs*third + &
                   four*b2*rs*third)/g
           else
              alpha = log(rs)
              ec = c1*alpha - c2 + (c3*alpha - c4)*rs
              vxc(i) = vxc(i) + ec - (c1 + (c3*alpha - c5)*rs)*third
           endif
           exc = exc + rho*ec*weight

           ! Compute Leeuwen-Baerends correction to the potential.
           abstemp = absgrad(i)
           asinhx = log(abstemp+sqrt(abstemp*abstemp+one))
           corrlb = two*betalb*(rho**third)/ac1*abstemp*abstemp/ &
                (one + three*betalb*abstemp*asinhx)

           ! If correcting for the ionization energy difference, shift the
           ! correlation energy.
           if (icorr == 'cs') then
              if (corrlb > dioniz) corrlb = dioniz
           endif
           vxc(i) = vxc(i) - corrlb
        endif
     enddo
  endif
  !
  ! Scale the total energy integral by h^3 (factor necessary due to the
  ! expression of an integral as a summation over grid points) and
  ! sum it up over processors.
  !
  excvec = exc
  call psum(excvec,1,parallel%group_size,parallel%group_comm)
  exc  = excvec(1) * (grid%hcub)

  deallocate(laplac)
  deallocate(absgrad)
  deallocate(gradient)
  deallocate(gradgrad)

end subroutine exc_nspn
!===============================================================
