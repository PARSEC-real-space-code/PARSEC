!===============================================================
!
! Copyright (C) 2005 Finite Difference Research Group
! This file is part of parsec, http://www.ices.utexas.edu/parsec/
!
! This subroutine computes the exchange correlation potential,
! Vxc, and the exchange correlation energy, Exc, both in Rydbergs,
! based on an input charge density given in e/(bohr^3).
!
! NOTE: Extra factor of two in formulas is due to conversion from
! Hartree to Rydberg!
!
! The following exchange-correlation functionals are implemented:
!
! Local density functionals:
! Ceperley-Alder - D. M. Ceperley and B. J. Alder, Phys. Rev.
! Lett., 45, 566 (1980).
!
! PL = PW92 = PWLDA - J. P. Perdew and Y. Wang, Phys. Rev. B 45, 13244 (1992). 
!
! Gradient-corrected functionals:
! PBE - Perdew, Burke, Ernzerhof - J. P. Perdew, K. Burke, 
! and M. Ernzerhof, Phys. Rev. Lett., 77, 3865 (1996).
!
! PW(91) = PWGGA : J. Perdew et. al., PRB 46, 6671 (1992)
!
!---------------------------------------------------------------
subroutine exc_spn(grid,parallel,nrep,icorr,dioniz,vxc,wvec,exc)

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
  ! Input: charge density. Output: exchange-correlation potential
  real(dp), intent(inout) :: vxc(parallel%mydim,2)
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

  ! should be changed to parameters
  integer, parameter :: UP = 1,DOWN = 2

  ! GGA variables:
  ! quantities fed to pbe.f:
  ! gradient - grad(n), where n is the charge density
  ! absgrad  - |grad(n)|
  ! grdgrd - components of grad(|grad(n)|
  ! laplac - laplacian(n)
  ! prod  -  grad(n) dot grad(|grad(n)|
  ! abs_total_grad - |grad(total rho)| 

  real(dp), dimension (:,:), allocatable :: &
       gradUP, gradDN, gradgradUP, gradgradDN, total_gradgrad
  real(dp),dimension (:,:), allocatable :: absgrad,laplac

  real(dp), allocatable :: abs_total_grad(:)
  real(dp),dimension (UP:DOWN):: prod
  real(dp) :: total_prod

  integer spin

  integer :: &
       lcor, &              ! flag to do correlation(=0=>don't)
       lpot                 ! flag to do potential(=0=>don't)


  real(dp) :: &
       exlsd, &             ! LSD exchange energy density, so that
                            !     ExLSD=int d^3r rho(r) exlsd(r)
       vxuplsd, &           ! up LSD exchange potential
       vxdnlsd, &           ! down LSD exchange potential
       eclsd, &             ! LSD correlation energy density
       vcuplsd, &           ! up LSD correlation potential
       vcdnlsd              ! down LSD correlation potential

  real(dp) :: &
       expw91,vxuppw91, &   !
       vxdnpw91,ecpw91, &   ! corresponding PW91 quantities
       vcuppw91,vcdnpw91    !

  real(dp) :: &
       expbe,vxuppbe, &     !
       vxdnpbe,vcuppbe, &   ! corresponding PBE quantities
       vcdnpbe,ecpbe        !

  ! LDA variables:
  ! rs is the local value of the Wigner-Seitz radius
  ! rho,vx,vc temporarily store the charge,exchange and correlation 
  ! potentials at a given grid point
  real(dp) :: rs, rh(2),vx(2),vc(2)

  ! High density (rs<1) contants, LDA
  real(dp), parameter :: ca_c1 = 0.0622d0
  real(dp), parameter :: ca_c2 = 0.0960d0
  real(dp), parameter :: ca_c3 = 0.0040d0
  real(dp), parameter :: ca_c4 = 0.0232d0
  real(dp), parameter :: ca_c5 = 0.0192d0

  real(dp), parameter :: ca_as =  0.0311d0
  real(dp), parameter :: ca_bs = -0.0538d0
  real(dp), parameter :: ca_cs =  0.0014d0
  real(dp), parameter :: ca_ds = -0.0096d0
  ! Low density (rs>1) constants, LDA
  real(dp), parameter :: ca_g = -0.2846d0
  real(dp), parameter :: ca_b1 = 1.0529d0
  real(dp), parameter :: ca_b2 = 0.3334d0
  real(dp), parameter :: ca_gs = -0.1686d0
  real(dp), parameter :: ca_b1s = 1.3981d0
  real(dp), parameter :: ca_b2s = 0.2611d0
  ! Other parameters
  real(dp), parameter :: ac2 = four/three
  real(dp), parameter :: p75vpi = 0.75d0/pi
  real(dp), parameter :: ax = -0.738558766382022405884230032680836d0
  real(dp), parameter :: gamma = 0.5198421d0
  real(dp), parameter :: fzz = 1.709921d0
  ! other temporary holders
  real(dp) :: ex(2),ec,eu,f,zet,z4,ecrs,eurs,eczet,comm,fz,d
  real(dp) :: lrs, srs, es, vs
  real(dp) :: ep,eprs,alfrsm,alfm
  real(dp) :: temp_val

  ! Asymptotic corrections parameters and temporary variables
  real(dp), parameter :: betalb = 0.05d0
  real(dp) :: abstemp1, abstemp2, corrlb1, corrlb2, asinhx

  ! counters
  integer i 
  ! allocation check
  integer alcstat

  real(dp), dimension (3) :: temp1,temp2,temp3,temp4,temp5,temp6

  ! local number of grid points
  integer ndim

  !---------------------------------------------------------------

  ndim = parallel%mydim
  weight = real(nrep,dp)
  !
  ! initialize the total exchange-correlation energy to zero
  !
  exc = zero

  if (icorr == 'pb' .or. icorr == 'pw') then

     ! using PBE - a gradient-corrected functional

     ! quantities fed to pbe.f:
     ! gradient - grad(n), where n is the charge density

     ! absgrad  - |grad(n)|
     ! gradgrad - components of grad(|grad(n)|
     ! laplac - laplacian(n)
     ! prod  -  grad(n) dot grad(|grad(n)|

     ! abs_total_grad - |grad(total rho)|

     ! initialize working arrays:

     allocate(abs_total_grad(ndim),stat=alcstat)
     call alccheck('abs_total_grad',ndim,alcstat)
     allocate(absgrad(ndim,UP:DOWN),stat=alcstat)
     call alccheck('absgrad',ndim*(DOWN-UP+1),alcstat)
     allocate(laplac(ndim,UP:DOWN),stat=alcstat)
     call alccheck('laplac',ndim*(DOWN-UP+1),alcstat)
     allocate(gradUP(3,ndim),stat=alcstat)
     call alccheck('gradUP',3*ndim,alcstat)
     allocate(gradDN(3,ndim),stat=alcstat)
     call alccheck('gradDN',3*ndim,alcstat)
     allocate(gradgradUP(3,ndim),stat=alcstat)
     call alccheck('gradgradUP',3*ndim,alcstat)
     allocate(gradgradDN(3,ndim),stat=alcstat)
     call alccheck('grdgrdDN',3*ndim,alcstat)
     allocate(total_gradgrad(3,ndim),stat=alcstat)
     call alccheck('total_gradgrad',3*ndim,alcstat)

     ! prepare vectors for invoking pbe
         
     ! calculate the gradient of the charge density (rho==exc) at each
     ! grid point: gradUP, gradDN
     call gradmvs(grid,parallel,vxc(1,UP),gradUP,grid%coe1 &
          ,grid%norder,wvec)
     call gradmvs(grid,parallel,vxc(1,DOWN),gradDN,grid%coe1 &
          ,grid%norder,wvec)

     ! calculate abs_total_grad
     ! (|grad(total rho)| = |grad(up)+grad(down)|)
     do i = 1,ndim
        temp1 = gradUP(:,i) + gradDN(:,i)
        abs_total_grad(i) = sqrt(dot_product(temp1,temp1))
     enddo

     ! calculate total_gradgrad(i) = grad|grad(rho(i))| 
     call gradmvs(grid,parallel,abs_total_grad,total_gradgrad, &
          grid%coe1,grid%norder,wvec)

     ! calculate absgrad(i,spin) -  the absolute value of the gradient
     ! in UP/DN: |grad(rho(i,spin)|
     do i = 1,ndim
        temp1 = gradUP(:,i)
        absgrad(i,UP)= sqrt(dot_product(temp1,temp1))
        temp2 = gradDN(:,i)
        absgrad(i,DOWN)= sqrt(dot_product(temp2,temp2))
     enddo

     ! calculate gradgrad_spin = grad|grad(rho(i,spin))| 
     call gradmvs(grid,parallel,absgrad(1,UP),gradgradUP,grid%coe1 &
          ,grid%norder,wvec)
     call gradmvs(grid,parallel,absgrad(1,DOWN),gradgradDN,grid%coe1 &
          ,grid%norder,wvec)

     ! Calculate laplac(i,spin) -  the Laplac of rho(spin) at each
     ! grid point.
     do spin = UP,DOWN
        call lapmvs(parallel,vxc(1,spin),laplac(1,spin),grid%coe2 &
             ,grid%norder,grid%lap_dir_num,wvec)
     enddo
     laplac = -laplac

     ! do correlation
     lcor = 1
     ! do potential
     lpot = 1

     do i = 1,ndim

        ! prod(spin) = grad(rho(spin))_dot_grad|grad(rho(spin))| 
        ! grad(rho(spin)) = gradUP , gradDN
        ! grad|grad(rho(spin))| = gradgradUP , gradgradDN
        ! UP
        temp1 = gradUP(:,i)
        temp2 = gradgradUP(:,i)
        prod(UP) = dot_product (temp1,temp2)
        ! DOWN
        temp3 = gradDN(:,i)
        temp4 = gradgradDN(:,i)
        prod(DOWN) = dot_product (temp3,temp4)

        ! Total Prod
        ! total_prod = delgr=(grad rho).(grad |grad rho|) 
        ! (grad rho) = grad(rho up) + grad (rho dn)
        ! (grad |grad rho|) = total_gradgrad
        temp5 = gradUP(:,i) + gradDN(:,i)
        temp6 = total_gradgrad(:,i)
        total_prod = dot_product(temp5,temp6)
        
        call spn_pbe(vxc(i,UP),absgrad(i,UP),prod(UP),laplac(i,UP), &
             vxc(i,DOWN),absgrad(i,DOWN),prod(DOWN),laplac(i,DOWN), &
             abs_total_grad(i),total_prod,lcor,lpot, &
             exlsd,vxuplsd,vxdnlsd,eclsd,vcuplsd,vcdnlsd, &
             expw91,vxuppw91,vxdnpw91,ecpw91,vcuppw91,vcdnpw91, &
             expbe,vxuppbe,vxdnpbe,ecpbe,vcuppbe,vcdnpbe)

        if (icorr == 'pb') then
           ! using PBE - a gradient-corrected functional
           exc=exc+( (vxc(i,UP)+vxc(i,DOWN))*two*(expbe+ecpbe) )*weight
           vxc(i,UP)=two*(vxuppbe+vcuppbe)
           vxc(i,DOWN)=two*(vxdnpbe+vcdnpbe)
           ! using non-local functional of rho, PW91
        else
           ! using GGA functional of rho, PW91
           exc=exc+( (vxc(i,UP)+vxc(i,DOWN))*two*(expw91+ecpw91) )*weight
           vxc(i,UP)=two*(vxuppw91+vcuppw91)
           vxc(i,DOWN)=two*(vxdnpw91+vcdnpw91)
        endif
     enddo

     ! free up memory
     deallocate(abs_total_grad)
     deallocate(absgrad)
     deallocate(laplac)
     deallocate(gradUP)
     deallocate(gradDN)
     deallocate(gradgradUP)
     deallocate(gradgradDN)
     deallocate(total_gradgrad)
  else if (icorr == 'pl') then
     ! Ceperly-Alder exchange correlation
     temp_val = ax*two**third
     do i=1,ndim
        rh(1)=vxc(i,1)
        rh(2)=vxc(i,2)
        if (minval(rh) > zero) then
           ! exchange
           d=two*rh(1)
           ex(1)= ax*d**third
           vx(1)=ex(1)*ac2
           d=two*rh(2)
           ex(2)= ax*d**third
           vx(2)=ex(2)*ac2
           ! local correlation
           d=sum(rh)
           rs=(p75vpi/d)**third
           call spn_gcor(0.0310907d0,0.21370d0,7.5957d0,3.5876d0, &
                1.6382d0,0.49294d0,1.00d0,rs,eu,eurs)
           call spn_gcor(0.01554535d0,0.20548d0,14.1189d0,6.1977d0, &
                3.3662d0,0.62517d0,1.00D0,rs,ep,eprs)
           call spn_gcor(0.0168869d0,0.11125d0,10.357d0,3.6231d0, &
                0.88026d0,0.49671d0,1.00d0,rs,alfm,alfrsm)
           zet=(rh(1)-rh(2))/d
           f = ((one+zet)**ac2+(one-zet)**ac2-two)/gamma
           z4 = zet**4
           ec = eu*(one-f*z4)+ep*f*z4-alfm*f*(one-z4)/fzz
           ecrs = eurs*(one-f*z4)+eprs*f*z4-alfrsm*f*(one-z4)/fzz
           fz = ac2*((one+zet)**third-(one-zet)**third)/gamma
           eczet = four*(zet**3)*f*(ep-eu+alfm/fzz)+fz*(z4*ep-z4*eu &
                -(one-z4)*alfm/fzz)
           comm = ec -rs*ecrs/three-zet*eczet
           vc(1) = comm + eczet
           vc(2) = comm - eczet

           vxc(i,1)=two*(vx(1)+vc(1))
           vxc(i,2)=two*(vx(2)+vc(2))
           ec=two*(rh(1)*ex(1)+rh(2)*ex(2)+ec*d)
           exc = exc + ec*weight
        endif
     enddo
  else if (icorr == 'ca') then
     ! Ceperly-Alder exchange correlation
     temp_val = ax*two**third
     do i=1,ndim
        rh(1)=vxc(i,1)
        rh(2)=vxc(i,2)
        if (minval(rh) > zero) then
           ! exchange
           d=two*rh(1)
           ex(1)= ax*d**third
           vx(1)=ex(1)*ac2
           d=two*rh(2)
           ex(2)= ax*d**third
           vx(2)=ex(2)*ac2
           ! local correlation
           d=sum(rh)
           rs=(p75vpi/d)**third

           zet=(rh(1)-rh(2))/d
           f = ((one+zet)**ac2+(one-zet)**ac2-two)/gamma
           fz = ac2*((one+zet)**third-(one-zet)**third)/gamma
           if (rs < one) then
              ! correlation: high density
              lrs = log(rs)
              eu = ca_c1 * lrs - ca_c2 + (ca_c3 * lrs - ca_c4) * rs  
              vc(1) = eu - (ca_c1 + (ca_c3 * lrs + ca_c3 - ca_c4) * rs) /3.d0
              es = ca_as * lrs + ca_bs + (ca_cs * lrs + ca_ds) * rs
              vs = es - (ca_as + (ca_cs * lrs + ca_cs + ca_ds) * rs) /3.d0
           else
              ! correlation: low density
              srs = sqrt(rs)
              eu = ca_g / (1.d0 + ca_b1 * srs + ca_b2 * rs)
              vc(1) = eu * eu * (1.d0 + 7.0d0 * ca_b1 * srs / 6.d0 + &
                   4.d0/3.d0 * ca_b2 * rs) / ca_g
              es = ca_gs / (1.d0 + ca_b1s * srs + ca_b2s * rs)
              vs = es * es * (1.d0 + 7.0d0 * ca_b1s * srs / 6.d0 + &
                   4.d0/3.d0 * ca_b2s * rs) / ca_gs
           endif
           vc(2) = vc(1) + f * (vs - vc(1)) - (es - eu) * &
                (1.d0 + zet) * fz
           vc(1) = vc(1) + f * (vs - vc(1)) + (es - eu) * &
                (1.d0 - zet) * fz

           ec = eu + f * (es - eu)

           vxc(i,1)=two*vx(1)+vc(1)
           vxc(i,2)=two*vx(2)+vc(2)
           ec=two*(rh(1)*ex(1)+rh(2)*ex(2))+ec*d
           exc = exc + ec*weight
        endif
     enddo

  else if (icorr == 'lb' .or. icorr == 'cs') then
     ! Asymptotically corrected Ceperly-Alder exchange correlation

     ! initialize working arrays:
     allocate(absgrad(ndim,UP:DOWN),stat=alcstat)
     call alccheck('absgrad',ndim*(DOWN-UP+1),alcstat)
     allocate(gradUP(3,ndim),stat=alcstat)
     call alccheck('gradUP',3*ndim,alcstat)
     allocate(gradDN(3,ndim),stat=alcstat)
     call alccheck('gradDN',3*ndim,alcstat)

     ! calculate |grad(rho)|/rho^(4/3) -- spinup & spindown
     call gradmvs(grid,parallel,vxc(1,UP),gradUP,grid%coe1 &
          ,grid%norder,wvec)
     call gradmvs(grid,parallel,vxc(1,DOWN),gradDN,grid%coe1 &
          ,grid%norder,wvec)
     do i = 1,ndim
        temp1 = gradUP(:,i)
        abstemp1 = sqrt(dot_product(temp1,temp1))
        temp2 = gradDN(:,i)
        abstemp2 = sqrt(dot_product(temp2,temp2))
        if (vxc(i,UP) > zero) then
           absgrad(i,UP)= abstemp1/(vxc(i,UP)**ac2)
        endif
        if (vxc(i,DOWN) > zero) then
           absgrad(i,DOWN)= abstemp2/(vxc(i,DOWN)**ac2)
        endif
     enddo

     ! compute exchange-correlation potential and energy
     temp_val = ax*two**third
     do i=1,ndim
        rh(1)=vxc(i,1)
        rh(2)=vxc(i,2)
        if (minval(rh) > zero) then

           ! compute exchange
           d=two*rh(1)
           ex(1)= ax*d**third
           vx(1)=ex(1)*ac2
           d=two*rh(2)
           ex(2)= ax*d**third
           vx(2)=ex(2)*ac2

           ! compute Ceperely-Alder term
           d=sum(rh)
           rs=(p75vpi/d)**third

           zet=(rh(1)-rh(2))/d
           f = ((one+zet)**ac2+(one-zet)**ac2-two)/gamma
           fz = ac2*((one+zet)**third-(one-zet)**third)/gamma
           if (rs < one) then
              ! correlation: high density
              lrs = log(rs)
              eu = ca_c1 * lrs - ca_c2 + (ca_c3 * lrs - ca_c4) * rs
              vc(1) = eu - (ca_c1 + (ca_c3 * lrs + ca_c3 - ca_c4) * rs) /3.d0
              es = ca_as * lrs + ca_bs + (ca_cs * lrs + ca_ds) * rs
              vs = es - (ca_as + (ca_cs * lrs + ca_cs + ca_ds) * rs) /3.d0
           else
              ! correlation: low density
              srs = sqrt(rs)
              eu = ca_g / (1.d0 + ca_b1 * srs + ca_b2 * rs)
              vc(1) = eu * eu * (1.d0 + 7.0d0 * ca_b1 * srs / 6.d0 + &
                   4.d0/3.d0 * ca_b2 * rs) / ca_g
              es = ca_gs / (1.d0 + ca_b1s * srs + ca_b2s * rs)
              vs = es * es * (1.d0 + 7.0d0 * ca_b1s * srs / 6.d0 + &
                   4.d0/3.d0 * ca_b2s * rs) / ca_gs
           endif
           vc(2) = vc(1) + f * (vs - vc(1)) - (es - eu) * &
                (1.d0 + zet) * fz
           vc(1) = vc(1) + f * (vs - vc(1)) + (es - eu) * &
                (1.d0 - zet) * fz

           ec = eu + f * (es - eu)

           vxc(i,1)=two*vx(1)+vc(1)
           vxc(i,2)=two*vx(2)+vc(2)
           ec=two*(rh(1)*ex(1)+rh(2)*ex(2))+ec*d
           exc = exc + ec*weight

           ! Compute Leeuwen-Baerends correction to the potential.
           abstemp1 = absgrad(i,UP)
           abstemp2 = absgrad(i,DOWN)
           asinhx = log(abstemp1+sqrt(abstemp1*abstemp1+one))
           corrlb1 = two*betalb*(rh(1)**third)*abstemp1*abstemp1/ &
                (one + three*betalb*abstemp1*asinhx)
           asinhx = log(abstemp2+sqrt(abstemp2*abstemp2+one))
           corrlb2 = two*betalb*(rh(2)**third)*abstemp2*abstemp2/ &
                (one + three*betalb*abstemp2*asinhx)

           ! If correcting for the ionization energy difference, shift the
           ! correlation energy.
           if (icorr == 'cs') then
              if (corrlb1 > dioniz) corrlb1 = dioniz
              if (corrlb2 > dioniz) corrlb2 = dioniz
           endif
           vxc(i,1) = vxc(i,1) - corrlb1
           vxc(i,2) = vxc(i,2) - corrlb2
        endif
    enddo

    ! free up memory
     deallocate(absgrad)
     deallocate(gradUP)
     deallocate(gradDN)

  endif
  !
  ! Scale the total energy integral by h^3 (factor necessary due to
  ! the expression of an integral as a summation over grid points) and
  ! sum it up over processors.
  !
  excvec = exc
  call psum(excvec,1,parallel%group_size,parallel%group_comm)
  exc  = excvec(1)  * (grid%hcub)

end subroutine exc_spn
!===============================================================
