!===============================================================
!
! Copyright (C) 2005 Finite Difference Research Group
! This file is part of parsec, http://www.ices.utexas.edu/parsec/
!
! Computes the Ewald contribution to the total energy and forces
! by performing two convergent summations, one over lattice
! vectors and the other over reciprocal lattice vectors.
!
! Adapted from plane-wave programs written by S. Froyen and
! J. L. Martins.
!
!---------------------------------------------------------------
subroutine ewald_sum(clust,pbc,eewald,ipr,zv)

  use constants
  use cluster_module
  use pbc_module
  implicit none
  !
  ! Input/Output variables:
  !
  ! the cluster
  type (cluster), intent(inout) :: clust
  ! pbc related data
  type (pbc_data), intent(inout) :: pbc
  ! print flag
  integer, intent(in) :: ipr
  ! Ewald summation
  real(dp), intent(out) :: eewald

  real(dp), intent(in) :: zv(clust%type_num)
  !
  ! Work variables:
  !
  ! allocation check
  integer alcstat

  integer :: i, j, k, n, imx, jmx, kmx, ii, jj, ngg, nss
  real(dp) :: vcell,arg,bdot(3,3),adot(3,3)

  real(dp) :: fsub(3),rp(3),tt(3)
  real(dp) :: eps,seps,sepi,rmax,fpiv,fpiv2
  real(dp) :: esumg,esumr,esum0,ztemp,rmod,exp1,exp2
  real(dp) :: esub,sumc,sumr,gdr,cosg
  real(dp) :: sing,expgc,expgr
  real(dp), dimension(3,pbc%atom_num) ::  &
       fsumg, fcg, fsg, fsumr, tau, rc

  real(dp), allocatable :: expg(:)
  real(dp) :: zz(pbc%atom_num*pbc%type_num)
  real(dp) :: pi_dev_eps_v
  ! constants
  real(dp), parameter :: twopio = 0.5d0/pi
  real(dp), parameter :: small = 1.d-20
  real(dp), parameter :: p25 = 0.25d0
  real(dp), parameter :: tol = 200.d0, small2 = 6.5d0
  !
  ! External functions:
  !
  real(dp), external :: erfc

  !---------------------------------------------------------------

  bdot = pbc%bdot

  adot = pbc%adot

  allocate(expg(pbc%nstar),stat=alcstat)
  call alccheck('expg',pbc%nstar,alcstat)

  vcell = pbc%vcell
  eps = p25*two*pbc%ek(pbc%nstar)/tol
  seps = sqrt(eps)
  sepi = two*seps/sqrt(pi)
  rmax = sqrt(tol/eps)
  fpiv = four*pi/vcell
  fpiv2 = two*fpiv

  imx = int(rmax*sqrt(bdot(1,1))*twopio) + 1
  jmx = int(rmax*sqrt(bdot(2,2))*twopio) + 1
  kmx = int(rmax*sqrt(bdot(3,3))*twopio) + 1

  ! Store data for atoms into new arrays and fold into first cell:
  ! prevents errors in sums from too small a real space sum.
  k = 0
  write(7,*) 'Reciprocal lattice atom data'
  do i = 1, clust%type_num
     do j = 1, clust%natmi(i)
        k = k + 1
        zz(k) = zv(i)
        rc(1:3,k) =  pbc%rat(1:3,j,i)*twopio
        write(7,'(3(3x,f11.6))') pbc%rat(:,j,i)
        do n = 1, 3
           if (rc(n,k) >= one) &
                rc(n,k) = modulo(rc(n,k) + eight,one)
           if (rc(n,k) < zero) &
                rc(n,k) = modulo(rc(n,k),one)
        enddo
     enddo
  enddo

  esumg = zero
  esumr = zero

  fsumg(:,:) = zero
  fsumr(:,:) = zero

  ! start sum in g space
  esum0 = zero
  ztemp = -p25/eps
  do i = 2, pbc%nstar
     expg(i) = exp(ztemp*two*pbc%ek(i))/(two*pbc%ek(i))
     if (expg(i) < small) exit
  enddo
  nss = i - 1
  ngg = pbc%ng 

  ngg = ngg - sum(pbc%mstar(nss+1:pbc%nstar))

  ztemp = half/eps
  do i = 2, ngg
     sumc=zero
     sumr=zero

     ! start loop over atoms in cell
     do j = 1, clust%atom_num
        gdr = twopi * dot_product(real(pbc%kgv(:,i),dp),rc(:,j))
        cosg = zz(j)*cos(gdr)
        sing = zz(j)*sin(gdr)
        sumc = sumc + cosg
        sumr = sumr + sing
        fcg(:,j) = -real(pbc%kgv(:,i),dp)*cosg
        fsg(:,j) =  real(pbc%kgv(:,i),dp)*sing
     enddo
     expgc = expg(pbc%inds(i))*sumc
     expgr = expg(pbc%inds(i))*sumr
     esumg = esumg + expgc*sumc+expgr*sumr
     do j = 1, clust%atom_num
        fsumg(:,j) = fsumg(:,j) + expgc*fsg(:,j) + expgr*fcg(:,j)
     enddo
  enddo

  esumg = fpiv*esumg
  fsumg = fpiv2*fsumg

  pi_dev_eps_v =  pi/(eps*vcell)

  ! end g sum
  !
  ! start sum in r space
  !
  ! AMIR: the different summations are repeated several times
  ! to avoid the situation of rmod=0. perhaps there is
  ! a problem here when working with a general metric. 
  !
  esum0 = zero
  do i = -imx, imx
     rp(1) = i
     do j = -jmx, jmx
        rp(2) = j
        do k = -kmx, kmx
           rp(3) = k
           ! computation of the metric.
           call matvec3('N',adot,rp,tt)
           rmod = sqrt(dot_product(rp,tt))
           ! skip point rp = (0,0,0)
           if (rp(1) == 0 .and. rp(2) == 0 .and. rp(3) == 0 ) cycle
           arg = seps*rmod
           if (arg < small2) then
              exp1 = erfc(arg) / rmod
              esum0 = esum0 + exp1
           endif
        enddo
     enddo
  enddo

  esum0 = esum0 - pi/(eps*vcell) - sepi     

  ! start loop over atoms in cell
  ! term with a=b (atom 1)
  ztemp = zz(1)*zz(1)
  esumr = esumr + ztemp*esum0

  do i = 2, clust%atom_num

     ! term with a=b (atom >1)
     ztemp = zz(i)*zz(i)
     esumr = esumr + ztemp*esum0

     ! terms with a#b
     do j=1,i-1

        ! loop over lattice points
        esub = zero
        fsub(:) = zero
        do ii = -imx, imx
           rp(1) = real(ii,dp) + rc(1,i) - rc(1,j)
           do jj = -jmx, jmx
              rp(2) = real(jj,dp) + rc(2,i) - rc(2,j)
              do k = -kmx, kmx
                 rp(3) = real(k,dp) + rc(3,i) - rc(3,j)
                 call matvec3('N',adot,rp,tt)
                 rmod = sqrt(dot_product(rp,tt))
                 arg = seps*rmod
                 if (arg < small2) then
                    exp1 = erfc(arg) / rmod
                    exp2 = (exp1+sepi*exp(-arg*arg))/(rmod*rmod)
                    esub = esub + exp1
                    tt = rp * exp2
                    fsub = fsub + tt
                 endif
              enddo
           enddo
        enddo
        esub = esub - pi_dev_eps_v
        ztemp = two*zz(i)*zz(j)
        esumr = esumr + ztemp*esub
        tt = ztemp*fsub
        fsumr(:,i) = fsumr(:,i) + tt
        fsumr(:,j) = fsumr(:,j) - tt
     enddo
  enddo

  ! end r sum
  eewald = esumg + esumr

  ! force
  ! note - returned force is in units of real space lattice vectors
  ! printed force is in cartesian coordinates
  do ii = 1, 3
     clust%force(ii,:) = fsumr(ii,:)
     do jj = 1, 3
        clust%force(ii,:) = clust%force(ii,:) + &
             pbc%bdot(ii,jj)*fsumg(jj,:)*twopio
     enddo
  enddo
  if (allocated(expg)) deallocate(expg)

  ! if print flag on, print the ion-ion forces
  do i = 1, clust%atom_num
     call matvec3('N',pbc%latt_vec,clust%force(1,i),tau(1,i))
     clust%force(:,i)=tau(:,i)
  enddo
  if(ipr >= 1) then
     write(7,*)
     write(7,*) 'Forces from ion-ion interaction:'
     write(7,*) '================================'
     write(7,22)
     write(7,*)
     do i = 1, clust%atom_num
        write(7,20) i,clust%xatm(i),clust%yatm(i),clust%zatm(i),tau(:,i)
     enddo
     write(7,*)
  endif

20 format(i4,2x,3(f11.6,1x),2x,3(f11.6,1x))
22 format('-atom-  ----x----   ----y----   ----z-----' &
        ,'   ----Fx----  ----Fy----   ----Fz---',/ &
        '                       [bohr]           ' &
        ,'                  [Ry/bohr]            ')

end subroutine ewald_sum
!===============================================================
