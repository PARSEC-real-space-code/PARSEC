!===============================================================
!
! Copyright (C) 2005 Finite Difference Research Group
! This file is part of parsec, http://www.ices.utexas.edu/parsec/
!
! Computes the Ewald contribution to the total energy and forces
! for slab geometry by performing two convergent summations, 
! one over lattice vectors and the other over reciprocal lattice vectors.
!
! Full description and equations can be found in the book 
! "First-Principles Calculations In Real-Space Formalism"
! by Tomoya Ono and Kikuji Hirose, Appendix A.
!
! Author: Ayelet Mor (2007) 
!
!---------------------------------------------------------------
subroutine ewald_slab(clust,pbc,eewald,ipr,zv)

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
  ! Ewald energy summation
  real(dp), intent(out) :: eewald
  ! print flag
  integer, intent(in) :: ipr
  ! Vector of Electric charge for each chemical element
  real(dp), intent(in) :: zv(clust%type_num)


  !
  ! Work variables:
  !
  
  ! number of G vectors in the reciprocal space
  integer :: ng
  ! area of one unit cell
  real(dp) :: latt_area
  ! pia = pi/latt_area
  real(dp) :: pia
  ! sqrtpi = sqrt(pi)
  real(dp) :: sqrtpi
  ! inveta = 1 / eta, 
  ! eta is the partitioning constant between the real and reciprocal space
  real(dp) :: inveta
  ! eta_sqrtpi = eta / sqrt(pi)
  real(dp) :: eta_sqrtpi
  ! rmax - the maximum for real space P vectors
  real(dp) :: rmax
  ! maximum P vector values in each direction (u & v)
  integer :: imx, jmx
  ! Charged-sheet energy
  real(dp) :: enucsh

  ! Vector quantities:
  real(dp) :: &
       ratom(3, pbc%atom_num*pbc%type_num), &   ! ratom(ia) = coordinates for atom no. ia
       zz(pbc%atom_num*pbc%type_num)            ! zz(ia) = the electric charge on atom ia

  ! Reciprocal space variables:
  real(dp) :: &
       gvec(3), &        ! the G vector
       gvec_size, &      ! gvec_size = |G|
       zdiff, &          ! Z(ja) - Z(ia)
       foutput1, &       ! foutput1 = fplus(G(ig),Z(ja)-Z(ia))
       foutput2, &       ! foutput2 = fminus(G(ig),Z(ja)-Z(ia))
       gmultri, &        ! gmultri = G*R(ia)
       cosi,    &        ! cos(gmultri), cos(gmultrj)
       sini              ! sin(gmultri), sin(gmultrj)
       
  ! General summation variables
  real(dp) :: &
       rrz , &           ! rrz = abs(ratom(3,ja) - ratom(3,ia))
       rrz_sign, &       ! rrz_sign = (Z(ja) - Z(ia))/|Z(ja) - Z(ia)|
       zdiff_eta, &      ! zdiff_eta = eta * rrz
       erf_output        ! erf_output = erf(zdiff_eta)

  ! Real space variables:
  real(dp) :: &
       adot(3,3), &      ! adot = pbc%adot
       bdot(3,3), &      ! bdot = pbc%bdot
       rp(3), &          ! real space coefficients of vector P
       pvec(3), &        ! the vector (P + R(ja) - R(ia))
       rmod, &           ! rmod = |P + R(ja) - R(ia)|
       arg, &            ! arg = eta * rmod
       tmp, tmp2  

  ! Counters:
  integer :: i, j, ia, ja, ig, natom, ii, jj, n

  ! Summing variables
  real(dp) :: &
       esumr, &          ! energy sum over real-space
       esumg, &          ! energy sum over the reciprocal space
       esum0, &          ! energy sum over ia == ja
       esumm, &          ! general energy sum over ia and ja
       esub              ! energy sub-summation in real-space

  real(dp), dimension(3,pbc%atom_num) :: &
       fsumr, &          ! force sum over real-space
       fsumg, &          ! force sum over reciprocal space
       fsumm, &          ! general force sum over ia and ja
       tau               ! total force in latt vec coordinates
  
  real(dp) :: fsub(3)    ! force sub-summationn in real-space

  real(dp) :: &
       sum_fact, &       ! multipliance factor
       sum_fact2, &      ! multipliance factor
       sum_fact3         ! multipliance factor

  ! Convergenace check
  real(dp) :: prevesumg, prevfsumg(3)

  ! Constants
  real(dp), parameter :: eta = 0.5d0     ! Determines the border between the two summations
  real(dp), parameter :: twopio = 0.5d0/pi
  real(dp), parameter :: small = 1.d-9
  real(dp), parameter :: small2 = 6.5d0
  real(dp), parameter :: tol = 200.d0

  ! External functions:
  real(dp), external :: erfc             ! Complementary error function

  ! Temporal variables
  real(dp) :: tmpvec(3), rtmp(3)
    
  !---------------------------------------------------------------

  !
  ! Initiate vsriables
  !
  latt_area = abs( pbc%latt_vec(1,1) * pbc%latt_vec(2,2) - &
       pbc%latt_vec(2,1) * pbc%latt_vec(1,2) )
  pia = pi / latt_area
  sqrtpi = sqrt(pi)
  inveta = one / eta
  eta_sqrtpi = eta / sqrtpi

  ng = pbc%ng
  adot = pbc%adot
  bdot = pbc%bdot

  gvec(:) = zero
  rp(:) = zero
  pvec(:) = zero
  ratom(:,:) = zero

  rmax = sqrt(tol/eta)
  imx = int(rmax*sqrt(bdot(1,1))*twopio) + 1
  jmx = int(rmax*sqrt(bdot(2,2))*twopio) + 1

  write(9,*) 'ewald real space limits ', imx, jmx

  ! Define coordinates of all atoms with respect to the corner of the
  ! periodic cell.


  natom = 0
  do i = 1, clust%type_num
     do j = 1, clust%natmi(i)
        natom = natom + 1

        tmpvec(1)=clust%xatm(natom)
        tmpvec(2)=clust%yatm(natom)
        tmpvec(3)=clust%zatm(natom)

        ! Multiplying by the bvec matrix to find
        ! the coordinates relative to lattice vec.
        ! Since (u,v,w) = (x,y,z) * A^-1 = (x,y,z) * B' /twopi
        call matvec3('T',pbc%bvec,tmpvec,pbc%rat(1:3,j,i))
     enddo
  enddo




  ! Build the Z table (electric charge) for all atoms
  ! and fill in their coordinates
  natom = 0
  do i = 1, clust%type_num
     do j = 1, clust%natmi(i)
        natom = natom + 1
        zz(natom) = zv(i)
        ratom(1:3,natom) = pbc%rat(1:3,j,i)*twopio
        write(9,*) 'ratom(',natom,')=',ratom(1:3,natom)
        do n = 1, 3
           if (ratom(n,natom) >= one) &
                ratom(n,natom) = modulo(ratom(n,natom) + eight,one)
           if (ratom(n,natom) < zero) &
                ratom(n,natom) = modulo(ratom(n,natom),one)
        enddo
     enddo
  enddo

  !
  ! Sum over the G-space (reciprocal space)
  ! stop when convergence is reached
  esumg = zero
  fsumg(:,:) = zero
  prevesumg = zero
  prevfsumg(:) = zero

  write(9,*) 'ng=',ng

  ng=min(ng,1000)

  do ig = 2, ng
     gvec = pbc%kgv(:,ig)
     gvec(3) = 0
     tmpvec=gvec(1)*pbc%bvec(:,1)+gvec(2)*pbc%bvec(:,2)
     gvec_size = sqrt(dot_product(tmpvec,tmpvec))
     sum_fact = one / gvec_size
     
     do ja = 1, natom

        do ia = 1, natom
           rtmp=ratom(:,ja)-ratom(:,ia)
           gmultri = twopi*dot_product(gvec,rtmp)
           cosi = cos(gmultri)
           sini = sin(gmultri)
           
           zdiff = clust%zatm(ja) - clust%zatm(ia)
           call fplus(gvec_size,zdiff,eta,foutput1)
           call fminus(gvec_size,zdiff,eta,foutput2)


           sum_fact2 = zz(ja)*zz(ia)*sum_fact*foutput1
           sum_fact3 = zz(ja)*zz(ia)*foutput2

           esumg = esumg + sum_fact2 * cosi
           fsumg(1,ia) = fsumg(1,ia) + tmpvec(1) * sum_fact2 * sini
           fsumg(2,ia) = fsumg(2,ia) + tmpvec(2) * sum_fact2 * sini
           fsumg(3,ia) = fsumg(3,ia) + sum_fact3 * cosi
        enddo
     enddo
     
     ! Check if convergence is reached.
     ! This part should be checked again for better convergence parameters. AYELET


  enddo

  esumg = esumg * pia
  fsumg = -fsumg * pia


  !
  ! Sum the general summation
  !
  esumm = zero
  fsumm(:,:) = zero

  do ja = 1, natom
     sum_fact = zz(ja)
     do ia = 1, natom
        rrz = abs(clust%zatm(ja) - clust%zatm(ia))
        if(rrz>0) then
           rrz_sign = (clust%zatm(ja) - clust%zatm(ia))/rrz
        else
           rrz_sign = 1
        endif
        zdiff_eta = rrz*eta
        erf_output=1-erfc(zdiff_eta)
!        call erf(zdiff_eta,erf_output)

        esumm = esumm + sum_fact*zz(ia) * (inveta*exp(-(zdiff_eta**2)) + sqrtpi*rrz*erf_output)
        fsumm(3,ia) = fsumm(3,ia) + sum_fact*zz(ia) * erf_output * rrz_sign
     enddo
  enddo

  esumm = esumm * (-two*sqrtpi/latt_area)
  fsumm = fsumm * (-two*pia)

  !
  ! Sum over the real space  
  ! the summation is devided into two to avoid the rmod==0 case 

  esum0 = zero
  esumr = zero
  fsumr(:,:) = zero

!  esum0 = esum0 - two*eta_sqrtpi

  do ja = 1, natom
     esumr = esumr - zz(ja)*zz(ja)*two*eta_sqrtpi
  enddo


  ! Sum over ia \= ja. 

  esub = zero

  do ja = 1, natom
     do ia = 1, natom
        esub = zero
        fsub = zero
        sum_fact = zz(ja) * zz(ia)
        rp(3) = clust%zatm(ja) - clust%zatm(ia)

        do ii = -imx, imx          
           rp(1) = real(ii,dp) + ratom(1,ja) - ratom(1,ia) 
           do jj = -jmx, jmx
              rp(2) = real(jj,dp) + ratom(2,ja) - ratom(2,ia)
                          
              call matvec3('N',adot,rp,pvec)
              call matvec3('N',pbc%latt_vec,rp,tmpvec)
              rmod = sqrt(dot_product(rp,pvec))
              
              ! avoid rmod = 0
              if (rmod == 0) cycle
              
              arg = rmod * eta
              if (arg < small2) then
                 tmp = erfc(arg) / rmod
                 tmp2 = (tmp + two*eta_sqrtpi*exp(-arg*arg))/(rmod*rmod)
                 esub = esub + tmp
                 fsub = fsub + tmpvec * tmp2
              endif
           enddo
        enddo

        esum0 = esum0 + sum_fact * esub
        fsumr(:,ia) = fsumr(:,ia) - sum_fact * fsub
!        fsumr(:,ja) = fsumr(:,ja) - sum_fact * fsub

     enddo
  enddo

  ! The Ewald energy is the sum of all the above summations

  write(9,*) 'ewald terms', esumg, esumm, esum0, esumr

  eewald = esumg + esumm + esum0 + esumr

  !eewald = eewald/2
 
  !
  ! Set the ionic forces and print if flaged
  ! note - returned force is in units of real space lattice vectors
  ! printed force is in cartesian coordinates

  do ii = 1, 3
     clust%force(ii,:) = 2*(fsumr(ii,:)+fsumm(ii,:)+fsumg(ii,:))
!     do jj = 1, 3
!        clust%force(ii,:) = clust%force(ii,:) + &
!             pbc%bdot(ii,jj)*fsumg(jj,:)*twopio
!     enddo
  enddo

!  do i = 1, clust%atom_num
!     call matvec3('N',pbc%latt_vec,clust%force(1,i),tau(1,i))
!     clust%force(:,i)=tau(:,i)
!  enddo


  ! Calculate forces and energy from sheet-nuclei interaction
  ! when charged sheets are present.
  if (clust%has_charged_sheet) then
     call force_charged_sheet(clust,pbc,zv,ipr,enucsh)
     eewald = eewald + enucsh
  endif


  ! if print flag on, print the ion-ion forces
  if(ipr >= 1) then
     write(7,*)
     write(7,*) 'Forces from ion-ion interaction:'
     write(7,*) '================================'
     write(7,22)
     write(7,*)
     do i = 1, clust%atom_num
        write(7,20) i,clust%xatm(i),clust%yatm(i),clust%zatm(i),clust%force(:,i)
     enddo
     write(7,*)
     write(7,*) 'fsumg'
     do i=1, clust%atom_num
        write(7,*) fsumg(:,i)
     enddo
     write(7,*) 'fsumr'
     do i=1, clust%atom_num
        write(7,*) fsumr(:,i)
     enddo
     write(7,*) 'fsumm'
     do i=1, clust%atom_num
        write(7,*) fsumm(:,i)
     enddo
  endif


20 format(i4,2x,3(f11.6,1x),2x,3(f11.6,1x))
22 format('-atom-  ----x----   ----y----   ----z-----' &
        ,'   ----Fx----  ----Fy----   ----Fz---',/ &
        '                       [bohr]           ' &
        ,'                  [Ry/bohr]            ')



end subroutine ewald_slab

!===================================================================

function erf(x)
  use constants
  implicit none

  real(dp), intent(in) :: x
  real(dp) :: erf
  real(dp), external :: erfc


  erf = 1 - erfc(x)

end function erf

!===================================================================

subroutine fplus(vecsize,zdiff,eta,output)
  use constants
  implicit none
  
  real(dp), intent(in) :: vecsize
  real(dp), intent(in) :: zdiff
  real(dp), intent(in) :: eta
  real(dp), intent(out) :: output
  real(dp), external :: erfc
  
  real(dp) :: mult     ! mult = vecsize*zdiff
  real(dp) :: etasqr   ! etasqr = 2*(eta**2)
  real(dp) :: inveta   ! inveta = 1/(2*eta)

  mult = vecsize * zdiff
  etasqr = 2*(eta**2)
  inveta = 0.5d0 / eta

  output = exp(-mult) * erfc(inveta*(vecsize - etasqr*zdiff)) + &
           exp(mult)  * erfc(inveta*(vecsize + etasqr*zdiff))

  !if(ISNAN(output)) output=0

  return
end subroutine fplus
!==================================================================

subroutine fminus(vecsize,zdiff,eta,output)
  use constants
  implicit none

  real(dp), intent(in) :: vecsize
  real(dp), intent(in) :: zdiff
  real(dp), intent(in) :: eta
  real(dp), intent(out) :: output
  real(dp), external :: erfc

  real(dp) :: mult     ! mult = vecsize*zdiff
  real(dp) :: etasqr   ! etasqr = 2*(eta**2)
  real(dp) :: inveta   ! inveta = 1/(2*eta)

  mult = vecsize * zdiff
  etasqr = 2*(eta**2)
  inveta = 0.5d0 / eta

  output = exp(-mult) * erfc(inveta*(vecsize - etasqr*zdiff)) - &
           exp(mult)  * erfc(inveta*(vecsize + etasqr*zdiff))

  !if(ISNAN(output)) output=0.0

  return
end subroutine fminus
!==================================================================
