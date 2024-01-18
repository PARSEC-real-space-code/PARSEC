! ===============================================================
!
! Calculates the local spin distribution by integrating the spin density
! (imbalance between spin-up electron density and spin-down electron
! density) over a Wigner-Seitz cell centered on each atom site
! In periodic systems, usual wrapping-around is done.
! Atomic coordinates are input in file angproj.in.
!
! Input files:
!     parsec.dat
!     angproj.in
! Output files:
!     wz_spin.dat
!
! Copyright (C) Murilo Tiago (Univ. of Texas, February 2007).
!          mtiago@ices.utexas.edu.
! This is free software; you can redistribute it and/or modify it
! under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 2 of the License, or
! (at your option) any later version.
! This code is distributed in the hope that it will be useful, but
! WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
! General Public License for more details.
! You can receive a copy of the GNU General Public License
! by writing to the Free Software Foundation,
! Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
!
! ---------------------------------------------------------------
!
! Constants
!
module constants
  integer, parameter :: dp = kind(1.0d0)
end module constants
!
! ---------------------------------------------------------------
program wz_spin

  use constants
  implicit none

  integer :: &
       nspin, &               ! number of spins
       ndim, &                ! number of grid points (irreducible wedge)
       ntrans, &              ! number of symmetry operations
       natom, &               ! number of atoms to project onto
       pbcflag                ! number of periodic directions (0, 1, 2, or 3)

  real(dp) :: &
       hh(3), &               ! grid spacing for x,y,z directions
       rsize(3), &            ! size of grid for x,y,z directions
       shift(3), &            ! grid shift from wfn.dat
       avec(3,3), &           ! unit lattice vectors
       avec_inv(3,3), &       ! inverse of unit lattice vectors
       avec_norm(3,3), &      ! normalized lattice vectors
       hcub                   ! unit volume (volume around each grid point)

  integer, allocatable :: &
       rgrid(:,:), &          ! coordinates of grid points, in units of hh
       indx(:)                ! ordered list of atoms around a grid point

  real(dp), allocatable :: &
       ratm(:,:), &           ! atomic coordinates
       trans(:,:,:), &        ! rotation matrices from symmetry operations
       rho(:,:), &            ! charge densities
       dist(:), &             ! distance array
       proj_at(:,:)           ! projection onto each atom

  ! tolerance in the distance between grid points
  real(dp), parameter :: dtol = 1.d-8
  ! counters
  integer ii, jj, isp, jr, iat, itrans, nat_wz, wz_bound, &
       wz_weight, stype(10)
  real(dp) :: rr(3), rt(3), tmp, dmin

  ! ---------------------------------------------------------------
  !
  ! Read atom coordinates from angproj.in.
  !
  open(10,file='angproj.in',status='old')
  read(10,*) natom

  allocate(ratm(3,natom))
  do iat = 1, natom
     read(10,*) (ratm(ii,iat),ii=1,3)
  enddo
  close(10)
  write(6,'(a,i8)') ' number of atoms = ',natom
  !
  ! Start reading parsec.dat.
  !
  open(20,file='parsec.dat',form='unformatted')

  read(20)
  read(20) stype
  nspin = stype(1)

  if (stype(3) == 0) then
     ! Confined system.
     pbcflag = 0
     write(6,*) ' Local spin distribution for confined system'
     read(20) hh(1), tmp
     hh(2:3) = hh(1)
     avec_norm = 0.d0
     ! Define unit volume and metric as identity matrices.
     do ii = 1, 3
        avec_norm(ii,ii) = 1.d0
     enddo
     avec = avec_norm
     avec_inv = avec_norm
     hcub = hh(1) * hh(2) * hh(3)
  else
     select case( stype(3) )
     case(1)
        ! Bulk (3-d) periodic system.
        pbcflag = 3
        write(6,*) ' Local spin distribution for bulk system'
     case(2)
        ! Wire (1-d) periodic system.
        pbcflag = 1
        write(6,*) ' Local spin distribution for wire (1-d) system'
     end select
     read(20) hh(:),rsize(:)
     read(20) ((avec(ii,jj),ii=1,3),jj=1,3)
     do ii = 1, 9
        read(20)
     enddo
     ! Define unit volume and metric.
     avec_inv = avec
     call mtrxin(avec_inv,tmp,tmp)
     do ii = 1, 3
        tmp = sqrt(sum(avec(:,ii)**2))
        avec(:,ii) = avec(:,ii)/tmp
     enddo
     avec_norm = avec
     call mtrxin(avec,hcub,tmp)
     hcub = hcub * hh(1) * hh(2) * hh(3)
     avec = avec_inv
     call mtrxin(avec,tmp,tmp)
  endif
  ! Set cell size to zero along non-periodic directions
  rsize(pbcflag+1:3) = 0.d0

  !
  ! Convert atomic coordinates to units of lattice vectors
  !
  do iat = 1, natom
     ratm(:,iat) = matmul(avec_inv,ratm(:,iat))
  enddo

  read(20)
  read(20) shift
  read(20) ndim, ntrans
  allocate(trans(3,3,ntrans))
  read(20)
  read(20) trans
  read(20)
  read(20)
  read(20)
  allocate(rgrid(3,ndim))
  read(20) ( (rgrid(ii,jr),ii=1,3), jr=1,ndim )
  !
  ! Initialize arrays and read electron densities for spin up/down.
  !
  allocate(rho(ndim,nspin))
  allocate(proj_at(nspin,natom))
  proj_at = 0.d0

  do isp = 1, nspin
     read(20) ii
     write(6,*) ' read spin ',isp,' number of orbitals = ',ii
     do ii = 1, 4 
        read(20)
     enddo
     read(20) (rho(jr,isp), jr=1,ndim)
     write(6,*) ' total charge = ', sum(rho(:,isp)) * hcub * ntrans
  enddo                     ! isp = 1, nspin
  close(20)
  !
  ! Start calculation.
  !
  allocate(dist(natom))
  allocate(indx(natom))
  wz_bound = 0
  wz_weight = 0

  do jr = 1, ndim
     do itrans = 1, ntrans
        do iat = 1, natom
           ! Get the position of current atom with respect to each grid point.
           rr = rgrid(:,jr) + shift
           rt = matmul(rr,trans(:,:,itrans))
           ! For periodic directions, ratm is given in coordinates of lattice
           ! vectors. Also, atom coordinates must be wrapped around the cell.
           do ii = 1, pbcflag
              rt(ii) = rt(ii) / rsize(ii)
              rr(ii) = ratm(ii,iat) - rt(ii)
              rr(ii) = mod( rr(ii) + 8.d0 + 0.5d0 , 1.d0 ) - 0.5d0
              rr(ii) = rr(ii) * rsize(ii)
              rt = matmul(avec,rr)
              dist(iat) = sqrt( DOT_PRODUCT(rt,rt) )
           enddo
           ! For non-periodic directions, ratm is given in Cartesian
           ! coordinates.
           do ii = pbcflag + 1, 3
              rr(ii) = ratm(ii,iat) - hh(ii) * rt(ii)
           enddo
           dist(iat) = sqrt( DOT_PRODUCT(rr,rr) )
        enddo

        call heapsort(natom,dist,indx)
        dmin = dist(indx(1))

        nat_wz = 0
        do iat = 1, natom
           if (abs(dist(indx(iat)) - dmin) < dtol) then
              nat_wz = nat_wz + 1
           endif
        enddo

        if (nat_wz > 1) then
           wz_bound = wz_bound + 1
           wz_weight = wz_weight + nat_wz
           do iat = 1, natom
              if (abs(dist(indx(iat)) - dmin) < dtol) then
                 do isp = 1, nspin
                    proj_at(isp,iat) = proj_at(isp,iat) + &
                         rho(jr,isp) * hcub/real(nat_wz,dp)
                 enddo
              endif
           enddo
        elseif (nat_wz < 1) then
           write(6,*) ' ERROR: could not find nearest atom to ' , &
                ' grid point # ',jr,rgrid(:,jr)
        else
           iat = indx(1)
           do isp = 1, nspin
              proj_at(isp,iat) = proj_at(isp,iat) + rho(jr,isp) * hcub
           enddo
        endif
        
     enddo                  ! itrans = 1, ntrans
  enddo                     ! jr = 1, ndim
  deallocate(rho)
  if (wz_bound > 0) then
     write(6,*) ' Number of grid points on the surface of a ', &
          'Wigner-Seitz cell = ',wz_bound
     write(6,*) ' Average weight of points on the surface = ', &
          real(wz_bound,dp)/real(wz_weight,dp)
  endif
  !
  ! Print out projections.
  ! (note: projections are not weighted by occupancy factor!)
  !
  open(30,file='wz_spin.dat')
  do iat = 1, natom
     write(30,'(i5,3x,4(g12.4))') iat, proj_at(:,iat), &
          proj_at(1,iat) - proj_at(2,iat)
  enddo
  do isp = 1, nspin
     write(6,*) ' spin ', isp, ' total charge = ', sum(proj_at(isp,:))
  enddo
  close(30)

end program wz_spin
! ===============================================================
!
!     SORTS AN ARRAY BY THE HEAPSORT METHOD
!     W. H. PREUSS ET AL. NUMERICAL RECIPES
!
!-------------------------------------------------------------------
subroutine heapsort(n, arrin, indx)
  use constants
  implicit none

  ! arguments
  integer, intent(in) :: n
  real(dp), intent(in) :: arrin(n)
  integer, intent(out) :: indx(n)

  ! local variables
  integer :: i, j, l, indxt, ir
  double precision q

  if (n == 0) return  
  if (n == 1) then  
     indx(1) = 1  
     return  
  end if

  do j = 1, n  
     indx(j) = j  
  end do

  l = n / 2 + 1  
  ir = n

  do
     if (l > 1) then  
        l = l - 1  
        indxt = indx(l)  
        q = arrin(indxt)  
     else  
        indxt = indx(ir)  
        q = arrin(indxt)  
        indx(ir) = indx(1)  
        ir = ir - 1  
        if (ir == 1) then  
           indx(1) = indxt
           return  
        end if
     end if
     i = l
     j = l + l
     do while (j <= ir)
        if (j < ir) then  
           if (arrin(indx(j)) < arrin(indx(j + 1))) j = j + 1
        end if
        if (q < arrin(indx(j))) then  
           indx(i) = indx(j)  
           i = j
           j = j + j  
        else  
           j = ir + 1  
        end if
     end do
     indx(i) = indxt
  end do

end subroutine heapsort
!===================================================================
!
! Inverts 3x3 matrix m corresponding to a symmetry operation,
! storing the result in m. It also calculates determinant and
! trace of the input matrix. Matrix inversion is aborted if
! det<del.
!
!---------------------------------------------------------------
subroutine mtrxin(m,det,tr)

  use constants
  implicit none
  !
  ! Input/Output variables:
  !
  real(dp), intent(inout) :: m(3,3)
  real(dp), intent(out) :: det,tr
  !
  ! Work variables:
  !
  real(dp) :: a(3,3),del,x
  integer i,j
  !---------------------------------------------------------------
  !
  ! compute matrix of cofactors
  !
  a(1,1) = m(2,2)*m(3,3) - m(2,3)*m(3,2)
  a(2,1) = -m(2,1)*m(3,3) + m(2,3)*m(3,1)
  a(3,1) = m(2,1)*m(3,2) - m(2,2)*m(3,1)
  a(1,2) = -m(1,2)*m(3,3) + m(1,3)*m(3,2)
  a(2,2) = m(1,1)*m(3,3) - m(1,3)*m(3,1)
  a(3,2) = -m(1,1)*m(3,2) + m(1,2)*m(3,1)
  a(1,3) = m(1,2)*m(2,3) - m(1,3)*m(2,2)
  a(2,3) = -m(1,1)*m(2,3) + m(1,3)*m(2,1)
  a(3,3) = m(1,1)*m(2,2) - m(1,2)*m(2,1)
  !
  ! compute determinant
  !
  det = m(1,1)*a(1,1) + m(1,2)*a(2,1) + m(1,3)*a(3,1)
  tr = m(1,1) + m(2,2) + m(3,3)
  del = 1.0d-05
  if (abs(det) < del) stop 501
  !
  ! form mi
  !
  do i=1,3
     do j=1,3
        x = a(i,j)/det
        m(i,j) = x
     enddo
  enddo

end subroutine mtrxin
!===============================================================
