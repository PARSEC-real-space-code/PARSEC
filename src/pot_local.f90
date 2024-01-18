!===============================================================
!
! Copyright (C) 2005 Finite Difference Research Group
! This file is part of parsec, http://www.ices.utexas.edu/parsec/
!
! Computes the local potential vscr on the real space using
! fast fourier transforms.
!
! Adapted from plane-wave programs written by S. Froyen and
! J. L. Martins.
!
!---------------------------------------------------------------
subroutine pot_local(pbc,ipr,veff)

  use constants
  use pbc_module
  implicit none
  !
  ! Input/Output variables:
  !
  ! pbc related data
  type (pbc_data), intent(inout) :: pbc
  ! local potential in reciprocal space
  complex(dpc), dimension (pbc%ng), intent(in) :: veff
  ! printout flag
  integer, intent(in) :: ipr
  !
  ! Work variables:
  !
  integer ntot,i,j,k,ijk2,k1,k2,k3
  integer iadd,ierr

  real(dp) :: dmax,dmin,cmax,abschd
  complex(dpc) :: phase

  real(dp), parameter :: small = 1.0d-9

  !---------------------------------------------------------------

  ! Print out local potential if requested.
  if (ipr >= 6) write(7,10) (real(veff(i)),aimag(veff(i)),i=1,pbc%ng)

  ntot = pbc%n1 * pbc%n2 * pbc%n3
  write(7,14) pbc%n1,pbc%n2,pbc%n3
  !
  ! Initialize charge density array and enter symmetrized
  ! charge. Must add phase from displacement of FFT origin:
  ! coordinates of corner of FFT box are not (0,0,0).
  !
  pbc%vscr4(:) = zero

  do i=1,pbc%ng
     k1 = 1 + pbc%kgv(1,i)
     if (k1 <= 0) k1 = k1 + pbc%n1
     k2 = 1 + pbc%kgv(2,i)
     if (k2 <= 0) k2 = k2 + pbc%n2
     k3 = 1 + pbc%kgv(3,i)
     if (k3 <= 0) k3 = k3 + pbc%n3
     call get_address(pbc,k1,k2,k3,iadd,phase)
     pbc%vscr4(iadd) = pbc%vscr4(iadd)  + veff(i) * phase
  enddo

  ! Fourier transform to real space.
  call cfftw (pbc%vscr4(:),pbc%n1,pbc%n2,pbc%n3,-1)

  dmax = real( pbc%vscr4(1) ,dp)
  dmin = dmax
  cmax = zero
  ierr = 0
  do k=1,pbc%n3
     do j=1,pbc%n2
        do i=1,pbc%n1
           ijk2 = ((k-1)*pbc%n2 + j-1)*pbc%n1 + i
           pbc%vscr2(ijk2) = real( pbc%vscr4(ijk2) ,dp)
           if (pbc%vscr2(ijk2) > dmax) dmax=pbc%vscr2(ijk2)
           if (pbc%vscr2(ijk2) < dmin) dmin=pbc%vscr2(ijk2)
           abschd = abs( aimag( pbc%vscr4(ijk2) ) )
           if (abschd > cmax) cmax=abschd
           if (abschd > small) ierr=ierr+1
        enddo
     enddo
  enddo

  write(7,16) dmax,dmin
  if (ierr /= 0) write(7,18) cmax

10 format(/,'  local potential in g-space',/,100(/,1x,20f6.2))
14 format(/,'  in fft for local potential n = ',3i4)
16 format(/,'  max and min of potential ',2f10.4)
18 format(/,'  complex potential: max value:',e12.4)

end subroutine pot_local
!===============================================================
