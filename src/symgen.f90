!===============================================================
!
! Copyright (C) 2005 Finite Difference Research Group
! This file is part of parsec, http://www.ices.utexas.edu/parsec/
!
! Adapted by J.L. Martins from the program GROUP
! written in 1974 by Warren and Worlton,
! Computer Physics Communications, vol 8, 71-74 (1974)
! incorporated by Eckold et al into UNISOFT
! Modified by Alberto Garcia (1990)
! Cleaned up by Bernd Pfrommer (1996)
! Abelian group routines added by M. Tiago (2005)
!
! Let's see if we can get the rotations both in real and
! in reciprocal space:
!
! gmtrx : g-space representation
! rmtrx: r-space representation
!
! Input:
!
! a(i,j) is the i-th cartesian component of the j-th primitive
! translation vector of the direct lattice and thus it is the
! transpose of the matrix A defined by Jones in 'The Theory
! of Brillouin Zones and Electronic States in Crystals'
!
! coorat(i,j,k) is the k-th component (lattice coordinates) of
! the position of the j-th atom of type i.
!
! natom(i) is the number of atoms of type i.
!
! ntype is the total number of types of atoms.
!
! nops indicates whether or not to count nonzero fractional
! translations
!
! internal:
!
! b contains the reciprocal lattice vectors in the
! crystallographic usage, that is, WITHOUT the 2pi factor.
! This matrix IS Jones'.
!
! na is the number of atoms.
!
! ity(i) is an integer distinguishing atoms of different type
! i.e. different atomic species.
!
! x(j,i) is the j-th cartesian component of the position vector for
! the i-th atom in the unit cell.
!
! output:
!
! ntrans is the number of point group operations.
!
! gmtrx is the matrix of rotations in g-lattice coordinates.
! rmtrx is the matrix of rotations in r-lattice coordinates.
!
! tnp(oper,ilat) is the fractional translation in latt. coor.
!
! invers_no is the operation number of the inversion (0 if not
! present). It is used to restore the inversion symmetry by
! a simple change of origin when possible
!
!---------------------------------------------------------------
subroutine symgen(ipr, a, coorat, natom, ntype, invers_no,nops, &
     norot, ntrans, gmtrx, tnp, mxdatm, rmtrx, trans, nab, indx, chi, &
     cell_symmetry)

  use constants
  implicit none  
  !
  ! Input/Output variables:
  !
  integer, intent(in) :: ipr, & ! ipr>=2 print
       mxdatm                   ! array dimensioning

  integer, intent(out) :: ntrans, cell_symmetry
  real(dp), intent(out) :: tnp(48, 3), trans(3, 3, 48)
  ! rotation matrices in cartesian coordinates
  integer, intent(out) ::   rmtrx(48, 3, 3), gmtrx(48, 3, 3)
  ! number of symmetry operations in Abelian subgroup
  integer, intent(out) :: nab
  ! index of symmetry operations in Abelian subgroup
  integer, intent(out) :: indx(12)
  ! character table of Abelian subgroup
  complex(dpc), intent(out) :: chi(12,12)

  integer, intent(in) :: ntype, nops, norot
  integer, intent(out) :: invers_no
  real(dp), intent(in) :: a(3, 3), coorat(ntype, mxdatm, 3)
  integer, intent(in) :: natom(ntype)
  !
  ! Work variables:
  !
  real(dp) :: xdum,ydum
  integer :: i, ihg, ind, ipm, j, k, l, li, m, na, nat, ierr, &
       include_fractional, ihc

  real(dp) :: b(3, 3), r(49, 3, 3), r1(3, 3), rlat(48, 3, 3)
  real(dp), allocatable :: x(:,:)
  integer, allocatable :: ity(:)
  integer :: ib(48), mult(48,48)
  character(len=3) :: id(48), label(12)
  character(len=5) :: name
  !
  ! External subroutines:
  !
  external atftmt, pgl, symm_ident, symchk

  !---------------------------------------------------------------

  ipm = 1
  include_fractional = nops

  allocate(x(3, ntype * mxdatm))
  allocate(ity(ntype * mxdatm))
  !
  ! Calculate cartesian coordinates, atom types, and number of atoms.
  !
  na = 0
  do i = 1, ntype
     nat = natom(i)
     do j = 1, nat
        ind = na + j
        do k = 1, 3
           x(k, ind) = zero
           do l = 1, 3
              x(k, ind) = x(k, ind) + a(k, l) * coorat(i, j, l)
           end do
        end do
        ity(ind) = i
     end do
     na = na + natom(i)
  end do
  !
  ! Determine reciprocal lattice basis vectors.
  ! We know that
  !
  !           T                             T -1         -1
  !          B A  = 2pi (1), so   B  = 2pi(A )   = 2pi(a )
  !
  ! and we can use the linpack (SCILIB version) routines
  ! sgefa and sgedi to compute the determinant (celvol) and the
  ! inverse.  But note that the 2pi factor is NOT needed here.
  !
  do i = 1, 3
     do j = 1, 3
        b(i, j) = a(i, j)
     end do
  end do

  call mtrxin(b,xdum,ydum)

  call pgl(a, b, r, ntrans, ib, ihg, ihc)
  cell_symmetry = 1 - ihc  ! 0 for cubic, 1 for hexagonal

  call atftmt(ipr, b, x, r, tnp, trans, ity, na, ib, ihg, &
       include_fractional, li, ntrans, invers_no, ntype, mxdatm)
  !
  ! We have the rotations in cartesian coordinates.
  ! Transform into lattice coordinates (r-space and g-space)
  !
  do l = 1, ntrans
     !
     ! In terms of the real-space basis vectors:
     !                  T
     !          y' = ( B R a ) y     ( y are the r-lattice coord.)
     !
     !                          T
     !          watch out: b = B
     !
     ! Trans * a ...
     !
     do j = 1, 3
        do k = 1, 3
           r1(j, k) = zero
           do m = 1, 3
              r1(j, k) = r1(j, k) + trans(m, j, l) * a(m, k)
           end do
        end do
     end do
     !
     ! B * Trans * a
     !
     do j = 1, 3
        do k = 1, 3
           rlat(l, j, k) = zero
           do m = 1, 3
              rlat(l, j, k) = rlat(l, j, k) + b(j, m) * r1(m, k)
           end do
           rmtrx(l, j, k) = nint(rlat(l, j, k))
        end do
     end do
  end do

  call symm_ident(ntrans, rmtrx, tnp, id)

  do l = 1, ntrans
     !
     ! In terms of the g-space basis vectors:
     !                  T  T
     !          z' = ( a  R  b ) z    ( z are the g-lattice coord.)
     !
     do k = 1, 3
        gmtrx(l, :, k) = rmtrx(l, k, :)
     end do
  end do
  !
  ! write the matrices and fractional translations
  !
  write(7, 900)
900 format(//4x,'Rotation matrices (r-lattice) and fractional', &
         ' translations (r-lattice)',/)

  do i = 1, ntrans
     write(7, 901) i, ((rmtrx(i, j, k), k = 1, 3), j = 1, 3), &
          (tnp(i, k), k = 1, 3), id(i)
     
  end do
901 format(i5,3(3x,3i3),4x,3f9.5,1x,a5)

  ierr = 0
  if (norot >= 0) call symchk(ipr,ierr,ntrans,rmtrx,tnp,mult)
  if (ierr /= 0) ntrans = 1
  if (norot < 0) then
     nab = 1
     indx(1) = 1
     chi(1,1) = zone
     goto 950
  endif
  call abelian(ntrans, norot, mult, rmtrx, tnp, id, nab, indx, chi, &
       name, label)
  !
  ! write symmetry operations for the Abelian subgroup
  !
  write(7,*)
  write(7,*) ' Found Abelian subgroup ',name,' with ',nab, &
       ' symmetry operations'
  write(7,*)
  write(7,*) ' Character table : character = exp(i*Pi*x/3)'
  write(7,'(100a)') ('-',i=1,18 + 5*nab)
  write(7,'(a,12(a5))') ' Representation |', (trim(id(indx(k))),k=1,nab)
  do i = 1, nab
     write(7,'(3x,a3,1x,i5,2x,a,12(2x,i3))') label(i),i,' | ', &
          (mod(nint(three*aimag(log(chi(i,k)))/pi),6), k = 1, nab)
  enddo
  write(7,'(100a)') ('-',i=1,18 + 5*nab)
  write(7,'(a,1x,i3,11(2x,i3))') ' Operation number:', &
       (indx(k), k = 1, nab)

950 continue
  write(7,*)
  !
  ! If identity is not a symmetry operation, then the symmetry group is wrong
  !
  k = 0
  do i = 1, ntrans
     if (trim(id(i)) == 'E') k = 1
  enddo
  if (k == 0) then
     write(7,*)
     write(7,*) ' WARNING! SYMMETRY GROUP DOES NOT HAVE IDENTITY'
     write(7,*) ' GROUP SYMMETRIES ARE INCORRECT, POSSIBLY BECAUSE'
     write(7,*) ' OF OVERLAPPING ATOMS. STOP...'
     write(7,*)
     ntrans = 0
  endif

  deallocate(x)
  deallocate(ity)

end subroutine symgen
!===============================================================
!
! Determines the point group of the crystal,
! the atom transformation table,f0, the fractional translations,
! tnp, associated with each rotation and finally the
! multiplication table mt, for the point group of the crystal. ib
! now contains operations in the p.g. of the crystal and ntrans
! is the order of this group.
!
! 1997 Bernd Pfrommer, based on a routine from Alberto Garcia's
! code. I cleaned up some more, put in dynamic memory allocation,
! explicit typing, fortran 90 style, and more comments. Also,
! some shortcuts were put in to speed up the case where the
! symmetry is low. The rdiff array precomputes the rotated
! vectors to speed up things.
!
!---------------------------------------------------------------
subroutine atftmt(ipr, ai, x, r, tnp, trans, ity, na, ib, &
     ihg, include_fractional, li, nc, invers_no, mxdtyp, mxdatm)

  use constants
  implicit none
  !
  ! Input/Output variables:
  !
  integer, intent(in) :: &
       ipr, &               ! print flag
       ihg, &               ! holohedral group number
       na, &                ! total number of atoms (of all kinds)
       mxdtyp, &            ! max number of atom types
       mxdatm, &            ! max number of atoms of each type
       include_fractional   ! -2 means don't include frac translations

  real(dp), intent(in) :: &
       x(3, mxdtyp * mxdatm), & ! compact list of all coordinates (cartesian)
       ai(3, 3)                 ! inverse of lattice vectors (no 2pi)

  integer, intent(inout) :: &
       ity(mxdtyp * mxdatm), & ! compact list of all the types of all atoms
       ib(48), &               ! index map for symmetry operations
       nc                      ! number of symm-ops without/with basis
  real(dp), intent(inout) :: &
       r(49, 3, 3)          ! the rotation matrices as they come
                            ! out of the pgl subroutine

  integer, intent(out) :: &
       li, &                ! something to do with inversion?
       invers_no            ! which operation is the inversion
  real(dp), intent(out) :: &
       trans(3, 3, 48), &   ! the cartesian realspace transf. matrices
       tnp(48, 3)           ! the nonprimitive translations
  !
  ! Work variables:
  !
  real(dp) :: v(3, 48), vr(3), vt(3), xb(3)
  integer :: ia(48), ic(48), mt(48, 48), &
       ipm                  ! print flag for multiplication table

  real(dp) :: da, dif, ts, vs, difmax
  real(dp), parameter :: eps = 1.0d-8
  integer :: i, il, is, isy, iu, j, k, k1, k2, k3, k4, ks, l, &
       m, n,n1, n2, n3, nca, ni
  real(dp), allocatable :: rx(:,:), rdiff(:,:,:)
  integer, allocatable :: if0(:,:)
  character(len=12), parameter :: cst(7) =(/ 'triclinic   ', &
       'monoclinic  ', 'orthorhombic','tetragonal  ', &
       'cubic       ', 'trigonal    ', 'hexagonal   ' /)

  !---------------------------------------------------------------

  invers_no = 0
  ipm = 1

  allocate(rx(3, mxdtyp * mxdatm))
  allocate(if0(48, mxdtyp * mxdatm))
  allocate(rdiff(3, na, na))
  !
  ! eps should be slightly larger than computer precision
  !
  nca = 0
  ni = 13
  if (ihg < 6) ni = 25
  li = 0
  op_loop : do n = 1, nc    ! loop over all lattice symmetry operations
     l = ib(n)              ! get index of symmetry operation
     ic(n) = ib(n)
     !
     ! Operate on all atoms with symmetry operation l, and store in
     ! list rx.
     !
     do k = 1, na
        do i = 1, 3
           rx(i, k) = zero
           do j = 1, 3
              rx(i, k) = rx(i, k) + r(l, i, j) * x(j, k)
           end do
        end do
     end do
     !
     ! This piece of code is pretty naive, and scales like order(n**3).
     ! It basically checks if for each
     !
     ! R*x_1 - x_2   there is a matching R*x_3-x4
     !
     ! excluding the trivial case of x_1 == x_3  and  x_2 == x_4
     !
     ! precompute the rdiffs first
     do k1 = 1, na
        do k2 = 1, na
           xb(:) = rx(:, k1) - x(:, k2)
           ! rdiff = R*x_1 -x_2
           call rlv(ai, xb, rdiff(1, k2, k1), il)
           ! Subroutine rlv removes a direct lattice vector from xb
           ! leaving the remainder in rdiff. if a nonzero lattice vector was
           ! removed, il is made nonzero.
        end do
     end do

     difmax = zero
     do k1 = 1, na          ! double loop over compact atom list
        do k2 = 1, na
           if (ity(k1) == ity(k2)) then ! same type atoms?
              vr = rdiff(:, k2, k1) !vr stands for v-reference.
              ks = 0
              do k3 = 1, na
                 do k4 = 1, na
                    if (ity(k3) == ity(k4)) then
                       vt = rdiff(:, k4, k3) ! vt stands for v-test
                       dif = zero
                       do i = 1, 3
                          da = abs(vr(i) - vt(i)) + eps
                          dif = dif + mod(da, one)
                       end do
                       if (dif <= 10.0d0 * eps) then
                          if0(l, k3) = k4
! if0 is the function defined in maradudin and vosko by eq.(2.35).
! it defines the atom transformation table
                          ks = ks + k4
! if (ks == na * (na + 1) / 2) goto 110       ! found
                          if (ks == na*(na+1)/2) then
! Reject all symops with non-zero fractional translations.  It
! should be helpful when you have weird fractional translations
! that confuse adjustfft.fp.  -- David Roundy
                             if (include_fractional == -1 .or. &
                                  vr(1)*vr(1)+vr(2)*vr(2)+  &
                                  vr(3)*vr(3) < eps * 0.1d0) &
                                  go to 110 ! found all
                          endif
                          goto 80
                       end if
                    end if
                 end do
                 exit
80               continue
              end do

! BP put in this shortcut (check carefully). If there
! is a single R*x_1- x_2 without match, then give up.
! this is not a symmetry operation
              if (ks == 0) cycle op_loop
           end if
        end do
     end do

     ! this was not a symmetry operation
     cycle op_loop

     ! we found a symmetry operation
110  continue

     nca = nca + 1
     !
     ! v(i,l) is the i-th cartesian component of the fractional
     ! translation associated with the rotation r(l).
     !
     v(1:3, l) = vr(1:3)

     ib(nca) = l
     if (l == ni) then
        li = l
        invers_no = nca
     end if

  end do op_loop
  deallocate(rdiff)
  !
  ! ------------ there are no gotos across this line -------------
  !
  !  if (ipr >= 2) then
  if ((ihg == 7 .and. nca == 24) .or. (ihg == 5 .and. nca == 48)) then
     write(7, 901) trim(cst(ihg))
901  format      ( /' The point group of the system is the full ', &
          a12,' group')
  else
     write(7, 900) trim(cst(ihg)), (ic(i), i = 1, nc)
900  format      (/' The crystallographic system is ',a, &
          ' with operations: ',4(/5x,16i3))
  end if
  !  end if
  vs = zero
  nc = nca
  do n = 1, nc
     l = ib(n)
     vs = sum(abs(v(1:3, l)))
  end do
  if (vs > eps) then
     write(7, 903)
903  format   (/' The space group is non-symmorphic',/, &
          ' (Or a non-standard origin of coordinates is used)',/)
     isy = 0
     is = 0
  else
     write(7, 902)
902  format   (/' The space group of the crystal is symmorphic',/, &
          ' (this usually means that symmetry operations do not',/ &
          ,' involve rotations and translations simultaneously)'/)
     isy = 1
     is = 1
  end if
  !
  ! Construct the multiplication table.
  !
  do n1 = 1, nc
     do n2 = 1, nc
        l = ib(n1)
        m = ib(n2)
        do i = 1, 3
           do j = 1, 3
              r(49, i, j) = zero
              do k = 1, 3
                 r(49, i, j) = r(49, i, j) + r(l, i, k) * r(m , k, j)
              end do
           end do
        end do
        do n3 = 1, nc
           n = ib(n3)
           ts = zero
           do i = 1, 3
              do j = 1, 3
                 ts = ts + abs(r(49, i, j) - r(n, i, j))
              end do
           end do
           if (ts > 1.0d2 * eps) cycle
           mt(l, m) = n
           exit
        end do
     end do
  end do

  il = 1
  iu = nc
  if (iu > 24) iu = 24
  write(7,*) 'Operation number: '
280 continue
  write(7, '(5x,12i3)') (ib(i), i = il, iu)
  do i = 1, na
     do j = 1, nc
        l = ib(j)
        ia(j) = if0(l, i)
     end do
  end do
  if (nc <= iu) goto 320
  il = 25
  iu = nc

  goto 280
  !
  ! Print multiplication table and fractional translations.
  !
320 continue
  if (ipm == 0) goto 410
  il = 1
  iu = nc
  if (nc > 24) iu = 24
  if (is > 0) goto 340
  if (ipr >= 2) then
     write(7, 906)
906  format   (/,3x,'Multiplication table',30x,'Fractional translations')
     write(7, 907) (ib(i), i = il, iu)
907  format   (5x,24i4)
     write(7, 908)
908  format   (10x,'v(1)      v(2)      v(3)')
  end if

  goto 360

340 continue
  if (ipr >= 2) then
     write(7, 909)
  end if
909 format(/,3x,'Multiplication table')
350 continue
  if (ipr >= 2) then
     write(7, 910) 'x',(ib(i), i = il, iu)
  end if
910 format(4x,a,24i4)
360 continue
  do j = 1, nc
     l = ib(j)
     do i = il, iu
        n = ib(i)
        ia(i) = mt(l, n)
     end do
     if (is > 0) goto 390
     if (ipr >= 2) then
        write(7, 912) ib(j), (ia(i), i = il, iu)
        write(7, 911) (v(i, l), i = 1, 3)
911     format      (10x,3f10.4)
     end if

     cycle

390  continue
     if (ipr >= 2) then
        write(7, 912) ib(j), (ia(i), i = il, iu)
     end if
  end do
912 format(i5,24i4)
  if (iu == nc) goto 410
  il = 25
  iu = nc
  is = 1

  goto 350

410 continue

  do i = 1, nc
     l = ib(i)
     do j = 1, 3
        tnp(i, j) = -v(j, l)
        do k = 1, 3
           trans(j, k, i) = r(l, j, k)
        end do
     end do
  end do

  deallocate(rx)
  deallocate(if0)

end subroutine atftmt
!===============================================================
!
! Determines the point group of the lattice and
! the crystal system. The array ib contains the locations of the
! group operations and nc is the order of the group.
!
!---------------------------------------------------------------
subroutine pgl(a, b, r, nc, ib, ihg, ihc)

  use constants
  implicit none
  !
  ! Input/Output variables:
  !
  integer, intent(out) :: ihg, ihc, nc
  real(dp), intent(in) :: a(3, 3), b(3, 3)
  real(dp), intent(inout) :: r(49, 3, 3)
  integer, intent(out) :: ib(48)
  !
  ! Work variables:
  !
  real(dp) :: tr
  !
  ! eps should be slightly larger than computer precision
  !
  real(dp), parameter :: eps = 1.0d-8
  integer :: i, j, k, lx, n, nr

  real(dp) :: vr(3), xa(3)

  !---------------------------------------------------------------

  ihc = 0
  !
  ! ihc is 0 for hexagonal groups and 1 for cubic groups.
  !
  nr = 24
10 continue
  nc = 0
  call rot(r, nr)
  do n = 1, nr
     ib(n) = 0
     tr = zero
     do k = 1, 3
        do i = 1, 3
           xa(i) = zero
           do j = 1, 3
              xa(i) = xa(i) + r(n, i, j) * a(j, k)
           end do
        end do
        call rlv(b, xa, vr, lx)
        do i = 1, 3
           tr = tr + abs(vr(i))
        end do
     end do
     if (tr <= 10.0d0 * eps) then
        nc = nc + 1
        ib(nc) = n
     end if
  end do
  if (ihc == 0) then
     if (nc == 12) then
        ihg = 6

        return

     end if
     if (nc > 12) then
        ihg = 7

        return

     end if
     if (nc < 12) then
        nr = 48
        ihc = 1

        goto 10

     end if
  else
     if (nc == 16) then
        ihg = 4

        return

     end if
     if (nc > 16) then
        ihg = 5

        return

     end if
     if (nc < 16) then
        if (nc == 4) then
           ihg = 2

           return

        end if
        if (nc > 4) then
           ihg = 3

           return

        end if
        if (nc < 4) then
           ihg = 1

           return

        end if
     end if
  end if
  !
  ! ihg stands for holohedral group number.
  !
end subroutine pgl
!===============================================================
!
! Define the generators for the rotation matrices
!
!---------------------------------------------------------------
subroutine rot(r, nr)

  use constants
  implicit none
  !
  ! Input/Output variables:
  !
  integer, intent(in) :: nr

  real(dp), intent(out) :: r(49, 3, 3)
  !
  ! Work variables:
  !
  real(dp) :: f
  integer :: i, j, k, n, nv
  !---------------------------------------------------------------
  do n = 1, nr
     do i = 1, 3
        do j = 1, 3
           r(n, i, j) = zero
        end do
     end do
  end do

  if (nr <= 24) then
     !
     !                             --hexagonal group
     f = root3 * half
     r(2, 1, 1) = half
     r(2, 1, 2) = -f
     r(2, 2, 1) = f
     r(2, 2, 2) = half
     r(7, 1, 1) = mhalf
     r(7, 1, 2) = -f
     r(7, 2, 1) = -f
     r(7, 2, 2) = half
     do n = 1, 6
        r(n, 3, 3) = one
        r(n + 18, 3, 3) = one
        r(n + 6, 3, 3) = mone
        r(n + 12, 3, 3) = mone
     end do
     !
     ! generate the rest of the rotation matrices
     !
     do i = 1, 2
        r(1, i, i) = one
        do j = 1, 2
           r(6, i, j) = r(2, j, i)
           do k = 1, 2
              r(3, i, j) = r(3, i, j) + r(2, i, k) * r(2, k, j)
              r(8, i, j) = r(8, i, j) + r(2, i, k) * r(7, k, j)
              r(12, i, j) = r(12, i, j) + r(7, i, k) * r(2, k, j)
           end do
        end do
     end do
     do i = 1, 2
        do j = 1, 2
           r(5, i, j) = r(3, j, i)
           do k = 1, 2
              r(4, i, j) = r(4, i, j) + r(2, i, k) * r(3, k, j)
              r(9, i, j) = r(9, i, j) + r(2, i, k) * r(8, k, j)
              r(10, i, j) = r(10, i, j) + r(12, i, k) * r(3, k, j)
              r(11, i, j) = r(11, i, j) + r(12, i, k) * r(2, k, j)
           end do
        end do
     end do

     do n = 1, 12
        nv = n + 12
        do i = 1, 2
           do j = 1, 2
              r(nv, i, j) = -r(n, i, j)
           end do
        end do
     end do
  else
     !
     !                                      --cubic group
     !
     r(9, 1, 3) = one
     r(9, 2, 1) = one
     r(9, 3, 2) = one
     r(19, 1, 1) = one
     r(19, 2, 3) = mone
     r(19, 3, 2) = one
     do i = 1, 3
        r(1, i, i) = one
        do j = 1, 3
           r(20, i, j) = r(19, j, i)
           r(5, i, j) = r(9, j, i)
           do k = 1, 3
              r(2, i, j) = r(2, i, j) + r(19, i, k) * r(19, k, j)
              r(16, i, j) = r(16, i, j) + r(9, i, k) * r(19, k, j)
              r(23, i, j) = r(23, i, j) + r(19, i, k) * r(9, k, j)
           end do
        end do
     end do
     do i = 1, 3
        do j = 1, 3
           do k = 1, 3
              r(6, i, j) = r(6, i, j) + r(2, i, k) * r(5, k, j)
              r(7, i, j) = r(7, i, j) + r(16, i, k) * r(23, k, j)
              r(8, i, j) = r(8, i, j) + r(5, i, k) * r(2, k, j)
              r(10, i, j) = r(10, i, j) + r(2, i, k) * r(9, k, j)
              r(11, i, j) = r(11, i, j) + r(9, i, k) * r(2, k, j)
              r(12, i, j) = r(12, i, j) + r(23, i, k) * r(16, k,j)
              r(14, i, j) = r(14, i, j) + r(16, i, k) * r(2, k, j)
              r(15, i, j) = r(15, i, j) + r(2, i, k) * r(16, k, j)
              r(22, i, j) = r(22, i, j) + r(23, i, k) * r(2, k, j)
              r(24, i, j) = r(24, i, j) + r(2, i, k) * r(23, k, j)
           end do
        end do
     end do
     do i = 1, 3
        do j = 1, 3
           do k = 1, 3
              r(3, i, j) = r(3, i, j) + r(5, i, k) * r(12, k, j)
              r(4, i, j) = r(4, i, j) + r(5, i, k) * r(10, k, j)
              r(13, i, j) = r(13, i, j) + r(23, i, k) * r(11, k,j)
              r(17, i, j) = r(17, i, j) + r(16, i, k) * r(12, k,j)
              r(18, i, j) = r(18, i, j) + r(16, i, k) * r(10, k,j)
              r(21, i, j) = r(21, i, j) + r(12, i, k) * r(15, k,j)
           end do
        end do
     end do
     do n = 1, 24
        nv = n + 24
        do i = 1, 3
           do j = 1, 3
              r(nv, i, j) = -r(n, i, j)
           end do
        end do
     end do
  end if

end subroutine rot
!===============================================================
!
! Removes a direct lattice vector from g by
! operation p, leaving the remainder in y. If a nonzero lattice
! vector was removed, l is made nonzero.
!
!---------------------------------------------------------------
subroutine rlv(p, g, y, l)

  use constants
  implicit none
  !
  ! Input/Output variables:
  !
  real(dp), intent(in) :: g(3), & ! vector to be multiplied
       p(3, 3)                    ! multiplication matrix, e.g. lattice vectors
  
  real(dp), intent(out) :: y(3) ! mod(multiplied vector,lattice vector)
  integer, intent(out) :: l     ! is nonzero if a lattice vector was removed
  !---------------------------------------------------------------
  call matvec3('N',p,g,y)
  l = sum(nint(abs(y(1:3))))
  y(1:3) = y(1:3) - one * nint(y(1:3))

end subroutine rlv
!===============================================================
!
! symchk checks if the symmetry operations defined by
! mtrx and tnp really form a group.
!
! Jan 7, 1990: AG/ 2pi factor removed ! /
!
!---------------------------------------------------------------
subroutine symchk(ipr, ierr, ntrans, mtrx, tnp, mult)

  use constants
  implicit none
  !
  ! Input/Output variables:
  !
  integer, intent(in) :: ntrans, ipr
  !
  ! ierr: returns 0 if no error, otherwise it returns
  ! the number of the suspected operation.
  integer, intent(out) :: ierr

  real(dp), intent(in) :: tnp(48, 3)
  integer, intent(in) :: mtrx(48, 3, 3)
  ! for operations i, j such that i * j = k, then mult(i,j) = k
  integer, intent(out) :: mult(48,48)
  !
  ! Work variables:
  !
  character(len=80), parameter :: not_a_grp = 'The symmetry'// &
       ' operations do not form a group. Check operation no%i4$\n'

  real(dp) :: ttest
  integer :: i, im, itest, j, k, l, m, maxerr

  integer :: nerr(48)
  integer :: mtest(3, 3)
  !---------------------------------------------------------------

  mult = 0
  ierr = 0
  !
  ! trivial group:
  !
  if (ntrans == 1) then
     mult(1,1) = 1
     return
  endif
  !
  ! check for duplicate operations
  !
  do i = 2, ntrans
     im = i - 1
     double_loop: do j = 1, im
        do k = 1, 3
           if (abs(tnp(i, k) - tnp(j, k)) > 1.0d-8) cycle double_loop
           do l = 1, 3
              if (mtrx(i, k, l) /= mtrx(j, k, l)) cycle double_loop
           end do
        end do
        if (ipr >= 2) then
           write(7, 900) j, i
        end if
900     format(/' symmetry operations',i3,' and',i3,' are equal')
        ierr = i
        write(7,*) ' ERROR: repeated symmetry operations'
        return
     end do double_loop
  end do
  !
  ! construct muliplication table
  !
  do i = 1, ntrans
     nerr(i) = 0
     do j = 1, ntrans
        mult(i, j) = 0

        ! multiply i and j

        do k = 1, 3
           do l = 1, 3
              mtest(k, l) = 0
              do m = 1, 3
                 mtest(k, l) = mtest(k, l) + mtrx(i, k, m) * mtrx(j, m, l)
              end do
           end do
        end do
        !
        ! check for match
        !
        match_loop: do k = 1, ntrans
           do l = 1, 3
              do m = 1, 3
                 if (mtest(l, m) /= mtrx(k, l, m)) cycle match_loop
              end do
           end do
           mult(i, j) = k
        end do match_loop
     end do
  end do
  !
  ! if translations not correct set mult(i,j) to -1
  !
  if (ipr >= 2) then
     write(7, *) 'Checking in r-space'
  end if
  do i = 1, ntrans
     trl_loop: do j = 1, ntrans
        k = mult(i, j)
        if (k == 0) cycle
        do l = 1, 3
           ttest = tnp(k, l) - tnp(j, l)
           do m = 1, 3
              ttest = ttest - mtrx(j, m, l) * tnp(i, m)
           end do
           itest = nint(ttest)
           if (abs(ttest - real(itest, dp)) < 1.0d-4) cycle
           write(7, *) i, j, k, l, abs(ttest - real(itest, dp))
           mult(i, j) = -1

           cycle trl_loop

        end do
     end do trl_loop
  end do
  !
  ! check multiplication table
  !
  do i = 1, ntrans
     do j = 1, ntrans
        if (mult(i, j) > 0) cycle
        nerr(i) = nerr(i) + 1
        nerr(j) = nerr(j) + 1
     end do
  end do
  !
  ! find element with max error
  !
  ierr = 0
  maxerr = 0
  do i = 1, ntrans
     if (nerr(i) <= maxerr) cycle
     maxerr = nerr(i)
     ierr = i
  end do
  if (ipr >= 2) then
     write(7, 901)
901  format(' Multiplication table',/)
903  format(1x,48i2)
     do i = 1, ntrans
        write(7, 903) (mult(i, j), j = 1, ntrans)
     end do
  end if
  if (ierr /= 0) write(7,*) ' ERROR: repeated symmetry operations'

end subroutine symchk
!===============================================================
!
! Identify the symmetry operations
!
!---------------------------------------------------------------
subroutine symm_ident(ntrans, mtrx, tnp, id)

  use constants
  implicit none
  !
  ! Input/Output variables:
  !
  integer, intent(in) :: ntrans
  integer, intent(in) :: mtrx(48, 3, 3)
  real(dp), intent(in) :: tnp(48, 3)
  character(len=3), intent(out) :: id(48)
  !
  ! Work variables:
  !
  integer :: i, j, oper
  logical :: proper
  real(dp) :: det, trace
  real(dp) :: a (3, 3), t (3)
  character(len=2) :: axes(-2:2)
  data axes / 'C2', 'C3', 'C4', 'C6', 'E ' /
  !---------------------------------------------------------------
  do oper = 1, ntrans
     !
     ! Copy the matrix to a double precision format for
     ! processing with EISPACK. Also, copy the non-primitive
     ! translation.
     !
     do i = 1, 3
        do j = 1, 3
           a(i, j) = mtrx(oper, i, j)
           t(j) = tnp(48, j)
        end do
     end do
     !
     ! Compute determinant  and trace
     !
     det = a(1, 1) * a(2, 2) * a(3, 3) + a(2, 1) * a(3, 2) * a(1, &
          3) + a(1, 2) * a(2, 3) * a(3, 1) - a(3, 1) * a(2, 2) &
          * a(1, 3) - a(2, 1) * a(1, 2) * a(3, 3) - a(1, 1) * a(3, &
          2) * a(2, 3)

     proper = (nint(det) == 1)

     trace = a(1, 1) + a(2, 2) + a(3, 3)

     if (proper) then
        !
        ! Proper operation
        !
        id(oper) = axes(nint(trace - 1))

     else
        !
        ! R = IS , where S is proper
        !
        a = a*mone
        id(oper) = 'I'//axes(nint(-trace - 1))

     end if
  end do

end subroutine symm_ident
!===============================================================
!
! For a given point group, search for the Abelian subgroup with
! maximum number of elements and REAL characters. This is a
! search by inspection, and selects only subgroups with real
! characters. These are the only possibilities: C_1, C_s, C_i,
! C_2, C_2h, C_2v, D_2, and D_2h. Output is in array indx, with
! the indices of symmetry operations in the subgroup with respect
! to the full point group.
!
!---------------------------------------------------------------
subroutine abelian(ntrans, norot, mult, rmtrx, tnp, id, nab, indx, &
     chi, name, label)

  use constants
  implicit none
  !
  ! Input/Output variables
  !
  integer, intent(in) :: &
       ntrans, &            ! number of operations in point group
       norot, &             ! if positive, Abelian group has only E, IE
       rmtrx(48,3,3), &     ! rotation matrix
       mult(48,48)          ! multiplication table of point group
  ! translation vectors
  real(dp), intent(in) :: tnp(48, 3)
  ! identity of symmetry operations, as determined by symm_ident
  character(len=3), intent(in) :: id(48)

  integer, intent(out) :: &
       nab, &               ! number of operations in Abelian group
       indx(12)             ! index of operations in Abelian group
  complex(dpc), intent(out) :: chi(12,12) ! character table
  ! name of Abelian subgroup and their representations
  character(len=5), intent(out) :: name
  character(len=3), intent(out) :: label(12)
  !
  ! Work variables
  !
  ! flags, indicate wether a particular symmetry is present or not
  logical :: &
       l_i, &               ! inversion symmetry
       l_c2, &              ! C2 rotation
       l_c22, &             ! multiple C2 rotations
       l_s2, &              ! IC2 mirror plane orthogonal to C2
       l_c3, &              ! C3 rotation
       l_s3, &              ! IC2 mirror plane orthogonal to C3
       l_c6                 ! C6 rotation
  !
  ! indices, locate flagged symmetries
  integer i_e, i_i, i_c2, i_c22, i_s2, i_c3, i_s3, i_c6
  ! counters
  integer it, jt, ii, jj, ntr, itr(ntrans), tmtrx(3,3)
  ! length of translation vectors
  real(dp) :: length
  ! tolerance
  real(dp), parameter :: tol = 1.d-10
  !
  ! constants for the phases
  complex(dpc) :: ph1,ph2,ph3,phm3,ph6,phm6

  !---------------------------------------------------------------

  l_i = .false.
  l_c2 = .false.
  l_c22 = .false.
  l_s2 = .false.
  l_c3 = .false.
  l_s3 = .false.
  l_c6 = .false.
  ph1 = exp(two*pi*zi)      ! 1
  ph2 = exp(pi*zi)          ! -1 = exp(i * Pi)
  ph3 = exp(two*pi*zi/three) ! exp(2 * i * Pi / 3)
  phm3 = exp(-two*pi*zi/three) ! exp(-2 * i * Pi / 3)
  ph6 = exp(pi*zi/three)    ! exp(i * Pi / 3)
  phm6 = exp(-pi*zi/three)  ! exp(-i * Pi / 3)
  indx = 0
  nab = 0
  !
  ! selected allowed transformations
  itr = 0
  ntr = 0
  do ii = 1, ntrans
     ! ignore operations with non-zero translation vector
     length = sqrt( dot_product(tnp(ii,:),tnp(ii,:)) )
     if (length > tol) cycle
     ! ignore operations where rotation matrix is not diagonal
     tmtrx = rmtrx(ii,:,:)
     do jj = 1, 3
        tmtrx(jj,jj) = 0
     enddo
     if (maxval(abs(tmtrx)) /= 0) cycle
     if (norot > 0) then
        if (trim(id(ii)) /= 'E' .and. trim(id(ii)) /= 'IE') cycle
     endif
     ntr = ntr + 1
     itr(ntr) = ii
  enddo

  do it = 1, ntr
     ii = itr(it)
     ! locate identity
     if (trim(id(ii)) == 'E') then
        nab = nab + 1
        i_e = ii
     endif
     ! is there inversion?
     if (trim(id(ii)) == 'IE') then
        l_i = .true.
        nab = nab + 1
        i_i = ii
     endif
     ! is there C2 rotation?
     if (.not. l_c2 .and. trim(id(ii)) == 'C2') then
        l_c2 = .true.
        nab = nab + 1
        i_c2 = ii
     endif
     ! is there C3 rotation?
     if (.not. l_c3 .and. trim(id(ii)) == 'C3') then
        l_c3 = .true.
        nab = nab + 1
        i_c3 = ii
     endif
     ! is there C6 rotation?
     if (.not. l_c6 .and. trim(id(ii)) == 'C6') then
        l_c6 = .true.
        nab = nab + 1
        i_c6 = ii
     endif
     ! is there IC2 mirror?
     if (.not. l_s2 .and. trim(id(ii)) == 'IC2') then
        l_s2 = .true.
        nab = nab + 1
        i_s2 = ii
     endif
  enddo

  ! is there a second C2 that commutes with the first C2?
  c2_loop: do it = 1, ntr
     ii = itr(it)
     if (trim(id(ii)) /= 'C2') cycle
     do jt = 1, ntr
        jj = itr(jt)
        if (trim(id(jj)) /= 'C2') cycle
        if (mult(ii,jj) == mult(jj,ii) .and. ii /= jj) then
           l_c22 = .true.
           i_c2 = ii
           i_c22 = jj
           nab = nab + 1
           exit c2_loop
        endif
     enddo
  enddo c2_loop
  !
  ! is there IC2 mirror? (if there is C2, it must commute with IC2)
  l_s2 = .false.
  do it = 1, ntr
     if (trim(id(itr(it))) == 'IC2') l_s2 = .true.
  enddo
  if (l_s2 .and. l_c2) then
     l_s2 = .false.
     nab = nab - 1
  endif
  mirror_loop: do it = 1, ntr
     ii = itr(it)
     if (trim(id(ii)) /= 'C2') cycle
     do jt = 1, ntr
        jj = itr(jt)
        length = sqrt( dot_product(tnp(jj,:),tnp(jj,:)) )
        if (length > tol) cycle
        if (trim(id(jj)) /= 'IC2') cycle
        if (mult(ii,jj) == mult(jj,ii)) then
           l_s2 = .true.
           i_c2 = ii
           i_s2 = jj
           nab = nab + 1
           exit mirror_loop
        endif
     enddo
  enddo mirror_loop
  !
  ! is there IC2 mirror? (if there is C3, it must commute with IC2)
  l_s3 = .false.
  do it = 1, ntr
     if (trim(id(itr(it))) == 'IC2') l_s3 = .true.
  enddo
  if (l_s3 .and. l_c3) then
     l_s3 = .false.
     nab = nab - 1
  endif
  mirror2_loop: do it = 1, ntr
     ii = itr(it)
     if (trim(id(ii)) /= 'C3') cycle
     do jt = 1, ntr
        jj = itr(jt)
        length = sqrt( dot_product(tnp(jj,:),tnp(jj,:)) )
        if (length > tol) cycle
        if (trim(id(jj)) /= 'IC2') cycle
        if (mult(ii,jj) == mult(jj,ii)) then
           l_s3 = .true.
           i_c3 = ii
           i_s3 = jj
           nab = nab + 1
           exit mirror2_loop
        endif
     enddo
  enddo mirror2_loop
  !
  ! Determine the point group based on available generators:
  !
  if (l_i) then
     if (l_c6) then
        name = 'C_6h'
     elseif (l_c2) then
        if (l_c22) then
           name = 'D_2h'
        else
           name = 'C_2h'
        endif
     else
        name = 'C_i'
     endif
  else
     if (l_c3) then
        if (l_s3) then
           name = 'C_3h'
        else
           name = 'C_3'
        endif
     elseif (l_c2) then
        if (l_c22) then
           name = 'D_2'
        elseif (l_s2) then
           name = 'C_2v'
        else
           name = 'C_2'
        endif
     else
        if (l_s2) then
           name = 'C_s'
        else
           name = 'C_1'
        endif
     endif
  endif
  !
  ! construct character table and index
  !
  select case(trim(name))
  case('C_6h')
     nab = 12
     indx(1) = i_e          ! E
     indx(2) = i_c6         ! C6
     indx(3) = mult(i_c6,i_c6) ! C3 = C6 * C6
     indx(4) = mult(i_c6,indx(3)) ! C2 = C6 * C6 * C6
     indx(5) = mult(i_c6,indx(4)) ! C3' = C6 * C6 * C6 * C6
     indx(6) = mult(i_c6,indx(5)) ! C6' = C6 * C6 * C6 * C6 * C6
     indx(7) = i_i          ! IE
     indx(8) = mult(indx(2),i_i) ! IC6 = C6 * IE
     indx(9) = mult(indx(3),i_i) ! IC3 = C3 * IE
     indx(10) = mult(indx(4),i_i) ! IC2 = C2 * IE
     indx(11) = mult(indx(5),i_i) ! IC3' = C3' * IE
     indx(12) = mult(indx(6),i_i) ! IC6' = C6' * IE
     chi( 1,1:12) = (/ ph1, ph1, ph1, ph1, ph1, ph1, &
          ph1, ph1, ph1, ph1, ph1, ph1 /)
     chi( 2,1:12) = (/ ph1,-ph1, ph1,-ph1, ph1,-ph1, &
          ph1,-ph1, ph1,-ph1, ph1,-ph1 /)
     chi( 3,1:12) = (/ ph1, ph6,-phm6,-ph1,-ph6, phm6, &
          ph1, ph6,-phm6,-ph1,-ph6, phm6 /)
     chi( 4,1:12) = conjg( chi(3,1:12) )
     chi( 5,1:12) = (/ ph1,-phm6,-ph6, ph1,-phm6,-ph6, &
          ph1,-phm6,-ph6, ph1,-phm6,-ph6 /)
     chi( 6,1:12) =  conjg( chi(5,1:12) )
     chi( 7,1:12) = (/ ph1, ph1, ph1, ph1, ph1, ph1, &
          -ph1,-ph1,-ph1,-ph1,-ph1,-ph1 /)
     chi( 8,1:12) = (/ ph1,-ph1, ph1,-ph1, ph1,-ph1, &
          -ph1, ph1,-ph1, ph1,-ph1, ph1 /)
     chi( 9,1:12) = (/ ph1, ph6,-phm6,-ph1,-ph6, phm6, &
          -ph1,-ph6, phm6, ph1, ph6,-phm6 /)
     chi(10,1:12) = conjg( chi(9,1:12) )
     chi(11,1:12) = (/ ph1,-phm6,-ph6, ph1,-phm6,-ph6, &
          -ph1, phm6, ph6,-ph1, phm6, ph6 /)
     chi(12,1:12) = conjg( chi(10,1:12) )
     label(1) = 'Au '
     label(2) = 'Bg '
     label(3) = 'E1g'
     label(4) = 'E1g'
     label(5) = 'E2g'
     label(6) = 'E2g'
     label(7) = 'Au '
     label(8) = 'Bu '
     label(9) = 'E1u'
     label(10) = 'E1u'
     label(11) = 'E2u'
     label(12) = 'E2u'
  case('C_3h')
     nab = 6
     indx(1) = i_e          ! E
     indx(2) = i_c3         ! C3
     indx(3) = mult(i_c3,i_c3) ! C3' = C3 * C3
     indx(4) = i_s3        ! IC3
     indx(5) = mult(indx(2),i_s3) ! IC3 = C3 * IE
     indx(6) = mult(indx(3),i_s3) ! IC3' = C3' * IE
     chi(1,1:6) = (/ ph1, ph1, ph1, ph1, ph1, ph1 /)
     chi(2,1:6) = (/ ph1, ph3, phm3, ph1, ph3, phm3 /)
     chi(3,1:6) = conjg( chi(2,1:6) )
     chi(4,1:6) = (/ ph1, ph1, ph1,-ph1,-ph1,-ph1 /)
     chi(5,1:6) = (/ ph1, ph3, phm3,-ph1,-ph3,-phm3 /)
     chi(6,1:6) = conjg( chi(5,1:6) )
     label(1) = 'Ap '
     label(2) = 'Ep '
     label(3) = 'Ep '
     label(4) = 'App'
     label(5) = 'Epp'
     label(6) = 'Epp'
  case('C_3')
     nab = 3
     indx(1) = i_e          ! E
     indx(2) = i_c3         ! C3
     indx(3) = mult(i_c3,i_c3) ! C3' = C3 * C3
     chi(1,1:3) = (/ ph1, ph1, ph1 /)
     chi(2,1:3) = (/ ph1, ph3, phm3 /)
     chi(3,1:3) = conjg( chi(3,1:3) )
     label(1) = 'A  '
     label(2) = 'E  '
     label(3) = 'E  '
  case('D_2h')
     nab = 8
     indx(1) = i_e          ! E
     indx(2) = i_c2         ! C2_x
     indx(3) = i_c22        ! C2_y
     indx(4) = mult(i_c22, i_c2) ! C2_z = C2_y * C2_x
     indx(5) = i_i          ! IE
     indx(6) = mult(indx(2), indx(5)) ! IC2_x
     indx(7) = mult(indx(3), indx(5)) ! IC2_y
     indx(8) = mult(indx(4), indx(5)) ! IC2_z
     chi(1,1:8) = (/ ph1, ph1, ph1, ph1, ph1, ph1, ph1, ph1/)
     chi(2,1:8) = (/ ph1, ph1,-ph1,-ph1, ph1, ph1,-ph1,-ph1/)
     chi(3,1:8) = (/ ph1,-ph1, ph1,-ph1, ph1,-ph1, ph1,-ph1/)
     chi(4,1:8) = (/ ph1,-ph1,-ph1, ph1, ph1,-ph1,-ph1, ph1/)
     chi(5,1:8) = (/ ph1, ph1, ph1, ph1,-ph1,-ph1,-ph1,-ph1/)
     chi(6,1:8) = (/ ph1, ph1,-ph1,-ph1,-ph1,-ph1, ph1, ph1/)
     chi(7,1:8) = (/ ph1,-ph1, ph1,-ph1,-ph1, ph1,-ph1, ph1/)
     chi(8,1:8) = (/ ph1,-ph1,-ph1, ph1,-ph1, ph1, ph1,-ph1/)
     label(1) = 'Ag '
     label(2) = 'B1g'
     label(3) = 'B2g'
     label(4) = 'B3g'
     label(5) = 'Au '
     label(6) = 'B1u'
     label(7) = 'B2u'
     label(8) = 'B3u'
  case('C_2h')
     nab = 4
     indx(1) = i_e          ! E
     indx(2) = i_c2         ! C2
     indx(3) = i_i          ! IE
     indx(4) = mult(i_c2, i_i) ! IC2 = C2 * IE
     chi(1,1:4) = (/ ph1, ph1, ph1, ph1/)
     chi(2,1:4) = (/ ph1,-ph1, ph1,-ph1/)
     chi(3,1:4) = (/ ph1, ph1,-ph1,-ph1/)
     chi(4,1:4) = (/ ph1,-ph1,-ph1, ph1/)
     label(1) = 'Ag '
     label(2) = 'Bg '
     label(3) = 'Au '
     label(4) = 'Bu '
  case('D_2')
     nab = 4
     indx(1) = i_e          ! E
     indx(2) = i_c2         ! C2_x
     indx(3) = i_c22        ! C2_y
     indx(4) = mult(i_c22, i_c2) ! C2_z = C2_y * C2_x
     chi(1,1:4) = (/ ph1, ph1, ph1, ph1/)
     chi(2,1:4) = (/ ph1, ph1,-ph1,-ph1/)
     chi(3,1:4) = (/ ph1,-ph1, ph1,-ph1/)
     chi(4,1:4) = (/ ph1,-ph1,-ph1, ph1/)
     label(1) = 'A  '
     label(2) = 'B1 '
     label(3) = 'B2 '
     label(4) = 'B3 '
  case('C_2v')
     nab = 4
     indx(1) = i_e          ! E
     indx(2) = i_c2         ! C2
     indx(3) = i_s2        ! IC2
     indx(4) = mult(i_c2, i_s2) ! IC2' = C2 * IC2
     chi(1,1:4) = (/ ph1, ph1, ph1, ph1/)
     chi(2,1:4) = (/ ph1, ph1,-ph1,-ph1/)
     chi(3,1:4) = (/ ph1,-ph1, ph1,-ph1/)
     chi(4,1:4) = (/ ph1,-ph1,-ph1, ph1/)
     label(1) = 'A1 '
     label(2) = 'A2 '
     label(3) = 'B1 '
     label(4) = 'B2 '
  case('C_2')
     nab = 2
     indx(1) = i_e          ! E
     indx(2) = i_c2         ! C2
     chi(1,1:2) = (/ ph1, ph1/)
     chi(2,1:2) = (/ ph1,-ph1/)
     label(1) = 'A  '
     label(2) = 'B  '
  case('C_s')
     nab = 2
     indx(1) = i_e          ! E
     indx(2) = i_s2        ! IC2
     chi(1,1:2) = (/ ph1, ph1/)
     chi(2,1:2) = (/ ph1,-ph1/)
     label(1) = 'A1 '
     label(2) = 'A2 '
  case('C_i')
     nab = 2
     indx(1) = i_e          ! E
     indx(2) = i_i          ! IE
     chi(1,1:2) = (/ ph1, ph1/)
     chi(2,1:2) = (/ ph1,-ph1/)
     label(1) = 'Ag '
     label(2) = 'Au '
  case('C_1')
     nab = 1
     indx(1) = i_e          ! E
     chi(1,1) =  ph1
     label(1) = 'A  '
  end select

end subroutine abelian
!===============================================================
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
