!===============================================================
!
! Copyright (C) 2005 Finite Difference Research Group
! This file is part of parsec, http://www.ices.utexas.edu/parsec/
!
! This is a driver routine for the symmetry operation routine
! symgen. It determines the set of existing symmetry operations
! based on crystallographic systems. For non-periodic systems
! (pbc%is_on = .false.), we need a periodic cell in order to find
! the crystallographic system. Since a regular grid is always
! going to impose cubic symmetry, even for highly symmetric
! systems like an isolated atom, choose a cubic supercell with
! side four times the boudary radius.
!
! Unit lattice vectors:
! A = 3x3 matrix that stores unit lattice vectors in column-wise
! format
!
! transformation of coordinates from cartesian to lattice
! coordinates is done in the usual way:
!    ( r )_cart =     A  * ( r )_latt
!    ( r )_latt = inv(A) * ( r )_cart
!
! A generic symmetry operation T is defined this way (see Jones,
! 'The Theory of Brillouin Zones and Electronic States in
! Crystals'):
!
!   T ( r ) -> r' = M * ( r - t)
!
! where:
!   r is original position (vector written in units of lattice vectors)
!   r' is a point related to r by operation T (same units)
!   M is a 3x3 matrix that executes rotation
!   t is a translation vector
!   operation T itself is defined by (M,t)
!
! Notice that M is not necessarily orthogonal 
! ( M * transpose(M) is not I ), but the matrix 
! U = A * M * inv(A) IS orthogonal!
!
! In structure symm:
!   M == rmtrx
!   inv(transpose(M)) = gmtrx (needed for rotations in Fourier space)
!   t == tnp
!   U == trans (rotation matrix in cartesian components)
!   A == alatt
!   transpose(inv(A)) == invlat
!
! Author: Murilo Tiago (2005)
!
!---------------------------------------------------------------
subroutine symmetries(clust,grid,pbc,symm,rsymm,ipr,info)

  use constants
  use cluster_module 
  use grid_module
  use pbc_module
  use symmetry_module

  implicit none
  !
  ! Input/Output variables:
  !
  ! the cluster
  type (cluster), intent(in) :: clust
  ! grid related data
  type (grid_data), intent(in) :: grid
  ! periodic boundary conditions data
  type (pbc_data), intent(in) :: pbc
  ! symmetry operations
  type (symmetry), intent(inout) :: symm,rsymm
  ! output flag
  integer, intent(in) :: ipr
  ! on input:
  !     info = 0 : build full symmetry group and Abelian subgroup
  !                for the system;
  !     info = 1 : build full symmetry group for the system but do
  !                not build Abelian subgroup; this will make sure
  !                that no symmetry operations are included in the
  !                real-space grid or in the Hamiltonian;
  !     info = 2 : do not build any symmetry group, and
  !                ignore symmetry operations throughout the
  !                calculation.
  ! on output:
  !     info = 0 : successfull exit;
  !     info /= 0 : failure.
  integer, intent(inout) :: info
  !
  ! Work variables:
  !
  ! flags passed to symgen
  integer invers_no, nops, norot
  ! variables to be saved to symm structure at the end
  ! (notice the reshape of arrays!)
  integer ntrans
  integer gmtrx(48, 3, 3)
  real(dp) :: tnp(48, 3), trans(3, 3, 48)
  integer rmtrx(48, 3, 3)
  ! number of elements in maximal Abelian subgroup
  integer nab
  ! index for the maximal Abelian subgroup and character table
  integer :: indx(12)
  complex(dpc) :: chi(12,12)
  ! counters, temporary variables
  integer mxdatm, ii, jj, ity
  real(dp) :: tmpvec(3)
  real(dp), allocatable :: coorat(:,:,:)

  !---------------------------------------------------------------
  !
  ! alatt: matrix A. If no PBC are used, construct an
  ! artificial cubic supercell that encloses the cluster and
  ! with unit vectors along cartesian axes x,y,z
  !
  symm%alatt(:,:) = zero
  if (pbc%is_on) then
     symm%alatt = pbc%latt_vec
     do ii = pbc%per + 1, 3
        symm%alatt(ii,ii) = four*grid%rmax
     enddo
  else
     symm%alatt(1,1) = four*grid%rmax
     symm%alatt(2,2) = four*grid%rmax
     symm%alatt(3,3) = four*grid%rmax
  endif

  symm%invlat = transpose(symm%alatt)
  call mtrxin(symm%invlat,tmpvec(1),tmpvec(2))

  mxdatm = maxval(clust%natmi)
  !
  ! coorat: coordinates of all atoms with respect to unit lattice
  ! vectors, in a layout proper for symgen
  !
  allocate(coorat(clust%type_num,mxdatm,3))
  coorat(:,:,:) = zero
  ii = 0
  do ity = 1, clust%type_num
     do jj = 1,clust%natmi(ity)
        ii = ii + 1
        tmpvec(1) = clust%xatm(ii)
        tmpvec(2) = clust%yatm(ii)
        tmpvec(3) = clust%zatm(ii)
        call matvec3('T',symm%invlat,tmpvec,tmpvec)
        coorat(ity,jj,:) = tmpvec
     enddo
  enddo

  if (info /= 2) then
     symm%use_symm = .true.
     !
     ! If symmetry operations are to be used, search for symmetry
     ! operations.
     !
     !     nops = -1  : include fractional translations
     !     nops = -2  : do not include fractional translations
     nops = -1
     !
     ! If lattice vectors are not orthogonal (pbc%latt_vec not diagonal),
     ! the Abelian group must not have rotation/mirror matrices.
     !
     norot = 0
     trans(:,:,1) = symm%alatt
     do ii = 1, 3
        trans(ii,ii,1) = zero
     enddo
     if (maxval(abs(trans(:,:,1))) > 1.d-8) norot = 1

     if (info == 1) norot = -1

     call symgen(ipr, symm%alatt, coorat, clust%natmi, &
          clust%type_num, invers_no, nops, norot, &
          ntrans, gmtrx, tnp, mxdatm, rmtrx, trans, nab, indx, chi, &
          symm%cell_symmetry)
  else
     symm%use_symm = .false.
     !
     ! Otherwise, create the structure but allow only Identity as
     ! symmetry operation.
     !
     ntrans = 1
     gmtrx = 0
     rmtrx = 0
     trans = 0
     tnp(:,:) = zero
     nab = 1
     indx(1) = 1
     chi(1,1) = zone
     do jj=1,3
        gmtrx(1,jj,jj) = 1
        rmtrx(1,jj,jj) = 1
        trans(jj,jj,1) = 1
     enddo
     write(7,'(/,a)') 'skipping non-trivial symmetry operations'
  endif
  !
  ! Save information about symmetry operations:
  ! define trans = transpose(symm%trans)
  !
  ! Full point group stored in structure symm.
  call create_symmetry(ntrans,symm)
  do jj = 1, ntrans
     symm%gmtrx(:,:,jj) = gmtrx(jj,:,:)
     symm%rmtrx(:,:,jj) = rmtrx(jj,:,:)
     symm%trans(:,:,jj) = transpose(trans(:,:,jj))
     symm%tnp(:,jj) = tnp(jj,:)
  enddo
  !
  ! Sanity check.
  !
  if (ntrans == 0) then
     info = 1
     return
  else
     info = 0
  endif
  !
  ! Abelian subgroup (reduced group) stored in rsymm.
  !
  call create_symmetry(nab,rsymm)
  rsymm%use_symm = symm%use_symm
  rsymm%alatt = symm%alatt
  rsymm%invlat = symm%invlat
  do ii = 1, nab
     jj = indx(ii)
     rsymm%gmtrx(:,:,ii) = gmtrx(jj,:,:)
     rsymm%rmtrx(:,:,ii) = rmtrx(jj,:,:)
     rsymm%trans(:,:,ii) = trans(:,:,jj)
     rsymm%tnp(:,ii) = tnp(jj,:)
     rsymm%chi(1:nab,ii) = anint(real(chi(1:nab,ii),dp))
  enddo
  deallocate(coorat)

  if (nab == 1) then
     write(7,'(a,/)') 'Abelian subgroup contains only the Identity operation'
  endif

end subroutine symmetries
!===============================================================
!
! For a given vector r in cartesian coordinates, apply it-th
! symmetry operation, T, according to the usual rule:
!              T ( r ) -> tr = M * ( r - t)
! Output is vector tr in cartesian coordinates.
!
!---------------------------------------------------------------
subroutine symop(symm,it,rin,tr)

  use constants
  use symmetry_module
  implicit none
  !
  ! Input/Output variables:
  !
  ! symmetry operations in reduced group:
  type (symmetry), intent(in) :: symm
  ! index of this operation
  integer, intent(in) :: it
  ! input vector
  real(dp), intent(in) :: rin(3)
  ! transformed vector, tr = T ( r )
  real(dp), intent(out) :: tr(3)
  !
  ! Work variables:
  !
  ! translation vectors in cartesian units (needed because the symm%tnp
  ! are defined in units of supercell unit vectors)
  real(dp) :: tnplat(3)
  !---------------------------------------------------------------
  call matvec3('N',symm%alatt,symm%tnp(:,it),tnplat)
  tnplat = rin - tnplat
  call matvec3('N',symm%trans(:,1,it),tnplat,tr)

end subroutine symop
!===============================================================
!
! For a given vector r in coordinates relative to lattice vectors, apply it-th
! symmetry operation, T, according to the usual rule:
!              T ( r ) -> tr = M * ( r - t)
! Output is vector tr in cartesian coordinates.
!
! WARNING: input vector r is really given in lattice vector
! coordinates! For example, vector (1,0,0) is identical to the first
! lattice vector, with no normalization.
!
!---------------------------------------------------------------
subroutine symop_relative(symm,it,rin,tr)

  use constants
  use symmetry_module
  implicit none
  !
  ! Input/Output variables:
  !
  ! symmetry operations in reduced group:
  type (symmetry), intent(in) :: symm
  ! index of this operation
  integer, intent(in) :: it
  ! input vector
  real(dp), intent(in) :: rin(3)
  ! transformed vector, tr = T ( r )
  real(dp), intent(out) :: tr(3)
  !
  ! Work variables:
  !
  ! temporary vectors and matrices
  real(dp) :: tnplat(3), rmtrx(3,3)
  !---------------------------------------------------------------
  tnplat = rin - symm%tnp(:,it)
  rmtrx = symm%rmtrx(:,:,it)
  call matvec3('N',rmtrx,tnplat,tr)

end subroutine symop_relative
!===============================================================
