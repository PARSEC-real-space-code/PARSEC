!===============================================================
!
! Copyright (C) 2005 Finite Difference Research Group
! This file is part of parsec, http://www.ices.utexas.edu/parsec/
!
! Builds the irreducible wedge of the Brillouin zone given the
! specification of the Monkhorst-Pack grid (order and shift).
!
! Author: Murilo Tiago (2006)
!
!---------------------------------------------------------------
subroutine irrbz(elec_st,symm,pbc,lpr,ierr)

  use constants
  use electronic_struct_module
  use symmetry_module
  use pbc_module
  implicit none
  !
  ! Input/Output variables:
  !
  ! electronic structure
  type (electronic_struct), intent(inout) :: elec_st
  ! symmetry operations structure
  type (symmetry), intent(in) :: symm
  ! periodic boundary conditions data
  type (pbc_data), intent(in) :: pbc
  ! printing flag
  logical, intent(inout) :: lpr
  ! error flag, 240 < ierr < 251
  integer, intent(out) :: ierr
  !
  ! Work variables:
  !
  ! counters
  integer :: ii, jj, ll, itran, ikpt, nkpt, itest, inv, a
  ! temporary vectors
  real(dp), dimension(3) :: kprot, krot, ktest
  ! maximum number of k-points in the IBZ and in the full BZ
  integer :: nkptmax
  ! Monkhorst-Pack shift
  real(dp) :: mpshift(3)
  ! Monkhorst-Pack order
  integer :: mpgrid(3)
  ! temporary logical flag
  logical :: new_kpt
  ! temporary scalars
  real(dp) :: norm2, norm1, gmtrx(3,3)
  ! temporary arrays for k-points in IBZ and full BZ
  real(dp), dimension(:,:), allocatable :: kibz
  ! temporary arrays for index of k-points
  integer, dimension(:), allocatable :: kptrep
  ! temporary array for kinetic energy
  real(dp), dimension(:), allocatable :: ekin
  !
  ! tolerance in the construction of the IBZ
  real(dp), parameter :: tol_bz = 1.d-8

  !---------------------------------------------------------------
  !
  ! Allocate and initialize data.
  !
  nkptmax = 2 * symm%ntrans * elec_st%mpgrid(1) * elec_st%mpgrid(2) * &
       elec_st%mpgrid(3)
  allocate(kibz(3,nkptmax))
  allocate(kptrep(nkptmax))
  allocate(ekin(nkptmax))
  allocate(elec_st%kptgrid(elec_st%mpgrid(1),elec_st%mpgrid(2),elec_st%mpgrid(3)))
  mpshift = elec_st%mpshift
  mpgrid = elec_st%mpgrid

  ! Search for a set of k-points which are not mutually equivalent.
  ! This will be the irreducible wedge of the Brillouin zone (IBZ). The
  ! first point in this set is the first point in the Monkhorst-Pack grid.

  ikpt = 0
  kptrep = 0

  do ii = 1, elec_st%mpgrid(1)
     do jj = 1, elec_st%mpgrid(2)
        do ll = 1, elec_st%mpgrid(3)
           kprot(1) = (real(ii-1,dp) + mpshift(1))/real(mpgrid(1),dp)
           kprot(2) = (real(jj-1,dp) + mpshift(2))/real(mpgrid(2),dp)
           kprot(3) = (real(ll-1,dp) + mpshift(3))/real(mpgrid(3),dp)
           new_kpt = .true.
           do itran = 1, symm%ntrans
              gmtrx = symm%gmtrx(:,:,itran)
              call matvec3('N',gmtrx,kprot,krot)
              call translate01(krot)
              if (ikpt == 0) then
                 ikpt = ikpt + 1
                 kibz(:,ikpt) = krot
                 ktest = krot
                 ktest = ktest - anint(ktest)
                 call matvec3('N',pbc%bvec,ktest,ktest)
                 ekin(ikpt) = dot_product(ktest,ktest)
                 kptrep(1) = 0
                 elec_st%kptgrid(ii,jj,ll)=ikpt
              endif
              do itest = 1, ikpt
                 ! Remove time-reversal symmetry if noncollinear calculation
                 a = 0
                 ! Impose inversion (time-reversal) symmetry explicitly.
                 do inv = 1, -1, -2
                    a = a + 1
                    if (elec_st%ncl .and. a /= 1) cycle
                    ktest = krot - inv*kibz(:,itest)
                    call translate01(ktest)
                    call matvec3('N',pbc%bvec,ktest,ktest)
                    norm2 = dot_product(ktest,ktest)
                    !
                    ! Is ktest equivalent to any kpoint in IBZ?
                    ! If so, redefine the kpoint in IBZ as the shortest between
                    ! it and ktest.
                    ! If not, add ktest to the IBZ.
                    !
                    if (norm2 < tol_bz .and. new_kpt) then
                       new_kpt = .false.
                       ktest = kprot
                       ktest = ktest - anint(ktest)
                       call matvec3('N',pbc%bvec,ktest,ktest)
                       norm2 = dot_product(ktest,ktest) + eight
                       ktest = kibz(:,itest)
                       elec_st%kptgrid(ii,jj,ll)=itest
                       call matvec3('N',pbc%bvec,ktest,ktest)
                       norm1 = dot_product(ktest,ktest)
                       if (norm2 < norm1) then
                          ktest = kprot
                          call translate01(ktest)
                          kibz(:,itest) = ktest
                          elec_st%kptgrid(ii,jj,ll)=itest
                          ekin(itest) = norm2
                       endif
                       kptrep(itest) = kptrep(itest) + 1
                    endif
                 enddo
              enddo
           enddo
           if (new_kpt) then
              ikpt = ikpt + 1
              kibz(:,ikpt) = kprot
              ktest = krot
              ktest = ktest - anint(ktest)
              call matvec3('N',pbc%bvec,ktest,ktest)
              ekin(ikpt) = dot_product(ktest,ktest)
              kptrep(ikpt) = kptrep(ikpt) + 1
              elec_st%kptgrid(ii,jj,ll)=ikpt
           endif
        enddo
     enddo
  enddo

  elec_st%nkpt = ikpt

  nkpt = elec_st%mpgrid(1)*elec_st%mpgrid(2)*elec_st%mpgrid(3)
  if (sum(kptrep) /= nkpt) then
     if (lpr) then
        write(7,*) ' ERROR: incorrect BZ ! ',nkpt, sum(kptrep)
        write(7,*) kptrep
     endif
     ierr = 241
     return
  endif

  ! Done. Save data and clean up.
  ! K-points in elec_st%kpts and elec_st%kptf must be stored in Cartesian
  ! coordinates, units of inverse Bohr radius.

  allocate(elec_st%kpts(3,elec_st%nkpt))
  allocate(elec_st%kpwt(elec_st%nkpt))
  do ikpt = 1, elec_st%nkpt
     call matvec3('N',pbc%bvec,kibz(1,ikpt),elec_st%kpts(:,ikpt))
     elec_st%kpts(:,ikpt) = elec_st%kpts(:,ikpt)
     elec_st%kpwt(ikpt) = real(kptrep(ikpt),dp)/real(nkpt,dp)
  enddo

  if (lpr) then
     write(7,'(/,a)') ' K-point generation:'
     write(7,'(a)') ' -------------------'
     write(7,'(1x,i5,a,/,4x,a,3i5,6x,a,3f8.4,/)') elec_st%nkpt, &
          ' K points generated from Monkhorst-Pack grid :','order = ', &
          mpgrid,'shift = ',mpshift  
     write(7,'(4x,a,5x,a,8x,a,8x,a,5x,a,/)') 'I','weight', &
          'relative coord.','Cartesian coord.','E_kin (Ry)'
     do ikpt = 1, elec_st%nkpt
        ktest = kibz(:,ikpt)
        do itest = 1, 3
           ktest(itest) = ktest(itest) - anint(ktest(itest))
        enddo
        call matvec3('N',pbc%bvec,ktest,kprot)
        norm2 = dot_product(kprot,kprot)
        write(7,'(1x,i4,1x,f11.7,1x,3f8.4,3x,3f7.3,3x,f7.3)') &
             ikpt, elec_st%kpwt(ikpt),kibz(:,ikpt), &
             elec_st%kpts(:,ikpt),ekin(ikpt)
     enddo
  endif

  if (allocated(kibz)) deallocate(kibz)
  if (allocated(kptrep)) deallocate(kptrep)
  if (allocated(ekin)) deallocate(ekin)

  return
end subroutine irrbz
!===============================================================
!
! Translates a generic 3-D vector to the rectangular [0,1) cell.
!
! Input is vector kvec = (k1,k2,k3).
! Output is vector kvec = (q1,q2,q3) such that
!        0 <= q1 < 1  and  (q1-k1) is integer
!        similar prescriptions for (q2,k2) and (q3,k3).
!
!---------------------------------------------------------------
subroutine translate01(kvec)

  use constants
  implicit none
  !
  ! Input/Output variables:
  !
  real(dp), intent(inout) :: kvec(3)
  !
  ! Work variables:
  !
  integer :: ii
  real(dp) :: k_tmp
  ! tolerance in the construction of the IBZ
  real(dp), parameter :: tol_bz = 1.d-8
  !---------------------------------------------------------------
  do ii = 1, 3
     k_tmp = kvec(ii)
     k_tmp = k_tmp + half
     k_tmp = mod(k_tmp,one)   !translate to (-1,1)
     k_tmp = k_tmp - half
     if (abs(k_tmp) < tol_bz) k_tmp = zero
     if (k_tmp < zero) then
        k_tmp = k_tmp + one    !translate to [0,1)
     endif
     kvec(ii) = k_tmp
  enddo

  return
end subroutine translate01
!===============================================================
