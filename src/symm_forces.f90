!===============================================================
!
! Copyright (C) 2005 Finite Difference Research Group
! This file is part of parsec, http://www.ices.utexas.edu/parsec/
!
! Symmetrize the set of atomic forces according to symmetry
! operations.
!
! Author: Murilo Tiago (2005)
!
!---------------------------------------------------------------
subroutine symm_forces(clust,symm,ierr)

  use constants
  use cluster_module
  use symmetry_module

  implicit none
  !
  ! Input/Output variables:
  !
  ! the cluster
  type (cluster), intent(inout) :: clust
  ! symmetry operations
  type (symmetry), intent(in) :: symm
  ! error flag, 410 < ierr < 421
  integer, intent(out) :: ierr
  !
  ! Work variables:
  !
  ! symmetrized forces
  real(dp) :: fsym(3,clust%atom_num)
  ! input forces on ions
  real(dp), allocatable :: cforce(:,:,:)
  ! temporary variables
  integer :: ii,jj,kk,iat,ity,itrans,mxdatm,ishift
  real(dp) :: cro(3),tmpvec(3),dist
  ! atomic coordinates
  real(dp), allocatable :: ratom(:,:,:)
  ! tolerance in the position of one atom and its equivalent by
  ! symmetry
  real(dp), parameter :: tol = 1.d-6

  !---------------------------------------------------------------

  fsym(:,:) = zero
  mxdatm = maxval(clust%natmi)
  allocate(ratom(clust%type_num,mxdatm,3))
  allocate(cforce(clust%type_num,mxdatm,3))
  ratom(:,:,:) = zero
  iat = 0
  do ity = 1, clust%type_num
     do jj = 1, clust%natmi(ity)
        iat = iat + 1
        ratom(ity,jj,1) = clust%xatm(iat)
        ratom(ity,jj,2) = clust%yatm(iat)
        ratom(ity,jj,3) = clust%zatm(iat)
        cforce(ity,jj,:) = clust%force(:,iat)
     enddo
  enddo

  iat = 0
  do ity = 1, clust%type_num
     do jj = 1, clust%natmi(ity)
        iat = iat + 1
        do itrans = 1, symm%ntrans
           !
           !                -1
           ! find cro = mtrx    * (ratom - tnp)
           !
           tmpvec = ratom(ity,jj,:)
           call symop(symm,itrans,tmpvec,cro)
           do ii = 1, clust%natmi(ity)
              ! make sure the rotated atom is inside the periodic cell
              tmpvec = cro - ratom(ity,ii,:)
              call matvec3('T',symm%invlat,tmpvec,tmpvec)
              do kk = 1, 3
                 ishift = nint(tmpvec(kk))
                 tmpvec(kk) = tmpvec(kk) - one*ishift
              enddo
              dist = sqrt( dot_product(tmpvec,tmpvec) )
              if (dist > tol) cycle
              goto 10
           enddo
           write(7,*) ' ERROR: unable to find equivalent of atom'
           write(7,*) ' at position ',ratom(ity,jj,:)
           ierr = 411
           return
10         continue
           !
           ! rotate force and add
           !
           tmpvec = cforce(ity,ii,:)
           call matvec3('T',symm%trans(:,1,itrans),tmpvec,tmpvec)
           fsym(:,iat) = fsym(:,iat) + tmpvec(:)
        enddo
     enddo
  enddo
  !
  ! update forces
  !
  clust%force = fsym/(one*symm%ntrans)

  deallocate(ratom, cforce)

end subroutine symm_forces
!===============================================================
