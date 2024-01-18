!===============================================================
!
! Copyright (C) 2005 Finite Difference Research Group
! This file is part of parsec, http://www.ices.utexas.edu/parsec/
!
! Computes the local potential contribution to the Hellmann-
! Feynman forces. It uses the previous charge density (vsc_chg)
! and potential (vsc). Calculations are done in reciprocal space.
!
! Adapted from plane-wave programs written by S. Froyen and
! J. L. Martins.
!
!---------------------------------------------------------------
subroutine force_loc_pb(clust,pbc,vsc,vsc_chg)

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
  ! charge density and potential
  complex(dpc), dimension(pbc%nstar), intent(in) :: vsc, vsc_chg
  !
  ! Work variables:
  !
  integer i,j,k,l,natomi,kadd,jj,mstarj,jaa,ntype
  integer natom(clust%type_num)
  real(dp) :: fi,exp1,vdgr,vdgi
  real(dp) :: floc(3,clust%atom_num,clust%type_num),ftmp(3)
  complex(dpc) :: vds

  !---------------------------------------------------------------
  natom = clust%natmi     
  ntype = clust%type_num

  floc(:,:,:) = zero

  ! Loop over atomic types.
  do i=1,ntype
     ! Loop over stars.
     kadd = pbc%ng + 1
     do jj=2,pbc%nstar
        j = pbc%nstar - jj + 2
        vds =  conjg( pbc%vql(i,j) * vsc_chg(j) + vsc(j) * pbc%dnc(i,j))
        ! Loop over g vectors in star.
        mstarj = pbc%mstar(j)
        do k=1,mstarj
           kadd = kadd - 1
           vdgr = real(vds) * real(pbc%phase(kadd),dp) - &
                aimag(vds) * aimag(pbc%phase(kadd))
           vdgi = pbc%conj(kadd) * (aimag(vds)* real(pbc%phase(kadd),dp) + &
                real(vds,dp) * aimag(pbc%phase(kadd)))
           ! Loop over atoms of same type.
           natomi = natom(i)
           do l=1,natomi
              fi = dot_product(real(pbc%kgv(:,kadd),dp) ,pbc%rat(:,l,i))
              exp1 = sin(fi)*vdgr - cos(fi)*vdgi
              ! Add to forces.
              floc(:,l,i) = floc(:,l,i) + real(pbc%kgv(:,kadd),dp) * exp1
           enddo
        enddo
     enddo
  enddo

  jaa = 0
  do i = 1, ntype
     natomi = natom(i)
     do j=1,natomi
        jaa = jaa + 1
        call matvec3('N',pbc%bvec,floc(1,j,i),ftmp)
        clust%force(:,jaa) = clust%force(:,jaa) + two * ftmp
     enddo
  enddo

end subroutine force_loc_pb
!===============================================================
