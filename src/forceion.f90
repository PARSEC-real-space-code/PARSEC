!===============================================================
!
! Copyright (C) 2005 Finite Difference Research Group
! This file is part of parsec, http://www.ices.utexas.edu/parsec/
!
! This subroutine calculates ionic force contributions and the
! total nuclear energy using classical formulas, based on the
! ionic charges and the distance between the atoms, as well as
! point charges if present.
!  
!---------------------------------------------------------------
subroutine forceion(clust,zion,ipr,enuc)

  use constants
  use cluster_module
  implicit none
  !
  ! Input/Output variables:
  !
  ! the cluster
  type (cluster), intent(inout) :: clust
  ! positive charge of pseudo ion for each ion type
  real(dp), intent(in) :: zion(clust%type_num)
  ! print flag
  integer, intent(in) :: ipr
  ! total nuclear energy
  real(dp), intent(out) :: enuc
  !
  ! Work variables:
  !
  ! distance between two atoms and their cube,
  ! multiplication of their charges
  real(dp) :: rij, rij2, qiqj, dxx, dyy, dzz,qdevr, enucpt
  ! counters
  integer ia,ja,ity,jty,i,j

  !---------------------------------------------------------------

  ! initialize nuclear energy
  enuc = zero

  ! initialize force on each atom
  clust%force(:,:) = zero

  ! initialize first atom counter ("i atoms")
  ia = 0
  ! for all "i" atoms (atom types and atom within type)
  do ity = 1,clust%type_num
     do i = 1, clust%natmi(ity)
        ia = ia + 1
        ! initialize second atom counter ("j atoms")
        ja = 0
        ! for all "j" atoms (atom types and atom within type)
        do jty = 1, clust%type_num
           do j = 1, clust%natmi(jty)
              ja = ja + 1
              ! for each i, j pair (disregard ia=ja because this is unphysical 
              ! self-interaction):
              ! calculate distance between atoms and its cube, and the
              ! multiplication of their ionic charges
              if (ia /= ja) then
                 dxx = clust%xatm(ia)-clust%xatm(ja)
                 dyy = clust%yatm(ia)-clust%yatm(ja)
                 dzz = clust%zatm(ia)-clust%zatm(ja)
                 rij2 = (dxx*dxx+dyy*dyy+dzz*dzz)
                 rij = sqrt(rij2)
                 qiqj = zion(ity)*zion(jty)
                 qdevr = qiqj/(rij2*rij)
                 ! update forces: F=qi*qj*r/|r|^3
                 ! the factor of two is because we want the force in Ryd/bohr
                 ! and not Hartree/bohr !
                 clust%force(1,ia) = clust%force(1,ia) + &
                      two*(clust%xatm(ia)-clust%xatm(ja))*qdevr
                 clust%force(2,ia) = clust%force(2,ia) + &
                      two*(clust%yatm(ia)-clust%yatm(ja))*qdevr
                 clust%force(3,ia) = clust%force(3,ia) + &
                      two*(clust%zatm(ia)-clust%zatm(ja))*qdevr

                 ! Compute total nuclear energy.
                 ! Here there is no unit problem because the two of
                 ! Hartree->Rydberg cancels with the half resulting from
                 ! summing over all atom pairs, so each pair is counted twice.
                 enuc = enuc + qiqj/rij
              endif
           enddo
        enddo
     enddo
  enddo

  ! Calculate forces and energy from N-Q and Q-Q
  ! when point charges Q are present.

  if (clust%has_ptchrg) then
     call forcept(clust,zion,ipr,enucpt)
     enuc = enuc + enucpt
  endif

  ! if print flag on, print the ion-ion forces
  if(ipr >= 1) then
     write(7,*)
     write(7,*) 'Forces from ion-ion interaction:'
     write(7,*) '================================'
     write(7,22)
     write(7,*)
     do ia = 1, clust%atom_num
        write(7,20) ia,clust%xatm(ia),clust%yatm(ia) &
             ,clust%zatm(ia),clust%force(:,ia)
     enddo
     write(7,*)
  endif

20 format(i4,2x,3(f11.6,1x),2x,3(f11.6,1x))
22 format('-atom-  ----x----   ----y----   ----z-----' &
        ,'   ----Fx----  ----Fy----   ----Fz---',/ &
        '                       [bohr]           ' &
        ,'                  [Ry/bohr]            ')

end subroutine forceion
!===============================================================
