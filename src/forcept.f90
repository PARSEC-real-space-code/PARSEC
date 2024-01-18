!===============================================================
!
! Copyright (C) 2005 Finite Difference Research Group
! This file is part of parsec, http://www.ices.utexas.edu/parsec/
!
! This subroutine calculates ionic force contributions and the
! total nuclear energy based on interactions (Q-Q,Q-N) between
! point charges Q and "real" nuclei N. Contributions are Q-Q
!
! authors: C. Buergel (Humboldt U) and O. Guliamov (Weizmann), 2006
!
!---------------------------------------------------------------
subroutine forcept(clust,zion,ipr,enucpt)
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
  real(dp), intent(out) :: enucpt
  !
  ! Work varipbles:
  !
  ! distance between Q,N and Q,Q and their cube,
  real(dp) :: rij, rij2, qiqj, dxx, dyy, dzz,qdevr
  ! counters
  integer ip,jp,ity,jty,i,j,jatty,ja

  !---------------------------------------------------------------

  ! Initialize nuclear energy.
  enucpt = zero

  ! Initialize first point charge counter ("i Q").
  ip = 0
  ! for all "i" Q (pc types and pc within type).
  do ity = 1,clust%npttyp
     do i = 1, clust%nptt(ity)
        ip = ip + 1
        ! Initialize second point charge counter ("j Q").
        jp = 0
        ! For all "j" Q (pc types and pc within type).
        do jty = 1, clust%npttyp
           do j = 1, clust%nptt(jty)
              jp = jp + 1
              ! For each i, j pair (disregard ip=jp because this is unphysical 
              ! self-interaction):
              ! Calculate distance between pcs and its cube, and the
              ! multiplication of their ionic charges.
              if (ip /= jp) then
                 dxx = clust%xpt(ip)-clust%xpt(jp)
                 dyy = clust%ypt(ip)-clust%ypt(jp)
                 dzz = clust%zpt(ip)-clust%zpt(jp)
                 rij2 = (dxx*dxx+dyy*dyy+dzz*dzz)
                 rij = sqrt(rij2)
                 qiqj = clust%qpt(ity)*clust%qpt(jty)

                 ! Compute total nuclear energy.
                 ! Here there is no unit problem because the two of
                 ! Hartree->Rydberg cancels with the half resulting from
                 ! summing over all pcs pairs, so each pair is counted twice.
                 ! Decomment if Q-Q should be added to energy.
                 !                 enucpt = enucpt + qiqj/rij
              endif
           enddo
        enddo

        ! Calculate forces and energy contribution due to Q-N interaction
        ! Initialize atom charge counter ("j" atoms).
        ja=0
        ! For all "j" atoms (atom types and within type).
        do jatty = 1, clust%type_num
           do j = 1, clust%natmi(jatty)
              ja = ja + 1
              ! Calculate distance between pc(i) and atom(j) and its cube,
              ! and the multiplication of their ionic charges.
              dxx = clust%xatm(ja)-clust%xpt(ip) 
              dyy = clust%yatm(ja)-clust%ypt(ip)
              dzz = clust%zatm(ja)-clust%zpt(ip)
              rij2 = (dxx*dxx+dyy*dyy+dzz*dzz)
              rij = sqrt(rij2)
              qiqj = clust%qpt(ity)*zion(jatty)
              qdevr = qiqj/(rij2*rij)

              ! Update forces: F=qi*qj*r/|r|^3.
              ! The factor of two is because we want the force in Rydberg/bohr
              ! and not Hartree/bohr !
              clust%force(1,ja) = clust%force(1,ja)+two*dxx*qdevr
              clust%force(2,ja) = clust%force(2,ja)+two*dyy*qdevr
              clust%force(3,ja) = clust%force(3,ja)+two*dzz*qdevr

              ! Compute total nuclear energy. Factor two because here no self
              ! cancellation and we want the energy in Rydberg/bohr.
              enucpt = enucpt + two*qiqj/rij
           enddo
        enddo
        ! End of Q-N part.
     enddo
  enddo

  ! If print flag on, print the ion-ion forces.
  if(ipr >= 1) then
     write(7,*)
     write(7,*) 'Forces from ion/point charge interaction:'
     write(7,*) '(combined forces)'
     write(7,*) '========================================='
     write(7,22)
     write(7,*)
     do ip = 1, clust%atom_num

        write(7,20) ip,clust%xatm(ip),clust%yatm(ip),clust%zatm(ip), &
             clust%force(:,ip)
     enddo
     write(7,*)
  endif

20 format(i4,2x,3(f11.6,1x),2x,3(f11.6,1x))
22 format('-atom-  ----x----   ----y----   ----z-----', &
        '   ----Fx----  ----Fy----   ----Fz---',/ &
        '                       [bohr]           ', &
        '                  [Ry/bohr]            ')

end subroutine forcept
!===============================================================
