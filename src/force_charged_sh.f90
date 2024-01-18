!===============================================================
!
! Copyright (C) 2015 Finite Difference Research Group
! This file is part of parsec, http://www.ices.utexas.edu/parsec/
!
! This subroutine is based on forcept.f90. It calculates changes
! forces and energy due to ion/charged sheet -- charged sheet
! (electrostatic) interaction.
!
!
! author: O. Sinai (Weizmann), 2015
!
!---------------------------------------------------------------
subroutine force_charged_sheet(clust,pbc,zion,ipr,enucsh)
  use constants
  use pbc_module
  use cluster_module
  implicit none
  !
  ! Input/Output variables:
  !
  ! the cluster
  type (cluster), intent(inout) :: clust
  !  periodic boundary conditions data
  type (pbc_data), intent(in) :: pbc
  ! positive charge of pseudo ion for each ion type
  real(dp), intent(in) :: zion(clust%type_num)
  ! print flag
  integer, intent(in) :: ipr
  ! Sheet-derived nuclear energy
  real(dp), intent(out) :: enucsh
  !
  ! Work variables:
  !
  ! distances and charges
  real(dp) :: dzabs, qjqi, dzz, qsarea
  ! counters
  integer ics,ity,i,j,jatty,ja,jcs,jty

  !---------------------------------------------------------------

  ! Initialize nuclear energy.
  enucsh = zero

  ! Initialize charged sheet counter ("i csheet").
  ics = 0
  ! for all "i" csheet (csheet types and csheet within type).
  do ity = 1, clust%number_charged_sheet_type
     do i = 1, clust%number_charged_sheets_per_type(ity)
        ics = ics + 1
        ! Calculate sheet charge per area
        qsarea = clust%sheet_charge_of_type(ity)/pbc%a_surface_cell

        ! Calculate forces and energy contribution due to 
        ! charged-sheet -- nuclei interaction
        ! Initialize atom charge counter ("j" atoms).
        ja=0
        ! For all "j" atoms (atom types and within type).
        do jatty = 1, clust%type_num
           do j = 1, clust%natmi(jatty)
              ja = ja + 1
              ! Calculate distance along z between csheet(i) and atom(j)
              dzz = clust%zatm(ja)-clust%z_charged_sheets(ics)
              dzabs = abs(dzz)
              ! Product of nuclear charge and sheet area charge
              qjqi = zion(jatty)*qsarea

              ! Update forces: Fz=qj(atom)*qi(sheet)*z/|z|.
              ! The factor of two is because we want the force in Rydberg/bohr
              ! and not Hartree/bohr!
              clust%force(3,ja) = clust%force(3,ja)+two*qjqi*dzz/dzabs
              
              ! Compute energy of ion in the potential of this sheet. Factor two
              ! because no self cancellation and energy should be in Rydberg,
              ! and a 2*pi factor comes from the potential expression.
              enucsh = enucsh - two*pi*two*qjqi*dzabs
           enddo
        enddo
        ! End charge-sheet -- nuclei part.

        ! Calculate energy contribution due to interaction between sheets
        ! Initialize additional charged sheet counter ("j csheet").
        jcs = 0
        ! for all "j" csheet (csheet types and csheet within type).
        do jty = 1, clust%number_charged_sheet_type
           do j = 1, clust%number_charged_sheets_per_type(jty)
              jcs = jcs + 1
              ! Ignore interaction of sheet with itself
              if (jcs /= ics) then
                 ! Distance along z between csheets (i) (j)
                 dzz = clust%z_charged_sheets(jcs)-clust%z_charged_sheets(ics)
                 dzabs = abs(dzz)
                 ! Product of sheet area charges. j'th sheet charge is not
                 ! divided by the surface cell area since this is cancelled by
                 ! multiplication by the same area to get energy per entire cell
                 qjqi = clust%sheet_charge_of_type(jty)*qsarea

                 ! Compute energy of j'th sheet in potential of i'th sheet. No factor
                 ! of two because the two of Hartree->Rydberg cancels with the
                 ! half resulting from summing over each sheet pair twice, but
                 ! certainly a factor of 2*pi to get the potential of the i'th sheet.
                 enucsh = enucsh - 2*pi*qjqi*dzabs
              endif
           enddo
        enddo
        ! End charge-sheet -- charged-sheet part.
     enddo
  enddo

  ! If print flag on, print the ion-ion forces.
  if(ipr >= 1) then
     write(7,*)
     write(7,*) 'Forces from ion/charged sheet interaction:'
     write(7,*) '(combined forces)'
     write(7,*) '========================================='
     write(7,22)
     write(7,*)
     do ja = 1, clust%atom_num

        write(7,20) ja,clust%xatm(ja),clust%yatm(ja),clust%zatm(ja), &
             clust%force(:,ja)
     enddo
     write(7,*)
  endif

20 format(i4,2x,3(f11.6,1x),2x,3(f11.6,1x))
22 format('-atom-  ----x----   ----y----   ----z-----', &
        '   ----Fx----  ----Fy----   ----Fz---',/ &
        '                       [bohr]           ', &
        '                  [Ry/bohr]            ')

end subroutine force_charged_sheet
!===============================================================
