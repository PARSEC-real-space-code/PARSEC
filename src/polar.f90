!===============================================================
!
! Copyright (C) 2005 Finite Difference Research Group
! This file is part of parsec, http://www.ices.utexas.edu/parsec/
!
! This subroutine prints out the electric field direction, total
! energy, and dipoles, for each of the seven polarizability runs.
! It then updates the run index, and changes the ionic potential
! in the Hamiltonian (pot%vion), FOR THE NEXT RUN, to contain
! a perturbing electric field in the new run index direction - 
! necessary for eventually calculating the polarizabilities 
! (after all runs were completed).
! For details see: Vasiliev et al., PRL, 78, 4805 (1997).
!
! ionic potential, pot%vion - 
! input:
! Initial run - superposition of local pseudpotentials, given on
!               the 3-d grid
! Following runs - same, with contribution of electric field of
!                  previous run
! output:
! same, with contribution of current electric field
!
!---------------------------------------------------------------
subroutine polar(rsymm,elec_st,grid,pot,parallel,natom,field,ifield)

  use constants
  use symmetry_module
  use electronic_struct_module
  use grid_module
  use potential_module
  use parallel_data_module
  implicit none
  !
  ! Input/Output variables:
  !
  ! symmetry operations in reduced group:
  type (symmetry), intent(in) :: rsymm
  ! electronic structure
  type (electronic_struct), intent(in) :: elec_st
  ! grid related data
  type (grid_data), intent(in) :: grid
  ! potential related data
  type (potential), intent(inout) :: pot
  ! parallel computation related data
  type (parallel_data), intent(in) :: parallel
  ! number of atoms in the system
  integer, intent(in) :: natom
  ! polarizability field (in a.u.)
  real(dp), intent(in) :: field
  ! polarizability field index (direction of field) - it being updated
  integer, intent(inout) :: ifield
  !
  ! Work variables:
  !
  ! electric field directions of next run
  real(dp), save :: xfield = zero, yfield = zero, zfield = zero
  ! electric field directions of current (completed) run
  real(dp) :: xfieldold, yfieldold, zfieldold
  ! inverse of the electric field
  real(dp) :: invfld
  ! Cartesian coordinates of each grid point
  real(dp) :: rw(3), rrw(3)
  ! arrays of energies, spin-corrected energies, and dipoles
  ! necessary for calculating the polarizability
  real(dp), save :: enrg(7), dplx(7), dply(7), dplz(7)
  ! polarizabilities based on energy - along each axis and average
  real(dp) :: xalfae, yalfae, zalfae, alfae
  ! polarizabilities based on dipole - along each axis and average
  real(dp) :: xalfad, yalfad, zalfad, alfad
  ! polariazbility per atom
  real(dp) :: alfatom
  ! counters
  integer j, ioffset, itrans
  !---------------------------------------------------------------
  !
  ! write out field direction, energy, and dipole
  !
  if (parallel%iammaster) then
     write(91,93) xfield,yfield,zfield, &
          elec_st%etot,elec_st%dipx,elec_st%dipy,elec_st%dipz
     !
     ! Update components of energy and dipole arrays - convert energies
     ! from Rydberg to Hartree and dipoles from Debye to atomic units.
     !
     enrg(ifield) = elec_st%etot/two
     dplx(ifield) = elec_st%dipx/debye
     dply(ifield) = elec_st%dipy/debye
     dplz(ifield) = elec_st%dipz/debye
  endif

  ifield = ifield+1
  !
  ! Define direction of electric field in x,y,z, directions
  ! according to the run number flag. Also, define the fields used
  ! in the previous run (with zero taken for the first run).
  !
  xfield = zero
  yfield = zero
  zfield = zero
  xfieldold = zero
  yfieldold = zero
  zfieldold = zero
  if(ifield == 1) then
     xfield = zero
     yfield = zero
     zfield = zero
     xfieldold = zero
     yfieldold = zero
     zfieldold = zero
  endif
  if(ifield == 2) then
     xfield = one
     yfield = zero
     zfield = zero
     xfieldold = zero
     yfieldold = zero
     zfieldold = zero
  endif
  if(ifield == 3) then
     xfield = -one
     yfield = zero
     zfield = zero
     xfieldold = one 
     yfieldold = zero
     zfieldold = zero
  endif
  if(ifield == 4) then
     xfield = zero
     yfield = one
     zfield = zero
     xfieldold = -one
     yfieldold = zero
     zfieldold = zero
  endif
  if(ifield == 5) then
     xfield = zero
     yfield = -one
     zfield = zero
     xfieldold = zero
     yfieldold = one
     zfieldold = zero
  endif
  if(ifield == 6) then
     xfield = zero
     yfield = zero
     zfield = one
     xfieldold = zero
     yfieldold = -one
     zfieldold = zero
  endif
  if(ifield == 7) then
     xfield = zero
     yfield = zero
     zfield = -one
     xfieldold = zero
     yfieldold = zero
     zfieldold = one
  endif

  pot%vnew = pot%vold
  !
  ! Subtract contribution of previous electric field and add
  ! contribution of present one to vion.
  !
  ioffset = parallel%irows(parallel%group_iam) - 1
  do j = parallel%irows(parallel%group_iam), &
       parallel%irows(parallel%group_iam+1)-1
     rw(1) = (grid%shift(1) + grid%kx(j)) * grid%step(1)
     rw(2) = (grid%shift(2) + grid%ky(j)) * grid%step(2)
     rw(3) = (grid%shift(3) + grid%kz(j)) * grid%step(3)
     do itrans = 1, rsymm%ntrans
        call symop(rsymm,itrans,rw,rrw)
        pot%vion(j-ioffset) = pot%vion(j-ioffset) &
             + two*field* ((xfield-xfieldold)*rrw(1) + &
                           (yfield-yfieldold)*rrw(2) +  &
                           (zfield-zfieldold)*rrw(3))
     enddo
  enddo

  ! If final run compute the polarizabilities and write results
  ! to polar.dat.

  if (ifield == 8 .and. parallel%iammaster) then
     invfld = one/field
     xalfae = invfld*invfld*(two*enrg(1)-enrg(2)-enrg(3))*cubangs
     yalfae = invfld*invfld*(two*enrg(1)-enrg(4)-enrg(5))*cubangs
     zalfae = invfld*invfld*(two*enrg(1)-enrg(6)-enrg(7))*cubangs

     xalfad = invfld*((dplx(2)-dplx(3))/two)*cubangs
     yalfad = invfld*((dply(4)-dply(5))/two)*cubangs
     zalfad = invfld*((dplz(6)-dplz(7))/two)*cubangs

     alfae = (xalfae + yalfae + zalfae)/three
     alfad = (xalfad + yalfad + zalfad)/three

     write(91,*)
     write(91,94) elec_st%charge
     write(91,*)
     write(91,*) 'Polarizability of the cluster (A**3):'
     write(91,*) '                alfaXX    alfaYY','    alfaZZ       average'
     write(91,95) xalfae, yalfae, zalfae, alfae
     write(91,96) xalfad, yalfad, zalfad, alfad

     write(91,*)
     alfatom = alfae/real(natom,dp)
     write(91,97) alfatom
     alfatom = alfad/real(natom,dp)
     write(91,98) alfatom
  endif

93 format(f4.0,f3.0,f3.0,2x,f17.12,2x,f11.6,f11.6,f11.6)
94 format(1x,'Tot. charge of the cluster: ',f7.3)
95 format(1x,'From energy: ',f10.4,f10.4,f10.4,3x,f10.4)
96 format(1x,'From dipole: ',f10.4,f10.4,f10.4,3x,f10.4)
97 format(1x,'Polarizability per atom =',f10.4,' A**3  (energy)')
98 format('                          ',f10.4,' A**3  (dipole)')

end subroutine polar
!===============================================================
