!===============================================================
!
! Copyright (C) 2005 Finite Difference Research Group
! This file is part of parsec, http://www.ices.utexas.edu/parsec/
!
! =======  Periodic Table =======
! Give the core charge and atomic mass for each element
!
! Atomic masses from NIST document SP966 (physics.nist.gov/data),
! September 2003.
!
! note: if the pseudopotential file already contains core
! charge, that will override the output from this subroutine.
!
!---------------------------------------------------------------
subroutine ptable(element,am,ierr)

  use constants
  implicit none
  !
  ! Input/Output variables:
  !
  ! name of the element
  character (len=2), intent(in) ::  element
  ! charge, mass
  real(dp), intent(out) :: am
  ! error flag
  integer, intent(out) :: ierr
  !---------------------------------------------------------------
  if (element == 'H ') then; am = 1.00794d0
  else if (element == 'H0') then; am = 1.00794d0
  else if (element == 'H1') then; am = 1.00794d0
  else if (element == 'H2') then; am = 1.00794d0
  else if (element == 'H3') then; am = 1.00794d0
  else if (element == 'H4') then; am = 1.00794d0
  else if (element == 'H5') then; am = 1.00794d0
  else if (element == 'H6') then; am = 1.00794d0
  else if (element == 'H7') then; am = 1.00794d0
  else if (element == 'H8') then; am = 1.00794d0
  else if (element == 'H9') then; am = 1.00794d0
  else if (element == 'F0') then; am = 1.00794d0
  else if (element == 'F1') then; am = 1.00794d0
  else if (element == 'F2') then; am = 1.00794d0
  else if (element == 'F3') then; am = 1.00794d0
  else if (element == 'F4') then; am = 1.00794d0
  else if (element == 'F5') then; am = 1.00794d0
  else if (element == 'F6') then; am = 1.00794d0
  else if (element == 'F7') then; am = 1.00794d0
  else if (element == 'F8') then; am = 1.00794d0
  else if (element == 'HI') then; am = 1.00794d0
  else if (element == 'HJ') then; am = 1.00794d0
  else if (element == 'HK') then; am = 1.00794d0
  else if (element == 'HL') then; am = 1.00794d0
  else if (element == 'HM') then; am = 1.00794d0
  else if (element == 'HN') then; am = 1.00794d0
  else if (element == 'HO') then; am = 1.00794d0
  else if (element == 'HP') then; am = 1.00794d0
  else if (element == 'HQ') then; am = 1.00794d0
  else if (element == 'HR') then; am = 1.00794d0
  else if (element == 'FI') then; am = 1.00794d0
  else if (element == 'FJ') then; am = 1.00794d0
  else if (element == 'FK') then; am = 1.00794d0
  else if (element == 'FL') then; am = 1.00794d0
  else if (element == 'FM') then; am = 1.00794d0
  else if (element == 'FN') then; am = 1.00794d0
  else if (element == 'FO') then; am = 1.00794d0
  else if (element == 'FP') then; am = 1.00794d0
  else if (element == 'FQ') then; am = 1.00794d0
  else if (element == 'Li') then; am = 6.941d0
  else if (element == 'Na') then; am = 22.98977d0
  else if (element == 'K ') then; am = 39.0983d0
  else if (element == 'Rb') then; am = 85.4678d0
  else if (element == 'Cs') then; am = 132.90545d0
  else if (element == 'Fr') then; am = 223.0d0
  else if (element == 'Be') then; am = 9.012182d0
  else if (element == 'Mg') then; am = 24.3050d0
  else if (element == 'Ca') then; am = 40.078d0
  else if (element == 'Sr') then; am = 87.62d0 
  else if (element == 'Ba') then; am = 137.327d0 
  else if (element == 'Ra') then; am = 226.0d0
  else if (element == 'B ') then; am = 10.811d0
  else if (element == 'Al') then; am = 26.981538d0
  else if (element == 'Ga') then; am = 69.723d0
  else if (element == 'In') then; am = 114.818d0
  else if (element == 'Tl') then; am = 204.3833d0
  else if (element == 'C ') then; am = 12.0107d0
  else if (element == 'Si') then; am = 28.0855d0
  else if (element == 'Ge') then; am = 72.64d0
  else if (element == 'Sn') then; am = 118.710d0
  else if (element == 'Pb') then; am = 207.2d0
  else if (element == 'N ') then; am = 14.0067d0
  else if (element == 'P ') then; am = 30.973761d0
  else if (element == 'As') then; am = 74.92160d0
  else if (element == 'Sb') then; am = 121.760d0
  else if (element == 'Bi') then; am = 208.98038d0
  else if (element == 'O ') then; am = 15.9994d0
  else if (element == 'S ') then; am = 32.065d0
  else if (element == 'Se') then; am = 78.96d0
  else if (element == 'Te') then; am = 127.60d0
  else if (element == 'Po') then; am = 209.00d0
  else if (element == 'F ') then; am = 18.9984032d0
  else if (element == 'Cl') then; am = 35.453d0
  else if (element == 'Br') then; am = 79.904d0
  else if (element == 'I ') then; am = 126.90447d0
  else if (element == 'At') then; am = 210d0
  else if (element == 'He') then; am = 4.002602d0
  else if (element == 'Ne') then; am = 20.1797d0
  else if (element == 'Ar') then; am = 39.948d0
  else if (element == 'Kr') then; am = 83.798d0
  else if (element == 'Xe') then; am = 131.293d0
  else if (element == 'Rn') then; am = 222.0d0
  else if (element == 'Sc') then; am = 44.955910d0
  else if (element == 'Y ') then; am = 88.90585d0
  else if (element == 'La') then; am = 138.9055d0
  else if (element == 'Ac') then; am = 227.0d0
  else if (element == 'Ti') then; am = 47.867d0
  else if (element == 'Zr') then; am = 91.224d0
  else if (element == 'Hf') then; am = 178.49d0
  else if (element == 'V ') then; am = 50.9415d0
  else if (element == 'Nb') then; am = 92.90638d0
  else if (element == 'Ta') then; am = 180.9479d0
  else if (element == 'Cr') then; am = 51.9961d0
  else if (element == 'Mo') then; am = 95.94d0
  else if (element == 'W ') then; am = 183.84d0
  else if (element == 'Mn') then; am = 54.938049d0
  else if (element == 'Tc') then; am = 98.0d0
  else if (element == 'Re') then; am = 186.207d0
  else if (element == 'Fe') then; am = 55.845d0
  else if (element == 'uF') then; am = 55.845d0
  else if (element == 'Ru') then; am = 101.07d0
  else if (element == 'Os') then; am = 190.23d0
  else if (element == 'Co') then; am = 58.933200d0
  else if (element == 'Rh') then; am = 102.90550d0
  else if (element == 'Ir') then; am = 192.217d0
  else if (element == 'Ni') then; am = 58.6934d0
  else if (element == 'Pd') then; am = 106.42d0
  else if (element == 'Pt') then; am = 195.078d0
  else if (element == 'Cu') then; am = 63.546d0
  else if (element == 'Ag') then; am = 107.8682d0
  else if (element == 'Au') then; am = 196.96655d0
  else if (element == 'Zn') then; am = 65.409d0
  else if (element == 'Cd') then; am = 112.411d0
  else if (element == 'Hg') then; am = 200.59d0
  else if (element == 'Ce') then; am = 140.116d0
  else if (element == 'Pr') then; am = 140.90765d0
  else if (element == 'Nd') then; am = 144.24d0
  else if (element == 'Pr') then; am = 145.0d0
  else if (element == 'Sm') then; am = 150.36d0
  else if (element == 'Eu') then; am = 151.964d0
  else if (element == 'Gd') then; am = 157.25d0
  else if (element == 'Tb') then; am = 158.92534d0
  else if (element == 'Dy') then; am = 162.500d0
  else if (element == 'Ho') then; am = 164.93032d0
  else if (element == 'Er') then; am = 167.259d0
  else if (element == 'Tm') then; am = 168.93421d0
  else if (element == 'Yb') then; am = 173.04d0
  else if (element == 'Lu') then; am = 174.967d0
  else if (element == 'Th') then; am = 232.0381d0
  else if (element == 'Pa') then; am = 231.03588d0
  else if (element == 'U ') then; am = 238.02891d0
  else
     write(7,*) 'ERROR: unknown element: ',element
     write(7,*) 'go to ptable.f and add element'
     write(7,*) 'STOP in ptable'
     ierr = 201
     return
  endif

end subroutine ptable
!===============================================================
!===============================================================
!
! Copyright (C) 2005 Finite Difference Research Group
! This file is part of parsec, http://www.ices.utexas.edu/parsec/
!
!
!---------------------------------------------------------------
subroutine rtable(element,am,ierr)

  use constants
  implicit none
  !
  ! Input/Output variables:
  !
  ! name of the element
  character (len=2), intent(in) ::  element
  ! Wigner-Seitz radii
  real(dp), intent(out) :: am
  ! error flag
  integer, intent(out) :: ierr
  !---------------------------------------------------------------
  if (element == 'H ') then; am = 0.70000d0
  else if (element == 'H0') then; am = 0.70000d0
  else if (element == 'H1') then; am = 0.70000d0
  else if (element == 'H2') then; am = 0.70000d0
  else if (element == 'H3') then; am = 0.70000d0
  else if (element == 'H4') then; am = 0.70000d0
  else if (element == 'H5') then; am = 0.70000d0
  else if (element == 'H6') then; am = 0.70000d0
  else if (element == 'H7') then; am = 0.70000d0
  else if (element == 'H8') then; am = 0.70000d0
  else if (element == 'H9') then; am = 0.70000d0
  else if (element == 'F0') then; am = 0.70000d0
  else if (element == 'F1') then; am = 0.70000d0
  else if (element == 'F2') then; am = 0.70000d0
  else if (element == 'F3') then; am = 0.70000d0
  else if (element == 'F4') then; am = 0.70000d0
  else if (element == 'F5') then; am = 0.70000d0
  else if (element == 'F6') then; am = 0.70000d0
  else if (element == 'F7') then; am = 0.70000d0
  else if (element == 'F8') then; am = 0.70000d0
  else if (element == 'HI') then; am = 0.70000d0
  else if (element == 'HJ') then; am = 0.70000d0
  else if (element == 'HK') then; am = 0.70000d0
  else if (element == 'HL') then; am = 0.70000d0
  else if (element == 'HM') then; am = 0.70000d0
  else if (element == 'HN') then; am = 0.70000d0
  else if (element == 'HO') then; am = 0.70000d0
  else if (element == 'HP') then; am = 0.70000d0
  else if (element == 'HQ') then; am = 0.70000d0
  else if (element == 'HR') then; am = 0.70000d0
  else if (element == 'FI') then; am = 0.70000d0
  else if (element == 'FJ') then; am = 0.70000d0
  else if (element == 'FK') then; am = 0.70000d0
  else if (element == 'FL') then; am = 0.70000d0
  else if (element == 'FM') then; am = 0.70000d0
  else if (element == 'FN') then; am = 0.70000d0
  else if (element == 'FO') then; am = 0.70000d0
  else if (element == 'FP') then; am = 0.70000d0
  else if (element == 'FQ') then; am = 0.70000d0 
  else if (element == 'Li') then; am = 2.60000d0
  else if (element == 'Na') then; am = 3.32000d0
  else if (element == 'K ') then; am = 1
  else if (element == 'Rb') then; am = 1
  else if (element == 'Cs') then; am = 1
  else if (element == 'Fr') then; am = 1
  else if (element == 'Be') then; am = 1
  else if (element == 'Mg') then; am = 1
  else if (element == 'Ca') then; am = 1
  else if (element == 'Sr') then; am = 1
  else if (element == 'Ba') then; am = 1
  else if (element == 'Ra') then; am = 1
  else if (element == 'B ') then; am = 1
  else if (element == 'Al') then; am = 1
  else if (element == 'Ga') then; am = 1
  else if (element == 'In') then; am = 1
  else if (element == 'Tl') then; am = 1
  else if (element == 'C ') then; am = 1
  else if (element == 'Si') then; am = 1
  else if (element == 'Ge') then; am = 1
  else if (element == 'Sn') then; am = 1
  else if (element == 'Pb') then; am = 1
  else if (element == 'N ') then; am = 1
  else if (element == 'P ') then; am = 1
  else if (element == 'As') then; am = 1
  else if (element == 'Sb') then; am = 1
  else if (element == 'Bi') then; am = 1
  else if (element == 'O ') then; am = 1
  else if (element == 'S ') then; am = 1
  else if (element == 'Se') then; am = 1
  else if (element == 'Te') then; am = 1
  else if (element == 'Po') then; am = 1
  else if (element == 'F ') then; am = 1
  else if (element == 'Cl') then; am = 1
  else if (element == 'Br') then; am = 1
  else if (element == 'I ') then; am = 1
  else if (element == 'At') then; am = 1
  else if (element == 'He') then; am = 1
  else if (element == 'Ne') then; am = 1
  else if (element == 'Ar') then; am = 1
  else if (element == 'Kr') then; am = 1
  else if (element == 'Xe') then; am = 1
  else if (element == 'Rn') then; am = 1
  else if (element == 'Sc') then; am = 1
  else if (element == 'Y ') then; am = 1
  else if (element == 'La') then; am = 1
  else if (element == 'Ac') then; am = 1
  else if (element == 'Ti') then; am = 1
  else if (element == 'Zr') then; am = 1
  else if (element == 'Hf') then; am = 1
  else if (element == 'V ') then; am = 1
  else if (element == 'Nb') then; am = 1
  else if (element == 'Ta') then; am = 1
  else if (element == 'Cr') then; am = 2.50000d0
  else if (element == 'Mo') then; am = 1
  else if (element == 'W ') then; am = 1
  else if (element == 'Mn') then; am = 1
  else if (element == 'Tc') then; am = 1
  else if (element == 'Re') then; am = 1
  else if (element == 'Fe') then; am = 2.46000d0
  else if (element == 'uF') then; am = 1
  else if (element == 'Ru') then; am = 1
  else if (element == 'Os') then; am = 1
  else if (element == 'Co') then; am = 1
  else if (element == 'Rh') then; am = 1
  else if (element == 'Ir') then; am = 1
  else if (element == 'Ni') then; am = 1
  else if (element == 'Pd') then; am = 1
  else if (element == 'Pt') then; am = 1
  else if (element == 'Cu') then; am = 1
  else if (element == 'Ag') then; am = 1
  else if (element == 'Au') then; am = 1
  else if (element == 'Zn') then; am = 1
  else if (element == 'Cd') then; am = 1
  else if (element == 'Hg') then; am = 1
  else if (element == 'Ce') then; am = 1
  else if (element == 'Pr') then; am = 1
  else if (element == 'Nd') then; am = 1
  else if (element == 'Pr') then; am = 1
  else if (element == 'Sm') then; am = 1
  else if (element == 'Eu') then; am = 1
  else if (element == 'Gd') then; am = 1
  else if (element == 'Tb') then; am = 1
  else if (element == 'Dy') then; am = 1
  else if (element == 'Ho') then; am = 1
  else if (element == 'Er') then; am = 1
  else if (element == 'Tm') then; am = 1
  else if (element == 'Yb') then; am = 1
  else if (element == 'Lu') then; am = 1
  else if (element == 'Th') then; am = 1
  else if (element == 'Pa') then; am = 1
  else if (element == 'U ') then; am = 1
  else
     write(7,*) 'ERROR: unknown element: ',element
     write(7,*) 'go to ptable.f and add element'
     write(7,*) 'STOP in ptable'
     ierr = 201
     return
  endif

end subroutine rtable
!===============================================================
