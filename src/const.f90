!===============================================================
!
!  Copyright (C) 2005 Finite Difference Research Group
!  This file is part of parsec, http://www.ices.utexas.edu/parsec/
!
!  NUMERICAL CONSTANTS
!  This file should NEVER be touched, barring new
!  and exciting developments in the value of, say, pi.
!
!---------------------------------------------------------------
module constants

  implicit none
  !
  !  Definitions of double precision real/complex types:
  !
  integer, parameter :: dp = kind(1.0d0)
  integer, parameter :: dpc = kind((1.0d0,1.0d0))
  !  Machine precision for 64bit 2**-52
  real(dp), parameter :: mach_eps    = 2.220446049250313080847263336182d-16
  real(dp), parameter :: twomach_eps = 4.440892098500626161694526672363d-16
  !
  !  Numerical constants:
  !
  real(dp), parameter :: zero = 0.d0, one = 1.d0, &
       two = 2.d0, three = 3.d0, four = 4.d0, five = 5.d0, &
       six = 6.d0, seven = 7.d0, eight = 8.d0, nine = 9.d0, mone = -1.d0, &
       half = 0.5d0, third = 1.d0/3.d0, mhalf = -0.5d0, &
       root3 = 1.73205080756887729352744634151d0, &
       pi = 3.1415926535897932384626433832795d0, &
       twopi = two * pi
  complex(dpc), parameter :: zzero = (0.d0,0.d0), &
       zone = (1.d0,0.d0), zi = (0.d0,1.d0)
  !
  !  Physical constants:
  !
  !  convert energy units from rydbergs to eV
  real(dp), parameter :: rydberg = 13.6058d0
  !  convert length units from atomic units (bohrs) to angstroms
  real(dp), parameter :: angs = 0.529177d0, cubangs = angs**3
  !  convert dipole units from au to debyes
  real(dp), parameter :: debye = 2.54d0
  !  convert temperature from kelvins to rydbergs
  real(dp), parameter :: tempkry = 6.33327186d-06
  !  convert time from a.u. to ps
  real(dp), parameter :: timeaups = 2.4188843d-05
  !
  !  Parsec constants:
  !
  !  constants for minimization scheme
  integer, parameter :: NONE = 0, STEEPDESC = 1, BFGS = 2, MANUAL = 3
  !  constants for mixing scheme
  integer, parameter :: ANDERSON = 0, BROYDEN = 1 &
       ,MSECANT1 = 2, MSECANT2 = 3, MSECANT3 = 4
  !  constants for diagonalization scheme
  integer, parameter :: TEST=-1, ARPACK=0, DIAGLA=1, TRLANC=2, CHEBDAV=3, CHEBFF=4 !!added CHEBFF
  !  constants for pseudopotential format
  integer, parameter :: MARTINS_OLD = 0, MARTINS_NEW = 1, &
       MARTINS_WANG = 2, FHIPP = 3
  !  constants for k-point sampling options
  integer, parameter :: KP_MANUAL = 1, MONKHORST_PACK = 2
  !  constants for molecular dynamics scheme
  integer, parameter :: MD_STAIR = 1, MD_LINEAR = 2, MD_LOG = 3
  !  constants for grid export routine. !! Relies also on NONE= 0 assigned above !!
  integer, parameter :: CHGDENS = 1, VEXT = 2, VHAR = 3, VEXTHAR = 23, MAX_EXPORT_OPTS = 10
  !  ratio between the dispersion of Gaussian distributions and
  !  pseudopotential cut-off radius.
  real(dp) :: sigmafac = half

end module constants
!===============================================================
