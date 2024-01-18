! ==============================================================
!
!     Copyright (C) 2005 Finite Difference Research Group
!     This file is part of parsec, http://www.ices.utexas.edu/parsec/
!
! ==============================================================
!                               FLAGS
! ==============================================================
! The following flags are provided by the user (in cluster.in):

! restart flag
      integer :: istart
! restart from wfn.dat
      logical :: readwfndat
! polarizability flag
      integer :: npolflg
! flag for printout of setup time
      logical :: is_setup
! true if ignoring symmetry ; defaults to false
      logical :: ignoresym
! true if symmetrizing charge density at each SCF iteration
      logical :: chsym
! extrapolation upon restart flag
      integer :: ixtrpflg
! printing flag
      integer :: ipr
! print/export flag
      integer :: export_griddata_flag(MAX_EXPORT_OPTS)
! error flag
      integer :: ierr
! allocation check
      integer :: alcstat
! output info from MPI communication
      integer :: mpinfo
! stop flag
      logical :: ioflag
!
!     output flag for writing parsec.dat file:
!     mod(outflag,2) = 1 : write all calculated wave functions 
!     mod(outflag,2) = 0 : write only the ones for occupied states
!     mod(outflag/2,2) = 1 : write wfn.dat only after SCF is finished
!     mod(outflag/2,2) = 0 : write wfn.dat between steps of SCF loop
!
      integer :: outflag
! additional output flag - outEvFlag
!     output flag for writing eigen.dat file:
!     mod(outevflag,2) = 0 : write all calculated wave functions
!     mod(outevflag,2) = 1 : write only the ones for occupied states
!     mod((outevflag-1)/2,2) = 0 : write eigen.dat only after SCF is finished
!     mod(outevflag/2,2) = 1 : write eigen.dat between steps of SCF loop
!
      integer :: outevflag
! interpolation flag
      logical :: oldinpformat
! flag for GW
      logical :: outputgw
! flag for bypassing/enabling data out
      logical :: enable_data_out 

! ==============================================================
!                   NUMERICAL ACCURACY PARAMETERS
! ==============================================================
! maximal # of iterations before declaring failure to converge
      integer :: mxiter
! convergence criterion
      real(dp) :: vconv
! indicator of convergence-approach
      real(dp) :: vconv_approach
! the variable which holds the SRE chosen to compare with vconv
      real(dp) :: current_sre

! ==============================================================
!                        ELECTRONIC STRUCTURE
! ==============================================================
! total energy in Rydberg, total energy per atom in eV
      real(dp) :: bdev
! exchange-correlation energy
      real(dp) :: exc
! total nuclear energy
      real(dp) :: enuc

! ==============================================================
!                             ATOM DATA
! ==============================================================
! elval - the # of valence electrons for a neutral cluster
      real(dp) :: elval

! ==============================================================
!                          POLARIZABILITY 
! ==============================================================
! polarizability field (in a.u.)
      real(dp) :: field
! polarizability run counter (indicator of field direction)
      integer :: ifield

! ==============================================================
!                           COUNTERS 
! ==============================================================
! counters
      integer :: i, ii, j, iter, irp, isp
! temporary variables
      real(dp) :: tdummy

! ==============================================================
!                       TIMING VARIABLES
! ==============================================================
! total timing
      real(dp) :: tstrt, tfinish, wstrt, wfinish
! Hartree potential timers
      real(dp) :: thart0, thart1, thart_sum
! Diagonalization timers
      real(dp) :: tdiag0, tdiag1, tdiag_sum
! Movement and self-consistent-loop timers
      real(dp) :: tmove0, tmove1, tscf
! Forces timers
      real(dp) :: tforce0, tforce1, tforce_sum
! Setup timers
      real(dp) :: tset0, tset1, tset_sum
! General purpose timer
      real(dp) :: ttimer0,ttimer1,ttimer_sum

! ==============================================================
!                    PERIODIC BOUNDARY CONDITIONS
! ==============================================================
! variable to specify which kpoint number we are looking at 
      integer :: kpnum, kplp

! ==============================================================
!                    LITERAL VARIABLES
! ==============================================================
! date label
      character (len=26) :: datelabel
! string with mixing history file name
      character(len=6) :: fnam
! string for storing the extension of the present node; this is
! concatenated to the end of stdout files to indicate messages
! from the nodes separately
      character(len=4) :: idstring
! ==============================================================
