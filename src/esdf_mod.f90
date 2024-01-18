!===============================================================
!
! e l e c t r o n i c   s t r u c t u r e   d a t a   f o r m a t
! ---------------------------------------------------------------
!
!                            e s d f
!                            =======
!
! Author: Chris J. Pickard (c)
! Email : cp@min.uni-kiel.de
! Place : kiel, Germany
! Date  : 5/6th august 1999
!
! Summary
! -------
!
! This module is designed to simplify and enhance the input of data into
! electronic structure codes (for example, castep). It works from a
! fashion, and input is independent of the ordering of the input
! file. An important feature is the requirement that most inputs require
! default settings to be supplied within the main program calling
! esdf. This means that rarely used variables will not clutter everyday
! input files, and, even more usefully, "intelligence" may be built into
! the main code as the defaults may be dependent of other set
! variables. Block data may also be read in. Another important feature
! is the ability to define "physical" values. This means that the input
! files need not depend on the internal physical units used by the main
! program.
!
! History
! -------
!
! Esdf has been written from scratch in f90, but is heavily based
! (especially for the concept) on the fdf package developed by alberto
! garcia and jose soler. It is not as "flexible" as fdf - there is no
! provision for branching to further input files. This simplifies the
! code, and I hope that it is still useful without this feature. Also,
! the input and defaults are not dumped to a output file currently. I've
! not found this a hindrance as of now.
!
! Future
! ------
!
! My intention is to make this release available to alberto garcia and
! jose soler for their comments. It might be a good idea to use this as
! a base for fully converting the fdf package to f90. Or it may remain
! as a cut down version of fdf. I certainly hope that a package of the
! fdf sort becomes widely used in the electronic structure community. My
! experience has been very positive.
!
! Usage
! -----
!
! First, "use esdf" wherever you wish to make use of its features. In
! the main program call the initialisation routine: call
! esdf_init('input.esdf'). "input.esdf" is the name of the input file -
! it could be anything. This routine opens the input file, and reads
! into a dynamically allocated storage array. The comments and blank
! lines are stripped out. You are now ready to use the
! esdf_functions. For example, if you want to read in the number of
! atoms in your calculation, you would use: natom =
! esdf_integer('numberofatoms',1), where 'numberofatoms' is the label to
! search for in the input file, and '1' is the default. Call esdf_close to
! deallocate the data arrays. You may then open another input file using
! esdf_init. It is not currently possible to open more that on input
! file at one time.
!
! Syntax
! ------
!
! The input file can contain comments. These are defined as anything to
! the right of, and including, '#', ';', or '!'. It is straightforward
! to modify the routine to accept further characters. Blank lines are
! ignored -- use comments and blank lines to make you input file
! readable.
!
! The "labels" are case insensitive (e.g. unitcell is equivalent to
! unitcell) and punctuation insensitive (unit.cell is equivalent to
! unit_cell is equivalent to unitcell). Punctuation characters are '.'
! and '-' at the moment. Again - use this feature to improve
! readability.
!
! The following are equivalent ways of defining a physical quantity:
!
! "Ageofuniverse = 24.d0 s" or "ageofuniverse : 24.d0 s" or
! "ageofuniverse 24.d0 s"
!
! It would be read in by the main program in the following way:
!
! Aou = esdf_physical('ageofuniverse',77.d0,'ns')
!
! "Aou" is the double precision variable, 77.d0 is the default number of
! "ns" or nanoseconds. 24s will be converted automatically to its
! equivalent number of nanoseconds.
!
! Block data should be placed in the input file as follows:
!
! Begin cellvectors
! 1.0 1.0 0.0
! 0.0 1.0 1.0
! 1.0 0.0 1.0
! end cellvectors
!
! And it may be read:
!
!   If(esdf_block('cellvectors',nlines))
!     if(nlines /= 3) then (... break out here if the incorrect number
! of lines)
!     do i=1,nlines
!       read(block_data(i),*) x,y,z
!     end do
!   endif
!
! List of functions
! -----------------
!
! Self explanatory:
!
! Esdf_string(label,default)
! esdf_integer(label,default)
! esdf_single(label,default)
! esdf_double(label,default)
! esdf_physical(label,default,unit)
!
! a little more explanation:
!
! Esdf_defined(label) is true if "label" found, false otherwise
!
! Esdf_boolean(label,default) is true if "label yes/true/t (case/punct.insens)
!                             is false if"label no/false/f (case/punct.insens)
!
! The help feature
! ----------------
!
! The routine "esdf_help(helpword,searchword)" can be used to access the
! information contained within the "esdf_key_mod" module.
!
! If "helpword" is "search" (case insensitive), then all variables whose
! description contains "searchword" will be output.
!
! If "helpword" is "basic", "inter", "expert" or "dummy" the varibles of
! that type will be displayed.
!
! If "helpword" is one of the valid labels, then a description of this
! label will be output.
!
! Finishing off
! -------------
!
! Two routines, "esdf_warnout" and "esdf_close", can be used to finish
! the use of esdf. "esdf_warnout" outputs esdf warnings to screen, and
! "esdf_close" deallocates the allocated esdf arrays.
!
! Contact the author
! ------------------
!
! This code is under development, and the author would be very happy to
! receive comments by email. Any use in a commercial software package is
! forbidden without prior arrangement with the author (Chris J. Pickard).
!---------------------------------------------------------------

module esdf

  use esdf_key

  implicit none

  ! Kind parameters

  integer, private, parameter :: i4b = kind (1), &
       dp  = kind (1.0d0), sp  = kind (1.0)

  ! Set the length of the lines

  integer (i4b), public, parameter  :: llength = 80
  integer (i4b), private, parameter :: nphys = 57

  integer (i4b), private :: nrecords, nwarns

  character (len=llength), private, allocatable :: &
       llist(:),warns(:), tlist(:,:)

  ! The public block data array

  character (len=llength), public, allocatable :: block_data(:)

  ! Set the physical units database

  character (len=10), private :: &
       phy_d(nphys), phy_n(nphys) ! D - dimension n - name
  real (dp), private :: phy_u(nphys) ! U - unit

  ! We allow case variations in the units. This could be dangerous
  ! (mev --> mev !!) in real life, but not in this restricted field.

  ! m - mass l - length t - time e - energy f - force p - pressure c - charge
  ! d - dipole mom - mom inert ef - efield

  data phy_d(1)       / 'm' /
  data phy_n(1)       / 'kg' /
  data phy_u(1)       / 1.0_dp /
  data phy_d(2)       / 'm' /
  data phy_n(2)       / 'g' /
  data phy_u(2)       / 1.0e-3_dp /
  data phy_d(3)       / 'm' /
  data phy_n(3)       / 'amu' /
  data phy_u(3)       / 1.66054e-27_dp /
  data phy_d(4)       / 'l' /
  data phy_n(4)       / 'm' /
  data phy_u(4)       / 1.0_dp /
  data phy_d(5)       / 'l' /
  data phy_n(5)       / 'nm' /
  data phy_u(5)       / 1.0e-9_dp /
  data phy_d(6)       / 'l' /
  data phy_n(6)       / 'ang' /
  data phy_u(6)       / 1.0e-10_dp /
  data phy_d(7)       / 'l' /
  data phy_n(7)       / 'bohr' /
  data phy_u(7)       / 0.529177e-10_dp /
  data phy_d(8)       / 't' /
  data phy_n(8)       / 's' /
  data phy_u(8)       / 1.0_dp /
  data phy_d(9)       / 't' /
  data phy_n(9)       / 'ns' /
  data phy_u(9)       / 1.0e-9_dp /
  data phy_d(10)      / 't' /
  data phy_n(10)      / 'ps' /
  data phy_u(10)      / 1.0e-12_dp /
  data phy_d(11)      / 't' /
  data phy_n(11)      / 'fs' /
  data phy_u(11)      / 1.0e-15_dp /
  data phy_d(12)      / 'e' /
  data phy_n(12)      / 'j' /
  data phy_u(12)      / 1.0_dp /
  data phy_d(13)      / 'e' /
  data phy_n(13)      / 'erg' /
  data phy_u(13)      / 1.0e-7_dp /
  data phy_d(14)      / 'e' /
  data phy_n(14)      / 'ev' /
  data phy_u(14)      / 1.60219e-19_dp /
  data phy_d(15)      / 'e' /
  data phy_n(15)      / 'mev' /
  data phy_u(15)      / 1.60219e-22_dp /
  data phy_d(16)      / 'e' /
  data phy_n(16)      / 'ry' /
  data phy_u(16)      / 2.17991e-18_dp /
  data phy_d(17)      / 'e' /
  data phy_n(17)      / 'mry' /
  data phy_u(17)      / 2.17991e-21_dp /
  data phy_d(18)      / 'e' /
  data phy_n(18)      / 'hartree' /
  data phy_u(18)      / 4.35982e-18_dp /
  data phy_d(19)      / 'e' /
  data phy_n(19)      / 'kcal/mol' /
  data phy_u(19)      / 6.94780e-21_dp /
  data phy_d(20)      / 'e' /
  data phy_n(20)      / 'mhartree' /
  data phy_u(20)      / 4.35982e-21_dp /
  data phy_d(21)      / 'e' /
  data phy_n(21)      / 'kj/mol' /
  data phy_u(21)      / 1.6606e-21_dp /
  data phy_d(22)      / 'e' /
  data phy_n(22)      / 'hz' /
  data phy_u(22)      / 6.6262e-34_dp /
  data phy_d(23)      / 'e' /
  data phy_n(23)      / 'thz' /
  data phy_u(23)      / 6.6262e-22_dp /
  data phy_d(24)      / 'e' /
  data phy_n(24)      / 'cm-1' /
  data phy_u(24)      / 1.986e-23_dp /
  data phy_d(25)      / 'e' /
  data phy_n(25)      / 'cm^-1' /
  data phy_u(25)      / 1.986e-23_dp /
  data phy_d(26)      / 'e' /
  data phy_n(26)      / 'cm**-1' /
  data phy_u(26)      / 1.986e-23_dp /
  data phy_d(27)      / 'f' /
  data phy_n(27)      / 'N' /
  data phy_u(27)      / 1.0_dp /
  data phy_d(28)      / 'f' /
  data phy_n(28)      / 'ev/ang' /
  data phy_u(28)      / 1.60219e-9_dp /
  data phy_d(29)      / 'f' /
  data phy_n(29)      / 'ry/bohr' /
  data phy_u(29)      / 4.11943e-8_dp /
  data phy_d(30)      / 'l' /
  data phy_n(30)      / 'cm' /
  data phy_u(30)      / 1.0e-2_dp /
  data phy_d(31)      / 'p' /
  data phy_n(31)      / 'pa' /
  data phy_u(31)      / 1.0_dp /
  data phy_d(32)      / 'p' /
  data phy_n(32)      / 'mpa' /
  data phy_u(32)      / 1.0e6_dp /
  data phy_d(33)      / 'p' /
  data phy_n(33)      / 'gpa' /
  data phy_u(33)      / 1.0e9_dp /
  data phy_d(34)      / 'p' /
  data phy_n(34)      / 'atm' /
  data phy_u(34)      / 1.01325e5_dp /
  data phy_d(35)      / 'p' /
  data phy_n(35)      / 'bar' /
  data phy_u(35)      / 1.0e5_dp /
  data phy_d(36)      / 'p' /
  data phy_n(36)      / 'mbar' /
  data phy_u(36)      / 1.0e11_dp /
  data phy_d(37)      / 'p' /
  data phy_n(37)      / 'ry/bohr**3' /
  data phy_u(37)      / 1.47108e13_dp /
  data phy_d(38)      / 'p' /
  data phy_n(38)      / 'ev/ang**3' /
  data phy_u(38)      / 1.60219e11_dp /
  data phy_d(39)      / 'c' /
  data phy_n(39)      / 'c' /
  data phy_u(39)      / 1.0_dp /
  data phy_d(40)      / 'c' /
  data phy_n(40)      / 'e' /
  data phy_u(40)      / 1.602177e-19_dp /
  data phy_d(41)      / 'd' /
  data phy_n(41)      / 'C*m' /
  data phy_u(41)      / 1.0_dp /
  data phy_d(42)      / 'd' /
  data phy_n(42)      / 'D' /
  data phy_u(42)      / 3.33564e-30_dp /
  data phy_d(43)      / 'd' /
  data phy_n(43)      / 'debye' /
  data phy_u(43)      / 3.33564e-30_dp /
  data phy_d(44)      / 'd' /
  data phy_n(44)      / 'e*bohr' /
  data phy_u(44)      / 8.47835e-30_dp /
  data phy_d(45)      / 'd' /
  data phy_n(45)      / 'e*ang' /
  data phy_u(45)      / 1.602177e-29_dp /
  data phy_d(46)      / 'mom' /
  data phy_n(46)      / 'kg*m**2' /
  data phy_u(46)      / 1.0_dp /
  data phy_d(47)      / 'mom' /
  data phy_n(47)      / 'ry*fs**2' /
  data phy_u(47)      / 2.1799e-48_dp /
  data phy_d(48)      / 'ef' /
  data phy_n(48)      / 'v/m' /
  data phy_u(48)      / 1.0_dp /
  data phy_d(49)      / 'ef' /
  data phy_n(49)      / 'v/nm' /
  data phy_u(49)      / 1.0e9_dp /
  data phy_d(50)      / 'ef' /
  data phy_n(50)      / 'v/ang' /
  data phy_u(50)      / 1.0e10_dp /
  data phy_d(51)      / 'ef' /
  data phy_n(51)      / 'v/bohr' /
  data phy_u(51)      / 1.8897268e10_dp /
  data phy_d(52)      / 'ef' /
  data phy_n(52)      / 'ry/bohr/e' /
  data phy_u(52)      / 2.5711273e11_dp /
  data phy_d(53)      / 'ef' /
  data phy_n(53)      / 'har/bohr/e' /
  data phy_u(53)      / 5.1422546e11_dp /
  data phy_d(54)      / 'e' /
  data phy_n(54)      / 'k' /
  data phy_u(54)      / 1.38066e-23_dp /
  data phy_d(55)      / 'b' /
  data phy_n(55)      / 't' /
  data phy_u(55)      / 1.0_dp /
  data phy_d(56)      / 'b' /
  data phy_n(56)      / 'ry/mu_bohr' /
  data phy_u(56)      / 2.350499e5_dp /
  data phy_d(57)      / 'b' /
  data phy_n(57)      / 'g' /
  data phy_u(57)      / 1.0e-4_dp /

contains

  !   --------------  esdf_init  ---------------------- 

  subroutine esdf_init (filename)

    character (len=*), intent (in) :: filename

  ! Local

    integer (i4b), parameter :: ncomm = 3, ndiv  = 3

    integer (i4b)       :: unit, ierr, i, j, ic, nt, ndef, nread, &
         itemp, itemp2
    character (len=llength) :: cjunk, ctemp
    character (len=1)   :: comment(ncomm), divide(ndiv)
    logical             :: inblock ! Integer :: temp1!, index 
                                 !  djr:  had to change this to get
                                 !  it to compile under linux.

  ! Define comment characters

    data comment        / '#', ';', '!' /
    data divide         / ' ', '=', ':' /


  ! "Reduce" the keyword list for comparison

    do i = 1, numkw
       ctemp       = kw_label(i)
       kw_label(i) = esdf_reduce(ctemp)
    enddo


  ! Initializing the array kw_index
    kw_index = 1
         

 ! Open the esdf file

    call esdf_file (unit, filename, ierr)
    cjunk = 'Unable to open main input file "'//trim (filename)// '"'

    if (ierr == 1) then
       write (*,*) 'ESDF WARNING: '//trim (cjunk)// ' - using defaults'
       write (*,*) 'I actually should abort now - but guess what?'
       nread = 0
    else
       nread = huge (1)
    endif

  ! Count the number of records (excluding blank lines and commented lines)

    nrecords = 0

    do i = 1, nread
       read (unit, '(a)', end = 10) cjunk
       do j = 1, ncomm
          ic = index(cjunk, comment(j))
          if (ic > 0) cjunk(ic:) = ' '
       enddo
       if (len_trim (cjunk) > 0) nrecords = nrecords + 1
    enddo
10  continue
    rewind (unit)

  ! Allocate the array to hold the records and tokens

    allocate(llist(nrecords), block_data(nrecords), &
         tlist(llength,nrecords), warns(nrecords))

  ! Set the number of warnings to zero

    nwarns = 0
    warns  = ' '

  ! Read in the records

    nrecords = 0
    do i = 1, nread
       read (unit, '(a)', end = 11) cjunk
       do j = 1, ncomm
          ic = index(cjunk, comment(j))
          if (ic > 0) cjunk(ic:) = ' '
       enddo
       if (len_trim (cjunk) > 0) then
          nrecords = nrecords + 1
          llist(nrecords) = adjustl(cjunk)
       endif
    enddo
11  continue
    close (unit)

 ! Now read in the tokens from llist

    tlist = ' '

    do i = 1, nrecords
       ctemp = llist(i)
       nt    = 0
       do while (len_trim (ctemp) > 0)

  ! Apparently this a hack to make it compile under linux
  !        temp1=index(ctemp,divide(1))
  !        ic = minval(temp1)
  !        ic = minval(index(ctemp,divide),mask=index(ctemp,divide)>0)

          ic = len_trim (ctemp) + 1
          do itemp = 1, size(divide)
             itemp2 = index(ctemp, divide(itemp))
             if (itemp2 == 0) itemp2 = len_trim (ctemp) + 1
             if (itemp2 < ic) ic = itemp2
          enddo
          if (ic > 1) then
             nt = nt + 1
             tlist(nt,i) = adjustl(ctemp(:ic-1))
          endif
          ctemp = adjustl(ctemp(ic+1:))
       enddo
    enddo

  ! Check if any of the "labels" in the input file are unrecognised

    inblock = .false.
    do i = 1, nrecords

  ! Check if we are in a block

       if (esdf_reduce(tlist(1,i)) == 'begin') then
          inblock = .true.

  ! Check if block label is recognised

          if ((count(esdf_reduce(tlist(2,i)) == kw_label) == 0)) then
             ctemp = 'Label "'//trim (esdf_reduce(tlist(2,i)))// &
                  '" not in keyword list'
             if (count(ctemp == warns) == 0) call esdf_warn (ctemp)
          endif

  ! Check if "label" is multiply defined in the input file

          ndef = 0
          do j = 1, nrecords
             if (esdf_reduce(tlist(2,i)) == &
                  esdf_reduce(tlist(2,j))) ndef = ndef + 1
          enddo
          ctemp = 'Label "'//trim (esdf_reduce(tlist(2,i)))// &
               '" is multiply defined in the input file. '
          if ((ndef > 2) .and. &
               (count(ctemp == warns) == 0)) call esdf_warn (ctemp)
       endif

  ! Check it is in the list of keywords

       if ((count(esdf_reduce(tlist(1,i)) == kw_label) == &
            0) .and. (.not. inblock)) then 
          ctemp = 'Label "'//trim (esdf_reduce(tlist(1,i)))// &
               '" not in keyword list'
          if (count(ctemp == warns) == 0) call esdf_warn (ctemp)
       endif
       if (.not. inblock) then

  ! Check if "label" is multiply defined in the input file

          ndef = 0
          do j = 1, nrecords
             if (esdf_reduce(tlist(1,i)) == esdf_reduce(tlist(1,j))) &
                  ndef = ndef + 1
          enddo
          ctemp = 'Label "'//trim (esdf_reduce(tlist(1,i)))// &
               '" is multiply defined in the input file. '
          if ((ndef > 1) .and. (count(ctemp == warns) == 0)) &
               call esdf_warn (ctemp)
       endif

  ! Check if we have left a block

       if (esdf_reduce(tlist(1,i)) == 'end') inblock = .false.

    enddo

    return
  end subroutine esdf_init

  !   --------------  esdf_string  ---------------------- 

  ! Return the string attached to the "label"

  function esdf_string (label, default)

    character(len=*), intent(in) :: label, default
    character(len=llength) :: esdf_string
    ! Local
    integer(i4b)      :: i
    integer(i4b)      :: kw_number

    ! Check "label" is defined

    call esdf_lablchk (label, 'T',kw_number)

    ! Set to default

    esdf_string = default

    do i = kw_index(kw_number), nrecords
  ! Search in the first token for "label"
  ! the first instance is returned
       if (esdf_reduce(tlist(1,i)) == esdf_reduce(label)) then
          esdf_string = llist(i)(index(llist(i),trim(tlist(2,i))):)
          kw_index(kw_number) = i+1
          exit
       endif

    enddo

    return
  end function esdf_string

  !   --------------  esdf_integer  ---------------------- 

  ! Return the integer attached to the "label"

  function esdf_integer (label, default)

    integer(i4b), intent(in)     :: default
    character(len=*), intent(in) :: label
    integer(i4b)      :: esdf_integer
  ! Local
    integer(i4b)       :: i
    character(len=llength) :: ctemp
    integer(i4b)      :: kw_number

  ! Check "label" is defined
    call esdf_lablchk (label, 'I',kw_number)
        
  ! Set to default
    esdf_integer = default

    do i = kw_index(kw_number) , nrecords
  ! Search in the first token for "label"
  ! the first instance is returned
       if (esdf_reduce(tlist(1,i)) == esdf_reduce(label)) then
          read (tlist(2,i), *, err = 10) esdf_integer
          kw_index(kw_number) = i+1
          exit
       endif
    enddo
    return

10  continue
    ctemp = 'Unable to parse "'//trim (esdf_reduce(label))// &
         '" in esdf_integer'
    call esdf_die (ctemp)

    return
  end function esdf_integer

  !   --------------  esdf_single  ---------------------- 

  ! Return the single precisioned value attached to the "label"

  function esdf_single (label, default)

    real(sp), intent(in) :: default
    character(len=*), intent(in) :: label
    real(sp)          :: esdf_single
  ! Local
    integer(i4b)       :: i
    character(len=llength) :: ctemp
    integer(i4b)      :: kw_number

  ! Check "label" is defined

    call esdf_lablchk (label, 'S',kw_number)

  ! Set to default

    esdf_single = default

    do i = kw_index(kw_number), nrecords
  ! Search in the first token for "label"
  ! the first instance is returned
       if (esdf_reduce(tlist(1,i)) == esdf_reduce(label)) then
          read (tlist(2,i), *, err = 10) esdf_single
          kw_index(kw_number) = i+1
          exit
       endif
    enddo
    return

10  continue
    ctemp = 'Unable to parse "'//trim (esdf_reduce(label))// &
                '" in esdf_single'
    call esdf_die (ctemp)

    return
  end function esdf_single

  !   --------------  esdf_double  ---------------------- 

  ! Return the double precisioned value attached to the "label"

  function esdf_double (label, default)

    real(dp), intent(in) :: default
    character(len=*), intent(in) :: label
    real(dp)          :: esdf_double
  ! Local
    integer(i4b)       :: i
    character(len=llength) :: ctemp
    integer(i4b)      :: kw_number

  ! Check "label" is defined

    call esdf_lablchk (label, 'D',kw_number)

  ! Set to default

    esdf_double = default

    do i = kw_index(kw_number), nrecords
  ! Search in the first token for "label"
  ! the first instance is returned
       if (esdf_reduce(tlist(1,i)) == esdf_reduce(label)) then
          read (tlist(2,i), *, err = 10) esdf_double
          kw_index(kw_number) = i+1
          exit
       endif
    enddo

    return

10  continue
    esdf_double = default
    ctemp       = 'Unable to parse "'//trim (esdf_reduce(label))// &
                      '" in esdf_double'
    call esdf_die (ctemp)

    return
  end function esdf_double

  !   --------------  esdf_physical  ---------------------- 

  ! Return the double precisioned physical value attached to the "label"
  ! units converted to "dunit"

  function esdf_physical (label, default, dunit)

    real(dp), intent(in) :: default
    character(len=*), intent(in) :: label, dunit
    real(dp)          :: esdf_physical
  ! Local
    integer(i4b)       :: i
    character(len=llength) :: ctemp, iunit
    integer(i4b)      :: kw_number

  ! Check "label" is defined

    call esdf_lablchk (label, 'P',kw_number)

  ! Set to default

    esdf_physical = default

    do i = kw_index(kw_number), nrecords

  ! Search in the first token for "label"
  ! the first instance is returned

       if (esdf_reduce(tlist(1,i)) == esdf_reduce(label)) then
          read (tlist(2,i), *, err = 10) esdf_physical
          iunit = dunit//repeat(' ', llength - len(dunit))
          read (tlist(3,i), *, err = 10, end = 13) iunit
          esdf_physical = esdf_convfac(iunit,dunit) * esdf_physical
          kw_index(kw_number) = i+1
          exit
       endif
    enddo
    return
13  continue
    iunit = dunit//repeat(' ', llength - len(dunit))
    esdf_physical = esdf_convfac(iunit,dunit) * esdf_physical
    kw_index(kw_number) = i+1

    return

10  continue
    esdf_physical = default
    ctemp = 'Unable to parse "'//trim (esdf_reduce(label))// &
         '" in esdf_physical'
    call esdf_die (ctemp)

    return
  end function esdf_physical

  !   --------------  esdf_defined  ---------------------- 

  ! Is the "label" defined in the input file

  function esdf_defined (label)

    character(len=*), intent(in) :: label
    logical            :: esdf_defined
  ! Local
    integer(i4b)      :: i
    integer(i4b)      :: kw_number

  ! Check "label" is defined
    call esdf_lablchk (label, 'E',kw_number)       

  ! Set to default
    esdf_defined = .false.

    do i = kw_index(kw_number), nrecords
  ! Search in the first token for "label"
  ! the first instance is returned
       if (esdf_reduce(tlist(1,i)) == esdf_reduce(label)) then
          esdf_defined = .true.
          kw_index(kw_number) = i+1
          exit
       endif

    enddo

    return
  end function esdf_defined

  !   --------------  esdf_boolean  ---------------------- 

  ! Is the "label" defined in the input file

  function esdf_boolean (label, default)

    character(len=*), intent(in) :: label
    logical, intent(in) :: default
    logical            :: esdf_boolean
  ! Local
    integer(i4b)       :: i
    character(len=llength) :: positive(3), negative(3)
    data positive       / 'yes', 'true', 't' /
    data negative       / 'no', 'false', 'f' /
    integer(i4b)      :: kw_number

  ! Check "label" is defined

    call esdf_lablchk (label, 'L',kw_number)

  ! Set to default
    esdf_boolean = default

    do i = kw_index(kw_number), nrecords
  ! Search in the first token for "label"
  ! the first instance is returned
       if (esdf_reduce(tlist(1,i)) == esdf_reduce(label)) then
          kw_index(kw_number) = i+1
          if (len_trim (tlist(2,i)) == 0) then
             esdf_boolean = .true.
             exit
          endif
          if (any(index(positive, esdf_reduce(tlist(2,i))) > 0)) then
             esdf_boolean = .true.
             exit
          endif
          if (any(index(negative, esdf_reduce(tlist(2,i))) > 0)) then
             esdf_boolean = .false.
             exit
          endif
          call esdf_die ('Unable to parse boolean value')

       endif

    enddo

    return
  end function esdf_boolean

  !   --------------  esdf_block  ---------------------- 
  function esdf_block (label, nlines)

    character(len=*), intent(in) :: label
    integer(i4b), intent(out)    :: nlines
    logical            :: esdf_block
  ! Local
    integer(i4b)       :: i
    character(len=llength) :: ctemp
    integer(i4b)      :: kw_number

  ! Check "label" is defined
    call esdf_lablchk (label, 'B',kw_number)

    ctemp = 'Block "'//trim (esdf_reduce(label))// &
         '" not closed correctly '

    esdf_block = .false.

    nlines = 0

    do i = kw_index(kw_number), nrecords
       if ((esdf_reduce(tlist(1,i)) == esdf_reduce('begin')) .and. &
            (esdf_reduce(tlist(2,i)) == esdf_reduce(label))) then
          esdf_block = .true.
          kw_index(kw_number) = i+1
          do while (esdf_reduce(tlist(1,i+nlines+1)) /= &
               esdf_reduce('end'))
             nlines = nlines + 1
             if (nlines + i > nrecords) call esdf_die (ctemp)
             block_data(nlines) = llist(i+nlines)
          enddo
          if (esdf_reduce(tlist(2,i+nlines+1)) /= &
               esdf_reduce(label)) call esdf_die (ctemp)
          exit
       endif
    enddo

    return
  end function esdf_block

  !   --------------  esdf_reduce  ---------------------- 

  ! Reduce the string to lower case and remove punctuation

  function esdf_reduce (string)

    character(len=*), intent(in) :: string
    character(len=llength) :: esdf_reduce
  ! Local
    integer(i4b), parameter :: npunct = 2
    integer(i4b)       :: ia, iz, ishift, ic, i, ln
    character(len=llength) :: ctemp
    character(len=1)   :: punct(npunct)
    logical             :: keep_going

  ! Define the punctuation to be removed
    data punct          / '.', '-' /

  ! Initialise system dependant bounds in collating sequence

    ia     = ichar('A')
    iz     = ichar('Z')
    ishift = ichar('a') - ia

  ! Initialise output

    ln = len(string)

    if (ln < llength) then
       esdf_reduce(1:ln)  = string(1:ln)
       esdf_reduce(ln+1:) = ' '
    else
       esdf_reduce(1:llength) = string(1:llength)
    endif

  ! Drop all upper case characters to lower case

    do i = 1, llength
       ic = ichar(esdf_reduce(i:i))
       if ((ic >= ia) .and. (ic <= iz)) &
            esdf_reduce(i:i) = char(ishift + ic)
    enddo

  ! Now remove punctuation

    do i = 1, npunct
       keep_going = .true.
       do while (keep_going)
          ic = index(esdf_reduce, punct(i))
          if (ic > 0) then
             ctemp = esdf_reduce
             esdf_reduce(ic:) = ctemp(ic+1:)
          else
             keep_going = .false.
          endif
       enddo
    enddo

    esdf_reduce = trim (adjustl(esdf_reduce))

    return
  end function esdf_reduce

  !   --------------  esdf_split  ---------------------- 

  ! Split a string according to commas or semicolons (Ofer Sinai)

  function esdf_split (string, nmaxwords, nwords)

    character(len=*), intent(in) :: string
    integer, intent(in) :: nmaxwords
    integer, intent(out) :: nwords
    character(len=llength) :: esdf_split(nmaxwords)

  ! Local
    integer(i4b), parameter :: npunct = 2
    integer(i4b)       :: ic, ic_new, ip, i
    character(len=llength) :: ctemp
    character(len=1)   :: punct(npunct)
    logical             :: keep_going

  ! Define the punctuation to split by (OS: ';' does not work, all after it is
  ! ignored in practice)
    data punct          / ',', ';' /

  ! Split into separate words, store data in flag_words
    ctemp = string
    nwords = 0
    keep_going = .true.
    do while ((keep_going) .and. (nwords < nmaxwords))
       nwords = nwords+1
       ic = index(ctemp, punct(1))
       do ip = 1, npunct
          ic_new = index(ctemp, punct(ip))
          if ((ic_new > 0) .and. (ic_new < ic)) ic = ic_new
       enddo
       if (ic > 0) then
          esdf_split(nwords) = trim(adjustl(ctemp(1:ic-1)))
          ctemp = ctemp(ic+1:)
       else
          esdf_split(nwords) = trim(adjustl(ctemp))
          keep_going = .false.
       endif
    enddo

    return
  end function esdf_split

  !   --------------  esdf_convfac  ---------------------- 

  ! Find the conversion factor between physical units

  function esdf_convfac (from, to)

    character(len=*), intent(in) :: from, to
    real(dp)          :: esdf_convfac
  ! Local
    integer(i4b)       :: i, ifrom, ito
    character(len=llength) :: ctemp

  ! Find the index numbers of the from and to units

    ifrom = 0
    ito   = 0
    do i = 1, nphys
       if (esdf_reduce(from) == esdf_reduce(phy_n(i))) ifrom = i
       if (esdf_reduce(to) == esdf_reduce(phy_n(i))) ito = i
    enddo

  ! Check that the units were recognised

    if (ifrom == 0) then
       ctemp = 'Units not recognised in input file : '// &
            trim (esdf_reduce(from))
       call esdf_die (ctemp)
    endif

    if (ito == 0) then
       ctemp = 'Units not recognised in Program : '// &
            trim (esdf_reduce(to))
       call esdf_die (ctemp)
    endif

  ! Check that from and to are of the same dimensions

    if (phy_d(ifrom) /= phy_d(ito)) then
       ctemp = 'Dimensions Do not match : '// &
            trim (esdf_reduce(from))//' vs '// trim (esdf_reduce(to))
       call esdf_die (ctemp)
    endif

  ! Set the conversion factor

    esdf_convfac = phy_u(ifrom) / phy_u(ito)

    return
  end function esdf_convfac

  !   --------------  esdf_unit  ---------------------- 

  ! Find an unused i/o unit

  function esdf_unit (ierr)

    integer(i4b), intent(out) :: ierr
    integer(i4b)      :: esdf_unit
  ! Local
    logical            :: op

    ierr = 0
    do esdf_unit = 10, 99
       inquire (unit = esdf_unit, opened = op, err = 10)
       if (.not. op) return
    enddo
    call esdf_warn ('Unable to find a free i/o unit using esdf_u'// 'nit')
    ierr = 1

    return

10  continue
    call esdf_die ('Error opening files by esdf_unit')

    return
  end function esdf_unit

  !   --------------  esdf_file  ---------------------- 

  ! Open an old file

  subroutine esdf_file (unit, filename, ierr)

    character(len=*), intent(in) :: filename
    integer(i4b), intent(out)    :: unit, ierr
    logical            :: ex

    unit = esdf_unit(ierr)
    if (ierr > 0) return
    inquire (file = trim (filename), exist = ex, err = 10)
    if (.not. ex) goto 10
    open (unit = unit, file = trim (filename), form = 'formatted', &
         status = 'old', err = 10)

    return

10  continue
    ierr = 1

    return
  end subroutine esdf_file

  subroutine esdf_lablchk (string, typ,index)
      
    character(len=*), intent(in) :: string
    character(len=1), intent(in) :: typ
  ! Local
    character(len=llength) :: ctemp
    character(len=1)   :: tp
    integer(i4b)       :: i,index

  ! Check if label is recognised
      
    i     = count(esdf_reduce(string) == kw_label)
    ctemp = 'Label "'//trim (esdf_reduce(string))// &
         '" not recognised in keyword list'
    if (i == 0) call esdf_die (ctemp)
    ctemp = 'Label "'//trim (esdf_reduce(string))// &
         '" is multiply defined'
    if (i > 1) call esdf_die (ctemp)
    ctemp = 'Label "'//trim (esdf_reduce(string))// &
         '" has been used with the wrong type'
    tp    = ' '
    i     = 0
    do while (tp == ' ')
       i = i + 1
       if (esdf_reduce(string) == kw_label(i)) tp = kw_typ(i)
    enddo
    index = i
    if (typ /= tp) call esdf_die (ctemp)
      
    return
  end subroutine esdf_lablchk
      
  !   --------------  esdf_die  ---------------------- 

  ! Stop execution due to an error cause by esdf

  subroutine esdf_die (string)

    character(len=*), intent(in) :: string

    write (*, '(a,a)') ' ESDF error: ', trim (string)
    write (*, '(a)') ' Stopping now'

    stop

    return
  end subroutine esdf_die

  !   --------------  esdf_warn  ---------------------- 

  ! Warning due to an error cause by esdf

  subroutine esdf_warn (string)

    character(len=*), intent(in) :: string

    nwarns = nwarns + 1
    warns(nwarns) = string

    return
  end subroutine esdf_warn

  !   --------------  esdf_close  ---------------------- 

  ! Deallocate the data arrays --- call this before re-initialising

  subroutine esdf_close
    deallocate(llist,tlist,block_data)
    if (allocated(warns)) deallocate(warns)
  end subroutine esdf_close

end module esdf
!  ===============================================================
