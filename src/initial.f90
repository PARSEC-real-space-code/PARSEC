!===============================================================
!
! Copyright (C) 2005 Finite Difference Research Group
! This file is part of parsec, http://www.ices.utexas.edu/parsec/
!
! Initializing various variables, arrays, and files.
!
!---------------------------------------------------------------
subroutine initial(clust,elec_st,pbc,mol_dynamic,move, &
     elval,npolflg,field,ifield,istart,OutEvFlag,domain_shape, &
     ierr)

  use constants
  use cluster_module
  use electronic_struct_module
  use pbc_module
  use molecular_dynamic_module
  use movement_module
  implicit none
  !
  ! Input/Output variables:
  !
  ! the cluster
  type (cluster), intent(inout) :: clust
  ! electronic structure
  type (electronic_struct), intent(inout) :: elec_st
  type (pbc_data), intent(in) :: pbc
  ! molecular dynamic related data
  type (molecular_dynamic), intent(inout) :: mol_dynamic
  ! movement data
  type (movement), intent(inout) :: move
  ! number of valence electrons in system if neutral
  real(dp), intent(in) :: elval
  ! polarizability flag
  integer, intent(in) :: npolflg
  ! polarizability field (in a.u.)
  real(dp), intent(in) :: field
  ! polarizability field index (direction of field)
  integer, intent(out) :: ifield
  ! restart parameter 
  integer, intent(in) :: istart
  ! output flag
  integer, intent(in) :: OutEvFlag
  ! the shape of the domain, for cluster (non-periodic) calculations
  integer, intent(in) :: domain_shape
  ! error flag, 230 < ierr < 241
  integer, intent(out) :: ierr
  !
  ! Work variables:
  !
  ! size of temperature step in K, calculated according to the
  ! cooling scheme
  real(dp) :: tscale
  ! atomic coordinates (changed only if cluster is not centered
  ! about the origin)
  real(dp), dimension(clust%atom_num) :: xatm, yatm, zatm
  ! number of atoms for each atom type
  integer natmi(clust%type_num)
  ! LMD parameters:
  ! initial and final temperatures (in kelvin), 
  ! step temperature for stair cooling
  real(dp) :: tempi,tempf,tstep
  ! total actual number of atoms (from all types combined)
  integer natom
  ! counters
  integer ity, i, j, l, m, ii, isp
  ! number of electrons calculated from pre-fixed occupation file, if used
  real(dp) :: elsum
  ! temporary storage for minimum, maximum of coordinates
  real(dp) :: xmin,ymin,zmin,xmax,ymax,zmax
  ! amount by which to shift the cluster in order to center it
  real(dp) :: xmov,ymov,zmov
  ! actual number of spins used in computation
  integer nspin
  ! user-provided occupation of each eigenvalue
  real(dp) ::  occup(elec_st%nstate,elec_st%nspin)
  ! atom names
  character(len=2) name(1:clust%type_num)
  ! numerically negligible difference
  real(dp), parameter :: small = 1.0d-5
  ! holding array for clust%mvat as read from mdinit.dat
  real(dp) :: amove(clust%atom_num)

  !---------------------------------------------------------------

  natom = clust%atom_num
  xatm = clust%xatm
  yatm = clust%yatm
  zatm = clust%zatm
  natmi = clust%natmi
  name = clust%name

  nspin = elec_st%nspin
  tempi = mol_dynamic%tempi
  tempf = mol_dynamic%tempf
  tstep = mol_dynamic%step_tmp
  !
  ! If restart of molecular dynamics then we need to read this file
  ! first. We get the coordinates for our run from this file.
  ! read in restart data. We also read the last time step of the 
  ! previous run (iframe).
  !
  if ((mol_dynamic%is_on) .and. (mol_dynamic%is_restart)) then
     open(76,file='mdinit.dat',status='old',form='formatted',iostat=ii)
     if (ii /= 0) then
        write(7,*) 'ERROR! restart file mdinit.dat not found.'
        write(7,*) 'STOP in initial.'
        ierr = 234
        return
     endif
     read(76,20) mol_dynamic%iframe
     do j = 1,natom
        read(76,22) mol_dynamic%xcur(j),mol_dynamic%ycur(j), &
             mol_dynamic%zcur(j),mol_dynamic%xold(j), &
             mol_dynamic%yold(j),mol_dynamic%zold(j), &
             mol_dynamic%vxold(j),mol_dynamic%vyold(j), &
             mol_dynamic%vzold(j),mol_dynamic%accxold(j), &
             mol_dynamic%accyold(j),mol_dynamic%acczold(j), &
             amove(j)
     enddo
     close(76)
     ! Assign xatm, yatm and zatm to be xcur, ycur and zcur 
     ! They are moved to their rightful places right after this.
     do j = 1, natom
        xatm(j) = mol_dynamic%xcur(j)  
        yatm(j) = mol_dynamic%ycur(j)  
        zatm(j) = mol_dynamic%zcur(j)  
        clust%mvat(j) = int(amove(j))
     enddo
  endif
20 format('Step #',i4)

  if ((.not. pbc%is_on) .and. &
       (istart == 0) .and. (.not.(move%name == MANUAL))) then
     ! If confined system:
     ! if the coordinates are not symmetrized about the origin, shift
     ! them so that they are, in order to allow the results to
     ! converge already for a smaller Rmax.
     !
     ! find minimum and maximum coordinate along each axis
     xmin=xatm(1)
     xmax=xatm(1)
     ymin=yatm(1)
     ymax=yatm(1)
     zmin=zatm(1)
     zmax=zatm(1)

     do j = 2, natom
        if (xmin > xatm(j)) xmin = xatm(j)
        if (xmax < xatm(j)) xmax = xatm(j)
        if (ymin > yatm(j)) ymin = yatm(j)
        if (ymax < yatm(j)) ymax = yatm(j)
        if (zmin > zatm(j)) zmin = zatm(j)
        if (zmax < zatm(j)) zmax = zatm(j)
     enddo
 
     ! Move the coordinates so that along each axis the farthest atom
     ! in the positive direction is as far as the farthest atom in
     ! the negative direction.
     xmov = half*(xmin+xmax)
     ymov = half*(ymin+ymax)
     zmov = half*(zmin+zmax)

     if (sqrt(xmov*xmov+ymov*ymov+zmov*zmov) > 0.1) then
        write(7,*)
        write(7,*) 'WARNING: atom coordinates were shifted so'
        write(7,*) 'that the cluster is centered at the origin'
        if((.not. pbc%is_on) .and. domain_shape > 0) then
           ! note, this recentering should not cause problems for spherical
           ! domains
           write(7,*)
           write(7,*) 'WARNING: This may cause atoms to be outside the domain.'
           write(7,*) 'Coordinates in the input file may need to be adjusted.'
        endif
        write(7,*)
        write(7,*) 'The new atom coordinates are:'

        do j=1,natom
           xatm(j) = xatm(j)-xmov
           yatm(j) = yatm(j)-ymov
           zatm(j) = zatm(j)-zmov
           write(7,18) xatm(j),yatm(j),zatm(j)
        enddo
        if (clust%has_ptchrg) then
           write(7,*)
           write(7,*)'WARNING: point charge coordinates were shifted'
           write(7,*)'in accordance with shift of atom coordinates!'
           write(7,*)'The new point charge coordinates are:'

           do j=1, clust%nptchrg
              clust%xpt(j) = clust%xpt(j)-xmov
              clust%ypt(j) = clust%ypt(j)-ymov
              clust%zpt(j) = clust%zpt(j)-zmov
              write(7,18) clust%xpt(j),clust%ypt(j),clust%zpt(j)
           enddo
        endif
18      format(3(f15.9,3x))
     endif
  elseif (pbc%is_on) then

     if ((pbc%per==2) .and. &
       (istart == 0) .and. (.not.(move%name == MANUAL))) then

       zmin=zatm(1)
       zmax=zatm(1)

       do j = 2, natom
        if (zmin > zatm(j)) zmin = zatm(j)
        if (zmax < zatm(j)) zmax = zatm(j)
       enddo

     ! Move the coordinates so that along z axis the farthest atom
     ! in the positive direction is as far as the farthest atom in
     ! the negative direction.
       
       zmov = half*(zmin+zmax)
       if (abs(zmov) > 0.1) then
        write(7,*)
        write(7,*) 'WARNING: atom coordinates were shifted so'
        write(7,*) 'that the slab is centered at the origin'
        write(7,*) 'The new atom coordinates are:'

        do j=1,natom
           zatm(j) = zatm(j)-zmov
           write(7,18) xatm(j),yatm(j),zatm(j)
        enddo
       endif ! abs(zmov)>0.1
     endif ! pbc%per==2
       
     ! If periodic system: make sure all atoms are inside the unit cell
!  AMIR - have to check this change with Doron !!!!!!
!     call back_to_cell(pbc,natom,xatm,yatm,zatm)

  endif

  ! compute actual number of electrons from # of valence electrons
  ! and cluster charge.
  elec_st%xele = elval - elec_st%ncharge
  !
  ! If polarizability flag is on -
  ! set ifield, field direction index, to 1 (initial direction)
  ! open polar.dat file (this is an ascii file!), into which 
  ! polarizability results are written. Write in some headers for
  ! the information will be fed later on in the calculation.
  !
  if (npolflg == 1) then
     ifield = 1
     open(91,file='polar.dat', status='unknown',form='formatted')
     write(91,*) 'Tot. number of atoms is: ', natom
     write(91,*) 'Number of different types of atoms is: ',clust%type_num 
     l = 0

     do i = 1, clust%type_num
        write(91,*)
        write(91,15) natmi(i),name(i)
        do j = 1, natmi(i)
           l = l + 1
           write(91,92) xatm(l), yatm(l), zatm(l)
        enddo
     enddo
92   format(3x, f10.6, 3x, f10.6, 3x, f10.6)
15   format(i3,1x,a2,1x,'atoms with coordinates:')
     write(91,*)
     write(91,85) field
85   format(1x,'Applied Electric field  = ',f8.5,1x,'[Ry/au]')
     write(91,*)
     write(91,86) 
     write(91,87)
86   format(3x,'Field',11x,'Tot.',14x,'Dipole  Moment [Deb.]')
87   format(1x,'Component',6x,'Energy [Ry]',10x,'X',10x,'Y',10x,'Z')
  end if
  !
  ! If atoms move, open atom coordinate file and write out initial
  ! coord. Initialize 'is minimum force satisfied' flag to zero.
  !
  if (move%is_on) then
     open(66,file='atom.cor',status='unknown',form='formatted')
     write(66,12) move%num
     j=0

     do ity = 1, clust%type_num
        write(66,*) name(ity),'    ',natmi(ity)
        do i=1,natmi(ity)
           j = j + 1
           write(66,14) xatm(j), yatm(j), zatm(j)
        enddo
     enddo
     call myflush(66)
     !
     ! If minimization is being done restarting from previous run,
     ! we should keep the user-provided geometry otherwise that will
     ! not be consistent with information from previous runs.
     !
     if (move%is_restart) then
        open(45,file='relax_restart.dat', &
             form='formatted',status='old',iostat=ii)
        if (ii /= 0) then
           write(7,*) 'WARNING! relax_restart.dat not found'
           write(7,*) 'start structural relaxation from scratch'
           move%is_restart = .false.
        else
           write(7,*) 'WARNING: atom coordinates in input were'
           write(7,*) 'replaced by coordinates in relax_restart.dat'
           write(7,*) 'The new atom coordinates are:'
           read(45,*)
           do j=1,natom
              read(45,*) xatm(j),yatm(j),zatm(j)
              write(7,18) xatm(j),yatm(j),zatm(j)
           enddo
           write(7,*)
           close(45)
        endif
     endif
     ! if manual mode open the file from which coordinates are read
     if (move%name == MANUAL) then
        open(67,file='manual.dat',status='old',form='formatted',iostat=ii)
        if (ii /= 0) then
           write(7,*) 'WARNING! manual.dat not found.'
           if (move%num > 0) then
              write(7,*) 'ERROR! Movement of atoms is not specified.'
              write(7,*) 'Missing manual.dat. STOP in initial.'
              ierr = 235
              return
           endif
        endif
     endif
     ! if BFGS minimization, open file for output of BFGS routines
     if (move%name == BFGS) then
        open(46,file='bfgs.dat',status='unknown',form='formatted')
     endif
  endif
12 format('Step #',1x,i3)
14 format(3(3x,f10.6))
  !
  ! If molecular dynamics - open files.
  ! Also write initial statements into these files,initialize
  ! counters, and initialize temperature step according to cooling
  ! scheme.
  !
  if (mol_dynamic%is_on) then
     !
     ! initialize frame counter and open atom coordinate file
     ! if we are restarting the molecular dynamics iframe is read
     ! in from the restart file (mdinit.dat)
     !
     if (.not. mol_dynamic%is_restart) mol_dynamic%iframe = 0 
     open(66,file='atom.cor',status='unknown',form='formatted')
     ! open file containing energy data from MD run
     open(78,file='md_nrg.dat',status='unknown',form='formatted')
     write(78,88)
     ! open file containing mechanical data from MD run
     open(79,file='md_mech.dat',status='unknown',form='formatted')
     write(79,81)
     write(79,82)
     write(79,83)
81   format('Each line contains: six coordinates, three', &
          ' velocities, and three accelerations',/,'for each', &
          ' molecular dynamics step')
82   format('The data are: next position, current position', &
          'current velocity',/,'current acceleration')
83   format('The last step of this file should be copied to', &
          ' mdinit.dat for a restart')
     ! Open file containing coordinates for each step. "movie.xyz"
     ! has coordinates for each step, in angstrom and in .xyz format
     ! (compatible with xmol and other molecular visualization codes).
     open(92,file='movie.xyz',status='unknown',form='formatted')
     write(92,90) natom
     write(92,*) 'Step # ', mol_dynamic%iframe
     ii = 0

     do i = 1, clust%type_num
        do j = 1, natmi(i)
           ii = ii + 1
           write(92,93) name(i), xatm(ii)*angs,yatm(ii)*angs, zatm(ii)*angs
        enddo
     enddo

     select case (mol_dynamic%cool_type)
        ! "Stairstep" cooling
     case (MD_STAIR)
        tscale = -tstep
        ! error checks: moving in right direction?
        if (((tempi > tempf) .and. (tscale >= 0.d0)) .or. &
             ((tempi < tempf) .and. (tscale <= 0.d0))) then
           write(7,*) 'ERROR: temperature step of wrong sign!'
           write(7,*) 'STOP in initial'
           ierr = 231
           return
        endif
        ! unless it is an isothermal run (in which case the temperature
        ! step is irrelevant), check if the temperature step is not
        ! too large. "Too large" means that it is more than twice the
        ! temperature gap (because then we won't be able to clamp it to
        ! the final temperature - see the temperature clamping in the 
        ! molecular dynamics subroutine).
        if ((abs(tempi-tempf) > 1.d-5) .and. &
             (abs(tscale/2.d0) >= abs(tempi-tempf))) then
           write(7,*) 'ERROR: temperature step too large!'
           write(7,*) 'it goes far beyond the final temperature!'
           write(7,*) 'STOP in initial'
           ierr = 232
           return
        endif
        ! linear cooling - substract a fixed temperature step
     case (MD_LINEAR)
        tscale = (tempf-tempi)/real(mol_dynamic%step_num,dp)
        ! logarithm cooling - multiply the temperature by a fixed 
        ! ratio at each step
     case(MD_LOG)
        tscale = (tempf/tempi)**(one/real(mol_dynamic%step_num,dp))
     end select
     mol_dynamic%tscale = tscale

     ! if restart of molecular dynamics
     if (mol_dynamic%is_restart) then

        ! Restart data has already been read. Just do some reporting.
        ! report initial coordinates to atom.cor.
        write(66,*) 'Coordinates for step #', mol_dynamic%iframe
        ii = 0
        do i = 1, clust%type_num
           write(66,*) name(i)
           do j = 1, natmi(i)
              ii = ii + 1
              write(66,89) xatm(ii),yatm(ii),zatm(ii)
           enddo
        enddo

        ! Report initial position,velicity, acceleration to md_mech.dat.
        write(79,*) 'Step #', mol_dynamic%iframe
        do j=1,natom
           write(79,22) mol_dynamic%xcur(j),mol_dynamic%ycur(j), &
                mol_dynamic%zcur(j),mol_dynamic%xold(j), &
                mol_dynamic%yold(j),mol_dynamic%zold(j), &
                mol_dynamic%vxold(j),mol_dynamic%vyold(j), &
                mol_dynamic%vzold(j),mol_dynamic%accxold(j), &
                mol_dynamic%accyold(j),mol_dynamic%acczold(j), &
                amove(j)
        enddo
     else
        ! if we are not restarting then we need to make xcur etc equal
        ! the xatm etc.
        do j = 1, natom
           mol_dynamic%xcur(j) = xatm(j)
           mol_dynamic%ycur(j) = yatm(j)
           mol_dynamic%zcur(j) = zatm(j)
        enddo
     endif
  end if

22 format(13(1x,e15.8))
88 format('Step',2x,'Time[ps]',2x,'Ekin[eV/at]',2x,'Epot[eV/at]', &
        2x,'Etot[eV/at]',5x,'T[K]',2x,'Tcalc[K]')
89 format(3(2x,f11.6))
90 format(i3)
93 format(a2, 3(3x, f10.6))
  !
  ! In the case of artificial occupations, denoted by a negative
  ! Fermi temperature, read occupations from occup.in.
  !
  if (elec_st%tfermi < zero) then
     open(90,file='occup.in',status='old',form='formatted',iostat=ii)
     if (ii /= 0) then
        write(7,*) 'ERROR! File occup.in not found.'
        write(7,*) 'STOP in initial.'
        ierr = 236
        return
     endif
     do isp = 1, nspin/elec_st%mxwd
        ! skip title line
        read(90,*) 
        ! read state up to which the states are full
        read(90,*) elec_st%ifmax(isp)
        ! occupy states
        do i = 1, elec_st%ifmax(isp)
           occup(i,isp)=one
        enddo
        do i = elec_st%ifmax(isp) + 1, elec_st%nstate
           occup(i,isp)=zero
        enddo
        ! read in manual over-ride for occupation of specific states
        read(90,*) m
        do i = 1, m
           read(90,*) j, occup(j,isp)
        enddo
     enddo
     close(90)

     ! check if occupation consistent with number of electrons
     elsum = zero
     do isp = 1, nspin/elec_st%mxwd
        do i = 1, elec_st%nstate
           elsum = elsum +occup(i,isp)
        enddo
     enddo
     !
     ! 3-nspin is 2 if nspin=1 and 1 if nspin=2. This is because if
     ! there's only one spin channel, our occupation count is defined
     ! as 1.0 for TWO electrons filling the state, whereas if there
     ! are two spin channels our occupation count is defined as 1.0
     ! for ONE electron filling the state.
     !
     if (abs(real(3-nspin,dp) * elsum-elec_st%xele) > small) then
        write(7,*)
        write(7,*) 'ERROR: fixed occupation inconsistent'
        write(7,*) '       with stated number of electrons!' 
        write(7,*) 'STOP in initial'
        ierr = 233
        return
     endif

     ! Find highest eigenvalue with non-zero occupation
     do isp=1, nspin/elec_st%mxwd
        do i=elec_st%nstate,1,-1
           if (occup(i,isp) /= zero) then
              elec_st%ifmax(isp) = i
              goto 25
           endif
        enddo
25      continue
     enddo

     ! end fixed occupation read-in
  endif

  clust%xatm = xatm
  clust%yatm = yatm
  clust%zatm = zatm

  elec_st%occ_in = occup
  ! output flag
  if (OutEvFlag > 0) then
     open(61,file='eigen.dat',status='unknown',form='formatted')
     if (elec_st%nkpt == 0) then
        write(61,*) ' State     Eigenvalue [Ry]    Eigenvalue [eV]', &
             '    Occupation   Repr.   Movement   spin'
     else
        write(61,*) ' State     Eigenvalue [Ry]    Eigenvalue [eV]', &
             '    Occupation   Repr.   k-point  Movement   spin'
     endif
  endif
  ! end  output flag read-in

end subroutine initial
!===============================================================
!
! Subroutine to set the various clm coefficients. Separated so as
! to not clutter up the above code. calculates clm = (l-m)!/(l+m)! 
! for l = 0,lpole, m=0,l. 
! Evaluating these to 20 decimal places is certainly as accurate 
! as doing the calculations in the code but less work.
!
!---------------------------------------------------------------
subroutine setclm(lpole, clm)
  use constants
  implicit none
  !
  ! Input/Output variables:
  !
  ! order of mutipole expansion
  integer, intent(in) :: lpole
  ! clm coefficients
  real(dp), intent(out) :: clm(0:lpole, 0:lpole)
  !---------------------------------------------------------------
  clm(0,0) = 1.00000000000000000000e+00
  if (lpole >= 1) then
     clm(1,0) = 1.00000000000000000000e+00
     clm(1,1) = 5.00000000000000000000e-01
  endif
  if (lpole >= 2) then
     clm(2,0) = 1.00000000000000000000e+00
     clm(2,1) = 1.66666666666666666670e-01
     clm(2,2) = 4.16666666666666666670e-02
  endif
  if (lpole >= 3) then
     clm(3,0) = 1.00000000000000000000e+00
     clm(3,1) = 8.33333333333333333330e-02
     clm(3,2) = 8.33333333333333333330e-03
     clm(3,3) = 1.38888888888888888890e-03
  endif
  if (lpole >= 4) then
     clm(4,0) = 1.00000000000000000000e+00
     clm(4,1) = 5.00000000000000000000e-02
     clm(4,2) = 2.77777777777777777780e-03
     clm(4,3) = 1.98412698412698412700e-04
     clm(4,4) = 2.48015873015873015870e-05
  endif
  if (lpole >= 5) then
     clm(5,0) = 1.00000000000000000000e+00
     clm(5,1) = 3.33333333333333333330e-02
     clm(5,2) = 1.19047619047619047620e-03
     clm(5,3) = 4.96031746031746031750e-05
     clm(5,4) = 2.75573192239858906530e-06
     clm(5,5) = 2.75573192239858906530e-07
  endif
  if (lpole >= 6) then
     clm(6,0) = 1.00000000000000000000e+00
     clm(6,1) = 2.38095238095238095240e-02
     clm(6,2) = 5.95238095238095238100e-04
     clm(6,3) = 1.65343915343915343920e-05
     clm(6,4) = 5.51146384479717813050e-07
     clm(6,5) = 2.50521083854417187750e-08
     clm(6,6) = 2.08767569878680989790e-09
  endif
  if (lpole >= 7) then
     clm(7,0) = 1.00000000000000000000e+00
     clm(7,1) = 1.78571428571428571430e-02
     clm(7,2) = 3.30687830687830687830e-04
     clm(7,3) = 6.61375661375661375660e-06
     clm(7,4) = 1.50312650312650312650e-07
     clm(7,5) = 4.17535139757361979580e-09
     clm(7,6) = 1.60590438368216145990e-10
     clm(7,7) = 1.14707455977297247140e-11
  endif
  if (lpole >= 8) then
     clm(8,0) = 1.00000000000000000000e+00
     clm(8,1) = 1.38888888888888888890e-02
     clm(8,2) = 1.98412698412698412700e-04
     clm(8,3) = 3.00625300625300625300e-06
     clm(8,4) = 5.01042167708834375500e-08
     clm(8,5) = 9.63542630209296875960e-10
     clm(8,6) = 2.29414911954594494280e-11
     clm(8,7) = 7.64716373181981647590e-13
     clm(8,8) = 4.77947733238738529740e-14
  endif
  if (lpole >= 9) then
     clm(9,0) = 1.00000000000000000000e+00
     clm(9,1) = 1.11111111111111111110e-02
     clm(9,2) = 1.26262626262626262630e-04
     clm(9,3) = 1.50312650312650312650e-06
     clm(9,4) = 1.92708526041859375190e-08
     clm(9,5) = 2.75297894345513393130e-10
     clm(9,6) = 4.58829823909188988550e-12
     clm(9,7) = 9.55895466477477059490e-14
     clm(9,8) = 2.81145725434552076320e-15
     clm(9,9) = 1.56192069685862264620e-16
  endif

end subroutine setclm
!===============================================================
!
! During movement of atoms, it may happen that an atom falls
! outside the unit cell. If that happens, the atom must be
! dragged back. We do it by replacing the atom with its periodic
! image inside the unit cell. This is how it is done:
! 1. Atomic coordinates are first translated to lattice coordinates.
! 2. The integer part of the lattice coordinates is subtracted.
! 3. lattice coordinates are translated back to atomic coordinates.
!
!---------------------------------------------------------------
subroutine back_to_cell(pbc,natom,xatm,yatm,zatm)
  use constants
  use pbc_module
  implicit none
  !
  ! Input/Output variables:
  !
  type (pbc_data), intent(in) :: pbc
  ! number of atoms
  integer, intent(in) :: natom
  ! position of atoms, in Cartesian coordinates
  real(dp), intent(inout), dimension(natom) :: xatm, yatm, zatm
  !
  ! Work variables:
  !
  ! quantities used to push atoms back to the cell
  real(dp) :: tmp_xyz_old(3), tmp_uvw(3), tmp_xyz_new(3), tmp1, tmp2
  ! inverse unit lattice vectors
  real(dp) :: alinv(3,3)
  ! counters
  integer :: jj, ii

  !---------------------------------------------------------------

  alinv = pbc%latt_vec
  call mtrxin(alinv,tmp1,tmp2)

  do jj = 1, natom
     tmp_xyz_old(1) = xatm(jj)
     tmp_xyz_old(2) = yatm(jj)
     tmp_xyz_old(3) = zatm(jj)
     call matvec3('N',alinv,tmp_xyz_old,tmp_uvw)

     do ii = 1, pbc%per
        tmp_uvw(ii) = tmp_uvw(ii) - nint(tmp_uvw(ii))
     enddo
     call matvec3('N',pbc%latt_vec,tmp_uvw,tmp_xyz_new)
     if (pbc%per <= 2) tmp_xyz_new(3) = zatm(jj)
     if (pbc%per <= 1) tmp_xyz_new(2) = yatm(jj)

     xatm(jj) = tmp_xyz_new(1)
     yatm(jj) = tmp_xyz_new(2)
     zatm(jj) = tmp_xyz_new(3)
  enddo

end subroutine back_to_cell
!===============================================================
