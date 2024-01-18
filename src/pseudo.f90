!===============================================================
!
! Copyright (C) 2005 Finite Difference Research Group
! This file is part of parsec, http://www.ices.utexas.edu/parsec/
!
! Sets up the pseudopotentials and computes associated quantities.
!
!---------------------------------------------------------------
subroutine pseudo(clust,p_pot,nstate,elval1,icorr,ierr)

  use constants
  use cluster_module
  use pseudo_potential_module
  implicit none
  !
  !  Input/Output variables:
  !
  type (cluster), intent(inout) :: clust
  type (pseudo_potential), intent(inout) :: p_pot

  !  User specified number of states
  integer, intent(in) :: nstate

  !  # of valence electrons for a neutral system
  real(dp), intent(out) :: elval1

  !  correlation type
  character(len=2), intent(in) :: icorr

  !  error flag, 200 < ierr < 231
  integer, intent(out) :: ierr
  !
  !  Work variables:
  !
  !  order for the pseudopotential derivative (number of points
  !  on each side used); maximum value is mor = 10
  integer, parameter :: mor = 9

  !  number of different atom types
  integer naty

  !  number of components in the pseudopotential (per given atom)
  !  read from pseudopotential file
  integer lll

  !  maximum number of grid points in given pseudopotential files
  !  maxpot = max (p_pot%ns)
  integer maxpot

  !  number of spin-polarized pseudopotential components
  integer nlso

  !  counters, temporary variables
  integer i, iorb, j, ity, i0, j0, kk, k, jrel, isy
  integer itmp(10)
  integer npoint
  real(dp) :: a, b, c, temp1, temp0, fspd, delr, rtmp(10)
  integer, allocatable :: nloc(:)

  !  temporary core cutoff parameters
  real(dp) :: rc, tcore,cfac,rcfac

  !  # electrons in each pseudopotential component (user input)
  real(dp) :: elval2

  !  temporary storage of the atomic mass of each atom type
  real(dp) :: am

  !  list of coefficients for numerical first derivative
  real(dp) ::  coe(-mor:mor)

  !  <psi_l|psi_l> - check that the input pseudowavefunction is
  !  properly normalized on the radial grid
  real(dp) :: psinorm(4)

  !  vspd - pseudopotential (up to four components)
  real(dp), dimension (:,:,:), allocatable :: vspd

  !  temporary arrays to calculate the spin-orbit normalization factors
  real(dp), dimension (:,:), allocatable :: clj

  ! temporary arrays to calculate the spin-orbit projectors and forces
  real(dp), dimension (:,:,:), allocatable :: dvlj

  !  vlj - spin-orbit pseudopotential (p and d components only)
  real(dp), dimension (:,:,:), allocatable :: vlj

  !  tolerance for identifying ekbi as zero)
  real(dp), parameter :: small = 1.0d-16

  !  position of fixed point for the calculation of logarithmic grid
  !  parameters (par_a, par_b, par_c)
  integer, parameter :: i_0 = 101

  !  occupancy of orbitals (from pseudopotential file)
  real(dp) :: eleatm_psp(4)

  !  extension that turns atom name into the pseudopotential file name
  character(len=90), dimension(:), allocatable :: fname,fnames

  !  atom names
  character(len=2) :: name(clust%type_num)
  !  format of pseudopotential files
  integer :: psformt

  !  strings read from pseudopotential file
  character(len=2) :: nameat,icorrt,icorrt2
  character(len=3) :: irel
  character(len=4) :: nicore
  character(len=200) :: linestr

  !---------------------------------------------------------------
  naty = p_pot%type_num
  name = clust%name

  write(7,'(/,/,a)') ' Pseudopotential messages'
  write(7,'(a,/)') ' ------------------------'
  !  initialize counters
  !  for # of valence electrons
  elval1 = zero
  !  for number of atoms
  k=1

  allocate(fname(naty))
  allocate(fnames(naty))
  !  go over all atom types
  do ity = 1, naty  

     !  Get atomic mass from periodic table.
     call ptable(name(ity),am,ierr)
     psformt = p_pot%format(ity)
     if (psformt == MARTINS_OLD) then
        elval2 = sum(p_pot%eleatm(ity,:))
        p_pot%zion(ity) = elval2
     endif
     if (ierr /= 0) return
     !  spread atomic mass value to all atoms of this type
     do kk=k,k+clust%natmi(ity)-1
        clust%amass(kk) = am*1836.d0
     enddo
     k = k + clust%natmi(ity)
     !  open potential & wave function input files 
     select case (psformt)
     case (MARTINS_OLD)
        fname(ity) = trim(name(ity))//'_POTR.DAT'
     case (MARTINS_NEW)
        fname(ity) = trim(name(ity))//'_POTRE.DAT'
     case (MARTINS_WANG)
        fname(ity) = trim(name(ity))//'_POTRW.DAT'
     case (FHIPP)
        fname(ity) = trim(name(ity))//'_FHIPP.DAT'
     end select

     open(unit=20,file=trim(fname(ity)), &
          status='old',form='formatted', err=50)
       !
       !  Prepare to read in  Pseudopotential data
       !  Define these variables: p_pot%ns, p_pot%icore, p_pot%nlocp
       !
     select case(psformt)

        !  Old Martins format
     case (MARTINS_OLD)

        !  pseudopotential parameters read from first line of
        !  pseudopotential file:
        !  (icore = 0 - no core, = 1 - core correction)
        read(20,*) p_pot%ns(ity), lll, p_pot%icore(ity)

        if (lll /= p_pot%nlocp(ity)) then
           write(7,*) 'ERROR: # of pseudopot. components from user'
           write(7,*) 'incompatible with # from pseudopot. file'
           write(7,*) ' atom type ',ity,lll,p_pot%nlocp(ity)
           write(7,*) 'STOP in pseudo '
           ierr = 202
           return
        endif

        !  New Martins format
     case (MARTINS_NEW)

        read(20,'(1x,a2,1x,a2,1x,a3,1x,a4)') nameat, icorrt, irel, nicore
        if (trim(nicore) == 'nc') then
           p_pot%icore(ity) = 0
        else
           p_pot%icore(ity) = 1
        endif
        read(20,'(a200)',end=68) linestr
68      continue
        read(20,'(a200)',end=69) linestr
69      continue
        read(20,*) p_pot%nlocp(ity), i, p_pot%ns(ity)

        !  Martins-Wang format
     case (MARTINS_WANG)

        read(20,*) p_pot%ns(ity), p_pot%icore(ity),i,temp0,p_pot%loc(ity)
        read(20,*) i, j, k
        p_pot%nlocp(ity) = 3
        if (i == 0) p_pot%nlocp(ity) = p_pot%nlocp(ity) - 1
        if (j == 0) p_pot%nlocp(ity) = p_pot%nlocp(ity) - 1
        if (k == 0) p_pot%nlocp(ity) = p_pot%nlocp(ity) - 1

        !  FHIPP format
     case (FHIPP)

        read(20,*) temp0, lll
        if (lll /= p_pot%nlocp(ity)) then
           write(7,*) 'ERROR: # of pseudopot. components from user'
           write(7,*) 'incompatible with # from pseudopot. file'
           write(7,*) ' atom type ',ity,lll,p_pot%nlocp(ity)
           write(7,*) 'STOP in pseudo '
           ierr = 203
           return
        endif

        do kk = 1, 10
           read(20,*)
        enddo
        read(20,*) p_pot%ns(ity)

     end select

     if (p_pot%loc(ity) > p_pot%nlocp(ity) .or. p_pot%loc(ity) < 0 ) then
        write(7,*)
        write(7,*) 'ERROR: problem with atom type', ity
        write(7,*) 'local pseudopot. beyond # of components: '  &
             ,p_pot%loc(ity),p_pot%nlocp(ity)
        write(7,*) 'STOP in pseudo'
        ierr = 204
        return
     endif
     if (p_pot%nlocp(ity) > 4) then
        write(7,*)'ERROR: More than 4 pseudopot. components'
        write(7,*) 'Do you really need to go beyond f orbitals?'
        write(7,*) 'STOP in pseudo'
        ierr = 205
        return
     endif

     !  the increase by one is in order to "push in" the point at r=0 -
     !  see below
     p_pot%ns(ity) = p_pot%ns(ity)+1
     close(20)
  enddo                     ! ity = 1, naty (first loop)

  !  getting the maximum number of lines in the given
  !  pseudopotential files
  maxpot = maxval(p_pot%ns)
  ! why was this reported as uninialized memory access??
  !from set_mxpot ->allocate (p_pot%vw (-norder:mxpot,type_num,4))
  call pseudo_potential_set_mxpot (maxpot,p_pot)
  isy = 0
  do ity = 1, naty
      if (p_pot%so(ity)) isy = isy + 1
  enddo
  !  allocate working arrays
  allocate (vspd (-p_pot%norder:maxpot, p_pot%type_num, 4))
  vspd(:,:,:) = zero

  allocate (dvlj(-p_pot%norder:maxpot,isy,4))
  dvlj(:,:,:) = zero
  p_pot%vionr(:,:,:) = zero
  p_pot%vsor(:,:,:) = zero


  allocate (clj (4,isy))
  clj(:,:) = zero

  allocate (vlj (-p_pot%norder:maxpot,isy,4))
  vlj(:,:,:) = zero
  p_pot%wfspd(:,:,:) = zero
  p_pot%dwfspd(:,:,:) = zero
  elval2 = zero
  isy = 0
  p_pot%cc(:,:) = zero
  
  do ity = 1, naty

     jrel = 0
     psformt = p_pot%format(ity)
     open(unit=20,file=trim(fname(ity)),  &
          status='old',form='formatted', err=50)
     write(7,55) ity, trim(fname(ity))
55   format(1x,'Pseudopotential file for atom type',1x,i2,1x,'is',1x,a)

     select case(psformt)

     case (MARTINS_OLD)
        write(7,'(a,/)') 'Old Martins generator'
        read(20,*)

        !  first, read in s,p,d,f pseudopotentials
        !  then, read in s,p,d pseudo wave functions
        !  all reading starts from j=2 to leave room for r=0 at j=1
        allocate(nloc(p_pot%nlocp(ity)))
        read(20,*) (nloc(j),j=1,p_pot%nlocp(ity))
        do j = 2, p_pot%ns(ity)
           read(20,*) p_pot%rs(j,ity),(vspd(j,ity,nloc(i)),i=1  &
                ,p_pot%nlocp(ity))
        enddo
        do  j = 2,p_pot%ns(ity)
           read(20,*) p_pot%rs(j,ity),(p_pot%wfspd(j,ity,nloc(i))  &
                ,i=1,p_pot%nlocp(ity))
        enddo
        deallocate(nloc)

        if (p_pot%so(ity)) then
           write(7,*) 'ERROR!  spin-orbit coupling calculation work with new format pseudopotentials only'
           ierr = 86
!           isy = isy + 1
!           read(20,*) j0
!           read(20,*)clj(3,isy),clj(1,isy),clj(4,isy),clj(2,isy)
!           do j = 1,p_pot%ns(ity)-1
!              read(20,*)vlj(j,isy,3),vlj(j,isy,1),vlj(j,isy,4),vlj(j,isy,2)
!           enddo
        endif

        !
        !  Calculate the radius step parameters - the pseudopotential is
        !  given on a log grid !
        !          r(i) = par_a*exp(par_b*(i-1)) - par_c
        !  The calculation is based on the results at two grid points that
        !  are quite far away, for accuracy.
        !
        j0 = 2*(i_0-1) + 1
        if(j0 > p_pot%ns(ity)) then
           write(7,*) 'ERROR: choose smaller i_0 to get par_a & par_b, i_0 = ',i_0
           write(7,*) 'STOP in pseudo'
           ierr = 206
           return
        endif
        p_pot%par_a(ity) = p_pot%rs(i_0,ity)**2/  &
             (p_pot%rs(j0,ity)-two*p_pot%rs(i_0,ity))
        p_pot%par_b(ity) = log(p_pot%rs(j0,ity)/  &
             p_pot%par_a(ity) + one)/real(j0-1,dp)
        p_pot%par_c(ity) = p_pot%par_a(ity)
        !
        !  Read in the Partial-Core Charge Density, if present.
        !
        do j=2,p_pot%ns(ity)
           p_pot%denc(j,ity) = zero
        enddo
        if (p_pot%icore(ity) > 0) then
           do j=2,p_pot%ns(ity)
              read(20,*,end=80) p_pot%rs(j,ity), p_pot%denc(j,ity)
           enddo
80         continue
        endif
        !
        !  Read spin-orbit potential from file xx_SO.DAT
        !  This option has been removed. 
        !  One must use new frmat for SO calculation
        !
!        if (p_pot%so(ity)) then
!           write(7,*)'Spin-Orbit works with old Martins pseudopotential ', &
!                'format here'
!        endif

     case (MARTINS_NEW)

        write(7,'(a,/)') 'New Martins generator'
        read(20,'(1x,a2,1x,a2,1x,a3,1x,a4)') nameat, icorrt, irel, nicore
        write(7,82) nameat,icorrt,irel,nicore
82      format(/,a2,2x,a2,2x,a3,2x,a4)
        read(20,'(a200)',end=84) linestr
84      continue
        write(7,*) trim(linestr)
        read(20,'(a200)',end=86) linestr
86      continue
        write(7,'(a,/)') trim(linestr)

        if (trim(nicore) == 'nc') then
           read(20,*) p_pot%nlocp(ity), nlso, i,  &
                p_pot%par_a(ity),p_pot%par_b(ity), p_pot%zion(ity)
        else
           read(20,*) p_pot%nlocp(ity), nlso, i,  &
                p_pot%par_a(ity),p_pot%par_b(ity), p_pot%zion(ity), &
                cfac,rcfac
        endif
        p_pot%par_c(ity) = p_pot%par_a(ity)

        if (icorrt == 'ca') icorrt2 = 'CA'
        if (icorrt == 'xa') icorrt2 = 'XA'
        if (icorrt == 'wi') icorrt2 = 'WI'
        if (icorrt == 'hl') icorrt2 = 'HL'
        if (icorrt == 'pb') icorrt2 = 'PB'
        if (icorrt == 'CA') icorrt2 = 'ca'
        if (icorrt == 'XA') icorrt2 = 'xz'
        if (icorrt == 'WI') icorrt2 = 'wi'
        if (icorrt == 'HL') icorrt2 = 'hl'
        if (icorrt == 'PB') icorrt2 = 'pb'
        if ((icorrt /= icorr).and.(icorrt2 /= icorr)) write(7,'(a)') &
             ' *** warning in pseudo : correlation potential does not match'

        write(7,'(a,g11.4)') ' ion charge = ',p_pot%zion(ity)

        if (irel == 'rel') jrel = 1
        if (irel == 'isp') jrel = 2
        if (p_pot%so(ity) .and. (irel /= 'rel' .and. irel /= 'isp')) then
           write(7,*) 'ERROR: file ',fname(ity),' does not contain ', &
                'spin-orbit potentials for atom type ',name(ity)
           write(7,*) ' STOP in pseudo'
           ierr = 1
        endif

        !  read in s,p,d,f pseudopotentials
        !  all reading starts from j=2 to leave room for r=0 at j=1
        read(20,*)
        read(20,*) (p_pot%rs(j,ity),j=2,p_pot%ns(ity))
        do i = 1,p_pot%nlocp(ity)
           read(20,*)
           read(20,*) i0
           if (i0 > 4) then
              write(7,*)'ERROR: More than 4 pseudopot. components'
              write(7,*) 'Do you really need to go beyond f orbitals?'
              write(7,*) 'STOP in pseudo'
              ierr = 207
              return
           endif
           read(20,*) (vspd(j,ity,i0+1),j=2,p_pot%ns(ity))
           !
           !  Pseudopotentials in this format are multiplied by r
           !
           do j=2,p_pot%ns(ity)
              vspd(j,ity,i0+1) = vspd(j,ity,i0+1)/p_pot%rs(j,ity)
           enddo

        enddo
        !  Keep spin-orbit potentials if needed.
        if (p_pot%so(ity)) then
           isy = isy + 1
           do i = 1,nlso
              read(20,*)
              read(20,*) i0
              read(20,*) (vlj(j,isy,i0),j=2,p_pot%ns(ity))
              !
              !  Pseudopotentials in this format are multiplied by r
              !
              do j=2,p_pot%ns(ity)
                 vlj(j,isy,i0) = vlj(j,isy,i0)/p_pot%rs(j,ity)
              enddo
           enddo ! do i = 1,nlso

           if (irel == 'rel') then
              do i0 = 1,2 ! no f orbitals
                 do j = 2, p_pot%ns(ity)
                    temp0 = vspd(j,ity,i0+1) + half*real(i0,dp)*vlj(j,isy,i0)
                    temp1 = vspd(j,ity,i0+1) - half*real(i0+1,dp)*vlj(j,isy,i0)
                    vlj(j,isy,i0) = temp0
                    vlj(j,isy,i0+2) = temp1
                 enddo
              enddo
           elseif (irel == 'isp') then
              if (p_pot%is_so) then
                 write(7,*)'ERROR! YOU USE SPIN POLARISED PSEUDO ', &
                      'POTENTIAL FOR SPIN-ORBIT CALCULATION'
                 ierr = 208
              else
                 do i0 = 1,2
                    do j = 2, p_pot%ns(ity)
                       temp0 = vspd(j,ity,i0+1) + half*real(i0,dp)*vlj(j,isy,i0)
                       temp1 = vspd(j,ity,i0+1) - half*real(i0+1,dp)*vlj(j,isy,i0)
                       vlj(j,isy,i0) = temp0
                       vlj(j,isy,i0+2) = temp1
                    enddo
                 enddo
              endif ! if (p_pot%is_so) then
           endif !  if (irel == 'rel')
        else
           do i = 1,nlso
              read(20,*)
              read(20,*) temp1
              read(20,*) (temp0,j=2,p_pot%ns(ity))
           enddo
        endif ! if (p_pot%so(ity))
        !
        !  Read in the Partial-Core Charge Density, if present.
        !
        read(20,*)
        p_pot%denc(:,ity) = zero
        p_pot%ddenc(:,ity) = zero
        if (trim(nicore) == 'nc') then
           read(20,*) (temp0,j=2,p_pot%ns(ity))
        else
           read(20,*) (p_pot%denc(j,ity),j=2,p_pot%ns(ity))
           do j=2,p_pot%ns(ity)
              p_pot%denc(j,ity) = p_pot%denc(j,ity)/(four*pi  &
                   *p_pot%rs(j,ity)*p_pot%rs(j,ity))
           enddo
        endif
        !
        !  Read valence charge density.
        !
        read(20,*)
        read(20,*) (p_pot%rho_r(j,ity),j=2,p_pot%ns(ity))
!        p_pot%rho_r(1,ity)=p_pot%rho_r(2,ity)

        do j=2,p_pot%ns(ity)
           p_pot%rho_r(j,ity) = p_pot%rho_r(j,ity)/(four*pi*p_pot%rs(j,ity)*p_pot%rs(j,ity))
        enddo

        !
        !  Read in s,p,d pseudo wave functions.
        !  All reading starts from j=2 to leave room for r=0 at j=1
        !  for cutoff radius, choose the maximum cutoff in file.
        !
        tcore = zero
        do i = 1,p_pot%nlocp(ity)
           read(20,*)
           read(20,*) i0,p_pot%eleatm(ity,i),rc
           read(20,*) (p_pot%wfspd(j,ity,i0+1),j=2,p_pot%ns(ity))
           if (tcore < rc) tcore = rc
        enddo
        !
        !  Move core radius to the next grid point.
        !
        do i = 2,p_pot%ns(ity)
           if (p_pot%rs(i,ity) > tcore) exit
        enddo
        tcore = p_pot%rs(i,ity)
        p_pot%rcore(ity) = tcore

     case (MARTINS_WANG)

        write(7,'(a,/)') 'Martins-Wang generator'
        !
        !  Pseudopotentials in this format are given in hartrees, 
        !  hence the factor of two ( hartrees -> rydbergs).
        !
        read(20,*) i, p_pot%icore(ity), i, p_pot%zion(ity),  &
             p_pot%loc(ity), (eleatm_psp(i), i = 1, 3), jrel

        if (p_pot%so(ity) .and. (jrel /= 1)) then
           write(7,*) 'ERROR: file ',fname(ity),' does not contain ', &
                'spin-orbit potentials for atom type ',name(ity)
           write(7,*) ' STOP in pseudo'
           ierr = 209
        endif

        allocate( nloc(p_pot%nlocp(ity)) )
        read(20,*) itmp(1:3)
        kk = 0
        do i = 1, 3
           if (itmp(i) /= 0) then
              kk = kk + 1
              nloc(kk) = i
           endif
        enddo

        do i = 1, p_pot%nlocp(ity)
           p_pot%eleatm(ity,nloc(i)) = eleatm_psp(i)
        enddo
        do i = 2, p_pot%ns(ity)
           read(20,*) p_pot%rs(i,ity),(rtmp(j),j=1,6)
           do j = 1, p_pot%nlocp(ity)
              vspd(i,ity,j) = rtmp(nloc(j)) * two
              p_pot%wfspd(i,ity,j) = rtmp(nloc(j) + 3) * p_pot%rs(i,ity)
           enddo
        enddo
        !  This file has non-linear core correction.
        if (p_pot%icore(ity) /= 0) then
           rewind(20)
           read(20,*)
           read(20,*)
           do i = 2, p_pot%ns(ity)
              read(20,*) p_pot%rs(i,ity),(rtmp(j),j=1,6),p_pot%denc(i,ity)
           enddo
        endif
        !  Keep spin-orbit potentials if needed.
        if (p_pot%so(ity)) then
           rewind(20)
           read(20,*)
           read(20,*)
           do i = 2, p_pot%ns(ity)
              read(20,*) p_pot%rs(i,ity), &
                   (rtmp(j),j=1,6+p_pot%icore(ity)),vlj(1:2,isy,i)
           enddo
           do i0 = 1,2
              do j = 2, p_pot%ns(ity)
                 vlj(j,isy,i0) = vlj(j,isy,i0) * two
                 temp0 = vspd(j,ity,i0+1) + half*real(i0,dp)*vlj(j,isy,i0)
                 temp1 = vspd(j,ity,i0+1) - half*real(i0+1,dp)*vlj(j,isy,i0)
                 vlj(j,isy,i0) = temp0
                 vlj(j,isy,i0+2) = temp1
              enddo
           enddo
        endif
        !
        !  Calculate the radius step parameters - the pseudopotential is
        !  given on a log grid !
        !          r(i) = par_a*exp(par_b*(i-1)) - par_c
        !  The calculation is based on the results at two grid points that
        !  are quite far away, for accuracy.
        !
        j0 = 2*(i_0-1) + 1
        if(j0 > p_pot%ns(ity)) then
           write(7,*) 'ERROR: choose smaller i_0 to get par_a & par_b, i_0 = ',i_0
           write(7,*) 'STOP in pseudo'
           ierr = 210
           return
        endif
        p_pot%par_a(ity) = p_pot%rs(i_0,ity)**2/  &
             (p_pot%rs(j0,ity)-two*p_pot%rs(i_0,ity))
        p_pot%par_b(ity) = log(p_pot%rs(j0,ity)/  &
             p_pot%par_a(ity) + one)/real(j0-1,dp)
        p_pot%par_c(ity) = p_pot%par_a(ity)

        deallocate(nloc)

     case (FHIPP)

        write(7,'(a,/)') 'FHIPP generator'
        read(20,*) p_pot%zion(ity)
        do kk = 1, 10
           read(20,*)
        enddo
        do i = 1, p_pot%nlocp(ity)
           read(20,*) kk, rtmp(1)
           if (kk /= p_pot%ns(ity)-1) then
              write(7,*) 'ERROR in file ',trim(fname(ity))
              write(7,*) 'Radial grid size for angular ',  &
                   'component l = ',i-1,' is not equal to the'
              write(7,*) 'size for other angular components! ',  &
                   kk, p_pot%ns(ity)-1
              write(7,*) 'STOP in pseudo'
              ierr = 211
              return
           endif
           do j = 2, p_pot%ns(ity)
              read(20,*) kk, p_pot%rs(j,ity),  &
                   p_pot%wfspd(j,ity,i), vspd(j,ity,i)
              vspd(j,ity,i) = vspd(j,ity,i) * two
           enddo
           p_pot%par_a(ity) = p_pot%rs(2,ity)/rtmp(1)
           p_pot%par_b(ity) = log(rtmp(1))
           p_pot%par_c(ity) = zero
           if (i > 1 .and. rtmp(1) /= exp(p_pot%par_b(ity))) then
              write(7,*) 'ERROR in file ',trim(fname(ity))
              write(7,*) 'Radial increment for angular ',  &
                   'component l = ',i-1,' is not equal to the'
              write(7,*) 'increment for other angular ',  &
                   'components! ',rtmp(1), exp(p_pot%par_b(ity))
              write(7,*) 'STOP in pseudo'
              ierr = 212
              return
           endif
        enddo

        if (p_pot%icore(ity) > 0) then
           do j=2,p_pot%ns(ity)
              read(20,*) p_pot%rs(j,ity), p_pot%denc(j,ity)
              p_pot%denc(j,ity) = p_pot%denc(j,ity)/(four*pi)
           enddo
        endif

     end select
     !
     !  Print out data and close unit 20
     !
     write(7,*)
     write(7,*) 'logarithmic parameters of pseudopot. radial grid'
     write(7,*) 'r(i) = a*Exp(b(i-1)) - c'
     write(7,93) name(ity), p_pot%par_a(ity), p_pot%par_b(ity),  &
          p_pot%par_c(ity)
93   format(1x,a2,':   a [bohr] = ',g13.6,'  b = ',g13.6,  &
          '  c [bohr] = ',g13.6)
     if (p_pot%icore(ity) > 0) then
        write(7,'(/,a,/)') ' NOTICE: Using Core-Correction !'
        if (psformt == MARTINS_NEW) &
             write(7,'(a,f10.6,a,f10.6,/)') ' cfac = ',cfac,' rcfac = ',rcfac
     else
        write(7,'(/,a,/)') ' NOTICE: No Core-Correction !'
     endif
     if (jrel == 1 .and. (.not. p_pot%so(ity))) &
          write(7,*) ' RELATIVISTIC POTENTIALS : ' &
          ,'IGNORING SPIN-ORBIT PSEUDOPOTENTIAL COMPONENT'
     if (jrel == 2 .and. (.not. p_pot%so(ity))) &
          write(7,*) ' SPIN POLARIZED POTENTIALS : ' &
          ,'IGNORING SECOND PSEUDOPOTENTIAL COMPONENT'

     write(7,'(a,/,4f8.3,/)') ' Occupancy of orbitals : '  &
          ,(p_pot%eleatm(ity,kk),kk=1,p_pot%nlocp(ity))

     write(7,14) p_pot%nlocp(ity), p_pot%loc(ity)-1,  &
          p_pot%rcore(ity)
14   format(' # of pseudopot= ',i1,', Local component is l = ',  &
          i1,', Core radius [bohr] = ',f5.2,/)

     close(unit=20)
     !
     !  Calculate number of valence electrons in two ways - once based
     !  on the number of valence electrons found in ptable, and once
     !  based on the number of electrons per pseudopotential orbital
     !  provided by the user.
     !
     elval1 = elval1 + real(clust%natmi(ity),dp)*p_pot%zion(ity)
     do j = 1, p_pot%nlocp(ity)
        elval2 = elval2 + real(clust%natmi(ity),dp)*  &
             real(p_pot%eleatm(ity,j),dp)
     end do
  enddo                     !  ity = 1, naty (second loop)

  !  the two values for number of valence electrons must agree!
  if (elval1 /= elval2) then
     write(7,*) 'WARNING: # of electrons from pseudopot. '  &
         ,'does not match with atomic configuration.'
         write(7,*) 'This may be a charged pseudopotential.'
  endif
  !
  !  Make sure the user specified enough states to contain all electrons
  !  in the system + two empty states to facilitate convergence.
  !
  if (real(nstate,dp) < elval1*half + two) then
     write(7,*) 'ERROR: Increase number of states to at least',  &
          idint(elval1*half + two) + 1
     write(7,*) 'Stop in pseudo'
     write(7,*) 'In the future this will only be a warning and the number of states will change automatically'
     ierr = 213
     return
  endif
  !
  !  Get coefficients for numerical first derivative based on the 
  !  order specified. They will be used for derivatives of potentials
  write(7,*) 'Using order',2*mor,'for calculating dV/dr in pseudo' 
  call fornberg(1,mor,coe,ierr)
  if (ierr /= 0) then
     write(7,*) ' ERROR: Wrong parameter in call of Fornberg'
     write(7,*) ' STOP in pseudo'
     ierr = 214
     return
  endif

  !  for all atom types...
  isy = 0
  do ity = 1, naty
     npoint = p_pot%ns(ity)
     if (p_pot%so(ity)) then
        isy = isy + 1
        if ((psformt == MARTINS_NEW).or.(psformt == MARTINS_WANG)) then

           do i0 = 1,4
              do j = -p_pot%norder,npoint
                 vlj(j,isy,i0) = vlj(j,isy,i0) - vspd(j,ity,p_pot%loc(ity))
              enddo
           enddo
          p_pot%rs(1,ity) = zero
           do i0 = 1,2
              j0 = i0 + 1
              delr = p_pot%rs(2,ity)-p_pot%rs(1,ity)
              fspd = p_pot%wfspd(1,ity,j0)*p_pot%wfspd(1,ity,j0)
              clj(i0,isy) = vlj(1,isy,i0)*fspd*delr
              clj(i0+2,isy) = vlj(1,isy,i0+2)*fspd*delr
              do j = 2,npoint-1
                 delr = half*(p_pot%rs(j+1,ity)-p_pot%rs(j-1,ity))
                 fspd = p_pot%wfspd(j,ity,j0)*p_pot%wfspd(j,ity,j0)
                 clj(i0,isy) = clj(i0,isy) + vlj(j,isy,i0)*fspd*delr
                 clj(i0+2,isy) = clj(i0+2,isy) + vlj(j,isy,i0+2)*fspd*delr
              enddo
              delr = p_pot%rs(npoint,ity) - p_pot%rs(npoint-1,ity)
              fspd = p_pot%wfspd(npoint,ity,j0)*p_pot%wfspd(npoint,ity,j0)
              clj(i0,isy) = clj(i0,isy) + vlj(npoint,isy,i0)*fspd*delr
              clj(i0+2,isy) = clj(i0+2,isy) + vlj(npoint,isy,i0+2)*fspd*delr
              do j=2,npoint
                 vlj(j,isy,i0)=vlj(j,isy,i0)*p_pot%wfspd(j,ity,j0)/ &
                      sqrt(abs(clj(i0,isy)))/p_pot%rs(j,ity)
                 vlj(j,isy,i0+2)=vlj(j,isy,i0+2)*p_pot%wfspd(j,ity,j0)/ &
                      sqrt(abs(clj(i0+2,isy)))/p_pot%rs(j,ity)
              enddo

              do j=-p_pot%norder,mor
                 vlj(j,isy,i0)=vlj(2*mor,isy,i0)
                 vlj(j,isy,i0+2)=vlj(2*mor,isy,i0+2)
              enddo
           enddo
        endif !if ((psformt == MARTINS_NEW).or.(psformt == MARTINS_WANG))

        do i0 = 1,2
           j0 = i0 + 1
           temp0 = p_pot%so_hcub/real(j0+i0,dp)
           
           do j = -p_pot%norder,npoint
              p_pot%vsor(j,isy,i0) = temp0*two*(vlj(j,isy,i0)-vlj(j,isy,i0+2))
              p_pot%vionr(j,isy,i0) = temp0*(real(j0,dp) * vlj(j,isy,i0) &
                   + real(i0,dp)*vlj(j,isy,i0+2))
           enddo
        enddo

!        p_pot%vsor(1,isy,1:2) = zero
!        p_pot%vionr(1,isy,1:2) = zero
        
        write(7,*)''
        write(7,*)'KB integrals for channels P_up P_dn D_up D_dn'
        write(7,87) name(ity),clj(1,isy),clj(3,isy),clj(2,isy),clj(4,isy)
87 format(2x,a2,2x,' SO K-B Int. [Ry]', 4(3x,f12.8))

        do i0 = 1,2
           j0 = i0 + 2

           if (abs(clj(i0,isy))<small .or. abs(clj(j0,isy))<small)  then
              write(7,*) 'WARNING: The K-B Integral is too small!!'
              !p_pot%cc(i0,isy) = zero
           endif

           if (clj(i0,isy)/abs(clj(i0,isy)) == clj(j0,isy)/abs(clj(j0,isy))) then
              p_pot%cc(i0,isy) = clj(i0,isy)/abs(clj(i0,isy))
           else
              p_pot%cc(i0,isy) = zero
              if(abs(clj(i0,isy))>abs(clj(j0,isy)))then
                 p_pot%cc(i0,isy) = clj(i0,isy)/abs(clj(i0,isy))
              else
                 p_pot%cc(i0,isy) = clj(j0,isy)/abs(clj(j0,isy))
              endif
              write(7,*)''
              write(7,*)'WARNING ! in pseudo: atom ',name(ity),' L=',i0
              write(7,*)'sign of K-B integral is not equal for the two components'
!              write(7,*)'This projector will not contribute to spin-orbit energy'
              write(7,*)'*******************************************************'
           endif
        enddo
     endif ! if (p_pot%so(ity))
     !
     !  temporary variables 
     npoint = p_pot%ns(ity)
     a = p_pot%par_a(ity)
     b = p_pot%par_b(ity)
     c = p_pot%par_c(ity)
     !
     !  Define local component (notice that pseudopotentials of Martins-Wang
     !  format with p_pot%loc(:) = 0 already have local component defined).
     !
     do  j =  2, npoint
        p_pot%vion(j,ity) = vspd(j,ity,p_pot%loc(ity)) 
     enddo
     !
     !  Filling in the value at r<=0 for the pseudopotentials and wave
     !  functions. Value at r=0 NOT provided by pseudopotential file!
     !  For negative r, set rs arbitrarily be j-1 for preventing
     !  division by zero. This only affects several points near the
     !  atom center, which is a negligible effect.
     !
     do j = -p_pot%norder,1
        p_pot%rs(j,ity) = real(j-1,dp)
        vspd(j,ity,:) = vspd(2,ity,:)
        p_pot%vion(j,ity) = p_pot%vion(2,ity)
        p_pot%wfspd(j,ity,:) = zero
        p_pot%denc(j,ity) = p_pot%denc(2,ity)
        if (p_pot%so(ity)) vlj(j,isy,:) = zero
     enddo
     !
     !  Change non-local components of the pseudopotential into the
     !  difference between them and the local component, i.e., give
     !  pseudopotentials to the main program already in Kleinman-
     !  Bylander form.
     !
     do j0 = 1, p_pot%nlocp(ity)
        do i0 = -p_pot%norder,p_pot%ns(ity)
           vspd(i0,ity,j0) = vspd(i0,ity,j0) - p_pot%vion(i0,ity)
        enddo
     enddo
     !
     !  Calculate the derivative of the ionic (=local) pseudopotential,
     !  wfspd*vspd/r, its derivative, and the derivative of the 
     !  core charge density.
     !
     !  Comment: the wavefunctions and pseudopotentials are specified
     !  on a logarithmic grid, c but the numerical derivative assumes a
     !  uniform grid, therefore: d/dr = di/dr*d/di. 
     !  The term 1/(b*(rs(j,ity)+c)) is di/dr, hence its appearance
     !  wherever a derivative is taken.
     do j = -p_pot%norder,mor-1
        p_pot%dvion(j,ity) = zero
     enddo

     do j = mor, npoint - mor
        temp0 = zero
        do kk = -mor,mor,1
           temp0 = temp0 + coe(kk)*p_pot%vion(j+kk,ity)
        enddo
        p_pot%dvion(j,ity) = temp0/(b*(p_pot%rs(j,ity)+c))
     enddo
     !  For the last few points, vion is just the Coulomb potential
     !  (for rs is beyond the core, so the derivative can be taken
     !  analytically. The factor of 2 is for conversion from hartrees to
     !  rydbergs.
     do j = npoint - mor + 1, npoint
        p_pot%dvion(j,ity) = two*p_pot%zion(ity)/  &
             (p_pot%rs(j,ity)*p_pot%rs(j,ity))
     enddo

     if (p_pot%uspline(ity)) then
        call spline(p_pot%rs(-p_pot%norder:npoint,ity), &
                    p_pot%vion(-p_pot%norder:npoint,ity), npoint+1+p_pot%norder, &
                    zero,zero,p_pot%d2vion(-p_pot%norder:npoint,ity))
     endif

     !  Set first few points where the derivative is mostly numerical
     !  noise anyway to zero.
     do j=-p_pot%norder,npoint
        if(p_pot%rs(j,ity) <= 3.d-4) p_pot%dvion(j,ity) = zero
     enddo

     !  For all non-local pseudopotential components...
     do i0 = 1, p_pot%nlocp(ity) 
        if (i0 /= p_pot%loc(ity)) then
        
           do j = mor, npoint
              p_pot%vw(j,ity,i0) = vspd(j,ity,i0)*  &
                   p_pot%wfspd(j,ity,i0)/p_pot%rs(j,ity)
           enddo
        
           do j = mor, npoint-mor
              temp0 = zero
              do kk = -mor,mor,1
                 temp0 = temp0 + coe(kk)*p_pot%vw(j+kk,ity,i0)
              enddo
              p_pot%dvw(j,ity,i0) = temp0/(b*(p_pot%rs(j,ity) + c))
           enddo
        endif

        if (p_pot%so(ity) .and. i0 < 3) then
           do j = mor, npoint-mor
              temp0 = zero
              temp1 = zero
              do kk = -mor,mor,1
                 temp0 = temp0 + coe(kk)*vlj(j+kk,isy,i0)
                 temp1 = temp1 + coe(kk)*vlj(j+kk,isy,i0+2)
              enddo

              dvlj(j,isy,i0) = temp0/(b*(p_pot%rs(j,isy) + c))

              dvlj(j,isy,i0+2) = temp1/(b*(p_pot%rs(j,isy) + c)) 
           enddo

           j0 = i0 + 1
           do j = mor, npoint-mor
              p_pot%dvr_so(j,isy,i0) = &
                   two * (dvlj(j,isy,i0) - dvlj(j,isy,i0+2))
              p_pot%dvr_ion(j,isy,i0) = &
                   real(j0,dp) * vlj(j,isy,i0) + real(i0,dp) &
                   * vlj(j,isy,i0+2)
           enddo

           do j = npoint - mor + 1, npoint
              p_pot%dvr_so(j,isy,i0) = zero
              p_pot%dvr_ion(j,isy,i0) = zero
           enddo

           do j = -p_pot%norder,2*mor
              p_pot%dvr_so(j,isy,i0) = p_pot%dvr_so(mor*2,isy,i0)
              p_pot%dvr_ion(j,isy,i0) = p_pot%dvr_ion(mor*2,isy,i0)
           enddo
           if (p_pot%uspline(ity)) then
              call spline(p_pot%rs(-p_pot%norder:npoint,ity), &
                       p_pot%vionr(-p_pot%norder:npoint,ity,i0), npoint+1+p_pot%norder, &
                       zero,zero,p_pot%d2vr_ion(-p_pot%norder:npoint,ity,i0))
              call spline(p_pot%rs(-p_pot%norder:npoint,ity), &
                       p_pot%vsor(-p_pot%norder:npoint,ity,i0), npoint+1+p_pot%norder, &
                       zero,zero,p_pot%d2vr_so(-p_pot%norder:npoint,ity,i0))
           endif

        endif    !if (p_pot%so(ity) .and. i0 < 3)

        if (i0 /= p_pot%loc(ity)) then

        !  The last few one may be set to zero, because rs is large
        !  enough to be outside the non-local core.
           do j = npoint - mor + 1, npoint
              p_pot%dvw(j,ity,i0) = zero
           enddo
           
           do j = -p_pot%norder,mor
              p_pot%vw(j,ity,i0) = p_pot%vw(mor*2,ity,i0)
           enddo
           do j = -p_pot%norder,2*mor
              p_pot%dvw(j,ity,i0) = p_pot%dvw(mor*2,ity,i0)
           enddo
           if (p_pot%uspline(ity)) then
              call spline(p_pot%rs(-p_pot%norder:npoint,ity), &
                       p_pot%vw(-p_pot%norder:npoint,ity,i0), npoint+1+p_pot%norder, &
                       zero,zero,p_pot%d2vw(-p_pot%norder:npoint,ity,i0))
           endif
        endif               !if(i0 /= p_pot%loc(ity))

     enddo                  ! i0 = 1, p_pot%nlocp(ity)
     !
     !  Calculate the radial derivative of the core charge density, if
     !  exists.
     !
     if (p_pot%icore(ity) > 0) then
        do j = mor, npoint - mor
           temp0 = zero
           do kk = -mor,mor,1
              temp0 = temp0 + coe(kk)*p_pot%denc(j+kk,ity)
           enddo
           p_pot%ddenc(j,ity) = temp0/(b*(p_pot%rs(j,ity) + c))
        enddo
        do j = -p_pot%norder,mor-1
           p_pot%ddenc(j,ity) = zero
        enddo
        do j = npoint - mor + 1, npoint
           p_pot%ddenc(j,ity) = zero
        enddo
        if (p_pot%uspline(ity)) then
           call spline(p_pot%rs(-p_pot%norder:npoint,ity), &
                    p_pot%denc(-p_pot%norder:npoint,ity), npoint+1+p_pot%norder, &
                    zero,zero,p_pot%d2denc(-p_pot%norder:npoint,ity))
        endif
     endif
     !
     !   Calculate the radial derivative of the pseudo wavefunctions.
     !   This is needed for LDA+U.
     !
     do i0 = 1, p_pot%nlocp(ity)
        if ( p_pot%uu(i0,ity) == zero .and. p_pot%jj(i0,ity) == zero ) cycle

        do j = mor, npoint - mor
           temp0 = zero
           do kk = -mor,mor,1
              temp0 = temp0 + coe(kk)*p_pot%wfspd(j+kk,ity,i0)
           enddo
           p_pot%dwfspd(j,ity,i0) = temp0/(b*(p_pot%rs(j,ity) + c))
        enddo
        do j = -p_pot%norder,mor-1
           p_pot%dwfspd(j,ity,i0) = zero
        enddo
        do j = npoint - mor + 1, npoint
           p_pot%dwfspd(j,ity,i0) = zero
        enddo

     enddo
  
  enddo                     ! ity = 1, naty (third loop)
  !
  !  Calculate the VOLUME atomic charge density and its derivative,
  !  for each atom type, on the logarithmic pseudopotential grid.
  !
  !  Go over all atom types.
  do ity = 1, naty
     npoint = p_pot%ns(ity)

     !  Get the total RADIAL charge density for each atom type, by
     !  summing over the number of atoms per orbital, eleatm(ity,iorb),
     !  times the square of that orbital's pseudo wavefunction,
     !  wfspd(i,ity,iorb). Then divide by 4*pi*r^2 to convert the
     !  RADIAL density into a VOLUME one.
     if(.not. p_pot%rvcd(ity)) then
        do i = 2, npoint
           temp0 = zero
           do iorb = 1, p_pot%nlocp(ity)
              temp0 = temp0 +  &
                   p_pot%eleatm(ity,iorb)*p_pot%wfspd(i,ity,iorb)**2
           enddo
           p_pot%rho_r(i,ity) = temp0/(four*pi*p_pot%rs(i,ity)*  &
                p_pot%rs(i,ity))
        enddo
     else
        if(psformt /= MARTINS_NEW) then
           write(7,*)'Error!!'
           write(7,*)' valence charge density can be read only in MARTINS_NEW format'
           write(7,*)'Stop in pseudo.f90'
           ierr = 23
        endif
     endif
           
     !  Set the "negative index" points (used for derivation near the
     !  origin) to zero.

     do i = -p_pot%norder,1
        p_pot%rho_r(i,ity) = p_pot%rho_r(2,ity)
     enddo

     !  Calculate the radial derivative of the atomic charge density.
     !  It is used for calculating forces.

     do i = -p_pot%norder, mor
        p_pot%drhodr(i,ity) = zero
     enddo
     do i = mor, npoint - mor
        temp0 = zero
        do kk = -mor,mor,1
           temp0 = temp0 + coe(kk)*p_pot%rho_r(i+kk,ity)
        enddo
        p_pot%drhodr(i,ity) = temp0/(b*(p_pot%rs(i,ity)+c))
     enddo
     do i = npoint - mor + 1, npoint, 1
        p_pot%drhodr(i,ity) = zero
     enddo
     if (p_pot%uspline(ity)) then
        call spline(p_pot%rs(-p_pot%norder:npoint,ity), &
                 p_pot%rho_r(-p_pot%norder:npoint,ity), npoint+1+p_pot%norder, &
                 zero,zero,p_pot%d2rhodr(-p_pot%norder:npoint,ity))
     endif
  enddo
  !
  !  Calculate (in ekbi) the Kleinman-Bylander normalization factor 
  !  <psi|(V_l-V_locaa)l|psi> (this is the denominator of Eq. (3) in
  !  Chelikowsky et al., PRL, 72, 1240 (1994)). As this is needed
  !  for the non-local pseudopotential only, the results for the
  !  local pseudopotential are junk and are never used.
  !
  !  Also calculate (in psinorm) the integral of psi* times psi over
  !  the radial grid (in psinorm) in order to make sure it is close
  !  to one, and write it out to parsec.out. This is for information
  !  only - psinorm is not used anywhere else. 
  !
  write(7,'(/,a,/)') ' Normalization check and Kleinman-Bylander Integral'
  write(7,'(33x,a,17x,a,17x,a)') 's','p','d'

  do ity = 1,naty
     npoint = p_pot%ns(ity)
     do j0 = 1, p_pot%nlocp(ity)
        delr = p_pot%rs(2,ity) - p_pot%rs(1,ity)
        fspd = p_pot%wfspd(1,ity,j0)*p_pot%wfspd(1,ity,j0)
        psinorm(j0) = fspd*delr
        p_pot%ekbi(j0,ity) = vspd(1,ity,j0)*fspd*delr
        do j = 2,npoint-1
           delr = half*(p_pot%rs(j+1,ity)-p_pot%rs(j-1,ity))
           fspd = p_pot%wfspd(j,ity,j0)*p_pot%wfspd(j,ity,j0)
           psinorm(j0) =  psinorm(j0) + fspd*delr
           p_pot%ekbi(j0,ity) = p_pot%ekbi(j0,ity) +   &
                vspd(j,ity,j0)*fspd*delr
        enddo
        delr = p_pot%rs(npoint,ity) - p_pot%rs(npoint-1,ity)
        fspd = p_pot%wfspd(npoint,ity,j0)*p_pot%wfspd(npoint,ity,j0)
        psinorm(j0) =  psinorm(j0) + fspd*delr
        p_pot%ekbi(j0,ity) = p_pot%ekbi(j0,ity) +  &
             vspd(npoint,ity,j0)*fspd*delr
     enddo
     write(7,77) name(ity), (psinorm(j), j=1,p_pot%nlocp(ity))
     write(7,78) name(ity), (p_pot%ekbi(j,ity),j=1,p_pot%nlocp(ity))
     write(7,*)
  enddo

77 format(2x,a2,2x,' Normalization', 3(3x,f15.8))
78 format(2x,a2,2x,' K-B Int. [Ry]', 3(3x,f15.8))
  !
  !  Because the Kleinman-Bylander normalization factor is eventually used
  !  in a denominator, if it is close to zero its inverse will explode.
  !  Therefore, check its value and, if it is too small for any non-local
  !  component of any atom, print some warning.
  !
  do ity  = 1, naty
     do j = 1, p_pot%nlocp(ity)
        if (j /= p_pot%loc(ity)) then
           if (abs(p_pot%ekbi(j,ity)) < small)  then
               ! this warning should be in bold :)
              write(7,*) 'WARNING: The K-B Integral is too small!!'
              write(7,79) j-1,ity,p_pot%ekbi(j,ity)
           endif
           !
           !  Renormalize K-B projectors with normalization factor
           !
           temp0 = sqrt(abs(p_pot%ekbi(j,ity)))
           p_pot%vw(:,ity,j) = p_pot%vw(:,ity,j)/temp0
           p_pot%dvw(:,ity,j) = p_pot%dvw(:,ity,j)/temp0
           p_pot%d2vw(:,ity,j) = p_pot%d2vw(:,ity,j)/temp0
           p_pot%ekbi(j,ity) = p_pot%ekbi(j,ity)/temp0/temp0
        endif
     enddo
  enddo
  write(7,*) 'K-B Energy Integrals renormalized to unity'
  write(7,*)

79 format(1x,'the K-B factor of L=',1x,i1,1x,'for atom type #', &
        1x,i2,1x,'is',1x,f15.8)
  !
  !  Include volume factor in wfspd
  !
  do ity = 1, naty
     do i = 2, p_pot%ns(ity)
        do iorb = 1, p_pot%nlocp(ity)
           p_pot%wfspd(i,ity,iorb) = p_pot%wfspd(i,ity,iorb)/p_pot%rs(i,ity)
           p_pot%dwfspd(i,ity,iorb) = ( p_pot%dwfspd(i,ity,iorb) - &
                p_pot%wfspd(i,ity,iorb) ) /p_pot%rs(i,ity)
        enddo
     enddo
  enddo

  !  deallocate the arrays
  if (allocated (vspd)) deallocate (vspd)
  if (allocated (fname)) deallocate(fname)
  if (allocated (fnames)) deallocate(fnames)
  if (allocated (dvlj)) deallocate(dvlj)
  if (allocated (clj)) deallocate(clj)
  if (allocated (vlj)) deallocate(vlj)
  return

  !  print out error message and exit if files do not exist
50 write(7,92) trim(fname(ity))
  write(7,*) ' STOP in pseudo'
  ierr = 216
  return
92 format('ERROR: The file',1x,a,1x,' does not seem to exist',/)

end subroutine pseudo
!===============================================================
