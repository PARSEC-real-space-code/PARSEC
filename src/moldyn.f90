!===============================================================
!
! Copyright (C) 2005 Finite Difference Research Group
! This file is part of parsec, http://www.ices.utexas.edu/parsec/
!
! modified by M. Alemany (2003)
!
! This subroutine computes new positions, velocities, and
! accelerations for the atoms given an equation of motion. 
!
! For canonical ensemble : it couples the system to a Langevin
! heat bath.
! It solves the equation: d^2(r)/dt^2=F/m-beta*v+G/m, where, for
! each atom, r is the position vector, F is the previously
! computed quantum force, beta is a viscosity parameter, v is the
! velocity, m is the nass, and G is a "white "noise" (Gaussian)
! force with an autocorrelation of sqrt(2*beta*mass*kT*delta(t-t')),
! delta() being the Dirac delta function.
! For a detailed discussion of the Langevin equation see:
! Reif, Statistical and Thermal Physics, section 15.5, pp. 560
! -567 (mcGraw-Hill, New York, 1965).
!
! The subroutine integrates the Langevin equation at finite time
! steps, using the method suggested by Biswas and Hamann [Phys.
! Rev. B, 34 895 (1986)] but using the integration method
! suggested by Beeman  [J. Comput. Phys., 20, 130 (1976)] eqns
! 10A and 10C in that paper. The previous estimate of velocity 
! is used for the friction force.
!
! For micro-canonical ensemble: it integrates the newtonian eqns
! of motion.
! It solves the equation: d^2(r)/dt^2=F/m, where, per each atom,
! r is the position vector, F is the previously computed quantum
! force and m is the mass.
!
! The subroutine integrates the above equation at finite time
! steps, using a modification of the Verlet algorithm suggested
! by Beeman [J. Comput. Phys., 20, 130 (1976)] This is same as
! verlet algorithm except the estimate of velocity is improved.
!
! IMPORTANT : To decide whether it is microcanonical or canonical
! equation of motion this routine just tests the value of beta.
! If beta < 10^-6 then it is microcanonical else it is canonical.
!
! NOTE: This subroutine uses Hartree as the units. The reason for
! this choice rather than Rydbergs is that unless one puts in
! factors of 2 and sqrt(2) in in force and velocity, the are
! directly not in the same units and so you cannot just add them
! without a problem. If you consistently use hartrees you 
! are always fine (they are atomic units and hence compatible
! with everything else).
!
!---------------------------------------------------------------
subroutine moldyn(clust,mol_dynamic,pbc,imove,bdev)

  use constants
  use cluster_module
  use molecular_dynamic_module
  use pbc_module
  implicit none
  !
  ! Input/Output variables:
  !
  ! the cluster
  type(cluster), intent(inout) :: clust
  ! molecular dynamic related data
  type (molecular_dynamic), intent(inout) :: mol_dynamic
  ! periodic boundary conditions
  type (pbc_data), intent(in) :: pbc
  ! # of time steps
  integer, intent(inout) :: imove
  ! total energy per atom in eV
  real(dp), intent(in) :: bdev
  !
  ! Work variables:
  !
  ! atomic mass of each ion type
  real(dp) :: amass(clust%atom_num)
  ! array to specify if atom is mobile   
  real(dp) :: amove(clust%atom_num)

  ! atomic coordinates
  real(dp), dimension(clust%atom_num) ::  xatm, yatm, zatm

  ! Langevin trajectory stuff
  ! current positions, old positions (from previous step)
  real(dp), dimension(clust%atom_num) ::  &
       xcur, ycur, zcur, xold, yold, zold

  ! velocities (from previous step and the current ones); they are
  ! stored in the same variable (although it is called old)
  real(dp), dimension(clust%atom_num) :: vxold, vyold, vzold

  ! old accelerations (from previous step and current ones); again,
  ! they are stored in the same variable (although it is called old)
  real(dp), dimension(clust%atom_num) :: accxold, accyold, acczold

  ! total actual number of atoms
  integer natom
  ! number of different atom types
  integer naty
  ! Temperature in hartree
  real(dp) :: tempeh
  ! accumulators for total mass, total mobile mass and center of mass velocity
  real(dp) :: vxav, vyav, vzav, tmass, tmobilemass
  ! thermal velocity along each axis
  real(dp) :: vel
  ! the standard deviation of the random force
  real(dp) :: sigm
  ! temporary storage of older position components
  real(dp) :: xoldr, yoldr, zoldr
  ! temporary storage of older velocity components
  real(dp) :: vxoldr,vyoldr,vzoldr
  ! temporary storage of older acceleration components
  real(dp) :: accxoldr, accyoldr, acczoldr
  ! minimum and maximum of coordinates along each axis
  real(dp) :: xmin, xmax, ymin, ymax, zmin, zmax
  ! amount by which to shift xcur,ycur,zcur to center coordinates
  ! about the origin
  real(dp) :: xmov, ymov, zmov
  ! total kinetic energy [Ry] and total kinetic energy per atom [eV].
  real(dp) :: tke, tkeev, tke2
  ! effective actual temperature 
  real(dp) :: tcalc
  ! factor with which to rescale the velocities
  ! (ratio of t_input /tcalc)
  real(dp) :: rescale, tcorrect

  ! current and final temperatures (in kelvin)
  real(dp) :: tempi, tempf
  ! time step, friction coefficient
  real(dp) :: deltat,beta

  ! number of mobile atoms (from clust%mvat array)
  integer :: nmobileatom

  ! Counter of atom movement cycles in molecular dynamics. It is
  ! different from imove only when using stairstep cooling, where
  ! imove is reset to zero for each "stair" and iframe is the total
  ! number of movements.
  integer iframe

  ! counters
  integer i,j,l,ii
  !
  ! External functions:
  !
  ! a random Gaussian variable with a mean of zero and a standard 
  ! deviation of one
  real(dp), external :: gdev

  real(dp), dimension(clust%atom_num) :: xold2, yold2, zold2
  real(dp), dimension(clust%atom_num) :: vxold2, vyold2, vzold2
  real(dp), dimension(clust%atom_num) :: accxold2, accyold2, acczold2
  integer, parameter :: cdflag = 0
  !Nose-Hoover Thermostat stuff
  real(dp), allocatable :: Jacob(:,:), deltak(:), h(:), DWORK(:)
  real(dp), allocatable :: vxcur(:), vycur(:), vzcur(:)
  integer,  allocatable :: IPIV(:)
  real(dp) :: qnose(mol_dynamic%nchain)
  real(dp) :: nose, nosez, nosezd, nosezdd, nosezdcur
  real(dp) :: nosezold, nosezdold, enose, thkh, nodof, nosezoldr
  real(dp) :: nosezcur,nosezddold
  real(dp) :: exn, evn
  logical :: finiter, beeman
  integer :: info, inr,k
  integer :: nchain
  !Artificial force stuff
  integer :: ni, nj, ia, ja
  real(dp) :: dxx, dyy, dzz, dcell, rij, rij2
  real(dp) :: frep
  real(dp), dimension(:,:), allocatable :: forc
  !---------------------------------------------------------------

  natom = clust%atom_num
  naty = clust%type_num
  xatm = clust%xatm
  yatm = clust%yatm
  zatm = clust%zatm
  amass = clust%amass

  tempi = mol_dynamic%tempi
  tempf = mol_dynamic%tempf
  beta = mol_dynamic%friction_coef
  deltat =  mol_dynamic%time_step
  iframe = mol_dynamic%iframe

  xcur = mol_dynamic%xcur
  ycur = mol_dynamic%ycur
  zcur = mol_dynamic%zcur

  xold = mol_dynamic%xold
  yold = mol_dynamic%yold
  zold = mol_dynamic%zold

  vxold = mol_dynamic%vxold
  vyold = mol_dynamic%vyold
  vzold = mol_dynamic%vzold

  accxold = mol_dynamic%accxold
  accyold = mol_dynamic%accyold
  acczold = mol_dynamic%acczold

  nose = mol_dynamic%nose 
  nosez = mol_dynamic%nosez
  nosezd = mol_dynamic%nosezd
  nosezdd = mol_dynamic%nosezdd
  nosezcur = mol_dynamic%nosezcur
  nchain = mol_dynamic%nchain
  beeman = .false.

  if(nose .ge. 0.d0) then
    allocate(Jacob(3*natom+1,3*natom+1))
    allocate(IPIV(3*natom+1))
    allocate(deltak(3*natom+1))
    allocate(DWORK(3*natom+1))
    allocate(h(3*natom+1))
    allocate(vxcur(natom))
    allocate(vycur(natom))
    allocate(vzcur(natom))
  endif


  if (imove == 0 ) then
    if(mol_dynamic%class_restart) then
      write(7,*) "r, v, and a are read from the last step of CMD"
      xatm = mol_dynamic%xcur2
      yatm = mol_dynamic%ycur2
      zatm = mol_dynamic%zcur2
      xold = mol_dynamic%xold2
      yold = mol_dynamic%yold2
      zold = mol_dynamic%zold2
      xcur = mol_dynamic%xcur2
      ycur = mol_dynamic%ycur2
      zcur = mol_dynamic%zcur2
      vxold = mol_dynamic%vxold2
      vyold = mol_dynamic%vyold2
      vzold = mol_dynamic%vzold2
      accxold = mol_dynamic%accxold2
      accyold = mol_dynamic%accyold2
      acczold = mol_dynamic%acczold2
    endif
  endif

  write(7,*) xold
  write(7,*) xcur
  write(7,*) vxold
  write(7,*) accxold

  !Classical Repulsive Forces if necessary
  allocate(forc(3,natom))
  forc(:,:) = 0.d0
  ni = 0
  do ia=1,naty
    do k=1,clust%natmi(ia) 
      ni=ni+1
      nj=0 
      if(clust%rpat(ni) == 1) then
        do ja=1,naty
          do l=1,clust%natmi(ja)
            nj =nj+1
            !write(7,*) "anums",ni,nj
            if(clust%rpat(nj) == 1) then
              if (ni/=nj) then
                dxx =xatm(ni)-xatm(nj)
                dyy =yatm(ni)-yatm(nj)
                dzz =zatm(ni)-zatm(nj)
                if(pbc%per==3) then
                  dcell =abs(pbc%latt_vec(1,1)+pbc%latt_vec(1,2)+pbc%latt_vec(1,3))
                  if(abs(dxx) > dcell*half) dxx = dxx -  sign(dcell,dxx) 
                  dcell =abs(pbc%latt_vec(2,1)+pbc%latt_vec(2,2)+pbc%latt_vec(2,3))
                  if(abs(dyy) > dcell*half) dyy = dyy -  sign(dcell,dyy) 
                  dcell =abs(pbc%latt_vec(3,1)+pbc%latt_vec(3,2)+pbc%latt_vec(3,3))
                  if(abs(dzz) > dcell*half) dzz = dzz -  sign(dcell,dzz) 
                endif
                rij2=(dxx*dxx+dyy*dyy+dzz*dzz)
                rij =sqrt(rij2)
                !write(7,*) "Distance", rij
                if(rij < 6.d0) then
                  !write(7,*) "Computing force"
                  frep = clust%cr6 /(rij2*rij2*rij2)
                  forc(1,ni)=forc(1,ni)+dxx*frep/rij
                  forc(2,ni)=forc(2,ni)+dyy*frep/rij
                  forc(3,ni)=forc(3,ni)+dzz*frep/rij
                endif
                if(rij <= clust%rex + 1.d0 ) then
                  frep = exp( - clust%cex * (rij - clust%rex))
                  forc(1,ni)=forc(1,ni)+dxx*frep/rij
                  forc(2,ni)=forc(2,ni)+dyy*frep/rij
                  forc(3,ni)=forc(3,ni)+dzz*frep/rij
                endif
              endif
            endif
          enddo
        enddo
      endif
      if(abs(forc(1,ni)) > 1.d0) forc(1,ni) =  sign(clust%rf, forc(1,ni)) 
      if(abs(forc(2,ni)) > 1.d0) forc(2,ni) =  sign(clust%rf, forc(2,ni)) 
      if(abs(forc(3,ni)) > 1.d0) forc(3,ni) =  sign(clust%rf, forc(3,ni)) 
    enddo
  enddo
  write(*,*) iframe, "Total force and force including art_rep (in Ry/au)"
  do j = 1,natom
     write(7,*) forc(:,j)
     write(7,41)  "F",j,xatm(j),yatm(j),zatm(j), &
                 clust%force(1,j),clust%force(2,j),clust%force(3,j)
     if(clust%rpat(j) == 1) then
       clust%force(:,j) = clust%force(:,j) + forc(:, j) 
       write(7,41)  "A",j,xatm(j),yatm(j),zatm(j), &
                   clust%force(1,j),clust%force(2,j),clust%force(3,j)
     endif
  end do

  amove = real(clust%mvat,dp)
  ! number of mobile atoms
  nmobileatom = sum (clust%mvat)
  !
  ! change temperature from kelvin (input) to atomic unit (Hartree)
  tempeh = tempi * tempkry * half
 
  !
  ! if initial run, there are no previous steps to use a random guess 
  ! for the velocity
  !
  if ((imove == 0) .and. (.not. mol_dynamic%is_restart) .and. (.not. mol_dynamic%class_restart)) then
     ! set accumulators for total mass and center of mass velocity
     mol_dynamic%initen = bdev * real(natom,dp) / rydberg
     if(mol_dynamic%hybrid) write(7,*) "Doing first initialization MD run"
     if(mol_dynamic%hybrid) write(7,*) "Initial total energy is", mol_dynamic%initen
     mol_dynamic%de = zero
     do  j=1, natom
        ! find thermal velocity in each direction, x, y, z
        !                              (0.5kT=0 .5m(vel)^2)
        vel = sqrt(tempeh/amass(j))
        ! assign velocity based on a maxwell distribution:
        !      p(v) ~ exp (-0 .5*v^2/vel^2)
        ! gdev returns a random number x with a propbability
        !      p(x) ~ exp ( -0.5*x^2)
        ! and therefore the needed v is simply vel*x
        vxold(j) = vel*gdev(mol_dynamic)*amove(j)
        write(7,*) vxold(j) / vel
        vyold(j) = vel*gdev(mol_dynamic)*amove(j)
        write(7,*) vyold(j) / vel
        vzold(j) = vel*gdev(mol_dynamic)*amove(j)
        write(7,*) vzold(j) / vel
     enddo
     ! determined to total mass and mass of mobile atoms
     tmass = sum (amass)
     tmobilemass = dot_product(amass,amove)

     ! accumulate velocities for center of mass velocity calculation
     vxav = dot_product (amass,vxold)
     vyav = dot_product (amass,vyold)
     vzav = dot_product (amass,vzold)

     ! find center of mass velocity
     vxav = vxav/tmobilemass
     vyav = vyav/tmobilemass
     vzav = vzav/tmobilemass
     ! Correct random velocities so that the center of mass does not
     ! move and also find the kinetic energy of the atoms for
     ! calculation of temperature. Note that the kinetic energy is in
     ! hartrees.
     tke = zero
     tke2= zero

     vxold = vxold-vxav*amove
     vyold = vyold-vyav*amove
     vzold = vzold-vzav*amove 


     write(777,*) iframe, mol_dynamic%mdtime, 0, 0 
     if(mol_dynamic%read_inivel) then
        write(7,*) "Reading initial velocities"
        do j = 1, natom
          read(999,*) vxold(j), &
                      vyold(j), &
                      vzold(j)
        enddo 
     endif
     do j = 1, natom
        tke = tke+half*amass(j)*(vxold(j)*vxold(j)+vyold(j)* &
             vyold(j)+vzold(j)*vzold(j))
        ! Report velocitiy
        write(777,'(3f21.17,i5)') vxold(j), &
                                  vyold(j), &
                                  vzold(j), 0 
        ! Report Initial Velocity
        write(999,'(3f21.17,i5)') vxold(j), &
                                  vyold(j), &
                                  vzold(j) 
     enddo
     ! Calculate the actual temperature and rescale the velocities so
     ! that you get the input temperature.
     if (nmobileatom > 1) then
        if (beta < 1.0d-6 .and. nose .lt. zero) then
           rescale = two*tke/real(3*nmobileatom,dp)
        else if (beta .ge. 1.0d-6 .or. nose .ge. zero) then
           rescale = two*tke/real(3*nmobileatom - 3,dp)
        endif
        rescale = sqrt(tempeh/rescale)
     else
        rescale = one
     endif
     ! rescale the velocities
     vxold=vxold*rescale; vyold=vyold*rescale; vzold=vzold*rescale

     ! rescale the kinetic energy as well in hartree
     tke = tke * rescale * rescale
     ! 
     ! If beta is zero we are doing a micro canonical ensemble (and
     ! the system is not coupled to a heat bath at all we use the
     ! Beeman algorithm to integrate the equations of motions). 
     ! If however beta is non zero ie we are in a canonical ensemble
     ! and will have some random force terms.
     !
     ! save the current coordinates into {x,y,z}old
     xold = xcur; yold = ycur; zold = zcur

     if (nose .lt. zero .and. beta < 1.0d-6) then
        do j = 1, natom
           ! calculate the accelarations (in hartree/mass-au)
           accxold(j) = half*clust%force(1,j)/amass(j)*amove(j)
           accyold(j) = half*clust%force(2,j)/amass(j)*amove(j)
           acczold(j) = half*clust%force(3,j)/amass(j)*amove(j)
           ! calculate the new positions based on the velocity and
           ! acceleration as x1=x0+v*t+1/2*a*t^2
           xcur(j) = xold(j) + deltat*(vxold(j) + &
                half*accxold(j)*deltat)
           ycur(j) = yold(j) + deltat*(vyold(j) + &
                half*accyold(j)*deltat)
           zcur(j) = zold(j) + deltat*(vzold(j) + &
                half*acczold(j)*deltat)
           ! also update the velocity based on v1=v0+a*t
           vxold(j) = vxold(j) + accxold(j)*deltat
           vyold(j) = vyold(j) + accyold(j)*deltat
           vzold(j) = vzold(j) + acczold(j)*deltat
        enddo
     else if(nose .lt. zero) then
        do j = 1, natom
           ! Calculate the standard deviation (i.e., the autocorrelation)
           ! sigm = sqrt(two*beta*amass(j)*tempeh/deltat) is the standard
           ! deviation of the force. The devision by deltat comes from
           ! approximating delta(t-t') of continuous time as the square
           ! pulse function: 1/(t-t')[u(t) -u(t')], u() being the Heavyside
           ! square function.
           ! Because we are interested in accelerations and not forces, we
           ! will further divide by amass(j), moving it from numerator to
           ! denominator.
           sigm = sqrt(two*beta*tempeh/(deltat*amass(j)))
           ! For the random acceleration, again gdev is used to get a
           ! standard deviation of /= and then mutiplying by sigm corrects
           ! the standard deviation.
           ! Basically a = F/2*m - beta*v + G/m Again the division by 2 in
           ! the first term F/2m is to take care of Ry to hartree conversion
           accxold(j) = ( half * clust%force(1,j)/amass(j) - &
                beta*vxold(j) + sigm*gdev(mol_dynamic) ) * amove(j)
           accyold(j) = ( half * clust%force(2,j)/amass(j) - &
                beta*vyold(j) + sigm*gdev(mol_dynamic) ) * amove(j)
           acczold(j) = ( half * clust%force(3,j)/amass(j) - &
                beta*vzold(j) + sigm*gdev(mol_dynamic) ) * amove(j)
           ! calculate the first new position : x1=x0+v0*t+1/2*a*t^2
           xcur(j) = xold(j) + deltat*(vxold(j) +  &
                half*deltat*accxold(j))
           ycur(j) = yold(j) + deltat*(vyold(j) +  &
                half*deltat*accyold(j))
           zcur(j) = zold(j) + deltat*(vzold(j) +  &
                half*deltat*acczold(j))
           ! now update the velocities v1=v0+at
           vxold(j) = vxold(j) + deltat*accxold(j)
           vyold(j) = vyold(j) + deltat*accyold(j)
           vzold(j) = vzold(j) + deltat*acczold(j)
        enddo
     else if(nose .ge. zero) then
        !debug 
       !do j = 1, natom
       !  vxold(j) = zero
       !  vyold(j) = zero
       !  vzold(j) = zero
       !enddo
       !tke = zero
        !Nose-Hoover
        ! Num degree of freedom 3(N-1)+1?
        nodof = dble(3*(nmobileatom-1))
        ! 300 K in Hartree
        thkh = 0.02586d0*half/rydberg
        ! Fictitious mass Q
        qnose(1) = nodof*thkh*nose*nose
        if(nchain > 1) then
          do j = 2, nchain
            qnose(j) = qnose(1)/nodof
          enddo
        endif
        enose = tke + dble(nmobileatom)*bdev*half/rydberg
        call nhcmkt(clust,mol_dynamic,tke,vxold,vyold,vzold,&
               mol_dynamic%xnose,mol_dynamic%vnose)
        exn = zero
        evn = zero
        do j = 1, nchain
          evn = evn + half*mol_dynamic%vnose(j)*mol_dynamic%vnose(j)/qnose(j) 
          if(j == 1) then
            exn = exn + nodof*mol_dynamic%xnose(j)*thkh*(tempi/300.d0)
          else
            exn = exn +       mol_dynamic%xnose(j)*thkh*(tempi/300.d0)
          endif
        enddo
        enose = enose + exn + evn
        do j =1, natom
          accxold(j) = amove(j)*(half*clust%force(1,j)/amass(j))
          accyold(j) = amove(j)*(half*clust%force(2,j)/amass(j))
          acczold(j) = amove(j)*(half*clust%force(3,j)/amass(j))
        enddo
       ! Upd velocity and positions
        xold = xcur; yold = ycur; zold = zcur
        do j =1, natom
           vxold(j) = vxold(j) + deltat*half*accxold(j)
           vyold(j) = vyold(j) + deltat*half*accyold(j)
           vzold(j) = vzold(j) + deltat*half*acczold(j)
            xcur(j) =  xold(j) + deltat     * vxold(j)
            ycur(j) =  yold(j) + deltat     * vyold(j)
            zcur(j) =  zold(j) + deltat     * vzold(j)
        enddo
        write(7,*) "Nose-Hoover conserved energy is:", enose
        write(7,*) "Kinetic   energy", tke
        write(7,*) "Potential energy", dble(nmobileatom)*bdev*half/rydberg
        write(7,*) "Xnose     energy", exn
        write(7,*) "Vnose     energy", evn
        write(7,*) "================================"
     endif
     !
     ! print initial coordinates to atom.cor
     write(66,*) 'Coordinates for step #', iframe
     ii = 0
     do i = 1, naty
        write(66,*) clust%name(i)
        do j = 1, clust%natmi(i)
           ii = ii + 1
           write(66,82) xold(ii),yold(ii),zold(ii)
        enddo
     enddo
     ! also, initialize random number generator
     if ( .not. mol_dynamic%use_test_rng) then
     call random_seed
     i = 1
     call random_seed(SIZE=i)
     endif
     !
     ! If restarted MD, or MD has already been used at least once
     ! before, then use the two previous sets of positions,
     ! velocities, accelerations to create new atomic position
     !
  else
     ! again there are two cases here: microcanonical and canonical
     tke = zero
     tke2= zero
     if (nose .lt. zero .and. beta < 1.0d-6) then
        ! In this case, we also have to make sure that we subtarct the
        ! center of mass velocity. Tha just complicates things a little 
        ! bit and some careful handling of old and new coordinates and 
        ! accelerations is required. This will not be that big a problem
        ! because the newer accelerations are just F/m.
        ! To begin with we will not move anything around. Will do that 
        ! as and when needed.

        ! Set accumulators for total mass and center of mass velocity.
        do  j= 1, natom
           ! First predict the velocities
           ! while in the actual formula it should be 
           ! v=(x1-x0)/dt +  dt*(2*F/m+accold)/6 
           ! but since the force is in Ry and has to 
           ! be converted to hartrees the factor of 2 is missing.
           vxold(j) = amove(j) * ( (xcur(j) - xold(j))/deltat + &
                deltat*(clust%force(1,j)/amass(j)+accxold(j))/six )
           vyold(j) = amove(j) * ( (ycur(j) - yold(j))/deltat + &
                deltat*(clust%force(2,j)/amass(j)+accyold(j))/six )
           vzold(j) = amove(j) * ( (zcur(j) - zold(j))/deltat + &
                deltat*(clust%force(3,j)/amass(j)+acczold(j))/six )
        enddo
        mol_dynamic%tke = 0.d0
        do j = 1,natom 
           mol_dynamic%tke = mol_dynamic%tke+half*clust%amass(j)*(vxold(j)*vxold(j)+&
            vyold(j)*vyold(j)+vzold(j)*vzold(j))
        enddo
        tcorrect = -(mol_dynamic%tke+bdev*real(natom,dp)*half/rydberg) + mol_dynamic%initen*half
        if(mol_dynamic%rescale_bo .and. abs(tcorrect) .gt. 1.d-4) then
            write(7,*) "Correcting velocity because the conservation of energy is violated:"
            write(7,*) "Current kinetic energy is", mol_dynamic%tke, "Hartree"
            write(7,*) "Current total energy is",bdev*real(natom,dp)*half/rydberg,"Hartree"
            write(7,*) "Initial total energy is", mol_dynamic%initen*half, "Hartree"
            write(7,*) "Energy difference is",tcorrect,"Hartree"
            write(7,*) ""
            ! Resacle
            vxold(:) = vxold(:) * sqrt((mol_dynamic%tke + tcorrect) / mol_dynamic%tke)
            vyold(:) = vyold(:) * sqrt((mol_dynamic%tke + tcorrect) / mol_dynamic%tke)
            vzold(:) = vzold(:) * sqrt((mol_dynamic%tke + tcorrect) / mol_dynamic%tke)
            !do j = 1, natom
            !   if(vxold(j)*vxold(j) + &
            !     (tcorrect*two/clust%amass(j)) * &
            !     (half*clust%amass(j)*vxold(j)**2)/(mol_dynamic%tke) .lt. 0.d0) then 
            !     vxold(j) = vxold(j)*0.9d0
            !   else
            !     vxold(j) = sign(sqrt(vxold(j)*vxold(j) + & 
            !     (tcorrect*two/clust%amass(j)) * &
            !     (half*clust%amass(j)*vxold(j)**2)/(mol_dynamic%tke)), vxold(j))
            !   endif
            !   if(vyold(j)*vyold(j) + &
            !     (tcorrect*two/clust%amass(j)) * &
            !     (half*clust%amass(j)*vyold(j)**2)/(mol_dynamic%tke) .lt. 0.d0) then 
            !     vyold(j) = vyold(j)*0.9d0
            !   else
            !     vyold(j) = sign(sqrt(vyold(j)*vyold(j) + & 
            !     (tcorrect*two/clust%amass(j)) * &
            !     (half*clust%amass(j)*vyold(j)**2)/(mol_dynamic%tke)), vyold(j)) 
            !   endif
            !   if(vzold(j)*vzold(j) + &
            !     (tcorrect*two/clust%amass(j)) * &
            !     (half*clust%amass(j)*vzold(j)**2)/(mol_dynamic%tke) .lt. 0.d0) then 
            !     vzold(j) = vzold(j)*0.9d0
            !   else
            !     vzold(j) = sign(sqrt(vzold(j)*vzold(j) + & 
            !     (tcorrect*two/clust%amass(j)) * &
            !     (half*clust%amass(j)*vzold(j)**2)/(mol_dynamic%tke)), vzold(j))
            !   endif
            !end do
            mol_dynamic%tke = 0.d0
            do j = 1,natom 
             mol_dynamic%tke = mol_dynamic%tke+half*clust%amass(j)*(vxold(j)*vxold(j)+&
                vyold(j)*vyold(j)+vzold(j)*vzold(j))
            enddo
            write(7,*) "Modified kinetic energy is", mol_dynamic%tke,"Hartree"
            tcorrect = -(mol_dynamic%tke+bdev*real(natom,dp)*half/rydberg) + mol_dynamic%initen*half
            write(7,*) "Corrected energy difference is",tcorrect,"Hartree"
            !mol_dynamic%cstep = 0
        endif
        vxav = dot_product(amass,vxold)
        vyav = dot_product(amass,vyold)
        vzav = dot_product(amass,vzold)
        tmass = sum (amass)
        tmobilemass = dot_product(amass,amove)
        !
        ! find center of mass velocity
        vxav = vxav/tmobilemass
        vyav = vyav/tmobilemass
        vzav = vzav/tmobilemass
        ! correct velocities so that the center of mass does not move and
        vxold = vxold-vxav*amove
        vyold = vyold-vyav*amove
        vzold = vzold-vzav*amove
        ! move the positions to the old
        xold = xcur; yold = ycur; zold = zcur
        write(777,*) iframe, mol_dynamic%mdtime, 0, 0 
        do j = 1, natom

           ! Calculate the updated positions
           ! again the actual expression should have been
           ! x1 = x0 + dt*(v0 + dt*(4*F/m-accold)/6) but again
           ! since the units are hartree we replace 4*F/m by 2*F/m
           xcur(j) = xold(j) + deltat*(vxold(j) + deltat* &
                (two*clust%force(1,j)*amove(j)/amass(j)-accxold(j))/six)
           ycur(j) = yold(j) + deltat*(vyold(j) + deltat* &
                (two*clust%force(2,j)*amove(j)/amass(j)-accyold(j))/six)
           zcur(j) = zold(j) + deltat*(vzold(j) + deltat* &
                (two*clust%force(3,j)*amove(j)/amass(j)-acczold(j))/six)
           ! Find the kinetic energy of the atoms for calculation of
           ! temperature.
           ! Report velocitiy
           write(777,'(3f21.17,i5)') vxold(j), &
                                     vyold(j), &
                                     vzold(j), 0 
           tke = tke + half*amass(j)*(vxold(j)*vxold(j)+vyold(j)* &
                vyold(j)+vzold(j)*vzold(j))
           ! Finally correct the velocity again the force term has the
           ! factor half to take care of rydbergs to hartrees.
           vxold(j) = vxold(j) +  &
                deltat*(three*half*clust%force(1,j)*amove(j)/ &
                amass(j)-accxold(j))*half
           vyold(j) = vyold(j) +  &
                deltat*(three*half*clust%force(2,j)*amove(j)/ &
                amass(j)-accyold(j))*half
           vzold(j) = vzold(j) +  &
                deltat*(three*half*clust%force(3,j)*amove(j)/ &
                amass(j)-acczold(j))*half
           ! And now that we finally dont need acc{x,y,z}old anymore replace
           ! them with F/m again F/2m because of Ry to hartree
           accxold(j) = half*clust%force(1,j)*amove(j)/amass(j)
           accyold(j) = half*clust%force(2,j)*amove(j)/amass(j)
           acczold(j) = half*clust%force(3,j)*amove(j)/amass(j)
        enddo
     else if(nose .lt. zero) then
        ! for all atoms...
        write(777,*) iframe, mol_dynamic%mdtime, 0, 0 
        do j = 1, natom
           ! Calculate the standard deviation (i.e., the autocorrelation)
           ! sigm = sqrt(two*beta*amass(j)*temper/deltat) is the standard
           ! deviation of the force. The devision by deltat comes from
           ! approximating delta(t-t') of continuous time as the square
           ! pulse function: 1/(t-t')[u(t) -u(t')],  u() being the Heavyside
           ! square function.
           ! Because we are interested in accelerations and not forces, we
           ! will further divide by amass(j), moving it from numerator to
           ! denominator.
           sigm = sqrt(two*beta*tempeh/(deltat*amass(j)))
           ! Put the old position coordinates into temporary position
           ! coordinates.
           xoldr = xold(j)
           yoldr = yold(j)
           zoldr = zold(j)
           ! Put the current position coordinates into old /= (We have not
           ! updated the position yet).
           xold(j) = xcur(j)
           yold(j) = ycur(j)
           zold(j) = zcur(j)
           ! Put the older accelerations in the temporary acceleration
           ! holders.
           accxoldr = accxold(j)
           accyoldr = accyold(j)
           acczoldr = acczold(j)
           ! Use old velocity guess in Langevin equation to get new
           ! accelerations -  the half  factor is there because of rydbergs to
           ! hartrees.
           accxold(j) = amove(j) * ( half * clust%force(1,j)/amass(j)- &
                beta * vxold(j) + sigm * gdev(mol_dynamic) )
           accyold(j) = amove(j) * ( half * clust%force(2,j)/amass(j)- &
                beta * vyold(j) + sigm * gdev(mol_dynamic) )
           acczold(j) = amove(j) * ( half * clust%force(3,j)/amass(j)- &
                beta * vzold(j) + sigm * gdev(mol_dynamic) )
           ! predict the velocity
           vxoldr=(xold(j)-xoldr)/deltat + &
                deltat*(two*accxold(j)+ accxoldr)/six
           vyoldr=(yold(j)-yoldr)/deltat + &
                deltat*(two*accyold(j)+ accyoldr)/six
           vzoldr=(zold(j)-zoldr)/deltat + &
                deltat*(two*acczold(j)+ acczoldr)/six
           ! Calculate the new position coordinates 
           ! and put them in xcur.
           xcur(j)=xold(j)+deltat*(vxoldr + deltat*  &
                (four*accxold(j)-accxoldr)/six)
           ycur(j)=yold(j)+deltat*(vyoldr + deltat* &
                (four*accyold(j)-accyoldr)/six)
           zcur(j)=zold(j)+deltat*(vzoldr + deltat* &
                (four*acczold(j)-acczoldr)/six)
           ! Report velocitiy
           write(777,'(3f21.17,i5)') vxold(j), &
                                     vyold(j), &
                                     vzold(j), 0 
           ! calculate the kinetic energy
           tke = tke+half*amass(j)*(vxoldr*vxoldr+vyoldr* &
                vyoldr+vzoldr*vzoldr)
           ! finally correct the velocity
           vxold(j) = vxoldr + deltat* &
                (three*accxold(j)-accxoldr)*half
           vyold(j) = vyoldr + deltat* &
                (three*accyold(j)-accyoldr)*half
           vzold(j) = vzoldr + deltat* &
                (three*acczold(j)-acczoldr)*half
        enddo
     !Nose Hoover
     else if (nose .ge. zero) then
        ! Num degree of freedom 3(N-1)+1?
        nodof = dble(3*(nmobileatom-1))
        ! 300 K in Hartree
        thkh = 0.02586d0*half/rydberg
        ! Fictitious mass Q
        qnose(1) = nodof*thkh*nose*nose
        if(nchain > 1) then
          do j = 2, nchain
          qnose(j) = qnose(1)/nodof
          end do
        endif

        tke =zero
        do j =1, natom
           tke = tke + half*amass(j)*(vxold(j)*vxold(j)+vyold(j)* &
              vyold(j)+vzold(j)*vzold(j))
        enddo

        call nhcmkt(clust,mol_dynamic,tke,vxold,vyold,vzold,&
               mol_dynamic%xnose,mol_dynamic%vnose)
        enose = tke + dble(nmobileatom)*bdev*half/rydberg
        exn = zero
        evn = zero
        do j = 1, nchain
          evn = evn + half*mol_dynamic%vnose(j)*mol_dynamic%vnose(j)/qnose(j) 
          if(j == 1) then
            exn = exn + nodof*mol_dynamic%xnose(j)*thkh*(tempi/300.d0)
          else
            exn = exn +       mol_dynamic%xnose(j)*thkh*(tempi/300.d0)
          endif
        enddo
        enose = enose + evn + exn

        do j = 1, natom
          accxold(j) = amove(j)*(half*clust%force(1,j)/amass(j))
          accyold(j) = amove(j)*(half*clust%force(2,j)/amass(j))
          acczold(j) = amove(j)*(half*clust%force(3,j)/amass(j))
        enddo

        tke =zero
        write(777,*) iframe, mol_dynamic%mdtime, 0, 0 
        do j =1, natom
           vxold(j) = vxold(j) + deltat*half*accxold(j)
           vyold(j) = vyold(j) + deltat*half*accyold(j)
           vzold(j) = vzold(j) + deltat*half*acczold(j)
           tke = tke + half*amass(j)*(vxold(j)*vxold(j)+vyold(j)* &
              vyold(j)+vzold(j)*vzold(j))
           ! Report velocitiy
           write(777,'(3f21.17,i5)') vxold(j), &
                                     vyold(j), &
                                     vzold(j), 0 
        enddo

        xold = xcur; yold = ycur; zold = zcur

        call nhcmkt(clust,mol_dynamic,tke,vxold,vyold,vzold,&
               mol_dynamic%xnose,mol_dynamic%vnose)
        do j =1, natom
           vxold(j) = vxold(j) + deltat*half*accxold(j)
           vyold(j) = vyold(j) + deltat*half*accyold(j)
           vzold(j) = vzold(j) + deltat*half*acczold(j)
            xcur(j) =  xold(j) + deltat     * vxold(j)
            ycur(j) =  yold(j) + deltat     * vyold(j)
            zcur(j) =  zold(j) + deltat     * vzold(j)
        enddo

        write(7,*) "Nose-Hoover conserved energy is:", enose
        write(7,*) "Kinetic   energy", tke
        write(7,*) "Potential energy", dble(nmobileatom)*bdev*half/rydberg
        write(7,*) "Xnose     energy", exn 
        write(7,*) "Vnose     energy", evn
        exn = zero
        evn = zero
        do j = 1, nchain
          evn = evn + half*mol_dynamic%vnose(j)*mol_dynamic%vnose(j)/qnose(j) 
          if(j == 1) then
            exn = exn + nodof*mol_dynamic%xnose(j)*thkh*(tempi/300.d0)
          else
            exn = exn +       mol_dynamic%xnose(j)*thkh*(tempi/300.d0)
          endif
        enddo
        write(7,*) "Xnose(+dt)  energy", exn
        write(7,*) "Vnose(+dt)  energy", evn 
        write(7,*) "================================"
    endif !Nose Hoover 
  endif
  !
  ! If this is a confined system:
  ! If the coordinates are not symmetrized about the origin, shift
  ! them so that they are, in order to allow the results to
  ! converge already for a smaller Rmax.
 
  ! find minimum and maximum coordinate along each axis
  xmin = minval(xcur); ymin = minval(ycur); zmin = minval(zcur)
  xmax = maxval(xcur); ymax = maxval(ycur); zmax = maxval(zcur)

  ! Move the coordinates so that along each axis the farthest atom
  ! in the positive direction is as far as the farthest atom in
  ! the negative direction.
  xmov = half*(xmin+xmax)
  ymov = half*(ymin+ymax)
  zmov = half*(zmin+zmax)
 
  do j=1,natom
     if (pbc%per < 1) xatm(j) = xcur(j)-xmov
     if (pbc%per < 2) yatm(j) = ycur(j)-ymov
     if (pbc%per < 3) zatm(j) = zcur(j)-zmov
  enddo

  !
  ! If this is a periodic system:
  ! It is possible that some of the atoms have "moved outside of
  ! the cell" because of the Langevin dynamics. For the electronic
  ! structure problem, we need to "wrap" them back into the
  ! supercell. The stored trajectories (those printed to
  ! md_mech .dat) DO NOT have this wrapping so that the integration
  ! of the positions remains correct.
  !
  xatm = xcur
  yatm = ycur
  zatm = zcur
  call back_to_cell(pbc,natom,xatm,yatm,zatm)
  !
  ! print coordinates, velocities, accelerations to md_mech.dat
  write(79,24) iframe+1
  do j=1,natom
     write(79,22) xcur(j),ycur(j),zcur(j),xold(j),yold(j),zold(j) &
          ,vxold(j),vyold(j),vzold(j),accxold(j),accyold(j),acczold(j), &
          amove(j)
  enddo
  if(nose .ge. zero) then
    write(7,*)iframe, "For restart run"
    do j = 1, mol_dynamic%nchain
      write(7,'(2f19.14)') mol_dynamic%xnose(j), mol_dynamic%vnose(j)
    enddo
  endif
22 format(13(1x,e15.8))
24 format('Step #',i9)
  !
  ! Just for the purpose of reporting, calculate the center of mass
  ! velocity in the end as well...
  ! Set accumulators for total mass and center of mass velocity.
  vxav = dot_product(amass,vxold)
  vyav = dot_product(amass,vyold)
  vzav = dot_product(amass,vzold)
  tmass = sum(amass)
  tmobilemass = dot_product(amass,amove)
  !
  ! Find center of mass velocity.
  !
  vxav = vxav/tmobilemass
  vyav = vyav/tmobilemass
  vzav = vzav/tmobilemass
   
  ! Calculate total kinetic energy per atom (tke was calculated
  ! earlier) since it is in hartrees we convert it to Ry/atom and eV.
  tke = tke*two/real(nmobileatom,dp)
  tkeev = tke*rydberg

  ! Calculate the effective temperature (as opposed to the input one).
  ! Note here the natom/(natom-1) has to do with the correct
  ! degrees of ereedom! for canonical ensemble. For micro-canonical 
  ! there is no such correction needed.
  if (beta < 1.0d-6 .and. nose .lt. zero) then
     tcalc = two/three*tke/tempkry
  else if (beta .ge. 1.0d-6 .or. nose .ge. zero) then
     tcalc = two/three*tke/tempkry* &
          real(nmobileatom,dp)/real(nmobileatom-1,dp)
  endif

  mol_dynamic%tav = mol_dynamic%tav + tcalc 
  write(7,*) "Average temperature is", mol_dynamic%tav / dble(iframe + 1)

  ! Report the kinetic, potential, and total energies per atom, in
  ! eV, input and actual temperature  to md_nrg.dat
if(.not. mol_dynamic%hybrid .OR. mol_dynamic%seq) then
  write(78,12) iframe,mol_dynamic%mdtime,tkeev,bdev, &
       tkeev+bdev,tempi,tcalc,cdflag
  if(imove /= 0) then
    if(abs(bdev - mol_dynamic%bdevprev) .gt. mol_dynamic%de) then 
      mol_dynamic%de = abs(bdev - mol_dynamic%bdevprev)
    endif
  else if(imove==0) then
    mol_dynamic%de = zero
  endif
    mol_dynamic%bdevprev = bdev
12 format(i4,3x,f9.5,2x,f11.5,2(2x,f11.5),2x,f7.2,3x,f9.4,i4)
  !
  ! Unless this was already the last step...
  if ((imove /= mol_dynamic%step_num) .or.  &
       (abs(tempi-tempf) > 1.d-5) .or. mol_dynamic%hybrid) then
     !
     ! Report coordinates to atom.cor
     write(66,*) 'Coordinates for step #', iframe+1
     ii = 0
     do i = 1, naty
        write(66,*) clust%name(i)
        do j = 1, clust%natmi(i)
           ii = ii + 1
           write(66,82) xatm(ii),yatm(ii),zatm(ii)
        enddo
     enddo

     ! report xyz formatted data to movie.xyz
       write(92,84) natom
       write(92,*) 'Step # ', iframe+1
       ii = 0
       do i = 1, naty
          do j = 1, clust%natmi(i)
             ii = ii + 1
             write(92,88) clust%name(i), xatm(ii)*angs, &
                  yatm(ii)*angs, zatm(ii)*angs
          enddo
       enddo
       if(pbc%per==0) then
         write(931,*) 'ATOM', mol_dynamic%iframe+2
       elseif(pbc%per==3) then
         write(931,*) 'PRIMCOORD', mol_dynamic%iframe+2
         write(931,'(2i6)') natom, 1
       endif
       ii = 0
       do i = 1, clust%type_num
          do j = 1, clust%natmi(i)
             ii = ii + 1
             write(931,'(1i5,3f15.10)') clust%natmi(i), xatm(ii)*angs,yatm(ii)*angs,zatm(ii)*angs
          enddo
       enddo

     ! If this was the last step, report coordinates for next run to atom.cor.
  else
     write(66,*)
     write(66,*) ' Coordinates for next Run '
     do j = 1,natom
        write(66,80) xatm(j),yatm(j),zatm(j)
     enddo
     open(221,file='mdinit.dat',status='unknown',form='formatted')
     write(221,24) iframe+1
     do j=1,natom
        write(221,22) xcur(j),ycur(j),zcur(j),xold(j),yold(j),zold(j) &
             ,vxold(j),vyold(j),vzold(j),accxold(j),accyold(j),acczold(j), &
             amove(j)
     enddo
     do i = 1,naty
        do j = 1,naty
           write(221,*) mol_dynamic%dist_min(i,j), mol_dynamic%dist_max(i,j)
        enddo
     enddo
     write(221,*) bdev
     write(7,*) "writing the maximum delta E in deltat t", mol_dynamic%de, "(meV/atom)"
     write(221,*) mol_dynamic%de
     write(221,*) mol_dynamic%initen
     close(221)
     !NOSE
     if(nose .ge. zero) then
       open(222,file='nhinit.dat',status='unknown',form='formatted')
       do j = 1, mol_dynamic%nchain
         write(222,'(2f19.14)') mol_dynamic%xnose(j), mol_dynamic%vnose(j)
       enddo
       close(222)
     endif
  endif
endif
41 format(A4,i4,2x,3(f11.6,1x),2x,3(f11.6,1x))
80 format(3(2x,f11.6))
82 format(3(2x,f14.9))
84 format(i3)
88 format(a2, 3(3x, f10.6))
  !
  ! Update the temperature (in K) for the next run, depending on
  ! the MD scheme.

  select case (mol_dynamic%cool_type)
     ! "stairstep" annealing
  case(MD_STAIR)
     ! if at last step per the given "stair"
     if (mod(imove+1,mol_dynamic%stride) == 0) then
        ! Calculate the temperature for the next stair - 
        ! move according to fixed step, unless close to 
        ! final temperature, in which case clamp the temperature 
        ! to the final one.
        if (abs(tempi-tempf) > 1.d-5) then
           tempi = tempi + mol_dynamic%tscale
           !
           ! Do not let current temperature go past the final temperature.
           if ( (tempi-tempf)*mol_dynamic%tscale > zero ) tempi = tempf

        endif
     endif
     ! linear annealing
  case (MD_LINEAR)
     tempi = tempi + mol_dynamic%tscale
     ! logrithmic annealing
  case (MD_LOG)
     tempi = tempi * mol_dynamic%tscale
  end select
  !
  ! update frame and step number for next movement
  iframe = iframe + 1
  imove = imove + 1

  if(mol_dynamic%hybrid .and. (.not.mol_dynamic%seq)) then
    continue
  else 
    mol_dynamic%mdtime = mol_dynamic%mdtime + deltat*timeaups 
  endif
! if (imove .gt. mol_dynamic%step_num) then
!    xatm = xold
!    yatm = yold
!    zatm = zold
!    vxold = vxold2
!    vyold = vyold2
!    vzold = vzold2
!  accxold = accxold2
!  accyold = accyold2
!  acczold = acczold2
! endif

  clust%xatm = xatm
  clust%yatm = yatm
  clust%zatm = zatm

  mol_dynamic%tempi = tempi
  mol_dynamic%tempf = tempf
  mol_dynamic%iframe = iframe

  mol_dynamic%xcur = xcur
  mol_dynamic%ycur = ycur
  mol_dynamic%zcur = zcur

  mol_dynamic%xold = xold
  mol_dynamic%yold = yold
  mol_dynamic%zold = zold

  mol_dynamic%vxold = vxold
  mol_dynamic%vyold = vyold
  mol_dynamic%vzold = vzold

  mol_dynamic%accxold = accxold
  mol_dynamic%accyold = accyold
  mol_dynamic%acczold = acczold

  mol_dynamic%nosez    = nosez
  mol_dynamic%nosezd   = nosezd
  mol_dynamic%nosezdd  = nosezdd
  mol_dynamic%nosezcur = nosezcur

  if (imove .gt. mol_dynamic%step_num .and. mol_dynamic%seq) then
    mol_dynamic%xold2= mol_dynamic%xold
    mol_dynamic%yold2= mol_dynamic%yold
    mol_dynamic%zold2= mol_dynamic%zold
    mol_dynamic%xcur2= mol_dynamic%xcur
    mol_dynamic%ycur2= mol_dynamic%ycur
    mol_dynamic%zcur2= mol_dynamic%zcur
    mol_dynamic%vxold2= mol_dynamic%vxold
    mol_dynamic%vyold2= mol_dynamic%vyold
    mol_dynamic%vzold2= mol_dynamic%vzold
    mol_dynamic%accxold2= mol_dynamic%accxold
    mol_dynamic%accyold2= mol_dynamic%accyold
    mol_dynamic%acczold2= mol_dynamic%acczold
  endif

  write(7,*) xold
  write(7,*) xcur
  write(7,*) vxold
  write(7,*) accxold

  ! flush buffer of md output files
  call myflush(66)
  call myflush(78)
  call myflush(79)
  call myflush(92)
 
end subroutine moldyn
!===============================================================
!
! This function selects one random number between 0 and 1 from a
! given set and generates a random Gaussian variable with a mean
! of zero and a standard deviation of one.
! On output, mol_dynamic%iset is incremented by one.
! When random set is exhausted (mol_dynamic%iset = mol_dynamic%nrandom),
! update seed and create a new set.
!
!---------------------------------------------------------------
function gdev(mol_dynamic)

  use constants
  use molecular_dynamic_module
  implicit none
  !
  ! Input/Output variables:
  !
  ! molecular dynamic related data
  type (molecular_dynamic), intent(inout) :: mol_dynamic
  ! selected random number
  real(dp) :: gdev
  !
  ! Work variables:
  !
  ! two random uniform variables between -1 and +1, and the sum of
  ! their squares
  real(dp) :: v1, v2, r
  ! Two random Guassian variable with 0 mean and unity standard
  ! deviation
  real(dp) :: gset
  ! temporary variable
  real(dp) :: fac
  !
  ! External functions:
  !
  real(dp), external :: getrandom
  real(dp), external :: get_rand_test
  !---------------------------------------------------------------

  ! Create two random variables distributed uniformly between -1
  ! and +1, from two calls of getrandom, which returns a uniform
  ! variable between 0 and 1.
  r = 1.1d0
  if (mol_dynamic%use_test_rng) then
      do while (r >= 1.d0)
      v1 = two * get_rand_test(mol_dynamic) - one
      v2 = two * get_rand_test(mol_dynamic) - one
      ! in the (v1,v2) plane, is the results in the unit circle?
      ! If not, try again.
      r= v1 * v1 + v2 * v2
      ! throw out small numbers
      if ( r < 1d-7) r = 1.1d0
      end do
  else
      do while (r >= 1.d0)
      v1 = two * getrandom(mol_dynamic) - one
  v2 = two * getrandom(mol_dynamic) - one
  ! in the (v1,v2) plane, is the results in the unit circle?
  ! If not, try again.
  r= v1 * v1 + v2 * v2
      end do
  endif
  ! Use the Box-muller transformation to create two independent
  ! Gaussian variables - gset and gdev - from the two uniform ones.
  ! gdev will be returned.
  fac = sqrt( -two * log(r) / r )
  gset = v1 * fac
  gdev = v2 * fac

end function gdev
!===============================================================
!
! Picks element of order "iset" from an array of random numbers,
! "vrandom".
! On output, iset is incremented by one.
! When random set is exhausted (iset = nrandom), update seed and
! create a /= set. In order to avoid the unlikely possibility
! of two successive calls to system clock in the same milisecond
! creating the same seed, accumulate system clock to seed.
!
!---------------------------------------------------------------
function getrandom(mol_dynamic)

  use constants
  use molecular_dynamic_module
  implicit none
  !
  ! Input/Output variables:
  !
  ! molecular dynamic related data
  type (molecular_dynamic), intent(inout) :: mol_dynamic
  ! selected random number
  real(dp) :: getrandom
  !
  ! Work variables:
  !
  ! seed (updated if vrandom is recreated)
  integer :: iseed, nr, ii
  integer, allocatable :: oldseed(:)
  !---------------------------------------------------------------
  if (mol_dynamic%iset == 0) then
     call random_seed(SIZE=nr)
     allocate(oldseed(nr))
     if(mol_dynamic%rng_seed > 0) then
         oldseed = mol_dynamic%rng_seed
     else
         call system_clock(COUNT=iseed)
     call random_seed(GET=oldseed)
     oldseed = oldseed + iseed
     endif

     call random_seed(PUT=oldseed)
     call random_number(mol_dynamic%vrandom)
     mol_dynamic%iset = 1
     deallocate(oldseed)
  endif
  getrandom = mol_dynamic%vrandom(mol_dynamic%iset)
  !write(7,*) "random number is",getrandom
  mol_dynamic%iset = mol_dynamic%iset + 1
  if (mol_dynamic%iset > mol_dynamic%nrandom) mol_dynamic%iset = 0

end function getrandom
!===============================================================
! Portable random number generator - Jose Martins version
! based on ran0 from NR2
! For testing purposes only! seems to be only good for 50k calls
function get_rand_test(mol_dynamic)
    use molecular_dynamic_module
    implicit none
  ! molecular dynamic related data
    type (molecular_dynamic), intent(inout) :: mol_dynamic
    real(dp) :: get_rand_test
    integer  :: iseed

    integer ia, im ,iq ,ir, mask ,kk
    parameter (ia=16807,im=2147483647,iq=127773,ir=2836)
    parameter (mask=123459876)

    iseed = mol_dynamic%rng_seed
    if (iseed > im) iseed = mod(iseed,im)
    if (iseed == 0) iseed = mask

    kk = iseed/iq
    iseed = ia*(iseed-kk*iq) -ir*kk
    if (iseed < 0) iseed = iseed + im
    mol_dynamic%rng_seed = iseed
    get_rand_test = dfloat(iseed)/dfloat(im)

    return
end function get_rand_test
