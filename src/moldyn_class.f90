!==============================================================
!
! Copyright (C) 2005 Finite Difference Research Group
! This file is part of parsec, http://www.ices.utexas.edu/parsec/
!
! modified by Y. Sakai (2016)
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
subroutine moldyn_class(clust,mol_dynamic,pbc,imove,bdev)

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
  ! einimum and maximum of coordinates along each axis
  real(dp) :: xmin, xmax, ymin, ymax, zmin, zmax
  ! amount by which to shift xcur,ycur,zcur to center coordinates
  ! about the origin
  real(dp) :: xmov, ymov, zmov
  ! total kinetic energy [Ry] and total kinetic energy per atom [eV].
  real(dp) :: tke, tkeev
  ! effective actual temperature 
  real(dp) :: tcalc
  ! factor with which to rescale the velocities
  ! (ratio of t_input /tcalc)
  real(dp) :: rescale

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
  integer i,j,ii
  !
  ! External functions:
  !
  ! a random Gaussian variable with a mean of zero and a standard 
  ! deviation of one
  real(dp), external :: gdev

  !******************Classical Models*****************
  !EAM parameters
  real(dp), dimension(:,:,:), allocatable :: D_morse, &
                                             a_morse, &
                                             r_morse
  real(dp), dimension(:,:,:), allocatable :: gmt, lmt
  real(dp), dimension(:,:,:), allocatable :: gpot2(:,:,:), gpot3(:,:,:)
  real(dp), dimension(:,:,:), allocatable :: gpota(:,:,:)
  real(dp), dimension(:,:), allocatable :: ftheta 
  real(dp) :: etheta, ftot(3)
  real(dp), dimension(:,:), allocatable :: del_morse, &
                                           epslj, siglj, hsc
  real(dp), allocatable :: a_rho(:),b_rho(:),F_emb(:),gamma_emb(:), n(:), rr(:)
  real(dp) :: f_morse, f_tot, rij, rij2, &
              dxx, dyy, dzz, nsum, dndr, &
              dndrsum, dndrsum2, e_morse, e_F, &
              ddn1, ddn2, delta, h_global, rcut, &
              dx1, sc1, dx2, sc2, xcut, sc, ssdvr, &
              d1, d2, sx, v_morse,psi,dpsi, phir, phia, &
              invc,fcut,dcell 
  real(dp), dimension(:,:), allocatable :: forc,delta_dist
  real(dp), parameter :: dh = 0.0001d0 
  integer :: ia,ja, ni, nj, k, l
  integer :: alcstat
  character :: dummy
  ! temporary storage of older position components
  real(dp), dimension(clust%atom_num) :: xold2, yold2, zold2
  real(dp), dimension(clust%atom_num) :: vxold2, vyold2, vzold2
  real(dp), dimension(clust%atom_num) :: accxold2, accyold2, acczold2
  integer, parameter :: cdflag = 1
  logical :: safeflag
  real(dp):: mindist, tcorrect
  integer :: ppp, kkk, intsf
  real(dp):: dthr, tthr
  !Nose-Hoover Thermostat stuff
  real(dp), allocatable :: Jacob(:,:), deltak(:), h(:), DWORK(:)
  real(dp), allocatable :: vxcur(:), vycur(:), vzcur(:)
  integer,  allocatable :: IPIV(:)
  real(dp) :: qnose(mol_dynamic%nchain)
  real(dp) :: nose, nosez, nosezd, nosezdd, nosezdcur
  real(dp) :: nosezold, nosezdold, enose, thkh, nodof, nosezoldr
  real(dp) :: nosezcur,nosezddold
  real(dp) :: exn, evn
  integer :: info, inr
  integer :: nchain
  !---------------------------------------------------------------

  natom = clust%atom_num
  naty = clust%type_num
  xatm = clust%xatm
  !write(7,*) "xxxxatm",xatm
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


  if (imove == 0 ) then
    if(mol_dynamic%class_restart) then
      write(7,*) "RESTARTING CLASSICAL DYNAMICS"
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
  ! else
  !   write(7,*) "INITIALIZING MD PARAMETERS"
  !   xatm = mol_dynamic%xcur2
  !   yatm = mol_dynamic%ycur2
  !   zatm = mol_dynamic%zcur2
  !   xold = 0.d0 
  !   yold = 0.d0
  !   zold = 0.d0
  !   xcur = mol_dynamic%xcur2
  !   ycur = mol_dynamic%ycur2
  !   zcur = mol_dynamic%zcur2
  !   vxold = 0.d0 
  !   vyold = 0.d0
  !   vzold = 0.d0
  !   accxold = 0.d0 
  !   accyold = 0.d0
  !   acczold = 0.d0
    endif
  endif

 !write(7,*) "xatm",xatm
 !write(7,*) "xold",xold
 !write(7,*) "xcur",xcur
 !write(7,*) "vxold",vxold
 !write(7,*) "accxold",accxold

  if(mol_dynamic%ptype ==2) then
    allocate(D_morse(naty,naty,2))
    allocate(a_morse(naty,naty,2)) 
    allocate(r_morse(naty,naty,2))
  else
    allocate(D_morse(naty,naty,1))
    allocate(a_morse(naty,naty,1)) 
    allocate(r_morse(naty,naty,1))
  endif
  D_morse(:,:,:) = 0.0d0
  a_morse(:,:,:) = 0.0d0
  r_morse(:,:,:) = 0.0d0
    
  allocate(del_morse(naty,naty))
  del_morse(:,:) = 0.0d0
  allocate(epslj(naty,naty))
  epslj(:,:) = 0.0d0
  allocate(siglj(naty,naty))
  siglj(:,:) = 0.0d0
  allocate(hsc(naty,naty))
  hsc(:,:) = 0.0d0
  allocate(a_rho(naty))
  !call alccheck(a_rho,naty,alcstat)
  allocate(b_rho(naty))
  !call alccheck(b_rho,naty,alcstat)
  allocate(F_emb(naty))
  !call alccheck(F_emb,naty,alcstat)
  allocate(gamma_emb(naty))
  !call alccheck(gamma_emb,naty,alcstat)
  allocate(forc(3,natom))
  !call alccheck(forc,3*natom,alcstat)
  forc(:,:) = 0.d0
  allocate(n(natom))
  !call alccheck(n,natom,alcstat)
  allocate(rr(naty))
  !call alccheck(rr,naty,alcstat)
  rr(:) = 0.d0
  allocate(delta_dist(naty,naty))
  !call alccheck(delta_dist,naty**2,alcstat)
  allocate(gmt(naty,naty,naty))
  gmt(:,:,:) = 0.d0
  allocate(lmt(naty,naty,naty))
  lmt(:,:,:) = 0.d0
  allocate(ftheta(3,natom))
  ftheta(:,:) = 0.d0
  if(mol_dynamic%ptype==7) then
    allocate(gpot2(naty,naty,6))
    allocate(gpot3(naty,naty,2))
    allocate(gpota(naty,naty,naty))
  endif


  open(18, file = 'md.pot_end')
  open(765, file = 'fin_coord')
  !call read_eam_param(naty,D_morse,a_morse,r_morse,a_rho,b_rho,&
  !F_emb,gamma_emb,h_global,rcut,rr,18)
  if(mol_dynamic%ptype==1) then
     call read_morse_param(mol_dynamic,naty,D_morse,a_morse,r_morse,rcut,hsc,18)
  elseif(mol_dynamic%ptype==2) then
     call read_doublemorse_param(mol_dynamic,naty,D_morse,a_morse,r_morse,del_morse,rcut,hsc,18)
  elseif(mol_dynamic%ptype==3) then
     call read_ljmorse_param(mol_dynamic,naty,D_morse,a_morse,r_morse,epslj,siglj,del_morse,rcut,hsc,18)
  elseif(mol_dynamic%ptype==4) then
     call read_mtheta_param(mol_dynamic,naty,D_morse,a_morse,r_morse,rcut,hsc,gmt,lmt,18)
  elseif(mol_dynamic%ptype==5) then
     call read_ljmt_param(mol_dynamic,naty,D_morse,a_morse,r_morse,epslj,siglj,rcut,hsc,gmt,lmt,18)
  elseif(mol_dynamic%ptype==6) then
     call read_eam_param(             naty,D_morse,a_morse,r_morse,a_rho,b_rho,F_emb,gamma_emb,hsc,rcut,rr,18)
  elseif(mol_dynamic%ptype==7) then
     call read_sw_param(mol_dynamic,naty,gpot2,gpot3,gpota,rcut,18)
  endif
  close(18)
if(mol_dynamic%verbose==1) then
  if ((imove == 0) .and. (.not. mol_dynamic%is_restart)) then
    write(*,*) "rc", rcut
    do i = 1, naty
      do j = 1, naty
        write(7,*) i, j
        write(7,*) "D", D_morse(i,j,:)
        write(7,*) "a", a_morse(i,j,:)
        write(7,*) "r", r_morse(i,j,:)
        write(7,*) "delta", del_morse(i,j)
        if(mol_dynamic%ptype==3) then
           write(7,*) "epsilon", epslj(i,j)
           write(7,*) "sigma", siglj(i,j)
        endif
        write(7,*) "h", hsc(i,j)
        write(7,*) ""
        write(7,*) "The safe regions for ", i, j
        do ppp = 1, clust%natmi(i)
           do kkk = 1, clust%natmi(j)
              write(7,*) "Atom#", ppp, "and #", kkk
              write(7,*) "are min:", mol_dynamic%dist_min(i,j)
              write(7,*) "and max:", mol_dynamic%dist_max(i,j)
              write(7,*) ""
           enddo
        enddo
      end do
     write(7,*) "alpha",  a_rho(i)
     write(7,*) "beta" ,  b_rho(i)
     write(7,*) "r_rho",  rr(i)
     write(7,*) "F0",     F_emb(i)
     write(7,*) "gamma",  gamma_emb(i) 
    end do
  endif
  if(mol_dynamic%ptype==7) then
    write(7,*) gpot2(:,:,:)
    write(7,*) gpot3(:,:,:)
    write(7,*) gpota(:,:,:)
  endif
endif     
  write(*,*) "---------------------------------------"
  e_morse = 0.0d0
  e_F = 0.0d0
 
if(mol_dynamic%ptype == 6) then
   ! Create density for each atom
  ni = 0
  n(:) = 0.0d0
  do ia=1,naty
    do k=1,clust%natmi(ia) 
      ni=ni+1
      nj=0 
      do ja=1,naty
        do l=1,clust%natmi(ja)
          nj=nj+1
          if(ni/=nj) then
             if(mol_dynamic%verbose==1) write(*,*) "Atom #",ni,"and #",nj
             dxx =xatm(ni)-xatm(nj)
             dyy =yatm(ni)-yatm(nj)
             dzz =zatm(ni)-zatm(nj)
             rij2=(dxx*dxx+dyy*dyy+dzz*dzz)
             rij =sqrt(rij2)
             if(rij.lt.rcut) then
                n(ni)=n(ni)+a_rho(ja)*exp(-b_rho(ja)*(rij-rr(ja)))
             endif
             if(mol_dynamic%verbose==1) write(*,*) "ni",ni,n(ni),ia,ja
          endif
        enddo
      enddo
      if(mol_dynamic%verbose==1) write(*,*) "n of", ni, n(ni)
    enddo
  enddo
  ! EAM energy
  ni = 0
  do ia=1,naty
    do k=1,clust%natmi(ia) 
      ni=ni+1
      e_F=e_F+F_emb(ia)*(one-gamma_emb(ia)*log(n(ni)))*(n(ni)**gamma_emb(ia))
    end do
  end do
  write(*,*) "E_F", e_F
endif

  ! pairwise force and energy 
  delta_dist(:,:) = 0.0d0
  ni = 0
  safeflag=.true.
  do ia=1,naty
    do k=1,clust%natmi(ia) 
      ni=ni+1
      nj=0 
      do ja=1,naty
        do l=1,clust%natmi(ja)
          nj =nj+1
          f_morse = zero
          if (ni/=nj) then
            if(mol_dynamic%verbose==1) write(7,*) "Atom #",ni,"and #",nj
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
            if(mol_dynamic%verbose==1) write(7,*) "The interatomic distance is",rij
            mindist=5000.d0
            dthr = abs((mol_dynamic%dist_max(ia,ja) - mol_dynamic%dist_min(ia,ja)) * mol_dynamic%dist_thr)
           !write(7,*) "Threshold for distance is", dthr
           !write(7,*) mol_dynamic%dist_min(ia,ja) - dthr, mol_dynamic%dist_max(ia,ja) + dthr
            if((rij .ge. (mol_dynamic%dist_min(ia,ja) - dthr)) &
              .and. (rij .le. (mol_dynamic%dist_max(ia,ja) + dthr))) then 
              continue
            else
              if(mindist .gt. abs(rij - mol_dynamic%dist_min(ia,ja))) &
                 mindist   =  abs(rij - mol_dynamic%dist_min(ia,ja))
              if(mindist .gt. abs(rij - mol_dynamic%dist_max(ia,ja))) &
                 mindist   =  abs(rij - mol_dynamic%dist_max(ia,ja))
              safeflag = .false.
            endif
            if(mol_dynamic%sc) then
              sx = (rij - rcut) / hsc(ia,ja)
              psi= (sx**4) / (1.d0 + (sx**4))
              dpsi = ((4.d0*sx**3)/((1.d0+sx**4)**2))/hsc(ia,ja)
            else
              psi = 1.d0
              dpsi = 0.d0
            endif
            if(mol_dynamic%verbose==1) write(*,*) "rij", rij
            if(mol_dynamic%verbose==1) write(*,*) "psi and dpsi", psi, dpsi
            if (rij .lt. rcut) then
              if(mol_dynamic%ptype==1 .OR. mol_dynamic%ptype==4 .OR. mol_dynamic%ptype==6) then
                 v_morse=D_morse(ia,ja,1)*(&
                         exp(-two* a_morse(ia,ja,1)*(rij-r_morse(ia,ja,1))) & 
                         -two*exp(-a_morse(ia,ja,1)*(rij-r_morse(ia,ja,1)))&
                         ) 

                 f_morse=D_morse(ia,ja,1)*(-two*a_morse(ia,ja,1)*exp(-two*a_morse(ia,ja,1) & 
                            *(rij-r_morse(ia,ja,1))) &
                            +two*a_morse(ia,ja,1)*exp(-a_morse(ia,ja,1)  &
                            *(rij-r_morse(ia,ja,1)))) 
                 e_morse = e_morse + psi*v_morse / two
                 if(mol_dynamic%verbose==1) write(*,*) "E_morse", e_morse
                 f_morse = psi*f_morse + dpsi*v_morse 
                 if(mol_dynamic%verbose==1) write(*,*) "f_morse", f_morse

                 forc(1,ni)=forc(1,ni)-dxx*f_morse/rij
                 forc(2,ni)=forc(2,ni)-dyy*f_morse/rij
                 forc(3,ni)=forc(3,ni)-dzz*f_morse/rij
              !Double morse
              elseif(mol_dynamic%ptype==2) then
                 do i = 1, 2
                    v_morse=e_morse+D_morse(ia,ja,i)*(exp(-two*a_morse(ia,ja,i) &
                            *(rij-r_morse(ia,ja,i))) & 
                            -two*exp(-a_morse(ia,ja,i) &
                            *(rij-r_morse(ia,ja,i)))) + del_morse(ia,ja)

                    f_morse=D_morse(ia,ja,i)*(-two*a_morse(ia,ja,i)*exp(-two*a_morse(ia,ja,i) & 
                               *(rij-r_morse(ia,ja,i))) &
                               +two*a_morse(ia,ja,i)*exp(-a_morse(ia,ja,i)  &
                               *(rij-r_morse(ia,ja,i)))) 

                    e_morse = e_morse + psi*v_morse / two
                    if(mol_dynamic%verbose==1) write(*,*) "E_dmorse", e_morse
                    f_morse = psi*f_morse + dpsi*v_morse 
                    if(mol_dynamic%verbose==1) write(*,*) "f_dmorse", f_morse

                 enddo
                    
                 forc(1,ni)=forc(1,ni)-dxx*f_morse/rij
                 forc(2,ni)=forc(2,ni)-dyy*f_morse/rij
                 forc(3,ni)=forc(3,ni)-dzz*f_morse/rij
              !Morse plus LJ potential
              elseif(mol_dynamic%ptype == 3 .or. mol_dynamic%ptype== 5) then
                 if(mol_dynamic%ptype ==5) del_morse(:,:) = 0.d0
                 v_morse = D_morse(ia,ja,1)*(&
                          exp(-two*a_morse(ia,ja,1)*(rij-r_morse(ia,ja,1))) &
                         -two*exp(-a_morse(ia,ja,1)*(rij-r_morse(ia,ja,1))))
                 v_morse = v_morse + 4.d0*epslj(ia,ja)*((siglj(ia,ja)/rij)**12 - (siglj(ia,ja)/rij)**6 ) + del_morse(ia,ja) 

                 e_morse = e_morse + psi*v_morse / two

                 f_morse=D_morse(ia,ja,1)*(-two*a_morse(ia,ja,1)*exp(-two*a_morse(ia,ja,1) & 
                            *(rij-r_morse(ia,ja,1))) &
                            +two*a_morse(ia,ja,1)*exp(-a_morse(ia,ja,1)  &
                            *(rij-r_morse(ia,ja,1)))) 
                 if(mol_dynamic%verbose==1) write(*,*) "F_morse part", f_morse
                 f_morse = f_morse - 24.d0 * epslj(ia,ja) *(2.d0 *((siglj(ia,ja) ** 12) / (rij **13)) &
                                                         -          (siglj(ia,ja) **  6) / (rij **7 ))
                 if(mol_dynamic%verbose==1) write(*,*) "psi * f", psi*f_morse
                 f_morse = psi*f_morse + dpsi*v_morse 
                 if(mol_dynamic%verbose==1) write(*,*) "f_ljmorse", f_morse

                 forc(1,ni)=forc(1,ni)-dxx*f_morse/rij
                 forc(2,ni)=forc(2,ni)-dyy*f_morse/rij
                 forc(3,ni)=forc(3,ni)-dzz*f_morse/rij
              ! stiweb
              elseif(mol_dynamic%ptype==7) then
                if(rij .lt. gpot2(ia,ja,6)) then
                  phir =  gpot2(ia,ja,1)*(rij**(-gpot2(ia,ja,3)))
                  phia = -gpot2(ia,ja,2)*(rij**(-gpot2(ia,ja,4)))
                  invc = one/(rij-gpot2(ia,ja,6))
                  fcut = exp(gpot2(ia,ja,5)*invc)
                  v_morse = (phir + phia) * fcut
                  if(invc .lt. zero) then
                    e_morse = e_morse + v_morse * half
                    f_morse = -v_morse * gpot2(ia,ja,5) * invc * invc - fcut*(gpot2(ia,ja,3)*phir+gpot2(ia,ja,4)*phia)/rij
                  else 
                    v_morse = zero
                    f_morse = zero
                  endif
                  if(mol_dynamic%verbose==1) write(7,*) "vmorse", v_morse

                  forc(1,ni)=forc(1,ni)-dxx*f_morse/rij
                  forc(2,ni)=forc(2,ni)-dyy*f_morse/rij
                  forc(3,ni)=forc(3,ni)-dzz*f_morse/rij
                  if(mol_dynamic%verbose==1) then 
                    write(7,*)-dxx*f_morse/rij
                    write(7,*)-dyy*f_morse/rij
                    write(7,*)-dzz*f_morse/rij
                  endif
                endif
              endif
              if(mol_dynamic%ptype==6) then
                ssdvr = zero
                dndr =-b_rho(ja)*a_rho(ja)*exp(-b_rho(ja)*(rij-rr(ja))) 
                ssdvr=ssdvr-dndr*F_emb(ia)*log(n(ni))*(gamma_emb(ia)**2)*n(ni)**(gamma_emb(ia)-one)

                dndr =-b_rho(ia)*a_rho(ia)*exp(-b_rho(ia)*(rij-rr(ia))) 
                ssdvr=ssdvr-dndr*F_emb(ja)*log(n(nj))*(gamma_emb(ja)**2)*n(nj)**(gamma_emb(ja)-one)
                write(*,*) "ssdvr",ssdvr

                !ssdvr = F_emb(ia)*(one-gamma_emb(ia)*log(n(ni+1d-4)))*(n(ni+1d-4)**gamma_emb(ia)) &
                !       -F_emb(ia)*(one-gamma_emb(ia)*log(n(ni-1d-4)))*(n(ni-1d-4)**gamma_emb(ia)) 
                !write(*,*) "ssdvr",ssdvr

                forc(1,ni)=forc(1,ni)-dxx*ssdvr/rij 
                forc(2,ni)=forc(2,ni)-dyy*ssdvr/rij 
                forc(3,ni)=forc(3,ni)-dzz*ssdvr/rij 
              endif
            end if
          end if
        end do
       end do
      !write(*,*) "************************************************"
     enddo
    !write(*,*) "------------------------------------------------"
  enddo

  intsf = 1
  if(mol_dynamic%ptype==4 .or. mol_dynamic%ptype == 5) then 
    !call get_range(mol_dynamic,clust,1,intsf)
    if(intsf == 0) safeflag = .false.
    call threebody(mol_dynamic,clust,naty,clust%natmi,natom,etheta,ftheta,gmt,lmt)
   
    write(7,*) "Energy from three body term (in Ry)"
    write(7,*) "E2 =", e_morse*two
    write(7,*) "E3 =", etheta*two
    write(7,*) "Force from three body term (in Ry/au)"
    do j = 1,natom
       ftheta(:,j) = ftheta(:,j) * two
       write(7,41)  "A",j,xatm(j),yatm(j),zatm(j), &
                   ftheta(1,j),ftheta(2,j),ftheta(3,j)
    end do
  elseif(mol_dynamic%ptype==7) then
    if(intsf == 0) safeflag = .false.
    call thbsw(pbc,mol_dynamic,clust,naty,clust%natmi,natom,etheta,ftheta,gpot3,gpota)
   
    write(7,*) "Energy from three body term (in Ry)"
    write(7,*) "E2 =", e_morse*two
    write(7,*) "E3 =", etheta*two
    write(7,*) "Force from three body term (in Ry/au)"
    do j = 1,natom
       ftheta(:,j) = ftheta(:,j) * two
       write(7,41)  "A",j,xatm(j),yatm(j),zatm(j), &
                   ftheta(1,j),ftheta(2,j),ftheta(3,j)
    end do
  endif

  do j = 1,natom
     forc(:,j) = forc(:,j) * two
     if(mol_dynamic%ptype==4.OR.mol_dynamic%ptype==5 .OR. mol_dynamic%ptype==7) &
             forc(:,j) = forc(:,j)+ftheta(:,j)
  end do

  if(mol_dynamic%hybrid .and. (imove .ge.mol_dynamic%step_num_class)) then
    continue
  else
    write(7,*) iframe, "Total force (in Ry/au)"
    do j = 1,natom
       write(7,41)  "F",j,xatm(j),yatm(j),zatm(j), &
                   forc(1,j),forc(2,j),forc(3,j)
    end do
  endif
  ftot(:) = zero
  do j = 1, natom
     ftot(1) = ftot(1) + forc(1,j)
     ftot(2) = ftot(2) + forc(2,j)
     ftot(3) = ftot(3) + forc(3,j)
  enddo
  write(7,*) "Total force along x, y, z directions:"
  write(7,'(3f15.11)') ftot(:)

  e_morse = e_morse + mol_dynamic%iso_energy*dble(natom)

  write(*,*) "E_morse_fin (Ry)", e_morse*two
  !write(7,*) "E_F_fin (Ry)", e_F*two
  write(7,*) "E_tot (Ry)= ", (e_morse+e_F+etheta)*two
  write(7,*) "E_tot (eV/atom)= ", (e_morse+e_F+etheta)*two*rydberg/dble(natom)
  write(7,*) "************************************************"
  write(7,*) ""
  mol_dynamic%eclass = (e_morse+e_F+etheta)*two
  
  41 format(A3,i4,2x,3(f11.6,1x),2x,3(f11.6,1x))

  if(imove /= 0 .and. (.not. mol_dynamic%class) .and. (abs(mol_dynamic%bdevprev - mol_dynamic%eclass/real(natom,dp)*rydberg) .gt. mol_dynamic%de)) then
    write(7,*) 'The energy difference from previous step is', & 
      abs(mol_dynamic%bdevprev - mol_dynamic%eclass/real(natom,dp)*rydberg), "eV/atom"
    write(7,*) "larger than the max value in FPMD", mol_dynamic%de, "eV/atom"
    write(7,*) "Finish the classical MD run!"
    imove = mol_dynamic%step_num_class 
    mol_dynamic%contn = .true.
    goto 999
  elseif(imove==0 .and. (.not. mol_dynamic%class) .and. (abs(bdev - mol_dynamic%eclass/real(natom,dp)*rydberg) .gt. mol_dynamic%energy_thr) ) then
    write(7,*) 'The energy difference from the final DFT energy is', & 
      abs(bdev - mol_dynamic%eclass/real(natom,dp)*rydberg), "eV/atom"
    write(7,*) "larger than the threshold value", mol_dynamic%energy_thr, "eV/atom"
    write(7,*) "Finish the classical MD run!"
    imove = mol_dynamic%step_num_class 
    goto 999
  endif
  mol_dynamic%bdevprev = mol_dynamic%eclass/real(natom,dp)*rydberg

  if(mol_dynamic%cstep .gt. mol_dynamic%cnmax) then
    write(7,*) mol_dynamic%cnmax, "steps of sequential CMD has been done"
    write(7,*) "Should be refitted"
    mol_dynamic%cstep = 0
    mol_dynamic%contn = .true. 
    imove = mol_dynamic%step_num_class 
    goto 999
  endif
  if(safeflag == .false.) then
    write(7,*) 'The distance variation is', mindist
    write(7,*) "larger than the threshold value:", mol_dynamic%dist_thr*100.d0, "%"
    write(7,*) "Finish the classical MD run!"
    if(imove /= 0) mol_dynamic%contn = .false.
    imove = mol_dynamic%step_num_class 
    goto 999
 !elseif(abs(mol_dynamic%eclass/real(natom,dp)*rydberg - bdev) .gt. mol_dynamic%energy_thr) then
 !  write(7,*) 'The energy variation is', mol_dynamic%eclass/real(natom,dp)*rydberg- bdev, "eV/atom"
 !  write(7,*) "larger than the threshold value:", mol_dynamic%energy_thr, "eV/atom"
 !  write(7,*) "Finish the classical MD run!"
 !  if(imove /= 0) mol_dynamic%contn = .false.
 !  imove = mol_dynamic%step_num_class 
 !  goto 999
  endif
  if(mol_dynamic%hybrid .and. (imove .ge. mol_dynamic%step_num_class)) then 
    if(imove /= 0) mol_dynamic%contn = .false.
    goto 999 
  endif
  mol_dynamic%cstep = mol_dynamic%cstep + 1

 !if(maxval(delta_dist) .gt. mol_dynamic%dist_thr) then 
 !  write(7,*) "The distance change is larger than the threshold value:"
 !  write(7,*) "Finish the classical MD run!"
 !  imove = mol_dynamic%step_num_class 
 !  goto 999
 !endif

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
  if ((imove == 0) .and.((.not. mol_dynamic%class_restart) .and. (.not. mol_dynamic%seq))) then
     write(7,*) "@@@@@@@@@@@@@Fisrt MD run@@@@@@@@@@@"
     ! set accumulators for total mass and center of mass velocity
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
        vyold(j) = vel*gdev(mol_dynamic)*amove(j)
        vzold(j) = vel*gdev(mol_dynamic)*amove(j)
     enddo

     ! determined total mass and mass of mobile atoms
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

     vxold = vxold-vxav*amove
     vyold = vyold-vyav*amove
     vzold = vzold-vzav*amove 

     write(777,*) iframe, mol_dynamic%mdtime,1,1 
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
        write(777,'(3f24.17,i5)')vxold(j) , &
                                 vyold(j) , &
                                 vzold(j) , 1 
     enddo
     ! Calculate the actual temperature and rescale the velocities so
     ! that you get the input temperature.
     if (nmobileatom > 1) then
        if (beta < 1.0d-6 .and. nose .lt. zero) then
           rescale = two*tke/real(3*nmobileatom,dp)
        else if(beta .ge. 1.0d-6 .or. nose .ge. zero) then
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
           vxold(j) = 0.d0
           vyold(j) = 0.d0
           vzold(j) = 0.d0
           ! calculate the accelarations (in hartree/mass-au)
           accxold(j) = half*forc(1,j)/amass(j)*amove(j)
           accyold(j) = half*forc(2,j)/amass(j)*amove(j)
           acczold(j) = half*forc(3,j)/amass(j)*amove(j)
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
           accxold(j) = ( half * forc(1,j)/amass(j) - &
                beta*vxold(j) + sigm*gdev(mol_dynamic) ) * amove(j)
           accyold(j) = ( half * forc(2,j)/amass(j) - &
                beta*vyold(j) + sigm*gdev(mol_dynamic) ) * amove(j)
           acczold(j) = ( half * forc(3,j)/amass(j) - &
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
        enose = tke + mol_dynamic%eclass / two
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
          accxold(j) = amove(j)*(half*forc(1,j)/amass(j))
          accyold(j) = amove(j)*(half*forc(2,j)/amass(j))
          acczold(j) = amove(j)*(half*forc(3,j)/amass(j))
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
        write(7,*) "Potential energy", mol_dynamic%eclass / two 
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
     if (nose .lt. zero .and. beta < 1.0d-6) then
        ! In this case, we also have to make sure that we subtarct the
        ! center of mass velocity. Tha just complicates things a little 
        ! bit and some careful handling of old and new coordinates and 
        ! accelerations is required. This will not be that big a problem
        ! because the newer accelerations are just F/m.
        ! To begin with we will not move anything around. Will do that 
        ! as and when needed.

        ! Set accumulators for total mass and center of mass velocity.
        write(777,*) iframe, mol_dynamic%mdtime,1,1 
        do  j= 1, natom
           ! First predict the velocities
           ! while in the actual formula it should be 
           ! v=(x1-x0)/dt +  dt*(2*F/m+accold)/6 
           ! but since the force is in Ry and has to 
           ! be converted to hartrees the factor of 2 is missing.
           vxold(j) = amove(j) * ( (xcur(j) - xold(j))/deltat + &
                deltat*(forc(1,j)/amass(j)+accxold(j))/six )
           vyold(j) = amove(j) * ( (ycur(j) - yold(j))/deltat + &
                deltat*(forc(2,j)/amass(j)+accyold(j))/six )
           vzold(j) = amove(j) * ( (zcur(j) - zold(j))/deltat + &
                deltat*(forc(3,j)/amass(j)+acczold(j))/six )
        enddo
        mol_dynamic%tke = 0.d0
        do j = 1,natom
           mol_dynamic%tke = mol_dynamic%tke+half*clust%amass(j)*(vxold(j)*vxold(j)+&
            vyold(j)*vyold(j)+vzold(j)*vzold(j))
        enddo
        !write(7,*) "Correcting velocity if the number of CLMD",mol_dynamic%cstep, "is >= 5"
        !if(mol_dynamic%cstep .ge. 4) then
        tcorrect = -(mol_dynamic%tke+mol_dynamic%eclass*half) + mol_dynamic%initen*half
        if(mol_dynamic%rescale_bo .and. mol_dynamic%hybrid) then
            write(7,*) "Current kinetic energy is", mol_dynamic%tke, "Hartree"
            write(7,*) "Current total energy is",mol_dynamic%eclass / two,"Hartree"
            write(7,*) "Initial total energy is", mol_dynamic%initen*half, "Hartree"
            write(7,*) "Energy difference is",tcorrect,"Hartree"
            write(7,*) ""
            ! Rescale
            vxold(:) = vxold(:) * sqrt((mol_dynamic%tke + tcorrect) / mol_dynamic%tke)
            vyold(:) = vyold(:) * sqrt((mol_dynamic%tke + tcorrect) / mol_dynamic%tke)
            vzold(:) = vzold(:) * sqrt((mol_dynamic%tke + tcorrect) / mol_dynamic%tke)
           ! do j = 1,natom
           !   if(vxold(j)*vxold(j) + &
           !     (tcorrect*two/clust%amass(j)) * &
           !     (half*clust%amass(j)*vxold(j)**2)/(mol_dynamic%tke) .lt. 0.d0) then
           !     vxold(j) = 0.d0
           !   else
           !     vxold(j) = sign(sqrt(vxold(j)*vxold(j) + &
           !     (tcorrect*two/clust%amass(j)) * &
           !     (half*clust%amass(j)*vxold(j)**2)/(mol_dynamic%tke)), vxold(j))
           !   endif
           !   if(vyold(j)*vyold(j) + &
           !     (tcorrect*two/clust%amass(j)) * &
           !     (half*clust%amass(j)*vyold(j)**2)/(mol_dynamic%tke) .lt. 0.d0) then
           !     vyold(j) = 0.d0
           !   else
           !     vyold(j) = sign(sqrt(vyold(j)*vyold(j) + &
           !     (tcorrect*two/clust%amass(j)) * &
           !     (half*clust%amass(j)*vyold(j)**2)/(mol_dynamic%tke)), vyold(j))
           !   endif
           !   if(vzold(j)*vzold(j) + &
           !     (tcorrect*two/clust%amass(j)) * &
           !     (half*clust%amass(j)*vzold(j)**2)/(mol_dynamic%tke) .lt. 0.d0) then
           !     vzold(j) = 0.d0
           !   else
           !     vzold(j) = sign(sqrt(vzold(j)*vzold(j) + &
           !     (tcorrect*two/clust%amass(j)) * &
           !     (half*clust%amass(j)*vzold(j)**2)/(mol_dynamic%tke)), vzold(j))
           !   endif
           ! end do
            
            mol_dynamic%tke = 0.d0
            do j = 1,natom
             mol_dynamic%tke = mol_dynamic%tke+half*clust%amass(j)*(vxold(j)*vxold(j)+&
                vyold(j)*vyold(j)+vzold(j)*vzold(j))
            enddo
            write(7,*) "Modified kinetic energy is", mol_dynamic%tke,"Hartree"
            tcorrect = -(mol_dynamic%tke+mol_dynamic%eclass/two) + mol_dynamic%initen*half
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
        do j = 1, natom
           tke = tke + half*amass(j)*(vxold(j)*vxold(j)+vyold(j)* &
                vyold(j)+vzold(j)*vzold(j))
           write(777,'(3f24.17,i5)')vxold(j) , &
                                    vyold(j) , &
                                    vzold(j) , 1 
           ! Calculate the updated positions
           ! again the actual expression should have been
           ! x1 = x0 + dt*(v0 + dt*(4*F/m-accold)/6) but again
           ! since the units are hartree we replace 4*F/m by 2*F/m
           xcur(j) = xold(j) + deltat*(vxold(j) + deltat* &
                (two*forc(1,j)*amove(j)/amass(j)-accxold(j))/six)
           ycur(j) = yold(j) + deltat*(vyold(j) + deltat* &
                (two*forc(2,j)*amove(j)/amass(j)-accyold(j))/six)
           zcur(j) = zold(j) + deltat*(vzold(j) + deltat* &
                (two*forc(3,j)*amove(j)/amass(j)-acczold(j))/six)
           ! Find the kinetic energy of the atoms for calculation of
           ! temperature.
           ! Finally correct the velocity again the force term has the
           ! factor half to take care of rydbergs to hartrees.
           vxold(j) = vxold(j) +  &
                deltat*(three*half*forc(1,j)*amove(j)/ &
                amass(j)-accxold(j))*half
           vyold(j) = vyold(j) +  &
                deltat*(three*half*forc(2,j)*amove(j)/ &
                amass(j)-accyold(j))*half
           vzold(j) = vzold(j) +  &
                deltat*(three*half*forc(3,j)*amove(j)/ &
                amass(j)-acczold(j))*half
           ! And now that we finally dont need acc{x,y,z}old anymore replace
           ! them with F/m again F/2m because of Ry to hartree
           accxold(j) = half*forc(1,j)*amove(j)/amass(j)
           accyold(j) = half*forc(2,j)*amove(j)/amass(j)
           acczold(j) = half*forc(3,j)*amove(j)/amass(j)
        enddo
     else if (nose .lt. zero)then
        !
        ! for all atoms...
        write(777,*) iframe, mol_dynamic%mdtime,1,1 
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
           accxold(j) = amove(j) * ( half * forc(1,j)/amass(j)- &
                beta * vxold(j) + sigm * gdev(mol_dynamic) )
           accyold(j) = amove(j) * ( half * forc(2,j)/amass(j)- &
                beta * vyold(j) + sigm * gdev(mol_dynamic) )
           acczold(j) = amove(j) * ( half * forc(3,j)/amass(j)- &
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
           ! calculate the kinetic energy
           tke = tke+half*amass(j)*(vxoldr*vxoldr+vyoldr* &
                vyoldr+vzoldr*vzoldr)
           write(777,'(3f24.17,i5)')vxold(j) , &
                                    vyold(j) , &
                                    vzold(j) , 1 
           ! finally correct the velocity
           vxold(j) = vxoldr + deltat* &
                (three*accxold(j)-accxoldr)*half
           vyold(j) = vyoldr + deltat* &
                (three*accyold(j)-accyoldr)*half
           vzold(j) = vzoldr + deltat* &
                (three*acczold(j)-acczoldr)*half
        enddo
     else if (nose .ge. zero) then
        ! Num degree of freedom 3(N-1)
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
        enose = tke + mol_dynamic%eclass /two 
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
          accxold(j) = amove(j)*(half*forc(1,j)/amass(j))
          accyold(j) = amove(j)*(half*forc(2,j)/amass(j))
          acczold(j) = amove(j)*(half*forc(3,j)/amass(j))
        enddo

        write(777,*) iframe, mol_dynamic%mdtime,1,1 
        tke =zero
        do j =1, natom
           vxold(j) = vxold(j) + deltat*half*accxold(j)
           vyold(j) = vyold(j) + deltat*half*accyold(j)
           vzold(j) = vzold(j) + deltat*half*acczold(j)
           tke = tke + half*amass(j)*(vxold(j)*vxold(j)+vyold(j)* &
              vyold(j)+vzold(j)*vzold(j))
           write(777,'(3f24.17,i5)')vxold(j) , &
                                    vyold(j) , &
                                    vzold(j) , 1 
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
        write(7,*) "Potential energy", mol_dynamic%eclass /two 
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
     endif
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
    write(7,*) iframe, "For restart run"
    do j = 1, mol_dynamic%nchain
      write(7,'(2f19.14)') mol_dynamic%xnose(j), mol_dynamic%vnose(j)
    enddo
  endif
22 format(13(1x,e15.8))
24 format('Step #',i4)

  mol_dynamic%tke = tke 

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
  ! degrees of freedom! for canonical ensemble. For micro-canonical 
  ! there is no such correction needed.
  if (beta < 1.0d-6 .and. nose .lt. zero) then
     tcalc = two/three*tke/tempkry
  else if(beta .ge. 1.0d-6 .OR. nose .ge. zero) then
     tcalc = two/three*tke/tempkry* &
          real(nmobileatom,dp)/real(nmobileatom-1,dp)
  endif

  mol_dynamic%tav = mol_dynamic%tav + tcalc 
  write(7,*) "Average temperature is", mol_dynamic%tav / dble(iframe + 1)
  ! Report the kinetic, potential, and total energies per atom, in
  ! eV, input and actual temperature  to md_nrg.dat
  write(78,12) iframe,mol_dynamic%mdtime,tkeev,mol_dynamic%eclass*rydberg/dble(natom), &
       tkeev+mol_dynamic%eclass*rydberg/dble(natom),tempi,tcalc,cdflag
12 format(i4,3x,f9.5,2x,f11.5,2(2x,f11.5),2x,f7.2,3x,f9.4,i4)
  !
  ! Unless this was already the last step...
  if ((imove /= mol_dynamic%step_num_class) .or.  &
       (abs(tempi-tempf) > 1.d-5)) then
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
!    write(66,*)
!    write(66,*) ' Coordinates for next Run '
!    do j = 1,natom
!       write(66,80) xatm(j),yatm(j),zatm(j)
!    enddo
     call write_config(clust,mol_dynamic,pbc,imove,(e_F+e_morse)*two,765)
     close(765)
  endif
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
  mol_dynamic%mdtime = mol_dynamic%mdtime + deltat*timeaups 

  iframe = iframe + 1
999 continue
  imove = imove + 1

  if (imove .gt. mol_dynamic%step_num_class) then
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

 !write(7,*) xold
 !write(7,*) xcur
 !write(7,*) vxold
 !write(7,*) accxold

  ! flush buffer of md output files
  call myflush(66)
  call myflush(78)
  call myflush(79)
  call myflush(92)

  deallocate(forc)
  deallocate(D_morse,a_morse,r_morse,a_rho,b_rho,F_emb,gamma_emb)
  
 
end subroutine moldyn_class
!===============================================================
!
! This subroutine creates a potfit input file: md.config.
! Note: here Hartree units for forces and energies are used!
!
!---------------------------------------------------------------
subroutine write_config(clust,mol_dynamic,pbc,imove,bdev,iunit)

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
  ! unit
  integer, intent(in) ::  iunit
  !
  ! Work variables:
  !
  integer i, j, cnt, nt
  real(dp) fsum(3), fatom(3)
  character cfmt*9
  nt = size(clust%name)
  cfmt = "(A3, xA2)"
  write(cfmt(6:6),'(I1)') nt

  !Write down the structural infromation
  if(mol_dynamic%disfor) then
    write(iunit,'(A3, I5, A3)') "#N ", clust%atom_num, " 0"
  else
    write(iunit,'(A3, I5, A3)') "#N ", clust%atom_num, " 1"
  endif
  write(iunit,cfmt) "#C ",clust%name(:)
  if(pbc%per ==0) then
     write(iunit,'(A3, 3f15.10)') "#X ", pbc%latt_vec(:,1) * 100.d0 
     write(iunit,'(A3, 3f15.10)') "#Y ", pbc%latt_vec(:,2) * 100.d0   
     write(iunit,'(A3, 3f15.10)') "#Z ", pbc%latt_vec(:,3) * 100.d0  
  else
     write(iunit,'(A3, 3f15.10)') "#X ", pbc%latt_vec(:,1) 
     write(iunit,'(A3, 3f15.10)') "#Y ", pbc%latt_vec(:,2)  
     write(iunit,'(A3, 3f15.10)') "#Z ", pbc%latt_vec(:,3) 
  endif
  ! Energy in Hartree/atom
  write(iunit,'(A3, f14.9)') "#E ", bdev / rydberg / two - mol_dynamic%iso_energy 
  write(iunit,'(A3)') "#F "
  !Write the current atomic positions and forces
  fsum = sum( clust%force, dim=2 )
  fatom = fsum/real(clust%atom_num,dp)
  cnt = 1
  do i = 1, clust%type_num
     do  j = 1,clust%natmi(i)
        write(iunit,40) i - 1,clust%xatm(cnt),clust%yatm(cnt),clust%zatm(cnt), &
            (clust%force(:,cnt) - fatom) / two 
     cnt = cnt + 1
     enddo
  enddo

40 format(i4,2x,3(f11.6,1x),2x,3(f11.6,1x))
end subroutine write_config
!===============================================================
!
! This subroutine computes three-body energy and force terms from
! the potfit output file. 
!
!---------------------------------------------------------------
subroutine threebody(mol_dynamic,clust,nty,ntat,nat,etheta,ftheta,gm,lm)
  use constants  
  use molecular_dynamic_module
  use cluster_module
implicit none

  type (molecular_dynamic), intent(in) :: mol_dynamic
  type (cluster), intent(in) :: clust

  integer, intent(in) :: nty,nat 
  integer, intent(in) :: ntat(nty)
  !Gamma and Lambda parameters of mtheta potential
  real(dp), intent(in) :: gm(nty,nty,nty), &
                           lm(nty,nty,nty)
  ! counter
  real(dp), intent(out) :: etheta, ftheta(3,nat)
  integer :: i, j, k
  integer :: ia, ja, ka
  integer :: ni, nj, nk
  real(dp) rij2
  real(dp) dxxij, dyyij, dzzij
  real(dp) dxxik, dyyik, dzzik
  real(dp) dxxjk, dyyjk, dzzjk
  real(dp) ctheta,c0
  real(dp) fjtmp(3), fktmp(3)
  real(dp) tmp1, tmp2, tmp_grad
  
  real(dp) rij, rik, rjk

  real(dp) xatm(nat), yatm(nat), zatm(nat)

  xatm = clust%xatm
  yatm = clust%yatm
  zatm = clust%zatm

  etheta = 0.d0
  ftheta(:,:) = 0.d0

  ni = 0
  do ia = 1, nty
    do i=1,ntat(ia) 
      ni=ni+1
      nj = 0
      do ja = 1, nty
        do j = 1, ntat(ja)
          nj=nj+1
          if (ni/=nj) then
            dxxij =xatm(nj)-xatm(ni)
            dyyij =yatm(nj)-yatm(ni)
            dzzij =zatm(nj)-zatm(ni)
            rij2=(dxxij*dxxij+dyyij*dyyij+dzzij*dzzij)
            rij =sqrt(rij2) 
            if(rij > 6.0d0) goto 111
            nk = 0
            do ka = 1, nty
              do k =1, ntat(ka)
                nk = nk + 1
                if (ni /= nk .and. nj .lt. nk) then
                  write(*,*) "Proccesing Atoms #"
                  write(*,*) ni,nj,nk
                  write(*,*) "Atoms spcs are"
                  write(*,*) ia,ja,ka
                  write(*,*) "Gamma and Lambda:"
                  write(*,*) gm(ia,ja,ka),lm(ia,ja,ka)

                  write(89,*) "Atom comb #"
                  write(89,*) ni,nj,nk 
                  dxxik = xatm(nk)-xatm(ni)
                  dyyik = yatm(nk)-yatm(ni)
                  dzzik = zatm(nk)-zatm(ni)
                  rij2=(dxxik*dxxik+dyyik*dyyik+dzzik*dzzik)
                  rik =sqrt(rij2) 

                  dxxjk =xatm(nk)-xatm(nj)
                  dyyjk =yatm(nk)-yatm(nj)
                  dzzjk =zatm(nk)-zatm(nj)
                  rij2=(dxxjk*dxxjk+dyyjk*dyyjk+dzzjk*dzzjk)
                  rjk =sqrt(rij2) 
                  write(*,*) "rij, rik, rjk are"
                  write(*,*) rij, rik, rjk
    
                  if(rik > 6.0d0) goto 112

                  ctheta = (rij**2 + rik**2 - rjk**2) / (2.d0*rij*rik)
                  c0 = cos(lm(ia,ja,ka))
                  write(*,*) "Cosines: theta and theta0 = ", ctheta, c0 
                  etheta = etheta+half*gm(ia,ja,ka)*(ctheta-c0)*(ctheta-c0)
                  write(*,*) "Etheta = ", half*gm(ia,ja,ka)*(ctheta-c0)*(ctheta-c0)
                  write(89,*) "3body energy ", half*gm(ia,ja,ka)*(ctheta-c0)*(ctheta-c0)
                   
                  tmp_grad = - gm(ia,ja,ka) * (ctheta-c0)

                  tmp1 = tmp_grad/(rij*rik) 
                  tmp2 = tmp_grad*ctheta/(rij*rij) 
                  fjtmp(1) = tmp1*dxxik + tmp2*dxxij
                  fjtmp(2) = tmp1*dyyik + tmp2*dyyij
                  fjtmp(3) = tmp1*dzzik + tmp2*dzzij

                  write(89,*) "Fjx",fjtmp(1)
                  write(89,*) "Fjy",fjtmp(2)
                  write(89,*) "Fjz",fjtmp(3)
    
                  tmp2 = tmp_grad *ctheta/(rik*rik) 
                  fktmp(1) = tmp1*dxxij + tmp2*dxxik
                  fktmp(2) = tmp1*dyyij + tmp2*dyyik
                  fktmp(3) = tmp1*dzzij + tmp2*dzzik

                  write(89,*) "Fkx",fktmp(1)
                  write(89,*) "Fky",fktmp(2)
                  write(89,*) "Fkz",fktmp(3)
    
                  !Force on ni-th atom
                  ftheta(1,ni) = ftheta(1,ni) - fjtmp(1) - fktmp(1) 
                  ftheta(2,ni) = ftheta(2,ni) - fjtmp(2) - fktmp(2)  
                  ftheta(3,ni) = ftheta(3,ni) - fjtmp(3) - fktmp(3)  
    
                  !Force on nj-th atom
                  ftheta(1,nj) = ftheta(1,nj) + fjtmp(1) 
                  ftheta(2,nj) = ftheta(2,nj) + fjtmp(2) 
                  ftheta(3,nj) = ftheta(3,nj) + fjtmp(3)   
    
                  !Force on nk-th atom
                  ftheta(1,nk) = ftheta(1,nk) + fktmp(1)   
                  ftheta(2,nk) = ftheta(2,nk) + fktmp(2)   
                  ftheta(3,nk) = ftheta(3,nk) + fktmp(3) 
                endif
112 continue
              enddo
            enddo
          endif
111 continue
        enddo
      enddo
    enddo
  enddo

end subroutine threebody

!===============================================================
!
! This subroutine computes the distance and angle range scanned 
! in FPMD run.
!
!---------------------------------------------------------------
subroutine get_range(mol_dynamic,clust, flag, safeflag)
  use constants  
  use molecular_dynamic_module
  use cluster_module
implicit none

  type (molecular_dynamic), intent(inout) :: mol_dynamic
  type (cluster), intent(in) :: clust
  integer, intent(in) :: flag
  integer, intent(out) :: safeflag
  integer :: i, j, k,nty
  integer :: ia, ja, ka
  integer :: ni, nj, nk
  real(dp) rij, rik, rjk, rij2
  real(dp) xatm(clust%atom_num), yatm(clust%atom_num), zatm(clust%atom_num)
  real(dp) dxxij, dyyij, dzzij
  real(dp) dxxik, dyyik, dzzik
  real(dp) dxxjk, dyyjk, dzzjk
  real(dp) ctheta, deltat

  xatm = clust%xatm
  yatm = clust%yatm
  zatm = clust%zatm
  nty = clust%type_num

  ni = 0
  do ia = 1, nty
    do i=1, clust%natmi(ia)
      ni=ni+1
      nj = 0
      do ja = 1, nty
        do j = 1, clust%natmi(ja)
          nj=nj+1
          if (ni/=nj) then
            dxxij =xatm(nj)-xatm(ni)
            dyyij =yatm(nj)-yatm(ni)
            dzzij =zatm(nj)-zatm(ni)
            rij2=(dxxij*dxxij+dyyij*dyyij+dzzij*dzzij)
            rij =sqrt(rij2) 
            if(flag == 0) then
              if(sqrt(rij2) .lt. mol_dynamic%dist_min(ia,ja)) mol_dynamic%dist_min(ia,ja)=sqrt(rij2)
              if(sqrt(rij2) .gt. mol_dynamic%dist_max(ia,ja)) mol_dynamic%dist_max(ia,ja)=sqrt(rij2)
            endif
            nk = 0
            if(mol_dynamic%ptype == 4) then
              do ka = 1, nty
                do k =1, clust%natmi(ka)
                  nk = nk + 1
                  if (ni /= nk .and. nj .lt. nk) then
                    dxxik = xatm(nk)-xatm(ni)
                    dyyik = yatm(nk)-yatm(ni)
                    dzzik = zatm(nk)-zatm(ni)
                    rij2=(dxxik*dxxik+dyyik*dyyik+dzzik*dzzik)
                    rik =sqrt(rij2) 

                    dxxjk =xatm(nk)-xatm(nj)
                    dyyjk =yatm(nk)-yatm(nj)
                    dzzjk =zatm(nk)-zatm(nj)
                    rij2=(dxxjk*dxxjk+dyyjk*dyyjk+dzzjk*dzzjk)
                    rjk =sqrt(rij2) 
    
                    deltat = (mol_dynamic%theta_max(ia,ja,ka) - mol_dynamic%theta_min(ia,ja,ka)) * mol_dynamic%dist_thr
                    ctheta = (rij**2 + rik**2 - rjk**2) / (2.d0*rij*rik)
                    if(flag == 0) then
                      if(ctheta .lt. mol_dynamic%theta_min(ia,ja,ka)) mol_dynamic%theta_min(ia,ja,ka)=ctheta
                      if(ctheta .gt. mol_dynamic%theta_max(ia,ja,ka)) mol_dynamic%theta_max(ia,ja,ka)=ctheta
                    else if(flag == 1) then
                      if(ctheta .lt. mol_dynamic%theta_min(ia,ja,ka) - deltat) then 
                        safeflag = 0 
                        write(7,*) "Cos theta is out of the safe range", ctheta 
                      else if(ctheta .gt. mol_dynamic%theta_max(ia,ja,ka) + deltat) then 
                        safeflag = 0 
                        write(7,*) "Cos theta is out of the safe range", ctheta 
                      end if
                    endif
                  endif
                enddo
              enddo
            endif
          endif
        enddo
      enddo
    enddo
  enddo
  do ia = 1, nty
     do ja = ia, nty
        mol_dynamic%dist_min(ja,ia) = mol_dynamic%dist_min(ja,ia) 
     end do
  end do
  do i = 1, nty
    do j = 1, nty
      do k = 1, nty
           mol_dynamic%theta_min(i,k,j) = mol_dynamic%theta_min(i,j,k)
           mol_dynamic%theta_max(i,k,j) = mol_dynamic%theta_max(i,j,k)
      enddo
    enddo
  enddo
end subroutine get_range
!===============================================================
!
! This subroutine computes the distance and angle range scanned 
! in FPMD run.
!
!---------------------------------------------------------------
subroutine def_3bcutoff

end subroutine def_3bcutoff
!===============================================================
!
! This subroutine computes three-body energy and force terms from
! the potfit output file. 
!
!---------------------------------------------------------------
subroutine thbsw(pbc,mol_dynamic,clust,nty,ntat,nat,etheta,ftheta,sw3,swa)
  use constants  
  use molecular_dynamic_module
  use cluster_module
  use pbc_module
implicit none

  type (molecular_dynamic), intent(in) :: mol_dynamic
  type (cluster), intent(in) :: clust
  type (pbc_data), intent(in) :: pbc 

  integer, intent(in) :: nty,nat 
  integer, intent(in) :: ntat(nty)
  !Gamma and Lambda parameters of mtheta potential
  real(dp), intent(in) :: sw3(nty,nty,2), &
                          swa(nty,nty,nty)
  ! counter
  real(dp), intent(out) :: etheta, ftheta(3,nat)
  integer :: i, j, k
  integer :: ia, ja, ka
  integer :: ni, nj, nk
  real(dp) rij2,dcell
  real(dp) dxxij, dyyij, dzzij
  real(dp) dxxik, dyyik, dzzik
  real(dp) dxxjk, dyyjk, dzzjk
  real(dp) ctheta,c0
  real(dp) fjtmp(3), fktmp(3)
  real(dp) tmp1, tmp2, tmp_grad,tmp_grad2,tmpr
  real(dp) exik,exij,fxij,fxik
  
  real(dp) rij, rik, rjk

  real(dp) xatm(nat), yatm(nat), zatm(nat)

  xatm = clust%xatm
  yatm = clust%yatm
  zatm = clust%zatm

  etheta = 0.d0
  ftheta(:,:) = 0.d0

  ni = 0
  do ia = 1, nty
    do i=1,ntat(ia) 
      ni=ni+1
      nj = 0
      do ja = 1, nty
        do j = 1, ntat(ja)
          nj=nj+1
          if (ni/=nj) then
            dxxij =xatm(nj)-xatm(ni)
            dyyij =yatm(nj)-yatm(ni)
            dzzij =zatm(nj)-zatm(ni)
            if(pbc%per==3) then
              dcell =abs(pbc%latt_vec(1,1)+pbc%latt_vec(1,2)+pbc%latt_vec(1,3))
              if(abs(dxxij) > dcell*half) dxxij = dxxij - sign(dcell,dxxij) 
              dcell =abs(pbc%latt_vec(2,1)+pbc%latt_vec(2,2)+pbc%latt_vec(2,3))
              if(abs(dyyij) > dcell*half) dyyij = dyyij - sign(dcell,dyyij) 
              dcell =abs(pbc%latt_vec(3,1)+pbc%latt_vec(3,2)+pbc%latt_vec(3,3))
              if(abs(dzzij) > dcell*half) dzzij = dzzij - sign(dcell,dzzij) 
            endif
            rij2=(dxxij*dxxij+dyyij*dyyij+dzzij*dzzij)
            rij =sqrt(rij2) 
            nk = 0
            do ka = 1, nty
              do k =1, ntat(ka)
                nk = nk + 1
                if (ni /= nk .and. nj .lt. nk) then
                  if(mol_dynamic%verbose == 1) then 
                    write(*,*) "Proccesing Atoms #"
                    write(*,*) ni,nj,nk
                    write(*,*) "Atoms spcs are"
                    write(*,*) ia,ja,ka
                    write(*,*) "Gamma and Lambda:"
                    write(*,*) swa(ia,ja,ka)

                    write(7,*) "Atom comb #"
                    write(7,*) ni,nj,nk 
                  endif
                  dxxik = xatm(nk)-xatm(ni)
                  dyyik = yatm(nk)-yatm(ni)
                  dzzik = zatm(nk)-zatm(ni)
                  if(pbc%per==3) then
                    dcell =abs(pbc%latt_vec(1,1)+pbc%latt_vec(1,2)+pbc%latt_vec(1,3))
                    if(abs(dxxik) > dcell*half) dxxik = dxxik -  sign(dcell,dxxik) 
                    dcell =abs(pbc%latt_vec(2,1)+pbc%latt_vec(2,2)+pbc%latt_vec(2,3))
                    if(abs(dyyik) > dcell*half) dyyik = dyyik -  sign(dcell,dyyik) 
                    dcell =abs(pbc%latt_vec(3,1)+pbc%latt_vec(3,2)+pbc%latt_vec(3,3))
                    if(abs(dzzik) > dcell*half) dzzik = dzzik -  sign(dcell,dzzik) 
                  endif
                  rij2=(dxxik*dxxik+dyyik*dyyik+dzzik*dzzik)
                  rik =sqrt(rij2) 

                  !dxxjk =xatm(nk)-xatm(nj)
                  !dyyjk =yatm(nk)-yatm(nj)
                  !dzzjk =zatm(nk)-zatm(nj)
                  dxxjk = dxxik-dxxij
                  dyyjk = dyyik-dyyij
                  dzzjk = dzzik-dzzij
                  !if(pbc%per==3) then
                  !  dcell =abs(pbc%latt_vec(1,1)+pbc%latt_vec(1,2)+pbc%latt_vec(1,3))
                  !  if(abs(dxxjk) > dcell*half) dxxjk = dxxjk -  sign(dcell,dxxjk) 
                  !  dcell =abs(pbc%latt_vec(2,1)+pbc%latt_vec(2,2)+pbc%latt_vec(2,3))
                  !  if(abs(dyyjk) > dcell*half) dyyjk = dyyjk -  sign(dcell,dyyjk) 
                  !  dcell =abs(pbc%latt_vec(3,1)+pbc%latt_vec(3,2)+pbc%latt_vec(3,3))
                  !  if(abs(dzzjk) > dcell*half) dzzjk = dzzjk -  sign(dcell,dzzjk) 
                  !endif
                  rij2=(dxxjk*dxxjk+dyyjk*dyyjk+dzzjk*dzzjk)
                  rjk =sqrt(rij2) 
                  if(mol_dynamic%verbose == 1) write(7,*) "rij, rik, rjk are"
                  if(mol_dynamic%verbose == 1) write(7,*) rij, rik, rjk
    
                  ctheta = (rij**2 + rik**2 - rjk**2) / (2.d0*rij*rik)
                  if(mol_dynamic%verbose == 1) write(7,*) "Cosines: theta = ", ctheta
                   
                  tmpr =  rij-sw3(ia,ja,2)
                  if(tmpr .lt. -0.01d0*sw3(ia,ja,1)) then
                    exij =  exp(sw3(ia,ja,1)/tmpr)
                    fxij = -exij*sw3(ia,ja,1) / (tmpr*tmpr*rij)
                  else
                    exij = 0.d0
                    fxij = 0.d0
                  endif
                  !write(7,*) "tempr", tmpr 
                  tmpr =  rik-sw3(ia,ka,2)
                  if(tmpr .lt. -0.01d0*sw3(ia,ka,1)) then
                    exik =  exp(sw3(ia,ka,1)/tmpr)
                    fxik = -exik * sw3(ia,ka,1) / (tmpr*tmpr*rik)
                  else
                    exik = 0.d0
                    fxik = 0.d0
                  endif
                  !write(7,*) "tempr", tmpr 
                  if(mol_dynamic%verbose == 1) write(7,*) "etheta", swa(ia,ja,ka)*exij*exik*(ctheta+third)*(ctheta+third)
                  etheta = etheta+swa(ia,ja,ka)*exij*exik*(ctheta+third)*(ctheta+third)
                   
                  tmp_grad = two*swa(ia,ja,ka)*exij*exik*(ctheta+third)
                  tmp_grad2= swa(ia,ja,ka)*(ctheta+third)*(ctheta+third)

                  tmp1 = tmp_grad2*fxij*exik-tmp_grad*ctheta/(rij*rij)
                  tmp2 = tmp_grad/(rij*rik) 

                  fjtmp(1) = tmp1*dxxij + tmp2*dxxik
                  fjtmp(2) = tmp1*dyyij + tmp2*dyyik
                  fjtmp(3) = tmp1*dzzij + tmp2*dzzik

                 if(mol_dynamic%verbose == 1) write(7,*) "Fjx",fjtmp(1)
                 if(mol_dynamic%verbose == 1) write(7,*) "Fjy",fjtmp(2)
                 if(mol_dynamic%verbose == 1) write(7,*) "Fjz",fjtmp(3)
    
                  tmp1 = tmp_grad2*fxik*exij-tmp_grad*ctheta/(rik*rik)
                  fktmp(1) = tmp1*dxxik + tmp2*dxxij
                  fktmp(2) = tmp1*dyyik + tmp2*dyyij
                  fktmp(3) = tmp1*dzzik + tmp2*dzzij

                  if(mol_dynamic%verbose == 1) write(7,*) "Fkx",fktmp(1)
                  if(mol_dynamic%verbose == 1) write(7,*) "Fky",fktmp(2)
                  if(mol_dynamic%verbose == 1) write(7,*) "Fkz",fktmp(3)
    
                  !Force on ni-th atom
                  ftheta(1,ni) = ftheta(1,ni) + fjtmp(1) + fktmp(1) 
                  ftheta(2,ni) = ftheta(2,ni) + fjtmp(2) + fktmp(2)  
                  ftheta(3,ni) = ftheta(3,ni) + fjtmp(3) + fktmp(3)  
    
                  !Force on nj-th atom
                  ftheta(1,nj) = ftheta(1,nj) - fjtmp(1) 
                  ftheta(2,nj) = ftheta(2,nj) - fjtmp(2) 
                  ftheta(3,nj) = ftheta(3,nj) - fjtmp(3) 
    
                  !Force on nk-th atom
                  ftheta(1,nk) = ftheta(1,nk) - fktmp(1)   
                  ftheta(2,nk) = ftheta(2,nk) - fktmp(2)   
                  ftheta(3,nk) = ftheta(3,nk) - fktmp(3) 
                endif
              enddo
            enddo
          endif
        enddo
      enddo
    enddo
  enddo

end subroutine thbsw
