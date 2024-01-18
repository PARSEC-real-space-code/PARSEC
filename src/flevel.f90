!===============================================================
!
! Copyright (C) 2005 Finite Difference Research Group
! This file is part of parsec, http://www.ices.utexas.edu/parsec/
!
! This subroutine determines states occupations, based on an
! input Fermi temperature, or based (if the input Fermi
! temperature is negative) on a user provided file - occup.in.
!
!---------------------------------------------------------------
subroutine flevel(elec_st,latt_vec,nscf,ierr)

!_MLT_HP this is a patch for the HP fortran compiler
  use constants, ONLY : dp,zero,half,one,two,twopi,tempkry,rydberg
  use electronic_struct_module, ONLY : electronic_struct
  use nscf_module
!  use constants
!  use electronic_struct_module
!_MLT_HP
  implicit none
  !
  ! Input/Output variables:
  !
  ! electronic structure
  type (electronic_struct), intent(inout) :: elec_st
  ! unit lattice vectors
  real(dp), intent(in) :: latt_vec(3,3)
  ! non self consistent calculation 
  type(nscf_data) :: nscf
  ! error flag, 420 < ierr < 431
  integer, intent(out) :: ierr
  !
  ! Work variables:
  !
  ! useful parameters:
  ! Numerically negligible difference.
  real(dp), parameter :: small=1.0d-12
  ! Value of beta|E-Ef| for which the occupation is exactly one or
  ! zero given the numerical accuracy of the computation. 
  real(dp), parameter :: difmax = 35.0d0
  ! Temperature for which kT is 1 Rydberg.
  real(dp), parameter :: tryd = one/tempkry
  ! Maximum number of iterations for Fermi level position before giving up.
  integer, parameter :: iterlim = 100

  ! Fermi energy (in Rydbergs).
  real(dp) :: efermi
  ! Total number of electrons, per given guess of E_f.
  real(dp) :: enum, tmp

  ! Total number of electrons, divided by two or as is, for
  ! non-spin-polarized and spin-polarized computation, respectively.
  real(dp) :: xeletemp
  ! Distance of an eigenvalue from Fermi level guess, normalized by
  ! temperature.
  real(dp) :: etemp
  ! Boltzmann and Fermi levels for an eigenvalue, per given guess
  ! of the Fermi level.
  real(dp) :: bol, fd
  ! counters
  integer :: jj, nn, nspin, isp, nnn, nj, kplp, kpnum
  real(dp) :: ktmp(3)

  ! (1/kT), in Rydberg
  real(dp) :: beta
  ! upper and lower bounds for E_f during the iterative search 
  real(dp) :: efl, efh
  ! spin identifier for printing
  character (len=2) :: idsp
  ! eigenvalues (electron states) from all representations
  real(dp), dimension(:,:), allocatable :: eigen
  ! occupation of each eigenvalue, all representations
  real(dp), dimension(:,:), allocatable :: occup
  ! counters for representations
  integer irp,jrep(elec_st%nrep)
  integer, allocatable :: irep(:,:)
  ! maximum number of computed eigenstates
  integer nmax
  ! maximum number of computed eigenstates for each spin
  integer nstate(elec_st%nspin)
  ! counter for magnetic moment
  real(dp) ::  magmom
  ! temporary array for k-point weights
  real(dp), dimension(:,:), allocatable :: weight
  ! temporary array for the sorting of eigenvalues
  integer, dimension(:,:), allocatable :: wtindex

  !---------------------------------------------------------------

  !
  ! Start by merging eigenvalues and occupancy factors from all
  ! representations. They will be stored in arrays eigen, occup.
  !
  kpnum = max(elec_st%nkpt,1)
  nspin = elec_st%nspin/elec_st%mxwd
  nstate = 0
  do isp = 1, nspin
     nstate(isp) = 0
     do kplp = 1, kpnum
        do irp = 1, elec_st%nrep
           nstate(isp) = nstate(isp) + elec_st%eig(irp,kplp,isp)%nec
        enddo
     enddo
  enddo
  nmax = maxval(nstate)
  allocate(irep(nmax,nspin))
  allocate(weight(nmax,nspin))
  allocate(wtindex(nmax,nspin))
  irep(:,:) = 0
  weight(:,:) = zero
  wtindex(:,:) = 0

  allocate(eigen(nmax,nspin))
  allocate(occup(nmax,nspin))
  eigen(:,:) = zero
  occup(:,:) = zero

  ! Get the list of lowest eigenvalues and occupancy factors by
  ! going through all representations. array elec_st%eig%irep has
  ! the ordering of states in ascending order of energy eigenvalue.

  ! If kpoints off then take the presorted list of eigenvectors from
  ! their representations in order. If not, then cram them all into one
  ! list and sort them here.  

  do isp = 1, nspin
     jj = 0
     do kplp = 1, kpnum
        jrep = 0
        nnn = 0
        do irp = 1, elec_st%nrep
           nnn = nnn + elec_st%eig(irp,kplp,isp)%nec
        enddo
        do nj = 1, nnn
           jj = jj + 1
           irp = elec_st%irep(nj,kplp,isp)
           jrep(irp) = jrep(irp) + 1
           nn = jrep(irp)
           irep(jj,isp) = irp
           eigen(jj,isp) = elec_st%eig(irp,kplp,isp)%en(nn)
           occup(jj,isp) = elec_st%eig(irp,kplp,isp)%occ(nn)
           weight(jj,isp) = elec_st%kpwt(kplp)
        enddo
     enddo
  enddo

  do isp = 1, nspin
     do nn = 1, nstate(isp)
        wtindex(nn,isp) = nn
     enddo
     call sortreal('SA',nstate(isp),eigen(1,isp),wtindex(1,isp),one)
  enddo

  ! If we are doing a non-self consistent calculation, efermi is already 
  ! given. So just use the value and go straight to the computation of 
  ! occupations.

  ! Compute beta = (1/kT), in Rydberg. initial guess is the highest
  ! eigenvalue the Fermi energy.

  if (nscf%nscf_on) then
     beta = tryd/elec_st%tfermi
     efermi = elec_st%efermi
     goto 20
  endif

  !
  ! If tfermi is negative, occupations are pre-read from a file -
  ! occup.in and stored in array occ_in. Transfer that to occupancy
  ! array and skip calculation of Fermi energy.
  !
  if (elec_st%tfermi < zero) then
     if(elec_st%nkpt > 0) then
        write (7,*) 'CANNOT READ OCC WHEN KPOINTS ON!'
        ierr = 421
     else
        occup(1:elec_st%nstate,:) = elec_st%occ_in(1:elec_st%nstate,:)
        do isp = 1, nspin
           do kplp = 1, kpnum
              jrep = 0
              nnn = 0
              do irp = 1, elec_st%nrep
                 nnn = nnn + elec_st%eig(irp,kplp,isp)%nec
              enddo
              do irp = 1, elec_st%nrep
                 if (associated(elec_st%eig(irp,kplp,isp)%occ)) &
                      elec_st%eig(irp,kplp,isp)%occ(:) = zero
              enddo
              do nn = 1, nstate(isp)
                 irp = elec_st%irep(nn,kplp,isp)
                 jrep(irp) = jrep(irp) + 1
                 elec_st%eig(irp,kplp,isp)%occ(jrep(irp)) = occup(nn,isp)
              enddo
           enddo
        enddo
        goto 26
     endif
  endif
  ! 
  ! If not spin-polarized, the sum of all weighted occupations
  ! should be equal to HALF the total number of electrons (occup=1
  ! is two electrons). If spin-polarized, the sum of weighted
  ! ocupations from the up and down channels is simply the total
  ! number of electrons.
  !
  xeletemp = elec_st%xele * real(elec_st%nspin,dp) / two

  ! Find largest and smallest eigenvalues. Set them as lower and
  ! upper bound for the Fermi energy.
  efl = minval(eigen)
  efh = maxval(eigen)

  ! Compute beta = (1/kT), in Rydberg. initial guess is the highest
  ! eigenvalue the Fermi energy.

  beta = tryd/elec_st%tfermi
  efermi = efh

  ! Iterative search for the Fermi level position starts here. For
  ! each guess of the Fermi level, the occupation of each
  ! eigenlevel is computed  based on the Fermi-Dirac factor. 
  ! The total number of electrons (let's call it n) is then the sum
  ! of all occupations. If n is equal (within numerical accuracy)
  ! to xeletemp, we're done. If not, our next guess for the Fermi
  ! level is 0.5*(efl+efh), and efh or efl (depending on whether n
  ! was larger or smaller than the true number of electrons,
  ! respectively) is set to the old guess value. 

  do jj = 1, iterlim
     enum = zero
     ! For each eigenvalue...
     do isp = 1, nspin
        do nn = 1, nstate(isp)
           ! Calculate distance from Fermi level guess, normalized by
           ! temperature.
           etemp = (eigen(nn,isp)-efermi)*beta
           ! If far below the Fermi level, the occupation is one
           if (etemp < -difmax) then
              enum = enum+one*weight(wtindex(nn,isp),isp)
              ! If not too far above the Fermi level, calculate
              ! fractional occupation.
           elseif (etemp < difmax) then
              bol = exp(etemp)
              fd = one/(bol+one)
              enum = enum + fd*weight(wtindex(nn,isp),isp)
              ! It far above the Fermi level, no contribution to occupation -
              ! go to next eigenvalue.
           else
              exit
           endif
        enddo
        ! If total number of electrons is right - we're done.
     enddo
     if (abs(xeletemp-enum) < small) goto 20
     ! Update this guess for E_f as lower or upper bound for future
     ! guesses, depending if the number of electrons was too low or
     ! too high, respectively.
     if (enum < xeletemp) then
        efl = efermi
     else
        efh = efermi
     endif
     ! Create new guess for E_f, as explained above.
     efermi = (efh+efl)*half
     ! Iterative search for the Fermi level position ends here.
  enddo

  ! If Fermi level not found even after "iterlim" iterations,
  ! complain and quit.
  write(7,'(/,a)') ' ERROR: Fermi level not found'
  write(7,'(a,f10.4,a)') ' Last Fermi level guess: ',efermi,' [Ry]'
  write(7,*) 'STOP in flevel'
  ierr = 422
  return
20 continue

  ! For the Fermi level found, determine occupancy array.
  ! Define ifmax as the order of the highest energy band (or level) with
  ! non-zero occupancy, for each spin channel.

  do isp = 1, nspin

     do kplp = 1, kpnum
        jrep = 0
        nnn = 0
        do irp = 1, elec_st%nrep
           nnn = nnn + elec_st%eig(irp,kplp,isp)%nec
        enddo
        do nj = 1, nnn
           irp = elec_st%irep(nj,kplp,isp)
           jrep(irp) = jrep(irp) + 1
           nn = jrep(irp)
           etemp = (elec_st%eig(irp,kplp,isp)%en(nn)-efermi)*beta
           if (etemp < -difmax) then
              elec_st%eig(irp,kplp,isp)%occ(nn) = one
           elseif (etemp < difmax) then
              elec_st%eig(irp,kplp,isp)%occ(nn) = one/(exp(etemp)+one)
           else
              elec_st%eig(irp,kplp,isp)%occ(nn) = zero
              cycle
           endif
           elec_st%ifmax(isp) = nj
        enddo
     enddo

     do jj = 1, nstate(isp)
        etemp = (eigen(jj,isp)-efermi)*beta
        if (etemp < -difmax) then
           occup(jj,isp) = one
        elseif (etemp < difmax) then
           occup(jj,isp) = one/(exp(etemp)+one)
        else
           occup(jj,isp) = zero
        endif
     enddo

  enddo

26 continue

  ! If computation not spin-polarized, the total number of
  ! electrons is simply xele.
  if (elec_st%nspin == 1) then
     elec_st%totel(1)=elec_st%xele
     ! If it is spin-polarized, find the number of electrons in each
     ! spin channel.
  else
     if(elec_st%mxwd == 1) then
        do isp = 1, nspin
           elec_st%totel(isp) = dot_product(occup(1:nstate(isp),isp), &
                weight(wtindex(1:nstate(isp),isp),isp))
        enddo
     endif
     ! Report the number of electrons and the magnetic moment.
     if (elec_st%mxwd == 2 .and. (.not. elec_st%ncl)) then
        magmom = zero
        do nn = 1, elec_st%eig(1,1,1)%nn
           do kplp = 1, kpnum
              magmom = magmom + elec_st%magmom(nn,kplp)* &
                   elec_st%eig(1,kplp,1)%occ(nn) * elec_st%kpwt(kplp)                   
           enddo
        enddo
        elec_st%net_magmom = magmom
        do isp = 1,2
           tmp = (-one)**isp
           elec_st%totel(isp) = half*(elec_st%xele - tmp*elec_st%net_magmom)
        enddo
 
        write(7,*)
        if(elec_st%is_so) then
           write(7,*) 'Spin-Orbit: each state has a magnetic moment!'
           write(7,*) '   Spin is not a good quantum number here!'
           write(7,*) '---------------------------------------------'
           write(7,52) magmom
        endif
        write(7,*)
     elseif(elec_st%mxwd == 1) then
        write(7,*)
        write(7,51) elec_st%totel(1), elec_st%totel(2)
        write(7,52) elec_st%totel(1) - elec_st%totel(2)
        write(7,*)
     endif
  endif

51 format(1x,'Electrons in up, down spins:',1x,f7.2,1x,f7.2)
52 format(1x,'Net magnetic moment (in Bohr magnetons):',1x,f7.2)
53 format(1x,'Expectation value of spin |<S>|:',1x,f7.2)
  if (elec_st%tfermi >= zero) elec_st%efermi = efermi

  ! Report Fermi level. 
  if (elec_st%tfermi >= zero) write(7,61) efermi   
  do isp = 1, nspin
     ! Set spin identifier.
     jj = isp - 1 + elec_st%nspin
     if (jj == 1) idsp = '  '
     if (jj == 2) idsp = 'up'
     if (jj == 3) idsp = 'dn'
     write(7,*)
     do kplp = 1,kpnum
        if (elec_st%nkpt /= 0) then
           call matvec3('T',latt_vec,elec_st%kpts(1,kplp),ktmp)
           ktmp = ktmp / twopi
           write(7,'(a,3f9.4,a)') ' k-point = ', &
                elec_st%kpts(:,kplp),' (a.u.)^-1'
           write(7,'(a,3f9.4,a)') ' k-point = ', &
                ktmp,' (recip. latt. vec.) '
        endif
        if (elec_st%mxwd == 2) write(7,66) idsp
           
        if (elec_st%mxwd == 1) write(7,64) idsp
        write(7,*)

        jrep = 0
        nnn = 0
        do irp = 1, elec_st%nrep
           nnn = nnn + elec_st%eig(irp,kplp,isp)%nec
        enddo
        do nn=1,nnn
           irp = irep(wtindex(nn,isp),isp)
           jrep(irp) = jrep(irp) + 1
           if(elec_st%mxwd == 2) then
              write(7,65) nn,elec_st%eig(irp,kplp,isp)%en(jrep(irp)), &
                   rydberg*elec_st%eig(irp,kplp,isp)%en(jrep(irp)), &
                   elec_st%eig(irp,kplp,isp)%occ(jrep(irp)),irp, &
                   elec_st%magmom(nn,kplp)
           else
              write(7,63) nn,elec_st%eig(irp,kplp,isp)%en(jrep(irp)), &
                   rydberg*elec_st%eig(irp,kplp,isp)%en(jrep(irp)), &
                   elec_st%eig(irp,kplp,isp)%occ(jrep(irp)),irp
           endif
        enddo
        write(7,*)
     enddo
  enddo

  call myflush(7)

  deallocate(irep, weight, wtindex, eigen, occup)

61 format(1x,'Fermi level at',1x,f10.4,' [Ry]')
63 format(i5,3x,f18.10,3x,f18.10,1x,f9.4,3x,i6)
64 format(a2,' State   Eigenvalue [Ry]      Eigenvalue [eV]', &
        '    Occup.     Repr.')
65 format(i5,3x,f14.6,3x,f14.6,1x,f9.4,3x,i6,3x,f9.4)
66 format(a2,' State   Eigenv.[Ry]      Eigenv.[eV]', &
        '    Occup.     Repr.   Spin_z [hbar]')

end subroutine flevel
!===============================================================
