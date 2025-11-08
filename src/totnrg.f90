!===============================================================
!
! Copyright (C) 2005 Finite Difference Research Group
! This file is part of parsec, http://www.ices.utexas.edu/parsec/
!
! This subroutine determines the total energy of the system,
! based on Eq. (17) of Chelikowsky and Louie, Phys. Rev. B 29,
! 3470 (1984). 
!
!---------------------------------------------------------------
subroutine totnrg(elec_st,pot,pbc,parallel,totexc,enuc,hcub,natom,bdev,vdw_flag,newk,ierr)

  use constants
  use electronic_struct_module
  use potential_module
  use pbc_module
  use parallel_data_module
  implicit none
  !
  ! Input/Output variables:
  !
  ! electronic structure
  type (electronic_struct), intent(inout) :: elec_st
  ! potential related data
  type (potential), intent(in) :: pot
  ! periodic boundary condition data
  type (pbc_data), intent(in) :: pbc
  ! parallel computation related data
  type (parallel_data), intent(in) :: parallel
  ! new total exchange-correlation energy
  real(dp), intent(in) :: totexc
  ! total nuclear energy
  real(dp), intent(in) :: enuc
  ! hcub = h**3
  real(dp), intent(in) :: hcub
  ! total actual number of atoms (from all types combined)
  integer, intent(in) :: natom

  ! total energy per atom in eV
  real(dp), intent(out) :: bdev
  !
  logical,intent(in),optional :: vdw_flag
  !
  real(dp), intent(out) :: newk
  ! Work variables:
  !
  ! actual number of grid points
  integer ndim
  ! energy terms due to eigenvalue sum, new Hartree potential, new
  ! exchange- correlation potential, old Hartree-Exchange
  ! -correlation potential, ionic potential
  real(dp) :: eeig,ehnew,excnew,ehxcold,eion, psumvec(1)
  ! electronic part of total energy and difference from previous iteration
  real(dp) :: eele,dele
  ! save old electronic energy (computed below) for next iteration
  real(dp), save :: eold = zero
  ! counters
  integer irp,isp,jj,nspin
  ! total energy in Rydberg
  real(dp) :: etot
  ! kpoint variables
  integer kpnum, kplp
  !
  logical add_vdw
  !
  real(dp), allocatable :: vtemp(:)
  integer :: i
  !Additional energy terms of FDE
  real(dp) :: esh, enadd, etionsel, esionsel
  integer, intent(out) :: ierr
  !--------------------------------------------------------------- 

  ierr = 0

  if(present(vdw_flag)) then
!      write(9,*) 'found the vdw flag (DEBUG)', vdw_flag
      add_vdw = vdw_flag
  else
      add_vdw = .FALSE.
  endif

  kpnum = max(elec_st%nkpt,1)

  ndim = parallel%mydim

  nspin = elec_st%nspin/elec_st%mxwd

  ! Sum over eigenvalues, weighted by the occupation:
  eeig = zero
  do isp = 1, nspin
     do kplp = 1, kpnum
        do irp = 1, elec_st%nrep
           if (elec_st%eig(irp,kplp,isp)%nn == 0) cycle
           eeig = eeig +  dot_product(elec_st%eig(irp,kplp,isp)%occ, &
                elec_st%eig(irp,kplp,isp)%en)*elec_st%kpwt(kplp)
        enddo
     enddo
  enddo

  ! If not spin-polarized, multiply by two to account for two
  ! electrons in each eigenvalue.
  if (elec_st%nspin == 1) eeig = two*eeig
 
  ! Energy terms corresponding to new Hartree, new exchange
  ! -correlation, and old Hartree+exchange-correlation. Each one is
  ! an integral over the grid of the relevant potential weighted by
  ! the charge density. 
  ehnew = dot_product(pot%vhart,elec_st%rho(1:ndim,1))
  psumvec = ehnew
  call psum(psumvec,1,parallel%group_size,parallel%group_comm)
  ehnew = psumvec(1)

  excnew = zero
  ehxcold = zero
  do isp=1, elec_st%nspin
     jj=isp-1+elec_st%nspin
     ehxcold = ehxcold + &
          dot_product( pot%vhxcold(1:ndim,isp),elec_st%rho(1:ndim,jj) )
     excnew = excnew + &
          dot_product( pot%vxc(1:ndim,isp), elec_st%rho(1:ndim,jj) )
  enddo
  
  psumvec = ehxcold
  call psum(psumvec,1,parallel%group_size,parallel%group_comm)
  ehxcold = psumvec(1)

  psumvec = excnew
  call psum(psumvec,1,parallel%group_size,parallel%group_comm)
  excnew = psumvec(1)

  ehnew = ehnew*hcub * real(elec_st%nrep)
  excnew = excnew*hcub* real(elec_st%nrep)
  ehxcold = ehxcold*hcub * real(elec_st%nrep)
 
  ! Energy term due to ionic potential (NOT used explicitly for the
  ! total energy, computed for completeness of energies reported in
  ! parsec.out).
  eion = dot_product (pot%vion,elec_st%rho(1:ndim,1))
  psumvec = eion
  call psum(psumvec,1,parallel%group_size,parallel%group_comm)
  eion = psumvec(1)
  eion = eion*hcub * real(elec_st%nrep)

  ! Compute total energy.
  ! The total electronic energy is - (sum of eigenvalues)-(old hxc
  ! energy) +1/2*(new Hartree term) + (total new xc energy).
  eele = eeig - ehxcold + ehnew/two + totexc
  ! Add contribution from LDA+U potential.
  eele = eele + elec_st%etot_plusu
  ! Compute difference from previous iteration (only the electronic
  ! part changes), per atom.
  dele = (eele - eold)*rydberg/real(natom,dp)
  ! Add nuclear energy to get total energy, in rydbergs.
  etot = eele + enuc
  if (pbc%per == 3) etot = etot + pbc%ealpha

  !Old kinetic energy for FDE run
  newk = eeig - ehxcold - eion
  ! Add dispersion corrections if necessary
  if (add_vdw) then
  etot = etot + elec_st%evdw
  end if
  ! Add frozen density embedding part
  if (pot%fex_is) then
    enadd = dot_product (pot%vnadd,elec_st%rho(1:ndim,1))
    psumvec = enadd
    call psum(psumvec, 1, parallel%group_size, parallel%group_comm)
    enadd = psumvec(1)
    enadd = enadd*hcub * real(elec_st%nrep)

    ! Tip-ion/Sample-electron interaction energy
    allocate (vtemp(parallel%mydim))
    do i = 1, parallel%mydim
       vtemp(i) = pot%vion(i) - pot%vsion(i) - pot%vshart(i)
    end do
    etionsel = dot_product(vtemp, pot%rho0)
    psumvec = etionsel
    call psum(psumvec, 1, parallel%group_size, parallel%group_comm)
    etionsel = psumvec(1)
    etionsel = etionsel * hcub * real(elec_st%nrep)

    !Sample electron-ion energy
    esionsel = dot_product(pot%vsion,pot%rho0)
    psumvec = esionsel
    call psum(psumvec, 1, parallel%group_size, parallel%group_comm)
    esionsel = psumvec(1)
    esionsel = esionsel*hcub * real(elec_st%nrep)

    !Sample Hartree energy
    esh = dot_product(pot%vshart, pot%rho0)
    psumvec = esh
    call psum(psumvec, 1, parallel%group_size, parallel%group_comm)
    esh = psumvec(1)
    esh = esh*hcub * real(elec_st%nrep)

    etot = etot + etionsel + esh/two + esionsel + pot%tnadd - enadd + pot%oldke
  end if
  ! Total energy per atom/ in eV.
  bdev = rydberg*etot/real(natom,dp)

  ! Report results.
  if (parallel%iammaster) then
     write(7,*)
     write(7,40) eeig
     write(7,42) ehnew/two
     write(7,44) excnew
     write(7,46) totexc
     write(7,47) eion 
     write(7,49) enuc
     if (add_vdw) write(7,56) elec_st%evdw

     if (pot%fex_is) then
        write(7, 61) pot%tnadd
        write(7, 62) enadd 
        write(7, 63) etionsel
        write(7, 64) esh / two
        write(7, 65) esionsel
        write(7, 67) pot%oldke
      end if
     ! write(7,48) eeig-ehnew/two+totexc !T+Vion?
     if (elec_st%etot_plusu /= zero) write(7,45) elec_st%etot_plusu
     if (pbc%per == 3) write(7,55) pbc%ealpha

     if (eold /= zero) write(7,50) dele 
     write(7,*)
     write(7,52) etot 
     write(7,54) bdev 
     write(7,*)
     if (isnan(excnew) .or. isnan(totexc)) then
        write (7,'(A)') " ERROR: NaNs were observed when using the PBE functional with a large radius."
        write (7,'(A)') "        If this is the case, try reducing the radius."
        call myflush(7)
        ierr = 600
     end if

!48   format(2x,' QES like? "one electron E"     = ',f20.8,' [Ry]')
40   format(2x,' Eigenvalue Energy             = ',f20.8,' [Ry]')
42   format(2x,' Hartree Energy                = ',f20.8,' [Ry]')
44   format(2x,' Integral_{Vxc*rho}            = ',f20.8,' [Ry]')
46   format(2x,' Exc = Integral{eps_xc*rho}    = ',f20.8,' [Ry]')
47   format(2x,' Electron-Ion energy           = ',f20.8,' [Ry]')
49   format(2x,' Ion-Ion Energy                = ',f20.8,' [Ry]')
45   format(2x,' LDA+U Energy                  = ',f20.8,' [Ry]')
50   format(2x,' (E(new)-E(old))/atom  = ',2x,f20.8,' [eV]')
55   format(2x,' Alpha Energy                  = ',f20.8,' [Ry]')
56   format(2x,' TS-vDW correction (included)  = ',f20.8,' [Ry]')
52   format(2x,' Total Energy = ',2x,f20.8,' [Ry]')
54   format(2x,' Energy/atom  = ',2x,f20.8,' [eV]')
61   format(2x,' Non-additive Kinetic Energy   = ',f20.8,' [Ry]')
62   format(2x,' Integral_{Vnadd*rho}          = ',f20.8,' [Ry]')
63   format(2x,' Tip-ion/Sample-electron Energy= ',f20.8,' [Ry]')
64   format(2x,' Sample Hartree Energy         = ',f20.8,' [Ry]')
65   format(2x,' Sample Electron-Ion Energy    = ',f20.8,' [Ry]')
67   format(2x,' Sample Kinetic Energy         = ',f20.8,' [Ry]')
     elec_st%etot = etot
  endif

  ! Save the electronic energy as the old one for the next run.
  eold = eele

end subroutine totnrg
!===============================================================
