!===============================================================
!
! Copyright (C) 2005 Finite Difference Research Group
! This file is part of parsec, http://www.ices.utexas.edu/parsec/
!
! This subroutine updates thermostat positions, velocities in 
! Nose-Hoover thermostat. 
!
!---------------------------------------------------------------
subroutine nhcmkt(clust,mol_dynamic,tke,vx,vy,vz,xnose,vnose)

  use constants
  use cluster_module
  use molecular_dynamic_module
  implicit none
  !
  ! Input/Output variables:
  !
  ! the cluster
  type(cluster), intent(inout) :: clust
  ! molecular dynamic related data
  type (molecular_dynamic), intent(inout) :: mol_dynamic
  !
  real(dp), intent(in) :: tke
  real(dp), intent(inout) :: vx(clust%atom_num),vy(clust%atom_num),vz(clust%atom_num)
  real(dp), intent(inout) :: vnose(mol_dynamic%nchain),xnose(mol_dynamic%nchain)
  ! Work variables:
  !
  ! atomic mass of each ion type
  real(dp) :: amass(clust%atom_num)

  ! total actual number of atoms
  integer natom

  ! current and final temperatures (in kelvin)
  real(dp) :: tempi
  ! time step
  real(dp) :: deltat

  ! number of mobile atoms (from clust%mvat array)
  integer :: nmobileatom

  ! counters
  integer i,j,nm
  !
  !Nose-Hoover Thermostat stuff
  real(dp), allocatable :: wj(:)
  real(dp) :: qnose(mol_dynamic%nchain), gnose(mol_dynamic%nchain)
  real(dp) :: qscale, dts,thkh,nodof
  real(dp) :: AA
  integer :: nc, nys
  integer :: k,nchain
  !---------------------------------------------------------------

  natom = clust%atom_num
  amass = clust%amass
  nchain = mol_dynamic%nchain

  tempi = mol_dynamic%tempi
  deltat =  mol_dynamic%time_step

  nmobileatom = sum (clust%mvat)

  !Updating
  write(7,*) "Using Nose-Hoover chain thermostat with Martina integrator"
  write(7,*) "with the number of chain is",nchain
  ! Num degree of freedom
  nodof = dble(3*(nmobileatom-1))
  ! 300 K in Hartree
  thkh = 0.02586d0*half/rydberg
  ! Fictitious mass Q
  qnose(1) = nodof*thkh*mol_dynamic%nose*mol_dynamic%nose
  if(mol_dynamic%nchain .ge. 2) then
    do j = 2, nchain
      qnose(j) = qnose(1) / nodof 
    enddo
  endif
  ! Initialize scale
  qscale = 1.d0

  nc = 1
  nys = 5
  allocate(wj(nys))
  if(nys == 5) then
    wj(1) = one / (four - four**third)
    wj(2) = wj(1) 
    wj(4) = wj(1) 
    wj(5) = wj(1) 
    wj(3) = one - four * wj(1)
  endif

  gnose(1) = (two*tke-nodof*thkh*(tempi/300.d0))/qnose(1)
  do k = 1, nc
    do j = 1, nys
      dts = wj(j) * deltat / dble(nc)
      if(nchain==1) then
        vnose(1) = vnose(1) + dts*gnose(1)/four 
      else
        gnose(nchain) = (qnose(nchain-1)*vnose(nchain-1)*vnose(nchain-1)-thkh*tempi/300.d0)/qnose(nchain)
        vnose(nchain) = vnose(nchain) + dts*gnose(nchain)/four 
        do nm = 1, nchain-1
          AA = exp(-dts*vnose(nchain-nm+1)/eight)
          vnose(nchain-nm) = vnose(nchain-nm)*AA
          if(nchain-nm .ne. 1) gnose(nchain-nm) = (qnose(nchain-nm-1)*vnose(nchain-nm-1)*vnose(nchain-nm-1)-&
                                                   thkh*(tempi/300.d0))/qnose(nchain-nm)
          vnose(nchain-nm) = vnose(nchain-nm) + dts*gnose(nchain-nm)/four 
          vnose(nchain-nm) = vnose(nchain-nm)*AA 
        enddo
      endif

      AA = exp(-dts*vnose(1)/two)
      qscale = qscale*AA 
      do nm = 1, nchain
        xnose(nm) = xnose(nm) + half*dts*vnose(nm)
      enddo

      if(nchain == 1) then 
        gnose(1) = (two*qscale*qscale*tke - nodof*thkh*(tempi/300.d0) )/qnose(1)
        vnose(1) = vnose(1) + dts*gnose(1)/four
      else
        do nm = 1, nchain-1
          AA = exp(-dts*vnose(nm+1)/eight)
          vnose(nm) = vnose(nm)*AA
          if(nm == 1) then
            gnose(1) = (two*qscale*qscale*tke - nodof*thkh*(tempi/300.d0) )/qnose(1)
          else
            gnose(nm)= (qnose(nm-1)*vnose(nm-1)*vnose(nm-1)-&
                                          thkh*(tempi/300.d0))/qnose(nm)
          endif
          vnose(nm) = vnose(nm) + dts*gnose(nm)/four
          vnose(nm) = vnose(nm)*AA
        end do
        gnose(nchain) = (qnose(nchain-1)*vnose(nchain-1)*vnose(nchain-1)-&
                                    thkh*(tempi/300.d0))/qnose(nchain)
        vnose(nchain) = vnose(nchain) + dts*gnose(nchain)/four
      endif
    enddo 
  enddo

  write(7,*) "Scaling Velocity"
  do j = 1, natom
     vx(j) = vx(j)*qscale
     vy(j) = vy(j)*qscale
     vz(j) = vz(j)*qscale
  enddo

end subroutine nhcmkt
