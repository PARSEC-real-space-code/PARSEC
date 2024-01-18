!===============================================================
!
!  Copyright (C) 2005 Finite Difference Research Group
!  This file is part of parsec, http://www.ices.utexas.edu/parsec/
!
!  This subroutine saves grid, potential, wave function, energy
!  levels, etc, in parsec.dat file for restarting purposes.
!
!---------------------------------------------------------------
subroutine gw_write(elec_st,clust,grid,pbc,symm,parallel,rho,vxc)

#ifdef GW
  use constants
  use electronic_struct_module
  use cluster_module
  use grid_module
  use pbc_module
  use symmetry_module
  use parallel_data_module
  use potential_module
  use wfn_rho_vxc_io_m, only : write_binary_header
  use gspace_module
#ifdef MPI
  !  include mpi definitions
  use mpi
#endif

  implicit none 
  !
  !  Input/Output variables:
  !
  !  electronic structure
  type (electronic_struct), intent(inout) :: elec_st
  !  cluster structure
  type (cluster), intent(in) :: clust
  !  grid related data
  type (grid_data), intent(in) :: grid
  !  periodic boundary conditions data
  type (pbc_data), intent(in) :: pbc
  !  symmetry operations:
  type (symmetry), intent(in) :: symm
  !  parallel computation related data
  type (parallel_data), intent(inout) :: parallel
  ! distributed electron density and exchange-correlation potential
  ! (passed outside the structure to overcome a bug with the IBM compiler)
  real(dp), intent(in) :: rho(parallel%mydim,2*elec_st%nspin-1)
  real(dp), intent(in) :: vxc(parallel%mydim,elec_st%nspin)
  !
  !  Work variables:
  !
  !  communicator
  integer comm
  !  exit code for mpi calls
  integer mpinfo
  !  counters
  integer isp, jj, kk, ii, dim, g1, g2, g3
  real(dp) :: tmpd

  integer msgtype, nstate
  !  counters for representations
  integer irp
  !  temporary arrays for eigenvalues/occupations, before reshaping
  !  according to symmetry operations
  real(dp), dimension(:), allocatable :: en_tmp, occ_tmp
  !  jrep keeps track of how many eigenstates are already in each
  !  representation
  integer,allocatable :: jrep(:,:), jrepinv(:,:)
  !  work array
  complex(dpc), allocatable :: zwftmp(:)
  real(dp), allocatable :: wftmp(:)
  !  date/time tag
  character (len=26) :: datelabel
  !
  integer stype(10)
  ! Structure for gspace (PARATEC) -- charge density and vxc
  type(gspace) :: pot_gspace
  type(gspace), allocatable :: k_gspace(:)
  complex(dpc), allocatable :: vxc_g_full(:,:), rho_g_full(:)
  complex(dpc), allocatable :: wfn_g_full(:,:)
  ! blas norm
  real(dp), external :: dznrm2

#ifdef MPI
  integer status(MPI_STATUS_SIZE)
#endif

  ! cutoff - wavefunction cutoff
  real(dp) :: cutoff
  !  kpoint variables
  integer ikp, kpnum
  !  weight of each kpoint being examined
  real(dp), dimension(:), pointer :: weight
  !  coordinates of each kpoint
  real(dp), dimension(:,:), pointer :: kpts
  ! FFT size
  integer :: fftsize(3)
  ! Number of g-vectors per k point
  integer, pointer :: ngk(:)
  ! Maximum and minimum number of occupied states at each k point
  integer, pointer :: ifmin(:,:), ifmax(:,:)
  ! Number of atoms of each type
  integer, pointer :: atyp(:)
  ! Positions  of all atoms
  real(dp), pointer :: apos(:,:)
  ! Occupation for each band at each  kpoint for each spin 
  real(dp), pointer :: occupations(:,:,:)
  ! Eigenvalue of each band at each  kpoint for each spin 
  real(dp), pointer ::  energies(:,:,:)
  ! Header 
  character*3 :: sheader
  ! Number of spins
  integer :: nspin
  ! maximum number of g-vectors at any kpoint
  integer :: ngkmax

  ! allocation check
  integer alcstat
  ! exchange-correlation potential collected to master processor
  real(dp), allocatable :: vxc_m(:,:)
  ! total charge density collected to master processor
  real(dp), allocatable :: rho_m(:)
  ! wavefunction on the master processor
  complex(dpc), allocatable :: wfn_m(:)
  real(dp), allocatable     :: dwfn_m(:)

  integer :: i,j,k,iadd,js,p,n,iflag,irk, iadd2, kmax(3)
  real(dp) :: tmpr,tmpi, norm, factor
  complex(dpc) :: ztmp,fi_s
  complex(dpc), dimension(:), allocatable :: rho_tmp, rho_g
  complex(dpc), dimension(:,:), allocatable :: vxc_tmp,vxc_g
  complex(dpc), dimension(:), allocatable :: wfn_tmp


  !---------------------------------------------------------------

  nspin = elec_st%nspin
  if (pbc%is_on) then
     kpnum = elec_st%nkpt
     allocate(kpts(3,kpnum))
     kpts = zero
     do irk = 1, kpnum
       call matvec3('T',pbc%latt_vec,elec_st%kpts(1,irk),kpts(1:3,irk))
       kpts(1:3,irk) = kpts(1:3,irk) / twopi
     enddo
     allocate(weight(kpnum))
     weight = elec_st%kpwt
     cutoff = 1000.0
     do i = 1, 3
        cutoff = min((pi/grid%step(i))**2 - 1.0D-4, cutoff)
     enddo
     ngkmax = 0
     allocate(ngk(1:kpnum))
     ngk(1:kpnum) = 0
  else
     kpnum = 1
     allocate(kpts(3,kpnum))
     kpts = zero
     allocate(weight(kpnum))
     weight = one
     !AJB why stepin?
     cutoff = (pi/grid%stepin)**2 
     fftsize(1) = (grid%nxmax - grid%nxmin +1)
     fftsize(2) = (grid%nymax - grid%nymin +1)
     fftsize(3) = (grid%nzmax - grid%nzmin +1)
     !!!!ngkmax = pbc%ng
     !!!!allocate(ngk(1:kpnum))
     !!!!ngk(1:kpnum) = ngkmax
  endif


  allocate(energies(1:elec_st%nstate, 1:kpnum, 1:nspin))
  allocate(occupations(1:elec_st%nstate, 1:kpnum, 1:nspin))
  allocate(ifmin(1:kpnum, 1:nspin))
  allocate(ifmax(1:kpnum, 1:nspin))
  allocate(jrep(elec_st%nspin,elec_st%nrep))
  allocate(jrepinv(elec_st%nspin,elec_st%nstate))
  ifmin = 0
  ifmax = 0
  energies = zero
  occupations = zero
  jrepinv = 0
  jrep = 0

  do isp = 1, nspin
     do ikp = 1, kpnum
        if (parallel%iammaster) then
           !write(7,('(2i4,3f12.7,i4,f12.7)')) isp, ikp, kpts(1:3,ikp), ngk(ikp), weight(ikp)
           jrep(isp,:) = 0
           do kk = 1, elec_st%nstate
              irp = elec_st%irep(kk,ikp,isp)
              jrep(isp,irp) = jrep(isp,irp) + 1
              energies(kk, ikp, isp) = elec_st%eig(irp,ikp,isp)%en(jrep(isp,irp))
              occupations(kk, ikp,isp) = elec_st%eig(irp,ikp,isp)%occ(jrep(isp,irp))
              if(occupations(kk,ikp,isp).gt.0.5d0) then
                  if (ifmin(ikp, isp) == 0) ifmin(ikp, isp) = kk
                  ifmax(ikp, isp) = kk
              end if
              !write(7,('(3i4,2f12.7)')) irp, jrep(isp,irp), kk, energies(kk,ikp,isp),occupations(kk,ikp,isp)
           enddo
           !write(7,('(2i4)')) ifmin(ikp,isp), ifmax(ikp,isp)
        endif
     enddo
  enddo
 
  allocate(atyp(1:clust%atom_num),stat=alcstat)
  allocate(apos(1:3, 1:clust%atom_num), stat=alcstat)
  n = 0
  do k = 1, clust%type_num
    do i = 1, clust%natmi(k)
      n = n + 1
      atyp(n) = k
      apos(1, n) = clust%xatm(n)
      apos(2, n) = clust%yatm(n)
      apos(3, n) = clust%zatm(n)
    enddo
  enddo
 
  if (parallel%iammaster) then

     write(7,*) '\nAssuming wave-function cutoff of ',cutoff,' Ry \n'

! First generate the g-space needed -- only need number of gvectors
     pot_gspace%imap  = .false.
     pot_gspace%igvec = .false.  
     call generate_potential_gspace(pot_gspace, pbc, symm, cutoff, 1)
     fftsize(1:3) = pot_gspace%fftsize(1:3)
! Allocate the g-space for the kpoints -- only need number of gvectors
     allocate(k_gspace(kpnum))
     do irk = 1, kpnum
        k_gspace(irk)%imap  = .false.
        k_gspace(irk)%igvec = .false.  
        call generate_k_gspaces(k_gspace(irk), ngkmax, irk, fftsize, &
         pbc, kpts(1:3,irk), weight(irk), symm, cutoff, 1)
        ngk(irk) = k_gspace(irk)%length
     enddo
     iflag = 2
     sheader = 'VXC'
     open(21, file = 'VXC', status = 'replace', form = 'unformatted')
     call write_binary_header(21,sheader,iflag,elec_st%nspin,pot_gspace%length, &
          symm%ntrans,symm%cell_symmetry, &
          clust%atom_num, kpnum, elec_st%nstate, ngkmax, &
          four*cutoff, cutoff,pot_gspace%fftsize(1:3), &
          elec_st%mpgrid,elec_st%mpshift,pbc%vcell, &
          one, pbc%latt_vec, pbc%adot,twopi**3/pbc%vcell,one,pbc%bvec,pbc%bdot, &
          symm%gmtrx, symm%tnp*twopi, &
          atyp, apos, ngk, weight, kpts, &
          ifmin, ifmax, energies, occupations, warn=.true.)
     sheader = 'RHO'
     open(22, file = 'RHO', status = 'replace', form = 'unformatted')
     call write_binary_header(22,sheader,iflag,elec_st%nspin,pot_gspace%length, &
          symm%ntrans,symm%cell_symmetry, &
          clust%atom_num, kpnum, elec_st%nstate, ngkmax, &
          four*cutoff,cutoff,pot_gspace%fftsize(1:3), &
          elec_st%mpgrid,elec_st%mpshift,pbc%vcell, &
          one, pbc%latt_vec, pbc%adot,twopi**3/pbc%vcell,one,pbc%bvec,pbc%bdot, &
          symm%gmtrx, symm%tnp*twopi, &
          atyp, apos, ngk, weight, kpts, &
          ifmin, ifmax, energies, occupations, warn=.true.)
! generate the g-space needed -- only need imap and igvec
     pot_gspace%imap  = .true.
     pot_gspace%igvec = .true.  
     call generate_potential_gspace(pot_gspace, pbc, symm, cutoff, 1)
   endif

! Now collect charge density and vxc

  if (parallel%iammaster) then
     allocate(vxc_m(grid%nwedge,elec_st%nspin),stat=alcstat)
     call alccheck('vxc_m',grid%nwedge*elec_st%nspin,alcstat)
     vxc_m(:,:) = zero
     allocate(rho_m(grid%nwedge),stat=alcstat)
     call alccheck('rho_m',grid%nwedge,alcstat)
     rho_m(:) = zero
  endif

  call collect_function(parallel,rho(1,1))
  if (parallel%iammaster) &
     call dcopy(grid%nwedge,parallel%ftmp,1,rho_m,1)
  do isp = 1, elec_st%nspin
     call collect_function(parallel,vxc(1:parallel%mydim,isp))
     if (parallel%iammaster) &
         call dcopy(grid%nwedge,parallel%ftmp,1,vxc_m(1,isp),1)
  enddo

! Now unfold it to the full grid

  if (parallel%iammaster) then

     allocate(vxc_tmp(pbc%maxdfft,elec_st%nspin),stat=alcstat)
     call alccheck('vxc_tmp',pbc%maxdfft*elec_st%nspin,alcstat)
     allocate(rho_tmp(pbc%maxdfft),stat=alcstat)
     call alccheck('rho_tmp',pbc%maxdfft,alcstat)

     vxc_tmp(:,:) = cmplx(zero,zero,dp)
     rho_tmp(:) = cmplx(zero,zero,dp)
     do k = pbc%mz, pbc%mz + pbc%n3 - 1
        do j = pbc%my, pbc%my + pbc%n2 - 1
           do i = pbc%mx, pbc%mx + pbc%n1 - 1
              iadd = ((k-pbc%mz)*pbc%n2 + j-pbc%my)*pbc%n1 + i-pbc%mx+1
              jj = grid%rindex(grid%indexg(i,j,k) )
              do isp = 1, elec_st%nspin
                 vxc_tmp(iadd,isp) = vxc_m(jj,isp)
              enddo
              rho_tmp(iadd) = rho_m(jj) * pbc%vcell
           enddo
        enddo
     enddo
  
     deallocate(vxc_m)
     deallocate(rho_m)

     ! Perform the FFTs.
     do isp = 1, elec_st%nspin
        call cfftw (vxc_tmp(:,isp),pbc%n1,pbc%n2,pbc%n3,1)
     enddo
     call cfftw (rho_tmp,pbc%n1,pbc%n2,pbc%n3,1)

     allocate(vxc_g(pbc%nstar,elec_st%nspin),stat=alcstat)
     call alccheck('vxc_g',pbc%nstar*elec_st%nspin,alcstat)
     allocate(rho_g(pbc%nstar),stat=alcstat)
     call alccheck('rho_g',pbc%nstar,alcstat)

     ! Symmetrization of vxc and rho
     vxc_g(:,:) = cmplx(zero,zero,dp)
     rho_g(:) = cmplx(zero,zero,dp)
     do jj=1,pbc%ng
        js = pbc%inds(jj)
        i = 1 + pbc%kgv(1,jj)
        if (i <= 0) i=pbc%n1+i
        j = 1 + pbc%kgv(2,jj)
        if (j <= 0) j=pbc%n2+j
        k = 1 + pbc%kgv(3,jj)
        if (k <= 0) k=pbc%n3+k
        call get_address(pbc,i,j,k,iadd,fi_s)

        do isp = 1, elec_st%nspin
           ztmp = vxc_tmp(iadd,isp) 
           tmpr = real(ztmp,dp)
           tmpi = pbc%conj(jj)*aimag(ztmp)
           ztmp = cmplx(tmpr,tmpi,dp)/real(pbc%mstar(js),dp)
           vxc_g(js,isp) = vxc_g(js,isp) + ztmp * pbc%phase(jj)
        enddo

        ztmp = rho_tmp(iadd) 
        tmpr = real(ztmp,dp)
        tmpi = pbc%conj(jj)*aimag(ztmp)
        ztmp = cmplx(tmpr,tmpi,dp)/real(pbc%mstar(js),dp)
        rho_g(js) = rho_g(js) + ztmp * pbc%phase(jj)
     enddo

     deallocate(vxc_tmp)
     deallocate(rho_tmp)
 
! Now finally put this in the appropriate container for writing out
     allocate(vxc_g_full(pot_gspace%length,elec_st%nspin))
     allocate(rho_g_full(pot_gspace%length))
     vxc_g_full(:,:) = cmplx(zero,zero,dp)
     rho_g_full(:) = cmplx(zero,zero,dp)
     do jj=1,pbc%ng
        js = pbc%inds(jj)
        i = pbc%kgv(1,jj)
        j = pbc%kgv(2,jj)
        k = pbc%kgv(3,jj)
        iadd = pot_gspace%indexg(i,j,k)
        if (iadd > 0) then
           do isp = 1, elec_st%nspin
              vxc_g_full(iadd,isp) = vxc_g(js,isp)*pbc%phase(jj)
           enddo
           rho_g_full(iadd) = rho_g(js)*pbc%phase(jj)
        endif
     enddo

     write(7,*) 'rho',rho_g_full(pot_gspace%ig0)
     write(7,*) 'vxc',vxc_g_full(pot_gspace%ig0,1)
   
     deallocate(rho_g)
     deallocate(vxc_g)

! Write the g-space information to file
     write(21) 1   ! How many chunks
     write(22) 1   ! How many chunks
     write(21) pot_gspace%length ! How many records
     write(22) pot_gspace%length ! How many records
     write(21) ((pot_gspace%gvec(j, i), j = 1, 3), i = 1, pot_gspace%length)
     write(22) ((pot_gspace%gvec(j, i), j = 1, 3), i = 1, pot_gspace%length)
     

! Now write the vxc and rho out!
     write(21) 1   ! How many chunks
     write(22) 1   ! How many chunks
     write(21) pot_gspace%length ! How many records
     write(22) pot_gspace%length ! How many records
     write(21) ((vxc_g_full(i,isp), i=1,pot_gspace%length),isp=1,elec_st%nspin)
     write(22) (rho_g_full(i), i = 1, pot_gspace%length)
     close(21)
     close(22)

     deallocate(vxc_g_full)
     deallocate(rho_g_full)

!Now for the wavefunction file WFN
     iflag = 2
     sheader = 'WFN'
     open(23, file = 'WFN', status = 'replace', form = 'unformatted')
     call write_binary_header(23,sheader,iflag,elec_st%nspin,pot_gspace%length, &
          symm%ntrans,symm%cell_symmetry, &
          clust%atom_num, kpnum, elec_st%nstate, ngkmax, &
          four*cutoff, cutoff,pot_gspace%fftsize(1:3), &
          elec_st%mpgrid,elec_st%mpshift,pbc%vcell, &
          one, pbc%latt_vec, pbc%adot,twopi**3/pbc%vcell,one,pbc%bvec,pbc%bdot, &
          symm%gmtrx, symm%tnp*twopi, &
          atyp, apos, ngk, weight, kpts, &
          ifmin, ifmax, energies, occupations, warn=.true.)
     write(23) 1   ! How many chunks
     write(23) pot_gspace%length ! How many records
     write(23) ((pot_gspace%gvec(j, i), j = 1, 3), i = 1, pot_gspace%length)

! At this point the pot_gspace is done and the memory can be released
     call destroy_gspace(pot_gspace)
   endif ! master

   if (elec_st%cplx) then
      allocate(zwftmp(parallel%mydim),stat=alcstat)
      call alccheck('zwftmp',parallel%mydim,alcstat)
   else
      allocate(wftmp(parallel%mydim),stat=alcstat)
      call alccheck('wftmp',parallel%mydim,alcstat)
   endif

! For all kpoints
   do irk = 1, kpnum
      if (parallel%iammaster) then
! Now allocate the gspace for each kpoint
          k_gspace(irk)%imap  = .true.
          k_gspace(irk)%igvec = .true.  
          call generate_k_gspaces(k_gspace(irk), ngkmax, irk, fftsize, &
           pbc, kpts(1:3,irk), weight(irk), symm, cutoff, 1)
          write(23) 1   ! How many chunks
          write(23) k_gspace(irk)%length ! How many records
          write(23) ((k_gspace(irk)%gvec(j, i), j = 1, 3), i = 1, k_gspace(irk)%length)
! allocate space to hold the collected wavefunctions

          if (elec_st%cplx) then
             allocate(wfn_m(1:grid%nwedge), stat=alcstat)
             call alccheck('wfn_m',grid%nwedge,alcstat)
          else
             allocate(dwfn_m(1:grid%nwedge), stat=alcstat)
             call alccheck('dwfn_m',grid%nwedge,alcstat)
          endif
          allocate(wfn_tmp(pbc%maxdfft),stat=alcstat)
          call alccheck('wfn_tmp',pbc%maxdfft,alcstat)
          allocate(wfn_g_full(k_gspace(irk)%length,elec_st%nspin),stat=alcstat)
          call alccheck('wfn_g_full',k_gspace(irk)%length*elec_st%nspin,alcstat)

          kmax(:) = k_gspace(irk)%kmax(:)
      endif

      if (irk == 1) then 
         jrep = 0
         jrepinv = 0
      endif
      do isp = 1, elec_st%nspin
        do kk = 1, elec_st%nstate
           irp = elec_st%irep(kk,irk,isp)
           jrep(isp,irp) = jrep(isp,irp) + 1
           if (irk == 1) jrepinv(isp,kk) = jrep(isp,irp)
           !if (parallel%iammaster) write(7,('(4i)')) isp, irk, kk, jrepinv(isp,kk)
        enddo
      enddo
      !if (parallel%iammaster) call myflush(7)

! For each state at each k-point
      do ii = 1, elec_st%nstate
! Collect the wavefunction and unfold it

         do isp = 1, elec_st%nspin
           irp = elec_st%irep(ii,irk,isp)
           if (elec_st%cplx) then
              if ( elec_st%eig(irp,irk,isp)%group == parallel%mygroup ) then
                  call zcopy(parallel%mydim,elec_st%eig(irp,irk,isp)% &
                      zwf(1,jrepinv(isp,ii)),1,zwftmp,1)
              else
                  zwftmp = zzero
              endif
              call collect_zfunction(parallel,zwftmp)
#ifdef MPI
              if (parallel%iamgmaster) then
                 jj = elec_st%eig(irp,irk,isp)%group
                 call MPI_BCAST(parallel%zftmp,grid%nwedge, &
                      MPI_DOUBLE_COMPLEX,jj,parallel%gmaster_comm,mpinfo)
              endif
#endif

              if (parallel%iammaster) then
                 call zcopy(grid%nwedge,parallel%zftmp,1,wfn_m(1),1)
                 wfn_tmp(:) = cmplx(zero,zero,dp)
                 do k = pbc%mz, pbc%mz + pbc%n3 - 1
                    do j = pbc%my, pbc%my + pbc%n2 - 1
                       do i = pbc%mx, pbc%mx + pbc%n1 - 1
                          iadd = ((k-pbc%mz)*pbc%n2 + j-pbc%my)*pbc%n1 + i-pbc%mx+1
                          jj = grid%rindex(grid%indexg(i,j,k) )
                          wfn_tmp(iadd) = wfn_m(jj) !*character(grid%rtrans(grid%indexg(i,j,k)))
                       enddo
                    enddo
                 enddo
              endif

            else

              if ( elec_st%eig(irp,irk,isp)%group == parallel%mygroup ) then
                 call dcopy(parallel%mydim,elec_st%eig(irp,irk,isp)% &
                      wf(1,jrepinv(isp,ii)),1,wftmp,1)
              else
                 wftmp = zzero
              endif

              call collect_function(parallel,wftmp)
#ifdef MPI
              if (parallel%iamgmaster) then
                 jj = elec_st%eig(irp,irk,isp)%group
                 call MPI_BCAST(parallel%ftmp,grid%nwedge, &
                      MPI_DOUBLE_PRECISION,jj,parallel%gmaster_comm,mpinfo)
              endif
#endif

              if (parallel%iammaster) then
                 call dcopy(grid%nwedge,parallel%ftmp,1,dwfn_m(1),1)
                 wfn_tmp(:) = cmplx(zero,zero,dp)
                 do k = pbc%mz, pbc%mz + pbc%n3 - 1
                    do j = pbc%my, pbc%my + pbc%n2 - 1
                       do i = pbc%mx, pbc%mx + pbc%n1 - 1
                          iadd = ((k-pbc%mz)*pbc%n2 + j-pbc%my)*pbc%n1 + i-pbc%mx+1
                          jj = grid%rindex(grid%indexg(i,j,k) )
                          wfn_tmp(iadd) = cmplx(dwfn_m(jj),zero,dpc) ! *character(grid%rtrans(grid%indexg(i,j,k)))
                       enddo
                    enddo
                 enddo
              endif

           endif

           if (parallel%iammaster) then

              ! calculate the vxc matrix elements (wfn in real space and vxc also in real space -- fft_convolution)

              ! Perform the FFTs.
              call cfftw (wfn_tmp,pbc%n1,pbc%n2,pbc%n3,1)

              wfn_g_full(:,isp) = cmplx(zero,zero,dpc)
! Now finally put this in the appropriate container for writing out
              factor = dsqrt(real(pbc%n1*pbc%n2*pbc%n3,dp))

              do jj=1,pbc%ng
                 g1 = pbc%kgv(1,jj)
                 g2 = pbc%kgv(2,jj)
                 g3 = pbc%kgv(3,jj)
                 if (g1 > kmax(1) .or. g1 < -kmax(1) .or. & 
                     g2 > kmax(2) .or. g2 < -kmax(2) .or. &
                     g3 > kmax(3) .or. g3 < -kmax(3)) cycle
                 i = 1 + g1
                 if (i <= 0) i=pbc%n1+i
                 j = 1 + g2
                 if (j <= 0) j=pbc%n2+j
                 k = 1 + g3
                 if (k <= 0) k=pbc%n3+k
                 call get_address(pbc,i,j,k,iadd2,fi_s)
                 iadd = k_gspace(irk)%indexg(g1,g2,g3)
                 if (iadd > 0) then
                    wfn_g_full(iadd,isp) = wfn_tmp(iadd2)*factor
                 endif
              enddo

!Normalize the wave function
             norm = dznrm2(k_gspace(irk)%length, wfn_g_full(:,isp), 1)

!Scale the function by the norm
             call zscal(k_gspace(irk)%length, ONE/norm, wfn_g_full(:,isp), 1)
 
           endif !master
         enddo  ! spin

         if (parallel%iammaster) then
            write(23) 1
            write(23) k_gspace(irk)%length*elec_st%nspin ! How many records
            write(23) ((wfn_g_full(i,isp), i = 1, k_gspace(irk)%length),isp=1,elec_st%nspin)
         endif
      enddo !nstate

! At this point the k_gspace(irk) is done and the memory can be released
      if (parallel%iammaster) then
          call destroy_gspace(k_gspace(irk))
          deallocate(wfn_m)
          deallocate(wfn_tmp)
          deallocate(wfn_g_full)
      endif
   enddo  ! k-points
   if (parallel%iammaster) then
      close(23)
   endif
     
  

  return

#else
  return

#endif

end subroutine gw_write
