!===============================================================
!
!  Copyright (C) 2005 Finite Difference Research Group
!  This file is part of parsec, http://www.ices.utexas.edu/parsec/
!
!  This subroutine saves grid, potential, wave function, energy
!  levels, etc, in parsec.dat file for restarting purposes.
!
!---------------------------------------------------------------
subroutine wfnsave(elec_st,grid,pbc,rsymm,parallel,vnew,rho,iunit,spin3d)

  use constants
  use electronic_struct_module
  use grid_module
  use pbc_module
  use symmetry_module
  use parallel_data_module

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
  !  grid related data
  type (grid_data), intent(in) :: grid
  !  periodic boundary conditions data
  type (pbc_data), intent(in) :: pbc
  !  symmetry operations:
  type (symmetry), intent(in) :: rsymm
  !  parallel computation related data
  type (parallel_data), intent(inout) :: parallel
  !  distributed self-consistent potential and electron density
  !  (passed outside the structure to overcome a bug with the IBM compiler)
  real(dp), intent(in) :: vnew(parallel%mydim,elec_st%nspin)
  real(dp), intent(in) :: rho(parallel%mydim,2*elec_st%nspin-1)
  real(dp), intent(in), optional :: spin3d(parallel%mydim,3)
  !  number of output unit
  integer, intent(in) :: iunit
  !
  !  Work variables:
  !
  !  communicator
  integer comm
  !  exit code for mpi calls
  integer mpinfo
  !  counters
  integer isp, jj, kk, ii, dim
  real(dp) :: tmpd

  integer msgtype, nstate
  !  counters for representations
  integer irp
  !  temporary arrays for eigenvalues/occupations, before reshaping
  !  according to symmetry operations
  real(dp), dimension(:), allocatable :: en_tmp, occ_tmp
  !  jrep keeps track of how many eigenstates are already in each
  !  representation
  integer :: jrep(elec_st%nrep),jrepinv(elec_st%nspin,elec_st%nstate)
  !  work array
  complex(dpc) :: zwftmp(parallel%mydim*parallel%mxwd)
  real(dp) :: wftmp(parallel%mydim)
  !  date/time tag
  character (len=26) :: datelabel
  !
  !  flags for type of system
  !  stype(1) = nspin
  !  stype(2) = 0 : real wavefunctions
  !  stype(2) = 1 : complex wavefunctions
  !  stype(2) = 2 : spin-orbit (complex doubled wavefunctions)
  !  stype(2) = 3 : noncollinear (spin3d present)
  !  stype(2) = 4 : spin-orbit and noncollinear (spin3d present)
  !  stype(3) = 0 : confined system (cluster)
  !  stype(3) = 1 : system periodic on all directions (bulk)
  !  stype(3) = 2 : system periodic on x direction only (wire)
  !
  integer stype(10)

#ifdef MPI
  integer status(MPI_STATUS_SIZE)
#endif

  !  kpoint variables
  integer ikp, kpnum
  !  weight of each kpoint being examined
  real(dp), dimension(:), allocatable :: weight
  !  coordinates of each kpoint
  real(dp), dimension(:,:), allocatable :: kpts
  integer tt
  !---------------------------------------------------------------

  if (elec_st%nkpt > 0) then 
     kpnum = elec_st%nkpt
     allocate(kpts(3,kpnum))
     kpts = elec_st%kpts
     allocate(weight(kpnum))
     weight = elec_st%kpwt
  else
     kpnum = 1
     allocate(kpts(3,kpnum))
     kpts = zero
     allocate(weight(kpnum))
     weight = one
  endif

  !  Open unit.
  if (parallel%iammaster) &
       open(iunit,file='parsec.dat',form='unformatted',status='unknown')
  comm = parallel%comm
  !
  !  Build system flags.
  !
  stype = 0
  stype(1) = elec_st%nspin
  if (elec_st%cplx) then
     if(elec_st%mxwd == 2) then
         if(elec_st%is_so) then 
           if(present(spin3d)) then 
              stype(2) = 4
           else
              stype(2) = 2
          endif
        else
         if(present(spin3d))stype(2) = 3
        endif
     else
         stype(2) = 1
     endif
  endif
  select case(pbc%per)
  case(0)
     stype(3) = 0
  case(1)
     stype(3) = 2
  case(2)
     stype(3) = 3
  case(3)
     stype(3) = 1
  end select
  dim = parallel%mydim*parallel%mxwd

  if (parallel%iammaster) then
     rewind(iunit)
     call custom_date_time(datelabel,tmpd)

     write(iunit) datelabel
     write(iunit) stype
     if (pbc%is_on) then
        write(iunit) (grid%step(kk),kk=1,3),(pbc%box_size(kk),kk=1,3), &
             grid%rmax
        write(iunit) ((pbc%latt_vec(jj,kk),jj=1,3),kk=1,3)
        write(iunit) grid%lap_dir_num
        write(iunit) grid%lap_dir
        write(iunit) grid%lap_dir_step
        write(iunit) kpnum
        write(iunit) elec_st%kptmethod
        write(iunit) (elec_st%mpgrid(kk),kk=1,3)
        write(iunit) (elec_st%mpshift(kk),kk=1,3)
        write(iunit) ((kpts(jj,kk),jj=1,3),kk=1,kpnum)
        write(iunit) (weight(kk),kk=1,kpnum)
     else
        write(iunit) grid%step(1), grid%rmax
     endif
     write(iunit) grid%ndim
     write(iunit) grid%shift
     write(iunit) grid%nwedge,rsymm%ntrans
     write(iunit) rsymm%rmtrx
     write(iunit) rsymm%trans
     write(iunit) rsymm%tnp
     write(iunit) rsymm%alatt, rsymm%invlat
     write(iunit) ((rsymm%chi(kk,jj),kk=1,rsymm%ntrans),jj=1,rsymm%ntrans)
     write(iunit) (grid%kx(kk),grid%ky(kk),grid%kz(kk),kk=1,grid%nwedge)
  endif
  jrepinv = 0
  do isp = 1, elec_st%nspin/elec_st%mxwd
     jj = isp-1+elec_st%nspin
     !
     !  Print out eigenvalues for all k-points.
     !
     do ikp = 1, kpnum
        jrep = 0
        do kk = 1, elec_st%nstate
           irp = elec_st%irep(kk,ikp,isp)
           jrep(irp) = jrep(irp) + 1
           if (ikp == 1) jrepinv(isp,kk) = jrep(irp)
        enddo
        if (parallel%iammaster) then
           allocate(en_tmp(elec_st%nstate))
           allocate(occ_tmp(elec_st%nstate))
           jrep = 0
           do kk = 1, elec_st%nstate
              irp = elec_st%irep(kk,ikp,isp)
              jrep(irp) = jrep(irp) + 1
              if (ikp == 1) jrepinv(isp,kk) = jrep(irp)
              en_tmp(kk) = elec_st%eig(irp,ikp,isp)%en(jrep(irp))
              occ_tmp(kk) = elec_st%eig(irp,ikp,isp)%occ(jrep(irp))
           enddo
           write(iunit) elec_st%nstate
           write(iunit) (elec_st%irep(kk,ikp,isp),kk=1,elec_st%nstate)
           write(iunit) (en_tmp(kk),kk=1,elec_st%nstate)
           write(iunit) (occ_tmp(kk),kk=1,elec_st%nstate)
           deallocate(en_tmp, occ_tmp)
        endif
     enddo
#ifdef MPI
  call MPI_Barrier(comm,mpinfo)
#endif
#ifdef AJB_DEBUG
      write(9,*) ' going to collect potential for spin component , isp=', isp
#endif
     call collect_function(parallel,vnew(1,isp))
     if (parallel%iammaster) write(iunit) &
          (parallel%ftmp(kk),kk=1,grid%nwedge)
#ifdef AJB_DEBUG
      write(9,*) ' going to collect charge for component , jj= ', jj
#endif
     call collect_function(parallel,rho(1,jj))
     if (parallel%iammaster) write(iunit) &
          (parallel%ftmp(kk),kk=1,grid%nwedge)
  enddo
  !
  !  If spin-orbit potential is present or noncollinear states, then the 
  !  potential and density for the second spin component must be printed out.
  !  The loop above goes only over the first spin component.
  !
  if (elec_st%mxwd == 2) then
#ifdef AJB_DEBUG
      write(9,*) ' going to collect potential for spin component 2'
#endif
     call collect_function(parallel,vnew(1,2))
     if (parallel%iammaster) write(iunit) &
          (parallel%ftmp(kk),kk=1,grid%nwedge)
#ifdef AJB_DEBUG
      write(9,*) ' going to collect charge for component 3 '
#endif
     call collect_function(parallel,rho(1,3))
     if (parallel%iammaster) write(iunit) &
          (parallel%ftmp(kk),kk=1,grid%nwedge)
    if (stype(2)>2) then
#ifdef AJB_DEBUG
      write(9,*) 'skipping collection of spin3d '
#endif
       ! do tt = 1,3
       !    call collect_function(parallel,spin3d(1,tt))
       !    if (parallel%iammaster) write(iunit) &
       !    (parallel%ftmp(kk),kk=1,grid%nwedge)
      ! enddo
   endif
  endif
  !
  !  Print out wave-functions. With many k-points, must do it for all
  !  k-points. Also, must print the irreducible wedge (different
  !  k-points may have different wedges).
  !
  msgtype = 1
  do isp = 1, elec_st%nspin/elec_st%mxwd
     do ikp = 1, kpnum

        !  If elec_st%nsave < 0, then print out only occupied levels.
        if (elec_st%nsave < 0) then
           nstate = elec_st%ifmax(isp)
           allocate(elec_st%indxsave(nstate))
           do ii = 1, nstate
              elec_st%indxsave(ii) = ii
           enddo
        else
           nstate = elec_st%nsave
        endif

        if (parallel%iammaster) then
           write(iunit) grid%nwedge,rsymm%ntrans
           write(iunit) rsymm%rmtrx
           write(iunit) rsymm%trans
           write(iunit) rsymm%tnp
           write(iunit) rsymm%alatt, rsymm%invlat
           write(iunit) ((rsymm%chi(kk,jj),kk=1,rsymm%ntrans),jj=1 &
                ,rsymm%ntrans)
           write(iunit) (grid%kx(kk),grid%ky(kk),grid%kz(kk), &
                kk=1,grid%nwedge)
           write(iunit) nstate
           if (nstate > 0) write(iunit) elec_st%indxsave(1:nstate)
        endif
        do ii = 1, nstate
           kk = elec_st%indxsave(ii)
           irp = elec_st%irep(kk,ikp,isp)
           if (elec_st%cplx) then
              if ( elec_st%eig(irp,ikp,isp)%group == parallel%mygroup ) then
                 call zcopy(dim,elec_st%eig(irp,ikp,isp)% &
                      zwf(1,jrepinv(isp,kk)),1,zwftmp,1)
              else
                 zwftmp = zzero
              endif
              call collect_zfunction(parallel,zwftmp)
#ifdef MPI
              if (parallel%iamgmaster) then
                 jj = elec_st%eig(irp,ikp,isp)%group
                 call MPI_BCAST(parallel%zftmp,grid%nwedge*parallel%mxwd, &
                      MPI_DOUBLE_COMPLEX,jj,parallel%gmaster_comm,mpinfo)
              endif
#endif
              if (parallel%iammaster) write(iunit) &
                   (parallel%zftmp(jj),jj=1,grid%nwedge*parallel%mxwd)
           else
              if ( elec_st%eig(irp,ikp,isp)%group == parallel%mygroup ) then
                 call dcopy(parallel%mydim,elec_st%eig(irp,ikp,isp)% &
                      wf(1,jrepinv(isp,kk)),1,wftmp,1)
              else
                 wftmp = zzero
              endif
              call collect_function(parallel,wftmp)
#ifdef MPI
              if (parallel%iamgmaster) then
                 jj = elec_st%eig(irp,ikp,isp)%group
                 call MPI_BCAST(parallel%ftmp,grid%nwedge, &
                      MPI_DOUBLE_PRECISION,jj,parallel%gmaster_comm,mpinfo)
              endif
#endif
              if (parallel%iammaster) write(iunit) &
                   (parallel%ftmp(jj),jj=1,grid%nwedge)
           endif
        enddo

        if (elec_st%nsave < 0) deallocate(elec_st%indxsave)
     enddo
  enddo
#ifdef MPI
  call MPI_Barrier(comm,mpinfo)
#endif
  if (parallel%iammaster) close(iunit)

end subroutine wfnsave
!===============================================================
!
!  Gets the date (day-month-year) and time, with time zone. It
!  also returns the wall-clock time in absolute seconds.
!
!---------------------------------------------------------------
subroutine custom_date_time(datelabel,wallclock)

  use constants
  implicit none
  !
  !  Input/Output variables:
  !
  character (len=26), intent(out) :: datelabel
  real(dp), intent(out) :: wallclock
  !
  !  Work variables:
  !
  integer values(8)
  integer lmonth
  integer idate (8)
  character (len=11) :: bdate
  character (len=14) :: btime
  character (len=10) :: atime
  character (len=8) :: adate
  character (len=5) :: azone
  character (len=4) :: year
  character (len=2) :: day
  character (len=2) :: hour,min,sec
  character (len=3), parameter :: month(12) = (/'JAN','FEB','MAR', &
       'APR','MAY','JUN','JUL','AUG','SEP','OCT','NOV','DEC'/)
  character (len=1), parameter ::  dash = '-'
  !---------------------------------------------------------------
  call date_and_time(adate,atime,azone,idate)
  read(adate,'(a4,i2,a2)') year,lmonth,day
  write(bdate,'(a2,a1,a3,a1,a4)') day,dash,month(lmonth),dash,year
  read(atime,'(a2,a2,a2,a4)') hour,min,sec
  write(btime,'(a2,a1,a2,a1,a2,a1,a5)') hour,':',min,':',sec,' ',azone
  write(datelabel,'(a,1x,a)') bdate,btime

  call date_and_time(values=values)
  wallclock = ((values(3)*24.0d0+values(5))*60.0d0+values(6)) &
       *60.0d0+values(7)+values(8)*0.001d0

end subroutine custom_date_time
!===============================================================
