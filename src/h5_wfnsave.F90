!===============================================================
!
!  Copyright (C) 2005 Finite Difference Research Group
!  This file is part of parsec, http://www.ices.utexas.edu/parsec/
!
!  This subroutine saves grid, potential, wave function, energy
!  levels, etc, in parsec.dat file for restarting purposes.
!
!---------------------------------------------------------------
subroutine h5_wfnsave(elec_st,grid,pbc,rsymm,parallel,vnew,rho,iunit,h5_file_id,&
spin3d)

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
#ifdef USEHDF5
  use hdf5
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
#ifdef USEHDF5
  integer(hid_t), intent(in) :: h5_file_id
  integer(hid_t) :: tmp_dataset_id, tmp_attr_id, tmp_aspace_id
  integer(hid_t) :: tmp_group
  integer(hid_t) :: atype_id
  integer(hsize_t), dimension(1) :: tmp_aspace_size
#else
  integer, intent(in) :: h5_file_id ! not being used, still an input argument
#endif
  integer        :: hdferr, tmp_rank
  !---------------------------------------------------------------
#ifdef USEHDF5

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
     ! create /datalabel attribute
     ! The steps to create an attribute are as follows:
     !     1. Obtain the object identifier that the attribute is to be attached
     !     to.
     !     2. Define the characteristics of the attribute and specify the
     !     attribute creation property list.
     !     - Define the datatype.
     !     - Define the dataspace.
     !     - Specify the attribute creation property
     !     list.
     !     3. Create the attribute.
     !     4. Close the attribute and
     !     datatype, dataspace, and attribute
     !     creation property list, if
     !     necessary.

     call h5gcreate_f(h5_file_id,"/details",tmp_group,hdferr)
     ! ========================================
     ! Write a single line string attribute
     ! ========================================
     ! create space size defintions and simple dataspace
     tmp_aspace_size(1) = 1
     tmp_rank = 1
     call h5screate_simple_f(tmp_rank,tmp_aspace_size, tmp_aspace_id, hdferr)
     ! create string-like variable
     call h5tcopy_f(H5T_NATIVE_CHARACTER, atype_id, hdferr)
     call h5tset_size_f(atype_id, int(len(datelabel),HSIZE_T), hdferr)
     ! create attribute attached to base group
     call h5acreate_f(tmp_group,"date", H5T_STRING, tmp_aspace_id,&
     tmp_attr_id, hdferr)
     ! write  /datalabel attribute
     call h5awrite_f(tmp_attr_id, atype_id, datelabel, tmp_aspace_size, hdferr)
     ! close  /datalabel attribute
     call h5aclose_f(tmp_attr_id, hdferr)


     !============
     ! write stype
     !============
     call h5acreate_f(tmp_group,"stype", H5T_NATIVE_INTEGER, tmp_aspace_id,&
     tmp_attr_id, hdferr)
     tmp_aspace_size = 1
     call h5awrite_f(tmp_attr_id, H5T_NATIVE_INTEGER, stype, tmp_aspace_size, hdferr)
     call h5aclose_f(tmp_attr_id, hdferr)

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
     call collect_function(parallel,vnew(1,isp))
     if (parallel%iammaster) write(iunit) &
          (parallel%ftmp(kk),kk=1,grid%nwedge)
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
     call collect_function(parallel,vnew(1,2))
     if (parallel%iammaster) write(iunit) &
          (parallel%ftmp(kk),kk=1,grid%nwedge)
     call collect_function(parallel,rho(1,3))
     if (parallel%iammaster) write(iunit) &
          (parallel%ftmp(kk),kk=1,grid%nwedge)
    if (stype(2)>2) then
       do tt = 1,3
          call collect_function(parallel,spin3d(1,tt))
          if (parallel%iammaster) write(iunit) &
          (parallel%ftmp(kk),kk=1,grid%nwedge)
      enddo
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
#endif
end subroutine h5_wfnsave
