!===============================================================
!
! Copyright (C) 2005 Finite Difference Research Group
! This file is part of parsec, http://www.ices.utexas.edu/parsec/
!
! This subroutine terminates the procs and closes the parsec.out
! file if exiting upon an error.
!
!---------------------------------------------------------------
subroutine exit_err (clust,elec_st,grid,pot,p_pot,nloc_p_pot,u_pot, &
     move, mol_dynamic,pbc,mixer,solver,symm,rsymm,parallel,&
     file_id,ierr_in)

  use constants
  use matvec_module
  use cluster_module
  use electronic_struct_module
  use grid_module
  use potential_module
  use pseudo_potential_module
  use non_local_psp_module
  use movement_module
  use molecular_dynamic_module
  use pbc_module
  use mixer_module
  use eigen_solver_module
  use symmetry_module
  use parallel_data_module
#ifdef USEHDF5
  use hdf5
#endif
#ifdef MPI
  use mpi
#endif
  implicit none
#ifdef ITAC
  ! inlcude trace-analyzer definitions
  include 'VT.inc'
#endif 
  !
  ! Input/Output variables:
  !
  ! Structures to destroy
  ! the cluster
  type(cluster), intent(inout) :: clust
  ! electronic structure
  type (electronic_struct), intent(inout) :: elec_st
  ! grid related data
  type (grid_data), intent(inout) :: grid
  ! potential related data
  type (potential), intent(inout) :: pot
  ! pseudopotential related data
  type (pseudo_potential), intent(inout) :: p_pot
  ! non local pseudopotential related data
  type (nonloc_pseudo_potential), intent(inout) :: nloc_p_pot
  ! on-site Coulomb interaction related data
  type (nonloc_pseudo_potential), intent(inout) :: u_pot
  ! molecular dynamic related data
  type (molecular_dynamic), intent(inout) :: mol_dynamic
  ! movement step for lbfgs
  type (movement), intent(inout) :: move
  ! periodic boundary conditions data
  type (pbc_data), intent(inout) :: pbc
  ! mixer related data
  type (mixer_data), intent(inout) :: mixer
  ! solver related data
  type (eigen_solver), intent(inout) :: solver
  ! symmetry operations
  type (symmetry), intent(inout) :: symm,rsymm
  ! parallel computation related data
  type (parallel_data), intent(inout) :: parallel

  ! error flag
  integer, intent(in) :: ierr_in
#ifdef USEHDF5
     INTEGER(HID_T), INTENT(IN) :: file_id  ! File identifier
     INTEGER                    :: ierr  ! Error flag
#else
     INTEGER, INTENT(IN) :: file_id  ! integer to replace it
#endif
  !
  ! Work variables:
  !
#ifdef MPI
  ! communicator to be used
  integer comm
  ! variable for the return code of mpi calls
  integer mpinfo
#endif
  ! temporary variable for error flag
  integer ::  ierr_tmp
  ! file unit numbers
  integer :: iunit, units(1:22) =(/ 4, 6, 7, 8, 10, 11,&
             13, 15, 16, 20, 45,&
             46, 61, 66, 67, 76,&
             78, 79, 88, 90, 91,&
             92 /)
         !close(unit=6)  ! parsec.out
         !close(unit=7)  ! fort.7 in most cases..?
         !close(unit=45) ! relax_restart.dat
         !close(unit=46) ! bfgs.dat
         !close(unit=61) ! eigen.dat
         !close(unit=66) ! atom.cor
         !close(unit=67) ! manual.dat
         !close(unit=76) ! mdinit.dat
         !close(unit=78) ! md_nrg.dat
         !close(unit=79) ! md_mech.dat
         !close(unit=90) ! occup.in
         !close(unit=91) ! polar.dat
         !close(unit=92) ! movie.xyz

  ! whether or not a file is open
  logical :: lopen

  !---------------------------------------------------------------

  ! share error flag
#ifdef MPI
  call MPI_allREDUCE(ierr_in,ierr_tmp,1, &
       MPI_INTEGER,MPI_MAX,parallel%comm,mpinfo) !AJB: mpi max?
#else
  ierr_tmp = ierr_in
#endif

#ifdef DEBUG
  write(9,*) ' Error check ',ierr_tmp,ierr_in
  call myflush(9)
#endif

if(ierr_tmp /= 0) then
  if (ierr_tmp > 0) then
     write(9,*) ' Error code ierr = ',ierr_tmp
  elseif (ierr_tmp == -1) then
     write(9,*) ' End of Program '
  endif
     call myflush(9)
! destroy all the structures
#ifndef MDSTR
        call matvec_destroy()
        call destroy_cluster (clust)
        call destroy_electronic_struct (elec_st)
        call destroy_grid (grid)
        call destroy_potential (pot)
        call destroy_pseudo_potential (p_pot)
        call destroy_movement(move)
        call destroy_molecular_dynamic (mol_dynamic)
        call destroy_nonloc_pseudo_pot (nloc_p_pot)
        call destroy_nonloc_pseudo_pot (u_pot)
        if (pbc%is_on) call destroy_pbc (pbc)
        call destroy_mixer (mixer)
        call destroy_eigen_solver (solver)
        call destroy_symmetry (symm)
        call destroy_symmetry (rsymm)
        call destroy_parallel_data (parallel)
#endif
     ! close all files
     ! close parsec.out
     if (parallel%iammaster) then
            do iunit=1,size(units)
                 inquire(unit=units(iunit),opened=lopen)
                 if (lopen) close(unit=units(iunit))
            enddo
     endif
#ifdef USEFFTW3
     ! release all allocated memory from fftw
     call dfftw_cleanup()
#endif
#ifdef USEHDF5 
     ! Terminate access to the file.
     CALL h5fclose_f(file_id, ierr)
!    Close FORTRAN interface.
     CALL h5close_f(ierr)
#endif

#ifdef ITAC
  ! Finalize tracer
  ! call VTFINI(mpinfo)
#endif 
#ifdef MPI
     !call MPI_BARRIER(parallel%comm,mpinfo) !bugfix UNNEEDED
     call MPI_FINALIZE(mpinfo)
#endif
     if (ierr_tmp /= 0) stop
  endif

end subroutine exit_err
!===============================================================
!
! Checks allocation of array with size size, name wst. If
! allocation fails, an error message is writen in out.* file and
! the calculation is aborted in an abrupt and possibly nasty way.
!
! All "big arrays" (i.e., with size proportional to grid size or
! square of number of eigenvalues) should have allocation checked.
!
!---------------------------------------------------------------
subroutine alccheck(wst,isize,istat)
#ifdef MPI
  use mpi
#endif
  implicit none
  !
  ! Input/Output variables:
  !
  integer, intent(in) :: isize, istat
  character (len=*), intent(in) :: wst
  !
  ! Work variables:
  !
  logical :: ltest          ! true if unit 7 is open
  integer :: aborterr       ! just for completness of mpi call
  if(istat /= 0) then
     write(9,'(a,a)') ' *** MEMORY ALLOCATION FAILED ON ARRAY ',wst
     write(9,'(a,i12)') '     ARRAY SIZE: ',isize
     call myflush(9)
     inquire(unit=7,opened=ltest)
     if (ltest) then
        write(7,'(a,a)') ' *** MEMORY ALLOCATION FAILED ON ARRAY ',wst
        write(7,'(a,i12)') '     ARRAY SIZE: ',isize
        call myflush(7)
     endif
#ifdef MPI
     call MPI_ABORT(MPI_COMM_WORLD,-1,aborterr)
#endif
     stop
  endif

end subroutine alccheck
!===============================================================
