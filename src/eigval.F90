!=====================================================================
!
! Copyright (C) 2005 Finite Difference Research Group
! This file is part of parsec, http://www.ices.utexas.edu/parsec/
!
! Call one of the iterative diagonalization subroutines 
! (arpack, diagla, chebdav, chebff) for diagonalization at 
! each SCF step; or call a diagonalization subroutine only for the
! first SCF step, then switch to the subspace filtering method
! using adaptive Chebyshev polynomials (the latter approach is
! significantly faster than iterative diagonalization at each 
! SCF step).
!
! Ver1.4 main updates by ykzhou:
! (1)  implemented chebff for the 1st diagonalization 
!      (chebff is a lot faster than diagonalization, and it uses less memory)
! AJB: 2/3 need -DBETA in the make file for now
! (2)  implemented newer bounds estimator for the chebyshev filters 
! (3)  implemented the scaled chebyshev filters (in earlier versions, the filters were non-scaled)
! (4)  started removing trlanc, no longer need it, mainly due to complexity in compilation/installation
!      (this removing is not complete yet)
! Feb--April 2014, Austin
!---------------------------------------------------------------------
subroutine eigval(elec_st,pot,u_pot,solver,parallel,ipr,ierr)
  use constants
  use electronic_struct_module
  use potential_module
  use non_local_psp_module
  use eigen_solver_module
  use parallel_data_module
#ifdef MPI
  use mpi
#endif
  implicit none
#ifdef ITAC
  include 'VT.inc'
  include 'vtcommon.inc'
#endif 
  !
  ! Input/Output variables:
  !
  ! electronic structure
  type (electronic_struct), intent(inout) :: elec_st
  ! potential related data
  type (potential), intent(in) :: pot
  ! on-site Coulomb interaction related data
  type (nonloc_pseudo_potential), intent(inout) :: u_pot
  ! solver related data
  type (eigen_solver), intent(inout) :: solver
  ! parallel computation related data
  type (parallel_data), intent(in) :: parallel

  ! printout flag
  integer, intent(in) :: ipr
  ! error flag, 400 < ierr < 421
  integer, intent(out) :: ierr
  !
  ! Work variables:
  !
#ifdef MPI
  integer comm, masterid, mpinfo
#endif
  integer node, nelems,i,irp,ntotal,isc,iit,mxstate
  ! number of k-points
  integer :: kpnum
  ! number of spins
  integer :: isp,spnum
  ! integer array used to save number of truly converged eigenpairs
  ! per representation !!replaced by solver%nconvt!!
  !  integer :: nconvt(elec_st%nrep)
  integer :: kplp
  real(dp) :: dnorm, magmom
  complex(dpc) :: znorm
  ! number of eigenpairs converged
  integer,allocatable :: nec(:,:,:)
  ! maximum number of matrix-vector multiplications performed.
  integer,allocatable :: maxmvs(:,:,:)
  !
  ! Number of states to be calculated for each representation so
  ! that they add up to at least elec_st%nstate, required on input.
  !
  integer,allocatable :: nrep(:,:,:)
  !
  ! Number of additional states to be calculated for each
  ! representation. This is used as safety window, to prevent
  ! calculating too few eigenstates on some representation.
  integer,allocatable :: nrep_add(:,:,:)
  ! temporary representation array
  integer,  allocatable :: irep_in(:)
  ! output error
  integer,allocatable :: out_info(:,:,:)
  ! temporary array for communication
  integer, allocatable :: itmp(:,:,:),mxstt(:,:,:)
  ! output string
  ! maximum number of diagonalization iterations before quitting
  integer, parameter :: maxit = 5
  ! flag for chebdav_diag
  character(len=4)   :: flag

  Interface
     subroutine chebff_diag(kplp,max_spdim, nwant, nconv, max_mv, eval, &
          v_basis, conv_tol, verbose, flag, nconvt,  info, blksize)
       use parsec_global_data
       implicit none
       integer  :: max_spdim, nwant, nconv, max_mv
       real(dp) :: eval(max_spdim)
       real(dp) :: v_basis(parallel%ldn, max_spdim)
       real(dp) :: conv_tol
       integer  :: verbose, info, nconvt
       integer  :: kplp
       character(len=4)   :: flag
       integer, optional  :: blksize
     end subroutine chebff_diag
     subroutine zchebff_diag(kplp, max_spdim, nwant, nconv, max_mv, eval, &
          v_basis, conv_tol, verbose, flag, nconvt, info, blksize)
       use parsec_global_data
       implicit none
       integer  :: max_spdim, nwant, nconv, max_mv
       real(dp) :: eval(max_spdim)
       complex(dpc) :: v_basis(parallel%ldn*parallel%mxwd,max_spdim)
       real(dp) :: conv_tol
       integer  :: verbose, info, nconvt
       integer  :: kplp
       character(len=4)   :: flag
       integer, optional  :: blksize
     end subroutine zchebff_diag

     subroutine chebdav_diag(kplp,max_spdim, nwant, nconv, max_mv, eval, &
          v_basis, conv_tol, verbose, flag, nconvt, info, blksize)
       use parsec_global_data
       implicit none
       integer  :: max_spdim, nwant, nconv, max_mv
       real(dp) :: eval(max_spdim)
       real(dp) :: v_basis(parallel%ldn, max_spdim)
       real(dp) :: conv_tol
       integer  :: verbose, info, nconvt
       integer  :: kplp
       character(len=4)   :: flag
       integer, optional  :: blksize
     end subroutine chebdav_diag
     subroutine zchebdav_diag(kplp, max_spdim, nwant, nconv, max_mv, eval, &
          v_basis, conv_tol, verbose, flag, nconvt,  info, blksize)
       use parsec_global_data
       implicit none
       integer  :: max_spdim, nwant, nconv, max_mv
       real(dp) :: eval(max_spdim)
       complex(dpc) :: v_basis(parallel%ldn*parallel%mxwd,max_spdim)
       real(dp) :: conv_tol
       integer  :: verbose, info, nconvt
       integer  :: kplp
       character(len=4)   :: flag
       integer, optional  :: blksize
     end subroutine zchebdav_diag
  end Interface

  !---------------------------------------------------------------

#ifdef MPI
  masterid = parallel%masterid
  comm = parallel%group_comm
#endif
  !
  ! Set nec = elec_st%eig%nec if this eigensolver can reuse old info
  ! having eigenvalues and eigenvectors only.
  !
  kpnum = max(elec_st%nkpt,1)
  spnum = elec_st%nspin/elec_st%mxwd
  allocate(nec(elec_st%nrep,kpnum,spnum), stat=ierr)
  call alccheck('nec',elec_st%nrep*kpnum*spnum,ierr)
  nec(:,:,:) = 0
  allocate(maxmvs(elec_st%nrep,kpnum,spnum), stat=ierr)
  call alccheck('maxvms',elec_st%nrep*kpnum*spnum,ierr)
  maxmvs(:,:,:) = 0
  allocate(nrep(elec_st%nrep,kpnum,spnum), stat=ierr)
  call alccheck('nrep',elec_st%nrep*kpnum*spnum,ierr)
  nrep(:,:,:) = 0
  allocate(nrep_add(elec_st%nrep,kpnum,spnum), stat=ierr)
  call alccheck('nrep_add',elec_st%nrep*kpnum*spnum,ierr)
  nrep_add(:,:,:) = 0
  allocate(out_info(elec_st%nrep,kpnum,spnum), stat=ierr)
  call alccheck('out_info',elec_st%nrep*kpnum*spnum,ierr)
  out_info(:,:,:) = 0
  allocate(itmp(elec_st%nrep,kpnum,spnum), stat=ierr)
  call alccheck('itmp',elec_st%nrep*kpnum*spnum,ierr)
  itmp(:,:,:) = 0
  allocate(mxstt(elec_st%nrep,kpnum,spnum), stat=ierr)
  call alccheck('mxstt',elec_st%nrep*kpnum*spnum,ierr)
  mxstt(:,:,:) = 0
  do isp = 1, spnum
     do kplp = 1, kpnum
        do irp = 1, elec_st%nrep
           if (solver%name == CHEBFF) then
              nec(irp,kplp,isp) = max(elec_st%eig(irp,kplp,isp)%nec,0)
           else
              nec(irp,kplp,isp) = 0
           endif
           if (solver%name == CHEBDAV) then
              nec(irp,kplp,isp) = max(elec_st%eig(irp,kplp,isp)%nec,0)
           else
              nec(irp,kplp,isp) = 0
           endif
        enddo
     enddo
  enddo
  out_info = 0
  mxstt = 0

  ! Define variables; oversize solver%resn since nn and mxstate do
  ! not have fixed values.
  if (solver%name == DIAGLA) then
     nelems = elec_st%nstate + (2/elec_st%mxwd-1)*solver%nadd + &
          solver%winsize
     allocate(solver%resn(nelems,elec_st%nrep,kpnum,spnum),stat=ierr)
     call alccheck('solver%resn',nelems*elec_st%nrep*kpnum*spnum,ierr)
     solver%resn(:,:,:,:) = zero
  endif

  !
  ! Diagonalization loop.
  !
  elec_st%eig(:,:,:)%nec = 0
  do iit = 1, maxit
     maxmvs = 0
     !
     ! Start loop over spin
     !
     do isp = 1, spnum
        do kplp = 1, kpnum
           u_pot%isp = isp
!#ifdef AJB_DEBUG
  if (parallel%iammaster) then
     if (kpnum>1) then
         write(7,*) ' ...Working on k-point ',kplp
     endif
  endif
!#endif
           ! Retrieve the current Hartree potential
           if (parallel%mxwd == 1) then

              if (elec_st%nkpt == 0) then

                 call dcopy(parallel%mydim,pot%vnew(1,isp),1,solver%adiag,1)
                 solver%adiag(:) = solver%adiag(:) + sum(solver%coe2(0,:))

              else

                 call dcopy(parallel%mydim,pot%vnew(1,isp),1,solver%adiag,1)

                 solver%adiag(:) = solver%adiag(:) + sum(solver%coe2(0,:)) + &
                      elec_st%kpts(1,kplp)*elec_st%kpts(1,kplp) + &
                      elec_st%kpts(2,kplp)*elec_st%kpts(2,kplp) + &
                      elec_st%kpts(3,kplp)*elec_st%kpts(3,kplp)

              endif

           endif
           !
           ! Renormalize calculated eigenvectors to 1 (outside the eigensolver,
           ! eigenvectors are assumed to be normalized over the *full* grid,
           ! dot_product(wf,wf) = elec_st%nrep).
           !
           do irp = 1, elec_st%nrep
              ! if no eigenstates are needed, skip diagonalization
              if (elec_st%eig(irp,kplp,isp)%nec /= 0) cycle
              !
              ! Only processors belonging to the group in charge of this
              ! representation should work.
              !
              if ( elec_st%eig(irp,kplp,isp)%group /= &
                   parallel%mygroup ) cycle
              if (elec_st%cplx) then
                 do i = 1, elec_st%eig(irp,kplp,isp)%nec
                    znorm = DOT_PRODUCT( &
                         elec_st%eig(irp,kplp,isp)%zwf &
                         (1:parallel%mydim*parallel%mxwd,i), &
                         elec_st%eig(irp,kplp,isp)%zwf &
                         (1:parallel%mydim*parallel%mxwd,i))
                    call zpsum(znorm,1,parallel%procs_num,parallel%comm)
                    znorm = zone/sqrt(znorm)
                    call zscal(parallel%mydim*parallel%mxwd,znorm, &
                         elec_st%eig(irp,kplp,isp)%zwf(1,i),1)
                 enddo
              else
                 do i = 1, elec_st%eig(irp,kplp,isp)%nec
                    dnorm = DOT_PRODUCT( &
                         elec_st%eig(irp,kplp,isp)%wf(1:parallel%mydim,i), &
                         elec_st%eig(irp,kplp,isp)%wf(1:parallel%mydim,i))
                    call psum(dnorm,1,parallel%procs_num,parallel%comm)
                    dnorm = one/sqrt(dnorm)
                    call dscal(parallel%mydim,dnorm, &
                         elec_st%eig(irp,kplp,isp)%wf(1,i),1)
                 enddo
              endif
              ! initialize characters and matrix-vector counter
              maxmvs(irp,kplp,isp) = solver%maxmv
              solver%chi(:) = elec_st%chi(irp,:)
              ! if mxstate changes, reset eig_init flag
              ! ntotal stores the new value of mxstate
              ntotal = elec_st%eig(irp,kplp,isp)%nn
              mxstate = elec_st%eig(irp,kplp,isp)%mm
              if (mxstate /= ntotal) then
                 if (solver%do_subsp .and. mxstate < ntotal) &
                      solver%eig_init(irp,kplp,isp) = .false.
                 mxstate = ntotal
              endif
              call eig_adj_size(elec_st%eig(irp,kplp,isp),parallel%ldn* &
                   parallel%mxwd,mxstate+solver%winsize,elec_st%cplx)
              !
              ! Do Chebyshev filtering only if the number of eigenvalues has
              ! not changed.
              !
              if (solver%do_subsp .and.  &
                   solver%eig_init(irp,kplp,isp)) then
                 if (ipr >= 2) then
                    write(9,*) ' subspace (irp,kplp,isp,mm,nn,nec) = ', &
                         irp,kplp,isp,mxstate, &
                         elec_st%eig(irp,kplp,isp)%nn,nec(irp,kplp,isp)
                    call myflush(9)
                 endif
                 if (elec_st%cplx) then
                    call zsubspace(kplp,elec_st%eig(irp,kplp,isp)%nn, &
                         nec(irp,kplp,isp),maxmvs(irp,kplp,isp), &
                         elec_st%eig(irp,kplp,isp)%en, &
                         elec_st%eig(irp,kplp,isp)%zwf, &
                         irp,isp,ipr,solver%mv_blksize,out_info(irp,kplp,isp))
                 else
                    call subspace(kplp,elec_st%eig(irp,kplp,isp)%nn, &
                         nec(irp,kplp,isp),maxmvs(irp,kplp,isp), &
                         elec_st%eig(irp,kplp,isp)%en, &
                         elec_st%eig(irp,kplp,isp)%wf, &
                         irp,isp,ipr,solver%mv_blksize,out_info(irp,kplp,isp))
                 endif
                 ! In case the "chebdav" eigensolver is called after a call
                 ! to "subspace" (because the value of "nec" has changed,
                 ! for example), then it should know that all eigenpairs up
                 ! to the old value of "nec" are converged.
                 ! (--ykz: why call chebdav after a call to "subspace"?? shouldn't do this)
                 solver%nconvt(irp,kplp,isp) = nec(irp,kplp,isp)
            else
               if (ipr >= 2) then
                  write(9,*) ' eigval (irp,kplp,isp,mm,nn,nec) = ', &
                       irp,kplp,isp,mxstate, &
                       elec_st%eig(irp,kplp,isp)%nn,nec(irp,kplp,isp)
                  call myflush(9)
               endif
               do isc = 1, maxit
                  if (elec_st%cplx) then
                     select case(solver%name)
                     case (DIAGLA)
                        call zdiagla_diag (kplp,parallel, &
                             elec_st%eig(irp,kplp,isp)%nn,nec(irp,kplp,isp), &
                             solver%kss0,maxmvs(irp,kplp,isp), &
                             elec_st%eig(irp,kplp,isp)%en, &
                             elec_st%eig(irp,kplp,isp)%zwf, &
                             solver%resn(1,irp,kplp,isp),solver%toler, &
                             solver%winsize,ipr,out_info(irp,kplp,isp))
                     case (ARPACK)
! Residual vector is not reused for the moment.
!                    out_info(irp,kplp,isp) = solver%info(irp,kplp,isp)
                        out_info(irp,kplp,isp) = 0
#ifdef ITAC
                         call VTBEGIN(vt_arpack,ierr)
#endif
                        call zarpk_diag (kplp,parallel, &
                             elec_st%eig(irp,kplp,isp)%nn,nec(irp,kplp,isp), &
                             maxmvs(irp,kplp,isp),elec_st%eig(irp,kplp,isp)%en, &
                             elec_st%eig(irp,kplp,isp)%zwf,solver%toler, &
                             solver%zres_res(1,irp,kplp,isp), &
                             ipr,out_info(irp,kplp,isp))
                        solver%info(irp,kplp,isp) = out_info(irp,kplp,isp)
#ifdef ITAC
                         call VTEND(vt_arpack,ierr)
#endif
                     case (CHEBFF)
#ifdef ITAC
                         call VTBEGIN(vt_chebff,ierr)
#endif
                        if (nec(irp,kplp,isp) <= 0) then
                           flag = "new";
                        else
                           flag = "cont"
                        end if
                        call zchebff_diag(kplp,mxstate+solver%winsize, &
                             elec_st%eig(irp,kplp,isp)%nn, nec(irp,kplp,isp), &
                             maxmvs(irp,kplp,isp),elec_st%eig(irp,kplp,isp)%en, &
                             elec_st%eig(irp,kplp,isp)%zwf,solver%toler, &
                             ipr, flag, solver%nconvt(irp,kplp,isp), &
                             out_info(irp,kplp,isp),solver%mv_blksize)
#ifdef ITAC
                         call VTEND(vt_chebff,ierr)
#endif
                     case (CHEBDAV)
#ifdef ITAC
                         call VTBEGIN(vt_chebdav,ierr)
#endif
!  chebdav cannot restart correctly?
                     if (nec(irp,kplp,isp) <= 0) then
                        flag = "new";
                        ! else
                        !    flag = "cont"
                        !    !nec(irp,kplp,isp) = min(nec(irp,kplp,isp), solver%nconvt(irp,kplp,isp))
                        end if
                        if (ipr >= 2) then
                           if (nec(irp,kplp,isp) <= 0) then
                              write(9,'(a, a, a, i2, a, i5)') "Chebdav flag=" &
                                   , flag,", nec(",irp,")=", nec(irp,kplp,isp)
                           else
                              write(9,'(a, a, 2(a, i2, a, i5))') "Chebdav flag=", &
                                   flag,", nec(",irp,")=", nec(irp,kplp,isp), &
                                   ", nconvt(",irp,")=",solver%nconvt(irp,kplp,isp)
                           end if
                        endif
                        call zchebdav_diag(kplp,mxstate+solver%winsize, &
                             elec_st%eig(irp,kplp,isp)%nn, nec(irp,kplp,isp), &
                             maxmvs(irp,kplp,isp),elec_st%eig(irp,kplp,isp)%en, &
                             elec_st%eig(irp,kplp,isp)%zwf,solver%toler, &
                             ipr, flag, solver%nconvt(irp,kplp,isp), &
                             out_info(irp,kplp,isp),solver%mv_blksize)
#ifdef ITAC
                         call VTEND(vt_chebdav,ierr)
#endif
                     end select
                  else
                     select case(solver%name)
                     case (DIAGLA)
                        call diagla_diag (kplp,parallel, &
                             elec_st%eig(irp,kplp,isp)%nn,nec(irp,kplp,isp), &
                             solver%kss0,maxmvs(irp,kplp,isp), &
                             elec_st%eig(irp,kplp,isp)%en, &
                             elec_st%eig(irp,kplp,isp)%wf, &
                             solver%resn(1,irp,kplp,isp),solver%toler, &
                             solver%winsize,ipr,out_info(irp,kplp,isp))
                     case (ARPACK)
!    Residual vector is not reused for the moment
!                    out_info(irp,kplp,isp) = solver%info(irp,kplp,isp)
                        out_info(irp,kplp,isp) = 0
#ifdef ITAC
                         call VTBEGIN(vt_arpack,ierr)
#endif
                        call arpk_diag (kplp,parallel, &
                             elec_st%eig(irp,kplp,isp)%nn,nec(irp,kplp,isp), &
                             maxmvs(irp,kplp,isp),elec_st%eig(irp,kplp,isp)%en, &
                             elec_st%eig(irp,kplp,isp)%wf,solver%toler, &
                             solver%res_res(1,irp,kplp,isp), &
                             ipr,out_info(irp,kplp,isp))
#ifdef ITAC
                         call VTEND(vt_arpack,ierr)
#endif
                     !! case (TRLANC)
                     !!    call trlanc_diag (elec_st%eig(irp,kplp,isp)%nn, &
                     !!         nec(irp,kplp,isp),solver%kss0,maxmvs(irp,kplp,isp), &
                     !!         elec_st%eig(irp,kplp,isp)%en, &
                     !!         elec_st%eig(irp,kplp,isp)%wf,solver%toler, &
                     !!         out_info(irp,kplp,isp))
                     case (CHEBFF)
                        if (nec(irp,kplp,isp) <= 0) then
                           flag = "new";
                        else
                           flag = "cont"
                        end if
#ifdef ITAC
                         call VTBEGIN(vt_chebff,ierr)
#endif
                        call chebff_diag(kplp,mxstate+solver%winsize,  &
                             elec_st%eig(irp,kplp,isp)%nn, nec(irp,kplp,isp), &
                             maxmvs(irp,kplp,isp),elec_st%eig(irp,kplp,isp)%en, &
                             elec_st%eig(irp,kplp,isp)%wf,solver%toler, &
                             ipr, flag, solver%nconvt(irp,kplp,isp),  &
                             out_info(irp,kplp,isp),solver%mv_blksize)
#ifdef ITAC
                         call VTEND(vt_chebff,ierr)
#endif
                     case (CHEBDAV)
#ifdef ITAC
                         call VTBEGIN(vt_chebdav,ierr)
#endif
!  chebdav cannot restart correctly?
                    if (nec(irp,kplp,isp) <= 0) then
                        flag = "new";
                        ! else
                        !    flag = "cont"
                        !    !nec(irp,kplp,isp) = min(nec(irp,kplp,isp), solver%nconvt(irp,kplp,isp))
                        end if
                        if (ipr >= 2) then
                           if (nec(irp,kplp,isp) <= 0) then
                              write(9, '(a, a, a, i2, a, i5)') "Chebdav flag=", &
                                   flag,", nec(",irp,")=", nec(irp,kplp,isp)
                           else
                              write(9,'(a, a, 2(a, i2, a, i5))') "Chebdav flag=", &
                                   flag,", nec(",irp,")=", nec(irp,kplp,isp), &
                                   ", nconvt(",irp,")=",  &
                                   solver%nconvt(irp,kplp,isp)
                           end if
                        endif
                        call chebdav_diag(kplp,mxstate+solver%winsize,  &
                             elec_st%eig(irp,kplp,isp)%nn, nec(irp,kplp,isp), &
                             maxmvs(irp,kplp,isp),elec_st%eig(irp,kplp,isp)%en, &
                             elec_st%eig(irp,kplp,isp)%wf,solver%toler, &
                             ipr, flag, solver%nconvt(irp,kplp,isp),  &
                             out_info(irp,kplp,isp),solver%mv_blksize)
#ifdef ITAC
                         call VTEND(vt_chebdav,ierr)
#endif
                     end select
                  endif
                  
#ifdef MPI
                  call MPI_BARRIER(comm,mpinfo)
#endif
                  write(9,*) ' Eigensolver, end of iteration = ', isc
                  write(9,*) ' repr. = ',irp, ' k-point = ',kplp,' spin = ',isp 
                  write(9,*) ' maxmvs = ',maxmvs(irp,kplp,isp), &
                       ' nec = ', nec(irp,kplp,isp)
                  write(9,*) ' nstate = ', elec_st%eig(irp,kplp,isp)%nn, &
                       ' info = ',out_info(irp,kplp,isp)
                  if (out_info(irp,kplp,isp) /= 0) exit
                  if (nec(irp,kplp,isp) >= elec_st%eig(irp,kplp,isp)%nn) then
                     nec(irp,kplp,isp) = elec_st%eig(irp,kplp,isp)%nn
                     exit
                  endif
               enddo
               !
               ! If using Chebyshev-filtered subspace method, save information 
               ! for the subsequent iterations.
               !
               if (solver%do_subsp) then
                  solver%eval_loc(1:mxstate,irp,isp) = &
                       elec_st%eig(irp,kplp,isp)%en(1:mxstate)
                  solver%eig_init(irp,kplp,isp) = .true.
               endif
            endif
            mxstt(irp,kplp,isp) = mxstate
         enddo
      enddo
   enddo

#ifdef ITAC
    call VTBEGIN( vt_eigencomm_state,ierr)
#endif
   ! Share information across groups
#ifdef MPI
   itmp = out_info
   call MPI_ALLREDUCE(itmp,out_info,elec_st%nrep* &
        kpnum*spnum,MPI_INTEGER,MPI_MAX,parallel%comm,mpinfo)
   itmp = nec
   call MPI_ALLREDUCE(itmp,nec,elec_st%nrep* &
        kpnum*spnum,MPI_INTEGER,MPI_MAX,parallel%comm,mpinfo)
   itmp = maxmvs
   call MPI_ALLREDUCE(itmp,maxmvs,elec_st%nrep* &
        kpnum*spnum,MPI_INTEGER,MPI_MAX,parallel%comm,mpinfo)
   itmp = mxstt
   call MPI_ALLREDUCE(itmp,mxstt,elec_st%nrep* &
        kpnum*spnum,MPI_INTEGER,MPI_MAX,parallel%comm,mpinfo)
#endif

   ! report number of eigenpairs and matrix-vector operations
   do isp = 1, spnum
      do kplp = 1, kpnum
         if (parallel%iammaster .and. elec_st%nspin == 2) then
            if (isp == 1) write(7,*) ' SPIN UP'
            if (isp == 2) write(7,*) ' SPIN DOWN'
         endif
         do irp = 1, elec_st%nrep
            elec_st%eig(irp,kplp,isp)%mm = mxstt(irp,kplp,isp)
            if ( elec_st%eig(irp,kplp,isp)%group /= &
                 parallel%mygroup ) then
               call eig_e_adj_size(elec_st%eig(irp,kplp,isp), &
                    nec(irp,kplp,isp))
            endif
#ifdef MPI
            node = parallel%gmap(1,elec_st%eig(irp,kplp,isp)%group + 1)
            call MPI_BCAST(elec_st%eig(irp,kplp,isp)%en, &
                 nec(irp,kplp,isp),MPI_DOUBLE_PRECISION, &
                 node,parallel%comm,mpinfo)
#endif
            if (out_info(irp,kplp,isp) /= 0) then
               if (parallel%iammaster) then
                  write(7,*) ' ERROR in eigval.F while working at'
                  write(7,*) ' representation ',irp,' spin ',isp
                  write(7,*) ' error flag = ',out_info(irp,kplp,isp)
               endif
               ierr = 401
               return
            endif
            if ((nec(irp,kplp,isp) < elec_st%eig(irp,kplp,isp)%nn .or. &
                 maxmvs(irp,kplp,isp) >= solver%maxmv) .and. &
                 solver%name /= CHEBDAV) then
               if (parallel%iammaster) then
                  write(7,*)
                  write(7,*) 'ERROR: ', &
                       'Iterative diagonalization has failed!'
                  write(7,*) ' representation ',irp,' spin ',isp
                  write(7,*) isc,' iterations performed'
                  write(7,*) ' Increase value of Maximum_matvec ', &
                       'in parsec.in'
                  write(7,*) ' error flag = ',out_info(irp,kplp,isp)
               endif
               ierr = 402
               return
            endif
         enddo

         if (parallel%iammaster) then
            solver%totalmv = solver%totalmv + sum(maxmvs(:,kplp,isp))
            if(kpnum > 1) then
                write(7,*) 'Kpoint (spin)', kplp,'(',isp,')'
            endif
            write(7,*) 'Number of converged eigen-pairs per representation:', nec(:,kplp,isp)
            write(7,*) 'Number of matrix-vector multiplications: ',  maxmvs(:,kplp,isp)
            write(7,*) 'Cumulated number of matrix-vector mult. = ', solver%totalmv
            !also appearing in old format:
            !write(7,'(a,i7)') ' nec  total = ',sum(nec(:,kplp,isp))
            !write(7,'(a,i7)') ' maxmvs total = ',sum(maxmvs(:,kplp,isp))
            call myflush(7)
         endif
         !
         ! Error messages from eigensolver.
         !
         do irp = 1, elec_st%nrep
            if (out_info(irp,kplp,isp) /= 0) then
               if (parallel%iammaster) then
                  write(7,*) ' ERROR in eigval.F while working at'
                  write(7,*) ' representation ',irp,' spin ',isp
                  write(7,*) ' error flag = ',out_info(irp,kplp,isp)
               endif
               ierr = 403
               return
            endif
         enddo
      enddo
   enddo
   !
   ! Sort eigenvalues in ascending order and define the number of
   ! states per representation, nrep. Update elec_st%irep and all
   ! elec_st%eig%nec accordingly.
   !
#ifdef ITAC
    call VTEND( vt_eigencomm_state, ierr )
    call VTBEGIN( vt_eigensort_state, ierr )
#endif
   do isp = 1, spnum
      do kplp = 1, kpnum
         ntotal = 0
         do irp = 1, elec_st%nrep
            ntotal = ntotal + elec_st%eig(irp,kplp,isp)%nn
         enddo
         allocate(irep_in(ntotal))
         call eigen_sort(elec_st,isp,kplp,ntotal,solver%nadd, &
              nrep(1,kplp,isp),nrep_add(1,kplp,isp),irep_in)
         ntotal = min(sum(nrep(:,kplp,isp) +  &
              nrep_add(:,kplp,isp)),elec_st%nstate)
         elec_st%irep(:,kplp,isp) = 0
         elec_st%irep(1:ntotal,kplp,isp) = irep_in(1:ntotal)
         do irp = 1, elec_st%nrep
            elec_st%eig(irp,kplp,isp)%nn = nrep(irp,kplp,isp) &
                 + nrep_add(irp,kplp,isp)
            elec_st%eig(irp,kplp,isp)%nec = nrep(irp,kplp,isp)
         enddo
         ! add the lowest eigenvalue in set nrep_add as well
         irp = irep_in(elec_st%nstate)
         elec_st%eig(irp,kplp,isp)%nec = nrep(irp,kplp,isp) + 1
         deallocate(irep_in)

         elec_st%ntotal(isp) = ntotal
         !
         ! If no additional states are calculated, then the set of
         ! eigenstates for this representation may be incomplete. Increase
         ! the number of states to be computed (nn) and set nec to zero,
         ! so that a new diagonalization will be performed.
         !
         do irp = 1, elec_st%nrep
            if (nrep_add(irp,kplp,isp) == 0) then
               elec_st%eig(irp,kplp,isp)%nn = &
                    elec_st%eig(irp,kplp,isp)%nn + &
                    max(2*solver%nadd,elec_st%eig(irp,kplp,isp)%nn/4)
               elec_st%eig(irp,kplp,isp)%nec = 0
            endif
         enddo
         if (solver%fix_neig) then
            do irp = 1, elec_st%nrep
               elec_st%eig(irp,kplp,isp)%nn = elec_st%eig(irp,kplp,isp)%mm
               elec_st%eig(irp,kplp,isp)%nec = elec_st%eig(irp,kplp,isp)%mm
            enddo
            exit
         endif
      enddo
   enddo
#ifdef MPI
   itmp = nrep_add
   call MPI_ALLREDUCE(itmp,nrep_add,elec_st%nrep*kpnum*spnum, &
        MPI_INTEGER,MPI_MAX,parallel%comm,mpinfo)
#endif
   if ( all(nrep_add /= 0) ) exit
enddo
#ifdef ITAC
    call VTEND( vt_eigensort_state,ierr)
#endif
!
! Print error message if there are not enough eigenstates.
!
do isp = 1, spnum
   do kplp = 1, kpnum
      ntotal = 0
      do irp = 1, elec_st%nrep
         ntotal = ntotal + elec_st%eig(irp,kplp,isp)%nec
      enddo
      if (ntotal < elec_st%nstate) then
         if (parallel%iammaster) then
            write(7,*) 'ERROR: not enough eigenstates are computed! '
            write(7,*) (elec_st%eig(irp,kplp,isp)%nec,irp=1,elec_st%nrep)
            write(7,*) elec_st%nstate,' eigenstates requested but', &
                 ' only ',ntotal,' found.'
            write(7,*) 'Increase initial_subspace_acc and/or ', &
                 'subspace_buffer_size in parsec.in'
            write(7,*) ' STOP'
         endif
         ierr = 404
         return
      endif
      elec_st%ntotal(isp) = ntotal
   enddo
enddo
if (solver%name == DIAGLA) deallocate(solver%resn)

! If spin-orbit, calculate the magnetic moment of each state.
! (Must assume that wave-functions are complex!)
#ifdef ITAC
    call VTBEGIN( vt_so,ierr)
#endif
do kplp = 1, kpnum
   if (parallel%mxwd == 2) then
      do i = 1, elec_st%nstate
         magmom = zero
         irp = elec_st%irep(i,kplp,1)
         if (elec_st%eig(irp,kplp,1)%group == parallel%mygroup) then
            do iit = 1, parallel%mydim
               if(elec_st%ncl) then
                  magmom = magmom + &
                    (conjg(elec_st%eig(irp,kplp,1)%zwf(iit,i))* &
                    elec_st%eig(irp,kplp,1)%zwf(iit+parallel%mydim,i)+ &
                    conjg(elec_st%eig(irp,kplp,1)%zwf(iit+parallel%mydim,i))* &
                    elec_st%eig(irp,kplp,1)%zwf(iit,i))**2 + &
                    (zi*(-conjg(elec_st%eig(irp,kplp,1)%zwf(iit,i))* &
                    elec_st%eig(irp,kplp,1)%zwf(iit+parallel%mydim,i) + &
                    conjg(elec_st%eig(irp,kplp,1)%zwf(iit+parallel%mydim,i))* &
                    elec_st%eig(irp,kplp,1)%zwf(iit,i)))**2 + &
                    ((abs(elec_st%eig(irp,kplp,1)%zwf(iit,i))**2 - &
                    abs(elec_st%eig(irp,kplp,1)%zwf(iit+parallel%mydim,i))**2))**2
                  else
                     magmom = magmom + &
                          abs(elec_st%eig(irp,kplp,1)%zwf(iit,i))**2 - &
                          abs(elec_st%eig(irp,kplp,1)% &
                          zwf(iit+parallel%mydim,i))**2
                  endif
               enddo
            endif
         call psum(magmom,1,parallel%procs_num,parallel%comm)
         elec_st%magmom(i,kplp) = magmom
         if (elec_st%ncl) magmom = sqrt(magmom)
      enddo
   endif
enddo
deallocate(nec,maxmvs,nrep,nrep_add,out_info,itmp,mxstt)

#ifdef MPI
call MPI_BARRIER(comm,mpinfo)
#endif
#ifdef ITAC
    call VTEND( vt_so,ierr)
#endif

end subroutine eigval
!===============================================================
!
! For a given set of eigenvalues (for all representations), sort
! them in ascending order. Use straight insertion, no need for
! fancy algorithms since the arrays are always "short".
!
!---------------------------------------------------------------
subroutine eigen_sort(elec_st,isp,kplp,nmax,nadd,nrep,nrep_add,irep)

  use constants
  use electronic_struct_module
  implicit none
  !
  ! Input/Output variables:
  !
  ! electronic structure
  type (electronic_struct), intent(in) :: elec_st
  ! spin orientation
  integer, intent(in) :: isp
  ! kpoint
  integer, intent(in) :: kplp
  ! size of representation array
  integer, intent(in) :: nmax
  ! maximum number of additional states per representation
  integer, intent(in) :: nadd
  ! lowest number of states for each representation, such that the
  ! total sums up to (elec_st%nstate-1) after eigenvalues are sorted
  integer, intent(out) :: nrep(elec_st%nrep)
  ! number of additional states on each representation, after
  ! removing the first (elec_st%nstate-1) ones
  integer, intent(out) :: nrep_add(elec_st%nrep)
  ! representation array, irep(i) = 0 for i > sum( nrep_add + nrep )
  integer, intent(out) :: irep(nmax)
  !
  ! Work variables:
  !
  ! counters, temporary variables
  integer i, irp, jj, l, ntotal, nc
  real(dp) :: e_tmp

  ! temporary arrays: representation, order index, eigenvalues
  integer, dimension(:) , allocatable :: rep_tmp, ord_tmp
  real(dp), dimension(:), allocatable :: eig_tmp
  !---------------------------------------------------------------

  ntotal = 0
  do jj = 1, size(elec_st%eig(:,kplp,isp))
     ntotal = ntotal + elec_st%eig(jj,kplp,isp)%nn
  enddo
  allocate(eig_tmp(ntotal))
  allocate(rep_tmp(ntotal))
  allocate(ord_tmp(ntotal))
  nc = 0
  do irp = 1, elec_st%nrep
     do i = 1, elec_st%eig(irp,kplp,isp)%nn
        nc = nc + 1
        eig_tmp(nc) = elec_st%eig(irp,kplp,isp)%en(i)
        rep_tmp(nc) = irp
        ord_tmp(nc) = i
     enddo
  enddo

  do l = 2, nc
     e_tmp = eig_tmp(l)
     irp = rep_tmp(l)
     i = ord_tmp(l)
     do jj = l-1, 1, -1
        if (eig_tmp(jj) <= e_tmp) goto 10
        eig_tmp(jj+1) = eig_tmp(jj)
        rep_tmp(jj+1) = rep_tmp(jj)
        ord_tmp(jj+1) = ord_tmp(jj)
     enddo
     jj = 0
10   continue
     eig_tmp(jj+1) = e_tmp
     rep_tmp(jj+1) = irp
     ord_tmp(jj+1) = i
  enddo

  nrep = 0
  irep = 0
  do jj = 1, elec_st%nstate-1
     irep(jj) = rep_tmp(jj)
     nrep(irep(jj)) = nrep(irep(jj)) + 1
  enddo

  nrep_add = 0
  i = elec_st%nstate-1
  do jj = elec_st%nstate , nc
     irp = rep_tmp(jj)
     if (nrep_add(irp) < nadd) then
        nrep_add(irp) = nrep_add(irp) + 1
        i = i + 1
        irep(i) = rep_tmp(jj)
     endif
  enddo

  deallocate(eig_tmp, rep_tmp, ord_tmp)

end subroutine eigen_sort
!===============================================================
!
! Adjust size of eigenstate arrays.
!
!---------------------------------------------------------------
subroutine eig_adj_size(eigen,ldn,new_size,cplx)

  use constants
  use electronic_struct_module
  implicit none
  !
  ! Input/Output variables:
  !
  ! eigenstates
  type (eigenstates), intent(inout) :: eigen
  ! size of grid local to this processor (i.e. local number of rows)
  integer, intent(in) :: ldn
  ! new size of arrays
  integer, intent(in) :: new_size
  ! real/complex flag
  logical, intent(in) :: cplx
  !
  ! Work variables:
  !
  ! allocation flag
  integer alcstat
  ! old size of arrays
  integer old_size
  ! number of saved states (useful if old_size > new_size)
  integer nn
  ! buffer arrays for eigenvectors, eigenvalues, occupancy factors
  real(dp), allocatable :: eig_buff(:),occ_buff(:)
  complex(dpc), allocatable :: zwf_buff(:,:)
  real(dp), allocatable :: wf_buff(:,:)
  logical :: save_flag

  !---------------------------------------------------------------

  old_size = 0
  if (associated(eigen%en)) old_size = size(eigen%en)

  ! If new_size is smaller than old_size, do nothing.
  if (old_size >= new_size) return
  !
  ! If the structure already has some data, save it in temporary arrays.
  !
  if (associated(eigen%en)) then
     nn = min(old_size,new_size)

     allocate(eig_buff(nn))
     eig_buff(1:nn) = eigen%en(1:nn)
     allocate(occ_buff(nn))
     occ_buff(1:nn) = eigen%occ(1:nn)

     if (cplx) then
        allocate(zwf_buff(ldn,nn),stat=alcstat)
        call alccheck('zwf_buff',ldn*nn,alcstat)
        zwf_buff(:,1:nn) = eigen%zwf(:,1:nn)
     else
        allocate(wf_buff(ldn,nn),stat=alcstat)
        call alccheck('wf_buff',ldn*nn,alcstat)
        wf_buff(:,1:nn) = eigen%wf(:,1:nn)
     endif
     save_flag = .true.
  else
     save_flag = .false.
  endif
  !
  ! Now, reallocate arrays.
  !
  call destroy_eigenstate (eigen)
  if (new_size == 0) then
     if (save_flag) deallocate(occ_buff, eig_buff)
     if (cplx) then
        if (save_flag) deallocate(zwf_buff)
     else
        if (save_flag) deallocate(wf_buff)
     endif
     return
  endif
  call create_eigenstate (eigen,ldn,new_size,cplx)
  !
  ! Transfer information from old arrays to resized arrays.
  !
  if (save_flag) then
     eigen%en(1:nn) = eig_buff(1:nn)
     eigen%occ(1:nn) = occ_buff(1:nn)
     if (cplx) then
        eigen%zwf(:,1:nn) = zwf_buff(:,1:nn)
        deallocate(zwf_buff)
     else
        eigen%wf(:,1:nn) = wf_buff(:,1:nn)
        deallocate(wf_buff)
     endif
     deallocate(occ_buff, eig_buff)
  endif

end subroutine eig_adj_size
!===============================================================
!
! Adjust size of eigenstate arrays.
!
!---------------------------------------------------------------
subroutine eig_e_adj_size(eigen,new_size)

  use constants
  use electronic_struct_module
  implicit none
  !
  ! Input/Output variables:
  !
  ! eigenstates
  type (eigenstates), intent(inout) :: eigen
  ! new size of arrays
  integer, intent(in) :: new_size
  !
  ! Work variables:
  !
  ! old size of arrays
  integer old_size
  ! number of saved states (useful if old_size > new_size)
  integer nn
  ! buffer arrays for eigenvectors, eigenvalues, occupancy factors
  real(dp), allocatable :: eig_buff(:),occ_buff(:)
  logical :: save_flag

  !---------------------------------------------------------------

  old_size = 0
  if (associated(eigen%en)) old_size = size(eigen%en)

  ! If new_size is smaller than old_size, do nothing.
  if (old_size >= new_size) return
  !
  ! If the structure already has some data, save it in temporary arrays.
  !
  if (associated(eigen%en)) then
     nn = min(old_size,new_size)

     allocate(eig_buff(nn))
     eig_buff(1:nn) = eigen%en(1:nn)
     allocate(occ_buff(nn))
     occ_buff(1:nn) = eigen%occ(1:nn)

     save_flag = .true.
  else
     save_flag = .false.
  endif
  !
  ! Now, reallocate arrays.
  !
  if (associated(eigen%en)) deallocate(eigen%en)
  if (associated(eigen%occ)) deallocate(eigen%occ)
  if (new_size == 0) then
     if (save_flag) deallocate(occ_buff, eig_buff)
     return
  endif
  allocate(eigen%en(new_size))
  eigen%en = zero
  allocate(eigen%occ(new_size))
  eigen%occ = zero
  !
  ! Transfer information from old arrays to resized arrays.
  !
  if (save_flag) then
     eigen%en(1:nn) = eig_buff(1:nn)
     eigen%occ(1:nn) = occ_buff(1:nn)
     deallocate(occ_buff, eig_buff)
  endif

end subroutine eig_e_adj_size
!===============================================================
!
! Initialize various arrays and variables before any eigensolver
! is called. This used to be performed at the beginning of
! subroutine eigval.F. (now called in calc_nscf.F90 and parsec.F90)?
!
!---------------------------------------------------------------
subroutine initeigval(elec_st,nadd)

  use constants
  use electronic_struct_module
  use eigen_solver_module
  use parallel_data_module
  implicit none
  !
  ! Input/Output variables:
  !
  ! electronic structure
  type (electronic_struct), intent(inout) :: elec_st
  ! number of additional eigenvalues per representation
  integer, intent(in) :: nadd
  !
  ! Work variables:
  !
  integer kpnum,kplp,iit,ntotal,ii,irp

  !---------------------------------------------------------------

  kpnum = max(elec_st%nkpt,1)

  ! assume a "blind guess" for the number of eigenstates to be
  ! computed; later, nn will be dynamically adjusted
#ifndef BETA
  ntotal = 0
#endif
  do iit = 1, elec_st%nspin/elec_st%mxwd
     do kplp = 1, kpnum               
#ifdef BETA
        ntotal = 0 ! I moved this line up .. (you better not)
#endif
        do irp = 1, elec_st%nrep
           ntotal = ntotal + elec_st%eig(irp,kplp,iit)%nec
        enddo
        if (ntotal < elec_st%nstate) then
           do irp = 1, elec_st%nrep
              ii = nadd + elec_st%nstate/elec_st%nrep
              if (elec_st%eig(irp,kplp,iit)%nn < ii) &
                   elec_st%eig(irp,kplp,iit)%nn = ii
           enddo
           elec_st%ntotal(iit) = elec_st%nstate
        endif
        ! calculate at least one eigenvector
        do irp = 1, elec_st%nrep
           if (elec_st%eig(irp,kplp,iit)%nn == 0) &
                elec_st%eig(irp,kplp,iit)%nn = 1
        enddo
     enddo
  enddo

  return

end subroutine initeigval
!===============================================================
