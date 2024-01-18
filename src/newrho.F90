!===============================================================
!
! Copyright (C) 2005 Finite Difference Research Group
! This file is part of parsec, http://www.ices.utexas.edu/parsec/
!
! Updates the charge density based on the new wave functions that
! has been found after a self-consistency iteration.
!
! AJB 2015: NEEDS WORK.
!---------------------------------------------------------------
subroutine newrho(elec_st,parallel,hcub2)

  use constants
  use electronic_struct_module
  use parallel_data_module
#ifdef MPI
  use mpi
#endif
  implicit none
  !
  ! Input/Output variables:
  !
  ! electronic structure
  type (electronic_struct), intent(inout) :: elec_st
  ! parallel computation related data
  type (parallel_data), intent(in) :: parallel

  ! 2/(grid spacing)**3
  real(dp), intent(in) :: hcub2
  !
  ! Work variables:
  !
  ! 2/(h^3) if nspin=1, 1/(h^3) if nspin=2 (to account for two
  ! electrons or one electron in a filled state, respectively)
  real(dp) :: factor
  ! charge density
  real(dp) ::  rho(parallel%ldn,2*elec_st%nspin-1)
  ! counters
  integer i,isp,j,jj,irp
#ifdef MPI
  integer mpinfo,nelems
#endif
  real(dp) :: dnorm,tmp
  complex(dpc) :: znorm
  ! spin-orbit variables
  integer j2
  integer kplp, kpnum

  !---------------------------------------------------------------

  kpnum = max(elec_st%nkpt,1)

! #ifdef MPI
!   call MPI_BARRIER(parallel%comm,mpinfo)
! #endif
#ifdef MPI
  do isp = 1,elec_st%nspin/elec_st%mxwd
     call MPI_BCAST(elec_st%totel(isp),1,MPI_DOUBLE_PRECISION, &
          parallel%masterid,parallel%comm,mpinfo)
     call MPI_BCAST(elec_st%ifmax(isp),1, &
          MPI_INTEGER,parallel%masterid,parallel%comm,mpinfo)
     if(.not.elec_st%ncl)then
       do kplp = 1, kpnum
          do irp = 1, elec_st%nrep
             nelems = elec_st%eig(irp,kplp,isp)%nn
             if (nelems == 0) cycle
             call MPI_BCAST(elec_st%eig(irp,kplp,isp)%occ,nelems, &
                MPI_DOUBLE_PRECISION,parallel%masterid, &
                parallel%comm,mpinfo)
          enddo
       enddo
     endif
  enddo
  call MPI_BCAST(elec_st%totel,elec_st%nspin, &
       MPI_DOUBLE_PRECISION,parallel%masterid, &
       parallel%comm,mpinfo)
#endif

  factor = hcub2/ real(elec_st%nspin,dp)
  rho(:,:) = zero

  ! Update charge density for both spin orientations.
  ! If elec_st%nspin = 1, update the total charge density.
  ! Also, rescale eigenvectors: norm of each vector over the *full*
  ! grid is 1, but eigensolvers use norm = 1 over the irreducible
  ! wedge. Must divide that by the ratio of volumes full/wedge.
  !
  do isp = 1, elec_st%nspin/elec_st%mxwd
     jj = isp - 1 + elec_st%nspin
     do kplp = 1, kpnum
        do irp = 1, elec_st%nrep
           !
           ! Only processors belonging to the group in charge of this
           ! representation/k-point/spin should work.
           !
           if ( elec_st%eig(irp,kplp,isp)%group /= &
                parallel%mygroup ) cycle
 
           if (elec_st%cplx) then
              znorm = zone/sqrt( real(elec_st%nrep,dp) )
              i = size(elec_st%eig(irp,kplp,isp)%zwf)
              call zscal(i,znorm,elec_st%eig(irp,kplp,isp)%zwf(1,1),1)

              if (elec_st%mxwd == 2) then

!             AMIR - have to check this change with Doron..

!              if (elec_st%is_so) then

                 do i = 1, elec_st%eig(irp,kplp,isp)%nec
                    do j = 1, parallel%mydim
                       rho(j,2) = rho(j,2) +  &
                            elec_st%eig(irp,kplp,1)%occ(i) &
                            *abs(elec_st%eig(irp,kplp,1)%zwf(j,i))**2 &
                            *elec_st%kpwt(kplp)
                       
                       j2 = j + parallel%mydim
                       
                       rho(j,3) = rho(j,3) +  &
                            elec_st%eig(irp,kplp,1)%occ(i) &
                            *abs(elec_st%eig(irp,kplp,1)%zwf(j2,i))**2 &
                            *elec_st%kpwt(kplp)
                    enddo
                 enddo
              else
                 do i = 1, elec_st%eig(irp,kplp,isp)%nec
                    do j=1,parallel%mydim
                       rho(j,jj) = rho(j,jj)+ &
                            elec_st%eig(irp,kplp,isp)%occ(i) &
                            *abs(elec_st%eig(irp,kplp,isp)%zwf(j,i))**2 &
                            *elec_st%kpwt(kplp)
                    enddo
                 enddo
              endif
           else
              dnorm = one/sqrt( real(elec_st%nrep,dp) )
              i = size(elec_st%eig(irp,kplp,isp)%wf)
              call dscal(i,dnorm,elec_st%eig(irp,kplp,isp)%wf(1,1),1)
              do i = 1, elec_st%eig(irp,kplp,isp)%nec
                 do j=1,parallel%mydim
                    rho(j,jj) = &
                         rho(j,jj) + elec_st%eig(irp,kplp,isp)%occ(i) &
                         *abs(elec_st%eig(irp,kplp,isp)%wf(j,i))**2
                 enddo
              enddo
           endif
        enddo
     enddo
  enddo
  if(elec_st%ncl)then
     do j=1,parallel%mydim
        rho(j,1) = rho(j,2) + rho(j,3)

        rho(j,2) = half*(rho(j,1)+elec_st%spin3dn(j))

        rho(j,3) = half*(rho(j,1)-elec_st%spin3dn(j))
     enddo
  endif

  do jj = 1, 2*elec_st%nspin - 1
     do j=1,parallel%mydim
        rho(j,jj) = factor*rho(j,jj)
     enddo
     !
     ! Perform reduction across groups
     !
     call group_reduction(parallel,rho(1,jj))
  enddo
  ! update the total charge density, if elec_st%nspin = 2
  if (elec_st%nspin == 2 .and. (.not. elec_st%ncl)) then
     do j=1,parallel%mydim
        rho(j,1) = rho(j,2) + rho(j,3)
     enddo
  endif
!Consider removing:
#ifdef MPI
  call MPI_BARRIER(parallel%comm,mpinfo)
#endif
  do jj = 1, 2*elec_st%nspin - 1
     call dcopy(parallel%mydim,rho(1,jj),1,elec_st%rho(1,jj),1)
  enddo

end subroutine newrho
!===============================================================
!
! Copyright (C) 2005 Finite Difference Research Group
! This file is part of parsec, http://www.ices.utexas.edu/parsec/
!
! Symmetrize a scalar function according to symmetry
! operations.
!
! Author: Murilo Tiago (2005)
!
!---------------------------------------------------------------
subroutine symm_scalar(grid,rsymm,symm,parallel,fr,is_pbc)

  use constants
  use grid_module
  use symmetry_module
  use parallel_data_module
#ifdef MPI
  use mpi
#endif
  implicit none
  ! include mpi definitions
  !
  ! Input/Output variables:
  !
  ! grid related data
  type (grid_data), intent(inout) :: grid
  ! symmetry operations
  type (symmetry), intent(in) :: rsymm,symm
  ! parallel computation related data
  type (parallel_data), intent(inout) :: parallel
  ! function to be symmetrized
  real(dp), intent(inout) :: fr(parallel%mydim) 
  ! true if this system is periodic
  logical, intent(in) :: is_pbc
  !
  ! Work variables:
  !
  ! tolerance in the position of one atom and its equivalent by
  ! symmetry
  real(dp), parameter :: tol = 1.d-6
  ! buffer function
  real(dp) :: ftmp(parallel%ldn)
  ! symmetrized function
  real(dp) :: fsym(parallel%ldn)
  ! counters
  integer :: ii,jj,igrid,isym,ipe,jpe,ngrid,ioffset,ir(3),joffset
  real(dp) :: fin,dmax,rw(3),rrw(3),favg,tmp
  ! timers
  real(dp) t0,t1,t2
  !
  ! External function:
  !
  integer, external :: Part
#ifdef MPI
  integer :: mpinfo
#endif

  !---------------------------------------------------------------

     call mysecond(t0)
  call grid_set_ist(grid,parallel%group_size)

  fsym(:) = zero
  do ipe = 0, parallel%group_size - 1
     ftmp(:) = zero
     if (ipe == parallel%group_iam) then
        call dcopy(parallel%mydim,fr,1,ftmp,1)
        ngrid = parallel%mydim
     endif
#ifdef MPI
     call MPI_BCAST(ngrid,1, &
          MPI_INTEGER,ipe,parallel%group_comm,mpinfo)
     call MPI_BCAST(ftmp,parallel%ldn, &
          MPI_DOUBLE_PRECISION,ipe,parallel%group_comm,mpinfo)
#endif
     ioffset = parallel%irows(ipe) - 1
     joffset = parallel%irows(parallel%group_iam) - 1
     do ii = 1, ngrid
        fin = ftmp(ii)

        if(is_pbc) then
           rw(1) = (grid%shift(1) + grid%kx(ii+ioffset)) &
                / real(grid%n1,dp)
           rw(2) = (grid%shift(2) + grid%ky(ii+ioffset)) &
                / real(grid%n2,dp)
           rw(3) = (grid%shift(3) + grid%kz(ii+ioffset)) &
                / real(grid%n3,dp)
        else
           rw(1) = (grid%shift(1) + grid%kx(ii+ioffset)) &
                * grid%step(1)
           rw(2) = (grid%shift(2) + grid%ky(ii+ioffset)) &
                * grid%step(2)
           rw(3) = (grid%shift(3) + grid%kz(ii+ioffset)) &
                * grid%step(3)
        endif
        do isym = 1, symm%ntrans
           if (is_pbc) then
              call symop_relative(symm,isym,rw,rrw)
              ! Must make sure that ir is in the right cell:
              !   -n/2 <= ir < n/2
              ir(1) = nint(rrw(1)*real(grid%n1,dp) - grid%shift(1))
              ir(2) = nint(rrw(2)*real(grid%n2,dp) - grid%shift(2))
              ir(3) = nint(rrw(3)*real(grid%n3,dp) - grid%shift(3))
              ir(1) = mod( ir(1) + grid%n1*15/2, grid%n1) - grid%n1/2
              ir(2) = mod( ir(2) + grid%n2*15/2, grid%n2) - grid%n2/2
              ir(3) = mod( ir(3) + grid%n3*15/2, grid%n3) - grid%n3/2
           else
              call symop(symm,isym,rw,rrw)
              do jj = 1, 3
                 ir(jj)=nint(rrw(jj)/grid%step(jj)-grid%shift(jj))
              enddo
           endif
           jpe = Part(grid,parallel%group_size,ir(1), &
                ir(2),ir(3),.false.)
           if (jpe == parallel%group_iam) then
              igrid = grid%indexw(ir(1),ir(2),ir(3))
              fsym(igrid-joffset) = fsym(igrid-joffset) + fin
           endif
        enddo
     enddo
! AJB: no need for a barrier here.
! #ifdef MPI
!      call MPI_Barrier(parallel%group_comm,mpinfo)
! #endif
  enddo
     call mysecond(t1)
     write(9,*) 
     write(9,101) t1-t0
101  format(' symm_scalar fsym collection [sec]:',1x,f10.2)

  fsym = fsym / real(symm%ntrans,dp) * real(rsymm%ntrans,dp)
  call dcopy(parallel%mydim,fr,1,ftmp,1)
  dmax = maxval(abs(fsym-ftmp))
  favg = sum(fr)
#ifdef MPI
  call MPI_ALLREDUCE(favg,tmp,1,MPI_DOUBLE_PRECISION,MPI_SUM, &
       parallel%group_comm,mpinfo)
  favg = tmp
  call MPI_ALLREDUCE(dmax,tmp,1,MPI_DOUBLE_PRECISION,MPI_MAX, &
       parallel%group_comm,mpinfo)
  dmax = tmp
#endif
  if (parallel%iammaster) then
     write(7,1) dmax/favg*100.0
     write(7,2) favg*grid%hcub*real(rsymm%ntrans,dp)
  endif
  call dcopy(parallel%mydim,fsym,1,fr,1)

  favg = sum(fr)
#ifdef MPI
  call MPI_ALLREDUCE(favg,tmp,1,MPI_DOUBLE_PRECISION,MPI_SUM, &
       parallel%group_comm,mpinfo)
  favg = tmp
#endif
  if (parallel%iammaster) then
     write(7,3) favg*grid%hcub*real(rsymm%ntrans,dp)
  endif
     call mysecond(t2)
     write(9,*) 
     write(9,102) t2-t0
102  format(' symm_scalar total time [sec]:',1x,f10.2)

1 format(' maximum deviation after symmetrization: ',f10.4,' %')
2 format(' electron charge before symmetrization: ',f10.4)
3 format(' electron charge after symmetrization: ',f10.4)

end subroutine symm_scalar
!===============================================================
