!===============================================================
!
!  Copyright (C) 2005 Finite Difference Research Group
!  This file is part of parsec, http://www.ices.utexas.edu/parsec/
!
!  Use Conjugate Gradient method to compute the Hartree Potential.
!  Input is the charge distribution rho(r) and an initial guess for
!  the output potential.
!  Output is the Hartree potential V(r), defined as:
!
!      Laplacian[ V(r) ] = -4*pi * 2 * rho(r)
!                        = integral_r' 2 * rho(r')/|r - r'|
!
!      4*pi factor comes from the Laplace equation.
!      Factor of 2 comes from the hartrees -> rydbergs unit conversion.
!
!---------------------------------------------------------------
subroutine hpotcg(parallel,rho,pot,norder,coeke,lap_dir_num,vect)

  use parallel_data_module

  use constants
#ifdef MPI
  !  include mpi definitions
  use mpi
#endif
  implicit none
  !
  !  Input/Output variables:
  !
  !  parallel computation related data
  type (parallel_data), intent(in) :: parallel

  !  order of expansion
  integer, intent(in) :: norder
  !  the input eletron distribution
  real(dp), intent(in) :: rho(parallel%mydim)
  !  lap_dir_num - an integer that is holding how much
  !  extra directions (besides the main 3) are there 
  !  for the laplacian.
  integer, intent(in)  :: lap_dir_num
  !  coefficients in finite difference expansion
  real(dp), intent(in) :: coeke(-norder:norder,lap_dir_num+3)
  !  pot == on entry it contains the initial guess of the potential;
  !         on output it contains the solution found.
  real(dp), intent(inout) :: pot(parallel%ldn)
  !  work array
  real(dp), intent(inout) :: vect(parallel%nwedge+1)
  !
  !  Work variables:
  !
  integer, parameter :: kss0 = 6

  integer i, ipar(16)
  !  allocation check
  integer alcstat

  real(dp), dimension(:), allocatable :: wk,wk2
  real(dp) :: fpar(16)
  !
  !  External functions:
  !
  real(dp), external :: dnrm2
  !---------------------------------------------------------------
  !
  !  Set up the parameter arrays.
  !
  ipar(1) = 0 !init solver 
  ipar(2) = 0 !no preconditioning
  ipar(3) = 1 !stopping criteria
  ipar(4) = kss0*parallel%mydim !number of elements in w array (CG=5*n)
  ipar(5) = 5 !size of krylov subspsace (not needed for CG)
  ipar(6) = 1600 !new max num matvec operations (was:800)
  !fpar(1) = 1.0d-6 !relative tolerance
  fpar(1) = 1.0d-7 !relative tolerance
  !fpar(2) = 1.0d-12 !absolute tolerance
  fpar(2) = 1.0d-13 !absolute tolerance
  fpar(11) = zero !initialization for correct flops count, I think it should be integer

  ipar(1) = 0
  fpar(11) = zero
  allocate(wk(kss0*parallel%mydim),stat=alcstat)
  call alccheck('hpotcg_wk',kss0*parallel%mydim,alcstat)
  allocate(wk2(parallel%ldn+1),stat=alcstat)
  call alccheck('hpotcg_wk2',parallel%ldn+1,alcstat)
  wk2(:) = zero
  wk(:) = zero

  do
  !calls cg from sparskit 
     call cg(parallel%mydim,kss0,rho,pot,ipar,fpar,wk &
          ,parallel%group_size,parallel%group_comm)
     if (ipar(1) == 1) then
        call dcopy(parallel%mydim, wk(ipar(8)), 1, wk2, 1)
        call lapmvs(parallel,wk2,wk(ipar(9)),coeke,norder,lap_dir_num,vect)
        fpar(11) = fpar(11) + 74*parallel%mydim !why 74 ?
        cycle
     elseif (ipar(1) < 0) then
        if (ipar(1) == 0) then
           write(9,*) 'Iterative possion solver has satisfied convergence test.'
        else if (ipar(1) == -1) then
           write(9,*) 'WARNING: possion solver exceeded max iterations'
        else if (ipar(1) == -2) then
           write(9,*) 'Oops: Iterative poisson solver was not given enough work space.'
           write(9,*) 'The work space should at least have ', ipar(4), &
                ' elements.'
        else if (ipar(1) == -3) then
           write(9,*) 'Iterative solver is facing a (nervous) break-down.'
        else
           write(9,*) 'Iterative solver was somehow terminated. code =', ipar(1)
        endif
     endif
     exit
  enddo

  write(9, '(a,i3,a,g12.4)') '# Poisson solver, return code = ', &
       ipar(1),'  convergence rate = ', fpar(7)
  !
  !  Check the error.
  !
  call lapmvs(parallel,pot,wk,coeke,norder,lap_dir_num,vect)
  do i = 1, parallel%mydim
     wk(i) = wk(i) - rho(i)
  enddo
  ! write (9, *) '# the residual norm ', &
  !      real(dnrm2(parallel%mydim,wk,1),dp),real(fpar(5),dp)
   write(9,*) '# reported residual norm: ',real(fpar(5),dp)
   write (9, *) '# euclidean norm of lap(pot): ', &
        real(dnrm2(parallel%mydim,wk,1),dp)
  if (allocated(wk)) deallocate(wk)
  if (allocated(wk2)) deallocate(wk2)
#ifdef MPI
  call MPI_Barrier(parallel%comm,i)
#endif

end subroutine hpotcg
!===============================================================
!
!  Dot product, MPI version of the BLAS ddot function.
!
!---------------------------------------------------------------
function distdot(n,x,ix,y,iy,nnodes,comm)

  use constants
#ifdef MPI
  !  include mpi definitions
  use mpi
#endif
  implicit none
  !
  !  Input/Output variables:
  !
  integer, intent(in) :: n, ix, iy, nnodes, comm
  real(dp), intent(in) :: x(n), y(n)
  !  dot product
  real(dp) :: distdot
  !
  !  Work variables:
  !
  real(dp) :: d
  !
  !  External functions:
  !
  real(dp), external :: ddot
  !---------------------------------------------------------------
  d = ddot(n,x,ix,y,iy)
  call psum(d,1,nnodes,comm)
  distdot = d

end function distdot
!===============================================================
!
!  Distributed Laplacian operator. This subroutine receives as
!  input a distributed function over the grid (array func) and
!  operates the laplacian over it, using finite differences.
!  Output is distributed function lapl = -Laplacian(func).
!
!  WARNING: notice that the output laplacian has a minus sign!
!
!  WARNING: input function func is assumed to have the full
!  symmetry of the system, so there are no representation
!  characters (cf. subroutine matvec).
!
!---------------------------------------------------------------
subroutine lapmvs(parallel, func, lapl, coeke, norder, lap_dir_num, vec)

  use parallel_data_module

  use constants
#ifdef MPI
  !  include mpi definitions
  use mpi
#endif
  implicit none
#ifdef ITAC
  include 'VT.inc'
  include 'vtcommon.inc'
#endif
  !
  !  Input/Output variables:
  !
  !  parallel computation related data
  type (parallel_data), intent(in) :: parallel
  !  order of finite-difference expansion
  integer, intent(in) :: norder
  !  lap_dir_num - number of extra directions for the laplacian.
  integer, intent(in) :: lap_dir_num
  !  finite-difference coefficients for a second derivative
  real(dp), intent(in) :: coeke(-norder:norder,lap_dir_num+3)
  !  input function
  real(dp), intent(in) :: func(parallel%mydim)
  !  output function, lapl = -Laplacian(func)
  real(dp), intent(out) :: lapl(parallel%mydim)
  !  work array
  real(dp), intent(out) :: vec(parallel%nwedge+1)
  !
  !  Work variables:
  !
  !  local scalars
  integer :: i, jj, ish, irow,  ioffset
  real(dp) :: tmp
  !  maxdim = maximum number of rows handled by each proc
  !  needed for MPI communication inside bdxc_as
  integer :: maxdim
#ifdef ITAC
  integer :: lapmverr
#endif

  !---------------------------------------------------------------

  maxdim = 0
  do i = 0, parallel%group_size - 1
     irow = parallel%irows(i+1)-parallel%irows(i)
     if (irow > maxdim) maxdim = irow
  enddo

  !  Global dimension and local rows' offset.
  ioffset = parallel%irows(parallel%group_iam) - 1
  !
  !  Initiate the boundary exchange information for laplacian.
  !
  ! AJB: switch with buffers, not because you need the buffer,
  ! but because it will be with the neighborhood collectives
  ! and you are copying to a buffer anyways
  call dcopy(parallel%mydim,func,1,vec,1)
#ifdef MPI
  call bdxc_as( vec, maxdim, parallel%nwedge, parallel%irecvp, &
       parallel%jsendp,parallel%senrows, parallel%group_comm, &
       parallel%group_size, parallel%group_iam)
#endif
  !
  !  Set wave function outside the domain be zero.
  !  The index for points outside the domain is ndim+1.
  !
  vec(parallel%nwedge+1) = zero 
  !
  !  Diagonal part. AJB: time for Openmp
  !
  tmp = sum(coeke(0,1:3+lap_dir_num))
  do i = 1, parallel%mydim
     lapl(i) = tmp * vec(i)
  enddo
  !
  !  Non-local part, symmetric function (no characters).
  ! 
  !  AJB: time to do this with openmp. 
  !  next - would be to use the same laplacian for
  !  matvec and lapmvs.
  !
#ifdef ITAC
    call VTBEGIN(vt_lapmvs , lapmverr)
#endif
  select case (lap_dir_num)
  case (0)
     do i = 1,parallel%mydim
        irow = parallel%pint(i)
        jj = 0
        do ish = 1,norder
           lapl(irow) = lapl(irow) + &
                coeke(ish,1)*(vec( parallel%neibs(jj+1,irow) ) + &
                              vec( parallel%neibs(jj+2,irow) )) + &
                coeke(ish,2)*(vec( parallel%neibs(jj+3,irow) ) + &
                              vec( parallel%neibs(jj+4,irow) )) + &
                coeke(ish,3)*(vec( parallel%neibs(jj+5,irow) ) + &
                              vec( parallel%neibs(jj+6,irow) ))
           jj = jj + 6
        enddo
     enddo
  case (1)
     do i = 1,parallel%mydim
        irow = parallel%pint(i)
        jj = 0
        do ish = 1,norder
           lapl(irow) = lapl(irow) + &
                coeke(ish,1)*(vec( parallel%neibs(jj+1,irow) ) + &
                              vec( parallel%neibs(jj+2,irow) )) + &
                coeke(ish,2)*(vec( parallel%neibs(jj+3,irow) ) + &
                              vec( parallel%neibs(jj+4,irow) )) + &
                coeke(ish,3)*(vec( parallel%neibs(jj+5,irow) ) + &
                              vec( parallel%neibs(jj+6,irow) )) + &
                coeke(ish,4)*(vec( parallel%neibs(jj+7,irow) ) + &
                              vec( parallel%neibs(jj+8,irow) ))
           jj = jj + 8
        enddo
     enddo
  case (2)
     do i = 1,parallel%mydim
        irow = parallel%pint(i)
        jj = 0
        do ish = 1,norder
           lapl(irow) = lapl(irow) + &
                coeke(ish,1)*(vec( parallel%neibs(jj+1,irow) ) + &
                              vec( parallel%neibs(jj+2,irow) )) + &
                coeke(ish,2)*(vec( parallel%neibs(jj+3,irow) ) + &
                              vec( parallel%neibs(jj+4,irow) )) + &
                coeke(ish,3)*(vec( parallel%neibs(jj+5,irow) ) + &
                              vec( parallel%neibs(jj+6,irow) )) + &
                coeke(ish,4)*(vec( parallel%neibs(jj+7,irow) ) + &
                              vec( parallel%neibs(jj+8,irow) )) + &
                coeke(ish,5)*(vec( parallel%neibs(jj+9,irow) ) + &
                              vec( parallel%neibs(jj+10,irow) ))
           jj = jj + 10
        enddo
     enddo
  case (3)
     do i = 1,parallel%mydim
        irow = parallel%pint(i)
        jj = 0
        do ish = 1,norder
           lapl(irow) = lapl(irow) + &
                coeke(ish,1)*(vec( parallel%neibs(jj+1,irow) ) + &
                              vec( parallel%neibs(jj+2,irow) )) + &
                coeke(ish,2)*(vec( parallel%neibs(jj+3,irow) ) + &
                              vec( parallel%neibs(jj+4,irow) )) + &
                coeke(ish,3)*(vec( parallel%neibs(jj+5,irow) ) + &
                              vec( parallel%neibs(jj+6,irow) )) + &
                coeke(ish,4)*(vec( parallel%neibs(jj+7,irow) ) + &
                              vec( parallel%neibs(jj+8,irow) )) + &
                coeke(ish,5)*(vec( parallel%neibs(jj+9,irow) ) + &
                              vec( parallel%neibs(jj+10,irow) )) + &
                coeke(ish,6)*(vec( parallel%neibs(jj+11,irow) ) + &
                              vec( parallel%neibs(jj+12,irow) ))
           jj = jj + 12
        enddo
     enddo
  end select
#ifdef ITAC
    call VTEND(vt_lapmvs , lapmverr)
#endif

end subroutine lapmvs
!===============================================================
!
!  Distributed Gradient operator. This subroutine receives as
!  input a distributed function over the grid (array func) and
!  operates the Gradient over it, using finite differences.
!
!  WARNING: input function is assumed to have the full symmetry
!  of the system, so there are no representation characters (cf.
!  subroutine matvec).
!
!---------------------------------------------------------------
subroutine gradmvs(grid,parallel, func, gradf, coeke, norder, vec)

  use grid_module
  use parallel_data_module

  use constants
#ifdef MPI
  !  include mpi definitions
  use mpi
#endif
  implicit none
  !
  !  Input/Output variables:
  !
  !  grid related data
  type (grid_data), intent(in) :: grid
  !  parallel computation related data
  type (parallel_data), intent(in) :: parallel
  !  order of finite-difference expansion
  integer, intent(in) :: norder
  !  finite-difference coefficients for a first derivative
  real(dp), intent(in) :: coeke(-norder:norder,3)
  !  input function
  real(dp), intent(in) :: func(parallel%mydim)
  !  output gradient, gradf(:,i) = Gradient [ func(i) ]
  real(dp), intent(out) :: gradf(3,parallel%mydim)
  !  work array
  real(dp), intent(out) :: vec(parallel%nwedge+1)
  !
  !  Work variables:
  !
  !  local scalars
  integer i, ish, irow,  ioffset, kk
  !  maxdim = maximum number of rows handled by each proc
  !  needed for MPI communication inside bdxc_as
  integer maxdim
  real(dp), dimension(3) :: tmpvec
  !  work variable for extra directions
  integer :: lap_dir_nump6

  !---------------------------------------------------------------

  maxdim = 0
  do i = 0, parallel%group_size - 1
     irow = parallel%irows(i+1)-parallel%irows(i)
     if (irow > maxdim) maxdim = irow
  enddo

  !  Global dimension and local rows' offset
  ioffset = parallel%irows(parallel%group_iam) - 1
  !
  !  Initiate the boundary exchange information
  !
  call dcopy(parallel%mydim,func,1,vec,1)
  ! this sends too many neighbors - since it expects the laplacian operator
#ifdef MPI
  call bdxc_as( vec, maxdim, parallel%nwedge, parallel%irecvp, &
       parallel%jsendp,parallel%senrows, parallel%group_comm, &
       parallel%group_size, parallel%group_iam)
#endif
  !
  !  Set wave function outside the domain be zero.
  !  The index for points outside the domain is ndim+1.
  !
  vec(parallel%nwedge+1) = zero 
  !
  !  Finite difference.
  !
  gradf(:,:) = zero
  if (grid%lap_dir_num > 0) then
     lap_dir_nump6 = 6 + grid%lap_dir_num * 2
     do i = 1,parallel%mydim
        tmpvec = 0
        irow = parallel%pint(i)
        kk = 0
        do ish = 1,norder
           tmpvec(1) = tmpvec(1) + &
              coeke(ish,1)*(vec( parallel%neibs(kk+2,irow) ) - &
                            vec( parallel%neibs(kk+1,irow) ))
           tmpvec(2) = tmpvec(2) + &
              coeke(ish,2)*(vec( parallel%neibs(kk+4,irow) ) - &
                            vec( parallel%neibs(kk+3,irow) ))
           tmpvec(3) = tmpvec(3) + &
              coeke(ish,3)*(vec( parallel%neibs(kk+6,irow) ) - &
                            vec( parallel%neibs(kk+5,irow) ))
           kk = kk + lap_dir_nump6
        enddo
        call matvec3('N',grid%grad_bvec_norm,tmpvec,gradf(1,irow))
     enddo
  else
     do i = 1,parallel%mydim
        tmpvec = 0
        irow = parallel%pint(i)
        kk = 0
        do ish = 1,norder
           tmpvec(1) = tmpvec(1) + &
              coeke(ish,1)*(vec( parallel%neibs(kk+2,irow) ) - &
                            vec( parallel%neibs(kk+1,irow) ))
           tmpvec(2) = tmpvec(2) + &
              coeke(ish,2)*(vec( parallel%neibs(kk+4,irow) ) - &
                            vec( parallel%neibs(kk+3,irow) ))
           tmpvec(3) = tmpvec(3) + &
              coeke(ish,3)*(vec( parallel%neibs(kk+6,irow) ) - &
                            vec( parallel%neibs(kk+5,irow) ))
           kk = kk + 6
        enddo
        gradf(1:3,irow) = tmpvec
     enddo
  endif

end subroutine gradmvs
!===============================================================
