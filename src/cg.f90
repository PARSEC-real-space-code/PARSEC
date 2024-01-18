!===============================================================
!
!  Copyright (C) 2005 Finite Difference Research Group
!  This file is part of parsec, http://www.ices.utexas.edu/parsec/
!
!                       S P A R S K I T
!  ---------------------------------------------------------------
!        Basic Iterative Solvers with Reverse Communication
!  ---------------------------------------------------------------
!  This file currently has several basic iterative linear system
!  solvers. They are:
!  CG       -- Conjugate Gradient Method
!  CGNR     -- Conjugate Gradient Method (Normal Residual equation)
!  BCG      -- Bi-Conjugate Gradient Method
!  DBCG     -- BCG with partial pivoting
!  BCGSTAB  -- BCG stabilized
!  TFQMR    -- Transpose-Free Quasi-Minimum Residual method
!  FOM      -- Full Orthogonalization Method
!  GMRES    -- Generalized Minimum RESidual method
!  FGMRES   -- Flexible version of Generalized Minimum
!              RESidual method
!  DQGMRES  -- Direct versions of Quasi Generalize Minimum
!              Residual method
!  ---------------------------------------------------------------
!  They all have the following calling sequence:
!     subroutine solver(n, rhs, sol, ipar, fpar, w)
!     integer n, ipar(16)
!     real*8 rhs(n), sol(n), fpar(16), w(*)
!  Where
!  (1) 'n' is the size of the linear system,
!  (2) 'rhs' is the right-hand side of the linear system,
!  (3) 'sol' is the solution to the linear system,
!  (4) 'ipar' is an integer parameter array for the reverse
!  communication protocol,
!  (5) 'fpar' is an floating-point parameter array storing
!  information to and from the iterative solvers.
!  (6) 'w' is the work space (size is specified in ipar)
!
!  They are preconditioned iterative solvers with reverse
!  communication. The preconditioners can be applied from either
!  from left or right or both (specified by ipar(2), see below).

!  Author: Kesheng John Wu (kewu@mail.cs.umn.edu) 1993
!
!  NOTES:
!
!  (1) Work space required by each of the iterative solver
!  routines is as follows:
!    CG      == 5 * n
!    CGNR    == 5 * n
!    BCG     == 7 * n
!    DBCG    == 11 * n
!    BCGSTAB == 8 * n
!    TFQMR   == 11 * n
!    FOM     == (n+3)*(m+2) + (m+1)*m/2 (m = ipar(5), default m=15)
!    GMRES   == (n+3)*(m+2) + (m+1)*m/2 (m = ipar(5), default m=15)
!    FGMRES  == n*(2m+1) + (m+1)*m/2 + 3*m + 2 (m = ipar(5),
!               default m=15)
!    DQGMRES == n + lb * (2*n+4) (lb=ipar(5)+1, default lb = 16)
!
!  (2) ALL iterative solvers require a user-supplied DOT-product
!  routine named DISTDOT. The prototype of DISTDOT is
!
!  real*8 function distdot(n,x,ix,y,iy)
!  integer n, ix, iy
!  real*8 x(1+(n-1)*ix), y(1+(n-1)*iy)
!
!  This interface of DISTDOT is exactly the same as that of
!  DDOT (or SDOT if real == real*8) from BLAS-1. It should have
!  same functionality as DDOT on a single processor machine. On a
!  parallel/distributed environment, each processor can perform
!  DDOT on the data it has, then perform a summation on all the
!  partial results.
!
!  (3) To use this set of routines under SPMD/MIMD program
!  paradigm, several things are to be noted: (a) 'n' should be the
!  number of vector elements of 'rhs' that is present on the local
!  processor .
!  (b) if RHS(i) is on processor j, it is expected that SOL(i)
!  will be on the same processor, i.e. the vectors are distributed
!  to each processor in the same way. (c) the preconditioning and
!  stopping criteria specifications have to be the same on all
!  processor involved, ipar and fpar have to be the same on each
!  processor. (d) DISTDOT should be replaced by a distributed
!  dot-product function.
!
!  ..............................................................
!  Reverse Communication Protocols
!
!  When a reverse-communication routine returns, it could be
!  either that the routine has terminated or it simply requires
!  the caller to perform one matrix-vector multiplication. The
!  possible matrices that involve in the matrix-vector
!  multiplications are:
!  A       (the matrix of the linear system),
!  A^T     (A transposed),
!  Ml^{-1} (inverse of the left preconditioner),
!  Ml^{-T} (inverse of the left preconditioner transposed),
!  Mr^{-1} (inverse of the right preconditioner),
!  Mr^{-T} (inverse of the right preconditioner transposed).
!  For all the matrix vector multiplication, v = A u. The input
!  and output vectors are supposed to be part of the work space 'w',
!  and the starting positions of them are stored in ipar(8:9), see
!  below.
!
!  The array 'ipar' is used to store the information about the
!  solver. Here is the list of what each element represents:
!
!  ipar(1) -- status of the call/return.
!  A call to the solver with ipar(1) == 0 will initialize the
!  iterative solver. On return from the iterative solver, ipar(1)
!  carries the status flag which indicates the condition of the
!  return. The status information is divided into two categories,
!  (1) a positive value indicates the solver requires a
!  matrix-vector multiplication,
!  (2) a non-positive value indicates termination of the solver.
!  Here is the current definition:
!    1 == request a matvec with A,
!    2 == request a matvec with A^T,
!    3 == request a left preconditioner solve (Ml^{-1}),
!    4 == request a left preconditioner transposed solve (Ml^{-T}),
!    5 == request a right preconditioner solve (Mr^{-1}),
!    6 == request a right preconditioner transposed solve (Mr^{-T}),
!   10 == request the caller to perform stopping test,
!    0 == normal termination of the solver, satisfied the stopping
!         criteria,
!   -1 == termination because iteration number is greater than the
!         preset limit,
!   -2 == return due to insufficient work space,
!   -3 == return due to anticipated break-down / divide by zero,
!         in the case where Arnoldi procedure is used, additional
!         error code can be found in ipar(12), where ipar(12) is
!         the error code of orthogonalization procedure MGSRO:
!            -1: zero input vector
!            -2: input vector contains abnormal numbers
!            -3: input vector is a linear combination of others
!            -4: trianguler system in GMRES/FOM/etc. has rank 0 (zero)
!   -4 == the values of fpar(1) and fpar(2) are both <= 0, the valid
!         ranges are 0 <= fpar(1) < 1, 0 <= fpar(2), and they can
!         not be zero at the same time
!   -9 == while trying to detect a break-down, an abnormal number is
!         detected.
!  -10 == return due to some non-numerical reasons, e.g. invalid
!         floating-point numbers etc.
!
!  ipar(2) -- status of the preconditioning:
!    0 == no preconditioning
!    1 == left preconditioning only
!    2 == right preconditioning only
!    3 == both left and right preconditioning
!
!  ipar(3) -- stopping criteria (details of this will be
!  discussed later).
!
!  ipar(4) -- number of elements in the array 'w'. if this is less
!  than the desired size, it will be over-written with the minimum
!  requirement. In which case the status flag ipar(1) = -2.
!
!  ipar(5) -- size of the Krylov subspace (used by GMRES and its
!  variants), e.g. GMRES(ipar(5)), FGMRES(ipar(5)),
!  DQGMRES(ipar(5)).
!
!  ipar(6) -- maximum number of matrix-vector multiplies, if not a
!  positive number the iterative solver will run till convergence
!  test is satisfied.
!
!  ipar(7) -- current number of matrix-vector multiplies. It is
!  incremented after each matrix-vector multiplication. If there
!  is preconditioning, the counter is incremented after the
!  preconditioning associated with each matrix-vector
!  multiplication.
!
!  ipar(8) -- pointer to the input vector to the requested
!  matrix-vector multiplication.
!
!  ipar(9) -- pointer to the output vector of the requested
!  matrix-vector multiplication.
!
!  To perform v = A * u, it is assumed that u is
!  w(ipar(8):ipar(8)+n-1) and v is stored as w(ipar(9):ipar(9)+n-1).
!
!  ipar(10) -- the return address (used to determine where to go
!  to inside the iterative solvers after the caller has performed
!  the requested services).
!
!  ipar(11) -- the result of the external convergence test
!  On final return from the iterative solvers, this value
!  will be reflected by ipar(1) = 0 (details discussed later)
!
!  ipar(12) -- error code of MGSRO, it is
!               1 if the input vector to MGSRO is linear combination
!                 of others,
!               0 if MGSRO was successful,
!              -1 if the input vector to MGSRO is zero,
!              -2 if the input vector contains invalid number.
!
!  ipar(13) -- number of initializations. During each initilization
!              residual norm is computed directly from M_l(b - A x).
!
!  ipar(14) to ipar(16) are NOT defined, they are NOT USED by
!  any iterative solver at this time.
!
!  Information about the error and tolerance are stored in the
!  array FPAR. So are some internal variables that need to be
!  saved from one iteration to the next one. Since the internal
!  variables are not the same for each routine, we only define the
!  common ones.
!
!  The first two are input parameters:
!  fpar(1) -- the relative tolerance,
!  fpar(2) -- the absolute tolerance (details discussed later),
!
!  When the iterative solver terminates,
!  fpar(3) -- initial residual/error norm,
!  fpar(4) -- target residual/error norm,
!  fpar(5) -- current residual norm (if available),
!  fpar(6) -- current residual/error norm,
!  fpar(7) -- convergence rate,
!
!  fpar(8:10) are used by some of the iterative solvers to save
!  some internal information.
!
!  fpar(11) -- number of floating-point operations. The iterative
!  solvers will add the number of FLOPS they used to this variable,
!  but they do NOT initialize it, nor add the number of FLOPS due
!  to matrix-vector multiplications (since matvec is outside of
!  the iterative solvers). To insure the correct FLOPS count, the
!  caller should set fpar(11) = 0 before invoking the iterative
!  solvers and account for the number of FLOPS from matrix-vector
!  multiplications and preconditioners.
!
!  fpar(12:16) are not used in current implementation.
!
!  Whether the content of fpar(3), fpar(4) and fpar(6) are residual
!  norms or error norms depends on ipar(3). If the requested
!  convergence test is based on the residual norm, they will be
!  residual norms. If the caller want to test convergence based the
!  error norms (estimated by the norm of the modifications applied
!  to the approximate solution), they will be error norms.
!  Convergence rate is defined by (Fortran 77 statement)
!  fpar(7) = log10(fpar(3) / fpar(6)) / (ipar(7)-ipar(13))
!  If fpar(7) = 0.5, it means that approximately every 2 (= 1/0.5)
!  steps the residual/error norm decrease by a factor of 10.
!
!  ............................................................
!  Stopping criteria,
!
!  An iterative solver may be terminated due to (1) satisfying
!  convergence test; (2) exceeding iteration limit; (3)
!  insufficient work space; (4) break-down. Checking of the work
!  space is only done in the initialization stage, i.e. when it is
!  called with ipar(1) == 0. A complete convergence test is done
!  after each update of the solutions. Other conditions are
!  monitored continuously.
!
!  With regard to the number of iteration, when ipar(6) is
!  positive, the current iteration number will be checked against
!  it. If current iteration number is greater the ipar(6) than the
!  solver will return with status -1. If ipar(6) is not positive,
!  the iteration will continue until convergence test is satisfied.
!
!  Two things may be used in the convergence tests, one is the
!  residual 2-norm, the other one is 2-norm of the change in the
!  approximate solution. The residual and the change in
!  approximate solution are from the preconditioned system (if
!  preconditioning is applied). The DQGMRES and TFQMR use two
!  estimates for the residual norms. The estimates are not
!  accurate, but they are acceptable in most of the cases.
!  Generally speaking, the error of the TFQMR's estimate is less
!  accurate.
!
!  The convergence test type is indicated by ipar(3). There are
!  four type convergence tests: (1) tests based on the residual
!  norm; (2) tests based on change in approximate solution; (3)
!  caller does not care, the solver choose one from above two on
!  its own; (4) caller will perform the test, the solver should
!  simply continue.
!  Here is the complete definition:
!   -2 == || dx(i) || <= rtol * || rhs || + atol
!   -1 == || dx(i) || <= rtol * || dx(1) || + atol
!    0 == solver will choose test 1 (next)
!    1 == || residual || <= rtol * || initial residual || + atol
!    2 == || residual || <= rtol * || rhs || + atol
!  999 == caller will perform the test
!  where dx(i) denote the change in the solution at the ith update.
!  ||.|| denotes 2-norm. rtol = fpar(1) and atol = fpar(2).
!
!  If the caller is to perform the convergence test, the outcome
!  should be stored in ipar(11).
!  ipar(11) = 0 -- failed the convergence test, iterative solver
!  should continue
!  ipar(11) = 1 -- satisfied convergence test, iterative solver
!  should perform the clean up job and stop.
!
!  Upon return with ipar(1) = 10,
!  ipar(8)  points to the starting position of the change in
!           solution Sx, where the actual solution of the step is
!           x_j = x_0 + M_r^{-1} Sx.
!           Exception: ipar(8) < 0, Sx = 0. It is mostly used by
!           GMRES and variants to indicate (1) Sx was not necessary,
!           (2) intermediate result of Sx is not computed.
!  ipar(9)  points to the starting position of a work vector that
!           can be used by the caller.
!
!  NOTE: the caller should allow the iterative solver to perform
!  clean up job after the external convergence test is satisfied,
!  since some of the iterative solvers do not directly update
!  the 'sol' array. A typical clean-up stage includes performing
!  the final update of the approximate solution and computing
!  the convergence information (e.g. values of fpar(3:7)).
!
!  NOTE: fpar(4) and fpar(6) are not set by the accelerators (the
!  routines implemented here) if ipar(3) = 999.
!
!  ..............................................................
!  Usage:
!
!  To start solving a linear system, the user needs to specify
!  first 6 elements of the ipar, and first 2 elements of fpar.
!  The user may optionally set fpar(11) = 0 if one wants to count
!  the number of floating-point operations. (Note: the iterative
!  solvers will only add the floating-point operations inside
!  themselves, the caller will have to add the FLOPS from the
!  matrix-vector multiplication routines and the preconditioning
!  routines in order to account for all the arithmetic operations.)
!
!  Here is an example:
!  ipar(1) = 0       ! always 0 to start an iterative solver
!  ipar(2) = 2       ! right preconditioning
!  ipar(3) = 1       ! use convergence test scheme 1
!  ipar(4) = 10000   ! the 'w' has 10,000 elements
!  ipar(5) = 10      ! use *GMRES(10) (e.g. FGMRES(10))
!  ipar(6) = 100     ! use at most 100 matvec's
!  fpar(1) = 1.0D-6  ! relative tolerance 1.0D-6
!  fpar(2) = 1.0D-10 ! absolute tolerance 1.0D-10
!  fpar(11) = 0.0    ! clearing the FLOPS counter
!
!  After the above specifications, one can start to call an
!  iterative solver, say BCG. Here is a piece of pseudo-code
!  showing how it can be done,
!
! 10   call bcg(n,rhs,sol,ipar,fpar,w)
!   if (ipar(1) == 1) then
!      call amux(n,w(ipar(8)),w(ipar(9)),a,ja,ia)
!      goto 10
!   else if (ipar(1) == 2) then
!      call atmux(n,w(ipar(8)),w(ipar(9)),a,ja,ia)
!      goto 10
!   else if (ipar(1) == 3) then
!      left preconditioner solver
!      goto 10
!   else if (ipar(1) == 4) then
!      left preconditioner transposed solve
!      goto 10
!   else if (ipar(1) == 5) then
!      right preconditioner solve
!      goto 10
!   else if (ipar(1) == 6) then
!      right preconditioner transposed solve
!      goto 10
!   else if (ipar(1) == 10) then
!      call my own stopping test routine
!      goto 10
!   else if (ipar(1) > 0) then
!      ipar(1) is an unspecified code
!   else
!      the iterative solver terminated with code = ipar(1)
!   endif
!
!  This segment of pseudo-code assumes the matrix is in CSR
!  format, AMUX and ATMUX are two routines from the SPARSKIT
!  MATVEC module. They perform matrix-vector multiplications for
!  CSR matrices, where w(ipar(8)) is the first element of the
!  input vectors to the two routines, and w(ipar(9)) is the first
!  element of the output vectors from them. For simplicity, we did
!  not show the name of the routine that performs the
!  preconditioning operations or the convergence tests.
!
!  This is a implementation of the Conjugate Gradient (CG) method
!  for solving linear system.
!
!  NOTE: This is not the PCG algorithm. It is a regular CG
!  algorithm. To be consistent with the other solvers, the
!  preconditioners are applied by performing 
!  Ml^{-1} A Mr^{-1} P in place of A P in the CG algorithm. The
!  PCG uses its preconditioners very differently.
!
!  fpar(7) is used here internally to store <r, r>.
!  w(:,1) -- residual vector
!  w(:,2) -- P, the conjugate direction
!  w(:,3) -- A P, matrix multiply the conjugate direction
!  w(:,4) -- temporary storage for results of preconditioning
!  w(:,5) -- change in the solution (sol) is stored here until
!            termination of this solver
!
!---------------------------------------------------------------
subroutine cg(n, kss0, rhs, sol, ipar, fpar, w, nnodes,comm)

  use constants
  implicit none
  !
  !  Input/Output variables:
  !
  integer, intent(in) :: n, kss0, nnodes, comm
  integer, intent(inout) :: ipar(16)
  real(dp), intent(in) :: rhs(n)
  real(dp), intent(inout) :: sol(n), fpar(16), w(n,kss0)
  !
  !  Work variables:
  !
  integer i
  real(dp) :: alpha
  logical lp,rp
  save
  !
  !  External functions:
  !
  real(dp), external :: distdot
  logical, external :: stopbis, brkdn
  !---------------------------------------------------------------
  !
  !  Check the status of the call
  !
  if (ipar(1) <= 0) ipar(10) = 0
  select case (ipar(10))
  case (1)
     goto 10
  case (2)
     goto 20
  case (3)
     goto 40
  case (4)
     goto 50
  case (5)
     goto 60
  case (6)
     goto 70
  case (7)
     goto 80
  end select
  !
  !  Initialization
  !
  call bisinit(ipar,fpar,5*n,1,kss0*n,lp,rp,w)
  if (ipar(1) < 0) return
  !
  !  Request for matrix vector multiplication A*x in the initialization
  !
  ipar(1) = 1
  ipar(8) = n+1
  ipar(9) = ipar(8) + n
  ipar(10) = 1
  do i = 1, n
     w(i,2) = sol(i)
  enddo
  return
10 ipar(7) = ipar(7) + 1
  ipar(13) = 1
  do i = 1, n
     w(i,2) = rhs(i) - w(i,3)
  enddo
  fpar(11) = fpar(11) + n
  !
  !  If left preconditioned
  !
  if (lp) then
     ipar(1) = 3
     ipar(9) = 1
     ipar(10) = 2
     return
  endif

20 if (lp) then
     do i = 1, n
        w(i,2) = w(i,1)
     enddo
  else
     do i = 1, n
        w(i,1) = w(i,2)
     enddo
  endif

  fpar(7) = distdot(n,w,1,w,1,nnodes,comm)
  fpar(11) = fpar(11) + 2 * n
  fpar(3) = sqrt(fpar(7))
  fpar(5) = fpar(3)
  if (ipar(3) == 2.or.ipar(3) == -2) then
     fpar(4) = fpar(1) * sqrt(distdot(n,rhs,1,rhs,1,nnodes,comm)) + fpar(2)
     fpar(11) = fpar(11) + 2 * n
  else if (ipar(3) /= 999) then
     fpar(4) = fpar(1) * fpar(3) + fpar(2)
  endif
  !
  !  Before iteration can continue, we need to compute A * p, which
  !  includes the preconditioning operations.
  !
30 if (rp) then
     ipar(1) = 5
     ipar(8) = n + 1
     if (lp) then
        ipar(9) = ipar(8) + n
     else
        ipar(9) = 3*n + 1
     endif
     ipar(10) = 3
     return
  endif

40 ipar(1) = 1
  if (rp) then
     ipar(8) = ipar(9)
  else
     ipar(8) = n + 1
  endif
  if (lp) then
     ipar(9) = 3*n+1
  else
     ipar(9) = n+n+1
  endif
  ipar(10) = 4
  return

50 if (lp) then
     ipar(1) = 3
     ipar(8) = ipar(9)
     ipar(9) = n+n+1
     ipar(10) = 5
     return
  endif
  !
  !  Continuing with the iterations.
  !
60 ipar(7) = ipar(7) + 1
  alpha = distdot(n,w(1,2),1,w(1,3),1,nnodes,comm)
  fpar(11) = fpar(11) + 2*n
  if (brkdn(alpha,ipar)) goto 90
  alpha = fpar(7) / alpha
  do i = 1, n
     w(i,5) = w(i,5) + alpha * w(i,2)
     w(i,1) = w(i,1) - alpha * w(i,3)
  enddo
  fpar(11) = fpar(11) + 4*n
  !
  !  Are we ready to terminate ?
  !
  if (ipar(3) == 999) then
     ipar(1) = 10
     ipar(8) = 4*n + 1
     ipar(9) = 3*n + 1
     ipar(10) = 6
     return
  endif
70 if (ipar(3) == 999) then
     if (ipar(11) == 1) goto 90
  else if (stopbis(n,ipar,1,fpar,w,w(1,2),alpha,nnodes,comm)) then
     goto 90
  endif
  !
  !  Continue the iterations.
  !
  alpha = fpar(5)*fpar(5) / fpar(7)
  fpar(7) = fpar(5)*fpar(5)
  do i = 1, n
     w(i,2) = w(i,1) + alpha * w(i,2)
  enddo
  fpar(11) = fpar(11) + 2*n
  goto 30
  !
  !  Clean up -- necessary to accommodate the right-preconditioning.
  !
90 if (rp) then
     if (ipar(1) < 0) ipar(12) = ipar(1)
     ipar(1) = 5
     ipar(8) = 4*n + 1
     ipar(9) = ipar(8) - n
     ipar(10) = 7
     return
  endif
80 if (rp) then
     call tidycg(n,ipar,fpar,sol,w(1,4))
  else
     call tidycg(n,ipar,fpar,sol,w(1,5))
  endif

end subroutine cg
!===============================================================
!
!  Inititialize parameters and work arrays
!
!---------------------------------------------------------------
subroutine bisinit(ipar,fpar,wksize,dsc,lwk,lp,rp,wk)

  use constants
  implicit none
  !
  !  Input/Output variables:
  !
  integer, intent(in) :: wksize, dsc, lwk
  integer, intent(inout) :: ipar(16)
  real(dp), intent(inout) :: fpar(16)
  logical, intent(out) :: lp,rp
  real(dp), intent(inout) :: wk(lwk)
  !
  !  Work variables:
  !
  integer i
  !---------------------------------------------------------------
  !
  !  ipar(1) = -2 inidcate that there are not enough space in the work
  !  array
  !
  if (ipar(4) < wksize) then
     ipar(1) = -2
     ipar(4) = wksize
     return
  endif

  if (ipar(2) > 2) then
     lp = .true.
     rp = .true.
  else if (ipar(2) == 2) then
     lp = .false.
     rp = .true.
  else if (ipar(2) == 1) then
     lp = .true.
     rp = .false.
  else
     lp = .false.
     rp = .false.
  endif
  if (ipar(3) == 0) ipar(3) = dsc
  !  .. clear the ipar elements used
  ipar(7) = 0
  ipar(8) = 0
  ipar(9) = 0
  ipar(10) = 0
  ipar(11) = 0
  ipar(12) = 0
  ipar(13) = 0
  !
  !  fpar(1) must be between (0, 1), fpar(2) must be positive,
  !  fpar(1) and fpar(2) can NOT both be zero
  !  Normally return ipar(1) = -4 to indicate any of above error
  !
  if (fpar(1) < zero .or. fpar(1) >= one .or. fpar(2) < zero .or. &
       (fpar(1) == zero .and. fpar(2) == zero)) then
     if (ipar(1) == 0) then
        ipar(1) = -4
        return
     else
        fpar(1) = 1.0D-6
        fpar(2) = 1.0D-16
     endif
  endif
  !  .. clear the fpar elements
  do i = 3, 10
     fpar(i) = zero
  enddo
  if (fpar(11) < zero) fpar(11) = zero
  !  .. clear the used portion of the work array to zero
  do i = 1, wksize
     wk(i) = zero
  enddo

end subroutine bisinit
!===============================================================
!
!  Test whether alpha is zero or an abnormal number, if yes,
!  this routine will return .true.
!
!  If alpha == 0, ipar(1) = -3,
!  if alpha is an abnormal number, ipar(1) = -9.
!
!---------------------------------------------------------------
logical function brkdn(alpha, ipar)

  use constants
  implicit none
  !
  !  Input/Output variables:
  !
  integer, intent(inout) :: ipar(16)
  real(dp), intent(in) :: alpha
  !
  !  Work variables:
  !
  real(dp) :: beta
  !---------------------------------------------------------------
  brkdn = .false.
  if (alpha > zero) then
     beta = one / alpha
     if (.not. beta > zero) then
        brkdn = .true.
        ipar(1) = -9
     endif
  else if (alpha < zero) then
     beta = one / alpha
     if (.not. beta < zero) then
        brkdn = .true.
        ipar(1) = -9
     endif
  else if (alpha == zero) then
     brkdn = .true.
     ipar(1) = -3
  else
     brkdn = .true.
     ipar(1) = -9
  endif

end function brkdn
!===============================================================
!
!  Function for determining the stopping criteria. return value of
!  true if the stopbis criteria is satisfied.
!
!  ---------------------------------------------------------------
logical function stopbis(n,ipar,mvpi,fpar,r,delx,sx,nnodes,comm)

  use constants
  implicit none
  !
  !  Input/Output variables:
  !
  integer, intent(in) :: n, mvpi, nnodes, comm
  integer, intent(inout) :: ipar(16)
  real(dp), intent(inout) :: fpar(16)
  real(dp), intent(in) :: r(n), delx(n), sx
  !
  !  External functions:
  !
  real(dp), external :: distdot
  !---------------------------------------------------------------
  if (ipar(11) == 1) then
     stopbis = .true.
  else
     stopbis = .false.
  endif
  if (ipar(6) > 0 .and. ipar(7) >= ipar(6)) then
     ipar(1) = -1
     stopbis = .true.
  endif
  if (stopbis) return
  !
  !  computes errors
  !
  fpar(5) = sqrt(distdot(n,r,1,r,1,nnodes, comm))
  fpar(11) = fpar(11) + 2 * n
  if (ipar(3) < 0) then
  !
  !  compute the change in the solution vector
  !
     fpar(6) = sx * sqrt(distdot(n,delx,1,delx,1,nnodes,comm))
     fpar(11) = fpar(11) + 2 * n
     if (ipar(7) < mvpi+mvpi+1) then
  !
  !  if this is the end of the first iteration, set fpar(3:4)
  !
        fpar(3) = fpar(6)
        if (ipar(3) == -1) then
           fpar(4) = fpar(1) * fpar(3) + fpar(2)
        endif
     endif
  else
     fpar(6) = fpar(5)
  endif
  !
  !  .. the test is struct this way so that when the value in
  !  fpar(6) is not a valid number, STOPBIS is set to .true.
  !
  if (fpar(6) > fpar(4)) then
     stopbis = .false.
     ipar(11) = 0
  else
     stopbis = .true.
     ipar(11) = 1
  endif

end function stopbis
!===============================================================
!
!  Some common operations required before terminating the CG
!  routines.
!
!---------------------------------------------------------------
subroutine tidycg(n,ipar,fpar,sol,delx)

  use constants
  implicit none
  !
  !  Input/Output variables:
  !
  integer, intent(in) :: n
  integer, intent(inout) :: ipar(16)
  real(dp), intent(inout) :: fpar(16)
  real(dp), intent(in) :: delx(n)
  real(dp), intent(inout) :: sol(n)
  !
  !  Work variables:
  !
  integer i
  !---------------------------------------------------------------
  if (ipar(12) /= 0) then
     ipar(1) = -3
  else if (ipar(1) > 0) then
     if ((ipar(3) == 999 .and. ipar(11) == 1) .or. fpar(6) <= fpar(4)) then
        ipar(1) = 0
     else if (ipar(7) >= ipar(6) .and. ipar(6) > 0) then
        ipar(1) = -1
     else
        ipar(1) = -10
     endif
  endif
  if (fpar(3) > zero .and. fpar(6) > zero .and. ipar(7) > ipar(13)) then
     fpar(7) = log10(fpar(3) / fpar(6)) / real(ipar(7)-ipar(13),dp)
  else
     fpar(7) = zero
  endif
  do i = 1, n
     sol(i) = sol(i) + delx(i)
  enddo

end subroutine tidycg
!===============================================================
