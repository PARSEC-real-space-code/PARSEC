!===============================================================
!
!  Copyright (C) 2005 Finite Difference Research Group
!  This file is part of parsec, http://www.ices.utexas.edu/parsec/
!
!  Define "reciprocal space" information and constructs the
!  reciprocal space.
!  AJB:
!  Call fourier_f to construct the relevant quanitites
!  (this breaks fourier filtering of pseudos for non pbc calculations)
!---------------------------------------------------------------
subroutine inipbc(clust,p_pot,grid,pbc,symm,ipr,ierr)

  use constants
  use cluster_module
  use pseudo_potential_module
  use grid_module
  use pbc_module
  use symmetry_module

  implicit none
  !
  !  Input/Output variables:
  !
  ! the cluster - needed for fourier_f
  type (cluster), intent(in) :: clust
  ! pseudo_potential structure - needed for fourier_f
  type (pseudo_potential), intent(inout) :: p_pot
  type (grid_data), intent(inout) :: grid
  type (pbc_data), intent(inout) :: pbc 
  type (symmetry), intent(in) :: symm
  !  print flag
  integer, intent(in) :: ipr
  !  error flag
  integer, intent(out) :: ierr
  !
  !  Work variables:
  !
  integer ii, kmax(3)
  real(dp) :: htemp
  real(dp) :: tmpmat(3,3), tmptr, tmpdet, tmpnorm

  !---------------------------------------------------------------

  ierr = 0

  write(7,*) 'Reciprocal Space Data:'
  write(7,*) '----------------------'

  tmpmat = transpose(pbc%latt_vec)
  pbc%adot = matmul(tmpmat, pbc%latt_vec)

  !  Using the following equation for the reciprocal
  !  lattice vectors:     
  !
  !            T                 -1
  !  [b1 b2 b3]  =  2pi [a1 a2 a3]  
  !
  !  note the definition of bvec here means that
  !  reciprocal lattice vector i is given by bvec(:,i) !
  !  this exactly as latt_vec arrangement. 
  !  while lattic vector i is given by avec(i,:).

  pbc%bvec = transpose(pbc%latt_vec)

  call mtrxin(pbc%bvec, tmpdet, tmptr)

  !  Now multiplying bvec by 2Pi to agree to definition of 
  !  reciprocal lattice vector.

  pbc%bvec = twopi * pbc%bvec

  !  Calculating the metric matrix.

  pbc%bdot = matmul(transpose(pbc%bvec), pbc%bvec)

  write(7,*) 'Reciprocal Lattice (a.u.):'

  do ii = 1, pbc%per
     write(7,'(3(f15.10, 1x))') pbc%bvec(1:pbc%per,ii)
  enddo
  
  !AJB: this is not needed because adjust_step will do this in a moment
  grid%step = grid%stepin 


  !AJB: anyways, these all seem to me as simply inits, since they are being
  !changed. Also idint is an old int function.
  grid%n1 = idint(two*grid%rmax/grid%stepin) + 2
  grid%n2 = grid%n1
  grid%n3 = grid%n1

  pbc%n1 = 1
  pbc%n2 = 1
  pbc%n3 = 1

  !  The new grid spacings must be commensurate with the lattice vectors as
  !  well as with all fractional translations present in non-symmorphic
  !  symmetry operations.
  call adjust_step(symm,pbc%per,pbc%box_size,grid%stepin,grid%step,ierr)
  !AJB: did anyone even think to change also hcub or hcub2? no?
  htemp = grid%stepin
  grid%stepin=minval(grid%step(1:pbc%per))

  if (ierr > 0) return

  if ( pbc%per >= 1 ) then
      !nint rounds to nearest whole number
     grid%n1 = 2*nint(pbc%box_size(1)/(two*grid%step(1)))
     pbc%n1 = grid%n1
  endif
  if ( pbc%per >= 2 ) then
      !nint rounds to nearest whole number
     grid%n2 = 2*nint(pbc%box_size(2)/(two*grid%step(2)))
     pbc%n2 = grid%n2
  endif
  if ( pbc%per >= 3 ) then
      !nint rounds to nearest whole number
     grid%n3 = 2*nint(pbc%box_size(3)/(two*grid%step(3)))
     pbc%n3 = grid%n3
  endif


  !  Calculation of inverse of normalized lattice vectors
  !  matrix - this is needed for the calculation in gradmvs.

  tmpmat = pbc%latt_vec

  do ii = 1, 3
     tmpnorm = sqrt(sum(tmpmat(:,ii)**2))
     tmpmat(:,ii) = tmpmat(:,ii)/tmpnorm
  enddo

  pbc%avec_norm = tmpmat

  call mtrxin(tmpmat,tmpdet,tmptr)

  grid%grad_bvec_norm = transpose(tmpmat)
      

  !  Calculation of volume element is done like this:
  ! AJB: umm, what?
  !  hcub=det(grid%step(1)*x_vec+grid%step(2)*y_vec*grid%step(2)*z_vec)
  !  x_vec,y_vec and z_vec are the unit vectors in direction of lattice
  !  vectors.

  tmpmat = pbc%latt_vec
      
  do ii = 1, 3
     tmpmat(:,ii) = grid%step(ii)*tmpmat(:,ii)/sqrt(sum(tmpmat(:,ii)**2))
  enddo

  !  Calling inversion routine because it computes also the 
  !  required determinant and puts it in grid%hcub.
  !  AJB: assuming the program wont crash because the matrix is singular
  call mtrxin(tmpmat,grid%hcub,tmptr)

  !  So since we changed grid%step we also change hcub
  grid%hcub=abs(grid%hcub)
  ! ...and hcub2
  grid%hcub2 = two/grid%hcub
  !
  !  kmax = n/2 if n is even; kmax = (n+1)/2 otherwise
  !
  !  kmax = 0 along non-periodic directions
  kmax = 0
  if (pbc%per >= 1) kmax(1) = nint( (real(grid%n1,dp) + one/four)/ two )
  if (pbc%per >= 2) kmax(2) = nint( (real(grid%n2,dp) + one/four)/ two )
  if (pbc%per >= 3) kmax(3) = nint( (real(grid%n3,dp) + one/four)/ two )

  pbc%vcell = grid%n1*grid%n2*grid%n3*grid%hcub


  ! O. SINAI: Calculate 2-D cell area directly from lattice vectors
  if(pbc%per==2) then
     pbc%a_surface_cell = abs(pbc%latt_vec(1,1)*pbc%latt_vec(2,2) - pbc%latt_vec(1,2)*pbc%latt_vec(2,1))
  endif
 
  !moved fourier_f from parsec.f90 to here, since the step is being adjusted!

  ! ===============================================================
  ! Calculate Fourier transform of pseudopotentials
  ! Performed by master PE only.
  ! ===============================================================
  ! Actually, this is a fourier filtering and stuff is being replaced.
  ! However, this is only done if alpha and acore are set to all atom types
  ! More importantly, this also defines the Fourier transforms of the local
  ! pseudopotential and core-charge-density, so maybe it is a good idea to split
  ! this into two functions!
  ! It is called here because the step size is important, but we can't 
  ! hide the fact that the fourier transforms in fourier_f are 1D and don't
  ! care about grid%step being a vector.

  ! Also,
  ! Because this is called inside inipbc, fourier filtering is not available
  ! to non pbc calculations...
      call fourier_f(clust,p_pot,pbc,grid%step(1),ipr) !instead of grid%stepin



  !  rescale alpha energy
  ! (meaning that pbc%ealpha has to be set beforehand)
  ! this can easily be moved into fourier_f 
  pbc%ealpha = pbc%ealpha/pbc%vcell

  if(pbc%per==2) then
      !AJB: why only if per==2? what's going on with per==1?
    pbc%ealpha=0
  endif
  !  construct reciprocal space
  call g_space(pbc,symm,kmax,ipr,ierr)

  write(7,*)
  write(7,21) htemp
  write(7,22) grid%step
  write(7,24) grid%n1,grid%n2,grid%n3
  write(7,*)

21 format(' Input grid spacing:', 1x, f6.3, 1x, 'bohr')
22 format(' Final grid spacings:', 1x, 3(f6.3, 1x), 'bohr')
24 format(' Number of grid points along each direction:',3(1x,i3))

end subroutine inipbc
!===============================================================
!
!  For input coordinates (i,j,k) in FFT grid, calculate the address
!  and phase associated to these coordinates. Input coordinates are
!  of the form 0 < i =< n1, and similarly for j, k. Output phase is
!  exp( -i G . V ) where V is the shift vector: coordinates of the
!  corner point in the real-space grid.
!
!---------------------------------------------------------------
subroutine get_address(pbc,i,j,k,iadd,phase)
  use constants
  use pbc_module
  implicit none
  !
  !  Input/Output variables:
  !
  !  pbc related data
  type (pbc_data), intent(inout) :: pbc
  !  coordinates of g-vector in FFT grid
  integer, intent(in) :: i, j, k
  !  address of g-vector in FFT mesh
  integer, intent(out) :: iadd
  !  phase associated to this g-vector
  complex(dpc), intent(out) :: phase
  !
  !  Work variables:
  !
  !  exponent of phase
  real(dp) :: phi
  !  coordinates of g-vector
  integer :: gv(3)
  !---------------------------------------------------------------
  gv(1) = i - 1 - pbc%n1
  gv(2) = j - 1 - pbc%n2
  gv(3) = k - 1 - pbc%n3
  if (gv(1) < pbc%mx) gv(1) = gv(1) + pbc%n1
  if (gv(2) < pbc%my) gv(2) = gv(2) + pbc%n2
  if (gv(3) < pbc%mz) gv(3) = gv(3) + pbc%n3
  iadd = ( (k-1)*pbc%n2 + j-1 )*pbc%n1 + i
  phi = dot_product(gv,pbc%shift)
  !phase = dcmplx(cos(phi),sin(phi))
  phase = cmplx(cos(phi), sin(phi), kind=dp)

end subroutine get_address
!===============================================================
!
!  Performs matrix-vector multiplications in 3-d arrays. This
!  should be used in place of matmul (which assumes both arguments
!  to be matrices). Input is 3x3 matrix M and vector V. Output is
!  vector MV.
!  If op='T', performs MV = transpose(M)*V
!  If op='N', performs MV = M*V
!
!---------------------------------------------------------------
subroutine matvec3(op,mat,vec,mvec)

  use constants
  implicit none
  !
  !  Input/Output variables:
  !
  !  operation: transpose(M) or M
  character (len=1), intent(in) :: op
  !  input matrix
  real(dp), intent(in) :: mat(3,3)
  !  input vector V
  real(dp), intent(in) :: vec(3)
  !  output vector MV = op(M)*V
  real(dp), intent(out) :: mvec(3)
  !
  !  Work variables:
  !
  integer :: ii
  real(dp) :: vtmp(3)
  !---------------------------------------------------------------
  do ii = 1, 3
     if (op == 'N') then
        vtmp(ii) = dot_product(mat(ii,:),vec)
     elseif (op == 'T') then
        vtmp(ii) = dot_product(mat(:,ii),vec)
     endif
  enddo
  mvec = vtmp

  return
end subroutine matvec3
!===============================================================
!
!  Adjusts the grid spacing according to unit cell length and
!  translation vectors. Along each direction, the grid spacing (h)
!  must satisfy the following commensurability requirements:
!    a) a/h must be an even number, where "a" is the length of
!       the unit lattice vector along that direction. This way,
!       the real-space grid will be centered on the origin.
!    b) tau/h must be integer, where "tau" is any non-zero
!       translation vector (stored in symm%tnp). This requirement
!       is relevant only in non-symmorphic symmetry operations.
!       It makes the easy symmetrization of charge density
!       coded in subroutine scalar_symm possible.
!
!  author: Murilo L Tiago, ORNL, July 2007.
!
!---------------------------------------------------------------
subroutine adjust_step(symm,nper,box_size,stepin,stepout,ierr)

  use constants
  use symmetry_module

  implicit none
  !
  !  Input/Output variables:
  !
  type (symmetry), intent(in) :: symm

  !  number of periodic directions
  integer, intent(in) :: nper
  !  length of lattice vectors, input grid spacing
  real(dp), intent(in) :: box_size(3), stepin
  !  output grid spacing
  real(dp), intent(out) :: stepout(3)
  ! error flag, 265 < ierr < 271
  integer, intent(out) :: ierr
  !
  !  Work variables:
  !
  !  counters
  integer :: itrans, idir, nn
  !  temporary arrays
  real(dp) :: xtest, ytest, tmpvec(3), tmpmat(3,3)
  !  maximum number of grid order trials
  integer, parameter :: nmax = 1000
  !  tolerance in displacements
  real(dp), parameter :: tol = 1.d-8
  !  number of non-zero translations for each periodic direction
  integer :: ntau(nper)
  !  non-zero translations (in units of lattice vectors)
  real(dp) :: tau(nmax,nper)
  !  shortest translation
  real(dp) :: tau_min
  !  true if the current grid spacing is commensurate with all translations
  logical :: lcomm

  !---------------------------------------------------------------

  ntau = 0
  tau = zero

  !  Go through all symmetry operations and search for the ones which have
  !  fractional translations. Store those translations in array tau. When a
  !  new translation is found, make it positive and increment the counter ntau.
  do itrans = 1, symm%ntrans
     tmpmat = real(symm%rmtrx(:,:,itrans),dp)
     call matvec3('N',tmpmat,symm%tnp(:,itrans),tmpvec)
     do idir = 1, nper
        if (tmpvec(idir) < zero) tmpvec(idir) = -tmpvec(idir)
        if (tmpvec(idir) > tol) then
           ntau(idir) = ntau(idir) + 1
           tau(ntau(idir),idir) = tmpvec(idir)
        endif
     enddo
  enddo

  do idir = 1, nper
     !  Search for shortest translation. Notice that tau has zeros beyond
     !  ntau. If all translations are zero, choose tau_min = one.
     if (ntau(idir) == 0) then
        tau_min = one
     else
        tau_min = minval(tau(1:ntau(idir),idir))
     endif

     do nn = 1, nmax
        lcomm = .true.
        if (tau_min/real(nn,dp) < stepin/box_size(idir)) exit
        xtest = tau_min/real(nn,dp)
        !  Is this grid spacing commensurate with all translations?
        do itrans = 1, ntau(idir)
           ytest = tau(itrans,idir)/xtest - nint(tau(itrans,idir)/xtest)
           if ( abs(ytest) > tol) then
              lcomm = .false.
              cycle
           endif
        enddo
        !  Is the ratio between this grid spacing and lattice vector an
        !  even number?
        ytest = half/xtest - nint(half/xtest)
        if ( abs(ytest) > tol) then
           lcomm = .false.
           cycle
        endif
        !  Save this grid spacing if commensurability conditions are satisfied.
        if (lcomm) then
           stepout(idir) = xtest * box_size(idir)
        endif
     enddo
     !  Save this grid spacing if commensurability conditions are satisfied.
     if (nn == nmax) then
        write(9,*) 'ERROR in adjust_step! Could not find a commensurate ', &
             'grid spacing.'
        ierr = 266
        return
     endif
  enddo

  return
end subroutine adjust_step
!===============================================================
