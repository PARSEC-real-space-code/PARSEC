!===============================================================
!
! Copyright (C) 2005 Finite Difference Research Group
! This file is part of parsec, http://www.ices.utexas.edu/parsec/
!
!              Anderson mixing scheme
!
! This suborutine mixes the old and new potential to create a new
! guess for the potential in the self-consistent cycle.
! D. G. Anderson, J. Assoc. Computing Machinery, Vol12, p547(1965).
! V. Eyert, J. Comput. Phys. Vol124,P271(1996)
!
! NOTE: the mixed potential is returned as xin, BUT to be
! consistent with the old anderson.f and parsec, xout is also
! updated as the new potential.
!
! Author: Lingzhu Kong(2005)
!    
!---------------------------------------------------------------
subroutine anderson_mix(mixer,parallel,iter,nrep,xin,xout,ierr)

  use constants
  use mixer_module
  use parallel_data_module
  implicit none
  !
  ! Input/Output variables:
  !
  ! mixer related data
  type (mixer_data), intent(inout) :: mixer
  ! parallel computation related data
  type (parallel_data), intent(in) :: parallel

  ! iteration number
  integer, intent(in) :: iter
  ! order of the reduced group
  integer,intent(in) :: nrep
  ! error flag, 440 < ierr < 451
  integer,intent(out) :: ierr

  ! old potential
  real(dp), intent(inout) :: xin(mixer%ndim * mixer%nspin)

  ! new potentials 
  real(dp), intent(inout) :: xout(mixer%ndim * mixer%nspin)
  !
  ! Work variables:
  !
  ! Residual vectors F[n] = ( xout - xin ) (Eq(2.1) in Eyert)
  real(dp), allocatable :: f(:,:)

  ! Mixed input (Eq(4.1) in Eyert)
  real(dp), allocatable :: aver(:)

  ! Mixed residual(Eq(4.2) in Eyert)
  real(dp), allocatable :: faver(:)

  ! Coefficient matrix in (4.3)(Eyert)
  real(dp)   b(mixer%memory,mixer%memory)

  ! Right hand side(RHS) vector in (4.3)(Eyert)
  real(dp)   b0(mixer%memory,1)
  ! 
  ! dimension of the potential vector
  integer   dumb
  ! temp argument for psum
  real(dp), dimension(1) :: xsumvec
  ! multiplicity of grid points
  real(dp)  weight
  ! counters, temporary variables
  integer    i,j,k,anderi,ipiv(mixer%memory),info,alcstat
  real(dp) :: xsum

  !---------------------------------------------------------------
  !
  dumb = mixer%ndim * mixer%nspin

  allocate(f (dumb,mixer%memory+1),stat=alcstat)
  call alccheck('f ',dumb * (mixer%memory+1), alcstat)

  allocate(aver (dumb),stat=alcstat)
  call alccheck('aver ',dumb , alcstat)

  allocate(faver (dumb),stat=alcstat)
  call alccheck('faver ',dumb , alcstat)

  anderi = mod(iter-1,mixer%restart) + 1
  anderi = min(anderi,mixer%memory+1)

  weight = real(nrep,dp)

  if(anderi == 1)then
     ! first iteration - previous information not available yet!
     do k = 1,dumb

        ! store for next iteration
        mixer%xinold(k,1) = xin(k)
        mixer%xoutold(k,1) = xout(k)

        ! next guess obtained by pure linear mixing
        xin(k) = mixer%param* xout(k) + (one - mixer%param) * xin(k)
        xout(k) = xin(k)
     enddo

  else

     ! residual from current iteration 
     do k = 1,dumb
        f(k,1) = xout(k) - xin(k)
     enddo

     ! residual from previous iterations
     do j = 2,anderi
        do k = 1,dumb
           f(k,j) = mixer%xoutold(k,j-1) - mixer%xinold(k,j-1)
        enddo
     enddo

     ! Coefficient matirx b and RHS vector b0 initialised as zero
     b(:,:) = zero
     b0(:,:) = zero

     ! update matrix b from redidual vectors f
     do j = 1,anderi-1
        do i = 1,j
           xsum = zero
           do k = 1,dumb
              xsum = xsum+(f(k,1)-f(k,i+1)) * (f(k,1)-f(k,j+1))
           enddo
           xsumvec = xsum
           call psum(xsumvec,1,parallel%group_size,parallel%group_comm)
           xsum = xsumvec(1)
           b(i,j) = xsum*weight

           if(i /= j)b(j,i) = b(i,j)
        enddo
     enddo
     ! End of the update of matrix b

     ! update vectors b0
     do i = 1,anderi-1
        xsum = zero
        do k = 1,dumb
           xsum = xsum+(f(k,1)-f(k,i+1)) * f(k,1)
        enddo
        xsumvec = xsum
        call psum(xsumvec,1,parallel%group_size,parallel%group_comm)
        xsum = xsumvec(1)
        b0(i,1) = xsum * weight
     enddo
     ! End  of the update of b0

     ! Solve b * X = b0 for X, overwrite b0. vector X is the
     ! coefficients of the previous iterations that enter the new
     ! guess potential.
     info = 0
     call dgesv( anderi-1,1,b,mixer%memory,ipiv,b0,mixer%memory, info )
 
     if( info /= 0)then
        write(9,*) ' ERROR: dgesv diagonalization failed in '
        write(9,*) ' subroutine anderson.f, info = ',info
        ierr = 441
        return
     endif

     ! mixed input and mixed residual 
     do k = 1,dumb
        aver(k) = xin(k)
        faver(k) = f(k,1)
        do i = 1,anderi-1
           aver(k) = aver(k) + b0(i,1)*( mixer%xinold(k,i)-xin(k) )
           faver(k) = faver(k) + b0(i,1)*( f(k,i+1) - f(k,1) )
        enddo
     enddo

     ! Store the old information for next mxing
     if( anderi == mixer%memory+1 )then
        if (anderi > 2)then
           do k = 1,dumb
              do i = anderi-2,1,-1
                 mixer%xinold(k,i+1) = mixer%xinold(k,i)
                 mixer%xoutold(k,i+1) = mixer%xoutold(k,i)
              enddo
           enddo
        endif

     else
        do k = 1,dumb
           do i = anderi-1,1,-1
              mixer%xinold(k,i+1) = mixer%xinold(k,i)
              mixer%xoutold(k,i+1) = mixer%xoutold(k,i)
           enddo
        enddo
     endif
     ! end of store

     ! Linear mix the mixed input and mixed residual for the new potential.
     do k = 1,dumb
        mixer%xinold(k,1) = xin(k)
        mixer%xoutold(k,1) = xout(k)
        xin(k) = aver(k) + mixer%param*faver(k)
        xout(k) = xin(k)
     enddo

  endif

  deallocate( f , aver , faver)

end subroutine anderson_mix
!===============================================================
