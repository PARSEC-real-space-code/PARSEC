!===============================================================
!
! Copyright (C) 2005 Finite Difference Research Group
! This file is part of parsec, http://www.ices.utexas.edu/parsec/
!
! Authors: M. Jain, L. Kronik
!
! This subroutine mixes the old and new potential according to
! the GENERALIZED Broyden method, using p memory steps. 
!
! NOTE: This is NOT the modified Broyden method of D. D. Johnson.
! This is an original, unpublished algorithm developed by M. Jain
! and L. Kronik. The gist of the algorithm is similar to the
! generalization of the Broyden suggested by Eyert (J. Comp. Phys.
! 124, 271, 1996), but it is based on considerations for updating
! the inverse Jacobian with respect to that of the previous
! iteration, not with respect to that of the first iteration. 
! Eyert's paper also explains why this approach is inherently
! superior to the modified Broyden algorithm.
!
! The equations implemented are:
!
! (a) B_n = B_1 + \Sum_{i=2}^n (u_i*v_i^t)
! (b) f_i=xout_i-xin_i ;
!     df_i = f_i-f_(i-1);
!     dxin_i = xin_i-xin_(i-1)
! (c) u_i = B_1*df(i)+dxin_i-\Sum_{k=2}^{i-1} [(u_k*v_k^t*df_i]
! (d) T_{i,j} = df_i^t*df_j
! (e) v_i^t = \Sum_{k=1}^p {[T^(-1)]_pk * df_(i-p+k)^t
! (f) x_new = x_in-B_n*f_n
!
! where:
!         B is the inverse Jacobian,
!         1 denotes thec first iteration,
!         n denotes the nth iteration,
!         t indicates transposition,
!         T is the "mixing matrix" of dimensions pxp.
!
! By combining Eqs. (a) and (f), there is never a need to store
! B_n!! We only need to store the history of the u, v
! intermediate vectors, and this is what we do. This trick is
! known as the Srivastava Scheme (G. P. Srivastava, J. Phys. A:
! Math. Gen. 17, L317, 1984). Even this storage becomes difficult
! for many iterations. This subroutine is therefore an
! "out-of-core" one (i.e., it communicates with the hard disk) in
! order to avoid prohibitive memory requirements.
!
!---------------------------------------------------------------
subroutine genbro(mixer,parallel,iter,nrep,xin,xout,ierr)

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
  integer, intent(in) :: nrep
  ! input potential
  real(dp), intent(in) :: xin(mixer%ndim * mixer%nspin)
  ! output potential on input, mixed potential on output
  real(dp), intent(inout) :: xout(mixer%ndim * mixer%nspin)
  ! error flag, 450 < ierr < 461
  integer, intent(out) :: ierr
  !
  ! Work variables:
  !
  ! Broyden iteration number (zeroed after each restart)
  integer :: bitr
  ! suffix for temporary files
  character(len=4) :: idstring
  ! multiplicity of grid points in irreducible wedge (equal to
  ! order of the reduced group)
  real(dp) :: weight
  ! counters
  integer i,j,m,nn
  ! temporary holder of resid_in
  real(dp) :: dvin
  ! accumulators for scalar multiplications
  real(dp) :: vtemp, vfac
  ! error flag for matrix inversion calls
  integer info 
  ! number of read blocks and size of each read block from
  ! bro.* files
  integer nblk, rsize
  ! string with mixing history file name
  character(len=6) fnam
  ! initial Jacobian
  real(dp) :: pmix
  ! size parameter for memory block allocated for the u,v update
  ! arrays
  integer blk
  ! allocation check
  integer alcstat
  ! system dimension
  integer dim
  ! Jacobian update history vectors, work arrays
  real(dp), dimension (:,:), allocatable :: u,v
  ! history of resid_out-resid_in, where resid_out is out-out_old, 
  ! resid_in is in-in_old
  real(dp), dimension (:,:), allocatable :: df
  ! work arrays for matrix inversion calls
  integer, dimension(mixer%memory) :: ipiv
  real(dp), dimension(mixer%memory) :: work
  !
  ! External functions:
  !
  real(dp), external :: distdot

  !---------------------------------------------------------------
  !
  ! initialize parameters and work arrays
  !
  pmix = mixer%param
  weight = real(nrep,dp)
  dim = mixer%ndim * mixer%nspin
  blk = mixer%block
  write(idstring,'(I4.4)') parallel%iam

  allocate(u(dim,2:mixer%block+1),stat=alcstat)
  call alccheck('genbro_u',dim*mixer%block,alcstat)
  allocate(v(dim,2:mixer%block+1),stat=alcstat)
  call alccheck('genbro_v',dim*mixer%block,alcstat)
  allocate(df(dim , 1:mixer%memory),stat=alcstat)
  call alccheck('genbro_df',dim*mixer%memory,alcstat)

  u(:,:) = zero
  v(:,:) = zero
  df(:,:) = zero
  ipiv(:) = 0
  work(:) = zero
  bitr = mod(iter-1,mixer%restart) + 1
  !
  ! if p=0 the algorithm reduces to trivial linear mixing
  !
  if (mixer%memory == 0) then
     xout = (1-pmix)*xin + pmix*xout
     goto 99
  endif

  ! if p is not 0, there's "real work" to do:
 
  ! The first iteration must be treated separately because various
  ! things need to be initialized.
 
  if (bitr == 1) then
     ! initialize inverse Jacobian matrix. Here it is simply a scalar,
     ! pmix, but has been given an array because a more sophisticated
     ! first guess may be desired in the future.

     ! Also, save xin, xout as xinold, xoutold, for the next iteration

     mixer%resid1 = -pmix
     mixer%old1 = xin
     mixer%old2 = xout

     ! Perform a Broyden update based on the initial guess for the 
     ! inverse Jacobian [Eq. (f) above]
 
     do i = 1, dim
        xout(i) = xin(i)-mixer%resid1(i)*(xout(i) - xin(i))
     enddo
     mixer%t(:,:) = zero
     mixer%tinv(:,:) = zero
     goto 99
  endif
 
  ! If this is not the first iteration, we must distinguish whether
  ! we are filling up the "memory buffer" of previous iterations
  ! until it occupies all p previous memory steps (bitr <= p+1)
  ! and the case where we are updating it on a FIFO
  ! basis (bitr > p+1).
 
  m = bitr - 1

  if (bitr <= mixer%memory+1)  then

     ! read previous values of u,v vectors from bro.1
     ! read previous values of df from bro.df

     open(10,file='bro.1.'//idstring,form='unformatted')
     open(11,file='bro.df.'//idstring,form='unformatted')
 
     do j = 2, m
        read(10) (u(i,j),i=1,dim)
     enddo
     do j = 2, m
        read(10) (v(i,j),i=1,dim)
     enddo
     do j = 1, m-1
        read(11) (df(i,j),i=1,dim)
     enddo
     close(10)
     close(11)
 
     ! Compute current df, u, v [Eqs. (c),(d)]
 
     m = bitr
 
     do i = 1, dim
        dvin = xin(i) - mixer%old1(i)
        df(i,m-1) = xout(i) - mixer%old2(i) - dvin
        u(i,m) = dvin - mixer%resid1(i)*df(i,m-1)
        v(i,m) = zero
     enddo

     do j = 2, m-1
        vfac = distdot(dim,v(1,j),1,df(1,m-1),1, &
             parallel%group_size,parallel%group_comm)*weight
        do i = 1, dim
           u(i,m) = u(i,m) - vfac*u(i,j)
        enddo
     enddo
 
     ! Update the T matrix (Eq. (d)]
     ! First the row and then the column
 
     do j= 2, m
        vtemp = distdot(dim,df(1,m-1),1,df(1,j-1),1, &
             parallel%group_size,parallel%group_comm)*weight
        mixer%t(m-1,j-1) = vtemp
        mixer%t(j-1,m-1) = vtemp
     enddo

     ! Invert the T matrix
     do j = 1, m-1
        do i = 1, m-1
           mixer%tinv(j,i) = mixer%t(j,i)
        enddo
     enddo
     do j = 1, m-1
        ipiv(j) = 0
     enddo

     call dgetrf(m-1,m-1,mixer%tinv,mixer%memory,ipiv,info)

     if (info < 0) then
        write(9,*) 'ERROR: problem with mixer matrix inversion'
        write(9,*) 'STOP in genbro, info = ',info
        ierr = 451
        call myflush(9)
        return
     endif
     if (info > 0) then
        write(9,*) 'Warning: T matrix in mixer is Singular'
     endif

     call dgetri(m-1,mixer%tinv,mixer%memory,ipiv,work,m-1,info)
     if (info < 0) then
        write(9,*) 'ERROR: problem with mixer matrix inversion'
        write(9,*) 'STOP in genbro, info = ',info
        ierr = 452
        call myflush(9)
        return
     endif
     if (info > 0) then
        write(9,*) 'Warning: T matrix in mixer is Singular'
     endif

     ! update v according to the T matrix [Eq. (e) above]
     do i = 1, dim
        do j = 2,m
           v(i,m) = v(i,m) + mixer%tinv(m-1,j-1)*df(i,j-1)
        enddo
     enddo

     ! save xin, xout as xinold, xoutold, for the next iteration

     mixer%old1 = xin
     mixer%old2 = xout

     ! perform first part of computation of x_new (put in xout). This
     ! is the "B_1*f_n" fragment of Eqs. (a)+(f).

     do i = 1, dim
        xout(i) =  mixer%resid1(i)*(xout(i) - xin(i))
     enddo

     ! update x_out according to the update of B_1 to B_n
     do j = 2, m
        vfac = (distdot(dim,v(1,j),1,mixer%old2,1, &
             parallel%group_size,parallel%group_comm) - &
             distdot(dim,v(1,j),1,xin,1, &
             parallel%group_size,parallel%group_comm) )*weight
        do i = 1, dim
           xout(i) = xout(i) + vfac*u(i,j)
        enddo
     enddo

     ! complete the update of x_out

     xout = xin - xout

     ! write u,v,df out to disk for the next iteration.
     open(10,file='bro.1.'//idstring,form='unformatted')
     open(11,file='bro.df.'//idstring,form='unformatted')
     rewind (10)
     rewind (11)
     do j = 2, m
        write(10) (u(i,j),i=1,dim)
     enddo
     do j = 2, m
        write(10) (v(i,j),i=1,dim)
     enddo
     if (m /= mixer%memory+1) then
        do j = 1, m-1
           write(11) (df(i,j),i=1,dim)
        enddo
     else 
        do j = 2, m-1
           write(11) (df(i,j),i=1,dim)
        enddo
     endif
     close(10)
     close(11)
 
  else
 
     ! read the df arrays
     open(11,file='bro.df.'//idstring,form='unformatted')
     do j = 1, mixer%memory-1
        read(11) (df(i,j),i=1,dim)
     enddo
     close(11)

     ! Update u,v,df for inverse Jacobian matrix
     do i = 1, dim
        dvin = xin(i) - mixer%old1(i)
        df(i,mixer%memory) = xout(i) - mixer%old2(i) - dvin
        u(i,blk+1) = dvin - mixer%resid1(i)*df(i,mixer%memory)
        v(i,blk+1) = zero
     enddo
     ! update xin, xout to xinold, xoutold for next iteration
     mixer%old1 = xin
     mixer%old2 = xout

     ! perform the B_1 fragment of Eq. (f) (see similar comment above)
     do i = 1, dim
        xout(i) =  mixer%resid1(i)*(xout(i) - xin(i))
     enddo

     ! Shift the T matrix to make room for the new information
 
     do i = 1, mixer%memory-1
        do j = 1,mixer%memory-1
           mixer%t(i,j) = mixer%t(i+1,j+1)
        enddo
     enddo
 
     ! Update the T matrix [Eq. (d)]
     ! First the row and then the column
 
     do j= 2,mixer%memory+1
        vtemp = distdot(dim,df(1,mixer%memory),1,df(1,j-1),1, &
             parallel%group_size,parallel%group_comm)*weight
        mixer%t(mixer%memory,j-1) = vtemp
        mixer%t(j-1,mixer%memory) = vtemp
     enddo

     ! find the inverse of T
     do j = 1,mixer%memory
        do i = 1,mixer%memory
           mixer%tinv(j,i) = mixer%t(j,i)
        enddo
     enddo
     ipiv = 0

     call dgetrf(mixer%memory,mixer%memory,mixer%tinv, &
          mixer%memory,ipiv,info)
     if (info < 0) then
        write(9,*) 'ERROR: problem with mixer dgetrf'
        write(9,*) 'STOP in genbro, info = ',info
        ierr = 453
        call myflush(9)
        return
     endif
     if (info > 0) then
        write(9,*) 'Warning: T matrix in mixer is Singular'
     endif

     call dgetri(mixer%memory,mixer%tinv,mixer%memory,ipiv,work, &
          mixer%memory,info)
     if (info < 0) then
        write(9,*) 'ERROR: problem with mixer dgetrf'
        write(9,*) 'STOP in genbro, info = ',info
        ierr = 454
        call myflush(9)
        return
     endif
     if (info > 0) then
        write(9,*) 'Warning: T matrix in mixer is Singular'
     endif

     ! update current v vector (stored in v(i,blk+1)
     do i = 1, dim
        do j = 1,mixer%memory
           v(i,blk+1) = v(i,blk+1) + mixer%tinv(mixer%memory,j)*df(i,j)
        enddo
     enddo

     ! find out how many times the u,v blocks of memory will be read,
     ! using an m value appropriate for READING (bitr-1). Each block
     ! is stored in a different file (bro.1, bro.2,...) and read
     ! separately. This conserved memory and minimizes writes to disk.

     nblk = (m-1)/(blk-1)+1
     ! prepare the m value appropriate for further COMPUTING, bitr
     m = bitr
 
     ! for each of time we access the block:
     do nn = 1, nblk
 
        ! create the file name and open it
        write(fnam,20) nn
20      format('bro.',i1,'.')
        open(10,file=fnam//idstring,form='unformatted')
        rewind(10)

        ! read old u.v vectors
        if (nn /= nblk) then
           rsize = blk-1
        else
           rsize = bitr-1-(nblk-1)*(blk-1)
        endif
        do j = 2, rsize
           read(10) (u(i,j),i=1,dim)
        enddo
        do j = 2, rsize
           read(10) (v(i,j),i=1,dim)
        enddo
        close(10)
 
        ! find current u (stored in u(i,blk+1)
        do j = 2, rsize
           vfac = distdot(dim,v(1,j),1,df(1,mixer%memory),1, &
                parallel%group_size,parallel%group_comm)*weight
           do i = 1, dim
              u(i,blk+1) = u(i,blk+1) - vfac*u(i,j)
           enddo
        enddo

        ! update x_out according to the update of B_1 to B_(n-1) 

        do j = 2, rsize
           vfac = (distdot(dim,v(1,j),1,mixer%old2,1, &
                parallel%group_size,parallel%group_comm) - &
                distdot(dim,v(1,j),1,xin,1, &
                parallel%group_size,parallel%group_comm) )*weight
           do i = 1, dim
              xout(i) = xout(i) + vfac*u(i,j)
           enddo
        enddo

        ! If this is the last block to be read, write last block of u, v
        ! for the next iteration.

        if (nn == nblk) then
           if (rsize /= blk-1) then
              open(10,file=fnam//idstring,form='unformatted')
              rewind (10)
              do j = 2, rsize
                 write(10) (u(i,j),i=1,dim)
              enddo
              write(10) (u(i,blk+1),i=1,dim)
              do j = 2, rsize
                 write(10) (v(i,j),i=1,dim)
              enddo
              write(10) (v(i,blk+1),i=1,dim)
              close(10)
           else
              write(fnam,20) nn+1
              open(10,file=fnam//idstring,form='unformatted')
              rewind(10)
              write(10) (u(i,blk+1),i=1,dim)
              write(10) (v(i,blk+1),i=1,dim)
              close(10)
           endif
        endif
     enddo

     ! write df vectors for next time
     open(11,file='bro.df.'//idstring,form='unformatted')
     do j = 2,mixer%memory
        write(11) (df(i,j),i=1,dim)
     enddo
     close(11)

     ! update x_out according to the update of B_(n-1) to B_n (i.e.,
     ! use the current u,v, which weren't in the loop before)

     vfac = (distdot(dim,v(1,blk+1),1,mixer%old2,1, &
          parallel%group_size,parallel%group_comm) - &
          distdot(dim,v(1,blk+1),1,xin,1, &
          parallel%group_size,parallel%group_comm) )*weight
     do i = 1, dim
        xout(i) = xout(i) + vfac*u(i,blk+1)
     enddo

     ! complete the last fragment of Eq. (f)
     do i = 1, dim
        xout(i) = xin(i) - xout(i)
     enddo
  endif

99 continue
  deallocate(u)
  deallocate(v)
  deallocate(df)

end subroutine genbro
!===============================================================
