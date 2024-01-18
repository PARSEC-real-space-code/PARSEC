!===============================================================
! Multisecant methods for mixing, Type-I update on approximate
! Jacobian.
!
! Reference:
! H.-r. Fang and Y. Saad, "Two classes of multisecant methods
! for nonlinear acceleration." Report umsi-2007-100, Minnesota
! Supercomputer Institute, University of Minnesota, Minneapolis,
! MN, 2007.
!
! This implementation was with technical help of M. L. Tiago.
!---------------------------------------------------------------

subroutine msecant1_mix(mixer,parallel,iter,nrep,xin,xout,ierr)

    use constants
    use mixer_module
    use parallel_data_module
    implicit none

    !---------------------------
    ! Input/Output variables:
    !---------------------------
    ! Mixer related data.
    type (mixer_data), intent(inout) :: mixer
    ! Parallel computation related data.
    type (parallel_data), intent(in) :: parallel

    ! Iteration number.
    integer, intent(in) :: iter
    ! Order of the reduced group (not used).
    integer, intent(in) :: nrep
    ! Error flag.
    integer, intent(out) :: ierr

    ! Old potential.
    real(dp), intent(inout) :: xin(mixer%ndim * mixer%nspin)
    ! New potential.
    real(dp), intent(inout) :: xout(mixer%ndim * mixer%nspin)

    !---------------------------
    ! Local variables:
    !---------------------------
    ! Dimension of the potential vector.
    integer :: dim
    ! Residual vector f = xout - xin.
    ! Secant vector xx = xin_new - xin_old.
    real(dp) :: f(mixer%ndim * mixer%nspin), xx(mixer%ndim * mixer%nspin)
    real(dp), allocatable :: tmp(:)
    ! Counters, temporary variables, etc.
    integer :: i, j, k, m
    integer :: sz, res, ngroup
    integer :: info, alcstat
    ! For restarting scheme.
    logical :: restart
    real(dp) :: lhs, rhs
    real(dp),dimension(1) :: lhsvec, rhsvec

    !---------------------------
    ! Initialization:
    !---------------------------
    dim = mixer%ndim * mixer%nspin
    f = xout - xin
    if (iter == 1) then
        mixer%nsecant = 0
    endif
    m = mixer%nsecant  ! Number of available secant equations.
    if (m == mixer%memory) then
        ! Expand storage dx and df.
        allocate (mixer%t(dim,m), stat=alcstat)
        call alccheck('mixer%t', dim*m, alcstat)
        mixer%memory = max(mixer%memory+1, int(mixer%memory*mixer%expand_factor))

        ! Expand dx.
        mixer%t = mixer%dx
        deallocate (mixer%dx)
        allocate (mixer%dx(dim,mixer%memory), stat=alcstat)
        call alccheck('mixer%dx', dim*mixer%memory, alcstat)
        mixer%dx(1:dim,1:m) = mixer%t(1:dim,1:m)

        ! Expand df.
        mixer%t = mixer%df
        deallocate (mixer%df)
        allocate (mixer%df(dim,mixer%memory), stat=alcstat)
        call alccheck('mixer%df', dim*mixer%memory, alcstat)
        mixer%df(1:dim,1:m) = mixer%t(1:dim,1:m)

        ! Expand N, if it is appropriate.
        if (mixer%group_size == 0) then
            mixer%t = mixer%n
            allocate (mixer%n(dim,mixer%memory), stat=alcstat)
            call alccheck('mixer%n', dim*mixer%memory, alcstat)
            mixer%n(1:dim,1:m) = mixer%t(1:dim,1:m)
        endif

        deallocate (mixer%t)
    endif

    if (mixer%group_size == 0) then
        sz = m+1  ! Take all available iterates in one group.
    else
        sz = mixer%group_size
    endif
    allocate (tmp(sz), stat=alcstat)
    call alccheck('tmp', sz, alcstat)

    !---------------------------
    ! Mixing:
    !---------------------------
    if (m == 0) then
        ! No secant information available yet (the first iteration).
        ! Store xin and f for next iteration.
        mixer%dx(:,1) = xin(:)
        mixer%df(:,1) = f(:)
        ! Next guess obtained by simple mixing.
        xin(:) = xin(:) + mixer%param * f(:)
        xout = xin
        mixer%nsecant = 1
        if (mixer%en_stage == 1) then
            mixer%en_stage = 2
        endif

    else
        ! When the current iterate is bad, perform restart.
        ! More precisely, if ||f0|| < restart_factor*||f1||, then restart,
        ! where f0/f1 are the previous/current function values, respectively.
        ! In particular, restart_factor = 0 implies never restart.
        restart = .false.
        if ((m==1 .and. mixer%en_stage==1) .or. m>=2) then
            lhs = dot_product(mixer%f0,mixer%f0)
            lhsvec = lhs
            call psum(lhsvec, 1, parallel%group_size, parallel%group_comm)
            lhs = lhsvec(1)
            rhs = dot_product(f,f)
            rhsvec = rhs
            call psum(rhsvec, 1, parallel%group_size, parallel%group_comm)
            rhs = rhsvec(1)
            if (lhs < mixer%restart_factor*mixer%restart_factor*rhs) then
               restart = .true.
            endif
        endif
        if (restart) then
            mixer%dx(:,1) = mixer%x0(:)
            mixer%df(:,1) = mixer%f0(:)
            ! Set the new estimate by simple mixing with the previous iterate.
            xin(:) = mixer%x0(:) + mixer%param * mixer%f0(:)
            xout = xin
            mixer%nsecant = 1
            if (mixer%en_stage == 1) then
                mixer%en_stage = 2
            endif

        else
            res = mod(m+sz-1,sz) + 1
            ! res is the size of the last group.
            ngroup = (m-res)/sz
            ! ngroup does not count the last group; the # of groups is ngroup+1.
            if (mixer%en_stage /= 1) then
                mixer%dx(:,m) = xin(:) - mixer%dx(:,m)
                mixer%df(:,m) = f(:) - mixer%df(:,m)
                xx(:) = mixer%dx(:,m)

                mixer%dx(:,m) = mixer%dx(:,m) + mixer%param * mixer%df(:,m)
                do i=1,ngroup
                    ! Compute tmp = df(:,(i-1)*sz+1:i*sz)'*df(:,m).
                    tmp(:) = zero
                    do j=1,sz
                        do k=1,dim
                            tmp(j) = tmp(j) + mixer%df(k,(i-1)*sz+j)*mixer%df(k,m)
                        enddo
                    enddo
                    call psum(tmp, sz, parallel%group_size, parallel%group_comm)
                    ! Compute dx(:,m) = dx(:,m) - dx(:,(i-1)*sz+1:i*sz)*tmp.
                    do j=1,dim
                        do k=1,sz
                            mixer%dx(j,m) = mixer%dx(j,m) - mixer%dx(j,(i-1)*sz+k)*tmp(k)
                        enddo
                    enddo
                enddo
                ! mixer%dx(:,m) is now E(:,m) in the note.
            endif

            ! Compute new x.
            xout(:) = xin(:) + mixer%param*f(:)
            do i=1,ngroup
                ! Compute tmp = df(:,(i-1)*sz+1:i*sz)'*f.
                tmp(:) = zero
                do j=1,sz
                    do k=1,dim
                        tmp(j) = tmp(j) + mixer%df(k,(i-1)*sz+j)*f(k)
                    enddo
                enddo
                call psum(tmp, sz, parallel%group_size, parallel%group_comm)
                ! Compute xout = xout - dx(:,(i-1)*sz+1:i*sz)*tmp.
                do j=1,dim
                    do k=1,sz
                        xout(j) = xout(j) - mixer%dx(j,(i-1)*sz+k)*tmp(k)
                    enddo
                enddo
                ! dx(:,...) is E(:,...) in the note.
                ! df(:,...) is V(:,...) in the note.
            enddo

            ! Now deal with the last group.
            if (mixer%en_stage /= 1) then
                mixer%n(:,res) = -mixer%param*xx
                do i=1,ngroup
                    ! Compute tmp = dx(:,(i-1)*sz+1:i*sz)'*xx.
                    tmp(:) = zero
                    do j=1,sz
                        do k=1,dim
                            tmp(j) = tmp(j) + mixer%dx(k,(i-1)*sz+j)*xx(k)
                        enddo
                    enddo
                    call psum(tmp, sz, parallel%group_size, parallel%group_comm)
                    ! Compute N(:,res) = N(:,res)+df(:,(i-1)*sz+1:i*sz)*tmp.
                    do j=1,dim
                        do k=1,sz
                            mixer%n(j,res) = mixer%n(j,res) + mixer%df(j,(i-1)*sz+k)*tmp(k)
                        enddo
                    enddo
                enddo
            endif
            if (mixer%en_stage==1 .and. res==sz) then
                ! Compute tmp = df(:,ngroup*sz+1:m)'*f.
                tmp(:) = zero
                do j=1,res
                    do k=1,dim
                        tmp(j) = tmp(j) + mixer%df(k,ngroup*sz+j)*f(k)
                    enddo
                enddo
                call psum(tmp(1:res), res, parallel%group_size, parallel%group_comm)
                ! Compute xout = xout - dx(:,ngroup*sz+1:m)*tmp.
                do j=1,dim
                    do k=1,res
                        xout(j) = xout(j) - mixer%dx(j,ngroup*sz+k)*tmp(k)
                    enddo
                enddo
            else
                ! Compute T = N(:,1:res)'*df(:,m-res+1:m).
                ! Note that m = ngroup*sz + res.
                allocate (mixer%t(res,res), stat=alcstat)
                call alccheck('mixer%t', res*res, alcstat)
                mixer%t(:,:) = zero
                do i=1,res
                    do j=1,res
                        do k=1,dim
                            mixer%t(i,j) = mixer%t(i,j) + mixer%n(k,i)*mixer%df(k,ngroup*sz+j)
                        enddo
                    enddo
                enddo
                call psum(mixer%t, res*res, parallel%group_size, parallel%group_comm)
                allocate (mixer%tinv(res,res), stat=alcstat)
                call alccheck('mixer%tinv', res*res, alcstat)
                call qinv(res, mixer%t, mixer%tinv, info)
                if (info /= 0) then
                    write(9,*) ' ERROR: dgeqp3 invoked by qinv failed in'
                    write(9,*) ' subroutine msecant1, info = ', info
                    ierr = 1
                else
                    ! Compute T = N*Tinv'.
                    deallocate (mixer%t)
                    allocate (mixer%t(dim,res), stat=alcstat)
                    call alccheck('mixer%t', dim*res, alcstat)
                    mixer%t(:,:) = zero
                    do i=1,dim
                        do j=1,res
                            do k=1,res
                                mixer%t(i,j) = mixer%t(i,j) + mixer%n(i,k)*mixer%tinv(j,k)
                            enddo
                        enddo
                    enddo
                    ! Compute tmp = T'*f.
                    tmp(:) = zero
                    do i=1,res
                        do j=1,dim
                            tmp(i) = tmp(i) + mixer%t(j,i)*f(j)
                        enddo
                    enddo
                    call psum(tmp(1:res), res, parallel%group_size, parallel%group_comm)
                    ! Compute xout = xout - dx(:,m-res+1:m)*tmp.
                    do i=1,dim
                        do j=1,res
                            xout(i) = xout(i) - mixer%dx(i,m-res+j)*tmp(j)
                        enddo
                    enddo
                    if (res == sz) then
                        mixer%df(:,m-sz+1:m) = mixer%t(:,:)
                        ! df(:,m-sz+1:m) is now V(:,m-sz+1:m) in the note.
                    endif
                endif
                deallocate (mixer%t)
                deallocate (mixer%tinv)
            endif

            mixer%x0 = xin
            mixer%f0 = f
            if (mixer%en_stage /= 2) then
                mixer%dx(:,m+1) = xin(:)
                mixer%df(:,m+1) = f(:)
                mixer%nsecant = m+1
            endif
            ! xin = xout
            if (mixer%en_stage == 1) then
                mixer%en_stage = 2  ! Next EN stage is 1.
            else if (mixer%en_stage == 2) then
                mixer%en_stage = 1  ! Next EN stage is 2.
            endif
        endif
    endif

    deallocate (tmp)

end subroutine msecant1_mix
!===============================================================
