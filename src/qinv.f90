!===============================================================
! Quasi-inverse by Householder QR factorization with
! column-pivoting.
!
! H.-r. Fang (2007).
!---------------------------------------------------------------

subroutine qinv(n, A, Ainv, info)

    use constants
    implicit none

    !---------------------------
    ! Input/Output variables:
    !---------------------------
    ! Dimensions of A.
    integer, intent(in) :: n, info
    ! A is an m-by-n real matrix.
    real(dp), intent(in) :: A(n,n)
    ! Ainv is an n-by-m real matrix.
    real(dp), intent(out) :: Ainv(n,n)

    !---------------------------
    ! Local variables:
    !---------------------------
    integer :: jpvt(n), i, j, k, r
    real(dp) :: tol, tmp, tau(n), work(4*n), Q(n,n), Rinv(n,n), u(n), v(n)

    jpvt = zero
    call dgeqp3(n, n, A, n, jpvt, tau, work, 4*n, info)
    if (info /= 0) then
        return
    endif
    tol = abs(A(1,1))*2.2204e-16
    ! 2.2204e-16 is eps using IEEE754 arithmetic.
    r = 0
    do while (abs(A(r+1,r+1))>tol .AND. r<n)
        r = r+1
    enddo
    ! r is the (numerical) rank of A.
    ! if (r<n) write(*,*) 'Ill-conditioned system is encountered (n=', n, 'r=', r, ').'
    Q = zero
    do i = 1,n
        Q(i,i) = 1.0
    enddo
    do j = r,1,-1
        u(j) = 1.0
        u(j+1:n) = A(j+1:n,j)
        do i = j,n
            v(i) = 0.0
            do k = j,n
                v(i) = v(i) + u(k)*Q(k,i)
            enddo
        enddo
        do i = j,n
           u(i) = tau(j)*u(i)
        enddo
        do i = j,n
            do k = j,n
                Q(i,k) = Q(i,k) - u(i)*v(k)
            enddo
        enddo
    enddo
    ! Let AP=QR be the Householder QR factorization of A.
    ! After calling dgeqp3, R is stored in the upper triangular part of A.
    ! The information of permutation matrix P for pivoting is stored in jpvt.
    ! Q is now reconstructed. P is the permutation matrix for pivoting,

    Rinv = zero
    do j=r,1,-1
        Rinv(j,j) = 1/A(j,j)
        do i=j-1,1,-1
            tmp = 0.0
            do k=i+1,r
                tmp = tmp + A(i,k)*Rinv(k,j)
            enddo
            Rinv(i,j) = -tmp/A(i,i)
        enddo
    enddo
    ! R is stored in the upper triangular part of A.
    ! Rinv is the inverse of R.

    Ainv = zero
    do i=1,r
        do j=1,n
            do k=i,r
                Ainv(jpvt(i),j) = Ainv(jpvt(i),j) + Rinv(i,k)*Q(j,k)
            enddo
        enddo
    enddo

end subroutine qinv

!===============================================================
