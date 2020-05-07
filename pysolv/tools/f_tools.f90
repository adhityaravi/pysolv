!~~~~~~~~~~~~~~~~~~~~~~~~!
! SUBROUTINE mat_vec_mul !
!~~~~~~~~~~~~~~~~~~~~~~~~!
subroutine mat_vec_mul(a, b, c, n)
! subroutine to compute matrix vector inner product

    implicit none

    ! local
    double precision :: alpha, beta
    integer :: incx, incy

    ! inputs
    double precision, intent(in), dimension(n, n) :: a
    double precision, intent(in), dimension(n) :: b
    integer, intent(in) :: n

    ! outputs
    double precision, intent(inout), dimension(n) :: c

    ! blas matrix-vector multiplication parameters
    alpha = 1
    beta = 1
    incx = 1
    incy = 1

    ! compute matrix vetor inner product (using blas)
    call dgemv('N', n, n, alpha, a, n, b, incx, beta, c, incy)

end subroutine mat_vec_mul


!~~~~~~~~~~~~~~~~~~~~~~~~!
! SUBROUTINE vec_vec_mul !
!~~~~~~~~~~~~~~~~~~~~~~~~!
subroutine vec_vec_mul(a, b, c, n)
! subroutine to compute matrix vector inner product

    implicit none

    ! local
    integer :: incx, incy

    ! inputs
    double precision, intent(in), dimension(n) :: a, b
    integer, intent(in) :: n

    ! outputs
    double precision, intent(inout) :: c

    ! import ddot from blas
    double precision, external :: ddot

    ! blas matrix-vector multiplication parameters
    incx = 1
    incy = 1

    ! compute vetor inner product
    c = ddot(n, a, incx, b, incy)

end subroutine vec_vec_mul