!~~~~~~~~~~~~~~~~~~~~~~!
! MODULE f_jacobisolve !
!~~~~~~~~~~~~~~~~~~~~~~!
module f_jacobisolve
! module to perform Jacobi solver related functions

    implicit none

    integer :: iter_
    double precision :: res_

    contains
        !~~~~~~~~~~~~~~~~~~!
        ! SUBROUTINE solve !
        !~~~~~~~~~~~~~~~~~~!
        subroutine solve(A, b, x0, x, itermax, tol, omega, n)
        ! subroutine to perform Jacobi iteration

            implicit none

            ! inputs
            double precision, intent(in), dimension(n, n) :: A
            double precision, intent(in), dimension(n) :: b, x0
            double precision, intent(in) :: tol, omega
            integer, intent(in) :: itermax

            ! outputs
            double precision, intent(inout), dimension(n) :: x

            ! local
            double precision :: sum_
            integer :: iter, i, j, n
            double precision, dimension(n) :: Ax, x_old, res

            ! f2py intent(hide), depend(A) :: n=shape(A, 0)

            ! initialize local variables
            iter = 0
            Ax = 0
            x_old = x0
            res = 0

            ! continue the iteration till the iteration counter reaches the maximum count if convergence is not
            ! obtained
            do while (iter <= itermax)
                do i = 1, n
                    ! (b_i - sum(a_i_j * x_j(k)))
                    sum_ = 0

                    do j = 1, n
                        if (.not. j == i) then
                            sum_ = sum_ + (A(i, j) * x_old(j))
                        end if
                    end do

                    sum_ = b(i) - sum_
                    ! x_i(k) = (omega/a_i_i) * sum
                    x(i) = (omega/A(i, i)) * sum_
                end do

                ! stopping criteria: ||b - A.x|| < TOL
                call mat_vec_mul(A, x, Ax, n)
                res = b - Ax
                if (norm2(res) <= tol) then
                    exit
                end if

                ! update A.x for the next iteration
                Ax = 0

                ! update the iteration counter and x(iter-1)
                x_old = x
                iter = iter + 1
            end do

            iter_ = iter
            res_ = norm2(res)

        end subroutine solve

end module f_jacobisolve