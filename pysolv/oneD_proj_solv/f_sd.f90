!~~~~~~~~~~~~~~~~~~!
! MODULE f_sdsolve !
!~~~~~~~~~~~~~~~~~~!
module f_sdsolve
! module that contains steepest descent solver related functions

    implicit none

    ! module variables
    integer :: iter_
    double precision :: res_

    contains
        !~~~~~~~~~~~~~~~~~~!
        ! SUBROUTINE solve !
        !~~~~~~~~~~~~~~~~~~!
        subroutine solve(A, b, x, itermax, tol, n)
        ! subroutine that performs the steepest descent iterations

            implicit none

            ! inputs
            double precision, intent(in), dimension(n, n) :: A
            double precision, intent(in), dimension(n) :: b
            double precision, intent(in) :: tol
            integer, intent(in) :: itermax

            ! outputs
            double precision, intent(inout), dimension(n) :: x

            ! local
            integer :: n, iter
            double precision, dimension(n) :: dir, res, Adx
            double precision :: rdr, ddr, alpha

            ! f2py intent(hide), depend(A) :: n=shape(A, 0)

            ! initalize local variables
            iter = 0
            Adx = 0
            rdr = 0
            ddr = 0
            dir = 0

            ! intial search direction and residual
            call mat_vec_mul(A, x, Adx, n)
            res = b - Adx
            call mat_vec_mul(A, res, dir, n)
            Adx = 0

            ! continue the iteration till the iteration counter reaches the maximum count if convergence is not
            ! obtained
            do while (iter <= itermax)
                ! compute the step size alpha
                call vec_vec_mul(res, res, rdr, n)
                call vec_vec_mul(dir, res, ddr, n)
                alpha = rdr / ddr

                ! update the solution
                x = x + (alpha*res)

                ! update the residual
                res = res - (alpha*dir)

                ! check for convergence
                if (norm2(res) <= tol) then
                    exit
                end if

                ! update the search direction
                dir = 0
                call mat_vec_mul(A, res, dir, n)

                ! update the iteration counter
                iter = iter + 1

                ! reset the local variables
                rdr = 0
                ddr = 0
            end do

            iter_ = iter
            res_ = norm2(res)

        end subroutine solve

end module f_sdsolve