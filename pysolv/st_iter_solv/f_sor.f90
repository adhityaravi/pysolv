!~~~~~~~~~~~~~~~~~~~!
! MODULE f_sorsolve !
!~~~~~~~~~~~~~~~~~~~!
module f_sorsolve
! module for performing successive over-relaxation solver related functions

    implicit none

    ! variables declared for host association with the contained subroutines
    double precision, allocatable, dimension(:, :) :: A
    double precision, allocatable, dimension(:) :: b, x0
    double precision :: tol, omega, h, c1, c2, lambda1, lambda2, rho1, res_
    integer :: itermax, adaptive_omega, omega_update_frequency, n, iter_

    contains
        !~~~~~~~~~~~~~~~~~!
        ! SUBROUTINE init !
        !~~~~~~~~~~~~~~~~~!
        subroutine init(A_in, b_in, x0_in, itermax_in, tol_in, omega_in, h_in, adaptive_omega_in, &
                        omega_update_frequency_in, c1_in, c2_in, lambda1_in, lambda2_in, rho1_in, n_in)
        !subroutine for initiating py2f inputs that will be shared by other subroutines in the module

            implicit none

            ! inputs
            double precision, intent(in), dimension(n_in, n_in) :: A_in
            double precision, intent(in), dimension(n_in) :: b_in, x0_in
            integer, intent(in) :: itermax_in, adaptive_omega_in, omega_update_frequency_in
            double precision, intent(in) :: tol_in, omega_in, h_in, c1_in, c2_in, lambda1_in, lambda2_in, rho1_in

            ! local
            integer :: n_in

            ! f2py intent(hide), depend(A_in) :: n_in=shape(A_in, 0)

            ! copy local inputs into corresponding the module variables
            ! ints and floats
            n = n_in
            itermax = itermax_in
            adaptive_omega = adaptive_omega_in
            omega_update_frequency = omega_update_frequency_in
            tol = tol_in
            omega = omega_in
            h = h_in
            c1 = c1_in
            c2 = c2_in
            lambda1 = lambda1_in
            lambda2 = lambda2_in
            rho1 = rho1_in

            ! allocate sizes for global arrays
            ! sanity check
            if (allocated(A)) then
                deallocate(A)
                deallocate(b)
                deallocate(x0)
            end if
            allocate(A(n, n))
            allocate(b(n))
            allocate(x0(n))
            ! arrays
            A = A_in
            b = b_in
            x0 = x0_in

        end subroutine init

        !~~~~~~~~~~~~~~~~~~!
        ! SUBROUTINE solve !
        !~~~~~~~~~~~~~~~~~~!
        subroutine solve(x)
        ! subroutine for performing the sor iteration

            implicit none

            ! outputs
            double precision, intent(inout), dimension(:) :: x

            ! local
            double precision :: sum_
            integer :: iter, i, j
            double precision, dimension(n) :: Ax, x_old, res, res_old

            ! Initilializing local variables
            Ax = 0
            x_old = x0
            iter = 0
            res = 0
            res_old = 0

            ! continue the iteration till the iteration counter reaches the maximum count if convergence is not
            ! obtained
            do while (iter <= itermax)
                do i = 1, n
                    sum_ = 0

                    ! sum_(j=1) ^ (i - 1)(a_i_j * x_j(iter))
                    do j = 1, i-1
                        sum_ = sum_ + (A(i, j)*x(j))
                    end do

                    ! sum_(j=i + 1) ^ (n)(A_i_j * x_j(iter - 1))
                    do j = i+1, n
                        sum_ = sum_ + (A(i, j)*x_old(j))
                    end do

                    ! x_i(iter) = (1-omega) * x_i(iter-1) + ((omega/a_i_i) * (b_i - sum_)
                    x(i) = ((1-omega)*x_old(i)) + ((omega / A(i, i)) * (b(i) - sum_))
                end do

                ! stopping criteria: ||b - A.x|| < TOL
                call mat_vec_mul(A, x, Ax, n)
                res = b - Ax
                if (norm2(res) <= tol) then
                    exit
                end if

                ! update relaxation parameter
                if ((mod(iter, omega_update_frequency) == 0) .and. (.not. iter == 0)) then
                    if (.not. adaptive_omega == 0) then
                        call update_omega(x, x_old, res, res_old)
                    end if
                end if

                ! update residual for the next iteration
                res_old = res
                Ax = 0

                ! update the iteration counter and x(iter-1)
                x_old = x
                iter = iter + 1
            end do

            iter_ = iter
            res_ = norm2(res)

        end subroutine solve

        !~~~~~~~~~~~~~~~~~~~!
        ! SUBROUTINE s_solve !
        !~~~~~~~~~~~~~~~~~~~!
        subroutine s_solve(x)
        ! subroutine to perform symmetric sor iteration

            implicit none

            ! outputs
            double precision, intent(inout), dimension(:) :: x

            ! local
            double precision :: sum_
            integer :: iter, i, j
            double precision, dimension(n) :: Ax, x_old, res, res_old, x_half

            ! Initilializing local variables
            Ax = 0
            x_old = x0
            iter = 0
            res = 0
            res_old = 0
            x_half = 0

            ! continue the iteration till the iteration counter reaches the maximum count if convergence is not
            ! obtained
            do while (iter <= itermax)
                ! Forward sweep
                do i = 1, n
                    sum_ = 0

                    ! sum_(j=1) ^ (i - 1)(a_i_j * x_j(iter))
                    do j = 1, i-1
                        sum_ = sum_ + (A(i, j)*x_half(j))
                    end do

                    ! sum_(j=i + 1) ^ (n)(A_i_j * x_j(iter - 1))
                    do j = i+1, n
                        sum_ = sum_ + (A(i, j)*x_old(j))
                    end do

                    ! x_i(half) = (1-omega) * x_i(iter-1) + ((omega/a_i_i) * (b_i - sum_)
                    x_half(i) = ((1-omega)*x_old(i)) + ((omega / A(i, i)) * (b(i) - sum_))
                end do

                ! Backward sweep
                do i = n, 1, -1
                    sum_ = 0

                    ! sum_(j=1) ^ (i - 1)(a_i_j * x_j(iter))
                    do j = 1, i-1
                        sum_ = sum_ + (A(i, j)*x_half(j))
                    end do

                    ! sum_(j=i + 1) ^ (n)(A_i_j * x_j(iter - 1))
                    do j = i+1, n
                        sum_ = sum_ + (A(i, j)*x(j))
                    end do

                    ! x_i(iter) = (1-omega) * x_i(half) + ((omega/a_i_i) * (b_i - sum_)
                    x(i) = ((1-omega)*x_half(i)) + ((omega / A(i, i)) * (b(i) - sum_))
                end do

                ! stopping criteria: ||b - A.x|| < TOL
                call mat_vec_mul(A, x, Ax, n)
                res = b - Ax
                if (norm2(res) <= tol) then
                    exit
                end if

                ! update relaxation parameter
                if ((mod(iter, omega_update_frequency) == 0) .and. (.not. iter == 0)) then
                    if (.not. adaptive_omega == 0) then
                        call update_omega(x, x_old, res, res_old)
                    end if
                end if

                ! update residual for the next iteration
                res_old = res
                Ax = 0

                ! update the iteration counter and x(iter-1)
                x_old = x
                iter = iter + 1
            end do

            iter_ = iter
            res_ = norm2(res)

        end subroutine s_solve

        !~~~~~~~~~~~~~~~~~~~~~~~~~!
        ! SUBROUTINE update_omega !
        !~~~~~~~~~~~~~~~~~~~~~~~~~!
        subroutine update_omega(x, x_old, res, res_old)
        ! subroutine to update the relaxation parameter adaptively

            implicit none

            ! inputs
            double precision, intent(in), dimension(:) :: x, x_old, res, res_old

            ! call the correct updating subroutine
            if (adaptive_omega == 1) then
                call steepest_descent_update(res)

            elseif (adaptive_omega == 2) then
                call armijo_update(x, x_old, res_old)

            elseif (adaptive_omega == 3) then
                call wolfe_update(x, x_old, res, res_old)
            end if

        end subroutine update_omega

        !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
        ! SUBROUTINE steepest_descent_update !
        !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
        subroutine steepest_descent_update(res)
        ! subroutine to compute omega based on line search algorithm

            implicit none

            ! inputs
            double precision, intent(in), dimension(:) :: res

            ! local
            double precision :: rdr, rdAdr, h
            double precision, dimension(n) :: Adr

            ! initialize local variables
            rdr = 0
            rdAdr = 0
            h = 0
            Adr = 0

            ! residual.residual
            call vec_vec_mul(res, res, rdr, n)
            ! A.residual
            call mat_vec_mul(A, res, Adr, n)
            ! residual.(A.residual)
            call vec_vec_mul(res, Adr, rdAdr, n)

            ! h = residual.residual / residual.(A.residual)
            h = rdr / rdAdr
            ! omega = 2*h / (2 + h)
            omega = (2 * h) / (2 + h)
            ! omega safety check
            if ((omega > 1.95) .or. (omega < 0.05)) then
                omega = 1
            end if

        end subroutine steepest_descent_update

        !~~~~~~~~~~~~~~~~~~~~~~~~~~!
        ! SUBROUTINE armijo_update !
        !~~~~~~~~~~~~~~~~~~~~~~~~~~!
        subroutine armijo_update(x, x_old, res_old)
        ! subroutine to compute relaxation parameter based on Armijo conditions

            implicit none

            ! inputs
            double precision, intent(in), dimension(:) :: x, x_old, res_old

            ! local
            double precision, dimension(n) :: diff_x
            double precision :: f_dash_old_dot_diff_x, fx, fx_old

            ! initialize local variables
            diff_x = 0
            f_dash_old_dot_diff_x = 0
            fx = 0
            fx_old = 0

            ! x - x_old
            diff_x = x - x_old
            ! f'(x(iter-1)) . (x(iter) - x(iter-1))
            call vec_vec_mul(-res_old, diff_x, f_dash_old_dot_diff_x, n)
            ! f(x)
            call f(x, fx)
            ! f(x_old)
            call f(x_old, fx_old)

            ! if f(x(iter)) <= f(x(iter-1)) + c1*(f'(x(iter-1)).(x(iter) - x(iter-1)))
            if (fx <= (fx_old + c1*f_dash_old_dot_diff_x)) then
                h = lambda1 * h
            else
                h = rho1 * h
            end if

            ! omega = 2*h / (2 + h)
            omega = (2 * h) / (2 + h)
            ! omega safety check
            if ((omega > 1.95) .or. (omega < 0.05)) then
                omega = 1
            end if

        end subroutine armijo_update

        !~~~~~~~~~~~~~~~~~~~~~~~~~!
        ! SUBROUTINE wolfe_update !
        !~~~~~~~~~~~~~~~~~~~~~~~~~!
        subroutine wolfe_update(x, x_old, res, res_old)
        ! subroutine to compute relaxation parameter based on Wolfe condition

            implicit none

            ! inputs
            double precision, intent(in), dimension(:) :: x, x_old, res, res_old

            ! local
            double precision, dimension(n) :: diff_x
            double precision :: f_dash_dot_diff_x, f_dash_old_dot_diff_x, fx, fx_old

            ! initialize local variables
            diff_x = 0
            f_dash_old_dot_diff_x = 0
            f_dash_dot_diff_x = 0
            fx = 0
            fx_old = 0

            ! x - x_old
            diff_x = x - x_old
            ! f'(x(iter-1)) . (x(iter) - x(iter-1))
            call vec_vec_mul(-res_old, diff_x, f_dash_old_dot_diff_x, n)
            ! f'(x(iter)).(x(iter) - x(iter-1))
            call vec_vec_mul(-res, diff_x, f_dash_dot_diff_x, n)
            ! f(x)
            call f(x, fx)
            ! f(x_old)
            call f(x_old, fx_old)

            ! if f(x(iter)) <= f(x(iter-1)) + c1*(f'(x(iter-1)).(x(iter) - x(iter-1)))
            if (fx <= (fx_old + c1*f_dash_old_dot_diff_x)) then
                ! # if c2*(f'(x(iter-1)).(x(iter) - x(iter-1))) <= f'(x(iter)).(x(iter) - x(iter-1))
                if (c2*f_dash_old_dot_diff_x <= f_dash_dot_diff_x) then
                    ! h_new = lambda1 * h_old
                    h = lambda1 * h
                else
                    ! h_new = lambda2 * h_old
                    h = lambda2 * h
                end if
            else
                ! h_new = rho1 * h_old
                h = rho1 * h
            end if

            ! omega = 2*h / (2 + h)
            omega = (2 * h) / (2 + h)
            ! omega safety check
            if ((omega > 1.95) .or. (omega < 0.05)) then
                omega = 1
            end if

        end subroutine wolfe_update

        !~~~~~~~~~~~~~~!
        ! SUBROUTINE f !
        !~~~~~~~~~~~~~~!
        subroutine f(x, out)
        ! subroutine to compute 0.5 * (x.(A.x)) - (b.x)

            implicit none

            ! inputs
            double precision, intent(in), dimension(:) :: x

            ! outputs
            double precision, intent(inout) :: out

            ! local
            double precision, dimension(n) :: Adx
            double precision :: xdAdx, bdx

            ! initialize local variables
            Adx = 0
            xdAdx = 0
            bdx = 0

            ! A.x
            call mat_vec_mul(A, x, Adx, n)
            ! x.(A.x)
            call vec_vec_mul(Adx, x, xdAdx, n)
            ! b.x
            call vec_vec_mul(b, x, bdx, n)
            ! 0.5*(x.(A.x)) - b.x
            out = (0.5 * xdAdx) - (bdx)

        end subroutine f

end module f_sorsolve
