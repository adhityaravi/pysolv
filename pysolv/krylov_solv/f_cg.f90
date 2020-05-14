!~~~~~~~~~~~~~~~~~~!
! MODULE f_cgsolve !
!~~~~~~~~~~~~~~~~~~!
module f_cgsolve
! module for performing conjugate-gradient related functions

    implicit none

    ! variables declared for host association with the contained subroutines
    double precision, allocatable, dimension(:, :) :: A
    double precision, allocatable, dimension(:) :: b
    double precision :: tol, res_
    integer :: itermax, n, iter_

    contains
        !~~~~~~~~~~~~~~~~~!
        ! SUBROUTINE init !
        !~~~~~~~~~~~~~~~~~!
        subroutine init(A_in, b_in, tol_in, itermax_in, n_in)
        ! subroutine for initiating py2f inputs that will be shared by other subroutines in the module

            implicit none

            ! inputs
            double precision, intent(in), dimension(n_in, n_in) :: A_in
            double precision, intent(in), dimension(n_in) :: b_in
            double precision, intent(in) :: tol_in
            integer, intent(in) :: itermax_in

            ! local
            integer :: n_in

            ! f2py intent(hide), depend(A_in) :: n_in=shape(A_in, 0)

            ! copy local inputs into corresponding the module variables
            ! ints and floats
            tol = tol_in
            itermax = itermax_in
            n = n_in

            ! allocate sizes for global arrays
            ! sanity check
            if (allocated(A)) then
                deallocate(A)
                deallocate(b)
            end if
            allocate(A(n, n))
            allocate(b(n))
            ! arrays
            A = A_in
            b = b_in

        end subroutine init

        !~~~~~~~~~~~~~~~~~~!
        ! SUBROUTINE solve !
        !~~~~~~~~~~~~~~~~~~!
        subroutine solve(x)
        ! subroutine to perform the conjugate gradient iterations

            implicit none

            ! ouputs
            double precision, intent(inout), dimension(:) :: x

            ! local
            double precision, dimension(n) :: res, dir, Adx, Add
            double precision :: alpha, beta, rdr, ddAdd, rndrn
            integer :: iter

            ! initializing local variables
            iter = 0
            Adx = 0
            Add = 0
            ddAdd = 0
            rdr = 0
            rndrn = 0


            ! initial residual and search direction
            ! d = r = b - A.x_0
            call mat_vec_mul(A, x, Adx, n)
            res = b - Adx
            dir = res

            ! continue the iteration till the iteration counter reaches the maximum count if convergence is not
            ! obtained
            do while (iter < itermax)
               ! compute the step size:
               ! A.dir
               call mat_vec_mul(A, dir, Add, n)
               ! dir.(A.dir)
               call vec_vec_mul(dir, Add, ddAdd, n)
               ! res.res
               call vec_vec_mul(res, res, rdr, n)
               ! alpha = res.res / dir.(A.dir)
               alpha = rdr / ddAdd

               ! compute the new solution
               x = x + (alpha*dir)

               ! update the residual
               res = res - (alpha*Add)
               ! check for convergence
               if (norm2(res) <= tol) then
                   exit
               end if

               ! compute beta coefficient (from Gram-Schmidt orthogonalization)
               ! res_nxt.res_nxt
               call vec_vec_mul(res, res, rndrn, n)
               ! beta = (res_nxt.res_nxt) / (res.res)
               beta = rndrn / rdr

               ! update the search direction
               dir = res + (beta*dir)

               ! update the iteration counter
               iter = iter + 1

               ! reset the local variables
               Adx = 0
               Add = 0
               ddAdd = 0
               rdr = 0
               rndrn = 0

            end do

            iter_ = iter
            res_ = norm2(res_nxt)

        end subroutine solve

end module f_cgsolve