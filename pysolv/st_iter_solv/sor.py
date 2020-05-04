"""Successive over-relaxation linear solver.

Successive over-relaxation is a type of stationary iteration scheme for solving a set of linear equations.
The general linear iteration scheme reads:

x(k+1) = x(k) + (1/C)(b - Ax(k))

where A is the coefficient matrix
      b is the RHS vector
      x is the solution vector
      C is the approximation of A that is easier to invert
      k is the iteration number

In the Successive over-relaxation scheme, C is approximated by combining the lower and upper triangular matrices of A.
Then, a relaxation parameter omega is introduced for C. That is,

A = L + D + U; C = (1/omega) * (L + U)

Initial guess for the solution vector, if not mentioned, is taken as an unity vector.

Properties:
1. SOR is aimed at improving the convergence of the Gauss-Seidel solver by introducing the relaxation factor.

2. An apriori convergence criteria for the SOR solver is

    i. The A matrix is strictly or irreducibly diagonally dominant. Or

    ii. The A matrix is symmetric positive-definite (symmetric --> A = transpose(A); positive definite --> all the
        eigen values of A are positive)

3. A posteriori convergence criteria for the SOR scheme is: spectral radius of (1 - (1/omega*D + L)*A) is < 1.
"""

# import the necessary packages
from data import *
import numpy as np
import time as ti
import re


class SORSolve(Data):

    def __init__(self):
        """Initialize the class with the necessary inheritance and solve the linear system
        """

        # data inheritance
        Data.__init__(self)

        # initialize solution vector
        self._init_x0()

        # check for relaxation parameter
        self._check_omega()

        # call the gauss seidel solver
        self._solve()

    def _init_x0(self):
        """Check if an initial value exists for the solution vector or initialize the solution vector
        """

        try:
            self.x0
        except AttributeError:
            self.x0 = np.ones(self.n)

        # initialize the solution vectors
        self.x = np.empty(self.n)   # solution vector at the i'th iteration
        self.x_old = self.x0  # solution vector at the i-1'th iteration

    def _check_omega(self):
        """Check if the user has prescribed a value for the relaxation parameter. If not, prescribe a value of 1.7 based
           on Dr. Rosics's lecture (Institute of Scientific Computing, TU Braunschweig):
           [https://www.tu-braunschweig.de/en/wire/teaching/previous-terms/winter-2016-17]
        """

        # default relaxation parameter
        try:
            if type(self.omega) == str:
                self.ADAPTIVE_OMEGA = self.omega
                self.h = 2
                self.omega = 1
        except AttributeError:
            self.omega = 1.7

        # default number of iterations after which relaxation parameter will be updated
        try:
            self.omega_update_frequency
        except AttributeError:
            self.omega_update_frequency = 1

    def _solve(self):
        """Serial python implementation of successive over-relaxation solver based on Dr. Rosic's lectures (Institute of
           Scientific Computing, TU Braunschweig):
           [https://www.tu-braunschweig.de/en/wire/teaching/previous-terms/winter-2016-17]
        """

        # currently the sor solver of pysolv cannot handle non-square systems
        if self.m != self.n:
            msg = 'pysolv does not include a sor solver for non-square systems yet'
            raise ValueError(msg)

        # SOR iteration
        # x_i(iter) = (1 - omega) * x_i(iter-1) +
        #             (omega/a_i_i)(b_i - sum_(j=1)^(i-1)(a_i_j * x_j(iter)) - sum_(j=i+1)^(n)(A_i_j * x_j(iter-1)))
        else:
            # time the solver
            start_time = ti.time()

            # initialize iteration counter
            self.iter = 0

            # continue the iteration till the iteration counter reaches the maximum count if convergence is not
            # obtained
            while self.iter < self.ITERMAX:
                print(self.iter)
                print(self.omega)
                for i in range(0, self.n):
                    sum_ = 0

                    # sum_(j=1) ^ (i - 1)(a_i_j * x_j(iter))
                    for j in range(0, i):
                        sum_ = sum_ + (self.A[i, j] * self.x[j])

                    # sum_(j=i + 1) ^ (n)(A_i_j * x_j(iter - 1))
                    for j in range(i+1, self.n):
                        sum_ = sum_ + (self.A[i, j] * self.x_old[j])

                    # x_i(iter) = (1-omega) * x_i(iter-1) + ((omega/a_i_i) * (b_i - sum_)
                    self.x[i] = ((1-self.omega)*self.x_old[i]) + ((self.omega / self.A[i, i]) * (self.b[i] - sum_))

                # stopping criteria: (||x(iter) - x(iter-1)|| / ||x(iter)||) < TOL
                self.r = self.b - self.A.dot(self.x)
                if np.linalg.norm(self.r) <= np.linalg.norm(self.TOL*self.b):
                    break

                # update relaxation parameter
                if (self.iter % self.omega_update_frequency) == 0:
                    self._update_omega()

                # update the iteration counter and x(iter-1)
                self.x_old = self.x.copy()
                self.iter += 1
                print(np.linalg.norm(self.r))

            # time taken for convergence
            time_taken = ti.time() - start_time

            # add the computed information to the Data class
            Data.add_data('x', self.x)
            Data.add_data('time_taken', time_taken)
            Data.add_data('residual', np.linalg.norm(self.r))
            Data.add_data('iterations', self.iter-1)

    def _update_omega(self):
        """Update the relaxation parameter based on the selected type of parameter optimisation technique
        """

        if self.ADAPTIVE_OMEGA is None:
            pass

        elif re.sub('\W+', '', self.ADAPTIVE_OMEGA).casefold() == 'steepestdescent':
            self._steepest_descent_update()

        elif re.sub('\W+', '', self.ADAPTIVE_OMEGA).casefold() == 'wolfe':
            self._wolfe_update()

        elif re.sub('\W+', '', self.ADAPTIVE_OMEGA).casefold() == 'armijo':
            self._armijo_update()

        elif re.sub('\W+', '', self.ADAPTIVE_OMEGA).casefold() == 'experimental':
            self._exp_update()

        else:
            msg = 'Please select from available adaptive omega techniques: {}'.\
                  format(self.ADAPTIVE_OMEGA_FLAVOR)
            raise ValueError(msg)

    def _steepest_descent_update(self):
        """Update relaxation parameter based on steepest descent method from Miyatake et al.
           [https://www.sciencedirect.com/science/article/abs/pii/S037704271830181X?via%3Dihub]

           Steepest descent method for omega optimisation performs well for Poisson equations, but, might not be the
           optimum approach for all the linear systems. No additional parameter optimisation is required by this
           method. But, can get computationally expensive as this method includes a matrix-vector product computation.
           For further details, refer to Miyatake's paper.
        """

        # r(iter) = b - Ax(iter)
        r = self.b - (self.A.dot(self.x))

        # h(iter) = (r(iter).r(iter))/(r(iter).Ar(iter))
        self.h = r.dot(r) / (r.dot(self.A.dot(r)))

        # omega(iter) = 2*h(iter) / (2 + h(iter))
        self.omega = (2*self.h) / (2+self.h)

        # safety checks of calculated omega (steepest descent might not work perfectly for all systems, hence, this
        # check is necessary to make sure the solver converges)
        if self.omega > 1.9 or self.omega < 0.1:
            self.omega = 1

    def _wolfe_update(self):
        """Update relaxation parameter based on Wolfe conditions from Miyatake et al.
           [https://arxiv.org/abs/1806.09922]
        """

        # f(x) = 1/2 x.Ax - b.x
        f = lambda x: (0.5*(x.dot(self.A.dot(x)))) - (self.b.dot(x))

        # f'(x) = Ax - b
        f_dash = lambda x: (self.A.dot(x)) - (self.b)

        # x(iter) - x(iter-1)
        diff_x = self.x - self.x_old
        # f'(x(iter-1)) . (x(iter) - x(iter-1))
        f_dash_dot_diff_x = f_dash(self.x_old).dot(diff_x)

        # if f(x(iter)) <= f(x(iter-1)) + c1*(f'(x(iter-1)).(x(iter) - x(iter-1)))
        if f(self.x) <= f(self.x_old) + self.c1*f_dash_dot_diff_x:

            # if c2*(f'(x(iter-1)).(x(iter) - x(iter-1))) <= f'(x(iter)).(x(iter) - x(iter-1))
            if self.c2*f_dash_dot_diff_x <= f_dash(self.x).dot(diff_x):
                # h_new = lambda1 * h_old
                self.h = self.lambda1 * self.h
            else:
                # h_new = lambda2 * h_old
                self.h = self.lambda2 * self.h

        else:
            # h_new = rho1 * h_old
            self.h = self.rho1 * self.h

        # omega = 2*h / (2 + h)
        self.omega = (2*self.h) / (2+self.h)

        # safety check for omega
        if self.omega>1.9 or self.omega<0.1:
            self.h = 2
            self.omega = 1

    def _armijo_update(self):
        """Update relaxation parameter based on Armijo conditions method from Miyatake et al.
           [https://arxiv.org/abs/1806.09922]
        """

        # f(x) = 1/2 x.Ax - b.x
        f = lambda x: (0.5 * (x.dot(self.A.dot(x)))) - (self.b.dot(x))

        # f'(x) = Ax - b
        f_dash = lambda x: (self.A.dot(x)) - (self.b)

        # x(iter) - x(iter-1)
        diff_x = self.x - self.x_old
        # f'(x(iter-1)) . (x(iter) - x(iter-1))
        f_dash_dot_diff_x = f_dash(self.x_old).dot(diff_x)

        # if f(x(iter)) <= f(x(iter-1)) + c1*(f'(x(iter-1)).(x(iter) - x(iter-1)))
        if f(self.x) <= f(self.x_old) + self.c1*f_dash_dot_diff_x:
            # h_new = lambda1 * h_old
            print('increasing omega')
            self.h = self.lambda1 * self.h

        else:
            # h_new = rho1 * h_old
            print('decreasing omega')
            self.h = self.rho1 * self.h

        # omega = 2*h / (2 + h)
        self.omega = (2 * self.h) / (2 + self.h)

        # safety check for omega
        if self.omega > 1.9 or self.omega < 0.1:
            self.h = 2
            self.omega = 1

    def _exp_update(self):
        """Update relaxation parameter based on Armijo conditions method from Miyatake et al.
           [https://arxiv.org/abs/1806.09922]
        """

        self.oper = 'incr'

        # f(x) = 1/2 x.Ax - b.x
        f = lambda x: self.b - (self.A.dot(x))

        # if f(x(iter)) <= f(x(iter-1)) + c1*(f'(x(iter-1)).(x(iter) - x(iter-1)))
        if np.linalg.norm(f(self.x)) <= np.linalg.norm(f(self.x_old)):
            self.oper = self.oper
        else:
            if self.oper == 'incr':
                self.oper = 'decr'
            else:
                self.oper = 'incr'

        if self.oper == 'incr':
            # h_new = lambda1 * h_old
            print('increasing omega')
            self.h = self.lambda1 * self.h

        else:
            # h_new = rho1 * h_old
            print('decreasing omega')
            self.h = self.rho1 * self.h

        # omega = 2*h / (2 + h)
        self.omega = (2 * self.h) / (2 + self.h)

        # safety check for omega
        if self.omega > 1.9 or self.omega < 0.1:
            self.h = 2
            self.omega = 1