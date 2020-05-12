"""Conjugate Gradient linear solver.

"""

# import the necessary packages
from data import *
import numpy as np
import time as ti


class CGSolve(Data):

    def __init__(self):
        """Initialize the class with the necessary inheritance and solve the linear system
           ToDO: Deprecate Conjugate Gradient _solve() and rather call SOR with omega=1
        """

        # data inheritance
        Data.__init__(self)

        # initialize solution vector
        self._init_x0()

        # call the Conjugate Gradient solver
        self._solve()

    def _init_x0(self):
        """Check if an initial value exists for the solution vector or initialize the solution vector"""

        try:
            self.x0
        except AttributeError:
            self.x0 = np.ones(self.n)

        # initialize the solution vectors
        self.x = np.empty(self.n)   # solution vector at the i'th iteration
        self.x_old = self.x0  # solution vector at the i-1'th iteration
        self.x = self.x0

    def _solve(self):
        """Serial python implementation of Conjugate Gradient solver."""
        try:
            if self.pc == 'jacobi':
                pc = 1
            elif self.pc == 'gaussseidel':
                pc = 2
            elif self.pc == 'sor':
                pc = 3
            elif self.pc == 'ssor':
                pc = 4
        except AttributeError:
            pc = 0

        if self.m != self.n:
            msg = 'pysolv does not include a Conjugate Gradient solver for non-square systems yet'
            raise ValueError(msg)

        # Conjugate Gradient iterations
        else:
            # time the solver
            start_time = ti.time()

            # Initialize the iteration counter and solution vector
            iter = 0

            r = self.b - np.dot(self.A, self.x)
            z = np.empty(self.n)
            if pc == 1:
                # Jacobi preconditioner C = D
                z = np.multiply(1/np.diag(self.A), r)
            elif pc == 2:
                # Guass-seidel preconditioner C = (L + D)
                for i in range(0, self.n):
                    sum = 0
                    for j in range(0, i):
                        sum += self.A[i, j]*z[j]
                    z[i] = (1/self.A[i, i])*(r[i] - sum)
                print(z)
            elif pc == 3:
                # SOR preconditioner C = (1/omega)*(D + omega*L)
                # To be implemented
                for i in range(0, self.n):
                    sum = 0
                    for j in range(0, i):
                        sum += self.A[i, j]*z[j]
                    z[i] = (self.omega/self.A[i, i])*(r[i] - self.omega*sum)
            elif pc == 4:
                # SSOR preconditioner C = (D/omega + L)*(omega/(2-omega))*(D^-1)(D/omega+L^T)
                # To be implemented
                msg = 'SSOR Preconditioner yet to be implemented'
                raise ValueError(msg)
            elif pc == 0:
                # No preconditioner
                z = r.copy()
            d = z.copy()
            alpha_num = np.dot(z, r)
            while iter < self.ITERMAX:
                print(iter)

                ad = np.dot(self.A, d)
                alpha = alpha_num/np.dot(ad, d)
                self.x = self.x + alpha*d
                r = r - alpha*ad
                if pc == 1:
                    # Jacobi preconditioner C = D
                    z = np.multiply(1/np.diag(self.A), r)
                elif pc == 2:
                    # Guass-seidel preconditioner C = (L + D)
                    for i in range(0, self.n):
                        sum = 0
                        for j in range(0, i):
                            sum += self.A[i, j]*z[j]
                        z[i] = (1/self.A[i, i])*(r[i] - sum)
                elif pc == 3:
                    # SOR preconditioner C = (1/omega)*(D + omega*L)
                    for i in range(0, self.n):
                        sum = 0
                        for j in range(0, i):
                            sum += self.A[i, j]*z[j]
                        z[i] = (self.omega/self.A[i, i])*(r[i] - self.omega*sum)
                elif pc == 4:
                    # SSOR preconditioner C = (D/omega + L)*(omega/(2-omega))*(D^-1)(D/omega+L^T)
                    # To be implemented
                    msg = 'SSOR Preconditioner yet to be implemented'
                    raise ValueError(msg)
                elif pc == 0:
                    # No preconditioner
                    z = r.copy()
                alpha_num_next = np.dot(z, r)
                beta = alpha_num_next/alpha_num
                d = z + beta*d
                alpha_num = alpha_num_next

                # stopping criteria: (||x(iter) - x(iter-1)|| / ||x(iter)||) < TOL
                # ToDo: Implement a pysolv native function to compute the residual
                res = np.linalg.norm(np.dot(self.A, self.x) - self.b)

                self.x_old = self.x.copy()
                if res < self.TOL:
                    break
                # Update the iteration counter and x(iter-1)
                iter += 1

            # time taken for convergence
            time_taken = ti.time() - start_time
            print(np.dot(self.A, self.x))
            # add the solution to the Data class
            Data.add_data('x', self.x)
            Data.add_data('time_taken', time_taken)
            Data.add_data('iterations', iter - 1)
            Data.add_data('residual', res)