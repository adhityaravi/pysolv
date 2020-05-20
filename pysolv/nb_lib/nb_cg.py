"""Conjugate gradient algorithm.

This is a version of the conjugate gradient algorithm which will be compiled with numba.
Check pysolv/krylov_solv/cg.py for details on conjugate gradient algorithm.
"""

# import the necessary packages
import numpy as np
import numba as nb
from numba.pycc import CC


# CC instance
nbc_cg = CC('nbc_cg')

nbc_cg.target_cpu = 'skylake'


# CG solve function
@nbc_cg.export('solve', 'f8[:](f8[:, :], f8[:], f8[:], f8, i8)')
@nb.njit(fastmath=True)
def solve(A, b, x, tol, itermax):

    # initialize the iteration counter
    iter = 0

    # initial residual and search direction
    # dir = res = b - A.x_0
    res = b - mat_vec_mul(A, x)
    dir_ = res.copy()

    while iter < itermax:
        # compute the step size: alpha = res.res / dir.(A.dir)
        rdr = vec_vec_mul(res, res)
        Add = mat_vec_mul(A, dir_)
        ddAdd = vec_vec_mul(dir_, Add)
        alpha = rdr / ddAdd

        # compute the new solution
        x = x + (alpha * dir_)

        # update the residual
        res = res - (alpha * Add)
        # check for convergence
        if np.linalg.norm(res) < tol:
            break

        # compute beta coefficient (from Gram-Schmidt orthogonalization)
        # beta = res_nxt.res_nxt / res.res
        beta = vec_vec_mul(res, res) / rdr

        # update the search direction
        dir_ = res + (beta * dir_)

        # update iteration counter
        iter += 1

    return x


@nb.njit(fastmath=True)
def mat_vec_mul(a, b):

    sum_ = 0

    c = np.empty_like(b)

    for i in range(a.shape[0]):
        for j in range(a.shape[1]):
            sum_ = sum_ + (a[i, j] * b[j])

        c[i] = sum_
        sum_ = 0

    return c


@nb.njit(fastmath=True)
def vec_vec_mul(a, b):

    sum_ = 0

    for i in range(a.shape[0]):
        sum_ += (a[i] * b[i])

    return sum_


# call the compiler
nbc_cg.compile()
