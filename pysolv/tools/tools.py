"""Tool functions.

A set of tool functions to assist the various function/classes of the pysolv package.
"""

# import the necessary packages
from data import Data
import re
from krylov_solv import *
from st_iter_solv import *
from exceptions import *
import numpy as np


def check_solver(solver):
    """Check if the selected solver scheme is available in pysolv
    """

    # remove unwanted characters from the solver string and make it case insensitive
    solver = re.sub('\W+', '', solver).casefold()

    # if the requested solver scheme is not available in pysolv, raise an error
    if solver not in Data.SOLVERS:
        msg = 'Requested solver scheme not available in pysolv. Please select from the available schemes: {}'\
              .format(Data.SOLVERS)
        raise ValueError(msg)

    # else, add it to the Data class
    else:
        Data.add_data('solver', solver)


def add_kwargs(kwargs):
    """Add the keyword arguments, passed to the solve function, to the Data class
    """

    for key, value in kwargs.items():
        Data.add_data(key, value)


def ext_lsprops():
    """Extract various properties of the linear system like the size, etc. and store it in the Data class
    """

    # sanity checks
    # check if the input values are of type numpy array
    if not isinstance(Data.A, np.ndarray):
        msg = 'Expected {}. Got {}'.format(type(np.ndarray(0)), type(Data.A))
        raise TypeError(msg)
    if not isinstance(Data.b, np.ndarray):
        msg = 'Expected {}. Got {}'.format(type(np.ndarray(0)), type(Data.b))
        raise TypeError(msg)

    # shape of the coefficient matrix
    try:
        m, n = Data.A.shape
        if m == 1:
            msg = 'The minimum size of A matrix is 2x2'
            raise ValueError(msg)
        elif n == 1:
            msg = 'The minimum size of A matrix is 2x2'
            raise ValueError(msg)
    except ValueError:
        msg = 'The minimum size of A matrix is 2x2'
        raise ValueError(msg)

    # check if the dimensions of the A matrix and the b vector match
    if m != Data.b.size:
        msg = 'Input array dimensions do not match'
        raise DimensionMismatchException(msg)

    # add the dimensions of the linear system to the Data class
    Data.add_data('m', m)
    Data.add_data('n', n)


def call_solver():
    """Call the requested solver to solve the linear system of equations
    """

    if Data.solver == 'jacobi':
        JacobiSolve()

    elif Data.solver == 'gaussseidel':
        GaussSeidelSolve()

    elif Data.solver == 'sor':
        SORSolve()

    elif Data.solver == 'ssor':
        SSORSolve()

    elif Data.solver == 'cg':
        CGSolve()

    else:
        print('wubba lubba dub dub')
