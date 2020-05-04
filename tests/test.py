import sys
import os
sys.path.append(os.path.join(os.path.dirname(__file__), '..'))

import pysolv
import numpy as np
import time as ti
import scipy.io as io

#print(pysolv.__version__)

# A = [[1, 0, 1], [0, -3, 1], [2, 1, 3]]
# A = np.array(A)
#
# b = [6, 7, 15]
# b = np.array(b)

A = io.mmread('nos6.mtx')
A = A.toarray()

n, *_ = A.shape
b = np.ones(n)

# x, *_ = pysolv.solve(A, b, 'Jacobi')
# print(x)
#
# x = pysolv.solve(A, b, 'Gauss-Seidel')
# print(x)
#
# x, res, it, time = pysolv.solve(A, b, 'SOR', omega='Steepest Descent')
# print(res, it, time)
#
# x, res, it, time = pysolv.solve(A, b, 'SOR', omega='Wolfe')
# print(res, it, time)

x2, res, it, time = pysolv.solve(A, b, 'SOR', omega='experimental', omega_update_frequency=20)
print(it)

# print(np.allclose(x1, x2, 1e-2, 1e-2))
#x = np.linalg.solve(A, b)

# x2, res, it, time = pysolv.solve(A, b, 'SOR')
# print(it)

#print(A.dot(x2))

# x = pysolv.solve(A, b, 'Gauss Seidel')
# print(x)

# x = pysolv.solve(A, b, 'SSOR')
# print(x)

