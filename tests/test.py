import sys
import os
sys.path.append(os.path.join(os.path.dirname(__file__), '..'))

import pysolv
import numpy as np
import time as ti
import scipy.io as io
from pysolv.f_lib import f_jacobisolve

#print(pysolv.__version__)

# A = io.mmread('bcsstk05.mtx')
# A = A.toarray()
#
# n, *_ = A.shape
# b = np.ones(n)
#
A = [[1, 0, 1], [0, -3, 1], [2, 1, 3]]
A = np.array(A)

b = [6, 7, 15]
b = np.array(b)

x0 = np.ones(3)
x = np.empty(3)

#f_jacobisolve.solve(A, b, x0, x, 1000, 1e-6, 1)
x = pysolv.solve(A, b)
print(b)
print(A.dot(x))

#
# x = np.empty(3)
#
# f_sorsolve.ssolve(x)
# print(f_sorsolve.iter_)
# print(b)
# print(A.dot(x))
# print(f_sorsolve.res_)

# A = io.mmread('bcsstk04.mtx')
# A = A.toarray()
#
# n, *_ = A.shape
# b = np.ones(n)
#
# x0 = np.ones(n)
#
# st = ti.time()
# x = pysolv.solve(A, b, 'jacobi', ITERMAX=3000, TOL=1e-6)
# print(ti.time()-st)
# print(A.dot(x))

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

# x2, res, it, time = pysolv.solve(A, b, 'SOR', omega=1.85, ITERMAX=504)
# print(it, time)
# print(A.dot(x2))
#print(it, time)
#print(A.dot(x2))

# print(np.allclose(x1, x2, 1e-2, 1e-2))
#x = np.linalg.solve(A, b)

# x2, res, it, time = pysolv.solve(A, b, 'SOR')
# print(it)

#print(A.dot(x2))

# x = pysolv.solve(A, b, 'Gauss Seidel')
# print(x)

# x = pysolv.solve(A, b, 'SSOR')
# print(x)

