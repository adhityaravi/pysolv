import sys
import os
sys.path.append(os.path.join(os.path.dirname(__file__), '..'))

import pysolv
import numpy as np
import scipy.io as io
import time as ti
import scipy.io as io
from pysolv.f_lib import f_sdsolve

#print(pysolv.__version__)

# A = [[1, 0, 1], [0, -3, 1], [2, 1, 3]]
# A = np.array(A)
#
# b = [6, 7, 15]
# b = np.array(b)
#
# A = [[0.93355, 0.35572, 0.36875], [0.35572, 0.18057, 0.73134], [0.36875, 0.73134, 0.05757]]
# A = np.array(A)
#
# b = [1, 1, 1]
# b = np.array(b)

A = io.mmread('bcsstk05.mtx')
A = A.toarray()

n, *_ = A.shape
b = np.ones(n)

#print(n)

# st = ti.time()
# D = np.diag(np.diag(A))
# L = np.tril(A, -1)
#
# M = L + D
#
# st = ti.time()
# x = pysolv.solve(M, b, 'Gaussseidel')
# print(M.dot(x))
# print('Done in {} s'.format(ti.time()-st))
# exit()

# x = pysolv.solve(A, b, 'Jacobi')
# print('Jacobi')
# print(x)
#
# x = pysolv.solve(A, b, 'Gauss-Seidel')
# print('Gauss-Seidel')
# print(x)
#
# x = pysolv.solve(A, b, 'SOR', omega=1.2)
# print('SOR')
# print(x)
#
# x = pysolv.solve(A, b, 'SSOR')
# print('SSOR')
# print(x)
# st = ti.time()
# x = pysolv.solve(A, b, 'CG')
# print('Done in {} s'.format(ti.time() - st))
#
# st = ti.time()
# x = pysolv.solve(A, b, 'SOR', itermax=6000)
# print('Done in {} s'.format(ti.time() - st))

st = ti.time()
x = pysolv.solve(A, b, 'SD', itermax=100000)
print('Done in {} s'.format(ti.time() - st))

st = ti.time()
x = pysolv.solve(A, b, 'GaussSeidel', itermax=100000)
print('Done in {} s'.format(ti.time() - st))



# st = ti.time()
# x = pysolv.solve(A, b, 'SOR')
# print('Done in {} s'.format(ti.time() - st))
#
# print(A.dot(x))
# print('cg')
# print(x)

