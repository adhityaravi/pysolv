import sys
import os
sys.path.append(os.path.join(os.path.dirname(__file__), '..'))

import pysolv
import numpy as np
import scipy.io as io
import time as ti

#print(pysolv.__version__)

# A = [[1, 0, 1], [0, -3, 1], [2, 1, 3]]
# A = np.array(A)
#
# b = [6, 7, 15]
# b = np.array(b)

A = io.mmread('bcsstk04.mtx')
A = A.toarray()

n, *_ = A.shape
b = np.ones(n)

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

x = pysolv.solve(A, b, 'CG')
print('cg')
# print(x)

