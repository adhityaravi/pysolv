import sys
import os
sys.path.append(os.path.join(os.path.dirname(__file__), '..'))

import pysolv
import numpy as np

#print(pysolv.__version__)

A = [[1, 0, 1], [0, -3, 1], [2, 1, 3]]
A = np.array(A)

b = [6, 7, 15]
b = np.array(b)

x = pysolv.solve(A, b, 'Jacobi')
print(x)

x = pysolv.solve(A, b, 'Gauss-Seidel')
print(x)

x = pysolv.solve(A, b, 'SOR', ITERMAX=1000)
print(x)
