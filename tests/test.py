import sys
sys.path.append('/home/aravinth/pysolv')

import pysolv
import numpy as np

#print(pysolv.__version__)

A = [[1, 0, 1], [0, -3, 1], [2, 1, 3]]
A = np.array(A)

b = [6, 7, 15]
b = np.array(b)

x = pysolv.solve(A, b, 'Jacobi')
print(x)

x = pysolv.solve(A, b, 'Gaussseidel')
print(x)
