import sys
sys.path.append('/home/adhitya/pysolv')

import pysolv
import numpy as np

#print(pysolv.__version__)

A = [[2, 1], [5, 7]]
A = np.array(A)

b = [11, 13]
b = np.array(b)

x = pysolv.solve(A, b, 'Jacobi')

print(x)
