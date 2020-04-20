import sys
sys.path.append('/home/adhitya/pysolv')

import pysolv
import numpy as np

A = np.arange(12).reshape((4, 3))
b = np.arange(4)

pysolv.solve(A, b, 'Jacobi')
#print(pysolv.__version__)