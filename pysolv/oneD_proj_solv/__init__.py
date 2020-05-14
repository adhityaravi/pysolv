# author: Adhitya-Sriram Ravi

"""Simple 1D projection techniques for solving linear system of equations.

Projection techniques for solving a set of linear equations of the form

Ax = b,

involves searching for an approximate solution from a subspace of R{n}. The subspace of R{n} where the solution is
searched for is known as subspace of candidate approximation or search subspace. The dimension of this subspace is m
such that

m < n

This subspace is usually represented as K{m}. Apart from the search subspace, the subspace of constraints or the left
subspace is also important in projection based techniques which is represented as L{m}. Further information on these
topics can be found from the literature: Iterative Methods for Sparse Linear Systems, Yousef Saad.

1D projection techniques are a set of simple projection based solvers where the dimension of the search subspace and
left subspace is limited to 1, i.e, m = 1 --> K{1}, L{1}.
"""

# import the necessary packages
from .sd import SDSolve