# PySolv

## Overview

<p>A collection of linear solvers including stationary iterative solvers like Gauss-Siedel (GS), Jacobi (JM), Successive 
over-relaxation (SOR) and Krylov subspace methods like Conjugate gradient (CG), generalized minimal residual (GMRES) 
etc for Python 3.</p>

<p>The motivation behind PySolv is to provide a complete collection of possible different ways to solve a system of 
linear equations starting from the most basic solver like the Jacobi to complex subspace methods. Therfore, this package 
could be used for benchmarking, testing and educational purposes. Even though, this solver is not aimed at being the 
fastest, to achieve decent run times, fortran bindings are used for most of the solvers</p>

## Getting started

See **Usage** notes on how to use the installed PySolv package. To install PySolv on your local system, the following 
steps are necessary:

### Requirements
PySolv requires a Python version of 3.6 or higher. To check the Python version, the following command can be used:

    python --version
    
Additional python packages required by PySolv are:<br>
1. [NumPy][1] (>=1.18)

In case of building the package from source, the following additional dependencies are necessary:
1. A fortran compiler, for example, [gfortran][5]  (>=7.5)
2. LAPACK and BLAS library from a vendor, for example [OpenBLAS][6]. [OpenBLAS][6] or [ATLAS][7] is recommended over 
[The Netlib][8] for performance

### Installation

Given that, all the requirements are satisfied, the PySolv package can be installed through `pip` using the following 
command:

    pip install pysolv
    
### Building from source

To-Do 
    
### Testing

To-Do.

## Usage

PySolv can be imported into Python using the following command

    >>> import pysolv

* PySolv provides the `solve()` function to solve a system of linear equations. The parameters and the return values of 
  the solve function are as follows: </br>
  
  **Arguments:**
  
    |Argument|Usage|
    |---------|-----|
    |A| Coefficient matrix. Type: 2D numpy array. Note: It must be a 2D square matrix of full rank.|
    |b| RHS vector (dependent variable values). Type: 1D numpy array|

  **Keyword Arguments:**
  
    |Argument|Usage|
    |---------|-----|
    |solver| The type of solver to be used to solve the linear system. Type: str. Possible options are: 'Jacobi', 'Gauss-Seidel', 'SOR', 'SSOR', 'Richardson'. Jacobi solver is used by default.|
    |ITERMAX| The maximum number of iterations that the linear solver is allowed to iterate. Type: int. The default value is set to 3000.|
    |TOL| The tolerance value for stopping the linear solver. Type: float. The default value is set to 1e-6.|
    |omega| Relaxation parameter. Type: float. Significant only for SOR, SSOR and Jacobi solvers. For SOR and SSOR, the default value is taken as 1.7. Jacobi iterations can also be damped with omega. But, by default, the omega value is taken as 1 for Jacobi.|
    
  **Returns:**
  
    |Return|Value|
    |---------|-----|
    |x| Approximated solution to the input linear system. Type: 1D numpy array (matches with the shape of b)|
    |res| Convergence residual of the linear solver. Type: float.|
    |it| Convergence iteration number of the linear solver. Type: int.|
    |time_taken| Time taken by the linear solver for convergence in seconds. Type: float.|

  **Examples**
  
  An example snippet for calling the solve function is as follows:
  
        >>> import numpy as np
        >>> import pysolv
        >>>
        >>> A = np.array([[1, 0, 1], [0, -3, 1], [2, 1, 3]])
        >>> b = np.array([6, 7, 15])
        >>>
        >>> [x, res, it, time_taken] = pysolv.solve(A, b)
        
  In order to specify the type of solver use the following:
  
        >>> [x, res, it, time_taken] = pysolv.solve(A, b, solver='Gauss-Seidel')

  For limiting the linear solver either by number of iterations or by residual or both, following can be done:
  
        >>> [x, res, it, time_taken] = pysolv.solve(A, b, solver='Gauss-Seidel', ITERMAX=1000, TOL=1e-4)
        
  A relaxation parameter for the solver can be introduced as:
  
        >>> [x, res, it, time_taken] = pysolv.solve(A, b, solver='SOR', omega=1.3)
        
## References
1. [Dr. Bojana Rosics's lectures, Institute of Scientific Computing, TU Braunschweig][2]
2. [Iterative solver for sparse linear systems, Yousef Saad][3]
3. [Adaptive SOR methods based on Wolfe conditions, Miyatake et al.][4]


[1]: https://numpy.org/
[2]: https://www.tu-braunschweig.de/en/wire/teaching/previous-terms/winter-2016-17
[3]: https://www-users.cs.umn.edu/~saad/IterMethBook_2ndEd.pdf
[4]: https://link.springer.com/article/10.1007/s11075-019-00748-0
[5]: https://gcc.gnu.org/wiki/GFortran
[6]: https://www.openblas.net/
[7]: http://math-atlas.sourceforge.net/
[8]: https://www.netlib.org/lapack/lug/node11.html
