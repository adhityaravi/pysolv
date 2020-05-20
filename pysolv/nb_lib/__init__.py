# author: Adhitya-Sriram Ravi

"""Numba compiled python codes for performing various pysolv related operations.

The numba AOT compilation is used to test if any speed up over the currently used f2py extensions.

Possible advantages:
* Native python code with Fortran/C like run times
* If auto-parallelization (which is currently an issue in the numba repo) is forced for the numba AOT compiler like in
  the JIT compiler, further speed up is achievable with minimum effort for parallelization

Possible disadvantages:
* Might not be best suited for packaging and distributing
* Currently not able to specify the target folder for the compiled code (uses the source code folder by default)
* Relatively a new concept, arises support and stability doubts.

Because of the possible disadvantages, this enhancement will remain experimental for the time being.
"""

# import the necessary packages
