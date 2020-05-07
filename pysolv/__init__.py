# author: Adhitya-Sriram Ravi

# resolve path conflicts (only for testing purpose)
# ToDo: remove before packaging
import sys
import os
sys.path.append(os.path.dirname(__file__))

# set the version number
__version__ = "0.0.1"

# import the necessary packages
from .solve import solve
