"""Data class.

Store all the necessary data from different functions and classes.

This class provides a neat and understandable way for accessing data generated by different functions within the
package. This class does not function as standalone and is majorly inherited by other classes within the package.

Currently data class is setup for handling arrays only of dense numpy type, but in the future will be optimized for
sparse scipy arrays and for handling huge pytables.
"""


class Data:

    # available solvers in pysolv
    SOLVERS = ['jacobi', 'gaussseidel', 'sor', 'ssor', 'cg']

    # available pre-conditioners for CG in pysolv
    PRECONDITIONER_FLAVOR = ['jacobi', 'gaussseidel', 'sor', 'ssor']

    # available techniques to adapt relaxation parameter (significant only for SOR)
    ADAPTIVE_OMEGA_FLAVOR = {'steepestdescent': 1, 'armijo': 2, 'wolfe': 3}
    # parameters for Wolfe and Armijo condition for relaxation parameter update. These values are chosen based on paper
    # from Miyatake et al. (can be over-ridden by user)
    c1 = 0.89
    c2 = 0.95
    lambda1 = 1.15
    lambda2 = 1.4
    rho1 = 0.85

    # general solver settings (can be over-ridden by user)
    itermax = 3000
    tol = 1e-6

    # attributes acquired from other classes by the Data class. Declared as private and for internal use only by the
    # Data class
    __acquired_attrs = []

    def __init__(self):
        """Initialize the class that inherits Data class with all the attributes of the Data class
        """

        self._get_attr()

    @staticmethod
    def add_data(name, data):
        """Add any required data to the data class to be shared with other functions/classes

        Parameters:
            name (string): name with which the data can be called
            data (dtype): data to be stored
        """

        if name in ['SOLVERS', 'ADAPTIVE_OMEGA', 'ADAPTIVE_OMEGA_FLAVOR']:
            pass
        else:
            setattr(Data, name, data)
            Data.__acquired_attrs.append(name)

    def _get_attr(self):
        """Pass the __dict__ of the Data class to any child class of Data class
        """

        for name, data in Data.__dict__.items():
            try:
                setattr(self, name, data)
            except (TypeError, AttributeError):
                continue

    @staticmethod
    def __reset_data():
        """Flush all the attributes from the Data class. This method is intended to be a private static method. Call
           to this function from outside this class is not recommended.
        """

        for attr in Data.__acquired_attrs:
            try:
                delattr(Data, attr)
            except (TypeError, AttributeError):
                continue
