import functools
from sage.arith.misc import is_prime_power


class BaseAlgorithm:
    def __init__(self, n, m, q=None, w=None):
        """
        Base class for algorithms complexity estimator

        INPUT:

        - ``n`` -- no. of variables
        - ``m`` -- no. of polynomials
        - ``q`` -- order of the field (default: None)
        - ``w`` -- linear algebra constant (default: None)

        TESTS::

            sage: from mpkc.algorithms.base import BaseAlgorithm
            sage: BaseAlgorithm(n=-1, m=5)
            Traceback (most recent call last):
            ...
            ValueError: n must be >= 1
            sage: BaseAlgorithm(n=5, m=0)
            Traceback (most recent call last):
            ...
            ValueError: m must be >= 1
            sage: BaseAlgorithm(n=5, m=10, q=6)
            Traceback (most recent call last):
            ...
            ValueError: q must be a prime power
            sage: BaseAlgorithm(n=5, m=10, w=1)
            Traceback (most recent call last):
            ...
            ValueError: w must be in the range 2 <= w <= 3
        """

        if n < 1:
            raise ValueError("n must be >= 1")

        if m < 1:
            raise ValueError("m must be >= 1")

        if q is not None and not is_prime_power(q):
            raise ValueError("q must be a prime power")

        if w is not None and not 2 <= w <= 3:
            raise ValueError("w must be in the range 2 <= w <= 3")

        self._n = n
        self._m = m
        self._q = q
        self._w = w
        self._optimal_parameters = dict()

    def nvariables(self):
        """
        Return the number of variables

        TESTS::

            sage: from mpkc.algorithms.base import BaseAlgorithm
            sage: BaseAlgorithm(n=10, m=5).nvariables()
            10
        """
        return self._n

    def npolynomials(self):
        """"
        Return the number of polynomials

        TESTS::

            sage: from mpkc.algorithms.base import BaseAlgorithm
            sage: BaseAlgorithm(n=10, m=5).npolynomials()
            5
        """
        return self._m

    def time_complexity(self):
        """
        Return the time complexity of the algorithm

        TESTS::

            sage: from mpkc.algorithms.base import BaseAlgorithm
            sage: BaseAlgorithm(n=10, m=5).time_complexity()
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError

    def memory_complexity(self):
        """
        Return the memory complexity of the algorithm

        TESTS::

            sage: from mpkc.algorithms.base import BaseAlgorithm
            sage: BaseAlgorithm(n=10, m=5).memory_complexity()
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError

    def order_of_the_field(self):
        """
        Return the order of the field

        TESTS::

            sage: from mpkc.algorithms.base import BaseAlgorithm
            sage: BaseAlgorithm(n=10, m=5).order_of_the_field()
            <BLANKLINE>
            sage: BaseAlgorithm(n=10, m=5, q=256).order_of_the_field()
            256
        """
        return self._q

    def linear_algebra_constant(self):
        """
        Return the linear algebra constant

        TESTS::

            sage: from mpkc.algorithms.base import BaseAlgorithm
            sage: BaseAlgorithm(n=10, m=5).linear_algebra_constant()
            <BLANKLINE>
            sage: BaseAlgorithm(n=10, m=5, w=2).linear_algebra_constant()
            2
        """
        return self._w

    def is_overdefined_system(self):
        """
        Return `True` if the system is overdefined

        EXAMPLES::

            sage: from mpkc.algorithms.base import BaseAlgorithm
            sage: BaseAlgorithm(n=5, m=10).is_overdefined_system()
            True
            sage: BaseAlgorithm(n=10, m=5).is_overdefined_system()
            False

        TESTS::

            sage: BaseAlgorithm(n=10, m=10).is_overdefined_system()
            False
        """
        return self.npolynomials() > self.nvariables()

    def is_underdefined_system(self):
        """
        Return `True` if the system is underdefined

        EXAMPLES::

            sage: from mpkc.algorithms.base import BaseAlgorithm
            sage: BaseAlgorithm(n=10, m=5).is_underdefined_system()
            True
            sage: BaseAlgorithm(n=5, m=10).is_underdefined_system()
            False

        TESTS::

            sage: BaseAlgorithm(n=10, m=10).is_underdefined_system()
            False
        """
        return self.nvariables() > self.npolynomials()

    def is_square_system(self):
        """
        Return `True` is the system is square, i.e. there are equal no. of variables and polynomials

        EXAMPLES::

            sage: from mpkc.algorithms.base import BaseAlgorithm
            sage: BaseAlgorithm(n=10, m=10).is_square_system()
            True
            sage: BaseAlgorithm(n=5, m=10).is_square_system()
            False
        """
        return self.nvariables() == self.npolynomials()

    def optimal_parameters(self):
        """
        Return a dictionary of optimal parameters

        EXAMPLES::

            sage: from mpkc.algorithms.base import BaseAlgorithm
            sage: BaseAlgorithm(n=10, m=10).optimal_parameters()
            {}
        """
        return self._optimal_parameters


def optimal_parameter(func):
    """
    Decorator to indicate optimization parameter in BaseAlgorithm

    INPUT:

    - ``f`` -- a method of a BaseAlgoritm subclass
    """
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        name = func.__name__
        self = args[0]

        if name not in self._optimal_parameters:
            self._optimal_parameters[name] = func(*args, **kwargs)
        return self._optimal_parameters[name]
    return wrapper