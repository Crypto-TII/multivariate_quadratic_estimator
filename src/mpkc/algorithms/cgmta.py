"""
Module to compute the time and memory complexity of CGMT-A algorithm

The CGMT-A is an algorithm to solve the MQ problem

[CGM+02] Courtois, N., Goubin, L., Meier, W.,  and Tacier, J.-D. Solving underdefined systems of multivariate  quadratic
equations. In  D.  Naccache  and  P.  Paillier,  editors,Public  KeyCryptography, pages 211–227, Berlin, Heidelberg,
2002. Springer Berlin Heidelberg.
"""
from sage.all import Integer
from sage.functions.other import sqrt, floor
from sage.misc.functional import numerical_approx
from .base import BaseAlgorithm


class CGMTA(BaseAlgorithm):
    """
    Construct an instance of CGMT-A estimator

    INPUT:

    - ``n`` -- no. of variables
    - ``m`` -- no. of polynomials
    - ``q`` -- order of the finite field

    EXAMPLES::

        sage: from mpkc.algorithms import CGMTA
        sage: E = CGMTA(n=41, m=10, q=3)
        sage: E
        CGMT-A estimator for the MQ problem

    TESTS::

        sage: E.nvariables() == E.nvariables_reduced()
        True
    """
    def __init__(self, n, m, q):
        if not isinstance(q, (int, Integer)):
            raise TypeError("q must be an integer")

        if not m <= n:
            raise ValueError("m must be <= n")

        super().__init__(n=n, m=m, q=q)
        self._k = min(m / 2, floor(sqrt(n / 2 - sqrt(n / 2))))

        if not 2 * self._k ** 2 <= n - 2 * self._k or not m - 2 * self._k < 2 * self._k ** 2:
            raise ValueError(f'The condition 2k^2 <= n - 2k must be satisfied')

        self._n_reduced = n

    def time_complexity(self):
        """
        Return the time complexity of CGMT-A algorithm

        EXAMPLES::

            sage: from mpkc.algorithms import CGMTA
            sage: E = CGMTA(n=41, m=10, q=3)
            sage: E.time_complexity()
            4374.00000000000
        """
        m = self.npolynomials()
        q = self.order_of_the_field()
        k = self._k
        return numerical_approx(2 * q ** (m - k))

    def memory_complexity(self):
        """
        Return the memory complexity of CGMT-A algorithm

        EXAMPLES::

            sage: from mpkc.algorithms import CGMTA
            sage: E = CGMTA(n=41, m=10, q=3)
            sage: E.memory_complexity()
            162.000000000000
        """
        q = self.order_of_the_field()
        k = self._k
        return numerical_approx(2 * k * q ** k)

    def tilde_o_time(self):
        """
        Return the Ō time complexity of of CGMT-A algorithm

        EXAMPLES::

            sage: from mpkc.algorithms import CGMTA
            sage: E = CGMTA(n=41, m=10, q=3)
            sage: E.tilde_o_time()
            2187.00000000000
        """
        m = self.npolynomials()
        q = self.order_of_the_field()
        k = self._k
        return numerical_approx(q ** (m - k))

    def __repr__(self):
        return f"CGMT-A estimator for the MQ problem"
