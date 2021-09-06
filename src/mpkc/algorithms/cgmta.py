"""
Module to compute the time and memory complexity of CGMT-A algorithm

The CGMT-A is an algorithm to solve the MQ problem

[CGM+02] Courtois, N., Goubin, L., Meier, W.,  and Tacier, J.-D. Solving underdefined systems of multivariate  quadratic
equations. In  D.  Naccache  and  P.  Paillier,  editors,Public  KeyCryptography, pages 211–227, Berlin, Heidelberg,
2002. Springer Berlin Heidelberg.
"""

from sage.functions.other import sqrt
from sage.misc.functional import numerical_approx
from .base import BaseAlgorithm


class CGMTA(BaseAlgorithm):
    def __init__(self, n, m, q):
        """
        Construct an instance of CGMT-A estimator

        INPUT:

        - ``q`` -- order of the finite field
        - ``n`` -- no. of variables
        - ``m`` -- no. of polynomials

        EXAMPLES::

            sage: from mpkc.algorithms import CGMTA
            sage: E = CGMTA(n=12, m=10, q=3)
            sage: E
            CGMT-A estimator for the MQ problem
        """
        if not m <= n:
            raise ValueError("m must be <= n")

        super().__init__(n=n, m=m, q=q)
        self._k = min(m / 2, sqrt(n / 2 - sqrt(n / 2)))

    def time_complexity(self):
        """
        Return the time complexity of CGMT-A algorithm

        EXAMPLES::

            sage: from mpkc.algorithms import CGMTA
            sage: E = CGMTA(n=12, m=10, q=3)
            sage: E.time_complexity()
            14900.9039071367
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
            sage: E = CGMTA(n=12, m=10, q=3)
            sage: E.memory_complexity()
            29.8679427555961
        """
        q = self.order_of_the_field()
        k = self._k
        return numerical_approx(2 * k * q ** k)

    def tilde_o_time(self):
        """
        Return the Ō time complexity of of CGMT-A algorithm

        EXAMPLES::

            sage: from mpkc.algorithms import CGMTA
            sage: E = CGMTA(n=12, m=10, q=3)
            sage: E.tilde_o_time()
            7450.45195356836
        """
        m = self.npolynomials()
        q = self.order_of_the_field()
        k = self._k
        return numerical_approx(q ** (m - k))

    def __repr__(self):
        return f"CGMT-A estimator for the MQ problem"
