"""
Module to compute the time and memory complexity of KPG (Kipnis, Patarin, Goubin) algorithm

The KPG is an algorithm to solve a quadratic systems of equations over fields of even characteristic

[KPG99] Kipnis, A., Patarin, J.,  and  Goubin, L. Unbalanced  oil  andvinegar signature schemes. In Advances in
Cryptology —EUROCRYPT99, pages 206–222, Berlin, Heidelberg, 1999. Springer BerlinHeidelberg.
"""
from sage.all import Integer
from sage.arith.misc import is_power_of_two
from .base import BaseAlgorithm


class KPG(BaseAlgorithm):
    def __init__(self, n, m, q, w=2):
        """
        Construct an instance of kpg estimator

        INPUT:

        - ``n`` -- no. of variables
        - ``m`` -- no. of polynomials
        - ``q`` -- order of the finite field
        - ``w`` -- linear algebra constant (2 <= w <= 3) (default: 2)

        EXAMPLES::

            sage: from mpkc.algorithms import KPG
            sage: E = KPG(n=183, m=12, q=4, w=2.8)
            sage: E
            KPG estimator for the MQ problem

        TESTS::

            sage: E.nvariables() == E.nvariables_reduced()
            True
        """
        if not isinstance(q, (int, Integer)):
            raise TypeError("q must be an integer")

        if not is_power_of_two(q):
            raise ValueError("the order of finite field q must be a power of 2")

        if not m * (m + 1) < n:
            raise ValueError(f'The condition m(m + 1) < n must be satisfied')

        super().__init__(n=n, m=m, q=q, w=w)
        self._n_reduced = n

    def time_complexity(self):
        """
        Return the time complexity of kpg algorithm

        EXAMPLES::

            sage: from mpkc.algorithms import KPG
            sage: E = KPG(n=183, m=12, q=4, w=2.8)
            sage: float(log(E.time_complexity(),2))
            24.628922047916475
        """
        n, m = self.nvariables(), self.npolynomials()
        w = self.linear_algebra_constant()
        return m * n ** w

    def memory_complexity(self):
        """
        Return the memory complexity of kpg algorithm

        EXAMPLES::

            sage: from mpkc.algorithms import KPG
            sage: E = KPG(n=183, m=12, q=4, w=2.8)
            sage: E.memory_complexity()
            401868
        """
        n, m = self.nvariables(), self.npolynomials()
        return m * n ** 2

    def tilde_o_time(self):
        """
        Return the Ō time complexity

        EXAMPLES::

            sage: from mpkc.algorithms import KPG
            sage: E = KPG(n=183, m=12, q=4, w=2.8)
            sage: E.tilde_o_time()
            1
        """
        return 1

    def __repr__(self):
        return f"KPG estimator for the MQ problem"
