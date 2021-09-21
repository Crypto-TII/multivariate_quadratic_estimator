"""
Module to compute the time and memory complexity of MHT (Miura-Hashimoto-Takagi) algorithm

The MHT is an algorithm to solve the MQ problem when  m * (m + 3) / 2 <= n

[MHT13] Miura, H., Hashimoto, Y., and Takagi, T. Extended algorithm for solving underdefined multivariate quadratic
equations. In Post-Quantum Cryptography, 2013. Springer Berlin Heidelberg.
"""

from sage.arith.misc import is_power_of_two
from .base import BaseAlgorithm


class MHT(BaseAlgorithm):
    def __init__(self, n, m, q, w=2):
        """
        Construct an instance of MHT estimator

        INPUT:

        - ``q`` -- order of the finite field
        - ``n`` -- no. of variables
        - ``m`` -- no. of polynomials
        - ``w`` -- linear algebra constant (2 <= w <= 3) (default: 2)

        EXAMPLES::

            sage: from mpkc.algorithms import MHT
            sage: E = MHT(n=183, m=12, q=4, w=2.8)
            sage: E
            MHT estimator for the MQ problem
        """
        if not m * (m + 3) / 2 <= n:
            raise ValueError(f'The parameter n should be grater than or equal to m * (m + 3) / 2')

        super().__init__(n=n, m=m, q=q, w=w)

    def time_complexity(self):
        """
        Return the time complexity of MHT algorithm

        EXAMPLES::

            sage: from mpkc.algorithms import MHT
            sage: E = MHT(n=183, m=12, q=4, w=2.8)
            sage: float(log(E.time_complexity(),2))
            24.628922047916475
        """
        n, m = self.nvariables(), self.npolynomials()
        w = self.linear_algebra_constant()

        if is_power_of_two(self.order_of_the_field()):
            time = n ** w * m
        else:
            time = 2 ** m * n ** w * m

        return time

    def memory_complexity(self):
        """
        Return the memory complexity of MHT algorithm

        EXAMPLES::

            sage: from mpkc.algorithms import MHT
            sage: E = MHT(n=183, m=12, q=4, w=2.8)
            sage: E.memory_complexity()
            401868
        """
        n, m = self.nvariables(), self.npolynomials()
        return m * n ** 2

    def tilde_o_time(self):
        """
        Return the Ō time complexity

        EXAMPLES::

            sage: from mpkc.algorithms import MHT
            sage: E = MHT(n=183, m=12, q=4, w=2.8)
            sage: E.tilde_o_time()
            1
        """
        return 1

    def __repr__(self):
        return f"MHT estimator for the MQ problem"
