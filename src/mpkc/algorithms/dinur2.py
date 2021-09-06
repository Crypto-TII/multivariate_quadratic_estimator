"""
Module to compute the time and memory complexity of Dinur's second algorithm

The Dinur's second is an algorithm to solve the MQ problem over F_2

[Din21]  Dinur, I. Cryptanalytic Applications of the Polynomial Method for Solving Multivariate Equation Systems over GF(2).
Springer-Verlag, 2021.
"""
from sage.functions.log import log
from sage.functions.other import ceil
from sage.rings.infinity import Infinity
from ..utils import sum_of_binomial_coefficients
from .base import BaseAlgorithm, optimal_parameter


class DinurSecond(BaseAlgorithm):
    def __init__(self, n, m):
        """
        Construct an instance of Dinur's second estimator

        INPUT:

        - ``n`` -- no. of variables
        - ``m`` -- no. of polynomials

        EXAMPLES::

            sage: from mpkc.algorithms import DinurSecond
            sage: E = DinurSecond(n=10, m=12)
            sage: E
            Dinur's second estimator for the MQ problem
        """
        super().__init__(n=n, m=m)

    def time_specific_parameter(self, n1):
        """
        Return the time complexity of Dinur's second algorithm for the given input parameters

        INPUT:

        - ``n1`` -- no. of variables z

        EXAMPLES::

            sage: from mpkc.algorithms import DinurSecond
            sage: E = DinurSecond(n=10, m=12)
            sage: float(log(E.time_specific_parameter(3),2))
            15.81981143708297
        """
        if n1 < 1:
            raise ValueError('n1 must be >= 1')
        n, m = self.nvariables(), self.npolynomials()
        upper_bound = ((m - 2) // 2) - 1
        if n1 > upper_bound:
            raise ValueError(f"n1 must be <= {upper_bound}")
        return 16 * log(n, 2) * 2 ** n1 * sum_of_binomial_coefficients(n - n1, n1 + 3) + n1 * n * 2 ** (n - n1)

    @optimal_parameter
    def n1(self):
        """
        Return the optimal parameter `n1`

        EXAMPLES::

            sage: from mpkc.algorithms import DinurSecond
            sage: E = DinurSecond(n=10, m=12)
            sage: E.n1()
            4
        """
        n, m = self.nvariables(), self.npolynomials()
        max_n1 = ((m - 2) // 2) - 1
        exp_of_n = ceil(log(n, 2))
        min_time_complexity = Infinity
        best_n1 = None

        if exp_of_n <= max_n1:
            for n1 in range(exp_of_n, max_n1 + 1):
                time_complexity = self.time_specific_parameter(n1)
                if time_complexity < min_time_complexity:
                    best_n1 = n1
                    min_time_complexity = time_complexity

        return best_n1

    def time_complexity(self):
        """
        Return the time complexity of the Dinur's second algorithm

        EXAMPLES::

            sage: from mpkc.algorithms import DinurSecond
            sage: E = DinurSecond(n=10, m=12)
            sage: E.time_complexity().numerical_approx()
            56986.4699066345
        """
        n1 = self.n1()
        return self.time_specific_parameter(n1)

    def memory_complexity(self):
        """
        Return the memory complexity of the Dinur's second algorithm

        EXAMPLES::

            sage: from mpkc.algorithms import DinurSecond
            sage: E = DinurSecond(n=10, m=12)
            sage: E.memory_complexity()
            1280
        """
        n1 = self.n1()
        n = self.nvariables()
        return 4 * (n1 + 1) * 2 ** (n - n1)

    def __repr__(self):
        return f"Dinur's second estimator for the MQ problem"
