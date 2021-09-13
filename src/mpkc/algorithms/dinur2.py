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
        self._n1 = None
        self._time_complexity = None
        self._memory_complexity = None

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
        if self._n1 is None:
            self._compute_time_complexity_()
        return self._n1

    def time_complexity(self):
        """
        Return the time complexity of the Dinur's second algorithm

        EXAMPLES::

            sage: from mpkc.algorithms import DinurSecond
            sage: E = DinurSecond(n=10, m=12)
            sage: E.time_complexity().numerical_approx()
            57434.4699066345
        """
        if self._time_complexity is None:
            self._compute_time_complexity_()
        return self._time_complexity

    def memory_complexity(self):
        """
        Return the memory complexity of the Dinur's second algorithm

        EXAMPLES::

            sage: from mpkc.algorithms import DinurSecond
            sage: E = DinurSecond(n=10, m=12)
            sage: E.memory_complexity()
            2560
        """
        if self._memory_complexity is not None:
            return self._memory_complexity

        n1 = self.n1()
        n = self.nvariables()
        self._memory_complexity = 8 * (n1 + 1) * sum_of_binomial_coefficients(n - n1, n1 + 3)
        return self._memory_complexity

    def _compute_time_complexity_(self):
        n, m = self.nvariables(), self.npolynomials()
        max_n1 = ((m - 2) // 2) - 1
        min_time_complexity = Infinity
        optimal_n1 = None

        for n1 in range(1, max_n1 + 1):
            time_complexity = 16 * log(n, 2) * 2 ** n1 * sum_of_binomial_coefficients(n - n1, n1 + 3) +\
                              n1 * n * 2 ** (n - n1) + \
                              2 ** (self._n - 2 * n1 + 1) * sum_of_binomial_coefficients(self._n, 2)
            if time_complexity < min_time_complexity:
                optimal_n1 = n1
                min_time_complexity = time_complexity

        self._n1 = optimal_n1
        self._time_complexity = min_time_complexity

    def __repr__(self):
        return f"Dinur's second estimator for the MQ problem"
