"""
Module to compute the time and memory complexity of Dinur's second algorithm

The Dinur's second is an algorithm to solve the MQ problem over F_2

[Din21]  Dinur, I. Cryptanalytic Applications of the Polynomial Method for Solving Multivariate Equation Systems over GF(2).
Springer-Verlag, 2021.
"""

from math import inf
from sage.functions.log import log
from sage.functions.other import ceil
from mpkc.utils import sum_of_binomial_coefficients


class DinurSecond(object):
    def __init__(self, n, m):
        """
        Construct an instance of Dinur's second estimator

        INPUT:

        - ``n`` -- no. of variables
        - ``m`` -- no. of polynomials

        EXAMPLES::

            sage: from mpkc.algorithms.dinur2.complexity import DinurSecond
            sage: E = DinurSecond(n=10, m=12)
            sage: E
            Dinur's second estimator for the MQ problem
        """
        self._n = n
        self._m = m
        self._optimal_parameter = None

    def time_specific_parameters(self, n1):
        """
        Return the time complexity of Dinur's second algorithm for the given input parameters

        INPUT:

        - ``n1`` -- no. of variables z

        EXAMPLES::

            sage: from mpkc.algorithms.dinur2.complexity import DinurSecond
            sage: E = DinurSecond(n=10, m=12)
            sage: float(log(E.time_specific_parameters(3),2))
            15.81981143708297
        """
        if n1 < 1:
            raise ValueError(f'The paramereter n1 should be at least 1')
        if n1 > int((self._m - 2) / 2) - 1:
            raise ValueError(f'The paramereter n1 should be smaller than {int((self._m - 2) / 2)}')
        return 16 * log(self._n, 2) * 2 ** n1 * sum_of_binomial_coefficients(self._n - n1, n1 + 3) + \
               n1 * self._n * 2 ** (self._n - n1)

    def optimal_parameter(self):
        """
        Return the parameter optimizing the time complexity of dinur's second algorithm

        EXAMPLES::

            sage: from mpkc.algorithms.dinur2.complexity import DinurSecond
            sage: E = DinurSecond(n=10, m=12)
            sage: E.optimal_parameter()
            4
        """

        if self._optimal_parameter is not None:
            return self._optimal_parameter

        time = inf
        max_n1 = int((self._m - 2) / 2) - 1
        if ceil(log(self._n, 2)) <= max_n1:
            for n1 in range(ceil(log(self._n, 2)), max_n1 + 1):
                temp_time = self.time_specific_parameters(n1)
                if temp_time < time:
                    self._optimal_parameter = n1
                    time = temp_time
        return self._optimal_parameter

    def time(self):
        """
        Return the time complexity of dinur's second algorithm

        EXAMPLES::

            sage: from mpkc.algorithms.dinur2.complexity import DinurSecond
            sage: E = DinurSecond(n=10, m=12)
            sage: E.time().numerical_approx()
            56986.4699066345
        """
        n1 = self.optimal_parameter()
        return self.time_specific_parameters(n1)

    def memory(self):
        """
        Return the memory complexity of dinur's second algorithm

        EXAMPLES::

            sage: from mpkc.algorithms.dinur2.complexity import DinurSecond
            sage: E = DinurSecond(n=10, m=12)
            sage: E.memory()
            1280
        """
        n1 = self.optimal_parameter()
        return 4 * (n1 + 1) * 2 ** (self._n - n1)

    def tilde_o_time(self):
        """
        Return the Ō time complexity of Dinur's second algorithm

        EXAMPLES::

            sage: from mpkc.algorithms.dinur2.complexity import DinurSecond
            sage: E = DinurSecond(n=10, m=12)
            sage: E.tilde_o_time()
            283.68541077888506
        """
        return 2 ** ((1 - 1./(2.7*2)) * self._n)

    def __repr__(self):
        return f"Dinur's second estimator for the MQ problem"
