"""
Module to compute the time and memory complexity of Dinur's first algorithm

The Dinur's first is an algorithm to solve the MQ problem over F_2

[Din21]  Dinur, I. Improved algorithms for solving polynomial systems over GF(2) by multiple parity-counting.
In Proceedings of the 2021 ACM-SIAM Symposium on Discrete Algorithms (SODA), pages 2550–2564.
"""

from sage.functions.log import log
from sage.functions.other import floor
from sage.rings.infinity import Infinity
from .base import BaseAlgorithm, optimal_parameter
from ..utils import sum_of_binomial_coefficients


class DinurFirst(BaseAlgorithm):
    """
    Construct an instance of Dinur's first estimator

    INPUT:

    - ``n`` -- no. of variables
    - ``m`` -- no. of polynomials
    - ``nsolutions`` -- number of solutions (default: 1)

    EXAMPLES::

        sage: from mpkc.algorithms import DinurFirst
        sage: E = DinurFirst(n=10, m=12)
        sage: E
        Dinur's first estimator for the MQ problem
    """
    def __init__(self, n, m, nsolutions=1):
        super().__init__(n=n, m=m, q=2)
        self._nsolutions = nsolutions

        self._k = floor(log(nsolutions + 1, 2))
        self._kappa = None
        self._lambda = None
        self._time_complexity = None
        self._memory_complexity = None

    @optimal_parameter
    def λ(self):
        """
        Return the optimal λ

        EXAMPLES::

            sage: from mpkc.algorithms import DinurFirst
            sage: E = DinurFirst(n=10, m=12)
            sage: E.λ()
            1/9
        """
        if self._lambda is None:
            self._compute_kappa_and_lambda_()
        return self._lambda

    @optimal_parameter
    def κ(self):
        """
        Return the optimal κ

        EXAMPLES::

            sage: from mpkc.algorithms import DinurFirst
            sage: E = DinurFirst(n=10, m=12)
            sage: E.κ()
            2/9
        """
        if self._kappa is None:
            self._compute_kappa_and_lambda_()
        return self._kappa

    def time_complexity(self, **kwargs):
        """
        Return the time complexity of Dinur's first algorithm

        INPUT:

        - ``κ`` -- the parameter `κ` (kappa)
        - ``λ`` -- the parameter `λ`

        If `κ` and `λ` are specified, the function returns the time complexity w.r.t. the given parameter

        EXAMPLES::

            sage: from mpkc.algorithms import DinurFirst
            sage: E = DinurFirst(n=10, m=12)
            sage: float(log(E.time_complexity(), 2))
            26.819919688075288
            sage: float(log(E.time_complexity(κ=0.9, λ=0.9), 2))
            16.73237302312492

        TESTS::

            sage: E0 = DinurFirst(n=15, m=12)
            sage: E1 = DinurFirst(n=17, m=12)
            sage: E0.time_complexity().numerical_approx() == E1.time_complexity().numerical_approx()
            True
        """
        n = self.nvariables_reduced()
        k = self._k
        lambda_ = kwargs.get('λ', self.λ())
        kappa = kwargs.get('κ', self.κ())

        def w(i, kappa):
            return floor((n - i) * (1 - kappa))

        def n1(i, kappa):
            return floor((n - i) * kappa)

        if lambda_ == self.λ() and kappa == self.κ():
            if self._time_complexity is None:
                self._time_complexity = 8 * k * log(n, 2) * \
                                        sum([self._T(n - i, n1(i, kappa), w(i, kappa), lambda_) for i in range(1, n)])
            time_complexity = self._time_complexity
        else:
            time_complexity = 8 * k * log(n, 2) * \
                              sum([self._T(n - i, n1(i, kappa), w(i, kappa), lambda_) for i in range(1, n)])

        return time_complexity

    def memory_complexity(self):
        """
        Return the memory complexity of Dinur's first algorithm

        EXAMPLES::

            sage: from mpkc.algorithms import DinurFirst
            sage: E = DinurFirst(n=10, m=12)
            sage: float(log(E.memory_complexity(), 2))
            15.909893083770042

        TESTS::

            sage: E0 = DinurFirst(n=15, m=12)
            sage: E1 = DinurFirst(n=17, m=12)
            sage: E0.memory_complexity().numerical_approx() == E1.memory_complexity().numerical_approx()
            True
        """
        if self._memory_complexity is None:
            kappa = self.κ()
            n = self.nvariables_reduced()
            self._memory_complexity = (48 * n + 1) * 2 ** (floor((1 - kappa) * n))

        return self._memory_complexity

    def _compute_kappa_and_lambda_(self):
        min_complexity = Infinity
        n, m = self.nvariables_reduced(), self.npolynomials()
        k = self._k
        optimal_kappa = None
        optimal_lambda = None

        for n1 in range(1, min(m + k, (n - 1) // 3)):
            kappa = n1 / (n - 1)
            for n2 in range(1, n1):
                lambda_ = (n1 - n2) / (n - 1)
                n1 = floor((n - 1) * (1 - kappa))
                w = floor((n - 1) * kappa)
                complexity = self._T(n - 1, n1, w, lambda_)
                if complexity < min_complexity:
                    min_complexity = complexity
                    optimal_kappa = kappa
                    optimal_lambda = lambda_

        self._kappa = optimal_kappa
        self._lambda = optimal_lambda

    def tilde_o_time(self):
        """
        Return the Ō time complexity of dinur's first algorithm

        EXAMPLES::

            sage: from mpkc.algorithms import DinurFirst
            sage: E = DinurFirst(n=10, m=12)
            sage: float(log(E.tilde_o_time(), 2))
            6.943
        """
        n = self.nvariables_reduced()
        return 2 ** (0.6943 * n)

    def _T(self, n, n1, w, lambda_):
        t = 48 * n + 1
        n2 = floor(n1 - lambda_ * n)
        l = n2 + 2
        m = self.npolynomials()
        k = self._k

        if n2 <= 0:
            return n * sum_of_binomial_coefficients(n - n1, w) * 2 ** n1
        else:
            temp1 = self._T(n, n2, n2 + 4, lambda_)
            temp2 = n * sum_of_binomial_coefficients(n - n1, w) * 2 ** (n1 - n2)
            temp3 = n * sum_of_binomial_coefficients(n - n2, n2 + 4)
            temp4 = l * (m + k + 2) * sum_of_binomial_coefficients(n, 2)
            return t * (temp1 + temp2 + temp3 + temp4)

    def __repr__(self):
        return f"Dinur's first estimator for the MQ problem"
