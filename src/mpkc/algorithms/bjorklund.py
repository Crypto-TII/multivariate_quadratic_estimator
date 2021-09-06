"""
Module to compute the time and memory complexity of Bjöklund et al.'s algorithm

The Bjöklund et al.'s is an algorithm to solve the MQ problem over F_2

[BKW19]  Björklund, A., Kaski, P.,  and  Williams, R. Solving  Systemsof Polynomial Equations over GF(2) by a
Parity-Counting Self-Reduction. In International Colloquium on Automata, Languages, and Programming (ICALP 2019),
volume 132 of Leibniz International Proceedings in Informatics (LIPIcs), pages 26:1–26:13, Dagstuhl, Germany, 2019.
Schloss Dagstuhl–Leibniz-Zentrum fuer Informatik.
"""
from sage.rings.infinity import Infinity
from sage.functions.log import log
from sage.functions.other import floor
from .base import BaseAlgorithm, optimal_parameter
from ..utils import sum_of_binomial_coefficients


class Bjorklund(BaseAlgorithm):
    def __init__(self, n, m, nsolutions=1):
        """
        Construct an instance of bjorklund et al.'s estimator

        INPUT:

        - ``n`` -- no. of variables
        - ``m`` -- no. of polynomials
        - ``nsolutions`` -- number of solutions (default: 1)

        EXAMPLES::

            sage: from mpkc.algorithms import Bjorklund
            sage: E = Bjorklund(n=10, m=12)
            sage: E
            Björklund et al.'s estimator for the MQ problem
        """
        super().__init__(n=n, m=m)
        self._nsolutions = nsolutions
        self._k = floor(log(nsolutions + 1, 2))

    @staticmethod
    def _T(n, m, λ):
        if n <= 1:
            return 1
        else:
            l = floor(λ * n)
            T1 = (n + (l + 2) * m * sum_of_binomial_coefficients(n, 2) + (n - l) * 2 ** (n - l))
            s = 48 * n + 1
            return s * sum_of_binomial_coefficients(n - l, l + 4) * (Bjorklund._T(l, l + 2, λ) + T1)

    def nsolutions(self):
        """
        Return the number of solutions

        EXAMPLES::

            sage: from mpkc.algorithms import Bjorklund
            sage: B = Bjorklund(n=10, m=12, nsolutions=3)
            sage: B.nsolutions()
            3
        """
        return self._nsolutions

    @optimal_parameter
    def λ(self):
        """
        Return the optimal λ

        EXAMPLES::

            sage: from mpkc.algorithms import Bjorklund
            sage: E = Bjorklund(n=10, m=12)
            sage: E.λ()
            3/10
        """
        n, m = self.nvariables(), self.npolynomials()
        k = self._k
        min_time_complexity = Infinity
        optimal_λ = None

        for l in range(3, min(m, n - 1)):
            λ_ = l / n
            time_complexity = Bjorklund._T(n, m + k + 2, λ_)
            if time_complexity < min_time_complexity:
                min_time_complexity = time_complexity
                optimal_λ = λ_

        return optimal_λ

    def time_complexity(self):
        """
        Return the time complexity of Bjorklund et al.'s algorithm

        EXAMPLES::

            sage: from mpkc.algorithms import Bjorklund
            sage: E = Bjorklund(n=10, m=12)
            sage: float(log(E.time_complexity(), 2))
            35.48523010807851
        """
        λ = self.λ()
        k = self._k
        n, m = self.nvariables(), self.npolynomials()

        temp = 8 * k * log(n, 2)
        return temp * sum([Bjorklund._T(n - i, m + k + 2, λ) for i in range(1, n)])

    def memory_complexity(self):
        """
        Return the memory complexity of Bjorklund et al.'s algorithm

        EXAMPLES::

            sage: from mpkc.algorithms import Bjorklund
            sage: E = Bjorklund(n=10, m=12)
            sage: float(log(E.memory_complexity(), 2))
            10.89550378006907
        """
        def S(_n, _m, _λ):
            if _n <= 1:
                return 0
            else:
                s = 48 * _n + 1
                l = floor(_λ * _n)
                return S(l, l + 2, _λ) + 2 ** (_n - l) * log(s, 2) + _m * sum_of_binomial_coefficients(_n, 2)

        n, m = self.nvariables(), self.npolynomials()
        λ = self.λ()
        return S(n, m, λ)

    def tilde_o_time(self):
        """
        Return the Ō time complexity of Bjorklund et al.'s algorithm

        EXAMPLES::

            sage: from mpkc.algorithms import Bjorklund
            sage: E = Bjorklund(n=10, m=12)
            sage: float(log(E.tilde_o_time(), 2))
            8.03225
        """
        n = self.nvariables()
        return 2 ** (0.803225 * n)

    def __repr__(self):
        return f"Björklund et al.'s estimator for the MQ problem"
