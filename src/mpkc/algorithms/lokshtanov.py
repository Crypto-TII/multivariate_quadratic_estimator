"""
Module to compute the time and memory complexity of Losktanov et al.'s algorithm

The Losktanov et al.'s is an algorithm to solve the MQ problem

[LPT+17]  Lokshtanov, D.,  Paturi, R., Tamaki, S., Williams, R., and Yu, H. Beating brute force for systems of
polynomial equation sover finite fields. In Proceedings of the Twenty-Eighth Annual ACM-SIAM Symposium on Discrete
Algorithms, SODA ’17, page 2190–2202, USA, 2017. Society for Industrial and Applied Mathematics.
"""
from sage.all import Integer
from sage.arith.misc import is_power_of_two
from sage.functions.log import log
from sage.functions.other import floor
from sage.rings.infinity import Infinity
from sage.rings.finite_rings.finite_field_constructor import GF
from .base import BaseAlgorithm, optimal_parameter
from ..series.nmonomial import NMonomialSeries


class Lokshtanov(BaseAlgorithm):
    def __init__(self, n, m, q):
        """
        Construct an instance of Lokshtanov et al.'s estimator

        INPUT:

        - ``n`` -- no. of variables
        - ``m`` -- no. of polynomials
        - ``q`` -- order of the finite field

        EXAMPLES::

            sage: from mpkc.algorithms import Lokshtanov
            sage: E = Lokshtanov(n=10, m=12, q=9)
            sage: E
            Lokshtanov et al.'s estimator for the MQ problem
        """
        if not isinstance(q, (int, Integer)):
            raise TypeError("q must be an integer")

        super().__init__(n=n, m=m, q=q)
        self._time_complexity = None
        self._memory_complexity = None

    @optimal_parameter
    def δ(self):
        """
        Return the optimal δ for Lokshtanov et al.'s algorithm

        EXAMPLES::

            sage: from mpkc.algorithms import Lokshtanov
            sage: E = Lokshtanov(n=10, m=12, q=9)
            sage: E.δ()
            1/10
        """
        min_complexity = Infinity
        optimal_δ = None
        n, m = self.nvariables(), self.npolynomials()

        for np in range(1, min(m - 2, n)):
            δ = np / n
            time_complexity = self._C(n - 1, δ)
            if time_complexity < min_complexity:
                min_complexity = time_complexity
                optimal_δ = δ

        return optimal_δ

    def time_complexity(self, **kwargs):
        """
        Return the time complexity of lokshtanov et al.'s algorithm

        INPUT:

        - ``δ`` -- the parameter `δ`

        If `δ` is specified, the function returns the time complexity w.r.t. the given parameter

        EXAMPLES::

            sage: from mpkc.algorithms import Lokshtanov
            sage: E = Lokshtanov(n=10, m=12, q=9)
            sage: float(log(E.time_complexity(), 2))
            212.576588724275
            sage: float(log(E.time_complexity(δ=2/10), 2))
            214.16804105519708
        """
        q = self.order_of_the_field()
        n = self.nvariables()
        δ = kwargs.get('δ', self.δ())

        if δ is None:
            return Infinity
        else:
            if not 0 < δ < 1:
                raise ValueError("δ must be in the range 0 < δ < 1")

            if δ == self.δ():
                if self._time_complexity is None:
                    self._time_complexity = 100 * log(q, 2) * (q - 1) * sum([self._C(n - i, δ) for i in range(1, n)])
                    time_complexity = self._time_complexity
                else:
                    time_complexity = self._time_complexity
            else:
                time_complexity = 100 * log(q, 2) * (q - 1) * sum([self._C(n - i, δ) for i in range(1, n)])

        return time_complexity

    def memory_complexity(self):
        """
        Return the memory complexity of Lokshtanov et al.'s algorithm

        EXAMPLES::

            sage: from mpkc.algorithms import Lokshtanov
            sage: E = Lokshtanov(n=10, m=12, q=9)
            sage: float(log(E.memory_complexity(), 2))
            30.622995719758727
        """
        if self._memory_complexity is None:
            δ = self.δ()
            if δ is None:
                return Infinity

            n = self.nvariables()
            np = floor(n * δ)
            q = self.order_of_the_field()
            resulting_degree = 2 * (q - 1) * (np + 2)
            M = NMonomialSeries(n=n - np, q=q, max_prec=resulting_degree + 1).nmonomials_up_to_degree(resulting_degree)
            self._memory_complexity =  M + log(n, 2) * q ** (n - np)
        return self._memory_complexity

    def tilde_o_time(self):
        """
        Return the Ō time complexity of Lokshtanov et al.'s algorithm

        EXAMPLES::

            sage: from mpkc.algorithms import Lokshtanov
            sage: E = Lokshtanov(n=10, m=12, q=9)
            sage: float(log(E.tilde_o_time(), 2))
            31.62000188938707
        """
        e = 2.718
        q = self.order_of_the_field()
        n = self.nvariables()
        if q == 2:
            time = q ** (0.8765 * n)
        elif is_power_of_two(q):
            time = q ** (0.9 * n)
        elif log(GF(q).characteristic(), 2) < 8 * e:
            time = q ** (0.9975 * n)
        else:
            d = GF(q).degree()
            time = q ** n * (log(q, 2) / (2 * e * d))

        return time

    def _C(self, n, delta):
        q = self.order_of_the_field()
        np = floor(delta * n)
        resulting_degree = 2 * (q - 1) * (np + 2)
        M = NMonomialSeries(n=n - np, q=q, max_prec=resulting_degree + 1).nmonomials_up_to_degree(resulting_degree)
        return n * (q ** (n - np) + M * q ** np * n ** (6 * q))

    def __repr__(self):
        return f"Lokshtanov et al.'s estimator for the MQ problem"
