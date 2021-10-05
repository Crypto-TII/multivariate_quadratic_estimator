"""
Module to compute the time and memory complexity of Crossbred algorithm

The Crossbred is an algorithm to solve the MQ problem

[JV18]  Joux, A., and Vitse, V.  A crossbred algorithm for solving boolean polynomial  systems.
In  Jerzy  Kaczorowski,  Josef  Pieprzyk,  and  JacekPomyka la,  editors,Number-Theoretic  Methods  in  Cryptology,
pages  3–21,Cham, 2018. Springer International Publishing.

[Dua20] Duarte, J. D. On the complexity of the crossbred algorithm.
Cryptologye Print Archive, Report 2020/1058, 2020.https://eprint.iacr.org/2020/1058.
"""

from sage.all import Integer
from sage.functions.log import log
from sage.rings.all import QQ
from sage.rings.infinity import Infinity
from sage.rings.power_series_ring import PowerSeriesRing
from .base import BaseAlgorithm, optimal_parameter
from ..series.hilbert import HilbertSeries
from ..series.nmonomial import NMonomialSeries
from ..utils import nmonomials_up_to_degree


class Crossbred(BaseAlgorithm):
    def __init__(self, n, m, q, w=2, max_D=10):
        """
        Construct an instance of crossbred estimator

        INPUT:

        - ``n`` -- no. of variables
        - ``m`` -- no. of polynomials
        - ``q`` -- order of the finite field
        - ``w`` -- linear algebra constant (2 <= w <= 3) (default: 2)
        - ``max_D`` -- upper bound to the parameter D (default: 10)

        EXAMPLES::

            sage: from mpkc.algorithms import Crossbred
            sage: E = Crossbred(n=10, m=12, q=5)
            sage: E
            Crossbred estimator for the MQ problem
        """
        if not isinstance(q, (int, Integer)):
            raise TypeError("q must be an integer")

        super().__init__(n=n, m=m, q=q, w=w)
        self._max_D = max_D
        self._k = None
        self._D = None
        self._d = None
        self._time_complexity = None
        self._memory_complexity = None

    def max_D(self):
        """
        Return the upper bound of the degree of the initial Macaulay matrix

        EXAMPLES::

            sage: from mpkc.algorithms import Crossbred
            sage: E = Crossbred(n=10, m=12, q=5, max_D=9)
            sage: E.max_D()
            9
        """
        return self._max_D

    def ncols_in_preprocessing_step(self, k, D, d):
        """
        Return the number of columns involve in the preprocessing step

        INPUT:

        - ``k`` -- no. variables in the resulting system
        - ``D`` -- degree of the initial Macaulay matrix
        - ``d`` -- degree resulting Macaulay matrix

        EXAMPLES::

            sage: from mpkc.algorithms import Crossbred
            sage: E = Crossbred(n=10, m=12, q=5)
            sage: E.ncols_in_preprocessing_step(4, 6, 3)
            297
        """
        if not d < D:
            raise ValueError("d must be smaller than D")

        n = self.nvariables()
        q = self.order_of_the_field()

        nms0 = NMonomialSeries(n=k, q=q, max_prec=D+1)
        nms1 = NMonomialSeries(n=n-k, q=q, max_prec=D+1)

        ncols = 0
        for dk in range(d + 1, D):
            ncols += sum([nms0.nmonomials_of_degree(dk) * nms1.nmonomials_of_degree(dp) for dp in range(D - dk)])

        return ncols

    def ncols_in_linearization_step(self, k, d):
        """
        Return the number of columns involve in the linearization step

        INPUT:

        - ``k`` -- no. variables in the resulting system
        - ``d`` -- degree resulting Macaulay matrix

        EXAMPLES::

            sage: from mpkc.algorithms import Crossbred
            sage: E = Crossbred(n=10, m=12, q=5)
            sage: E.ncols_in_linearization_step(4, 3)
            35
        """
        return nmonomials_up_to_degree(d, k, q=self.order_of_the_field())

    @optimal_parameter
    def k(self):
        """
        Return the optimal `k`, i.e. no. of variables in the resulting system

        EXAMPLES::

            sage: from mpkc.algorithms import Crossbred
            sage: E = Crossbred(n=10, m=12, q=5)
            sage: E.k()
            5
        """
        if self._k is None:
            _ = self.time_complexity()
        return self._k

    @optimal_parameter
    def D(self):
        """
        Return the optimal `D`, i.e. degree of the initial Macaulay matrix

        EXAMPLES::

            sage: from mpkc.algorithms import Crossbred
            sage: E = Crossbred(n=10, m=12, q=5)
            sage: E.D()
            3
        """
        if self._D is None:
            _ = self.time_complexity()
        return self._D

    @optimal_parameter
    def d(self):
        """
        Return the optimal `d`, i.e. degree resulting Macaulay matrix

        EXAMPLES::

            sage: from mpkc.algorithms import Crossbred
            sage: E = Crossbred(n=10, m=12, q=5)
            sage: E.d()
            1
        """
        if self._d is None:
            _ = self.time_complexity()
        return self._d

    def admissible_parameter_series(self, k):
        """
        Return a the series S_k of admissible parameters

        INPUT:

        - ``k`` -- no. variables in the resulting system

        EXAMPLES::

            sage: from mpkc.algorithms import Crossbred
            sage: E = Crossbred(n=10, m=12, q=5, max_D=2)
            sage: E.admissible_parameter_series(2)
            -1 - 3*x - 3*y - 10*x^2 - 3*x*y + 6*y^2 + O(x, y)^3
        """
        n, m = self.nvariables(), self.npolynomials()
        q = self.order_of_the_field()
        max_D = self.max_D()

        R = PowerSeriesRing(QQ, names=['x', 'y'], default_prec=max_D + 1)
        x, y = R.gens()

        Hk = HilbertSeries(n=k, degrees=[2]*m, q=q)
        k_y, k_xy = Hk.series(y), Hk.series(x * y)

        Hn = HilbertSeries(n=n, degrees=[2]*m, q=q)
        n_x = Hn.series(x)

        N = NMonomialSeries(n=n - k, q=q, max_prec=max_D + 1)
        nk_x = N.series_monomials_of_degree()(x)

        return (k_xy * nk_x - n_x - k_y) / ((1 - x) * (1 - y))

    def admissible_parameters(self):
        """
        Return a list of admissible parameters `(k, D, d)`

        EXAMPLES::

            sage: from mpkc.algorithms import Crossbred
            sage: E = Crossbred(n=10, m=12, q=5)
            sage: E.admissible_parameters()[:5]
            [(1, 2, 1), (1, 3, 1), (1, 4, 1), (1, 3, 2), (1, 5, 1)]
        """
        n = self.nvariables()
        max_D = self.max_D()

        admissible_parameters = []
        for k in range(1, n):
            Sk = self.admissible_parameter_series(k)
            possibles_D_d = [monomial.exponents()[0]
                             for (monomial, coefficient) in Sk.coefficients().items()
                             if (0 <= coefficient)
                             and (monomial.exponents()[0][0] > monomial.exponents()[0][1])
                             and (monomial.exponents()[0][0] <= max_D)
                             and (1 <= monomial.exponents()[0][1])]
            admissible_parameters.extend([(k, D, d) for D, d in possibles_D_d])

        return admissible_parameters

    def time_complexity(self, **kwargs):
        """
        Return the time complexity

        INPUT:

        - ``k`` -- no. of variables in the resulting system (default: None)
        - ``D`` -- degree of the initial Macaulay matrix (default: None)
        - ``d`` -- degree resulting Macaulay matrix (default: None)

        If `k`, `D`, and `d` are specified, the function returns the time complexity w.r.t to the given parameters

        EXAMPLES::

            sage: from mpkc.algorithms import Crossbred
            sage: E = Crossbred(n=10, m=12, q=5)
            sage: float(log(E.time_complexity(), 2))
            20.3663736923649
            sage: float(log(E.time_complexity(k=4, D=6, d=4), 2))
            29.775157881382952
        """
        k = kwargs.get('k', None)
        D = kwargs.get('D', None)
        d = kwargs.get('d', None)

        if all(var is not None for var in (k, D, d)):
            return self._time_complexity_(k, D, d)

        min_time_complexity = Infinity
        if self._time_complexity is None:
            for (k, D, d) in self.admissible_parameters():
                time_complexity = self._time_complexity_(k, D, d)
                if time_complexity < min_time_complexity:
                    min_time_complexity = time_complexity
                    self._k = k
                    self._D = D
                    self._d = d
            self._time_complexity = min_time_complexity

        return self._time_complexity

    def memory_complexity(self):
        """
        Return the memory complexity

        EXAMPLES::

            sage: from mpkc.algorithms import Crossbred
            sage: E = Crossbred(n=10, m=12, q=5)
            sage: float(log(E.memory_complexity(), 2))
            8.027905996569885
        """
        if self._memory_complexity is None:
            D = self.D()
            k = self.k()
            d = self.d()
            ncols_pre_step = self.ncols_in_preprocessing_step(k, D, d)
            ncols_lin_step = self.ncols_in_linearization_step(k, d)
            self._memory_complexity = ncols_pre_step ** 2 + ncols_lin_step ** 2

        return self._memory_complexity

    def tilde_o_time(self):
        """
        Return the Ō time complexity of crossbred algorithm

        EXAMPLES::

            sage: from mpkc.algorithms import Crossbred
            sage: E = Crossbred(n=10, m=12, q=5)
            sage: float(log(E.tilde_o_time(), 2))
            16.782447984412244
        """
        k, D, d = self.k(), self.D(), self.d()
        np = self.ncols_in_preprocessing_step(k=k, D=D, d=d)
        nl = self.ncols_in_linearization_step(k=k, d=d)
        q = self.order_of_the_field()
        n = self.nvariables()
        w = self.linear_algebra_constant()

        return np ** 2 + q ** (n - k) * nl ** w

    def _time_complexity_(self, k, D, d):
        n, m = self.nvariables(), self.npolynomials()
        w = self.linear_algebra_constant()
        q = self.order_of_the_field()
        np = self.ncols_in_preprocessing_step(k=k, D=D, d=d)
        nl = self.ncols_in_linearization_step(k=k, d=d)
        complexity = Infinity
        if np > 1 and log(np, 2) > 1:
            complexity = (np ** 2 * log(np, 2) * log(log(np, 2), 2)) + (m * q ** (n - k) * nl ** w)
        return complexity

    def __repr__(self):
        return "Crossbred estimator for the MQ problem"
