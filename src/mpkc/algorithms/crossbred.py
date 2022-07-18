# *****************************************************************************
# Multivariate Quadratic (MQ) Estimator
# Copyright (C) 2021-2022 Emanuele Bellini, Rusydi H. Makarim, Javier Verbel
# Cryptography Research Centre, Technology Innovation Institute LLC
#
# This file is part of MQ Estimator
#
# MQ Estimator is free software: you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.
#
# MQ Estimator is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with
# MQ Estimator. If not, see <https://www.gnu.org/licenses/>.
# *****************************************************************************


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
    r"""
    Construct an instance of crossbred estimator

    The Crossbred is an algorithm to solve the MQ problem [JV18]_. This algorithm consists of two steps, named the
    preprocessing step and the linearization step. In the preprocessing step, we find a set $S$ of degree-$D$
    polynomials in the ideal generated by the initial set of polynomials. Every specialization  of the first $n-k$
    variables of the polynomials in $S$ results in a set $S'$ of degree-$d$ polynomials in $k$ variables. Finally, in
    the linearization step, a solution to $S'$ is found by direct linearization.

    .. NOTE::

        Our complexity estimates are a generalization over any field of size `q` of the complexity formulas given in
        [Dua20]_, which are given either for `q=2` or generic fields.

    INPUT:

    - ``n`` -- no. of variables
    - ``m`` -- no. of polynomials
    - ``q`` -- order of the finite field
    - ``w`` -- linear algebra constant (2 <= w <= 3) (default: 2)
    - ``max_D`` -- upper bound to the parameter D (default: 10)
    - ``h`` -- external hybridization parameter (default: 0)

    EXAMPLES::

        sage: from mpkc.algorithms import Crossbred
        sage: E = Crossbred(n=10, m=12, q=5)
        sage: E
        Crossbred estimator for the MQ problem
    """
    def __init__(self, n, m, q, w=2, max_D=10, h=0):
        if not isinstance(q, (int, Integer)):
            raise TypeError("q must be an integer")

        super().__init__(n=n, m=m, q=q, w=w, h=h)
        self._max_D = max_D
        self._k = None
        self._D = None
        self._d = None
        self._time_complexity = None
        self._memory_complexity = None

    @property
    def max_D(self):
        """
        Return the upper bound of the degree of the initial Macaulay matrix

        EXAMPLES::

            sage: from mpkc.algorithms import Crossbred
            sage: E = Crossbred(n=10, m=12, q=5, max_D=9)
            sage: E.max_D
            9
        """
        return self._max_D

    @max_D.setter
    def max_D(self, value):
        """
        Set new upper bound of the degree of the initial Macaulay matrix

        TESTS::

            sage: from mpkc.algorithms import Crossbred
            sage: E = Crossbred(n=10, m=12, q=5, max_D=6)
            sage: E.admissible_parameter_series(1)
            -1 - 2*x - 2*y - 2*x*y + 9*y^2 + 65*x^3 + 9*x^2*y + 9*x*y^2 + 20*y^3 + 439*x^4 + 119*x^3*y + 9*x^2*y^2 +
            20*x*y^3 - 35*y^4 + 1705*x^5 + 658*x^4*y + 20*x^3*y^2 + 20*x^2*y^3 - 35*x*y^4 - 89*y^5 + 4892*x^6 +
            2419*x^5*y + 64*x^4*y^2 + 20*x^3*y^3 - 35*x^2*y^4 - 89*x*y^5 + 77*y^6 + O(x, y)^7
            sage: E.max_D = 5
            sage: E.admissible_parameter_series(1)
            -1 - 2*x - 2*y - 2*x*y + 9*y^2 + 65*x^3 + 9*x^2*y + 9*x*y^2 + 20*y^3 + 439*x^4 + 119*x^3*y + 9*x^2*y^2 +
            20*x*y^3 - 35*y^4 + 1705*x^5 + 658*x^4*y + 20*x^3*y^2 + 20*x^2*y^3 - 35*x*y^4 - 89*y^5 + O(x, y)^6
        """
        self._max_D = value
        self._k = None
        self._D = None
        self._d = None
        self._time_complexity = None
        self._memory_complexity = None

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

        n = self.nvariables_reduced()
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
        Return a the series $S_k$ of admissible parameters

        INPUT:

        - ``k`` -- no. variables in the resulting system

        EXAMPLES::

            sage: from mpkc.algorithms import Crossbred
            sage: E = Crossbred(n=10, m=12, q=5, max_D=2)
            sage: E.admissible_parameter_series(2)
            -1 - 3*x - 3*y - 10*x^2 - 3*x*y + 6*y^2 + O(x, y)^3
        """
        n, m = self.nvariables_reduced(), self.npolynomials_reduced()
        q = self.order_of_the_field()
        max_D = self.max_D

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
        n = self.nvariables_reduced()
        max_D = self.max_D

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

        TESTS::

            sage: E0 = Crossbred(n=15, m=12, q=5)
            sage: E1 = Crossbred(n=16, m=12, q=5)
            sage: E0.time_complexity().numerical_approx() == E1.time_complexity().numerical_approx()
            True
        """
        k = kwargs.get('k', None)
        D = kwargs.get('D', None)
        d = kwargs.get('d', None)

        h = self._h
        if all(var is not None for var in (k, D, d)):
            return 2 ** h * self._time_complexity_(k, D, d)

        min_time_complexity = Infinity
        if self._time_complexity is None:
            for (k, D, d) in self.admissible_parameters():
                time_complexity = self._time_complexity_(k, D, d)
                if time_complexity < min_time_complexity:
                    min_time_complexity = time_complexity
                    self._k = k
                    self._D = D
                    self._d = d
            self._time_complexity = 2 ** h * min_time_complexity

        return self._time_complexity

    def memory_complexity(self, **kwargs):
        """
        Return the memory complexity

        INPUT:

        - ``k`` -- no. of variables in the resulting system (default: None)
        - ``D`` -- degree of the initial Macaulay matrix (default: None)
        - ``d`` -- degree resulting Macaulay matrix (default: None)

        If `k`, `D`, and `d` are specified, the function returns the memory complexity w.r.t to the given parameters

        EXAMPLES::

            sage: from mpkc.algorithms import Crossbred
            sage: E = Crossbred(n=10, m=12, q=5)
            sage: float(log(E.memory_complexity(), 2))
            8.027905996569885
            sage: float(log(E.memory_complexity(k=4, D=6, d=4), 2))
            12.892542816648552

        TESTS::

            sage: E0 = Crossbred(n=15, m=12, q=5)
            sage: E1 = Crossbred(n=16, m=12, q=5)
            sage: E0.memory_complexity().numerical_approx() == E1.memory_complexity().numerical_approx()
            True
        """

        k = kwargs.get('k', None)
        D = kwargs.get('D', None)
        d = kwargs.get('d', None)

        if all(var is not None for var in (k, D, d)):
            ncols_pre_step = self.ncols_in_preprocessing_step(k, D, d)
            ncols_lin_step = self.ncols_in_linearization_step(k, d)
            return ncols_pre_step ** 2 + ncols_lin_step ** 2

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
        n = self.nvariables_reduced()
        w = self.linear_algebra_constant()
        h = self._h
        return 2 ** h * np ** 2 + q ** (n - k) * nl ** w

    def _time_complexity_(self, k, D, d):
        n, m = self.nvariables_reduced(), self.npolynomials_reduced()
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

    # all methods below are implemented to overwrite the parent's docstring while keeping the implementation

    def has_optimal_parameter(self):
        """
        Return `True` if the algorithm has optimal parameter

        EXAMPLES::

            sage: from mpkc.algorithms import Crossbred
            sage: H = Crossbred(q=256, n=5, m=10)
            sage: H.has_optimal_parameter()
            True
        """
        return super().has_optimal_parameter()

    def is_defined_over_finite_field(self):
        """
        Return `True` if the algorithm is defined over a finite field

        EXAMPLES::

            sage: from mpkc.algorithms import Crossbred
            sage: H = Crossbred(q=256, n=5, m=10)
            sage: H.is_defined_over_finite_field()
            True
        """
        return super().is_defined_over_finite_field()

    def is_overdefined_system(self):
        """
        Return `True` if the system is overdefined

        EXAMPLES::

            sage: from mpkc.algorithms import Crossbred
            sage: H = Crossbred(q=256, n=5, m=10)
            sage: H.is_overdefined_system()
            True
            sage: E = Crossbred(q=256, n=10, m=10)
            sage: E.is_overdefined_system()
            False
        """
        return super().is_overdefined_system()

    def is_square_system(self):
        """
        Return `True` if the system is square, there are equal no. of variables and polynomials

        EXAMPLES::

            sage: from mpkc.algorithms import Crossbred
            sage: H = Crossbred(q=256, n=5, m=10)
            sage: H.is_square_system()
            False
            sage: E = Crossbred(q=256, n=10, m=10)
            sage: E.is_square_system()
            True
        """
        return super().is_square_system()

    def is_underdefined_system(self):
        """
        Return `True` if the system is underdefined

        EXAMPLES::

            sage: from mpkc.algorithms import Crossbred
            sage: H = Crossbred(q=256, n=5, m=10)
            sage: H.is_underdefined_system()
            False
            sage: E = Crossbred(q=256, n=10, m=5)
            sage: E.is_underdefined_system()
            True
        """
        return super().is_underdefined_system()

    def linear_algebra_constant(self):
        """
        Return the linear algebra constant

        EXAMPLES::

            sage: from mpkc.algorithms import Crossbred
            sage: H = Crossbred(q=256, n=5, m=10, w=2)
            sage: H.linear_algebra_constant()
            2
        """
        return super().linear_algebra_constant()

    def npolynomials(self):
        """
        Return the number of polynomials

        EXAMPLES::

            sage: from mpkc.algorithms import Crossbred
            sage: H = Crossbred(q=256, n=5, m=10)
            sage: H.npolynomials()
            10
        """
        return super().npolynomials()

    def nvariables(self):
        """
        Return the number of variables

        EXAMPLES::

            sage: from mpkc.algorithms import Crossbred
            sage: H = Crossbred(q=256, n=5, m=10)
            sage: H.nvariables()
            5
        """
        return super().nvariables()

    def nvariables_reduced(self):
        """
        Return the no. of variables after fixing some values

        EXAMPLES::

            sage: from mpkc.algorithms import Crossbred
            sage: H = Crossbred(q=256, n=5, m=10)
            sage: H.nvariables_reduced()
            5
            sage: E = Crossbred(q=256, n=12, m=10)
            sage: E.nvariables_reduced()
            9
        """
        return super().nvariables_reduced()

    def npolynomials_reduced(self):
        """
        Return the no. of polynomials after applying the Thomae and Wolf strategy

        EXAMPLES::

            sage: from mpkc.algorithms import Crossbred
            sage: H = Crossbred(q=256, n=5, m=10)
            sage: H.npolynomials_reduced()
            10
            sage: E = Crossbred(q=256, n=12, m=10)
            sage: E.npolynomials_reduced()
            9
        """
        return super().npolynomials_reduced()

    def optimal_parameters(self):
        """
        Return a dictionary of optimal parameters

        EXAMPLES::

            sage: from mpkc.algorithms import Crossbred
            sage: H = Crossbred(q=256, n=15, m=10)
            sage: H.optimal_parameters()
            {'D': 7, 'd': 2, 'k': 7}
        """
        return super().optimal_parameters()

    def order_of_the_field(self):
        """
        Return the order of the field

        EXAMPLES::

            sage: from mpkc.algorithms import Crossbred
            sage: H = Crossbred(q=256, n=15, m=10)
            sage: H.order_of_the_field()
            256
        """
        return super().order_of_the_field()