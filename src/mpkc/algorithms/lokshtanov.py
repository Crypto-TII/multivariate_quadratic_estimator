
from sage.all import Integer
from sage.arith.misc import is_power_of_two
from sage.functions.log import log
from sage.functions.other import floor
from sage.rings.infinity import Infinity
from sage.rings.finite_rings.finite_field_constructor import GF
from .base import BaseAlgorithm, optimal_parameter
from ..series.nmonomial import NMonomialSeries


class Lokshtanov(BaseAlgorithm):
    r"""
    Construct an instance of Lokshtanov et al.'s estimator

    Lokshtanov et al.'s is a probabilistic algorithm to solve the MQ problem over GF(q) [LPTWY17]_. It describes an
    algorithm to determine the consistency of a given system of polynomial equations.


    INPUT:

    - ``n`` -- no. of variables
    - ``m`` -- no. of polynomials
    - ``q`` -- order of the finite field
    - ``h`` -- external hybridization parameter (default: 0)

    EXAMPLES::

        sage: from mpkc.algorithms import Lokshtanov
        sage: E = Lokshtanov(n=10, m=12, q=9)
        sage: E
        Lokshtanov et al.'s estimator for the MQ problem
    """
    def __init__(self, n, m, q, h=0):
        if not isinstance(q, (int, Integer)):
            raise TypeError("q must be an integer")

        super().__init__(n=n, m=m, q=q, h=h)
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
        n, m = self.nvariables_reduced(), self.npolynomials_reduced()

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

        TESTS::

            sage: E0 = Lokshtanov(n=15, m=12, q=9)
            sage: E1 = Lokshtanov(n=17, m=12, q=9)
            sage: E0.time_complexity().numerical_approx() == E1.time_complexity().numerical_approx()
            True
        """
        q = self.order_of_the_field()
        n = self.nvariables_reduced()
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

        h = self._h
        time_complexity *= 2 ** h
        return time_complexity

    def memory_complexity(self):
        """
        Return the memory complexity of Lokshtanov et al.'s algorithm

        EXAMPLES::

            sage: from mpkc.algorithms import Lokshtanov
            sage: E = Lokshtanov(n=10, m=12, q=9)
            sage: float(log(E.memory_complexity(), 2))
            30.622995719758727

        TESTS::

            sage: E0 = Lokshtanov(n=15, m=12, q=9)
            sage: E1 = Lokshtanov(n=17, m=12, q=9)
            sage: E0.memory_complexity().numerical_approx() == E1.memory_complexity().numerical_approx()
            True
        """
        if self._memory_complexity is None:
            δ = self.δ()
            if δ is None:
                return Infinity

            n = self.nvariables_reduced()
            np = floor(n * δ)
            q = self.order_of_the_field()
            resulting_degree = 2 * (q - 1) * (np + 2)
            M = NMonomialSeries(n=n - np, q=q, max_prec=resulting_degree + 1).nmonomials_up_to_degree(resulting_degree)
            self._memory_complexity = M + log(n, 2) * q ** (n - np)
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
        n = self.nvariables_reduced()
        if q == 2:
            time = q ** (0.8765 * n)
        elif is_power_of_two(q):
            time = q ** (0.9 * n)
        elif log(GF(q).characteristic(), 2) < 8 * e:
            time = q ** (0.9975 * n)
        else:
            d = GF(q).degree()
            time = q ** n * (log(q, 2) / (2 * e * d))

        h = self._h
        return 2 ** h * time

    def _C(self, n, delta):
        q = self.order_of_the_field()
        np = floor(delta * n)
        resulting_degree = 2 * (q - 1) * (np + 2)
        M = NMonomialSeries(n=n - np, q=q, max_prec=resulting_degree + 1).nmonomials_up_to_degree(resulting_degree)
        return n * (q ** (n - np) + M * q ** np * n ** (6 * q))

    def __repr__(self):
        return f"Lokshtanov et al.'s estimator for the MQ problem"

    # all methods below are implemented to overwrite the parent's docstring while keeping the implementation

    def has_optimal_parameter(self):
        """
        Return `True` if the algorithm has optimal parameter

        EXAMPLES::

            sage: from mpkc.algorithms import Lokshtanov
            sage: H = Lokshtanov(q=256, n=5, m=10)
            sage: H.has_optimal_parameter()
            True
        """
        return super().has_optimal_parameter()

    def is_defined_over_finite_field(self):
        """
        Return `True` if the algorithm is defined over a finite field

        EXAMPLES::

            sage: from mpkc.algorithms import Lokshtanov
            sage: H = Lokshtanov(q=256, n=5, m=10)
            sage: H.is_defined_over_finite_field()
            True
        """
        return super().is_defined_over_finite_field()

    def is_overdefined_system(self):
        """
        Return `True` if the system is overdefined

        EXAMPLES::

            sage: from mpkc.algorithms import Lokshtanov
            sage: H = Lokshtanov(q=256, n=5, m=10)
            sage: H.is_overdefined_system()
            True
            sage: E = Lokshtanov(q=256, n=10, m=10)
            sage: E.is_overdefined_system()
            False
        """
        return super().is_overdefined_system()

    def is_square_system(self):
        """
        Return `True` if the system is square, there are equal no. of variables and polynomials

        EXAMPLES::

            sage: from mpkc.algorithms import Lokshtanov
            sage: H = Lokshtanov(q=256, n=5, m=10)
            sage: H.is_square_system()
            False
            sage: E = Lokshtanov(q=256, n=10, m=10)
            sage: E.is_square_system()
            True
        """
        return super().is_square_system()

    def is_underdefined_system(self):
        """
        Return `True` if the system is underdefined

        EXAMPLES::

            sage: from mpkc.algorithms import Lokshtanov
            sage: H = Lokshtanov(q=256, n=5, m=10)
            sage: H.is_underdefined_system()
            False
            sage: E = Lokshtanov(q=256, n=10, m=5)
            sage: E.is_underdefined_system()
            True
        """
        return super().is_underdefined_system()

    def linear_algebra_constant(self):
        """
        Return the linear algebra constant

        EXAMPLES::

            sage: from mpkc.algorithms import Lokshtanov
            sage: H = Lokshtanov(q=256, n=5, m=10)
            sage: H.linear_algebra_constant()
            <BLANKLINE>
        """
        return super().linear_algebra_constant()

    def npolynomials(self):
        """
        Return the number of polynomials

        EXAMPLES::

            sage: from mpkc.algorithms import Lokshtanov
            sage: H = Lokshtanov(q=256, n=5, m=10)
            sage: H.npolynomials()
            10
        """
        return super().npolynomials()

    def nvariables(self):
        """
        Return the number of variables

        EXAMPLES::

            sage: from mpkc.algorithms import Lokshtanov
            sage: H = Lokshtanov(q=256, n=5, m=10)
            sage: H.nvariables()
            5
        """
        return super().nvariables()

    def nvariables_reduced(self):
        """
        Return the no. of variables after fixing some values

        EXAMPLES::

            sage: from mpkc.algorithms import Lokshtanov
            sage: H = Lokshtanov(q=256, n=5, m=10)
            sage: H.nvariables_reduced()
            5
            sage: E = Lokshtanov(q=256, n=12, m=10)
            sage: E.nvariables_reduced()
            9
        """
        return super().nvariables_reduced()

    def npolynomials_reduced(self):
        """
        Return the no. of polynomials after applying the Thomae and Wolf strategy

        EXAMPLES::

            sage: from mpkc.algorithms import Lokshtanov
            sage: H = Lokshtanov(q=256, n=5, m=10)
            sage: H.npolynomials_reduced()
            10
        """
        return super().npolynomials_reduced()

    def optimal_parameters(self):
        """
        Return a dictionary of optimal parameters

        EXAMPLES::

            sage: from mpkc.algorithms import Lokshtanov
            sage: H = Lokshtanov(q=256, n=15, m=10)
            sage: H.optimal_parameters()
            {'δ': 1/9}
        """
        return super().optimal_parameters()

    def order_of_the_field(self):
        """
        Return the order of the field

        EXAMPLES::

            sage: from mpkc.algorithms import Lokshtanov
            sage: H = Lokshtanov(q=256, n=15, m=10)
            sage: H.order_of_the_field()
            256
        """
        return super().order_of_the_field()