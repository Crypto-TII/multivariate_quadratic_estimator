from sage.functions.log import log
from sage.functions.other import floor
from sage.rings.infinity import Infinity
from .base import BaseAlgorithm, optimal_parameter
from ..utils import sum_of_binomial_coefficients


class DinurFirst(BaseAlgorithm):
    r"""
    Construct an instance of Dinur's first estimator

    The Dinur's first is a probabilistic algorithm to solve the MQ problem over GF(2) [Din21a]_. It computes the parity
    of the number of solutions of many quadratic polynomial systems. These systems come from the specialization, in the
    original system, of the values in a fixed set of variables.


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

    # all methods below are implemented to overwrite the parent's docstring while keeping the implementation

    def has_optimal_parameter(self):
        """
        Return `True` if the algorithm has optimal parameter

        EXAMPLES::

            sage: from mpkc.algorithms import DinurFirst
            sage: H = DinurFirst(n=5, m=10)
            sage: H.has_optimal_parameter()
            True
        """
        return super().has_optimal_parameter()

    def is_defined_over_finite_field(self):
        """
        Return `True` if the algorithm is defined over a finite field

        EXAMPLES::

            sage: from mpkc.algorithms import DinurFirst
            sage: H = DinurFirst(n=5, m=10)
            sage: H.is_defined_over_finite_field()
            True
        """
        return super().is_defined_over_finite_field()

    def is_overdefined_system(self):
        """
        Return `True` if the system is overdefined

        EXAMPLES::

            sage: from mpkc.algorithms import DinurFirst
            sage: H = DinurFirst(n=5, m=10)
            sage: H.is_overdefined_system()
            True
            sage: E = DinurFirst(n=10, m=10)
            sage: E.is_overdefined_system()
            False
        """
        return super().is_overdefined_system()

    def is_square_system(self):
        """
        Return `True` if the system is square, there are equal no. of variables and polynomials

        EXAMPLES::

            sage: from mpkc.algorithms import DinurFirst
            sage: H = DinurFirst(n=5, m=10)
            sage: H.is_square_system()
            False
            sage: E = DinurFirst(n=10, m=10)
            sage: E.is_square_system()
            True
        """
        return super().is_square_system()

    def is_underdefined_system(self):
        """
        Return `True` if the system is underdefined

        EXAMPLES::

            sage: from mpkc.algorithms import DinurFirst
            sage: H = DinurFirst(n=5, m=10)
            sage: H.is_underdefined_system()
            False
            sage: E = DinurFirst(n=10, m=5)
            sage: E.is_underdefined_system()
            True
        """
        return super().is_underdefined_system()

    def linear_algebra_constant(self):
        """
        Return the linear algebra constant

        EXAMPLES::

            sage: from mpkc.algorithms import DinurFirst
            sage: H = DinurFirst(n=5, m=10)
            sage: H.linear_algebra_constant()
            <BLANKLINE>
        """
        return super().linear_algebra_constant()

    def npolynomials(self):
        """
        Return the number of polynomials

        EXAMPLES::

            sage: from mpkc.algorithms import DinurFirst
            sage: H = DinurFirst(n=5, m=10)
            sage: H.npolynomials()
            10
        """
        return super().npolynomials()

    def nvariables(self):
        """
        Return the number of variables

        EXAMPLES::

            sage: from mpkc.algorithms import DinurFirst
            sage: H = DinurFirst(n=5, m=10)
            sage: H.nvariables()
            5
        """
        return super().nvariables()

    def nvariables_reduced(self):
        """
        Return the no. of variables after fixing some values

        EXAMPLES::

            sage: from mpkc.algorithms import DinurFirst
            sage: H = DinurFirst(n=5, m=10)
            sage: H.nvariables_reduced()
            5
            sage: E = DinurFirst(n=12, m=10)
            sage: E.nvariables_reduced()
            10
        """
        return super().nvariables_reduced()

    def optimal_parameters(self):
        """
        Return a dictionary of optimal parameters

        EXAMPLES::

            sage: from mpkc.algorithms import DinurFirst
            sage: H = DinurFirst(n=15, m=10)
            sage: H.optimal_parameters()
            {'κ': 2/9, 'λ': 1/9}
        """
        return super().optimal_parameters()

    def order_of_the_field(self):
        """
        Return the order of the field

        EXAMPLES::

            sage: from mpkc.algorithms import DinurFirst
            sage: H = DinurFirst(n=15, m=10)
            sage: H.order_of_the_field()
            2
        """
        return super().order_of_the_field()