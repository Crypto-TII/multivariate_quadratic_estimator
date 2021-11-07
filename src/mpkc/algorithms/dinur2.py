
from sage.functions.log import log
from sage.rings.infinity import Infinity
from ..utils import sum_of_binomial_coefficients
from .base import BaseAlgorithm, optimal_parameter


class DinurSecond(BaseAlgorithm):
    """
    Construct an instance of Dinur's second estimator

    Dinur's second is a probabilistic algorithm to solve the MQ problem over GF(2) [Din21b]_. It is based on ideas from
    [Din21a]_.

    INPUT:

    - ``n`` -- no. of variables
    - ``m`` -- no. of polynomials

    EXAMPLES::

        sage: from mpkc.algorithms import DinurSecond
        sage: E = DinurSecond(n=10, m=12)
        sage: E
        Dinur's second estimator for the MQ problem
    """
    def __init__(self, n, m):
        super().__init__(n=n, m=m, q=2)
        self._n1 = None
        self._time_complexity = None
        self._memory_complexity = None

    @optimal_parameter
    def n1(self):
        """
        Return the optimal parameter $n_1$

        EXAMPLES::

            sage: from mpkc.algorithms import DinurSecond
            sage: E = DinurSecond(n=10, m=12)
            sage: E.n1()
            4
        """
        if self._n1 is None:
            self._compute_time_complexity_()
        return self._n1

    def time_complexity(self, **kwargs):
        """
        Return the time complexity of the Dinur's second algorithm

        INPUT:

        - ``n1`` -- the parameter $n_1$ (default: None)

        If $n_1$ is provided, the function returns the time complexity w.r.t. the given parameter

        EXAMPLES::

            sage: from mpkc.algorithms import DinurSecond
            sage: E = DinurSecond(n=10, m=12)
            sage: E.time_complexity().numerical_approx()
            57434.4699066345
            sage: E.time_complexity(n1=2).numerical_approx()
            58848.1441779413

        TESTS::

            sage: E0 = DinurSecond(n=15, m=12)
            sage: E1 = DinurSecond(n=17, m=12)
            sage: E0.time_complexity().numerical_approx() == E1.time_complexity().numerical_approx()
            True
        """

        n1 = kwargs.get("n1", None)

        if n1 is not None:
            time_complexity = self._time_complexity_(n1)
        else:
            if self._time_complexity is None:
                self._compute_time_complexity_()
            time_complexity = self._time_complexity

        return time_complexity

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
        n = self.nvariables_reduced()
        self._memory_complexity = 8 * (n1 + 1) * sum_of_binomial_coefficients(n - n1, n1 + 3)
        return self._memory_complexity

    def _compute_time_complexity_(self):
        n, m = self.nvariables_reduced(), self.npolynomials()
        max_n1 = ((m - 2) // 2) - 1
        min_time_complexity = Infinity
        optimal_n1 = None

        for n1 in range(1, max_n1 + 1):
            time_complexity = self._time_complexity_(n1)
            if time_complexity < min_time_complexity:
                optimal_n1 = n1
                min_time_complexity = time_complexity

        self._n1 = optimal_n1
        self._time_complexity = min_time_complexity

    def _time_complexity_(self, n1):
        """
        Return the time complexity for the given parameter

        INPUT:

        - ``n1`` -- the parameter $n_1$
        """
        n = self.nvariables_reduced()

        return 16 * log(n, 2) * 2 ** n1 * sum_of_binomial_coefficients(n - n1, n1 + 3) + \
               n1 * n * 2 ** (n - n1) + \
               2 ** (n - 2 * n1 + 1) * sum_of_binomial_coefficients(n, 2)

    def tilde_o_time(self):
        """
        Return the Ō time complexity of dinur's second algorithm

        EXAMPLES::

            sage: from mpkc.algorithms.dinur2 import DinurSecond
            sage: E = DinurSecond(n=10, m=12)
            sage: E.tilde_o_time()
            283.68541077888506
        """
        return 2 ** ((1 - 1./(2.7*2)) * self._n)

    def __repr__(self):
        return f"Dinur's second estimator for the MQ problem"

    # all methods below are implemented to overwrite the parent's docstring while keeping the implementation

    def has_optimal_parameter(self):
        """
        Return `True` if the algorithm has optimal parameter

        EXAMPLES::

            sage: from mpkc.algorithms import DinurSecond
            sage: H = DinurSecond(n=5, m=10)
            sage: H.has_optimal_parameter()
            True
        """
        return super().has_optimal_parameter()

    def is_defined_over_finite_field(self):
        """
        Return `True` if the algorithm is defined over a finite field

        EXAMPLES::

            sage: from mpkc.algorithms import DinurSecond
            sage: H = DinurSecond(n=5, m=10)
            sage: H.is_defined_over_finite_field()
            True
        """
        return super().is_defined_over_finite_field()

    def is_overdefined_system(self):
        """
        Return `True` if the system is overdefined

        EXAMPLES::

            sage: from mpkc.algorithms import DinurSecond
            sage: H = DinurSecond(n=5, m=10)
            sage: H.is_overdefined_system()
            True
            sage: E = DinurSecond(n=10, m=10)
            sage: E.is_overdefined_system()
            False
        """
        return super().is_overdefined_system()

    def is_square_system(self):
        """
        Return `True` if the system is square, there are equal no. of variables and polynomials

        EXAMPLES::

            sage: from mpkc.algorithms import DinurSecond
            sage: H = DinurSecond(n=5, m=10)
            sage: H.is_square_system()
            False
            sage: E = DinurSecond(n=10, m=10)
            sage: E.is_square_system()
            True
        """
        return super().is_square_system()

    def is_underdefined_system(self):
        """
        Return `True` if the system is underdefined

        EXAMPLES::

            sage: from mpkc.algorithms import DinurSecond
            sage: H = DinurSecond(n=5, m=10)
            sage: H.is_underdefined_system()
            False
            sage: E = DinurSecond(n=10, m=5)
            sage: E.is_underdefined_system()
            True
        """
        return super().is_underdefined_system()

    def linear_algebra_constant(self):
        """
        Return the linear algebra constant

        EXAMPLES::

            sage: from mpkc.algorithms import DinurSecond
            sage: H = DinurSecond(n=5, m=10)
            sage: H.linear_algebra_constant()
            <BLANKLINE>
        """

    def npolynomials(self):
        """
        Return the number of polynomials

        EXAMPLES::

            sage: from mpkc.algorithms import DinurSecond
            sage: H = DinurSecond(n=5, m=10)
            sage: H.npolynomials()
            10
        """
        return super().npolynomials()

    def nvariables(self):
        """
        Return the number of variables

        EXAMPLES::

            sage: from mpkc.algorithms import DinurSecond
            sage: H = DinurSecond(n=5, m=10)
            sage: H.nvariables()
            5
        """
        return super().nvariables()

    def nvariables_reduced(self):
        """
        Return the no. of variables after fixing some values

        EXAMPLES::

            sage: from mpkc.algorithms import DinurSecond
            sage: H = DinurSecond(n=5, m=10)
            sage: H.nvariables_reduced()
            5
            sage: E = DinurSecond(n=12, m=10)
            sage: E.nvariables_reduced()
            10
        """
        return super().nvariables_reduced()

    def optimal_parameters(self):
        """
        Return a dictionary of optimal parameters

        EXAMPLES::

            sage: from mpkc.algorithms import DinurSecond
            sage: H = DinurSecond(n=15, m=10)
            sage: H.optimal_parameters()
            {'n1': 2}
        """
        return super().optimal_parameters()

    def order_of_the_field(self):
        """
        Return the order of the field

        EXAMPLES::

            sage: from mpkc.algorithms import DinurSecond
            sage: H = DinurSecond(n=15, m=10)
            sage: H.order_of_the_field()
            2
        """
        return super().order_of_the_field()
