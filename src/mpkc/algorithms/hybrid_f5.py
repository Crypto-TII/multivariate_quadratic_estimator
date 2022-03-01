from sage.all import Integer
from .base import BaseAlgorithm, optimal_parameter
from .f5 import F5


class HybridF5(BaseAlgorithm):
    r"""
    Construct an instance of HybridF5

    HybridF5 is an algorithm to solve systems of polynomials over a finite field proposed in [BFP09]_, [BFP12]_. The
    algorithm is a tradeoff between exhaustive search and Groebner bases computation. The idea is to fix the value of,
    say, $k$ variables and compute the Groebner bases of $q^{k}$ subsystems, where $q$ is the order of the finite
    field. The Grobner bases computation is done using F5 algorithm.

    .. SEEALSO::

        :class:`mpkc.algorithms.f5.F5` -- class to compute the complexity of F5 algorithm.

    INPUT:

    - ``n`` -- no. of variables
    - ``m`` -- no. of polynomials
    - ``q`` -- order of the finite field
    - ``w`` -- linear algebra constant (2 <= w <= 3) (default: 2)
    - ``use_quantum`` -- return the complexity using quantum computer (default: False)
    - ``degrees`` -- a list/tuple of degree of the polynomials (default: [2]*m, i.e. quadratic system)

    EXAMPLES::

        sage: from mpkc.algorithms import HybridF5
        sage: H = HybridF5(q=256, n=5, m=10)
        sage: H
        Complexity estimator for hybrid approach with 5 variables and 10 polynomials
    """
    def __init__(self, n, m, q, w=2, use_quantum=False, **kwargs):
        if not isinstance(q, (int, Integer)):
            raise TypeError("q must be an integer")

        degrees = kwargs.get('degrees', [2] * m)
        if len(degrees) != m:
            raise ValueError(f"len(degrees) must be equal to {m}")

        super().__init__(n, m, q=q, w=w)
        self._use_quantum = use_quantum
        if degrees == [2] * m:
            self._degrees = [2] * self.npolynomials_reduced()
        else:
            self._degrees = degrees
        self._time_complexity = None
        self._memory_complexity = None

    def degree_of_polynomials(self):
        """
        Return a list of degree of the polynomials

        EXAMPLES::

            sage: from mpkc.algorithms import HybridF5
            sage: H = HybridF5(q=31, n=5, m=5, degrees=[3]*5)
            sage: H.degree_of_polynomials()
            [3, 3, 3, 3, 3]
        """
        return self._degrees

    def use_quantum(self):
        """
        Return `True` if the complexity computation is in quantum model

        EXAMPLES::

            sage: from mpkc.algorithms import HybridF5
            sage: H = HybridF5(q=31, n=5, m=10, use_quantum=False)
            sage: H.use_quantum()
            False
            sage: H = HybridF5(q=31, n=5, m=10, use_quantum=True)
            sage: H.use_quantum()
            True
        """
        return self._use_quantum

    @optimal_parameter
    def k(self):
        """
        Return `k`, i.e. the optimal no. of fixed variables

        EXAMPLES::

            sage: from mpkc.algorithms import HybridF5
            sage: H = HybridF5(q=31, n=23, m=23)
            sage: H.k()
            2

        TESTS::

            sage: H = HybridF5(q=256, n=10, m=10)
            sage: H.k()
            1
            sage: H = HybridF5(q=256, n=20, m=10)
            sage: H.k()
            1
        """
        n = self.nvariables_reduced()
        complexities = [self._time_complexity_(k) for k in range(n)]
        return min(range(len(complexities)), key=complexities.__getitem__)

    def time_complexity(self, **kwargs):
        """
        Return the complexity of HybridF5

        INPUT:

        - ``k`` -- no. of fixed variables (default: None)

        .. NOTE::

            If ``k`` is specified, the function returns the time complexity w.r.t the given parameter

        .. SEEALSO::

            :meth:`mpkc.algorithms.f5.F5.time_complexity`

        EXAMPLES::

            sage: from mpkc.algorithms import HybridF5
            sage: H = HybridF5(q=256, n=10, m=10)
            sage: H.time_complexity()
            6412806400
            sage: H.time_complexity(k=2)
            108551798784

        TESTS::

            sage: H = HybridF5(q=256, n=10, m=15)
            sage: H.time_complexity()
            1002001
            sage: H.time_complexity(k=2)
            1784217600
        """
        k = kwargs.get('k', self.k())

        if k == self.k():
            if self._time_complexity is None:
                self._time_complexity = self._time_complexity_(k)
                time_complexity = self._time_complexity
            else:
                time_complexity = self._time_complexity
        else:
            n = self.nvariables_reduced()
            if not 0 <= k <= n:
                raise ValueError(f'k must be in the range 0 <= k <= {n}')
            else:
                time_complexity = self._time_complexity_(k)

        return time_complexity

    def memory_complexity(self):
        """
        Return the memory complexity of HybridF5

        .. SEEALSO::

            :meth:`mpkc.algorithms.f5.F5.memory_complexity`

        EXAMPLES::

            sage: from mpkc.algorithms import HybridF5
            sage: E = HybridF5(n=10, m=12, q=7)
            sage: E.memory_complexity()
            7056
        """
        if self._memory_complexity is None:
            n, m = self.nvariables_reduced(), self.npolynomials_reduced()
            q = self.order_of_the_field()
            w = self.linear_algebra_constant()
            degrees = self.degree_of_polynomials()
            k = self.k()
            self._memory_complexity = F5(n=n-k, m=m, q=q, w=w, degrees=degrees).memory_complexity()

        return self._memory_complexity

    def tilde_o_time(self):
        """
        Return the Ō time complexity of HybridF5

        .. SEEALSO::

            :meth:`mpkc.algorithms.f5.F5.tilde_o_time`

        EXAMPLES::

            sage: from mpkc.algorithms import HybridF5
            sage: E = HybridF5(n=10, m=12, q=7)
            sage: E.tilde_o_time()
            4939200
        """
        return self.time_complexity()

    def _time_complexity_(self, k):
        """
        Return the time complexity w.r.t. `k`.

        INPUT:

        - ``k`` -- no. of fixed variables
        """
        n, m = self.nvariables_reduced(), self.npolynomials_reduced()
        q = self.order_of_the_field()
        w = self.linear_algebra_constant()
        degrees = self.degree_of_polynomials()

        return q ** (k / 2 if self.use_quantum() else k) * F5(n=n-k, m=m, q=q, w=w, degrees=degrees).time_complexity()

    def __repr__(self):
        n, m = self.nvariables(), self.npolynomials_reduced()
        return f"Complexity estimator for hybrid approach with {n} variables and {m} polynomials"

    # all methods below are implemented to overwrite the parent's docstring while keeping the implementation

    def has_optimal_parameter(self):
        """
        Return `True` if the algorithm has optimal parameter

        EXAMPLES::

            sage: from mpkc.algorithms import HybridF5
            sage: H = HybridF5(q=256, n=5, m=10)
            sage: H.has_optimal_parameter()
            True
        """
        return super().has_optimal_parameter()

    def is_defined_over_finite_field(self):
        """
        Return `True` if the algorithm is defined over a finite field

        EXAMPLES::

            sage: from mpkc.algorithms import HybridF5
            sage: H = HybridF5(q=256, n=5, m=10)
            sage: H.is_defined_over_finite_field()
            True
        """
        return super().is_defined_over_finite_field()

    def is_overdefined_system(self):
        """
        Return `True` if the system is overdefined

        EXAMPLES::

            sage: from mpkc.algorithms import HybridF5
            sage: H = HybridF5(q=256, n=5, m=10)
            sage: H.is_overdefined_system()
            True
            sage: E = HybridF5(q=256, n=10, m=10)
            sage: E.is_overdefined_system()
            False
        """
        return super().is_overdefined_system()

    def is_square_system(self):
        """
        Return `True` if the system is square, there are equal no. of variables and polynomials

        EXAMPLES::

            sage: from mpkc.algorithms import HybridF5
            sage: H = HybridF5(q=256, n=5, m=10)
            sage: H.is_square_system()
            False
            sage: E = HybridF5(q=256, n=10, m=10)
            sage: E.is_square_system()
            True
        """
        return super().is_square_system()

    def is_underdefined_system(self):
        """
        Return `True` if the system is underdefined

        EXAMPLES::

            sage: from mpkc.algorithms import HybridF5
            sage: H = HybridF5(q=256, n=5, m=10)
            sage: H.is_underdefined_system()
            False
            sage: E = HybridF5(q=256, n=10, m=5)
            sage: E.is_underdefined_system()
            True
        """
        return super().is_underdefined_system()

    def linear_algebra_constant(self):
        """
        Return the linear algebra constant

        EXAMPLES::

            sage: from mpkc.algorithms import HybridF5
            sage: H = HybridF5(q=256, n=5, m=10, w=2)
            sage: H.linear_algebra_constant()
            2
        """
        return super().linear_algebra_constant()

    def npolynomials(self):
        """
        Return the number of polynomials

        EXAMPLES::

            sage: from mpkc.algorithms import HybridF5
            sage: H = HybridF5(q=256, n=5, m=10)
            sage: H.npolynomials()
            10
        """
        return super().npolynomials()

    def nvariables(self):
        """
        Return the number of variables

        EXAMPLES::

            sage: from mpkc.algorithms import HybridF5
            sage: H = HybridF5(q=256, n=5, m=10)
            sage: H.nvariables()
            5
        """
        return super().nvariables()

    def nvariables_reduced(self):
        """
        Return the no. of variables after fixing some values

        EXAMPLES::

            sage: from mpkc.algorithms import HybridF5
            sage: H = HybridF5(q=256, n=5, m=10)
            sage: H.nvariables_reduced()
            5
            sage: E = HybridF5(q=256, n=12, m=10)
            sage: E.nvariables_reduced()
            9
        """
        return super().nvariables_reduced()

    def npolynomials_reduced(self):
        """
        Return the no. of polynomials after applying the Thomae and Wolf strategy

        EXAMPLES::

            sage: from mpkc.algorithms import HybridF5
            sage: H = HybridF5(q=256, n=5, m=10)
            sage: H.npolynomials_reduced()
            10
        """
        return super().npolynomials_reduced()

    def optimal_parameters(self):
        """
        Return a dictionary of optimal parameters

        EXAMPLES::

            sage: from mpkc.algorithms import HybridF5
            sage: H = HybridF5(q=256, n=15, m=10)
            sage: H.optimal_parameters()
            {'k': 1}
        """
        return super().optimal_parameters()

    def order_of_the_field(self):
        """
        Return the order of the field

        EXAMPLES::

            sage: from mpkc.algorithms import HybridF5
            sage: H = HybridF5(q=256, n=15, m=10)
            sage: H.order_of_the_field()
            256
        """
        return super().order_of_the_field()
