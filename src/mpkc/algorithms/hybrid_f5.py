from .base import BaseAlgorithm, optimal_parameter
from .f5 import F5


class HybridF5(BaseAlgorithm):
    def __init__(self, q, n, m, w=2, use_quantum=False, **kwargs):
        """
        Return an instance of hybrid approach complexity estimator

        INPUT:

        - ``q`` -- order of the finite field
        - ``n`` -- no. of variables
        - ``m`` -- no. of polynomials
        - ``w`` -- linear algebra constant (2 <= w <= 3) (default: 2)
        - ``use_quantum`` -- return the complexity using quantum computer (default: False)
        - ``degrees`` -- a list/tuple of degree of the polynomials (default: [2]*m)

        EXAMPLES::

            sage: from mpkc.algorithms import HybridF5
            sage: H = HybridF5(q=256, n=5, m=10)
            sage: H
            Complexity estimator for hybrid approach with 5 variables and 10 polynomials
        """
        degrees = kwargs.get('degrees', [2] * m)
        if len(degrees) != m:
            raise ValueError(f"len(degrees) must be equal to {m}")

        super().__init__(n, m, q=q, w=w)
        self._use_quantum = use_quantum
        self._degrees = degrees

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
        """
        n = self.nvariables()
        complexities = [self._time_complexity_(k) for k in range(n)]
        return min(range(len(complexities)), key=complexities.__getitem__)

    def time_complexity(self, **kwargs):
        """
        Return the complexity of hybrid approach

        INPUT:

        - ``k`` -- no. of fixed variables

        If `k` is specified, the function returns the time complexity w.r.t the given parameter

        EXAMPLES::

            sage: from mpkc.algorithms import HybridF5
            sage: H = HybridF5(q=256, n=10, m=10)
            sage: H.time_complexity()
            6412806400

        TESTS::

            sage: H = HybridF5(q=256, n=10, m=15)
            sage: H.time_complexity()
            1002001
            sage: H.time_complexity(k=2)
            1784217600
        """
        n = self.nvariables()
        k = kwargs.get('k', self.k())

        if not 0 <= k <= n:
            raise ValueError(f'k must be in the range 0 <= k <= {n}')

        return self._time_complexity_(k)

    def memory_complexity(self):
        """
        Return the memory complexity

        EXAMPLES::

            sage: from mpkc.algorithms import HybridF5
            sage: E = HybridF5(n=10, m=12, q=7)
            sage: E.memory_complexity()
            14400
        """
        n, m = self.nvariables(), self.npolynomials()
        q = self.order_of_the_field()
        w = self.linear_algebra_constant()
        degrees = self.degree_of_polynomials()
        k = self.k()

        return F5(n=n-k, m=m, q=q, w=w, degrees=degrees).memory_complexity()

    def tilde_o_time(self):
        """
        Return the Ō time complexity of hybrid-F5 algorithm for quadratic system

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
        n, m = self.nvariables(), self.npolynomials()
        if self.is_underdefined_system():
            n -= (n - m)

        q = self.order_of_the_field()
        w = self.linear_algebra_constant()
        degrees = self.degree_of_polynomials()

        return q ** (k / 2 if self.use_quantum() else k) * F5(n=n-k, m=m, q=q, w=w, degrees=degrees).time_complexity()

    def __repr__(self):
        n, m = self.nvariables(), self.npolynomials()
        return f"Complexity estimator for hybrid approach with {n} variables and {m} polynomials"
