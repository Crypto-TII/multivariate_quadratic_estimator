from .base import BaseAlgorithm


class HybridApproach(BaseAlgorithm):
    def __init__(self, q, n, m, w=2, use_quantum=False):
        """
        Return an instance of hybrid approach complexity estimator

        INPUT:

        - ``q`` -- order of the finite field
        - ``n`` -- no. of variables
        - ``m`` -- no. of polynomials
        - ``w`` -- linear algebra constant (2 <= w <= 3) (default: 2)
        - ``use_quantum`` -- return the complexity using quantum computer (default: False)

        EXAMPLES::

            sage: from mpkc.algorithms import HybridApproach
            sage: H = HybridApproach(q=256, n=5, m=10)
            sage: H
            Complexity estimator for hybrid approach with 5 variables and 10 polynomials
        """
        super().__init__(n, m, q=q, w=w)
        self._use_quantum = use_quantum

    def use_quantum(self):
        """
        Return `True` if the complexity computation is in quantum model

        EXAMPLES::

            sage: from mpkc.algorithms import HybridApproach
            sage: H = HybridApproach(q=31, n=5, m=10, use_quantum=False)
            sage: H.use_quantum()
            False
            sage: H = HybridApproach(q=31, n=5, m=10, use_quantum=True)
            sage: H.use_quantum()
            True
        """
        return self._use_quantum

    def optimal_nfixed_vars(self, degrees):
        """
        Return the optimal no. of fixed variables

        INPUT:

        - ``degrees`` -- a list of integers representing the degree of the polynomials

        EXAMPLES::

            sage: from mpkc.algorithms import HybridApproach
            sage: H = HybridApproach(q=31, n=23, m=23)
            sage: H.optimal_nfixed_vars([2]*23)
            2

        TESTS::

            sage: H = HybridApproach(q=256, n=10, m=10)
            sage: H.optimal_nfixed_vars([2]*10)
            1
        """
        min_finder = lambda iterable: min(range(len(iterable)), key=iterable.__getitem__)
        return self._hybrid_approach_(degrees, min_finder)

    def optimal_nfixed_vars_quadratic_system(self):
        """
        Return the optimal no. of fixed variables for quadratic system

        EXAMPLES::

            sage: from mpkc.algorithms import HybridApproach
            sage: H = HybridApproach(q=31, n=23, m=23)
            sage: H.optimal_nfixed_vars_quadratic_system()
            2
        """
        return self.optimal_nfixed_vars([2]*self.npolynomials())

    def time_complexity(self, degrees):
        """
        Return the complexity of hybrid approach

        INPUT:

        - ``degrees`` -- a list of integers representing the degree of the polynomials

        EXAMPLES::

            sage: from mpkc.algorithms import HybridApproach
            sage: H = HybridApproach(q=256, n=10, m=10)
            sage: H.time_complexity([2]*10)
            6412806400

        TESTS::

            sage: H = HybridApproach(q=256, n=10, m=15)
            sage: H.time_complexity([2]*15)
            1002001
        """
        min_finder = lambda iterable: min(iterable)
        return self._hybrid_approach_(degrees, min_finder)

    def time_complexity_quadratic_system(self):
        """
        Return the time complexity for quadratic system

        EXAMPLES::

            sage: from mpkc.algorithms import HybridApproach
            sage: H = HybridApproach(q=256, n=10, m=15)
            sage: H.time_complexity_quadratic_system()
            1002001
        """
        return self.time_complexity([2]*self.npolynomials())

    def _hybrid_approach_(self, degrees, min_finder):
        if len(degrees) != self.npolynomials():
            raise ValueError(f"len(degrees) must be equal to {self.npolynomials()}")

        n, m = self.nvariables(), self.npolynomials()
        if self.is_underdefined_system():
            n -= (n - m)

        if self.use_quantum():
            def multiplier(x, e):
                return x ** (e / 2)
        else:
            def multiplier(x, e):
                return x ** e

        from .f5 import F5

        q = self.order_of_the_field()
        w = self.linear_algebra_constant()

        complexities = [
            multiplier(q, k) * F5(n=n-k, m=m, q=q, w=w).time_complexity(degrees) for k in range(n)
        ]
        return min_finder(complexities)

    def __repr__(self):
        n, m = self.nvariables(), self.npolynomials()
        return f"Complexity estimator for hybrid approach with {n} variables and {m} polynomials"
