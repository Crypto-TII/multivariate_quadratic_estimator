from .base import BaseAlgorithm


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

    def optimal_nfixed_vars(self):
        """
        Return the optimal no. of fixed variables

        EXAMPLES::

            sage: from mpkc.algorithms import HybridF5
            sage: H = HybridF5(q=31, n=23, m=23)
            sage: H.optimal_nfixed_vars()
            2

        TESTS::

            sage: H = HybridF5(q=256, n=10, m=10)
            sage: H.optimal_nfixed_vars()
            1
        """
        min_finder = lambda iterable: min(range(len(iterable)), key=iterable.__getitem__)
        return self._hybrid_approach_(min_finder)

    def time_complexity(self, **kwargs):
        """
        Return the complexity of hybrid approach

        EXAMPLES::

            sage: from mpkc.algorithms import HybridF5
            sage: H = HybridF5(q=256, n=10, m=10)
            sage: H.time_complexity()
            6412806400

        TESTS::

            sage: H = HybridF5(q=256, n=10, m=15)
            sage: H.time_complexity()
            1002001
        """
        min_finder = lambda iterable: min(iterable)
        return self._hybrid_approach_(min_finder)

    def _hybrid_approach_(self, min_finder):
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
        degrees = self.degree_of_polynomials()

        complexities = [
            multiplier(q, k) * F5(n=n-k, m=m, q=q, w=w, degrees=degrees).time_complexity() for k in range(n)
        ]
        return min_finder(complexities)

    def __repr__(self):
        n, m = self.nvariables(), self.npolynomials()
        return f"Complexity estimator for hybrid approach with {n} variables and {m} polynomials"
