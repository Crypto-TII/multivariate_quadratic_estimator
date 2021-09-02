from sage.functions.other import binomial
from mpkc import degree_of_regularity
from .base import BaseAlgorithm


class F5(BaseAlgorithm):
    def __init__(self, n, m, q=None, w=2):
        """
        Return an instance of F5 complexity estimator

        INPUT:

        - ``n`` -- no. of variables
        - ``m`` -- no. of polynomials
        - ``q`` -- order of the base field (default: Infinity)
        - ``w`` -- linear algebra constant (default: 2)

        EXAMPLES::

            sage: from mpkc.algorithms import F5
            sage: F5_ = F5(n=10, m=5)
            sage: F5_
            Complexity estimator for F5 with 10 variables and 5 polynomials
        """
        super().__init__(n, m, q=q, w=w)

    def time_complexity(self, degrees):
        """
        Return the time complexity of the F5 algorithm

        INPUT:

        - ``degrees`` -- a list/tuple of the degree of the polynomials

        EXAMPLES::

            sage: from mpkc.algorithms import F5
            sage: F5_ = F5(n=10, m=15)
            sage: F5_.time_complexity(degrees=[3]*15)
            8533694884
        """
        if self.is_overdefined_system():
            complexity = self.time_complexity_semi_regular_system(degrees)
        else:
            complexity = self.time_complexity_regular_system(degrees)

        return complexity

    def time_complexity_regular_system(self, degrees):
        """
        Return the time complexity for regular system

        INPUT:

        - ``degrees`` -- a list of integers representing the degree of the polynomials

        EXAMPLES::

            sage: from mpkc.algorithms import F5
            sage: F5_ = F5(n=10, m=5)
            sage: F5_.time_complexity_regular_system(degrees=[2]*5)
            626250625

        TESTS::

            sage: F5(n=15, m=5).time_complexity_regular_system(degrees=[2]*5)
            37558440000
        """
        if not (self.is_square_system() or self.is_underdefined_system()):
            raise ValueError("regularity assumption is valid only on square or underdefined system")

        if len(degrees) != self.npolynomials():
            raise ValueError(f"len(degrees) must be equal to {self.npolynomials()}")

        n, m = self.nvariables(), self.npolynomials()
        w = self.linear_algebra_constant()
        dreg = degree_of_regularity.regular_system(n, degrees)
        return (m * binomial(n + dreg - 1, dreg)) ** w

    def time_complexity_semi_regular_system(self, degrees):
        """
        Return the time complexity for semi-regular system

        INPUT:

        - ``degrees`` -- a list of integers representing the degree of the polynomials

        EXAMPLES::

            sage: from mpkc.algorithms import F5
            sage: F5_ = F5(n=5, m=10)
            sage: F5_.time_complexity_semi_regular_system([2]*10)
            3136

        TESTS::

            sage: F5(n=5, m=15).time_complexity_semi_regular_system([2]*15)
            441
        """
        if not self.is_overdefined_system():
            raise ValueError("semi regularity assumption is valid only on overdefined system")

        if len(degrees) != self.npolynomials():
            raise ValueError(f"len(degrees) must be equal to {self.npolynomials()}")

        n = self.nvariables()
        w = self.linear_algebra_constant()
        q = self.order_of_the_field()
        dreg = degree_of_regularity.semi_regular_system(n, degrees, q)
        return binomial(n + dreg, dreg) ** w

    def time_complexity_quadratic_system(self):
        """
        Return the time complexity for quadratic system

        EXAMPLES::

            sage: from mpkc.algorithms import F5
            sage: F5_ = F5(n=5, m=10)
            sage: F5_.time_complexity_quadratic_system()
            3136

        TESTS::

            sage: F5(n=10, m=5).time_complexity_quadratic_system()  # underdefined system
            626250625
            sage: F5(n=10, m=10).time_complexity_quadratic_system()  # square system
            2821056160000
        """
        return self.time_complexity([2]*self.npolynomials())

    def __repr__(self):
        n, m = self.nvariables(), self.npolynomials()
        return f"Complexity estimator for F5 with {n} variables and {m} polynomials"
