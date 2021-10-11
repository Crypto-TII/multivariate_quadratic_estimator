from sage.all import Integer
from sage.functions.log import log
from sage.misc.functional import numerical_approx
from .base import BaseAlgorithm


class ExhaustiveSearch(BaseAlgorithm):
    def __init__(self, n, m, q, nsolutions=1):
        """
        Construct an instance of Exhaustive Search estimator

        INPUT:

        - ``n`` -- no. of variables
        - ``m`` -- no. of polynomials
        - ``q`` -- order of the finite field
        - ``nsolutions`` -- number of solutions (default: 1)

        EXAMPLES::

            sage: from mpkc.algorithms import ExhaustiveSearch
            sage: E = ExhaustiveSearch(q=3, n=10, m=12)
            sage: E
            Exhaustive search estimator for the MQ problem
        """
        if not isinstance(q, (int, Integer)):
            raise TypeError("q must be an integer")

        super().__init__(n=n, m=m, q=q)
        self._nsolutions = nsolutions

    def nsolutions(self):
        """
        Return the number of solutions

        EXAMPLES::

            sage: from mpkc.algorithms import ExhaustiveSearch
            sage: E = ExhaustiveSearch(q=3, n=10, m=12, nsolutions=3)
            sage: E.nsolutions()
            3
        """
        return self._nsolutions

    def time_complexity(self):
        """
        Return the time complexity of the exhaustive search algorithm

        EXAMPLES::

            sage: from mpkc.algorithms import ExhaustiveSearch
            sage: E = ExhaustiveSearch(q=3, n=10, m=12)
            sage: E.time_complexity()
            61880.4962217569

        TESTS::

            sage: E0 = ExhaustiveSearch(n=15, m=12, q=3)
            sage: E1 = ExhaustiveSearch(n=17, m=12, q=3)
            sage: E0.time_complexity() == E1.time_complexity()
            True
        """
        n = self.nvariables_reduced()
        nsolutions = self.nsolutions()
        q = self.order_of_the_field()
        if q == 2:
            complexity = 4 * log(n, 2) * (2 ** n / (nsolutions + 1))
        else:
            complexity = log(n, q) * (q ** n / (nsolutions + 1))

        return numerical_approx(complexity)

    def memory_complexity(self):
        """
        Return the memory complexity of the exhaustive search algorithm

        EXAMPLES::

            sage: from mpkc.algorithms import ExhaustiveSearch
            sage: E = ExhaustiveSearch(q=3, n=10, m=12)
            sage: E.memory_complexity()
            1200

        TESTS::

            sage: E0 = ExhaustiveSearch(n=15, m=12, q=3)
            sage: E1 = ExhaustiveSearch(n=17, m=12, q=3)
            sage: E0.memory_complexity() == E1.memory_complexity()
            True
        """
        n, m = self.nvariables_reduced(), self.npolynomials()
        return m * n ** 2

    def tilde_o_time(self):
        """
        Return the Ō time complexity of the exhaustive search algorithm

        EXAMPLES::

            sage: from mpkc.algorithms import ExhaustiveSearch
            sage: E = ExhaustiveSearch(q=3, n=10, m=12)
            sage: E.tilde_o_time()
            59049
        """
        q = self.order_of_the_field()
        n = self.nvariables_reduced()
        return q ** n

    def __repr__(self):
        return f"Exhaustive search estimator for the MQ problem"
