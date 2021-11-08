
from sage.all import Integer
from sage.arith.misc import is_power_of_two
from .base import BaseAlgorithm


class KPG(BaseAlgorithm):
    r"""
    Construct an instance of KPG estimator

    The KPG is an algorithm to solve a quadratic systems of equations over fields of even characteristic [KPG99]_.

    INPUT:

    - ``n`` -- no. of variables
    - ``m`` -- no. of polynomials
    - ``q`` -- order of the finite field
    - ``w`` -- linear algebra constant (2 <= w <= 3) (default: 2)

    EXAMPLES::

        sage: from mpkc.algorithms import KPG
        sage: E = KPG(n=183, m=12, q=4, w=2.8)
        sage: E
        KPG estimator for the MQ problem

    TESTS::

        sage: E.nvariables() == E.nvariables_reduced()
        True
    """
    def __init__(self, n, m, q, w=2):
        if not isinstance(q, (int, Integer)):
            raise TypeError("q must be an integer")

        if not is_power_of_two(q):
            raise ValueError("the order of finite field q must be a power of 2")

        if not m * (m + 1) < n:
            raise ValueError(f'The condition m(m + 1) < n must be satisfied')

        super().__init__(n=n, m=m, q=q, w=w)
        self._n_reduced = n

    def time_complexity(self):
        """
        Return the time complexity of kpg algorithm

        EXAMPLES::

            sage: from mpkc.algorithms import KPG
            sage: E = KPG(n=183, m=12, q=4, w=2.8)
            sage: float(log(E.time_complexity(),2))
            24.628922047916475
        """
        n, m = self.nvariables(), self.npolynomials()
        w = self.linear_algebra_constant()
        return m * n ** w

    def memory_complexity(self):
        """
        Return the memory complexity of kpg algorithm

        EXAMPLES::

            sage: from mpkc.algorithms import KPG
            sage: E = KPG(n=183, m=12, q=4, w=2.8)
            sage: E.memory_complexity()
            401868
        """
        n, m = self.nvariables(), self.npolynomials()
        return m * n ** 2

    def tilde_o_time(self):
        """
        Return the Ō time complexity

        EXAMPLES::

            sage: from mpkc.algorithms import KPG
            sage: E = KPG(n=183, m=12, q=4, w=2.8)
            sage: E.tilde_o_time()
            1
        """
        return 1

    def __repr__(self):
        return f"KPG estimator for the MQ problem"

    # all methods below are implemented to overwrite the parent's docstring while keeping the implementation

    def has_optimal_parameter(self):
        """
        Return `True` if the algorithm has optimal parameter

        EXAMPLES::

            sage: from mpkc.algorithms import KPG
            sage: H = KPG(n=183, m=12, q=4, w=2.8)
            sage: H.has_optimal_parameter()
            False
        """
        return super().has_optimal_parameter()

    def is_defined_over_finite_field(self):
        """
        Return `True` if the algorithm is defined over a finite field

        EXAMPLES::

            sage: from mpkc.algorithms import KPG
            sage: H = KPG(n=183, m=12, q=4, w=2.8)
            sage: H.is_defined_over_finite_field()
            True
        """
        return super().is_defined_over_finite_field()

    def is_overdefined_system(self):
        """
        Return `True` if the system is overdefined

        EXAMPLES::

            sage: from mpkc.algorithms import KPG
            sage: H = KPG(n=183, m=12, q=4, w=2.8)
            sage: H.is_overdefined_system()
            False
        """
        return super().is_overdefined_system()

    def is_square_system(self):
        """
        Return `True` if the system is square, there are equal no. of variables and polynomials

        EXAMPLES::

            sage: from mpkc.algorithms import KPG
            sage: H = KPG(n=183, m=12, q=4, w=2.8)
            sage: H.is_square_system()
            False
        """
        return super().is_square_system()

    def is_underdefined_system(self):
        """
        Return `True` if the system is underdefined

        EXAMPLES::

            sage: from mpkc.algorithms import KPG
            sage: E = KPG(n=183, m=12, q=4, w=2.8)
            sage: E.is_underdefined_system()
            True
        """
        return super().is_underdefined_system()

    def linear_algebra_constant(self):
        """
        Return the linear algebra constant

        EXAMPLES::

            sage: from mpkc.algorithms import KPG
            sage: H = KPG(n=183, m=12, q=4, w=3)
            sage: H.linear_algebra_constant()
            3
        """
        return super().linear_algebra_constant()

    def npolynomials(self):
        """
        Return the number of polynomials

        EXAMPLES::

            sage: from mpkc.algorithms import KPG
            sage: H = KPG(n=183, m=12, q=4, w=2.8)
            sage: H.npolynomials()
            12
        """
        return super().npolynomials()

    def nvariables(self):
        """
        Return the number of variables

        EXAMPLES::

            sage: from mpkc.algorithms import KPG
            sage: H = KPG(n=183, m=12, q=4, w=2.8)
            sage: H.nvariables()
            183
        """
        return super().nvariables()

    def nvariables_reduced(self):
        """
        Return the no. of variables after fixing some values

        EXAMPLES::

            sage: from mpkc.algorithms import KPG
            sage: H = KPG(n=183, m=12, q=4, w=2.8)
            sage: H.nvariables_reduced()
            183
        """
        return super().nvariables_reduced()

    def optimal_parameters(self):
        """
        Return a dictionary of optimal parameters

        EXAMPLES::

            sage: from mpkc.algorithms import KPG
            sage: H = KPG(n=183, m=12, q=4, w=2.8)
            sage: H.optimal_parameters()
            {}
        """
        return super().optimal_parameters()

    def order_of_the_field(self):
        """
        Return the order of the field

        EXAMPLES::

            sage: from mpkc.algorithms import KPG
            sage: H = KPG(n=183, m=12, q=4, w=2.8)
            sage: H.order_of_the_field()
            4
        """
        return super().order_of_the_field()