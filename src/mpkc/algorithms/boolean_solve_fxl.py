"""
Module to compute the time and memory complexity of the algorithms BooleanSolve and FXL

The BooleanSolve and the FXL are algorithms to solve the MQ problem

[BFS+11] Bardet, M., Faugère, J.-C., Salvy, B., and Spaenlehauer, P.-J. On the complexity of solving quadratic
boolean systems. CoRR,abs/1112.6263, 2011.

[YC04]  Courtois, N., and Klimov, A., and Patarin, J., and Shamir, A. Efficient  algorithms  for  solving overdefined systems of multivariate polynomial
equations, In B. Preneel, editor,Advancesin Cryptology — EUROCRYPT 2000, pages 392–407, Berlin, Heidelberg, 2000.
SpringerBerlin Heidelberg.
"""
from sage.arith.misc import binomial
from sage.rings.infinity import Infinity
from .base import BaseAlgorithm, optimal_parameter
from .. import witness_degree


class BooleanSolveFXL(BaseAlgorithm):
    variants = ("las_vegas", "deterministic")

    def __init__(self, n, m, q=None, w=2):
        """
        Construct an instance of BooleanSolve and FXL estimator

        INPUT:

        - ``n`` -- no. of variables
        - ``m`` -- no. of polynomials
        - ``q`` -- order of the finite field
        - ``w`` -- linear algebra constant (2 <= w <= 3) (default: 2)

        EXAMPLES::

            sage: from mpkc.algorithms import BooleanSolveFXL
            sage: E = BooleanSolveFXL(n=10, m=12, q=7)
            sage: E
            BooleanSolve and FXL estimators for the MQ problem
        """
        super().__init__(n=n, m=m, q=q, w=w)

        if self.is_defined_over_finite_field():
            if not self.is_overdefined_system() and not self.is_square_system():
                raise ValueError("the no. of polynomials must be >= than the no. of variables")
        else:
            if not self.is_overdefined_system():
                raise ValueError("the no. of polynomials must be > than the no. of variables")

        self._k = None
        self._variant = None
        self._time_complexity = None
        self._memory_complexity = None
        self._compute_optimal_k_ = self._compute_time_complexity_
        self._compute_optimal_variant_ = self._compute_time_complexity_

    @optimal_parameter
    def k(self):
        """
        Return the optimal `k`

        EXAMPLES::

            sage: from mpkc.algorithms import BooleanSolveFXL
            sage: E = BooleanSolveFXL(n=10, m=12, q=7)
            sage: E.k()
            4
        """
        if self._k is None:
            self._compute_optimal_k_()
        return self._k

    @optimal_parameter
    def variant(self):
        """
        Return the optimal variant

        EXAMPLES::

            sage: from mpkc.algorithms import BooleanSolveFXL
            sage: E = BooleanSolveFXL(n=10, m=12, q=7)
            sage: E.variant()
            'deterministic'
        """
        if self._variant is None:
            self._compute_optimal_variant_()
        return self._variant

    def time_complexity(self, **kwargs):
        """
        Return the time complexity of BooleanSolve and FXL algorithms

        INPUT:

        - ``k`` -- the optimal `k` (default: None)
        - ``variant`` -- the selected variant (default: None)

        If `k` and `variant` are specified, the function returns the time complexity w.r.t the given parameters

        EXAMPLES::

            sage: from mpkc.algorithms import BooleanSolveFXL
            sage: E = BooleanSolveFXL(n=10, m=12, q=7)
            sage: float(log(E.time_complexity(), 2))
            27.599017034509096
            sage: float(log(E.time_complexity(k=2, variant="las_vegas"), 2))
            33.35111811760744
        """

        k = kwargs.get('k', None)
        variant = kwargs.get('variant', None)

        if k is not None and variant is not None:
            time_complexity = self._time_complexity_(k, variant)
        else:
            self._compute_time_complexity_()
            time_complexity = self._time_complexity

        return time_complexity

    def memory_complexity(self):
        """
        Return the memory complexity of BooleanSolve and FXL algorithms

        EXAMPLES::

            sage: from mpkc.algorithms import BooleanSolveFXL
            sage: E = BooleanSolveFXL(n=10, m=12, q=7)
            sage: E.memory_complexity()
            7056
        """
        if self._memory_complexity is None:
            n, m = self.nvariables(), self.npolynomials()
            q = self.order_of_the_field()
            k = self.k()
            wit_deg = witness_degree.quadratic_system(n=n - k, m=m, q=q)
            self._memory_complexity = max(binomial(n - k + wit_deg, wit_deg) ** 2, m * n ** 2)

        return self._memory_complexity

    def tilde_o_time(self):
        """
        Return the Ō time complexity of BooleanSolve and FXL algorithms

        EXAMPLES::

            sage: from mpkc.algorithms import BooleanSolveFXL
            sage: E = BooleanSolveFXL(n=10, m=12, q=7)
            sage: float(log(E.tilde_o_time(), 2))
            24.014054533787938
        """
        n, m = self.nvariables(), self.npolynomials()
        q = self.order_of_the_field()
        w = self.linear_algebra_constant()
        k = self.k()
        variant = self.variant()
        wit_deg = witness_degree.quadratic_system(n=n - k, m=m, q=q)

        if variant == 'las_vegas':
            complexity = q ** k * binomial(n - k + wit_deg, wit_deg) ** 2
        else:
            complexity = q ** k * binomial(n - k + wit_deg, wit_deg) ** w

        return complexity

    def _compute_time_complexity_(self):
        min_time_complexity = Infinity

        n, m = self.nvariables(), self.npolynomials()

        optimal_k = optimal_variant = None

        for variant in BooleanSolveFXL.variants:
            a = 0 if self.is_overdefined_system() else 1
            for k in range(a, n):

                time_complexity = self._time_complexity_(k, variant)

                if time_complexity < min_time_complexity:
                    min_time_complexity = time_complexity
                    optimal_k = k
                    optimal_variant = variant

        self._time_complexity = min_time_complexity
        self._k = optimal_k
        self._variant = optimal_variant

    def _time_complexity_(self, k, variant):
        """
        Return the time complexity for the given parameter

        INPUT:

        - ``k`` -- the value `k`
        - ``variant`` -- the variant of the algorithm
        """
        n, m = self.nvariables(), self.npolynomials()
        q = self.order_of_the_field()
        w = self.linear_algebra_constant()

        wit_deg = witness_degree.quadratic_system(n=n - k, m=m, q=q)

        if variant == "las_vegas":
            time_complexity = 3 * binomial(n - k + 2, 2) * q ** k * binomial(n - k + wit_deg, wit_deg) ** 2
        elif variant == "deterministic":
            time_complexity = q ** k * m * binomial(n - k + wit_deg, wit_deg) ** w
        else:
            raise ValueError("variant must either be las_vegas or deterministic")

        return time_complexity

    def __repr__(self):
        return f"BooleanSolve and FXL estimators for the MQ problem"
