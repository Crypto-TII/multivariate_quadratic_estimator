# *****************************************************************************
# Multivariate Quadratic (MQ) Estimator
# Copyright (C) 2021-2022 Emanuele Bellini, Rusydi H. Makarim, Javier Verbel
# Cryptography Research Centre, Technology Innovation Institute LLC
#
# This file is part of MQ Estimator
#
# MQ Estimator is free software: you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.
#
# MQ Estimator is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with
# MQ Estimator. If not, see <https://www.gnu.org/licenses/>.
# *****************************************************************************


from sage.functions.log import log
from sage.functions.other import floor, ceil, binomial
from sage.misc.misc_c import prod
from sage.all import vector, ZZ
from sage.rings.finite_rings.finite_field_constructor import FiniteField
try:
    from sage.rings.polynomial.pbori.pbori import BooleanPolynomialRing
except ImportError:
    from sage.rings.polynomial.pbori import BooleanPolynomialRing
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.power_series_ring import PowerSeriesRing
from sage.rings.rational_field import QQ
from sage.structure.sequence import Sequence
from ..utils import random_affine_map


class GeMSS:
    """
    Construct an instance of GeMSS

    INPUT:

    - ``D`` -- maximum degree of the HFE polynomial
    - ``n`` -- degree of the extension field
    - ``delta`` -- no. of minus
    - ``v`` -- no. of vinegar variables in the HFEv polynomial

    EXAMPLES::

        sage: from mpkc.schemes import GeMSS
        sage: G = GeMSS(D=17, n=11, delta=4, v=3)
        sage: G
        GeMSS with D=17, n=11, Δ=4, v=3
        sage: msg = G.random_message()
        sage: signature = G.sign(msg)
        sage: G.is_valid_signature(signature, msg)
        True
    """
    def __init__(self, D, n, delta, v):
        self._base_field = FiniteField(2)

        self._max_deg_of_hfe_polynomial = D
        self._nminus = delta
        self._nvinegar_vars = v
        self._extension_field = self._base_field.extension(n)

        self._inner_affine_map = None
        self._outer_affine_map = None
        self._hfev_ring = None
        self._hfev_polynomial = None
        self._ring = None
        self._public_key = None

    @property
    def max_deg_of_hfe_polynomial(self):
        """
        Return the degree of the secret polynomial

        EXAMPLES::

            sage: from mpkc.schemes import GeMSS
            sage: G = GeMSS(D=17, n=15, delta=4, v=2)
            sage: G.max_deg_of_hfe_polynomial
            17
        """
        return self._max_deg_of_hfe_polynomial

    @property
    def nminus(self):
        """
        Return the number of minues, i.e. no. of polynomials removed from the public-key

        EXAMPLES::

            sage: from mpkc.schemes import GeMSS
            sage: G = GeMSS(D=17, n=15, delta=4, v=2)
            sage: G.nminus
            4
        """
        return self._nminus

    @property
    def nvinegar_vars(self):
        """
        Return the number of vinegar variables

        EXAMPLES::

            sage: from mpkc.schemes import GeMSS
            sage: G = GeMSS(D=17, n=15, delta=4, v=2)
            sage: G.nvinegar_vars
            2
        """
        return self._nvinegar_vars

    @property
    def base_field(self):
        """
        Return the base field

        EXAMPLES::

            sage: from mpkc.schemes import GeMSS
            sage: G = GeMSS(D=17, n=15, delta=4, v=2)
            sage: G.base_field
            Finite Field of size 2
        """
        return self._base_field

    def nvariables(self):
        """
        Return the number of variables in the system of equations, i.e. the public-key

        EXAMPLES::

            sage: from mpkc.schemes import GeMSS
            sage: G = GeMSS(D=17, n=15, delta=4, v=2)
            sage: G.nvariables()
            17
        """
        return self.extension_field.degree() + self.nvinegar_vars

    def npolynomials(self):
        """
        Return the number of polynomials in the system of equations, i.e. the public-key

        EXAMPLES::

            sage: from mpkc.schemes import GeMSS
            sage: G = GeMSS(D=17, n=15, delta=4, v=2)
            sage: G.npolynomials()
            11
        """
        return self.extension_field.degree() - self.nminus

    def inner_affine_map(self):
        """
        Return the inner affine map

        EXAMPLES::

            sage: from mpkc.schemes import GeMSS
            sage: G = GeMSS(D=17, n=3, delta=4, v=2)
            sage: M, v = G.inner_affine_map()
            sage: M  # random
            [1 1 0 0 0]
            [0 0 1 1 0]
            [1 0 1 0 1]
            [1 1 0 0 1]
            [1 0 0 0 1]
            sage: v  # random
            (0, 1, 1, 0, 1)

        TESTS::

            sage: M.is_invertible()
            True
            sage: M.dimensions() == (G.nvariables(), G.nvariables())
            True
            sage: v.length() == G.nvariables()
            True
        """
        if self._inner_affine_map is not None:
            return self._inner_affine_map

        self._inner_affine_map = random_affine_map(self.base_field, self.nvariables())
        return self._inner_affine_map

    def outer_affine_map(self):
        """
        Return the outer affine map

        EXAMPLES::

            sage: from mpkc.schemes import GeMSS
            sage: G = GeMSS(D=17, n=11, delta=4, v=2)
            sage: M, v = G.outer_affine_map()
            sage: M  # random
            [1 0 0 0 0 0 0]
            [0 1 0 1 0 0 0]
            [0 1 1 0 1 0 0]
            [1 1 0 1 0 0 1]
            [0 1 0 1 1 0 1]
            [1 1 1 1 1 1 1]
            [1 1 0 0 0 0 0]
            sage: v  # random
            (1, 1, 1, 0, 0, 0, 1)

        TESTS::

            sage: n = G.extension_field.degree()
            sage: M.is_invertible()
            True
            sage: M.dimensions() == (n, n)
            True
            sage: v.length() == n
            True
        """
        if self._outer_affine_map is not None:
            return self._outer_affine_map

        self._outer_affine_map = random_affine_map(self.base_field, self.extension_field.degree())
        return self._outer_affine_map

    @property
    def extension_field(self):
        """
        Return the extension field where the central map is defined

        TESTS::

            sage: from mpkc.schemes import GeMSS
            sage: G = GeMSS(D=17, n=11, delta=4, v=2)
            sage: G.extension_field.base_ring() == G.base_field
            True
        """
        return self._extension_field

    def hfev_ring(self):
        """
        Return the polynomial ring for the HFEv polynomial

        EXAMPLES::

            sage: from mpkc.schemes import GeMSS
            sage: G = GeMSS(D=17, n=11, delta=4, v=3)
            sage: R = G.hfev_ring()
            sage: R.gens()
            (X, v0, v1, v2)

        TESTS::

            sage: R.ngens() - 1 == G.nvinegar_vars
            True
        """
        if self._hfev_ring is not None:
            return self._hfev_ring

        var_names = ['X'] + [f'v{i}' for i in range(self.nvinegar_vars)]
        self._hfev_ring = PolynomialRing(self.extension_field, var_names)

        return self._hfev_ring

    def hfe_var(self):
        """
        Return the variable for the HFE polynomial

        EXAMPLES::

            sage: from mpkc.schemes import GeMSS
            sage: G = GeMSS(D=17, n=11, delta=4, v=3)
            sage: G.hfe_var()
            X
        """
        return self.hfev_ring().gens()[0]

    def vinegar_vars(self):
        """
        Return a tuple of vinegar variables

        EXAMPLES::

            sage: from mpkc.schemes import GeMSS
            sage: G = GeMSS(D=17, n=11, delta=4, v=3)
            sage: G.vinegar_vars()
            (v0, v1, v2)
        """
        return self.hfev_ring().gens()[1:]

    def hfev_polynomial(self):
        """
        Return the HFEv polynomial

        TESTS::

            sage: from mpkc.schemes import GeMSS
            sage: G = GeMSS(D=17, n=11, delta=4, v=3)
            sage: F = G.hfev_polynomial()
            sage: v = G.vinegar_vars()
            sage: HFE_polynomial = F.subs({v[i] : G.base_field.random_element() for i in range(G.nvinegar_vars)})
            sage: HFE_polynomial.degree() <= G.max_deg_of_hfe_polynomial
            True
        """
        if self._hfev_polynomial is not None:
            return self._hfev_polynomial

        n = self.extension_field.degree()
        D = self.max_deg_of_hfe_polynomial
        F = self.hfev_ring().zero()
        E = self.extension_field
        X = self.hfe_var()

        # "quadratic" term for the HFEv polynomial
        for i in range(n):
            for j in range(i):
                e = 2**i + 2**j
                if e > D:
                    break
                F += E.random_element() * X**e

        # "linear" term for the HFEv polynomial
        v = self.vinegar_vars()
        for i in range(n):
            e = 2**i
            if e > D:
                break

            coeffs = [E.random_element() for _ in range(self.nvinegar_vars)]
            beta = sum([c*m for (c, m) in zip(coeffs, v)])
            F += beta * X**e

        # "constant" term for the HFEv polynomial
        from itertools import combinations
        F += sum([E.random_element() * prod(vars_) for vars_ in combinations(v, 2)])

        self._hfev_polynomial = F
        return self._hfev_polynomial

    def hfev_evaluation(self, u):
        """
        Return the result of evaluation map from HFEv polynomial

        INPUT:

        - ``u`` -- a input vector over `self.base_field`
        """
        F = self.central_map()
        x = self.ring().gens()
        return vector(self.base_field, F.subs(dict(zip(x, u))))

    def ring(self):
        """
        Return the Boolean polynomial ring for the public-key

        EXAMPLES::

            sage: from mpkc.schemes import GeMSS
            sage: G = GeMSS(D=17, n=11, delta=4, v=3)
            sage: G.ring()
            Boolean PolynomialRing in x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x13
        """
        if self._ring is not None:
            return self._ring

        self._ring = BooleanPolynomialRing(self.nvariables(), [f"x{i}" for i in range(self.nvariables())], order="deglex")
        return self._ring

    def vars(self):
        """
        Return a tuple of variables for the public-key

        EXAMPLES::

            sage: from mpkc.schemes import GeMSS
            sage: G = GeMSS(D=17, n=11, delta=4, v=3)
            sage: G.vars()
            (x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x13)
        """
        return self.ring().gens()

    def central_map(self):
        """
        Return a list of multivariate polynomials representing the central map

        TESTS::

            sage: from mpkc.schemes import GeMSS
            sage: G = GeMSS(D=17, n=11, delta=4, v=3)
            sage: F = G.central_map()
            sage: x = F.ring().gens()
            sage: values = [F.ring().base_ring().random_element() for _ in range(len(x))]
            sage: y = F.subs(dict(zip(x, values)))
            sage: HFEv_F = G.hfev_polynomial()
            sage: hfev_vars = HFEv_F.variables()
            sage: E = G.extension_field
            sage: n = E.degree()
            sage: e = E.fetch_int(ZZ(values[:n], base=2))
            sage: Y = HFEv_F.subs(dict(zip(hfev_vars, [e] + values[n:]))).constant_coefficient()._vector_().list()
            sage: y == Y
            True
        """
        R = self.ring().change_ring(base_ring=self.extension_field)
        P = self.hfev_ring().change_ring(base_ring=R)
        F = self.hfev_polynomial().change_ring(P)

        x = R.gens()
        X, v = P.gens()[0], P.gens()[1:]
        g = self.extension_field.gen()

        n = self.extension_field.degree()
        values = {X: sum([g ** i * x[i] for i in range(n)])}
        values.update({v[i]: x[n + i] for i in range(self.nvinegar_vars)})

        F_subs = F.subs(values).coefficients()[0]

        polynomials = [sum([c._vector_()[i]*m for (c, m) in F_subs]) for i in range(n)]

        return Sequence(polynomials, self.ring())

    def inner_affine_polynomials(self):
        """
        Return a list of polynomials representing the inner affine map

        EXAMPLES::

            sage: from mpkc.schemes import GeMSS
            sage: G = GeMSS(D=17, n=11, delta=4, v=3)
            sage: G.inner_affine_polynomials()  # random
            [x3 + x4 + x5 + x7 + x8 + x10 + x11 + x12 + 1,
            x1 + x3 + x6 + x10 + x11 + x12,
            x0 + x1 + x2 + x3 + x4 + x5 + x6 + x8 + x10 + x12 + 1,
            x6 + x7 + x9 + x11,
            x0 + x3 + x10 + x13,
            x2 + x3 + x5 + x11 + x12 + x13 + 1,
            x0 + x1 + x3 + x9 + x11 + x12,
            x1 + x3 + x4 + x5 + x7 + x11 + x13,
            x0 + x2 + x6 + x8 + x10 + 1,
            x0 + x1 + x2 + x10 + x13 + 1,
            x1 + x3 + x10,
            x1 + x2 + x3 + x4 + x5 + x8 + x9 + x13,
            x1 + x4 + x6 + x9 + x10 + x12,
            x0 + x1 + x3 + x4 + x6 + x9 + x10 + x13]
        """
        R = self.ring()
        M, v = self.inner_affine_map()
        x = vector(R, R.gens())

        return (M*x + v).list()

    def random_signature(self):
        """
        Return a random signature

        EXAMPLES::

            sage: from mpkc.schemes import GeMSS
            sage: G = GeMSS(D=17, n=11, delta=4, v=3)
            sage: G.random_signature()  # random
            (0, 1, 0, 0, 1, 0, 1, 0, 1, 1, 1, 1, 0, 1)
        """
        return vector(self.base_field, [self.base_field.random_element() for _ in range(self.nvariables())])

    def public_key(self):
        """
        Return a list of polynomials for the public key

        TESTS::

            sage: from mpkc.schemes import GeMSS
            sage: G = GeMSS(D=17, n=11, delta=4, v=3)
            sage: T, t = G.inner_affine_map()
            sage: S, s = G.outer_affine_map()
            sage: signature = G.random_signature()
            sage: y = (S * G.hfev_evaluation(T * signature + t) + s).list()[:G.npolynomials()]
            sage: P = G.public_key()
            sage: x = G.ring().gens()
            sage: y == P.subs({x[i] : signature[i] for i in range(G.nvariables())})
            True
        """
        if self._public_key is not None:
            return self._public_key

        T = self.inner_affine_polynomials()
        F = self.central_map()

        FT = []
        R = self.ring()
        for f in F:
            p = R.zero()

            for monomial in f:
                var_indices = [int(str(v)[1:]) for v in monomial.variables()]
                p += prod(T[var_index] for var_index in var_indices)

            FT.append(p)

        S, s = self.outer_affine_map()

        self._public_key = Sequence(S * vector(FT) + s)[:self.npolynomials()]
        return self._public_key

    def is_valid_signature(self, signature, msg):
        """
        Return whether the signature is valid for the message

        INPUT:

        - ``signature`` -- a list of `self.base_field` elements
        - ``msg`` -- a list of `self.base_field` elements

        EXAMPLES::

            sage: from mpkc.schemes import GeMSS
            sage: G = GeMSS(D=17, n=11, delta=4, v=3)
            sage: msg = [G.base_field.random_element() for _ in range(G.npolynomials())]
            sage: signature = [0, 1, 1, 0, 1, 1, 0, 1, 0, 0, 1, 1, 1, 0]
            sage: G.is_valid_signature(signature, msg)  # random
            False
        """
        if len(signature) != self.nvariables():
            raise ValueError(f"signature length must be equal to {self.nvariables()}")

        if len(msg) != self.npolynomials():
            raise ValueError(f"message length must be equal to {self.npolynomials()}")

        s = list(signature)
        m = list(msg)

        P = self.public_key()
        x = self.vars()

        w = P.subs(dict(zip(x, s)))

        return m == w

    def inverse_outer_affine_map(self):
        """
        Return the inverse of the outer affine map

        EXAMPLES::

            sage: from mpkc.schemes import GeMSS
            sage: G = GeMSS(D=17, n=11, delta=4, v=3)
            sage: Si, si = G.inverse_outer_affine_map()
            sage: Si  # random
            [0 1 1 0 1 1 0 0 0 1 1]
            [1 0 1 0 1 1 1 0 1 0 0]
            [1 0 0 0 1 0 1 0 0 1 0]
            [1 0 1 1 0 1 1 0 0 1 1]
            [1 1 1 1 1 0 1 0 0 1 0]
            [0 1 1 0 0 1 0 0 0 1 1]
            [1 1 1 0 1 1 1 1 0 1 0]
            [0 1 1 0 1 1 0 1 1 1 1]
            [1 1 1 1 1 0 1 0 1 0 0]
            [0 0 1 0 1 1 1 0 0 0 1]
            [1 1 1 0 1 0 1 0 0 0 0]
            sage: si  # random
            (1, 1, 0, 1, 1, 1, 1, 1, 0, 0, 1)

        TESTS::

            sage: S, s = G.outer_affine_map()
            sage: v = VectorSpace(G.base_field, G.extension_field.degree()).random_element()
            sage: Si*(S*v + s) + si == v
            True
        """
        S, s = self.outer_affine_map()
        return S.inverse(), -S.inverse()*s

    def inverse_inner_affine_map(self):
        """
        Return the inverse of the inner affine map

        EXAMPLES::

            sage: from mpkc.schemes import GeMSS
            sage: G = GeMSS(D=17, n=11, delta=4, v=3)
            sage: Ti, ti = G.inverse_inner_affine_map()
            sage: Ti  # random
            [1 1 1 0 0 0 0 0 0 1 1 1 1 0]
            [0 1 0 0 1 0 0 0 0 1 1 1 1 1]
            [1 1 1 0 1 0 0 1 1 0 1 1 1 0]
            [0 0 0 0 1 1 0 0 1 1 0 0 0 0]
            [1 1 1 1 1 0 1 1 0 1 1 1 0 0]
            [0 0 0 1 1 0 0 1 0 0 0 0 1 1]
            [1 0 1 0 0 0 0 0 1 0 1 1 1 1]
            [0 0 0 1 1 1 0 0 0 0 0 1 1 1]
            [0 1 0 1 1 1 1 1 0 1 1 1 0 1]
            [0 0 0 0 1 0 0 1 1 0 1 1 1 1]
            [0 0 0 1 1 0 1 1 0 1 1 0 1 0]
            [1 0 1 0 1 0 1 1 1 1 1 0 1 0]
            [1 1 1 0 1 1 0 0 0 1 1 0 0 0]
            [1 1 0 1 1 0 1 0 1 1 0 1 0 0]
            sage: ti  # random
            (1, 1, 1, 1, 1, 0, 1, 1, 1, 0, 1, 0, 0, 1)

        TESTS::

            sage: T, t = G.inner_affine_map()
            sage: v = VectorSpace(G.base_field, G.nvariables()).random_element()
            sage: Ti*(T*v + t) + ti == v
            True
        """
        T, t = self.inner_affine_map()
        return T.inverse(), -T.inverse()*t

    def sign(self, msg):
        """
        Return a signature for the given message

        INPUT:

        - ``msg`` -- a list of `self.base_field` elements

        TESTS::

            sage: from mpkc.schemes import GeMSS
            sage: G = GeMSS(D=17, n=11, delta=4, v=3)
            sage: msg = G.random_message()
            sage: signature = G.sign(msg)
            sage: G.is_valid_signature(signature, msg)
            True
        """

        if len(msg) != self.npolynomials():
            raise ValueError(f"length of msg must be equal to {self.npolynomials()}")

        E = self.extension_field
        Si, si = self.inverse_outer_affine_map()
        F = self.hfev_polynomial()
        vinegar_vars = self.vinegar_vars()
        Ti, ti = self.inverse_inner_affine_map()

        roots = []
        v = []
        while len(roots) == 0:
            appended_msg = list(msg) + self.random_vector(self.nminus)
            y = vector(self.base_field, appended_msg)
            Y = E.fetch_int(ZZ((Si * y + si).list(), base=2))
            v = self.random_vector(self.nvinegar_vars)
            Fs = F.subs(dict(zip(vinegar_vars, v))).univariate_polynomial() - Y
            roots = Fs.roots(multiplicities=False)

        x = vector(self.base_field, roots[0]._vector_().list() + v)
        return Ti * x + ti

    def random_vector(self, n):
        """
        Return a random vector over `self.base_field` of length `n`

        INPUT:

        - ``n`` -- a positive integer
        """
        if n < 1:
            raise ValueError("n must be a positive integer")
        return [self.base_field.random_element() for _ in range(n)]

    def random_message(self):
        """
        Return a random message
        """
        return self.random_vector(self.npolynomials())

    def complexity_classical_exhaustive_search(self):
        """
        Return the complexity of exhaustive search in classical setting

        TESTS::

            sage: from mpkc.schemes.gemss import GeMSS128, GeMSS192, GeMSS256
            sage: G_I = GeMSS128()
            sage: G_I.complexity_classical_exhaustive_search()
            166
            sage: G_III = GeMSS192()
            sage: G_III.complexity_classical_exhaustive_search()
            247
            sage: G_V = GeMSS256()
            sage: G_V.complexity_classical_exhaustive_search()
            329
        """
        m = self.npolynomials()
        complexity = 4 * log(m, base=2) * 2**m

        return floor(log(complexity, base=2))

    def complexity_quantum_exhaustive_search_ngates(self, use_less_qubits=False):
        """
        Return the complexity of quantum exhaustive search in terms of no. of quantum gates

        INPUT:

        - ``use_less_qubits`` -- whether use the variant with less qubits but higher no. of gates (default: False)

        TESTS::

            sage: from mpkc.schemes.gemss import GeMSS128, GeMSS192, GeMSS256
            sage: G_I = GeMSS128()
            sage: G_I.complexity_quantum_exhaustive_search_ngates(use_less_qubits=False)
            104
            sage: G_I.complexity_quantum_exhaustive_search_ngates(use_less_qubits=True)
            105
            sage: G_III = GeMSS192()
            sage: G_III.complexity_quantum_exhaustive_search_ngates(use_less_qubits=False)
            146
            sage: G_III.complexity_quantum_exhaustive_search_ngates(use_less_qubits=True)
            147
            sage: G_V = GeMSS256()
            sage: G_V.complexity_quantum_exhaustive_search_ngates(use_less_qubits=False)
            188
            sage: G_V.complexity_quantum_exhaustive_search_ngates(use_less_qubits=True)
            189
        """
        n = m = self.npolynomials()

        ngates = 2**((n + 1) / 2) * (2*(m + 1) * ((n + 1)**2 + 2*(n + 1)) + 1)
        if use_less_qubits:
            ngates *= 2

        return floor(log(ngates, base=2))

    def complexity_quantum_exhaustive_search_nqubits(self, use_less_qubits=False):
        """
        Return the complexity of quantum exhaustive search in terms of no. of qubits

        INPUT:

        - ``use_less_qubits`` -- whether use the variant with less qubits but higher no. of gates (default: False)

        TESTS::

            sage: from mpkc.schemes.gemss import GeMSS128, GeMSS192, GeMSS256
            sage: G_I = GeMSS128()
            sage: G_I.complexity_quantum_exhaustive_search_nqubits(use_less_qubits=False)
            328
            sage: G_I.complexity_quantum_exhaustive_search_nqubits(use_less_qubits=True)
            174
            sage: G_III = GeMSS192()
            sage: G_III.complexity_quantum_exhaustive_search_nqubits(use_less_qubits=False)
            490
            sage: G_III.complexity_quantum_exhaustive_search_nqubits(use_less_qubits=True)
            255
            sage: G_V = GeMSS256()
            sage: G_V.complexity_quantum_exhaustive_search_nqubits(use_less_qubits=False)
            652
            sage: G_V.complexity_quantum_exhaustive_search_nqubits(use_less_qubits=True)
            337
        """
        n = m = self.npolynomials()

        if use_less_qubits:
            nqubits = 3 + (n + 1) + ceil(log(m + 1, base=2))
        else:
            nqubits = (m + 1) + (n + 1) + 2

        return nqubits

    def complexity_approximation_algorithm(self):
        """
        Return the complexity of approximation algorithm (Lokshtanov et.al.)

        TESTS::

            sage: from mpkc.schemes.gemss import GeMSS128, GeMSS192, GeMSS256
            sage: G_I = GeMSS128()
            sage: G_I.complexity_approximation_algorithm()
            141
            sage: G_III = GeMSS192()
            sage: G_III.complexity_approximation_algorithm()
            212
            sage: G_V = GeMSS256()
            sage: G_V.complexity_approximation_algorithm()
            283
        """
        m = self.npolynomials()
        return floor(0.8765 * m)

    def complexity_classical_boolean_solve(self):
        """
        Return the complexity of Boolean Solve algorithm in classical setting

        TESTS::

            sage: from mpkc.schemes.gemss import GeMSS128, GeMSS192, GeMSS256
            sage: G_I = GeMSS128()
            sage: G_I.complexity_classical_boolean_solve()
            128
            sage: G_III = GeMSS192()
            sage: G_III.complexity_classical_boolean_solve()
            192
            sage: G_V = GeMSS256()
            sage: G_V.complexity_classical_boolean_solve()
            256
        """
        from mpkc.algorithms.boolean_solve_fxl import BooleanSolveFXL
        m = self.npolynomials()
        E = BooleanSolveFXL(n=m, m=m, q=2)
        return floor(log(E.tilde_o_time(),2))

    def complexity_quantum_boolean_solve(self):
        """
        Return the complexity of Boolean Solve algorithm in quantum settings (no. of quantum gates)

        TESTS::

            sage: from mpkc.schemes.gemss import GeMSS128, GeMSS192, GeMSS256
            sage: G_I = GeMSS128()
            sage: G_I.complexity_quantum_boolean_solve()
            74
            sage: G_III = GeMSS192()
            sage: G_III.complexity_quantum_boolean_solve()
            112
            sage: G_V = GeMSS256()
            sage: G_V.complexity_quantum_boolean_solve()
            149
        """
        m = self.npolynomials()
        return floor(0.462 * m)

    def complexity_minrank_kipnis_shamir(self):
        """
        Return the complexity of minrank attack using Kipnis-Shamir approach

        TESTS::

            sage: from mpkc.schemes.gemss import GeMSS128, GeMSS192, GeMSS256
            sage: G_I = GeMSS128()
            sage: G_I.complexity_minrank_kipnis_shamir()
            521
            sage: G_III = GeMSS192()
            sage: G_III.complexity_minrank_kipnis_shamir()
            853
            sage: G_V = GeMSS256()
            sage: G_V.complexity_minrank_kipnis_shamir()
            1253
        """
        omega = 2
        v = self.nvinegar_vars
        delta = self.nminus
        D = self.max_deg_of_hfe_polynomial
        n = self.extension_field.degree()

        complexity = n**(omega * (ceil(log(D, base=2)) + v + delta + 1))

        return floor(log(complexity, base=2))

    def complexity_grobner_bases(self):
        """
        Return the complexity of direct attack using Grobner bases

        TESTS::

            sage: from mpkc.schemes.gemss import GeMSS128, GeMSS192, GeMSS256
            sage: G_I = GeMSS128()
            sage: G_I.complexity_grobner_bases()  # official result: 131
            124
            sage: G_III = GeMSS192()
            sage: G_III.complexity_grobner_bases()  # official result: 192
            185
            sage: G_V = GeMSS256()
            sage: G_V.complexity_grobner_bases()  # official result: 260
            253
        """
        m = self.npolynomials()
        dreg = self.hfev_minus_dreg()
        complexity = binomial(m, dreg)**2

        return floor(log(complexity, base=2))

    def hfev_minus_dreg(self):
        """
        Return the lower bound of the degree of regularity for the HFEv- polynomial
        """
        D = self.max_deg_of_hfe_polynomial
        delta = self.nminus
        v = self.nvinegar_vars

        R = floor(log(D - 1, base=2)) + 1
        dreg = floor((R + delta + v + 7) / 3)

        return dreg

    def complexity_minrank_with_projections(self):
        """
        Return the complexity of minrank with projection

        TESTS::

            sage: from mpkc.schemes.gemss import GeMSS128, GeMSS192, GeMSS256
            sage: G_I = GeMSS128()
            sage: G_I.complexity_minrank_with_projections()
            295
            sage: G_III = GeMSS192()
            sage: G_III.complexity_minrank_with_projections()
            444
            sage: G_V = GeMSS256()
            sage: G_V.complexity_minrank_with_projections()
            605
        """
        from sage.functions.other import sqrt

        D = self.max_deg_of_hfe_polynomial
        r = ceil(log(D, base=2))
        v = self.nvinegar_vars
        delta = self.nminus
        n = self.extension_field.degree()

        complexity = min(
            binomial(n + v + r - c, delta + v + r - c)**2 *
            binomial(n - delta, 2) *
            2**(c * (r + delta + sqrt(n - delta)) - binomial(c + 1, 2))
            for c in range(1, v + 1)
        )

        return floor(log(complexity, base=2))

    def complexity_minrank_with_support_minors(self):
        """
        Return the complexity of minrank with support minors modelling

        TESTS::

            sage: from mpkc.schemes.gemss import GeMSS128, GeMSS192, GeMSS256
            sage: G_I = GeMSS128()
            sage: G_I.complexity_minrank_with_support_minors()  # long time ; official result: 158
            153
            sage: G_III = GeMSS192()
            sage: G_III.complexity_minrank_with_support_minors()  # long time ; official result: 224
            217
            sage: G_V = GeMSS256()
            sage: G_V.complexity_minrank_with_support_minors()  # long time ; official result: 304
            290
        """
        D = self.max_deg_of_hfe_polynomial
        delta = self.nminus
        v = self.nvinegar_vars
        K = self.npolynomials()
        m = self.nvariables()
        r = ceil(log(D, base=2)) + delta + v + 1

        def is_condition_satisfied(b, n):
            return binomial(n, r) * binomial(K + b - 1, b) - 1 <= \
                   sum(
                       (-1)**(i + 1) * binomial(n, r + i) * binomial(m + i - 1, i) * binomial(K + b - i - 1, b - i)
                       for i in range(1, b + 1)
                   )

        min_complexity = 2**512
        for b in range(1, r + 2):
            for n in range(r + b, self.nvariables()):
                if not is_condition_satisfied(b, n):
                    continue
                complexity = K * (r + 1) * (binomial(n, r) * binomial(K + b - 1, b))**2
                min_complexity = min(complexity, min_complexity)

        return floor(log(min_complexity, base=2))

    def complexity_classical_distinguishing_attack(self):
        """
        Return the complexity of distinguishing attack in classical setting

        TESTS::

            sage: from mpkc.schemes.gemss import GeMSS128, GeMSS192, GeMSS256
            sage: G_I = GeMSS128()
            sage: G_I.complexity_classical_distinguishing_attack()
            226
            sage: G_III = GeMSS192()
            sage: G_III.complexity_classical_distinguishing_attack()
            346
            sage: G_V = GeMSS256()
            sage: G_V.complexity_classical_distinguishing_attack()
            476
        """
        return self.complexity_distinguishing_attack(use_quantum=False)

    def complexity_quantum_distinguishing_attack(self):
        """
        Return the complexity of distinguishing attack in quantum setting

        TESTS::

            sage: from mpkc.schemes.gemss import GeMSS128, GeMSS192, GeMSS256
            sage: G_I = GeMSS128()
            sage: G_I.complexity_quantum_distinguishing_attack()
            175
            sage: G_III = GeMSS192()
            sage: G_III.complexity_quantum_distinguishing_attack()
            265
            sage: G_V = GeMSS256()
            sage: G_V.complexity_quantum_distinguishing_attack()
            364
        """
        return self.complexity_distinguishing_attack(use_quantum=True)

    def complexity_distinguishing_attack(self, use_quantum):
        """
        Return the complexity of distinguishing attack

        INPUT:

        - `use_quantum`` -- whether to return the complexity in quantum settings (True/False)
        """
        m = self.npolynomials()
        R = PowerSeriesRing(QQ, 'z', default_prec=2*m)
        z = R.gen()

        def dreg(np):
            G = (1 + z) ** np / (1 + z ** 2) ** m
            for dreg in range(2*m):
                if G[dreg] <= 0:
                    return ZZ(dreg)

        d = self.hfev_minus_dreg()
        n = self.extension_field.degree()
        v = self.nvinegar_vars
        kbar = 0
        for kbar in range(n + v, 0, -1):
            np = n + v - kbar
            if d <= dreg(np):
                break

        min_complexity = 2**512
        for k in range(kbar):
            complexity = 3 * binomial(n + v - k, d) ** 2 * binomial(n + v - k, 2)
            if use_quantum:
                complexity *= 2**((n - k) / 2)
            else:
                complexity *= 2**(n - k)

            min_complexity = min(complexity, min_complexity)

        return floor(log(min_complexity, base=2))

    def security_level_classical(self):
        """
        Return the security level of GeMSS in classical setting

        TESTS::

            sage: from mpkc.schemes.gemss import GeMSS128, GeMSS192, GeMSS256
            sage: G_I = GeMSS128()
            sage: G_I.security_level_classical()  # long time
            124
            sage: G_III = GeMSS192()
            sage: G_III.security_level_classical()  # long time
            185
            sage: G_V = GeMSS256()
            sage: G_V.security_level_classical()  # long time
            253
        """
        sec_level = min(
            self.complexity_classical_exhaustive_search(),
            self.complexity_approximation_algorithm(),
            self.complexity_classical_boolean_solve(),
            self.complexity_minrank_kipnis_shamir(),
            self.complexity_minrank_with_projections(),
            self.complexity_minrank_with_support_minors(),
            self.complexity_grobner_bases()
        )

        return sec_level

    def security_level_quantum(self):
        """
        Return the security level of GeMSS in quantum setting

        TESTS::

            sage: from mpkc.schemes.gemss import GeMSS128, GeMSS192, GeMSS256
            sage: G_I = GeMSS128()
            sage: G_I.security_level_quantum()
            74
            sage: G_III = GeMSS192()
            sage: G_III.security_level_quantum()
            112
            sage: G_V = GeMSS256()
            sage: G_V.security_level_quantum()
            149
        """
        sec_level = min(
            self.complexity_quantum_boolean_solve(),
            self.complexity_quantum_exhaustive_search_ngates()
        )

        return sec_level

    def __repr__(self):
        D = self.max_deg_of_hfe_polynomial
        n = self.extension_field.degree()
        delta = self.nminus
        v = self.nvinegar_vars
        return f"GeMSS with D={D}, n={n}, Δ={delta}, v={v}"


def GeMSS128():
    """
    Return an instance of GeMSS with 128-bit security

    EXAMPLES::

        sage: from mpkc.schemes.gemss import GeMSS128
        sage: GeMSS128()
        GeMSS with D=513, n=174, Δ=12, v=12
    """
    return GeMSS(D=513, n=174, delta=12, v=12)


def GeMSS192():
    """
    Return an instance of GeMSS with 192-bit security

    EXAMPLES::

        sage: from mpkc.schemes.gemss import GeMSS192
        sage: GeMSS192()
        GeMSS with D=513, n=265, Δ=22, v=20
    """
    return GeMSS(D=513, n=265, delta=22, v=20)


def GeMSS256():
    """
    Return an instance of GeMSS with 256-bit security

    EXAMPLES::

        sage: from mpkc.schemes.gemss import GeMSS256
        sage: GeMSS256()
        GeMSS with D=513, n=354, Δ=30, v=33
    """
    return GeMSS(D=513, n=354, delta=30, v=33)


def BlueGeMSS128():
    """
    Return an instance of BlueGeMSS with 128-bit security

    EXAMPLES::

        sage: from mpkc.schemes.gemss import BlueGeMSS128
        sage: BlueGeMSS128()
        GeMSS with D=129, n=175, Δ=13, v=14
    """
    return GeMSS(D=129, n=175, delta=13, v=14)


def BlueGeMSS192():
    """
    Return an instance of BlueGeMSS with 192-bit security

    EXAMPLES::

        sage: from mpkc.schemes.gemss import BlueGeMSS192
        sage: BlueGeMSS192()  # long time
        GeMSS with D=129, n=265, Δ=22, v=23
    """
    return GeMSS(D=129, n=265, delta=22, v=23)


def BlueGeMSS256():
    """
    Return an instance of BlueGeMSS with 256-bit security

    EXAMPLES::

        sage: from mpkc.schemes.gemss import BlueGeMSS256
        sage: BlueGeMSS256()  # long time
        GeMSS with D=129, n=358, Δ=34, v=32
    """
    return GeMSS(D=129, n=358, delta=34, v=32)


def RedGeMSS128():
    """
    Return an instance of RedGeMSS with 128-bit security

    EXAMPLES::

        sage: from mpkc.schemes.gemss import RedGeMSS128
        sage: RedGeMSS128()
        GeMSS with D=17, n=177, Δ=15, v=15
    """
    return GeMSS(D=17, n=177, delta=15, v=15)


def RedGeMSS192():
    """
    Return an instance of RedGeMSS with 192-bit security

    EXAMPLES::

        sage: from mpkc.schemes.gemss import RedGeMSS192
        sage: RedGeMSS192()  # long time
        GeMSS with D=17, n=266, Δ=23, v=25
    """
    return GeMSS(D=17, n=266, delta=23, v=25)


def RedGeMSS256():
    """
    Return an instance of RedGeMSS with 256-bit security

    EXAMPLES::

        sage: from mpkc.schemes.gemss import RedGeMSS256
        sage: RedGeMSS256()  # long time
        GeMSS with D=17, n=358, Δ=34, v=35
    """
    return GeMSS(D=17, n=358, delta=34, v=35)


def generate_instances(sec_level_classical, sec_level_quantum,
                       min_D=1, max_D=513, min_n=1, max_n=358, min_delta=1, max_delta=34, min_v=1, max_v=35):
    """
    Return instances of GeMSS satisfying the security level

    INPUT:

    - ``sec_level_classical`` -- the classical security level
    - ``sec_level_quantum`` -- the quantum security level
    - ``min_D`` -- minimum degree of HFE polynomial
    - ``max_D`` -- maximum degree of HFE polynomial
    - ``min_n`` -- minimum degree of the extension field
    - ``max_n`` -- maximum degree of the extension field
    - ``min_delta`` -- minimum no. of minus variables
    - ``max_delta`` -- maximum no. of minus variables
    - ``min_v`` -- minimum no. of vinegar variables
    - ``max_v`` -- maximum no. of vinegar variables

    EXAMPLES::

        sage: from mpkc.schemes.gemss import generate_instances
        sage: instances = generate_instances(sec_level_classical=25, sec_level_quantum=15,
        ....:                                min_D=35, max_D=40, min_n=30, max_n=35,
        ....:                                min_delta=2, max_delta=3, min_v=1, max_v=3)  # long time
        sage: len(instances)
        18
        sage: instances[0]
        GeMSS with D=35, n=35, Δ=2, v=1
    """
    gemss_instances = []

    for D in range(min_D, max_D + 1):
        for n in range(min_n, max_n + 1):
            for delta in range(min_delta, max_delta + 1):
                for v in range(min_v, max_v + 1):
                    G = GeMSS(D=D, n=n, delta=delta, v=v)
                    if G.security_level_classical() >= sec_level_classical and \
                            G.security_level_quantum() >= sec_level_quantum:
                        gemss_instances.append(G)

    return gemss_instances
