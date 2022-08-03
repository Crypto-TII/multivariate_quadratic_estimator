from .rainbow import Rainbow


class UOV:
    def __init__(self, q, n, m):
        """
        Construct an instance of UOV

        - ``q`` -- order of the finite field
        - ``n`` -- no. of variables
        - ``m`` -- no. of polynomials (also equal to the no. of oil variables)

        EXAMPLES::

            sage: from mpkc.schemes import UOV
            sage: U = UOV(q=16, n=25, m=10)
            sage: U
            UOV signature over GF(16) with 25 variables and 10 polynomials
        """
        noil_vars = m
        nvinegar_vars = n - noil_vars

        if not nvinegar_vars > noil_vars:
            raise ValueError("the no. of vinegar variables must be greater than the no. of oil variables")

        self._rainbow = Rainbow(q=q, n=n, v=[nvinegar_vars])

    @property
    def base_field(self):
        """
        Return the base field

        EXAMPLES::

            sage: from mpkc.schemes import UOV
            sage: U = UOV(q=31, n=25, m=10)
            sage: U.base_field
            Finite Field of size 31
        """
        return self._rainbow.base_field

    @property
    def nvariables(self):
        """
        Return the no. of variables

        EXAMPLES::

            sage: from mpkc.schemes import UOV
            sage: U = UOV(q=31, n=25, m=10)
            sage: U.nvariables
            25
        """
        return self._rainbow.nvariables

    @property
    def npolynomials(self):
        """
        Return the no. of polynomials

        EXAMPLES::

            sage: from mpkc.schemes import UOV
            sage: U = UOV(q=31, n=25, m=10)
            sage: U.npolynomials
            10
        """
        return self._rainbow.npolynomials

    def inner_affine_map(self):
        """
        Return the inner affine map

        EXAMPLES::

            sage: from mpkc.schemes import UOV
            sage: U = UOV(q=31, n=5, m=2)
            sage: M, v = U.inner_affine_map()
            sage: M  # random
            [30 10  8  4 11]
            [22 20  5 30 19]
            [27 24 24 18 27]
            [ 3  5 21  0 13]
            [ 8 12  4  0 13]
            sage: v  # random
            (2, 12, 7, 13, 0)
        """
        return self._rainbow.inner_affine_map()

    def outer_affine_map(self):
        """
        Return the outer affine map

        EXAMPLES::

            sage: from mpkc.schemes import UOV
            sage: U = UOV(q=31, n=5, m=2)
            sage: M, v = U.outer_affine_map()
            sage: M  # random
            [17  8]
            [12 15]
            sage: v  # random
            (7, 14)
        """
        return self._rainbow.outer_affine_map()

    def ring(self):
        """
        Return the polynomial ring

        EXAMPLES::

            sage: from mpkc.schemes import UOV
            sage: U = UOV(q=31, n=5, m=2)
            sage: U.ring()
            Multivariate Polynomial Ring in x0, x1, x2, x3, x4 over Finite Field of size 31
        """
        return self._rainbow.ring()

    def vars(self):
        """
        Return a tuple of variables in the polynomial ring

        EXAMPLES::

            sage: from mpkc.schemes import UOV
            sage: U = UOV(q=31, n=5, m=2)
            sage: U.vars()
            (x0, x1, x2, x3, x4)
        """
        return self._rainbow.vars()

    def nvinegar_vars(self):
        """
        Return the number of vinegar variables

        EXAMPLES::

            sage: from mpkc.schemes import UOV
            sage: U = UOV(q=31, n=5, m=2)
            sage: U.nvinegar_vars()
            3
        """
        return self._rainbow.nvinegar_vars_at_layer(0)

    def noil_vars(self):
        """
        Return the number of oil variables

        EXAMPLES::

            sage: from mpkc.schemes import UOV
            sage: U = UOV(q=31, n=5, m=2)
            sage: U.noil_vars()
            2
        """
        return self._rainbow.noil_vars_at_layer(0)

    def vinegar_vars(self):
        """
        Return a list of vinegar variables

        EXAMPLES::

            sage: from mpkc.schemes import UOV
            sage: U = UOV(q=31, n=5, m=2)
            sage: U.vinegar_vars()
            [x0, x1, x2]
        """
        return self._rainbow.vinegar_vars_at_layer(0)

    def oil_vars(self):
        """
        Return a list of oil variables

        EXAMPLES::

            sage: from mpkc.schemes import UOV
            sage: U = UOV(q=31, n=5, m=2)
            sage: U.oil_vars()
            [x3, x4]
        """
        return self._rainbow.oil_vars_at_layer(0)

    def central_map(self):
        """
        Return a list of polynomials representing the central map

        EXAMPLES::

            sage: from mpkc.schemes import UOV
            sage: U = UOV(q=5, n=5, m=2)
            sage: U.central_map()  # random
            [-2*x0*x2 - x1*x2 + 2*x2^2 + x0*x3 + 2*x2*x3 + x0*x4 - 2*x1*x4 + 2*x2*x4 - x0 + 2*x1 - 2*x2 - x3 - 2*x4 + 1,
            x0^2 + x1^2 + x0*x2 - 2*x1*x2 - 2*x2^2 - x0*x3 - x1*x3 + 2*x2*x3 - x1*x4 + x2*x4 - x0 - x1 + x2 + x4 - 2]
        """
        return self._rainbow.central_map()

    def eval_central_map(self, v):
        """
        Return the output of evaluation of the central map on `v`

        INPUT:

        - ``v`` -- a list of `self.base_field` elements

        EXAMPLES::

            sage: from mpkc.schemes import UOV
            sage: U = UOV(q=31, n=7, m=3)
            sage: v = U.random_vector(7)
            sage: U.eval_central_map(v)  # random
            (18, 13, 11)
        """
        return self._rainbow.eval_central_map(v)

    def preimage_central_map(self, y):
        """
        Return a preimage of vector `y`

        INPUT:

        - ``y`` -- a list of `self.base_field` elements

        EXAMPLES::

            sage: from mpkc.schemes import UOV
            sage: U = UOV(q=31, n=7, m=3)
            sage: y = vector(U.base_field, [2, 1, 0])
            sage: x = U.preimage_central_map(y)
            sage: x  # random
            (12, 5, 24, 0, 6, 29, 27)
        """
        return self._rainbow.preimage_central_map(y)

    def inner_affine_polynomials(self):
        """
        Return a list of polynomials representing the inner affine map

        EXAMPLES::

            sage: from mpkc.schemes import UOV
            sage: U = UOV(q=5, n=7, m=3)
            sage: U.inner_affine_polynomials()  # random
            [x0 + x1 - 2*x3 - 2*x4 - x5 - 1,
            -x0 + x1 - x2 + x4 + 2*x5 + x6 - 2,
            -2*x1 - 2*x2 - 2*x3 + x4 + 2*x6,
            x2 - x3 + x4 - x6 - 1,
            -x0 - x1 + 2*x3 + x5 - 1,
            x0 - 2*x1 + 2*x3 - x4 - 2*x6 - 2,
            -x1 + 2*x2 - 2*x3 + 2*x4 - x5 + 2*x6 - 1]
        """
        return self._rainbow.inner_affine_polynomials()

    def random_vector(self, n):
        """
        Return a random vector of length `n` over `self.base_field`

        INPUT:

        - ``n`` -- a positive integer

        EXAMPLES::

            sage: from mpkc.schemes import UOV
            sage: U = UOV(q=5, n=7, m=3)
            sage: U.random_vector(3)  # random
            (1, 4, 2)
        """
        return self._rainbow.random_vector(n)

    def random_msg(self):
        """
        Return a random vector representing a message

        EXAMPLES::

            sage: from mpkc.schemes import UOV
            sage: U = UOV(q=5, n=7, m=3)
            sage: U.random_msg()  # random
            (4, 1, 1)
        """
        return self._rainbow.random_msg()

    def random_signature(self):
        """
        Return a random vector representing a signature

        EXAMPLES::

            sage: from mpkc.schemes import UOV
            sage: U = UOV(q=5, n=7, m=3)
            sage: U.random_signature()  # random
            (1, 2, 0, 4, 3, 2, 1)
        """
        return self._rainbow.random_signature()

    def sign(self, msg):
        """
        Return a signature for the given message

        INPUT:

        - ``msg`` -- a list of `self.base_field` elements

        EXAMPLES::

            sage: from mpkc.schemes import UOV
            sage: U = UOV(q=5, n=7, m=3)
            sage: msg = U.random_msg()
            sage: signature = U.sign(msg)
            sage: signature  # random
            (3, 3, 0, 2, 3, 4, 4)
            sage: U.is_valid_signature(signature, msg)
            True
        """
        return self._rainbow.sign(msg)

    def public_key(self):
        """
        Return a list of polynomials for the public key

        EXAMPLES::

            sage: from mpkc.schemes import UOV
            sage: U = UOV(q=5, n=5, m=2)
            sage: U.public_key()  # random
            [-x0^2 + x0*x1 - x1*x2 + x2^2 + 2*x0*x3 - x1*x3 + x3^2 - x0*x4 - 2*x2*x4 + x3*x4 - x4^2 - x0 - x1 + x4 + 2,
            2*x0^2 + x0*x1 + x1^2 - x2^2 - x0*x3 - x2*x3 - x0*x4 - 2*x2*x4 + x3*x4 + 2*x4^2 + 2*x0 + 2*x1 + 2*x2 - 2*x4]
        """
        return self._rainbow.public_key()

    def inverse_inner_affine_map(self):
        """
        Return the inverse of the inner affine map

        EXAMPLES::

            sage: from mpkc.schemes import UOV
            sage: U = UOV(q=5, n=5, m=2)
            sage: Ti, ti = U.inverse_inner_affine_map()
            sage: Ti  # random
            [0 2 4 1 3]
            [0 2 3 1 3]
            [3 2 3 0 4]
            [0 2 0 3 4]
            [4 4 1 1 1]
            sage: ti  # random
            (1, 0, 1, 3, 2)
        """
        return self._rainbow.inverse_inner_affine_map()

    def inverse_outer_affine_map(self):
        """
        Return the inverse of the outer affine map

        EXAMPLES::

            sage: from mpkc.schemes import UOV
            sage: U = UOV(q=5, n=5, m=2)
            sage: Si, si = U.inverse_outer_affine_map()
            sage: Si  # random
            [3 1]
            [2 2]
            sage: si  # random
            (3, 2)
        """
        return self._rainbow.inverse_outer_affine_map()

    def is_valid_signature(self, signature, msg):
        """
        Return whether the signature is valid for the given message

        INPUT:

        - ``signature`` -- a list of `self.base_field` elements
        - ``msg`` -- a list of `self.base_field` elements

        EXAMPLES::

            sage: from mpkc.schemes import UOV
            sage: U = UOV(q=5, n=5, m=2)
            sage: msg = U.random_msg()
            sage: signature = U.sign(msg)
            sage: U.is_valid_signature(signature, msg)
            True
        """
        return self._rainbow.is_valid_signature(signature, msg)

    def __repr__(self):
        q = self.base_field.order()
        n = self.nvariables
        m = self.npolynomials
        return f"UOV signature over GF({q}) with {n} variables and {m} polynomials"
