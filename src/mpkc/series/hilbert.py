from sage.all import ZZ, QQ
from sage.misc.misc_c import prod
from sage.rings.power_series_ring import PowerSeriesRing


class HilbertSeries(object):
    def __init__(self, n, degrees):
        """
        Construct an instance of Hilbert series

        INPUT:

        - ``n`` -- no of variables
        - ``degrees`` -- a list of integers representing the degree of the polynomials

        EXAMPLES::

            sage: from mpkc.series.hilbert import HilbertSeries
            sage: H = HilbertSeries(10, [2]*15)
            sage: H
            Hilbert series for system with 10 variables and 15 polynomials
        """
        self._nvariables = n
        self._degrees = degrees
        self._ring = PowerSeriesRing(QQ, 'z', default_prec=2*len(self._degrees))
        z = self._ring.gen()
        self._series = prod([1 - z ** d for d in degrees]) / (1 - z) ** n

    @property
    def nvariables(self):
        """
        Return the no. of variables

        EXAMPLES::

            sage: from mpkc.series.hilbert import HilbertSeries
            sage: H = HilbertSeries(5, [2]*7)
            sage: H.nvariables
            5
        """
        return self._nvariables

    @property
    def degrees(self):
        """
        Return a list of degrees of the polynomials

        EXAMPLES::

            sage: from mpkc.series.hilbert import HilbertSeries
            sage: H = HilbertSeries(5, [2]*7)
            sage: H.degrees
            [2, 2, 2, 2, 2, 2, 2]
        """
        return self._degrees

    @property
    def precision(self):
        """
        Return the default precision of the series

        EXAMPLES::

            sage: from mpkc.series.hilbert import HilbertSeries
            sage: H = HilbertSeries(5, [2]*7)
            sage: H.precision
            14
        """
        return self.ring.default_prec()

    @property
    def ring(self):
        """
        Return the power series ring

        EXAMPLES::

            sage: from mpkc.series.hilbert import HilbertSeries
            sage: H = HilbertSeries(5, [2]*7)
            sage: H.ring
            Power Series Ring in z over Rational Field
        """
        return self._ring

    @property
    def series(self):
        """
        Return the series

        EXAMPLES::

            sage: from mpkc.series.hilbert import HilbertSeries
            sage: H = HilbertSeries(5, [2]*7)
            sage: H.series
            1 + 5*z + 8*z^2 - 14*z^4 - 14*z^5 + 8*z^7 + 5*z^8 + z^9 + O(z^14)
        """
        return self._series

    @property
    def npolynomials(self):
        """
        Return the no. of polynomials

        EXAMPLES::

            sage: from mpkc.series.hilbert import HilbertSeries
            sage: H = HilbertSeries(10, [2]*15)
            sage: H.npolynomials
            15
        """
        return len(self._degrees)

    def first_nonpositive_integer(self):
        """
        Return the first non-positive integer of the series

        EXAMPLES::

            sage: from mpkc.series.hilbert import HilbertSeries
            sage: H = HilbertSeries(10, [2]*15)
            sage: H.first_nonpositive_integer()
            4
        """
        s = self.series()
        for d in range(self.precision):
            if s[d] <= 0:
                return ZZ(d)
        else:
            raise ValueError("unable to find a nonpositive coefficient in the series")

    def __repr__(self):
        return f"Hilbert series for system with {self.nvariables} variables and {self.npolynomials} polynomials"
