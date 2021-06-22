from sage.arith.misc import is_prime_power
from sage.functions.log import log
from sage.functions.other import ceil


def ngates(q, n):
    """
    Return the number of gates for the given number of multiplications in a finite field

    INPUT:

    - ``q`` -- order of the finite field
    - ``n`` -- no. of multiplications

    EXAMPLES::

        sage: from mpkc.utils import ngates
        sage: ngates(16, 2**16)
        2359296

    TESTS::

        sage: ngates(6, 2**16)
        Traceback (most recent call last):
        ...
        ValueError: q must be a prime power
    """
    if not is_prime_power(q):
        raise ValueError("q must be a prime power")
    return n * (2 * log(q, 2) ** 2 + log(q, 2))


def nbits(q, n):
    """
    Return the number of bits required to store `n` elements of a finite field

    - ``q`` -- order of the finite field
    - ``n`` -- no. of field elements

    EXAMPLES::

        sage: from mpkc.utils import nbits
        sage: nbits(4, 256)
        512
    """
    return ceil(log(q, 2)) * n
