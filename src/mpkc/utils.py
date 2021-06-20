from sage.arith.misc import is_prime_power
from sage.functions.log import log


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
