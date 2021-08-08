"""
Module to compute the time and memory complexity of MHT (Miura-Hashimoto-Takagi) algorithm

The MHT is an algorithm to solve the MQ problem when  m * (m + 3) / 2 <= n

[MHT13] Miura, H., Hashimoto, Y., and Takagi, T. Extended algorithm for solving underdefined multivariate quadratic
equations. In Post-Quantum Cryptography, 2013. Springer Berlin Heidelberg.
"""

from sage.arith.misc import is_prime_power, is_power_of_two
from math import log2


def time(q, n, m, w=2):
    """
    Return the time complexity of the MHT algorithm

    INPUT:

    - ``q`` -- order of the finite field
    - ``n`` -- no. of variables
    - ``m`` -- no. of polynomials
    - ``w`` -- linear algebra constant (2 <= w <= 3) (default: 2)

    EXAMPLES::

        sage: from mpkc.algorithms import mht
        sage: complexity = mht.complexity.time(q=2, n=183, m=12, w=2.8)
        sage: log(complexity, 2)
        24.6289220479165
    """
    if not is_prime_power(q):
        raise ValueError("the order of finite field q must be a prime power")

    if not m * (m + 3) / 2 <= n:
        raise ValueError(f'The parameter n should be grater than or equal to m * (m + 3) / 2')

    if is_power_of_two(q):
        time = n ** w * m
    else:
        time = 2 ** m * n ** w * m

    return time


def memory(n, m):
    """
    Return the memory complexity of the MHT algorithm

    INPUT:

    - ``n`` -- no. of variables
    - ``m`` -- no. of polynomials

    EXAMPLES::

        sage: from mpkc.algorithms import mht
        sage: mht.complexity.memory(183, 12)
        401868
    """
    if not (m * (m + 3)) / 2 <= n:
        raise ValueError(f'The parameter n should be grater than or equal to m * (m + 3) / 2')

    return m * n**2
