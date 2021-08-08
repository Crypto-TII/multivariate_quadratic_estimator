"""
Module to compute the time and memory complexity of kpg (Kipnis, Patarin, Goubin) algorithm

The kpg is an algorithm to solve a quadratic systems of equations over fields of even characteristic

[KPG99] Kipnis, A., Patarin, J.,  and  Goubin, L. Unbalanced  oil  andvinegar signature schemes. In Advances in
Cryptology —EUROCRYPT99, pages 206–222, Berlin, Heidelberg, 1999. Springer BerlinHeidelberg.
"""

from sage.arith.misc import is_power_of_two


def time(q, n, m, w=2):
    """
    Return the time complexity of the kpg algorithm

    INPUT:

    - ``q`` -- order of the finite field
    - ``n`` -- no. of variables
    - ``m`` -- no. of polynomials
    - ``w`` -- linear algebra constant (2 <= w <= 3) (default: 2)

    EXAMPLES::

        sage: from mpkc.algorithms import kpg
        sage: complexity = kpg.complexity.time(q=4, n=183, m=12, w=2.8)
        sage: log(complexity, 2)
        24.6289220479165
    """
    if not is_power_of_two(q):
        raise ValueError("the order of finite field q must be a power of 2")

    if not m * (m + 1) < n:
        raise ValueError(f'The condition m(m + 1) < n must be satisfied')

    return m * n**w


def memory(q, n, m):
    """
    Return the memory complexity of the kpg algorithm

    INPUT:

    - ``q`` -- order of the finite field
    - ``n`` -- no. of variables
    - ``m`` -- no. of polynomials

    EXAMPLES::

        sage: from mpkc.algorithms import kpg
        sage: kpg.complexity.memory(q=4, n=183, m=12)
        401868
    """
    if not is_power_of_two(q):
        raise ValueError("the order of finite field q must be a power of 2")

    if not m * (m + 1) < n:
        raise ValueError(f'The condition m(m + 1) < n must be satisfied')

    return m * n**2
