from sage.arith.misc import is_prime_power


def time(q, n, degrees, w=2, use_quantum=False):
    """
    Return the complexity of hybrid approach using F5 algorithm

    INPUT:

    - ``q`` -- order of the finite field
    - ``n`` -- no. of variables
    - ``degrees`` -- a list of integers representing the degree of the polynomials
    - ``w`` -- linear algebra constant (2 <= w <= 3) (default: 2)
    - ``use_quantum`` -- return the complexity using quantum computer (default: False)

    EXAMPLES::

        sage: from mpkc.algorithms import F5
        sage: F5.generic.hybrid_approach.complexity.time(256, 10, [2]*10)
        6412806400
        sage: F5.generic.hybrid_approach.complexity.time(256, 10, [2]*15)
        1002001
    """
    min_finder = lambda iterable: min(iterable)
    return _hybrid_approach_(q, n, degrees, w, use_quantum, min_finder)


def _hybrid_approach_(q, n, degrees, w, use_quantum, min_finder):
    """
    Return either the complexity or the optimal no. of fixed variables for the hybrid approach

    INPUT:

    - ``q`` -- order of the finite field
    - ``n`` -- no. of variables
    - ``degrees`` -- a list of integers representing the degree of the polynomials
    - ``w`` -- linear algebra constant (2 <= w <= 3)
    - ``use_quantum`` -- return the complexity using quantum computer
    - ``min_finder`` -- a function to find the "minimum" element
    """
    if not is_prime_power(q):
        raise ValueError("q must be a prime power")

    m = len(degrees)
    if n > m:
        n -= (n - m)

    if use_quantum:
        def c(x, e):
            return x ** (e / 2)
    else:
        def c(x, e):
            return x ** e

    from mpkc.algorithms import F5

    complexities = [c(q, k) * F5.generic.complexity.time(n - k, degrees, w) for k in range(n)]
    return min_finder(complexities)
