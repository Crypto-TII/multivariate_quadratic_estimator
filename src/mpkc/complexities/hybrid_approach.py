from sage.arith.misc import is_prime_power


def hybrid_approach(q, n, degrees, w=2, use_quantum=False):
    """
    Return the complexity of hybrid approach using F5 algorithm

    INPUT:

    - ``q`` -- order of the finite field
    - ``n`` -- no. of variables
    - ``degrees`` -- a list of integers representing the degree of the polynomials
    - ``w`` -- linear algebra constant (2 <= w <= 3) (default: 2)
    - ``use_quantum`` -- return the complexity using quantum computer (default: False)

    EXAMPLES::

        sage: from mpkc.complexities.hybrid_approach import hybrid_approach
        sage: hybrid_approach(256, 10, [2]*10)
        6412806400
        sage: hybrid_approach(256, 10, [2]*15)
        1002001
    """
    return _hybrid_approach_(q, n, degrees, w, use_quantum, get_tradeoff=False)


def quadratic_system(q, n, m, w=2, use_quantum=False):
    """
    Return the complexity of hybrid approach on quadratic system of equations using F5 algorithm

    INPUT:

    - ``q`` -- order of the finite field
    - ``n`` -- no. of variables
    - ``m`` -- no. of polynomials
    - ``w`` -- linear algebra constant (2 <= w <= 3) (default: 2)
    - ``use_quantum`` -- return the complexity using quantum computer (default: False)

    EXAMPLES::

        sage: from mpkc.complexities.hybrid_approach import quadratic_system
        sage: quadratic_system(256, 10, 15)
        1002001
    """
    return hybrid_approach(q, n, [2]*m, w=w, use_quantum=use_quantum)


def best_tradeoff(q, n, degrees, w=2, use_quantum=False):
    """
    Return the best tradeoff, i.e. optimal no. of fixed variables

    INPUT:

    - ``q`` -- order of the finite field
    - ``n`` -- no. of variables
    - ``degrees`` -- a list of integers representing the degree of the polynomials
    - ``w`` -- linear algebra constant (2 <= w <= 3) (default: 2)
    - ``use_quantum`` -- return the complexity using quantum computer (default: False)

    EXAMPLES::

        sage: from mpkc.complexities.hybrid_approach import best_tradeoff
        sage: best_tradeoff(31, 23, [2]*23)
        2
        sage: best_tradeoff(256, 10, [2]*10)
        1
    """
    return _hybrid_approach_(q, n, degrees, w, use_quantum, get_tradeoff=True)


def _hybrid_approach_(q, n, degrees, w, use_quantum, get_tradeoff):
    """
    Return either the complexity or the best tradeoff of the hybrid approach

    INPUT:

    - ``q`` -- order of the finite field
    - ``n`` -- no. of variables
    - ``degrees`` -- a list of integers representing the degree of the polynomials
    - ``w`` -- linear algebra constant (2 <= w <= 3)
    - ``use_quantum`` -- return the complexity using quantum computer
    - ``get_tradeoff`` -- whether returning the best tradeoff of the hybrid approach
    """
    if not is_prime_power(q):
        raise ValueError("q must be a prime power")

    m = len(degrees)
    if n > m:
        n -= (n - m)

    if use_quantum:
        def c(x, e): return x ** (e / 2)
    else:
        def c(x, e): return x ** e

    from .F5 import F5

    min_complexity = float('inf')
    nfixed_variables = 0
    for k in range(n):
        complexity = c(q, k) * F5(n - k, degrees, w)
        if complexity < min_complexity:
            min_complexity = complexity
            nfixed_variables = k

    if get_tradeoff:
        ret = nfixed_variables
    else:
        ret = min_complexity

    return ret
