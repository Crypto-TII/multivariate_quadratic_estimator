from .complexity import _hybrid_approach_


def nfixed_variables(q, n, degrees, w=2, use_quantum=False):
    """
    Return the best tradeoff, i.e. optimal no. of fixed variables

    INPUT:

    - ``q`` -- order of the finite field
    - ``n`` -- no. of variables
    - ``degrees`` -- a list of integers representing the degree of the polynomials
    - ``w`` -- linear algebra constant (2 <= w <= 3) (default: 2)
    - ``use_quantum`` -- return the complexity using quantum computer (default: False)

    EXAMPLES::

        sage: from mpkc.algorithms.F5.generic.hybrid_approach import nfixed_variables
        sage: nfixed_variables(31, 23, [2]*23)
        2
        sage: nfixed_variables(256, 10, [2]*10)
        1
    """
    min_finder = lambda iterable: min(range(len(iterable)), key=iterable.__getitem__)
    return _hybrid_approach_(q, n, degrees, w, use_quantum, min_finder)
