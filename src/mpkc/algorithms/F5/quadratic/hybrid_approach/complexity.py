from mpkc.algorithms import F5


def time(q, n, m, w=2, use_quantum=False):
    """
    Return the complexity of hybrid approach on quadratic system of equations using F5 algorithm

    INPUT:

    - ``q`` -- order of the finite field
    - ``n`` -- no. of variables
    - ``m`` -- no. of polynomials
    - ``w`` -- linear algebra constant (2 <= w <= 3) (default: 2)
    - ``use_quantum`` -- return the complexity using quantum computer (default: False)

    EXAMPLES::

        sage: from mpkc.algorithms import F5
        sage: F5.quadratic.hybrid_approach.complexity.time(256, 10, 15)
        1002001
    """
    return F5.generic.hybrid_approach.complexity.time(q, n, [2]*m, w=w, use_quantum=use_quantum)
