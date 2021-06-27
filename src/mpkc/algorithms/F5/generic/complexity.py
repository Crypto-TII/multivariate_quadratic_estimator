from sage.functions.other import binomial
from mpkc import degree_of_regularity


def time(n, degrees, w=2):
    """
    Return the time complexity of the F5 algorithm

    INPUT:

    - ``n`` -- no. of variables
    - ``degrees`` -- a list of integers representing the degree of the polynomials
    - ``w`` -- linear algebra constant (2 <= w <= 3) (default: 2)

    EXAMPLES::

        sage: from mpkc.algorithms import F5
        sage: F5.generic.complexity.time(10, [3]*15)
        8533694884
    """
    m = len(degrees)
    if n >= m:
        complexity = regular_system(n, degrees, w)
    else:
        complexity = semi_regular_system(n, degrees, w)

    return complexity


def regular_system(n, degrees, w=2):
    """
    Return the complexity of F5 algorithm for regular system

    INPUT:

    - ``n`` -- no. of variables
    - ``degrees`` -- a list of integers representing the degree of the polynomials
    - ``w`` -- linear algebra constant (2 <= w <= 3) (default: 2)

    EXAMPLES::

        sage: from mpkc.algorithms import F5
        sage: F5.generic.complexity.regular_system(n=10, degrees=[2]*5)
        626250625
        sage: F5.generic.complexity.regular_system(n=15, degrees=[2]*5)
        37558440000
    """
    dreg = degree_of_regularity.regular_system(n, degrees)
    m = len(degrees)
    return (m * binomial(n + dreg - 1, dreg)) ** w


def semi_regular_system(n, degrees, w=2):
    """
    Return the complexity of F5 algorithm for semi-regular system

    INPUT:

    - ``n`` -- no. of variables
    - ``degrees`` -- a list of integers representing the degree of the polynomials
    - ``w`` -- linear algebra constant (2 <= w <= 3) (default: 2)

    EXAMPLES::

        sage: from mpkc.algorithms import F5
        sage: F5.generic.complexity.semi_regular_system(n=5, degrees=[2]*10)
        3136
        sage: F5.generic.complexity.semi_regular_system(n=5, degrees=[2]*15)
        441
    """
    dreg = degree_of_regularity.semi_regular_system(n, degrees)
    return binomial(n + dreg, dreg) ** w
