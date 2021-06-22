from sage.functions.other import binomial
from .. import degree_of_regularities


def F5(n, degrees, w=2):
    """
    Return the complexity of F5 algorithm

    INPUT:

    - ``n`` -- no. of variables
    - ``degrees`` -- a list of integers representing the degree of the polynomials
    - ``w`` -- linear algebra constant (2 <= w <= 3) (default: 2)

    EXAMPLES::

        sage: from mpkc.complexities.F5 import F5
        sage: F5(10, [3]*15)
        8533694884
    """
    m = len(degrees)
    if n >= m:
        complexity = regular_system(n, degrees, w)
    else:
        complexity = semi_regular_system(n, degrees, w)

    return complexity


def quadratic_system(n, m, w=2):
    """
    Return the complexity of F5 algorithm for quadratic system

    INPUT:

    - ``n`` -- no. of variables
    - ``m`` -- no. of polynomials
    - ``w`` -- linear algebra constant (2 <= w <= 3) (default: 2)

    EXAMPLES::

        sage: from mpkc.complexities import F5
        sage: F5.quadratic_system(5, 10)
        3136
        sage: F5.quadratic_system(10, 5)
        626250625
    """
    return F5(n, [2]*m, w)


def regular_system(n, degrees, w=2):
    """
    Return the complexity of F5 algorithm for regular system

    INPUT:

    - ``n`` -- no. of variables
    - ``degrees`` -- a list of integers representing the degree of the polynomials
    - ``w`` -- linear algebra constant (2 <= w <= 3) (default: 2)

    EXAMPLES::

        sage: from mpkc.complexities import F5
        sage: F5.regular_system(n=10, degrees=[2]*5)
        626250625
        sage: F5.regular_system(n=15, degrees=[2]*5)
        37558440000
    """
    dreg = degree_of_regularities.regular_system(n, degrees)
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

        sage: from mpkc.complexities import F5
        sage: F5.semi_regular_system(n=5, degrees=[2]*10)
        3136
        sage: F5.semi_regular_system(n=5, degrees=[2]*15)
        441
    """
    dreg = degree_of_regularities.semi_regular_system(n, degrees)
    return binomial(n + dreg, dreg) ** w
