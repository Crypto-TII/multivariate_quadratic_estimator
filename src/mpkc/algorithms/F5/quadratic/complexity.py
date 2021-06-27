from mpkc.algorithms import F5


def time(n, m, w=2):
    """
    Return the complexity of F5 algorithm for quadratic system

    INPUT:

    - ``n`` -- no. of variables
    - ``m`` -- no. of polynomials
    - ``w`` -- linear algebra constant (2 <= w <= 3) (default: 2)

    EXAMPLES::

        sage: from mpkc.algorithms import F5
        sage: F5.quadratic.complexity.time(5, 10)
        3136
        sage: F5.quadratic.complexity.time(10, 5)
        626250625
    """
    return F5.generic.complexity.time(n, [2]*m, w)
