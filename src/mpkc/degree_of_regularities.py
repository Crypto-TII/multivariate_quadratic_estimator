from .hilbert_series import HilbertSeries


def degree_of_regularity(n, degrees):
    """
    Return the degree of regularity for the system of polynomial equations

    INPUT:

    - ``n`` -- no. of variables
    - ``degrees`` -- a list of integers representing the degree of the polynomials

    EXAMPLES::

        sage: from mpkc import degree_of_regularity
        sage: degree_of_regularity(5, [2]*10)
        3
        sage: degree_of_regularity(10, [3]*5)
        11
    """
    m = len(degrees)
    if n >= m:
        dreg = regular_system(n, degrees)
    else:
        dreg = semi_regular_system(n, degrees)
    return dreg


def regular_system(n, degrees):
    """
    Return the degree of regularity for regular system

    INPUT:

    - ``n`` -- no. of variables
    - ``degree`` -- a list of integers representing the degree of the polynomials

    EXAMPLES::

        sage: from mpkc.degree_of_regularities import regular_system as dreg_regular
        sage: dreg_regular(15, [2]*10)
        11
    """
    m = len(degrees)
    if n < m:
        raise ValueError("the number of variables must be greater than or equal to the number of polynomials")
    return sum((d - 1) for d in degrees) + 1


def semi_regular_system(n, degrees):
    """
    Return the degree of regularity for semi-regular system

    INPUT:

    - ``n`` -- no. of variables
    - ``degrees`` -- a list of integers representing the degree of the polynomials

    EXAMPLES::

        sage: from mpkc.degree_of_regularities import semi_regular_system as dreg_semi_regular
        sage: dreg_semi_regular(10, [2]*15)
        4
    """
    m = len(degrees)
    if m <= n:
        raise ValueError("the number of polynomials must be strictly greater than the number of variables")

    s = HilbertSeries(n, degrees)
    return s.first_nonpositive_integer()


def quadratic_system(n, m):
    """
    Return the degree of regularity for quadratic system

    INPUT:

    - ``n`` -- no. of variables
    - ``m`` -- no. of polynomials

    EXAMPLES::

        sage: from mpkc.degree_of_regularities import quadratic_system
        sage: quadratic_system(10, 15)
        4
        sage: quadratic_system(15, 15)
        16
    """
    return degree_of_regularity(n, [2]*m)
