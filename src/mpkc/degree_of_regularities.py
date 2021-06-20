from .hilbert_series import HilbertSeries


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

    S = HilbertSeries(n, degrees)
    return S.first_nonpositive_integer()


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
    if n >= m:
        dreg = regular_system(n, [2]*m)
    else:
        dreg = semi_regular_system(n, [2]*m)

    return dreg
