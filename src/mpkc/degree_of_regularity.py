def regular_system(n, degrees):
    """
    Return the degree of regularity for regular system

    INPUT:

    - ``n`` -- no. of variables
    - ``degree`` -- a list of integers representing the degree of the polynomials

    EXAMPLES::

        sage: from mpkc.degree_of_regularity import regular_system as dreg_regular
        sage: dreg_regular(15, [2]*10)
        11
    """
    m = len(degrees)
    if n < m:
        raise ValueError("the number of variables must be greater than or equal to the number of polynomials")
    return sum((d - 1) for d in degrees) + 1
