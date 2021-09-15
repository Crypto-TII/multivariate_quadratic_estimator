from sage.functions.log import log
from ..algorithms import HybridF5


def min_npolynomials(security_level, q, w=2):
    """
    Return a minimum number of equations in a determined system that satisfies the given security level

    INPUT:

    - ``security_level`` -- the intended security level (in bits) (80/100/128/192/256)
    - ``q`` -- order of the finite field
    - ``w`` -- linear algebra constant (2 <= w <= 3) (default: 2)

    EXAMPLES::

        sage: from mpkc.estimators.mq import min_npolynomials
        sage: min_npolynomials(security_level=80, q=16)
        33

    TESTS::

        sage: min_npolynomials(security_level=80, q=31)
        32
        sage: min_npolynomials(security_level=80, q=256)
        28
        sage: min_npolynomials(security_level=100, q=16)
        43
        sage: min_npolynomials(security_level=100, q=31)
        40
        sage: min_npolynomials(security_level=100, q=256)
        36
        sage: min_npolynomials(security_level=128, q=16)
        56
        sage: min_npolynomials(security_level=128, q=31)
        52
        sage: min_npolynomials(security_level=128, q=256)
        47
        sage: min_npolynomials(security_level=192, q=16)
        86
        sage: min_npolynomials(security_level=192, q=31)
        80
        sage: min_npolynomials(security_level=192, q=256)
        72
        sage: min_npolynomials(security_level=256, q=16)  # long time
        116
        sage: min_npolynomials(security_level=256, q=31)  # long time
        109
        sage: min_npolynomials(security_level=256, q=256)  # long time
        98
    """
    if security_level not in (80, 100, 128, 192, 256):
        raise ValueError("the valid parameter for security_level is {80, 100, 128, 192, 256}")

    m = 1
    while log(HybridF5(n=m, m=m, q=q, w=w).time_complexity(), 2) < security_level:
        m += 1

    return m


def min_nvariables(security_level, q, w=2):
    """
    Return a minimum number of variables in a determined system that satisfies the given security level

    INPUT:

    - ``security_level`` -- the intended security level (in bits) (80/100/128/192/256)
    - ``q`` -- order of the finite field
    - ``w`` -- linear algebra constant (2 <= w <= 3) (default: 2)

    EXAMPLES::

        sage: from mpkc.estimators.mq import min_nvariables
        sage: min_nvariables(security_level=80, q=16)
        33
    """
    return min_npolynomials(security_level, q, w)
