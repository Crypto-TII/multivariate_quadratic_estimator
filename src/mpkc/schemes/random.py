from sage.rings.finite_rings.all import FiniteField
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.structure.sequence import Sequence
from ..utils import nmonomials_up_to_degree


def random_posso(q, n, m, d):
    """
    Return a random instance of PoSSo (Polynomial System Solving) problem

    INPUT:

    - ``q`` -- order of the finite field
    - ``n`` -- no. of variables
    - ``m`` -- no. of polynomials
    - ``d`` -- degree of the system

    EXAMPLES::

        sage: from mpkc.schemes import random_posso
        sage: F = random_posso(q=31, n=5, m=20, d=3)
        sage: F
        Polynomial Sequence with 20 Polynomials in 5 Variables
        sage: max(f.degree() for f in F)
        3

    TESTS::

        sage: F.groebner_basis() != [1]
        True
    """
    R = PolynomialRing(FiniteField(q), n, 'x', order="degrevlex")
    max_nterms = nmonomials_up_to_degree(q=q, n=n, d=d)
    F = Sequence([R.random_element(degree=d, terms=max_nterms) for _ in range(m)])

    if m >= n:
        # ensure a solution exists for (over)-determined system
        solution = [R.base_ring().random_element() for _ in range(n)]
        x = R.gens()
        y = F.subs(dict(zip(x, solution)))
        for i in range(len(F)):
            F[i] -= y[i]

    return F


def random_mq(q, n, m):
    """
    Return a random instance of random MQ problem

    INPUT:

    - ``q`` -- order of the finite field
    - ``n`` -- no. of variables
    - ``m`` -- no. of polynomials

    EXAMPLES::

        sage: from mpkc.schemes import random_mq
        sage: F = random_mq(q=7, n=10, m=20)
        sage: F
        Polynomial Sequence with 20 Polynomials in 10 Variables
        sage: max(f.degree() for f in F)
        2
    """
    return random_posso(q=q, n=n, m=m, d=2)
