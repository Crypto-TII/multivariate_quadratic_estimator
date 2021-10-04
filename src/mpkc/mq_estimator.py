import inspect
from prettytable import PrettyTable
from sage.functions.log import log
from .algorithms.base import BaseAlgorithm
from .algorithms import HybridF5


class MQEstimator(object):
    def __init__(self, n, m, q=None, w=2, nsolutions=1):
        """
        Construct an instance of MQ Estimator

        INPUT:

        - ``n`` -- no. of variables
        - ``m`` -- no. of polynomials
        - ``q`` -- order of the finite field
        - ``w`` -- linear algebra constant (default: 2)
        - ``nsolutions`` -- no. of solutions (default: 1)

        EXAMPLES::

            sage: from mpkc import MQEstimator
            sage: MQEstimator(n=10, m=5)
            MQ Estimator for system with 10 variables and 5 equations
        """
        constructor_args = locals()

        self._algorithms = []
        for Algorithm in BaseAlgorithm.__subclasses__():
            alg_constructor_args = inspect.getargs(Algorithm.__init__.__code__).args
            arg_and_values = {arg: constructor_args[arg] for arg in alg_constructor_args
                              if arg in constructor_args and arg != 'self'}

            try:
                algorithm = Algorithm(**arg_and_values)
            except ValueError:
                continue

            self._algorithms.append(algorithm)
            setattr(self, algorithm.__module__.split('.')[-1], algorithm)

    def algorithms(self):
        """
        Return a list of considered algorithms

        EXAMPLES::

            sage: from mpkc import MQEstimator
            sage: E = MQEstimator(n=10, m=15)
            sage: E.algorithms()
            [Complexity estimator for F5 with 10 variables and 15 polynomials,
             Complexity estimator for hybrid approach with 10 variables and 15 polynomials,
             Dinur's first estimator for the MQ problem,
             Dinur's second estimator for the MQ problem,
             Exhaustive search estimator for the MQ problem,
             Björklund et al.'s estimator for the MQ problem,
             Lokshtanov et al.'s estimator for the MQ problem,
             BooleanSolve and FXL estimators for the MQ problem,
             Crossbred estimator for the MQ problem]
        """
        return self._algorithms

    def algorithm_names(self):
        """
        Return a list of the name of considered algorithms

        EXAMPLES::

            sage: from mpkc import MQEstimator
            sage: E = MQEstimator(n=10, m=15)
            sage: E.algorithm_names()
            ['F5',
             'HybridF5',
             'DinurFirst',
             'DinurSecond',
             'ExhaustiveSearch',
             'Bjorklund',
             'Lokshtanov',
             'BooleanSolveFXL',
             'Crossbred']
        """
        return [algorithm.__class__.__name__ for algorithm in self.algorithms()]

    def nalgorithms(self):
        """
        Return the number of considered algorithms

        EXAMPLES::

            sage: from mpkc import MQEstimator
            sage: E0 = MQEstimator(n=10, m=15)
            sage: E0.nalgorithms()
            9
            sage: E1 = MQEstimator(n=183, m=12, q=4)
            sage: E1.nalgorithms()
            11
        """
        return len(self.algorithms())

    def table(self, use_tilde_o_time=False):
        """
        Return the table describing the complexity of each algorithm and its optimal parameters

        INPUT:

        - ``use_tilde_o_time`` -- use Ō time complexity (default: False)

        EXAMPLES::

            sage: from mpkc import MQEstimator
            sage: E = MQEstimator(n=15, m=15, q=2)
            sage: table = E.table()
            sage: print(table)
            +------------------+------------------+------------------+---------------------------+
            |    algorithm     |       time       |      memory      |         parameters        |
            +------------------+------------------+------------------+---------------------------+
            |        F5        | 62.0451351868504 | 23.1586318751600 |                           |
            |     HybridF5     | 17.1699250014423 | 3.90689059560852 |           k: 14           |
            |    DinurFirst    | 23.6655729769094 | 20.4938554492408 |      λ: 9/14, κ: 3/14     |
            |   DinurSecond    | 20.3499998825216 | 15.8017083589165 |           n1: 2           |
            | ExhaustiveSearch | 17.9660208563962 | 11.7206717868256 |                           |
            |    Bjorklund     | 42.4516669331353 | 15.3160789459123 |           λ: 1/5          |
            |    Lokshtanov    | 67.1234362737997 | 16.1050592581276 |          δ: 1/15          |
            | BooleanSolveFXL  | 20.3398500028846 | 5.82580271452019 | k: 14, variant: las_vegas |
            |    Crossbred     | 17.9309653178356 | 8.98013957763916 |      D: 3, k: 7, d: 1     |
            +------------------+------------------+------------------+---------------------------+
        """
        table = PrettyTable()
        table.field_names = ['algorithm', 'time', 'memory', 'parameters']

        for algorithm in self.algorithms():
            name = algorithm.__class__.__name__
            time_complexity = algorithm.tilde_o_time() if use_tilde_o_time else algorithm.time_complexity()
            memory_complexity = algorithm.memory_complexity()
            optimal_parameters = ', '.join([f"{k}: {v}" for k, v in algorithm.optimal_parameters().items()])

            table.add_row([name,
                           log(time_complexity, 2).numerical_approx(),
                           log(memory_complexity, 2).numerical_approx(),
                           optimal_parameters])

        return table

    def fastest_algorithm(self, use_tilde_o_time=False):
        """
         Return the algorithm with the smallest time complexity

         INPUT:

         - ``use_tilde_o_time`` -- use Ō time complexity (default: False)

         EXAMPLES::

             sage: from mpkc import MQEstimator
             sage: E = MQEstimator(n=15, m=15, q=2)
             sage: E.fastest_algorithm()
             Complexity estimator for hybrid approach with 15 variables and 15 polynomials
         """
        key = lambda algorithm: algorithm.tilde_o_time() if use_tilde_o_time else algorithm.time_complexity()
        return min(self.algorithms(), key=key)

    def __repr__(self):
        algorithm = self.algorithms()[0]
        n = algorithm.nvariables()
        m = algorithm.npolynomials()
        return f"MQ Estimator for system with {n} variables and {m} equations"


def min_npolynomials(security_level, q, w=2):
    """
    Return a minimum number of equations in a determined system that satisfies the given security level

    INPUT:

    - ``security_level`` -- the intended security level (in bits) (80/100/128/192/256)
    - ``q`` -- order of the finite field
    - ``w`` -- linear algebra constant (2 <= w <= 3) (default: 2)

    EXAMPLES::

        sage: from mpkc.mq_estimator import min_npolynomials
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

        sage: from mpkc.mq_estimator import min_nvariables
        sage: min_nvariables(security_level=80, q=16)
        33
    """
    return min_npolynomials(security_level, q, w)
