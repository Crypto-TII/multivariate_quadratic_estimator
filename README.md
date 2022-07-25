# Multivariate Quadratic (MQ) Complexity Estimator

This repository is a SageMath package providing functions to estimate to complexity of MQ problem. The result of this
estimator is published in the following paper

> Emanuele Bellini, Rusydi H. Makarim, Carlo Sanna and Javier Verbel. *An Estimator for the Hardness of the MQ
> Problem*. International Conference on Cryptology in Africa (AFRICACRYPT) 2022.

An eprint version is available as

> Cryptology ePrint Archive, Report 2022/708, 2022. https://eprint.iacr.org/2022/708.pdf.


## Quick Start

The easiest way to install the MQ estimator is by using the configured docker image under the folder `docker/`. You
need to first download [Docker Desktop](https://www.docker.com/get-started/) (see
https://docs.docker.com/get-started/overview/ to get an overview of what docker is).

We configure the Docker image with a pre-installed sage (v9.5). To build and run the docker image, simply execute

    $ make rundocker

Once completed, you will be able to use the estimator inside the sage shell.


## Local Installation

If you already have sage installed in your local machine, make sure that the `sage` binary is in your `PATH` environment
variable. To do so, run the following command

    $ export PATH=$PATH:/home/user/sage

by replacing `/home/user/sage` with the path where sage binary is located. After that, simply run

    $ make


## Usage Examples

    sage: from mpkc import MQEstimator
    sage: E = MQEstimator(q=31, m=10, n=15)
    sage: print(E.table())
    +------------------+---------+--------+------------------------------+
    |    algorithm     |   time  | memory |          parameters          |
    +------------------+---------+--------+------------------------------+
    |        F5        |  41.359 | 34.715 |                              |
    |     HybridF5     |  29.532 | 23.104 |             k: 1             |
    | ExhaustiveSearch |  47.965 | 9.965  |                              |
    |      CGMTA       |  47.918 | 11.908 |                              |
    |    Lokshtanov    | 649.473 | 46.684 |           δ: 1/10            |
    | BooleanSolveFXL  |  36.334 | 21.489 | k: 2, variant: deterministic |
    |    Crossbred     |  27.98  | 21.61  |       D: 6, k: 8, d: 3       |
    +------------------+---------+--------+------------------------------+


## Directory Structures

- `docker/` -- Dockerfile to generate the docker image
- `docs/` -- Sphinx-generated documentation
- `scripts/` -- scripts to generate the tables in the paper (requires the installation of the estimator)
- `src/mpkc/` -- source codes of the estimator


## Bugs and Contributions

Please reports bugs or new ideas/features through 
[GitHub issue tracker](https://github.com/Crypto-TII/multivariate_quadratic_estimator/issues). 


## Contributors

- Emanuele Bellini - emanuele.bellini@tii.ae
- Rusydi H. Makarim - rusydi.makarim@tii.ae
- Javier Verbel - javier.verbel@tii.ae


## License

This estimator is licensed under the [GPLv3](https://www.gnu.org/licenses/gpl-3.0.en.html) license.


## Estimators for Other Hard Problems

- Lattice Estimator -- https://github.com/malb/lattice-estimator
- Syndrome Decoding Estimator -- https://github.com/Crypto-TII/syndrome_decoding_estimator