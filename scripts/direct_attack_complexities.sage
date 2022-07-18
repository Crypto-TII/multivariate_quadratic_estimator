# *****************************************************************************
# Multivariate Quadratic (MQ) Estimator
# Copyright (C) 2021-2022 Emanuele Bellini, Rusydi H. Makarim, Javier Verbel
# Cryptography Research Centre, Technology Innovation Institute LLC
#
# This file is part of MQ Estimator
#
# MQ Estimator is free software: you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.
#
# MQ Estimator is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with
# MQ Estimator. If not, see <https://www.gnu.org/licenses/>.
# *****************************************************************************


from mpkc import MQEstimator
from mpkc.utils import ngates, nbits

rainbow_parameters = {
    'Rainbow-I': (16, 100, 64),
    'Rainbow-III': (256, 148, 80),
    'Rainbow-V': (256, 196, 100),
}

mqdss_parameters = {
    'MQDSS-I': (4, 88, 88),
    'MQDSS-III': (4, 128, 128),
    'MQDSS-V': (4, 160, 160),
}

gemss_parameters = {
    'GeMSS-I': (2, 186, 162),
    'GeMSS-III': (2, 258, 243),
    'GeMSS-V': (2, 387, 324),
}

mayo_parameters = {
    'MAYO-Ia': (7, 962, 76),
    'MAYO-Ib': (7, 1140, 78),
    'MAYO-IIIa': (7, 2220, 117),
    'MAYO-IIIb':  (7, 2240, 117),
    'MAYO-Va':  (7, 2960, 152),
    'MAYO-Vb': (7, 3874, 157),
}

theta=2
def security_direct_attack(scheme_parameters, maxD):
    print(70 * '-')
    print(f' Scheme \t Security \t Algorithm \t Optimal parameters')
    print(70 * '-')
    for name in scheme_parameters:
        q, n, m = scheme_parameters[name]
        E = MQEstimator(q=q, n=n, m=m, w=2.81)
        E.crossbred.max_D = maxD
        F = E.fastest_algorithm()
        alg_name = F.__class__.__name__
        time= ngates(q=q, n=F.time_complexity(), theta=theta)
        if F.has_optimal_parameter():
            print('{} \t {:.3f} \t {} \t {}'.format(name, float(log(time,2)), alg_name, F.optimal_parameters()))
        else:
            print('{} \t {:.3f} {} '.format(name, float(log(time,2)), alg_name))
    print(70 * '-')
    print('')

security_direct_attack(rainbow_parameters, 30)
security_direct_attack(mqdss_parameters, 35)
security_direct_attack(mayo_parameters, 45)
security_direct_attack(gemss_parameters, 40)


