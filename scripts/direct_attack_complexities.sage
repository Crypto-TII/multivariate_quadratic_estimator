from mpkc import MQEstimator

rainbow_parameters = {
    'Rainbow-I': (16, 100, 64),
    'Rainbow-III': (256, 148, 80),
    'Rainbow-V': (256, 148, 80),
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


def security_direct_attack(scheme_parameters):
    print(70 * '-')
    print(f' Scheme \t Security \t Algorithm \t Optimal parameters')
    print(70 * '-')
    for name in scheme_parameters:
        q, n, m = scheme_parameters[name]
        E = MQEstimator(q=q, n=n, m=m, w=2.81)
        F = E.fastest_algorithm()
        alg_name = F.__class__.__name__
        if F.has_optimal_parameter():
            print('{} \t {:.3f} \t {} \t {}'.format(name, float(log(F.time_complexity(),2)), alg_name, F.optimal_parameters()))
        else:
            print('{} \t {:.3f} {} '.format(name, float(log(F.time_complexity(),2)), alg_name))
    print(70 * '-')
    print('')

security_direct_attack(rainbow_parameters)
security_direct_attack(gemss_parameters)
security_direct_attack(mayo_parameters)


