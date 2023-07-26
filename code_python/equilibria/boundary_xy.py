"""
Author: Rayla Kurosaki

GitHub: https://github.com/rkp1503
"""

import numpy as np


def get_equilibria(parameter_values):
    equilibria = get_analytical_solutions(parameter_values)
    r_yx, r_zx, p, phi_xy, phi_yx, phi_xz, u_1, u_2, u_3, u_4 = parameter_values
    beta = 11
    for (i, equilibrium) in enumerate(equilibria):
        x, y, z = equilibrium
        # if not (x > 0):
        #     del equilibria[i]
        #     pass
        cond = phi_yx < (beta - 1) / ((phi_xy * beta ** 2 + 1) ** 2)
        if not all([cond]):
            del equilibria[i]
            pass
        pass
    return equilibria


def get_analytical_solutions(parameter_values):
    sols_real = []
    r_yx, r_zx, p, phi_xy, phi_yx, phi_xz, u_1, u_2, u_3, u_4 = parameter_values
    Y_0 = phi_yx + 1
    Y_1 = -1
    Y_2 = 2 * phi_xy * phi_yx
    Y_3 = 0
    Y_4 = phi_xy ** 2 * phi_yx
    for root in np.roots([Y_4, Y_3, Y_2, Y_1, Y_0]):
        if (root.real > 0) and (root.imag == 0):
            y = root.real
            x = 1 + phi_xy * y ** 2
            sols_real.append((x, y, 0))
            pass
        pass
    return sols_real


def get_stable_equilibria(equilibria, parameter_values):
    r_yx, r_zx, p, phi_xy, phi_yx, phi_xz, u_1, u_2, u_3, u_4 = parameter_values
    for (i, equilibrium) in enumerate(equilibria):
        x, y, z = equilibrium
        j_11 = 1 - 2 * x + phi_xy * y ** 2
        j_12 = 2 * phi_xy * x * y
        j_21 = 2 * r_yx * phi_yx * x * y
        j_22 = r_yx * (1 - 2 * y + phi_yx * x ** 2)
        j_33 = r_zx - u_4 + ((u_3 * (1 - p) * y) / (u_2 + (1 - p) * y))
        C_2 = -(j_11 + j_22 + j_33)
        C_1 = j_11 * j_22 + j_11 * j_33 + j_22 * j_33 - j_12 * j_21
        C_0 = j_33*(j_12*j_21-j_11*j_22)
        cond1 = C_2 > 0
        cond2 = C_1 > 0
        cond3 = C_0 > 0
        cond4 = C_2 * C_1 > C_0
        if not all([cond1, cond2, cond3, cond4]):
            del equilibria[i]
            pass
        pass
    return equilibria
