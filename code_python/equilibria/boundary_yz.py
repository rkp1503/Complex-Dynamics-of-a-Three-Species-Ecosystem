"""
Author: Rayla Kurosaki

GitHub: https://github.com/rkp1503
"""

import numpy as np


def get_equilibria(parameter_values):
    equilibria = get_analytical_solutions(parameter_values)
    r_yx, r_zx, p, phi_xy, phi_yx, phi_xz, u_1, u_2, u_3, u_4 = parameter_values
    for (i, equilibrium) in enumerate(equilibria):
        x, y, z = equilibrium
        # if not (z > 0):
        #     del equilibria[i]
        #     pass
        cond1 = y > (u_2 * (u_4 - r_zx)) / ((u_3 - u_4 + r_zx) * (1 - p))
        cond2 = (r_yx * r_zx * u_2) / (u_1 * (r_zx - u_4) * (1 - p)) > 1
        if not all([cond1, cond2]):
            del equilibria[i]
            pass
        pass
    return equilibria


def get_analytical_solutions(parameter_values):
    sols_real = []
    r_yx, r_zx, p, phi_xy, phi_yx, phi_xz, u_1, u_2, u_3, u_4 = parameter_values
    Y_0 = u_2 * (r_yx * r_zx * u_2 + u_1 * (u_4 - r_zx) * (1 - p))
    Y_1 = r_yx * r_zx * u_2 * (2 * (1 - p) - u_2) + u_1 * (
                u_4 - u_3 - r_zx) * (1 - p) ** 2
    Y_2 = r_yx * r_zx * ((1 - p) - 2 * u_2) * (1 - p)
    Y_3 = -r_yx * r_zx * (1 - p) ** 2
    if Y_0 > 0:
        for root in np.roots([Y_3, Y_2, Y_1, Y_0]):
            if (root.real > 0) and (root.imag == 0):
                y = root.real
                z = 1 + ((((u_3 * (1 - p) * y) / (
                            u_2 + (1 - p) * y)) - u_4) / r_zx)
                sols_real.append((0, y, z))
                pass
            pass
        pass
    return sols_real


def get_stable_equilibria(equilibria, parameter_values):
    r_yx, r_zx, p, phi_xy, phi_yx, phi_xz, u_1, u_2, u_3, u_4 = parameter_values
    for (i, equilibrium) in enumerate(equilibria):
        x, y, z = equilibrium
        j_11 = 1 + phi_xy * y ** 2 - phi_yx * z
        j_22 = r_yx * (1 - 2 * y) - (
                    (u_1 * u_2 * (1 - p) * z) / (u_2 + (1 - p) * y) ** 2)
        j_23 = -((u_1 * (1 - p) * y) / (u_2 + (1 - p) * y))
        j_32 = ((u_2 * u_3 * (1 - p) * z) / (u_2 + (1 - p) * y) ** 2)
        j_33 = r_zx * (1 - 2 * z) - u_4 + (
                    (u_3 * (1 - p) * y) / (u_2 + (1 - p) * y))
        C_2 = -j_11 - j_22 - j_33
        C_1 = j_11 * j_22 + j_11 * j_33 + j_22 * j_33 - j_23 * j_32
        C_0 = j_11 * (j_23 * j_32 - j_22 * j_33)
        cond1 = C_2 > 0
        cond2 = C_1 > 0
        cond3 = C_0 > 0
        cond4 = C_2 * C_1 > C_0
        if not all([cond1, cond2, cond3, cond4]):
            del equilibria[i]
            pass
        pass
    return equilibria
