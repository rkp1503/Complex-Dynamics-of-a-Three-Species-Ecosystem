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
        # if not ((x > 0) and (z > 0)):
        #     del equilibria[i]
        #     pass
        cond1 = (1 + phi_xy * y ** 2) / phi_xz > z
        cond2 = y > (u_2 * (u_4 - r_zx)) / ((u_3 - u_4 + r_zx) * (1 - p))
        if not all([cond1, cond2]):
            del equilibria[i]
            pass
        pass
    return equilibria


def get_analytical_solutions(parameter_values):
    sols_real = []
    r_yx, r_zx, p, phi_xy, phi_yx, phi_xz, u_1, u_2, u_3, u_4 = parameter_values
    Y_6 = r_yx * r_zx ** 2 * phi_xy ** 2 * phi_yx * (1 - p) ** 2
    Y_5 = 2 * r_yx * r_zx ** 2 * u_2 * phi_xy ** 2 * phi_yx * (1 - p)
    Y_4 = -2 * p ** 2 * r_yx * r_zx ** 2 * phi_xy * phi_xz * phi_yx + 2 * p ** 2 * r_yx * r_zx ** 2 * phi_xy * phi_yx - 2 * p ** 2 * r_yx * r_zx * u_3 * phi_xy * phi_xz * phi_yx + 2 * p ** 2 * r_yx * r_zx * u_4 * phi_xy * phi_xz * phi_yx + 4 * p * r_yx * r_zx ** 2 * phi_xy * phi_xz * phi_yx - 4 * p * r_yx * r_zx ** 2 * phi_xy * phi_yx + 4 * p * r_yx * r_zx * u_3 * phi_xy * phi_xz * phi_yx - 4 * p * r_yx * r_zx * u_4 * phi_xy * phi_xz * phi_yx + r_yx * r_zx ** 2 * u_2 ** 2 * phi_xy ** 2 * phi_yx - 2 * r_yx * r_zx ** 2 * phi_xy * phi_xz * phi_yx + 2 * r_yx * r_zx ** 2 * phi_xy * phi_yx - 2 * r_yx * r_zx * u_3 * phi_xy * phi_xz * phi_yx + 2 * r_yx * r_zx * u_4 * phi_xy * phi_xz * phi_yx
    Y_3 = -p ** 2 * r_yx * r_zx ** 2 + 4 * p * r_yx * r_zx ** 2 * u_2 * phi_xy * phi_xz * phi_yx - 4 * p * r_yx * r_zx ** 2 * u_2 * phi_xy * phi_yx + 2 * p * r_yx * r_zx ** 2 + 2 * p * r_yx * r_zx * u_2 * u_3 * phi_xy * phi_xz * phi_yx - 4 * p * r_yx * r_zx * u_2 * u_4 * phi_xy * phi_xz * phi_yx - 4 * r_yx * r_zx ** 2 * u_2 * phi_xy * phi_xz * phi_yx + 4 * r_yx * r_zx ** 2 * u_2 * phi_xy * phi_yx - r_yx * r_zx ** 2 - 2 * r_yx * r_zx * u_2 * u_3 * phi_xy * phi_xz * phi_yx + 4 * r_yx * r_zx * u_2 * u_4 * phi_xy * phi_xz * phi_yx
    Y_2 = p ** 2 * r_yx * r_zx ** 2 * phi_xz ** 2 * phi_yx - 2 * p ** 2 * r_yx * r_zx ** 2 * phi_xz * phi_yx + p ** 2 * r_yx * r_zx ** 2 * phi_yx + p ** 2 * r_yx * r_zx ** 2 + 2 * p ** 2 * r_yx * r_zx * u_3 * phi_xz ** 2 * phi_yx - 2 * p ** 2 * r_yx * r_zx * u_3 * phi_xz * phi_yx - 2 * p ** 2 * r_yx * r_zx * u_4 * phi_xz ** 2 * phi_yx + 2 * p ** 2 * r_yx * r_zx * u_4 * phi_xz * phi_yx + p ** 2 * r_yx * u_3 ** 2 * phi_xz ** 2 * phi_yx - 2 * p ** 2 * r_yx * u_3 * u_4 * phi_xz ** 2 * phi_yx + p ** 2 * r_yx * u_4 ** 2 * phi_xz ** 2 * phi_yx + 2 * p * r_yx * r_zx ** 2 * u_2 - 2 * p * r_yx * r_zx ** 2 * phi_xz ** 2 * phi_yx + 4 * p * r_yx * r_zx ** 2 * phi_xz * phi_yx - 2 * p * r_yx * r_zx ** 2 * phi_yx - 2 * p * r_yx * r_zx ** 2 - 4 * p * r_yx * r_zx * u_3 * phi_xz ** 2 * phi_yx + 4 * p * r_yx * r_zx * u_3 * phi_xz * phi_yx + 4 * p * r_yx * r_zx * u_4 * phi_xz ** 2 * phi_yx - 4 * p * r_yx * r_zx * u_4 * phi_xz * phi_yx - 2 * p * r_yx * u_3 ** 2 * phi_xz ** 2 * phi_yx + 4 * p * r_yx * u_3 * u_4 * phi_xz ** 2 * phi_yx - 2 * p * r_yx * u_4 ** 2 * phi_xz ** 2 * phi_yx - 2 * r_yx * r_zx ** 2 * u_2 ** 2 * phi_xy * phi_xz * phi_yx + 2 * r_yx * r_zx ** 2 * u_2 ** 2 * phi_xy * phi_yx - 2 * r_yx * r_zx ** 2 * u_2 + r_yx * r_zx ** 2 * phi_xz ** 2 * phi_yx - 2 * r_yx * r_zx ** 2 * phi_xz * phi_yx + r_yx * r_zx ** 2 * phi_yx + r_yx * r_zx ** 2 + 2 * r_yx * r_zx * u_2 ** 2 * u_4 * phi_xy * phi_xz * phi_yx + 2 * r_yx * r_zx * u_3 * phi_xz ** 2 * phi_yx - 2 * r_yx * r_zx * u_3 * phi_xz * phi_yx - 2 * r_yx * r_zx * u_4 * phi_xz ** 2 * phi_yx + 2 * r_yx * r_zx * u_4 * phi_xz * phi_yx + r_yx * u_3 ** 2 * phi_xz ** 2 * phi_yx - 2 * r_yx * u_3 * u_4 * phi_xz ** 2 * phi_yx + r_yx * u_4 ** 2 * phi_xz ** 2 * phi_yx
    Y_1 = -p ** 2 * r_zx ** 2 * u_1 - p ** 2 * r_zx * u_1 * u_3 + p ** 2 * r_zx * u_1 * u_4 - 2 * p * r_yx * r_zx ** 2 * u_2 * phi_xz ** 2 * phi_yx + 4 * p * r_yx * r_zx ** 2 * u_2 * phi_xz * phi_yx - 2 * p * r_yx * r_zx ** 2 * u_2 * phi_yx - 2 * p * r_yx * r_zx ** 2 * u_2 - 2 * p * r_yx * r_zx * u_2 * u_3 * phi_xz ** 2 * phi_yx + 2 * p * r_yx * r_zx * u_2 * u_3 * phi_xz * phi_yx + 4 * p * r_yx * r_zx * u_2 * u_4 * phi_xz ** 2 * phi_yx - 4 * p * r_yx * r_zx * u_2 * u_4 * phi_xz * phi_yx + 2 * p * r_yx * u_2 * u_3 * u_4 * phi_xz ** 2 * phi_yx - 2 * p * r_yx * u_2 * u_4 ** 2 * phi_xz ** 2 * phi_yx + 2 * p * r_zx ** 2 * u_1 + 2 * p * r_zx * u_1 * u_3 - 2 * p * r_zx * u_1 * u_4 - r_yx * r_zx ** 2 * u_2 ** 2 + 2 * r_yx * r_zx ** 2 * u_2 * phi_xz ** 2 * phi_yx - 4 * r_yx * r_zx ** 2 * u_2 * phi_xz * phi_yx + 2 * r_yx * r_zx ** 2 * u_2 * phi_yx + 2 * r_yx * r_zx ** 2 * u_2 + 2 * r_yx * r_zx * u_2 * u_3 * phi_xz ** 2 * phi_yx - 2 * r_yx * r_zx * u_2 * u_3 * phi_xz * phi_yx - 4 * r_yx * r_zx * u_2 * u_4 * phi_xz ** 2 * phi_yx + 4 * r_yx * r_zx * u_2 * u_4 * phi_xz * phi_yx - 2 * r_yx * u_2 * u_3 * u_4 * phi_xz ** 2 * phi_yx + 2 * r_yx * u_2 * u_4 ** 2 * phi_xz ** 2 * phi_yx - r_zx ** 2 * u_1 - r_zx * u_1 * u_3 + r_zx * u_1 * u_4
    Y_0 = p * r_zx ** 2 * u_1 * u_2 - p * r_zx * u_1 * u_2 * u_4 + r_yx * r_zx ** 2 * u_2 ** 2 * phi_xz ** 2 * phi_yx - 2 * r_yx * r_zx ** 2 * u_2 ** 2 * phi_xz * phi_yx + r_yx * r_zx ** 2 * u_2 ** 2 * phi_yx + r_yx * r_zx ** 2 * u_2 ** 2 - 2 * r_yx * r_zx * u_2 ** 2 * u_4 * phi_xz ** 2 * phi_yx + 2 * r_yx * r_zx * u_2 ** 2 * u_4 * phi_xz * phi_yx + r_yx * u_2 ** 2 * u_4 ** 2 * phi_xz ** 2 * phi_yx - r_zx ** 2 * u_1 * u_2 + r_zx * u_1 * u_2 * u_4
    if Y_0 < 0:
        for root in np.roots([Y_6, Y_5, Y_4, Y_3, Y_2, Y_1, Y_0]):
            if (root.real > 0) and (root.imag == 0):
                y = root.real
                z = 1 + ((((u_3 * (1 - p) * y) / (
                        u_2 + (1 - p) * y)) - u_4) / r_zx)
                x = 1 + phi_xy * y ** 2 - phi_xz * z
                sols_real.append((x, y, z))
                pass
            pass
        pass
    return sols_real


def get_stable_equilibria(equilibria, parameter_values):
    r_yx, r_zx, p, phi_xy, phi_yx, phi_xz, u_1, u_2, u_3, u_4 = parameter_values
    for (i, equilibrium) in enumerate(equilibria):
        x, y, z = equilibrium
        j_11 = 1 - 2 * x + phi_xy * y ** 2 - phi_xz * z
        j_12 = 2 * phi_xy * x * y
        j_13 = -phi_xz * x
        j_21 = 2 * phi_yx * r_yx * x * y
        j_22 = r_yx * (1 - 2 * y + phi_yx * x ** 2) - (
                (u_1 * u_2 * (1 - p) * z) / (u_2 + (1 - p) * y) ** 2)
        j_23 = -((u_1 * (1 - p) * y) / (u_2 + (1 - p) * y))
        j_32 = (u_2 * u_3 * (1 - p) * z) / ((u_2 + (1 - p) * y) ** 2)
        j_33 = r_zx * (1 - 2 * z) - u_4 + (
                (u_3 * (1 - p) * y) / (u_2 + (1 - p) * y))
        C_2 = -j_11 - j_22 - j_33
        C_1 = j_11 * j_22 + j_11 * j_33 + j_22 * j_33 - j_12 * j_21 - j_23 * j_32
        C_0 = (j_11 * j_23 * j_32 + j_12 * j_21 * j_33) - (
                j_11 * j_22 * j_33 + j_13 * j_21 * j_32)
        cond1 = C_2 > 0
        cond2 = C_1 > 0
        cond3 = C_0 > 0
        cond4 = C_2 * C_1 > C_0
        if not all([cond1, cond2, cond3, cond4]):
            del equilibria[i]
            pass
        pass
    return equilibria
