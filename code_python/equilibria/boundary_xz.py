"""
Author: Rayla Kurosaki

GitHub: https://github.com/rkp1503
"""


def get_equilibria(parameter_values):
    equilibria = get_analytical_solutions(parameter_values)
    r_yx, r_zx, p, phi_xy, phi_yx, phi_xz, u_1, u_2, u_3, u_4 = parameter_values
    for (i, equilibrium) in enumerate(equilibria):
        x, y, z = equilibrium
        # if not ((x > 0) and (z > 0)):
        #     del equilibria[i]
        #     pass
        cond1 = (u_4 / r_zx) + (1 / phi_xz) > 1
        cond2 = r_zx > u_4
        if not all([cond1, cond2]):
            del equilibria[i]
            pass
        pass
    return equilibria


def get_analytical_solutions(parameter_values):
    r_yx, r_zx, p, phi_xy, phi_yx, phi_xz, u_1, u_2, u_3, u_4 = parameter_values
    x = 1 - phi_xz * (1 - (u_4 / r_zx))
    z = 1 - (u_4 / r_zx)
    return [(x, 0, z)]


def get_stable_equilibria(equilibria, parameter_values):
    r_yx, r_zx, p, phi_xy, phi_yx, phi_xz, u_1, u_2, u_3, u_4 = parameter_values
    for (i, equilibrium) in enumerate(equilibria):
        x, y, z = equilibrium
        cond1 = 1 - (u_4 / r_zx) > (
                1 - 2 * (1 - phi_xz * (1 - (u_4 / r_zx)))) / phi_xz
        cond2 = 1 - (u_4 / r_zx) > r_yx * u_2 * (
                1 + phi_yx * (1 - phi_xz * (1 - (u_4 / r_zx))) ** 2) / (
                        u_1 * (1 - p))
        cond3 = 1 - (u_4 / r_zx) > (r_zx - u_4) / (2 * r_zx)
        if not all([cond1, cond2, cond3]):
            del equilibria[i]
            pass
        pass
    return equilibria
