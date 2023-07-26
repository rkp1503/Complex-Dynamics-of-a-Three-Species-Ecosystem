"""
Author: Rayla Kurosaki

GitHub: https://github.com/rkp1503
"""

import numpy as np

from utils import compute_data
from equilibria import axial_z
from equilibria import boundary_xy
from equilibria import boundary_xz
from equilibria import boundary_yz
from equilibria import interior

MAX_ITER = 1_000_000
ROUND = 3
TOL = 1e-3


def generate_parameters(model, variables_dict, parameters_dict, t_max,
                        equilibrium_type, debug=False):
    for i in range(1, MAX_ITER + 1):
        parameter_values = generate_parameters_helper(parameters_dict,
                                                      equilibrium_type, debug)
        if not parameter_values:
            parameters_dict_temp = dict(
                zip(parameters_dict.keys(), parameter_values))
            ts, sol = compute_data.solve_model(model, variables_dict,
                                               parameters_dict_temp, t_max)
            sol_x, sol_y, sol_z = sol[-1]
            if compare_results(sol, parameter_values, equilibrium_type):
                print(f"Sutible parameters found at iteration {i}:")
                print(f"Parameter Values: {parameter_values}")
                print("Equilibrium:")
                print(f"\tX: {np.round(sol_x, decimals=ROUND)}")
                print(f"\tY: {np.round(sol_y, decimals=ROUND)}")
                print(f"\tZ: {np.round(sol_z, decimals=ROUND)}")
                return parameter_values
            else:
                if debug:
                    print(
                        "True and approximated solutions are not close enough!\n")
                    pass
                pass
            pass
        pass
    print(f"Failed to find sutible parameters after {i} iterations.")
    return None


def generate_parameters_helper(parameters_dict, equilibrium_type, debug=False):
    while True:
        parameter_values = np.random.rand(len(parameters_dict))
        parameter_values = np.round(parameter_values, decimals=ROUND)
        if (0 not in parameter_values) and (parameter_values[2] < 1):
            equilibria = get_equilibria(parameter_values, equilibrium_type)
            if debug:
                print(f"Number of Equilibria: {len(equilibria)}")
                pass
            if not equilibria:
                equilibria = get_stable_equilibria(equilibria,
                                                   parameter_values,
                                                   equilibrium_type)
                if debug:
                    print(f"Number of Stable Equilibria: {len(equilibria)}")
                    pass
                if not equilibria:
                    return parameter_values
                    pass
                else:
                    if debug:
                        print(
                            f"The {equilibrium_type} equilibrium is not stable for the parameters {parameter_values}\n")
                        pass
                    pass
                pass
            else:
                if debug:
                    print(
                        f"The {equilibrium_type} equilibrium does not exist for the parameters {parameter_values}\n")
                    pass
                pass
            return []
        pass
    pass


def get_equilibria(parameter_values, equilibrium_type):
    if equilibrium_type == "z-axial":
        return axial_z.get_equilibria(parameter_values)
    elif equilibrium_type == "xy-boundary":
        return boundary_xy.get_equilibria(parameter_values)
    elif equilibrium_type == "xz-boundary":
        return boundary_xz.get_equilibria(parameter_values)
    elif equilibrium_type == "yz-boundary":
        return boundary_yz.get_equilibria(parameter_values)
    else:
        return interior.get_equilibria(parameter_values)
    pass


def get_stable_equilibria(equilibria_lst, parameter_values, equilibrium_type):
    if equilibrium_type == "z-axial":
        return axial_z.get_stable_equilibria(equilibria_lst, parameter_values)
    elif equilibrium_type == "xy-boundary":
        return boundary_xy.get_stable_equilibria(equilibria_lst,
                                                 parameter_values)
    elif equilibrium_type == "xz-boundary":
        return boundary_xz.get_stable_equilibria(equilibria_lst,
                                                 parameter_values)
    elif equilibrium_type == "yz-boundary":
        return boundary_yz.get_stable_equilibria(equilibria_lst,
                                                 parameter_values)
    else:
        return interior.get_stable_equilibria(equilibria_lst, parameter_values)
    pass


def compare_results(solutions, parameter_values, equilibrium_type):
    sol_x, sol_y, sol_z = solutions[-1]
    equilibria = get_equilibria(parameter_values, equilibrium_type)
    equilibria = get_stable_equilibria(equilibria, parameter_values,
                                       equilibrium_type)
    E_x, E_y, E_z = equilibria[1]
    x_err = np.round(abs(E_x - sol_x), decimals=ROUND)
    y_err = np.round(abs(E_y - sol_y), decimals=ROUND)
    z_err = np.round(abs(E_z - sol_z), decimals=ROUND)
    return (x_err <= TOL) and (y_err <= TOL) and (z_err <= TOL)


def validate_parameters(parameter_values, equilibrium_type, debug=False):
    equilibria = get_equilibria(parameter_values, equilibrium_type)
    if debug:
        print(f"Number of Equilibria: {len(equilibria)}")
        pass
    if not equilibria:
        print(
            f"The {equilibrium_type} equilibrium does not exist for the parameters {parameter_values}")
        pass
    equilibria = get_stable_equilibria(equilibria, parameter_values,
                                       equilibrium_type)
    if debug:
        print(f"Number of Stable Equilibria: {len(equilibria)}")
        pass
    if not equilibria:
        print(
            f"The {equilibrium_type} equilibrium is not stable for the parameters {parameter_values}")
        pass
    print(
        f"The {equilibrium_type} equilibrium exists and is stable for the parameters {parameter_values}")
    return None
