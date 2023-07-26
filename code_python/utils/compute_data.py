"""
Author: Rayla Kurosaki

GitHub: https://github.com/rkp1503
"""

import numpy as np
from scipy.integrate import odeint

SOL_MAX = 10


def solve_model(model, variables_dict, parameters_dict, t_max):
    ts = np.linspace(0, t_max, 10 * t_max + 1)
    variables_vals = list(variables_dict.values())
    parameter_vals = tuple(parameters_dict.values())
    return ts, odeint(model, variables_vals, ts, args=parameter_vals)


def get_bifurcation_data(model, variables_dict, parameters, t_max,
                         bifurcation_parameter, parameter_bounds,
                         parameter_values, data_points=10_000, debug=False):
    parameters_dict = dict(zip(parameters, parameter_values))
    x_t, x_s, x_u = [], [], []
    y_t, y_s, y_u = [], [], []
    z_t, z_s, z_u = [], [], []
    prev_x, prev_y, prev_z = "", "", ""
    curr_x, curr_y, curr_z = "", "", ""
    a, b = parameter_bounds[0], parameter_bounds[1]
    bifurcation_values = np.linspace(a, b, data_points + 1)
    for bifurcation_value in bifurcation_values:
        parameters_dict[bifurcation_parameter] = bifurcation_value
        ts, sol = solve_model(model, variables_dict, parameters_dict, t_max)
        n = int(len(ts) - np.floor(len(ts) * 0.95))
        x_min, x_max = min(sol[:, 0][-n:]), max(sol[:, 0][-n:])
        y_min, y_max = min(sol[:, 1][-n:]), max(sol[:, 1][-n:])
        z_min, z_max = min(sol[:, 2][-n:]), max(sol[:, 2][-n:])
        if val_in_bounds(x_min, x_max, data_points):
            x_t.append([bifurcation_value, x_min])
            x_t.append([bifurcation_value, x_max])
            if abs(x_min - x_max) <= 1 / data_points:
                curr_x = "Stable"
                x_s.append([bifurcation_value, x_min])
                x_s.append([bifurcation_value, x_max])
                pass
            else:
                curr_x = "Unstable"
                x_u.append([bifurcation_value, x_min])
                x_u.append([bifurcation_value, x_max])
                pass
            pass
        if val_in_bounds(y_min, y_max, data_points):
            y_t.append([bifurcation_value, y_min])
            y_t.append([bifurcation_value, y_max])
            if abs(y_min - y_max) <= 1 / data_points:
                curr_y = "Stable"
                y_s.append([bifurcation_value, y_min])
                y_s.append([bifurcation_value, y_max])
                pass
            else:
                curr_y = "Unstable"
                y_u.append([bifurcation_value, y_min])
                y_u.append([bifurcation_value, y_max])
                pass
            pass
        if val_in_bounds(z_min, z_max, data_points):
            z_t.append([bifurcation_value, z_min])
            z_t.append([bifurcation_value, z_max])
            if abs(z_min - z_max) <= 1 / data_points:
                curr_z = "Stable"
                z_s.append([bifurcation_value, z_min])
                z_s.append([bifurcation_value, z_max])
                pass
            else:
                curr_z = "Unstable"
                z_u.append([bifurcation_value, z_min])
                z_u.append([bifurcation_value, z_max])
                pass
            pass
        if (curr_x == curr_y == curr_z) and (curr_x != prev_x):
            if debug:
                if prev_x == "":
                    print(f"Start -> {curr_x}: {bifurcation_value}")
                    pass
                else:
                    print(f"{prev_x} -> {curr_x}: {bifurcation_value}")
                    pass
                pass
            prev_x = curr_x
            prev_y = curr_y
            prev_z = curr_z
            pass
        pass
    x_t = reformat_data(x_t)
    x_s = reformat_data(x_s)
    x_u = reformat_data(x_u)
    y_t = reformat_data(y_t)
    y_s = reformat_data(y_s)
    y_u = reformat_data(y_u)
    z_t = reformat_data(z_t)
    z_s = reformat_data(z_s)
    z_u = reformat_data(z_u)
    return x_t, x_s, x_u, y_t, y_s, y_u, z_t, z_s, z_u


def val_in_bounds(val_min, val_max, data_points):
    cond1 = 1 / data_points <= val_min < SOL_MAX
    cond2 = 1 / data_points <= val_max < SOL_MAX
    return cond1 and cond2


def reformat_data(data):
    if not data:
        return []
    return np.transpose(np.array(data))


def approximate_bounds():
    return None


def get_common_bounds(data):
    x_s, x_u = data[1][:, 0], data[2][:, 0]
    y_s, y_u = data[4][:, 0], data[5][:, 0]
    z_s, z_u = data[7][:, 0], data[8][:, 0]
    lst_s = list(x_s.intersection(y_s).intersection(z_s))
    lst_u = list(x_u.intersection(y_u).intersection(z_u))
    lst_t = list(set(lst_s) | set(lst_u))
    print(f"Full Bounds: [{min(lst_t)}, {max(lst_t)}]")
    print(f"Stable Bounds: [{min(lst_s)}, {max(lst_s)}]")
    print(f"Untable Bounds: [{min(lst_u)}, {max(lst_u)}]")
    return None
