"""
Author: Rayla Kurosaki

GitHub: https://github.com/rkp1503
"""

import sympy as sym

from utils import compute_data, generate_figures

x, X, y, Y, z, Z, t = sym.symbols("x, X, y, Y, z, Z, t")

r_yx, r_zx, p = sym.symbols("r_yx, r_zx, p")
phi_xy, phi_yx, phi_xz = sym.symbols("phi_xy, phi_yx, phi_xz")
u_1, u_2, u_3, u_4 = sym.symbols("u_1, u_2, u_3, u_4")

filename = ""


def main(variables_dict, parameters_dict, parameters_dict_2,
         t_max, t_max_b, model, labels_colors_dict):
    params_symb = parameters_dict.keys()
    ###########################################################################
    # The z-axial equilibrium
    ###########################################################################
    param_values = [0.007, 1.136, 0.874, 0.318, 0.416, 1.59, 1.655, 0.791,
                    0.994,
                    0.356]
    params_dict_temp = dict(zip(params_symb, param_values))
    ts, sol = compute_data.solve_model(model, variables_dict, params_dict_temp,
                                       t_max)
    title = "Time Evolution of Each Species"
    xaxis = "Time in Days"
    yaxis = "Population Density"
    generate_figures.my_plot(ts, sol, title, xaxis, yaxis, labels_colors_dict)
    ###########################################################################
    # The xy-boundary equilibrium
    ###########################################################################
    param_values = [0.049, 0.467, 0.645, 0.024, 0.163, 0.031, 0.31, 0.978, 0.9,
                    1.004]
    params_dict_temp = dict(zip(params_symb, param_values))
    ts, sol = compute_data.solve_model(model, variables_dict, params_dict_temp,
                                       t_max)
    generate_figures.my_plot(ts, sol, title, xaxis, yaxis, labels_colors_dict)
    ###########################################################################
    # The xz-boundary equilibrium
    ###########################################################################
    param_values = [0.199, 1.494, 0.482, 0.449, 0.993, 1.152, 1.671, 0.663,
                    1.556,
                    1.04]
    params_dict_temp = dict(zip(params_symb, param_values))
    ts, sol = compute_data.solve_model(model, variables_dict, params_dict_temp,
                                       t_max)
    generate_figures.my_plot(ts, sol, title, xaxis, yaxis, labels_colors_dict)
    ###########################################################################
    # The yz-boundary equilibrium
    ###########################################################################
    param_values = [1.219, 0.452, 0.589, 0.047, 1.587, 1.908, 1.658, 1.812,
                    1.473,
                    0.289]
    params_dict_temp = dict(zip(params_symb, param_values))
    ts, sol = compute_data.solve_model(model, variables_dict, params_dict_temp,
                                       t_max)
    generate_figures.my_plot(ts, sol, title, xaxis, yaxis, labels_colors_dict)
    ###########################################################################
    # The interior equilibrium
    ###########################################################################
    ts, sol = compute_data.solve_model(model, vars_dict, params_dict, t_max)
    generate_figures.my_plot(ts, sol, title, xaxis, yaxis, labels_colors_dict)
    ###########################################################################
    # 3D Phase Portrait for the interior equilibrium
    ###########################################################################
    sol_x, sol_y, sol_z = [sol[:, i] for i in range(len(vars_dict))]
    title = "3D Phase Portrait"
    xaxis = "Species $X$"
    yaxis = "Species $Y$"
    zaxis = "Species $Z$"
    generate_figures.my_3D_phase_portrait(sol, title, xaxis, yaxis, zaxis)
    ###########################################################################
    # xy-Phase Portrait for the interior equilibrium
    ###########################################################################
    title = f"2D Phase Portrait: {xaxis} and {yaxis}"
    generate_figures.my_2D_phase_portrait(sol_x, sol_y, title, xaxis, yaxis)
    ###########################################################################
    # xz-Phase Portrait for the interior equilibrium
    ###########################################################################
    title = f"2D Phase Portrait: {xaxis} and {zaxis}"
    generate_figures.my_2D_phase_portrait(sol_x, sol_z, title, xaxis, zaxis)
    ###########################################################################
    # yz-Phase Portrait for the interior equilibrium
    ###########################################################################
    title = f"2D Phase Portrait: {yaxis} and {zaxis}"
    generate_figures.my_2D_phase_portrait(sol_y, sol_z, title, yaxis, zaxis)
    ###########################################################################
    # Bifurcation Diagrams for r_yx
    ###########################################################################
    bifurcation_parameter = r_yx
    bifurcation_parameter_str = "r_{yx}"
    bifurcation_parameter_bounds = [0, 1]
    bifurcation_parameter_values = parameters_dict_2.values()
    bifurcation_parameter_values[
        list(params_symb).index(bifurcation_parameter)] = 0.5
    data = compute_data.get_bifurcation_data(model, vars_dict, params_symb,
                                             t_max_b, bifurcation_parameter,
                                             bifurcation_parameter_bounds,
                                             bifurcation_parameter_values)
    params_dict_temp = dict(zip(params_symb, bifurcation_parameter_values))
    ts, sol = compute_data.solve_model(model, vars_dict, params_dict_temp,
                                       t_max)
    title = "Time Evolution of Each Species"
    xaxis = "Time in Days"
    yaxis = "Population Density"
    generate_figures.my_plot(ts, sol, title, xaxis, yaxis, labels_colors_dict)
    generate_figures.my_bifurcation_diagram(data[1], data[2], X,
                                            bifurcation_parameter_str, "black")
    generate_figures.my_bifurcation_diagram(data[4], data[5], Y,
                                            bifurcation_parameter_str, "red")
    generate_figures.my_bifurcation_diagram(data[7], data[8], Z,
                                            bifurcation_parameter_str, "blue")
    ###########################################################################
    # Bifurcation Diagrams for r_zx
    ###########################################################################
    bifurcation_parameter = r_zx
    bifurcation_parameter_str = "r_{zx}"
    bifurcation_parameter_bounds = [0.133, 0.6155]
    bifurcation_parameter_values = parameters_dict.values()
    bifurcation_parameter_values[
        list(params_symb).index(bifurcation_parameter)] = 0.35
    data = compute_data.get_bifurcation_data(model, vars_dict, params_symb,
                                             t_max_b, bifurcation_parameter,
                                             bifurcation_parameter_bounds,
                                             bifurcation_parameter_values)
    params_dict_temp = dict(zip(params_symb, bifurcation_parameter_values))
    ts, sol = compute_data.solve_model(model, vars_dict, params_dict_temp,
                                       t_max)
    generate_figures.my_plot(ts, sol, title, xaxis, yaxis, labels_colors_dict)
    generate_figures.my_bifurcation_diagram(data[1], data[2], X,
                                            bifurcation_parameter_str, "black")
    generate_figures.my_bifurcation_diagram(data[4], data[5], Y,
                                            bifurcation_parameter_str, "red")
    generate_figures.my_bifurcation_diagram(data[7], data[8], Z,
                                            bifurcation_parameter_str, "blue")
    ###########################################################################
    # Bifurcation Diagrams for p
    ###########################################################################
    bifurcation_parameter = p
    bifurcation_parameter_str = "p"
    bifurcation_parameter_bounds = [0, 0.9476]
    bifurcation_parameter_values = parameters_dict.values()
    bifurcation_parameter_values[
        list(params_symb).index(bifurcation_parameter)] = 0.1
    data = compute_data.get_bifurcation_data(model, vars_dict, params_symb,
                                             t_max_b, bifurcation_parameter,
                                             bifurcation_parameter_bounds,
                                             bifurcation_parameter_values)
    params_dict_temp = dict(zip(params_symb, bifurcation_parameter_values))
    ts, sol = compute_data.solve_model(model, vars_dict, params_dict_temp,
                                       t_max)
    generate_figures.my_plot(ts, sol, title, xaxis, yaxis, labels_colors_dict)
    generate_figures.my_bifurcation_diagram(data[1], data[2], X,
                                            bifurcation_parameter_str, "black")
    generate_figures.my_bifurcation_diagram(data[4], data[5], Y,
                                            bifurcation_parameter_str, "red")
    generate_figures.my_bifurcation_diagram(data[7], data[8], Z,
                                            bifurcation_parameter_str, "blue")
    ###########################################################################
    # Bifurcation Diagrams for phi_xy
    ###########################################################################
    bifurcation_parameter = phi_xy
    bifurcation_parameter_str = "\\varphi_{xy}"
    bifurcation_parameter_bounds = [0, 0.9476]
    bifurcation_parameter_values = parameters_dict_2.values()
    bifurcation_parameter_values[
        list(params_symb).index(bifurcation_parameter)] = 0.15
    data = compute_data.get_bifurcation_data(model, vars_dict, params_symb,
                                             t_max_b, bifurcation_parameter,
                                             bifurcation_parameter_bounds,
                                             bifurcation_parameter_values)
    params_dict_temp = dict(zip(params_symb, bifurcation_parameter_values))
    ts, sol = compute_data.solve_model(model, vars_dict, params_dict_temp,
                                       t_max)
    generate_figures.my_plot(ts, sol, title, xaxis, yaxis, labels_colors_dict)
    generate_figures.my_bifurcation_diagram(data[1], data[2], X,
                                            bifurcation_parameter_str, "black")
    generate_figures.my_bifurcation_diagram(data[4], data[5], Y,
                                            bifurcation_parameter_str, "red")
    generate_figures.my_bifurcation_diagram(data[7], data[8], Z,
                                            bifurcation_parameter_str, "blue")
    ###########################################################################
    # Bifurcation Diagrams for phi_yx
    ###########################################################################
    bifurcation_parameter = phi_yx
    bifurcation_parameter_str = "\\varphi_{yx}"
    bifurcation_parameter_bounds = [0.0, 0.444]
    bifurcation_parameter_values = parameters_dict.values()
    bifurcation_parameter_values[
        list(params_symb).index(bifurcation_parameter)] = 0.43
    data = compute_data.get_bifurcation_data(model, vars_dict, params_symb,
                                             t_max_b, bifurcation_parameter,
                                             bifurcation_parameter_bounds,
                                             bifurcation_parameter_values)
    params_dict_temp = dict(zip(params_symb, bifurcation_parameter_values))
    ts, sol = compute_data.solve_model(model, vars_dict, params_dict_temp,
                                       t_max)
    generate_figures.my_plot(ts, sol, title, xaxis, yaxis, labels_colors_dict)
    generate_figures.my_bifurcation_diagram(data[1], data[2], X,
                                            bifurcation_parameter_str, "black")
    generate_figures.my_bifurcation_diagram(data[4], data[5], Y,
                                            bifurcation_parameter_str, "red")
    generate_figures.my_bifurcation_diagram(data[7], data[8], Z,
                                            bifurcation_parameter_str, "blue")
    ###########################################################################
    # Bifurcation Diagrams for phi_xz
    ###########################################################################
    bifurcation_parameter = phi_xz
    bifurcation_parameter_str = "\\varphi_{xz}"
    bifurcation_parameter_bounds = [0, 2.252]
    bifurcation_parameter_values = parameters_dict_2.values()
    bifurcation_parameter_values[
        list(params_symb).index(bifurcation_parameter)] = 0.5
    data = compute_data.get_bifurcation_data(model, vars_dict, params_symb,
                                             t_max_b, bifurcation_parameter,
                                             bifurcation_parameter_bounds,
                                             bifurcation_parameter_values)
    params_dict_temp = dict(zip(params_symb, bifurcation_parameter_values))
    ts, sol = compute_data.solve_model(model, vars_dict, params_dict_temp,
                                       t_max)
    generate_figures.my_plot(ts, sol, title, xaxis, yaxis, labels_colors_dict)
    generate_figures.my_bifurcation_diagram(data[1], data[2], X,
                                            bifurcation_parameter_str, "black")
    generate_figures.my_bifurcation_diagram(data[4], data[5], Y,
                                            bifurcation_parameter_str, "red")
    generate_figures.my_bifurcation_diagram(data[7], data[8], Z,
                                            bifurcation_parameter_str, "blue")
    ###########################################################################
    # Bifurcation Diagrams for u_1
    ###########################################################################
    bifurcation_parameter = u_1
    bifurcation_parameter_str = "u_1"
    bifurcation_parameter_bounds = [0.0, 1.0]
    bifurcation_parameter_values = parameters_dict_2.values()
    bifurcation_parameter_values[
        list(params_symb).index(bifurcation_parameter)] = 0.8
    data = compute_data.get_bifurcation_data(model, vars_dict, params_symb,
                                             t_max_b, bifurcation_parameter,
                                             bifurcation_parameter_bounds,
                                             bifurcation_parameter_values)
    params_dict_temp = dict(zip(params_symb, bifurcation_parameter_values))
    ts, sol = compute_data.solve_model(model, vars_dict, params_dict_temp,
                                       t_max)
    generate_figures.my_plot(ts, sol, title, xaxis, yaxis, labels_colors_dict)
    generate_figures.my_bifurcation_diagram(data[1], data[2], X,
                                            bifurcation_parameter_str, "black")
    generate_figures.my_bifurcation_diagram(data[4], data[5], Y,
                                            bifurcation_parameter_str, "red")
    generate_figures.my_bifurcation_diagram(data[7], data[8], Z,
                                            bifurcation_parameter_str, "blue")
    ###########################################################################
    # Bifurcation Diagrams for u_2
    ###########################################################################
    bifurcation_parameter = u_2
    bifurcation_parameter_str = "u_2"
    bifurcation_parameter_bounds = [0.0354, 0.1]
    bifurcation_parameter_values = parameters_dict.values()
    bifurcation_parameter_values[
        list(params_symb).index(bifurcation_parameter)] = 0.04
    data = compute_data.get_bifurcation_data(model, vars_dict, params_symb,
                                             t_max_b, bifurcation_parameter,
                                             bifurcation_parameter_bounds,
                                             bifurcation_parameter_values)
    params_dict_temp = dict(zip(params_symb, bifurcation_parameter_values))
    ts, sol = compute_data.solve_model(model, vars_dict, params_dict_temp,
                                       t_max)
    generate_figures.my_plot(ts, sol, title, xaxis, yaxis, labels_colors_dict)
    generate_figures.my_bifurcation_diagram(data[1], data[2], X,
                                            bifurcation_parameter_str, "black")
    generate_figures.my_bifurcation_diagram(data[4], data[5], Y,
                                            bifurcation_parameter_str, "red")
    generate_figures.my_bifurcation_diagram(data[7], data[8], Z,
                                            bifurcation_parameter_str, "blue")
    ###########################################################################
    # Bifurcation Diagrams for u_3
    ###########################################################################
    bifurcation_parameter = u_3
    bifurcation_parameter_str = "u_3"
    bifurcation_parameter_bounds = [0.208, 3]
    bifurcation_parameter_values = parameters_dict_2.values()
    bifurcation_parameter_values[
        list(params_symb).index(bifurcation_parameter)] = 0.75
    data = compute_data.get_bifurcation_data(model, vars_dict, params_symb,
                                             t_max_b, bifurcation_parameter,
                                             bifurcation_parameter_bounds,
                                             bifurcation_parameter_values)
    params_dict_temp = dict(zip(params_symb, bifurcation_parameter_values))
    ts, sol = compute_data.solve_model(model, vars_dict, params_dict_temp,
                                       t_max)
    generate_figures.my_plot(ts, sol, title, xaxis, yaxis, labels_colors_dict)
    generate_figures.my_bifurcation_diagram(data[1], data[2], X,
                                            bifurcation_parameter_str, "black")
    generate_figures.my_bifurcation_diagram(data[4], data[5], Y,
                                            bifurcation_parameter_str, "red")
    generate_figures.my_bifurcation_diagram(data[7], data[8], Z,
                                            bifurcation_parameter_str, "blue")
    ###########################################################################
    # Bifurcation Diagrams for u_4
    ###########################################################################
    bifurcation_parameter = u_4
    bifurcation_parameter_str = "u_4"
    bifurcation_parameter_bounds = [0.078, 0.55]
    bifurcation_parameter_values = parameters_dict_2.values()
    bifurcation_parameter_values[
        list(params_symb).index(bifurcation_parameter)] = 0.3
    data = compute_data.get_bifurcation_data(model, vars_dict, params_symb,
                                             t_max_b, bifurcation_parameter,
                                             bifurcation_parameter_bounds,
                                             bifurcation_parameter_values)
    params_dict_temp = dict(zip(params_symb, bifurcation_parameter_values))
    ts, sol = compute_data.solve_model(model, vars_dict, params_dict_temp,
                                       t_max)
    generate_figures.my_plot(ts, sol, title, xaxis, yaxis, labels_colors_dict)
    generate_figures.my_bifurcation_diagram(data[1], data[2], X,
                                            bifurcation_parameter_str, "black")
    generate_figures.my_bifurcation_diagram(data[4], data[5], Y,
                                            bifurcation_parameter_str, "red")
    generate_figures.my_bifurcation_diagram(data[7], data[8], Z,
                                            bifurcation_parameter_str, "blue")
    return None


if __name__ == '__main__':
    # Initializing variable values
    vars_dict = {
        # Density of species X
        x: 0.7,
        # Density of species Y
        y: 0.4,
        # Density of species Z
        z: 0.2,
    }
    # Initializing graph labels and colors
    labels_colors_dict = {
        x: ["X", "Black"],
        y: ["Y", "Red"],
        z: ["Z", "Blue"],
    }
    # Initializing parameter values
    params_dict = {
        # Ratio of intrinsic growth rate of species Y to species X
        r_yx: 0.5,
        # Ratio of intrinsic growth rate of species Z to species X
        r_zx: 0.5,
        # Refuge rate of species Y
        p: 0.6,
        # Scaled interspecies mutualism coefficient of species Y on species X
        phi_xy: 0.6,
        # Scaled interspecies competition coefficient of species X on species Y
        phi_yx: 0.15,
        # Scaled commensal coefficient of species Z on species X
        phi_xz: 0.4,
        # Scaled attack rate of species Z on species Y
        u_1: 0.6,
        # Scaled half saturation constant for Holling type II function
        u_2: 0.08,
        # Scaled conservation rate of species Y
        u_3: 0.5,
        # Scaled death rate of species Z
        u_4: 0.5,
    }
    # Creating a dictionary for a second set of parameter values
    params_vals_2 = [0.7, 0.15, 0.4, 0.05, 0.5, 0.04, 0.7, 0.2, 0.5,
                     0.32]
    params_dict_2 = dict(zip(params_dict.keys(), params_vals_2))
    # Number of days to model.
    t_max = 1_000
    t_max_b = 10_000


    def model(variables, t, *parameters):
        x, y, z = variables
        r_yx, r_zx, p, phi_xy, phi_yx, phi_xz, u_1, u_2, u_3, u_4 = parameters[
            0] if len(parameters) == 1 else parameters
        dx = x * (1 - x + phi_xy * y ** 2) - phi_xz * x * z
        dy = r_yx * y * (1 - y + phi_yx * x ** 2) - (
                (u_1 * (1 - p) * y * z) / (u_2 + (1 - p) * y))
        dz = r_zx * z * (1 - z) + z * (
                ((u_3 * (1 - p) * y) / (u_2 + (1 - p) * y)) - u_4)
        return [dx, dy, dz]


    main(vars_dict, params_dict, params_dict_2, t_max, t_max_b, model,
         labels_colors_dict)
    pass
