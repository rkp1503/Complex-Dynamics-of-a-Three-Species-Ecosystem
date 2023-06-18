"""
Author: Ramsey (Rayla) Phuc

Alias: Rayla Kurosaki

GitHub: https://github.com/rkp1503

Co-author: Ephraim Agyingi
"""

import numpy as np
import sympy as sym

from utils import utils


def main(model: any, t_max: int, vars_dict: dict[sym.core.Symbol, float], params_dict: dict[sym.core.Symbol, float], labels_colors_dict: dict[sym.core.Symbol, list[str]]) -> None:
    # ------------------------------------------------------------------------
    # Plotting time evolution for the z-axial equilibrium
    # ------------------------------------------------------------------------
    param_values: list[float] = [0.404, 0.903, 0.182, 0.639, 0.283, 0.301, 0.110, 0.645, 0.175, 0.145]
    ts, sol = utils.solve_model(model, t_max, vars_dict, param_values)
    title: str = "Time Evolution of Each Species"
    xaxis: str = "Time in Days"
    yaxis: str = "Population size"
    utils.my_plot(ts, sol, title, xaxis, yaxis, labels_colors_dict)
    # ------------------------------------------------------------------------
    # Plotting time evolution for the xy-boundary equilibrium
    # ------------------------------------------------------------------------
    param_values: list[float] = [0.978, 0.613, 0.326, 0.245, 0.015, 0.920, 0.696, 0.523, 0.951, 0.570]
    ts, sol = utils.solve_model(model, t_max, vars_dict, param_values)
    utils.my_plot(ts, sol, title, xaxis, yaxis, labels_colors_dict)
    # ------------------------------------------------------------------------
    # Plotting time evolution for the xz-boundary equilibrium
    # ------------------------------------------------------------------------
    param_values: list[float] = [0.102, 0.763, 0.271, 0.182, 0.301, 0.109, 0.198, 0.983, 0.186, 0.113]
    ts, sol = utils.solve_model(model, t_max, vars_dict, param_values)
    utils.my_plot(ts, sol, title, xaxis, yaxis, labels_colors_dict)
    # ------------------------------------------------------------------------
    # Plotting time evolution for the yz-boundary equilibrium
    # ------------------------------------------------------------------------
    param_values: list[float] = [0.978, 0.310, 0.843, 0.002, 0.407, 0.859, 0.446, 0.872, 0.201, 0.959]
    ts, sol = utils.solve_model(model, t_max, vars_dict, param_values)
    utils.my_plot(ts, sol, title, xaxis, yaxis, labels_colors_dict)
    # ------------------------------------------------------------------------
    # Plotting time evolution for the interior equilibrium
    # ------------------------------------------------------------------------
    param_values: list[float] = list(params_dict.values())
    ts, sol = utils.solve_model(model, t_max, vars_dict, param_values)
    utils.my_plot(ts, sol, title, xaxis, yaxis, labels_colors_dict)
    # ------------------------------------------------------------------------
    # Plotting a 3D phase portrait for the interior equilibrium
    # ------------------------------------------------------------------------
    xaxis: str = "Species X"
    yaxis: str = "Species Y"
    zaxis: str = "Species Z"
    sol_x, sol_y, sol_z = [sol[:, i] for i in range(len(vars_dict))]
    title: str = "3D Phase Portrait"
    utils.my_3D_phase_portrait(sol, title, xaxis, yaxis, zaxis)
    # ------------------------------------------------------------------------
    # Plotting a 2D phase portrait of Species X and Y for the interior equilibrium
    # ------------------------------------------------------------------------
    title: str = f"2D Phase Portrait: {xaxis} and {yaxis}"
    utils.my_2D_phase_portrait(sol_x, sol_y, title, xaxis, yaxis)
    # ------------------------------------------------------------------------
    # Plotting a 2D phase portrait of Species X and Z for the interior equilibrium
    # ------------------------------------------------------------------------
    title: str = f"2D Phase Portrait: {xaxis} and {zaxis}"
    utils.my_2D_phase_portrait(sol_x, sol_z, title, xaxis, zaxis)
    # ------------------------------------------------------------------------
    # Plotting a 2D phase portrait of Species Y and Z for the interior equilibrium
    # ------------------------------------------------------------------------
    title: str = f"2D Phase Portrait: {yaxis} and {zaxis}"
    utils.my_2D_phase_portrait(sol_y, sol_z, title, yaxis, zaxis)
    return None


def model(variables: np.ndarray, t: float, *parameters: tuple) -> list[np.float64]:
    """
    This function constructs the system of equations which resembles the equations in the proposed model.

    :param variables: The dependent variables.
    :param t: The independent variable.
    :param parameters: The parameters.
    :return: The system of equations.
    """
    # Unpacking the variables and parameters
    x, y, z = variables
    r1, r2, p, gamma12, gamma21, gamma13, gamma31, v1, v2, v3 = parameters[0] if len(parameters) == 1 else parameters
    # Constructing each equation
    dx = x * (1 - x + gamma12 * y ** 2) - gamma13 * x * z
    dy = r1 * y * (1 - y + gamma21 * x ** 2) - (((1 - p) * y * z) / (v1 + (1 - p) * y))
    dz = r2 * z * (1 - gamma31 * z) + z * (((v3 * (1 - p) * y) / (v1 + (1 - p) * y)) - v2)
    # Return the model as a list of equations
    return [dx, dy, dz]


if __name__ == "__main__":
    # Initializing variables
    x, y, z, t = sym.symbols("x, y, z, t")
    vars_dict: dict[sym.core.Symbol, float] = {
        # Densities of two logistically growing competing species
        x: 0.7,
        y: 0.4,
        # Density of a species which is a predator of 𝑌
        z: 0.2,
    }
    # Initializing parameters
    r1, r2, p, gamma12, gamma21, gamma13, gamma31, v1, v2, v3 = sym.symbols("r1 r2 p gamma12 gamma21 gamma13 gamma31 v1 v2 v3")
    params_dict: dict[sym.core.Symbol, float] = {
        # Ratio of intrinsic growth rate of Y to X
        r1: 0.6353,
        # Ratio of intrinsic growth rate of Z to X
        r2: 0.742,
        # Refuge rate of Y
        p: 0.852,
        # Interspecies competition coefficient of Y on X
        gamma12: 0.142,
        # Interspecies competition coefficient of X on Y
        gamma21: 0.002,
        #
        gamma13: 0.148,
        #
        gamma31: 0.215,
        # Half saturation constant for Holling type II function
        v1: 0.090,
        # Natural death rate of Z
        v2: 0.891,
        # Conservation rate of Y
        v3: 0.980,
    }
    # Initializing graph labels and colors for each variable
    labels_colors_dict: dict[sym.core.Symbol, list[str]] = {
        x: ["X", "Black"],
        y: ["Y", "Red"],
        z: ["Z", "Blue"],
    }
    # Number of days to model
    t_max: int = 500
    # Displaying indicators
    display_equilibria: bool = False
    display_stability: bool = False
    # Call the function to plot graphs
    main(model, t_max, vars_dict, params_dict, labels_colors_dict)
    pass
