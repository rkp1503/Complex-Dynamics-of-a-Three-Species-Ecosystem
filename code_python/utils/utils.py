"""
Author: Rayla Kurosaki

GitHub: https://github.com/rkp1503
"""


import os

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
import sympy as sym


def compute_equilibria(model, vars_dict, params_dict, equilibria, equilibrium_vals_temp):
    vars_symbols = list(vars_dict.keys())
    params_symbols = list(params_dict.keys())
    variables_dict = dict(zip(vars_symbols, equilibrium_vals_temp))
    system_model = model(vars_symbols, None, tuple(params_symbols))
    for i, eq in enumerate(system_model):
        system_model[i] = eq.subs(variables_dict)
        pass
    if system_model == [0, 0, 0]:
        equilibrium = dict(zip(vars_symbols, list((0, 0, 0))))
        if equilibrium not in equilibria:
            return equilibrium
        pass
    else:
        solutions = sym.solve(system_model, vars_symbols)
        equilibrium = dict(zip(vars_symbols, equilibrium_vals_temp))
        if type(solutions) == dict:
            for (k, v) in solutions.items():
                equilibrium[k] = v
                pass
            if equilibrium not in equilibria:
                return equilibrium
            pass
        else:
            for sol in solutions:
                equilibrium = dict(zip(vars_symbols, sol))
                for i, (k, v) in enumerate(equilibrium.items()):
                    if k == v:
                        equilibrium[k] = equilibrium_vals_temp[i]
                        pass
                    pass
                if equilibrium not in equilibria:
                    return equilibrium
                pass
            pass
        pass
    return None


def construct_jacobian(model, vars_dict, params_dict, equilibrium):
    vars_symbols = list(vars_dict.keys())
    params_symbols = list(params_dict.keys())
    system_model = sym.Matrix(model(vars_symbols, None, params_symbols))
    model_j = system_model.jacobian(vars_symbols)
    J = model_j.subs(equilibrium)
    (n, m) = sym.shape(J)
    J_s = sym.zeros(n, m)
    for i in range(n):
        for j in range(m):
            J[i, j] = sym.simplify(J[i,j])
            if J[i, j] != 0:
                J_s[i, j] = sym.symbols(f"j{i + 1}{j + 1}")
                pass
            pass
        pass
    return J, J_s


def stability_analysis(jacobian, jacobian_s, compact=False, method="Linearization"):
    lmbd = sym.symbols("lambda")
    display(jacobian)
    if method == "RHC":
        if compact:
            display(jacobian_s)
            char_poly = jacobian_s.charpoly(lmbd).as_expr()
            pass
        else:
            char_poly = (jacobian.charpoly(lmbd)).as_expr()
            pass
        display(char_poly)
        C_3 = char_poly.coeff(lmbd**3)
        C_2 = char_poly.coeff(lmbd**2)
        C_1 = char_poly.coeff(lmbd**1)
        C_0 = char_poly-(C_3*lmbd**3+C_2*lmbd**2+C_1*lmbd**1)
        print("C_3:")
        display(sym.simplify(C_3))
        print("C_2:")
        display(sym.simplify(C_2))
        print("C_1:")
        display(sym.simplify(C_1))
        print("C_0:")
        display(sym.simplify(C_0))
        print("C_2C_1-C_0:")
        display(sym.simplify(C_2*C_1-C_0))
        pass
    else:
        if compact:
            display(jacobian_s)
            char_poly = sym.simplify((jacobian_s.charpoly(lmbd)).as_expr())
            pass
        else:
            char_poly = sym.simplify((jacobian.charpoly(lmbd)).as_expr())
            pass
        display(char_poly)
        eigenvalues = sym.solve(char_poly, lmbd)
        display(eigenvalues)
        pass
    return None


def solve_model(model, t_max, vars_dict, param_values):
    ts = np.linspace(0, t_max, 10*t_max+1)
    return ts, odeint(model, list(vars_dict.values()), ts, args=tuple(param_values))


def save_figure(filename: str):
    path_to_project_dir: str = os.path.dirname(os.getcwd())
    path_to_report_dir: str = os.path.join(path_to_project_dir, "report")
    path_to_figs_dir: str = os.path.join(path_to_report_dir, "figures")
    if os.path.exists(path_to_figs_dir):
        path_to_png_dir: str = os.path.join(path_to_figs_dir, "png")
        if not os.path.exists(path_to_png_dir):
            os.makedirs(path_to_png_dir)
            pass
        path_to_png: str = os.path.join(path_to_png_dir, filename)
        plt.savefig(f"{path_to_png}.png", bbox_inches="tight")
        path_to_eps_dir: str = os.path.join(path_to_figs_dir, "eps")
        if not os.path.exists(path_to_eps_dir):
            os.makedirs(path_to_eps_dir)
            pass
        path_to_eps: str = os.path.join(path_to_eps_dir, filename)
        plt.savefig(f"{path_to_eps}.eps", bbox_inches="tight")
        pass
    return None


def my_plot(ts, sol, title, xaxis, yaxis, labels_colors_dict, filename=None):
    graphs = []
    for (i, label_color) in enumerate(labels_colors_dict.values()):
        label, color = label_color
        temp_graph = plt.plot(ts, sol[:, i], label=label, color=color, linewidth=1)
        graphs += temp_graph
        pass
    labels = [graph.get_label() for graph in graphs]
    plt.legend(graphs, labels)
    # plt.legend(graphs, labels, bbox_to_anchor=(1.01, 1.015), loc='upper left')
    plt.title(title)
    plt.xlabel(xaxis)
    plt.ylabel(yaxis)
    # plt.minorticks_on()
    # plt.grid(which='major', linestyle='-', linewidth='0.5', color='#000000')
    # plt.grid(which='minor', linestyle=':', linewidth='0.5', color='#FF0000')
    if filename is not None:
        save_figure(filename)
        pass
    plt.show()
    return None


def my_3D_phase_portrait(sol, title, xaxis, yaxis, zaxis, filename=None):
    ax = plt.axes(projection='3d')
    ax.plot3D(sol[:, 0], sol[:, 1], sol[:, 2])
    plt.title(title)
    ax.set_xlabel(xaxis)
    ax.set_ylabel(yaxis)
    ax.set_zlabel(zaxis)
    # plt.minorticks_on()
    # plt.grid(which='major', linestyle='-', linewidth='0.5', color='#000000')
    # plt.grid(which='minor', linestyle=':', linewidth='0.5', color='#FF0000')
    if filename is not None:
        save_figure(filename)
        pass
    plt.show()
    return None


def my_2D_phase_portrait(sol_x, sol_y, title, xaxis, yaxis, filename=None):
    plt.plot(sol_x, sol_y, linewidth=1)
    plt.title(title)
    plt.xlabel(xaxis)
    plt.ylabel(yaxis)
    # plt.minorticks_on()
    # plt.grid(which='major', linestyle='-', linewidth='0.5', color='#000000')
    # plt.grid(which='minor', linestyle=':', linewidth='0.5', color='#FF0000')
    if filename is not None:
        save_figure(filename)
        pass
    plt.show()
    return None
