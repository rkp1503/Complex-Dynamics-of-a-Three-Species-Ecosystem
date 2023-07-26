"""
Author: Rayla Kurosaki

GitHub: https://github.com/rkp1503
"""

import os

import matplotlib.pyplot as plt


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


def my_plot(ts, solutions, title, xaxis, yaxis, legend_data, norm=1,
            detailed_fig=False, filename=""):
    graphs = []
    for (i, style) in enumerate(legend_data.values()):
        label, color = style
        temp_graph = plt.plot(ts, solutions[:, i] / norm, label=label,
                              color=color, linewidth=1)
        graphs += temp_graph
        pass
    labels = [graph.get_label() for graph in graphs]
    plt.title(title)
    plt.xlabel(xaxis)
    plt.ylabel(yaxis)
    if detailed_fig:
        plt.legend(graphs, labels, bbox_to_anchor=(1.01, 1.015),
                   loc='upper left')
        plt.minorticks_on()
        plt.grid(which='major', linestyle='-', linewidth='0.5',
                 color='#000000')
        plt.grid(which='minor', linestyle=':', linewidth='0.5',
                 color='#FF0000')
        pass
    else:
        plt.legend(graphs, labels)
        pass
    if filename != "":
        save_figure(filename)
        pass
    plt.show()
    return None


def my_3D_phase_portrait(sol, title, xaxis, yaxis, zaxis, detailed_fig=False,
                         filename=""):
    ax = plt.axes(projection='3d')
    ax.plot3D(sol[:, 0], sol[:, 1], sol[:, 2])
    plt.title(title)
    ax.set_xlabel(xaxis)
    ax.set_ylabel(yaxis)
    ax.set_zlabel(zaxis)
    if detailed_fig:
        plt.minorticks_on()
        plt.grid(which='major', linestyle='-', linewidth='0.5',
                 color='#000000')
        plt.grid(which='minor', linestyle=':', linewidth='0.5',
                 color='#FF0000')
        pass
    if filename != "":
        save_figure(filename)
        pass
    plt.show()
    return None


def my_2D_phase_portrait(solution_x, solution_y, title, xaxis, yaxis,
                         detailed_fig=False, filename=""):
    plt.plot(solution_x, solution_y, linewidth=1)
    plt.title(title)
    plt.xlabel(xaxis)
    plt.ylabel(yaxis)
    if detailed_fig:
        plt.minorticks_on()
        plt.grid(which='major', linestyle='-', linewidth='0.5',
                 color='#000000')
        plt.grid(which='minor', linestyle=':', linewidth='0.5',
                 color='#FF0000')
        pass
    if filename is not None:
        save_figure(filename)
        pass
    plt.show()
    return None


def my_bifurcation_diagram(data_stable, data_unstable, variable, parameter,
                           color, detailed_fig=False, filename=""):
    if data_stable == []:
        parameters_stable, solution_stable = [], []
        pass
    else:
        parameters_stable, solution_stable = data_stable[0], data_stable[1]
        pass
    if data_stable == []:
        parameters_unstable, solution_unstable = [], []
        pass
    else:
        parameters_unstable, solution_unstable = data_unstable[0], data_unstable[1]
        pass
    plt.title(f"Density of Species ${variable}$ vs ${parameter}$")
    plt.xlabel(f"${parameter}$")
    plt.ylabel(f"Density of Species ${variable}$")
    plt.scatter(parameters_stable, solution_stable, s=1, color=color,
                label="Stable")
    plt.scatter(parameters_unstable, solution_unstable, s=1, color="g",
                label="Unstable")
    plt.legend()
    if filename is not None:
        save_figure(filename)
        pass
    plt.show()
    return None
