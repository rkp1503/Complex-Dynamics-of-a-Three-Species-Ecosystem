#=============================================================================
Author: Rayla Kurosaki
GitHub: https://github.com/rkp1503
=============================================================================#

module GenerateFigures

import Plots: @layout, plot, plot!, scatter, scatter!

function my_plot(solutions, title, xaxis, yaxis, legend_data; norm=1, detailed_fig=false, export_fig=false)
    figure = plot(title=title, xaxis=xaxis, yaxis=yaxis, gridalpha=0)
    if detailed_fig
        plot!(legend=:outertopright, background_color_inside="#DDDDDD", gridalpha=1, gridlinewidth=:1, foreground_color_grid="#000000", minorgrid=true, minorgridalpha=0.5, foreground_color_minor_grid=:"#FF0000")
    end
    for (i, style) in enumerate(collect(values(legend_data)))
        label, color = style
        plot!(figure, solutions.t, solutions'[:, i]/norm, label=label, linecolor=color, linewidth=1.5)
    end
    if export_fig
    end
    return figure
end;

function my_3D_phase_portrait(solutions, title, xaxis, yaxis, zaxis; detailed_fig=false, export_fig=false)
    figure = plot(title=title, xaxis=xaxis, yaxis=yaxis, zaxis=zaxis, legend=false, gridalpha=0)
    if detailed_fig
        plot!(background_color_inside="#DDDDDD", gridalpha=1, gridlinewidth=:1, foreground_color_grid="#000000", minorgrid=true, minorgridalpha=0.5, foreground_color_minor_grid=:"#FF0000")
    end
    plot!(figure, solutions'[:, 1], solutions'[:, 2], solutions'[:, 3], linewidth=1.5)
    if export_fig
    end
    return figure
end;

function my_2D_phase_portrait(solution_x, solution_y, title, xaxis, yaxis; detailed_fig=false, export_fig=false)
    figure = plot(title=title, xaxis=xaxis, yaxis=yaxis, legend=false, gridalpha=0)
    if detailed_fig
        plot!(background_color_inside="#DDDDDD", gridalpha=1, gridlinewidth=:1, foreground_color_grid="#000000", minorgrid=true, minorgridalpha=0.5, foreground_color_minor_grid=:"#FF0000")
    end
    plot!(figure, solution_x, solution_y, linewidth=1.5)
    if export_fig
    end
    return figure
end;

function my_phase_portraits(solutions, xaxis, yaxis, zaxis)
    solution_x = solutions'[:, 1]
    solution_y = solutions'[:, 2]
    solution_z = solutions'[:, 3]
    fig_a = my_3D_phase_portrait(solutions, "3D Phase Portrait", xaxis, yaxis, zaxis; detailed_fig=true)
    fig_b = my_2D_phase_portrait(solution_x, solution_y, "Species X and Y", xaxis, yaxis; detailed_fig=true)
    fig_c = my_2D_phase_portrait(solution_x, solution_z, "Species X and Z", xaxis, zaxis; detailed_fig=true)
    fig_d = my_2D_phase_portrait(solution_y, solution_z, "Species Y and Z", yaxis, zaxis; detailed_fig=true)
    l = @layout [a b ; c d]
    return plot(fig_a, fig_b, fig_c, fig_d, layout=l)
end;

function my_bifurcation_diagram(data_stable, data_unstable, variable, parameter, color; detailed_fig=false, export_fig=false)
    if data_stable == []
        parameters_stable, solution_stable = [], []
    else
        parameters_stable, solution_stable = data_stable[:, 1], data_stable[:, 2]
    end
    if data_unstable == []
        parameters_unstable, solution_unstable = [], []
    else
        parameters_unstable, solution_unstable = data_unstable[:, 1], data_unstable[:, 2]
    end
    figure = scatter(title="Density of Species $(variable) vs $(parameter)", xaxis="$(parameter)", yaxis="Density of Species $(variable)")
    if detailed_fig
        scatter!(legend=:outertopright, background_color_inside="#DDDDDD", gridalpha=1, gridlinewidth=:1, foreground_color_grid="#000000", minorgrid=true, minorgridalpha=0.5, foreground_color_minor_grid=:"#FF0000")
    end
    scatter!(figure, parameters_stable, solution_stable, markersize=1, markercolor=color, markerstrokecolor=color, label="Stable")
    scatter!(figure, parameters_unstable, solution_unstable, markersize=1, markercolor=:lime, markerstrokecolor=:lime, label="Unstable")
    if export_fig
    end
    return figure
end;

function my_bifurcation_diagrams(solutions, data, variables, bifurcartion_parameter, legend_data)
    x_arr_t, x_arr_s, x_arr_u, y_arr_t, y_arr_s, y_arr_u, z_arr_t, z_arr_s, z_arr_u = data
    x, y, z = variables
    title = "Time Evolution of Each Species"
	xaxis = "Time in Days"
	yaxis = "Population size"
    fig_a = my_plot(solutions, title, xaxis, yaxis, legend_data; detailed_fig=true)
    fig_b = my_bifurcation_diagram(x_arr_s, x_arr_u, x, bifurcartion_parameter, "black"; detailed_fig=true)
    fig_c = my_bifurcation_diagram(y_arr_s, y_arr_u, y, bifurcartion_parameter, "red"; detailed_fig=true)
    fig_d = my_bifurcation_diagram(z_arr_s, z_arr_u, z, bifurcartion_parameter, "blue"; detailed_fig=true)
    l = @layout [a b ; c d]
    return plot(fig_a, fig_b, fig_c, fig_d, layout=l)
end;
    
end