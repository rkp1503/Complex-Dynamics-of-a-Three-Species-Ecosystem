#=============================================================================
Author: Ramsey (Rayla) Phuc
Alias: Rayla Kurosaki
GitHub: https://github.com/rkp1503
Co-author: Ephraim Agyingi
=============================================================================#

module utils

import DifferentialEquations: ODEProblem, solve, Rodas5P
import ModelingToolkit: @named, ODESystem, structural_simplify
import OrderedCollections: OrderedDict
import Plots: plot, plot!

include("equilibria/axial_z.jl")
import .axial_z
include("equilibria/boundary_xy.jl")
import .boundary_xy
include("equilibria/boundary_xz.jl")
import .boundary_xz
include("equilibria/boundary_yz.jl")
import .boundary_yz
include("equilibria/interior.jl")
import .interior

const ITER_PARAMS = 1_000_000
const ITER_DEBUG = 1_000_000
const ITER = 1_000
const MAX_VAL = 1e6
const ROUND = 3
const TOL_1 = 1e-12
const TOL_2 = 1e-3

function solve_model(model, D, vars_dict, params_dict, t, tₘₐₓ)
    vars_keys = collect(keys(vars_dict))
    vars_vals = collect(values(vars_dict))
    params_keys = collect(keys(params_dict))
    params_vals = collect(values(params_dict))
    @named sys = ODESystem(model(D, vars_keys, params_keys), t, vars_keys, params_keys)
    sys = structural_simplify(sys)
    prob = ODEProblem(sys, vars_vals, (0.0, tₘₐₓ), params_vals, jac = true)
    return solve(prob, Rodas5P(), reltol=TOL_1, abstol=TOL_1; verbose=false)
end;

function my_plot(sol, title, xaxis, yaxis, legend_dict; detailed=false, param_vals=[], normalize=1)
    figure = plot(title=title, xaxis=xaxis, yaxis=yaxis,
        # legend=:outertopright,
        # background_color_inside="#DDDDDD",
        # gridalpha=1,
        # gridlinewidth=:1,
        # foreground_color_grid="#000000",
        # minorgrid=true,
        # minorgridalpha=0.5,
        # foreground_color_minor_grid=:"#FF0000",
        gridalpha=0,
        )
    for (i, style) in enumerate(collect(values(legend_dict)))
        label, color = style
        plot!(figure, sol.t, sol'[:, i]/normalize, 
            label=label, linecolor=color, linewidth=1.5
            )
        if detailed
            upper_bounds = get_upper_bounds(param_vals)
            plot!(figure, sol.t, [upper_bounds[i] for _ in sol.t]/normalize, label="", linecolor=color, linewidth=1.5, linestyle=:dash)
        end
    end
    return figure
end;

function my_3D_phase_portrait(sol, title, xaxis, yaxis, zaxis)
    figure = plot(title=title, xaxis=xaxis, yaxis=yaxis, zaxis=zaxis, legend=false,
        # background_color_inside="#DDDDDD",
        # gridalpha=1,
        # gridlinewidth=:1,
        # foreground_color_grid="#000000",
        # minorgrid=true,
        # minorgridalpha=0.5,
        # foreground_color_minor_grid=:"#FF0000",
        gridalpha=0,
        )
    plot!(figure, sol'[:, 1], sol'[:, 2], sol'[:, 3], 
        linewidth=1.5
        )
    return figure
end;

function my_2D_phase_portrait(sol_x, sol_y, title, xaxis, yaxis)
    figure = plot(title=title, xaxis=xaxis, yaxis=yaxis, legend=false,
        # background_color_inside="#DDDDDD",
        # gridalpha=1,
        # gridlinewidth=:1,
        # foreground_color_grid="#000000",
        # minorgrid=true,
        # minorgridalpha=0.5,
        # foreground_color_minor_grid=:"#FF0000",
        gridalpha=0,
        )
    plot!(figure, sol_x, sol_y, 
        linewidth=1.5
        )
    return figure
end;

function my_phase_portrait(sol, title, xaxis, yaxis, variables; zaxis="")
    figure = plot(title=title, xaxis=xaxis, yaxis=yaxis, zaxis=zaxis, legend=false,
        # background_color_inside="#DDDDDD",
        # gridalpha=1,
        # gridlinewidth=:1,
        # foreground_color_grid="#000000",
        # minorgrid=true,
        # minorgridalpha=0.5,
        # foreground_color_minor_grid=:"#FF0000",
        gridalpha=0,
        )
    plot!(figure, sol, idxs=variables, 
        linewidth=1.5
        )
    return figure
end;

function get_upper_bounds(param_vals)
    r₂₁, r₃₁, p, ϕ₁₂, ϕ₂₁, ϕ₁₃, u₁, u₂, u₃, u₄ = param_vals
    x_ub = 0
    y_ub = 0
    z_ub = 1+((u₃-u₄)/r₃₁)
    return [x_ub, y_ub, z_ub]
end;


function generate_parameters(model, D, vars_dict, params_dict, t, tₘₐₓ, eq_type; debug=false)
    sol_too_large = 0
    numerical_differ_analytical = 0
    iteration = 1
    while true
        println("Iteration: $(iteration)")
        params_keys = collect(keys(params_dict))
        params_vals = generate_parameters_helper(params_dict, eq_type)
        params_dict_temp = OrderedDict(zip(params_keys, params_vals))
        sol = solve_model(model, D, vars_dict, params_dict_temp, t, tₘₐₓ)
        xₛₒₗ, yₛₒₗ, zₛₒₗ = last(sol.u)
        if !((xₛₒₗ > MAX_VAL) || (yₛₒₗ > MAX_VAL) || (zₛₒₗ > MAX_VAL))
            if compare_results(sol, params_dict_temp, eq_type)
                if !debug
                    println("Iteration: $(iteration)")
                    println("Parameters: $(collect(keys(params_dict)))")
                    println("Parameter values: $(params_vals)")
                    print_error(sol, params_dict_temp, eq_type)
                    return params_vals
                end
            else
                numerical_differ_analytical += 1
            end
        else
            sol_too_large += 1
        end
        iteration += 1
        if debug
            if iteration == ITER_DEBUG
                break 
            end
        else
            if iteration == ITER_PARAMS
                break 
            end
        end
    end
    if debug
        return sol_too_large, numerical_differ_analytical
    else
        println("Failed to find sutible parameters after $(iteration) iterations.")
        return nothing
    end
end;


function generate_parameters_helper(params_dict, eq_type)
    while true
        param_vals = [round(2*rand(), digits=ROUND) for _ in 1:length(params_dict)]
        if (0 ∉ param_vals) && (param_vals[3] < 1)
            if eq_exist(param_vals, eq_type)
                if eq_stable(param_vals, eq_type)
                    return param_vals
                end
            end
        end
    end
end;

function eq_exist(parameter_values, equilibrium_type)
    if equilibrium_type == "z axial"
        return axial_z.equilibrium_exist(parameter_values)
    elseif equilibrium_type == "xy boundary"
        return boundary_xy.equilibrium_exist(parameter_values)
    elseif equilibrium_type == "xz boundary"
        return boundary_xz.equilibrium_exist(parameter_values)
    elseif equilibrium_type == "yz boundary"
        return boundary_yz.equilibrium_exist(parameter_values)
    else
        return interior.equilibrium_exist(parameter_values)
    end
end;

function eq_stable(parameter_values, equilibrium_type)
    if equilibrium_type == "z axial"
        return axial_z.equilibrium_stable(parameter_values)
    elseif equilibrium_type == "xy boundary"
        return boundary_xy.equilibrium_stable(parameter_values)
    elseif equilibrium_type == "xz boundary"
        return boundary_xz.equilibrium_stable(parameter_values)
    elseif equilibrium_type == "yz boundary"
        return boundary_yz.equilibrium_stable(parameter_values)
    else
        return interior.equilibrium_stable(parameter_values)
    end
end;

function compare_results(solution, parameters_dict, equilibrium_type)
    xₛₒₗ, yₛₒₗ, zₛₒₗ = last(solution.u)
    E_x, E_y, E_z = get_analytical_solution(parameters_dict, equilibrium_type)
    x_err = round(abs(E_x-xₛₒₗ), digits=ROUND)
    y_err = round(abs(E_y-yₛₒₗ), digits=ROUND)
    z_err = round(abs(E_z-zₛₒₗ), digits=ROUND)
    return (x_err <= TOL_2) && (y_err <= TOL_2) && (z_err <= TOL_2)
end;

function get_analytical_solution(parameters_dict, equilibrium_type)
    param_vals = collect(values(parameters_dict))
    if equilibrium_type == "z axial"
        return axial_z.analytical_solution(param_vals)
    elseif equilibrium_type == "xy boundary"
        return boundary_xy.analytical_solution(param_vals)
    elseif equilibrium_type == "xz boundary"
        return boundary_xz.analytical_solution(param_vals)
    elseif equilibrium_type == "yz boundary"
        return boundary_yz.analytical_solution(param_vals)
    else
        return interior.analytical_solution(param_vals)
    end
end;
		
function print_error(solution, parameters_dict, equilibrium_type)
    xₛₒₗ, yₛₒₗ, zₛₒₗ = last(solution.u)
    println("Numerical: ($(round(xₛₒₗ, digits=ROUND)), $(round(yₛₒₗ, digits=ROUND)), $(round(zₛₒₗ, digits=ROUND)))")
    E_x, E_y, E_z = get_analytical_solution(parameters_dict, equilibrium_type)
    println("Analytical: ($(round(E_x, digits=ROUND)), $(round(E_y, digits=ROUND)), $(round(E_z, digits=ROUND)))")
end;

function test_accuracy(model, D, vars_dict, params_dict, t, tₘₐₓ, eq_type)
    debug_stuff = utils.generate_parameters(model, D, vars_dict, params_dict, t, tₘₐₓ, eq_type, debug=true)
    sol_too_large, numerical_differ_analytical = debug_stuff
    println("Equilibrium: $(eq_type)")
    println("Number of Iterations: $(ITER_DEBUG)")
    println("sol_too_large: $(sol_too_large)")
    println("numerical_differ_analytical: $(numerical_differ_analytical)")
    success_rate = 1-(sol_too_large+numerical_differ_analytical)/ITER_DEBUG
    println("Success Rate: $(success_rate)")
end

end
