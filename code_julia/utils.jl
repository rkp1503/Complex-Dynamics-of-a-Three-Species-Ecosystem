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

# const ITER_PARAMS = 1_000_000
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
    # return solve(prob, reltol=TOL_1, abstol=TOL_1)
    return solve(prob, Rodas5P(), reltol=TOL_1, abstol=TOL_1)
end;

function my_plot(sol, title, xaxis, yaxis, legend_dict; normalize=1)
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

function success_rate(model, D, vars_dict, params_dict, t, tₘₐₓ, eq_type; method="eig")
    sol_large_c = 0
    err_c = 0
    for i in 1:ITER
        params_keys = collect(keys(params_dict))
        params_vals = generate_sutible_parameters(params_dict, eq_type, method)
        params_dict_temp = OrderedDict(zip(params_keys, params_vals))
        sol = solve_model(model, D, vars_dict, params_dict_temp, t, tₘₐₓ)
        if !compare_results(sol, params_dict_temp, eq_type)
            # print_error(sol, params_dict_temp, eq_type)
            xₛₒₗ, yₛₒₗ, zₛₒₗ = last(sol.u)
            if (xₛₒₗ > MAX_VAL) || (yₛₒₗ > MAX_VAL) || (zₛₒₗ > MAX_VAL)
                sol_large_c += 1
            else
                err_c += 1
            end
        end
    end
    # Probability of true solutions. (Considers very large solutions as not true solutions)
    proba_a = 1-((err_c+sol_large_c)/ITER)
    # Probability of valid solutions. (Considers very large solutions as true solutions)
    proba_b = 1-(err_c/ITER)
    println("Success rate (True solutions): $(100*proba_a)%")
    println("Success rate (Valid solutions): $(100*proba_b)%")
    println("Success rate error: $(100*abs(proba_a-proba_b))%")
    println("")
end;

function generate_sutible_parameters(params_dict, eq_type, method="eig")
    while true
        param_vals = [round(rand(), digits=ROUND) for x in 1:length(params_dict)]
        if (0 ∉ param_vals) && (1 ∉ param_vals)
            if eq_exist(param_vals, eq_type)
                if eq_stable(param_vals, eq_type, method)
                    return param_vals
                end
            end
        end
    end
end;

function generate_valid_parameters(model, D, vars_dict, params_dict, t, tₘₐₓ, eq_type; method="eig")
    iteration_counter = 1
    while true
        params_keys = collect(keys(params_dict))
        params_vals = generate_sutible_parameters(params_dict, eq_type, method)
        params_dict_temp = OrderedDict(zip(params_keys, params_vals))
        sol = solve_model(model, D, vars_dict, params_dict_temp, t, tₘₐₓ)
        xₛₒₗ, yₛₒₗ, zₛₒₗ = last(sol.u)
        if !((xₛₒₗ > MAX_VAL) || (yₛₒₗ > MAX_VAL) || (zₛₒₗ > MAX_VAL))
            if compare_results(sol, params_dict_temp, eq_type)
                # println("Valid parameters found at iteration $(iteration_counter)")
                return params_vals
            end
        end
        iteration_counter += 1
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

function eq_stable(parameter_values, equilibrium_type, method)
    if equilibrium_type == "z axial"
        return axial_z.equilibrium_stable(parameter_values, method)
    elseif equilibrium_type == "xy boundary"
        return boundary_xy.equilibrium_stable(parameter_values, method)
    elseif equilibrium_type == "xz boundary"
        return boundary_xz.equilibrium_stable(parameter_values, method)
    elseif equilibrium_type == "yz boundary"
        return boundary_yz.equilibrium_stable(parameter_values, "RHC")
    else
        return interior.equilibrium_stable(parameter_values, "RHC")
    end
end;

function compare_results(solution, parameters_dict, equilibrium_type)
    xₛₒₗ, yₛₒₗ, zₛₒₗ = last(solution.u)
    E_x, E_y, E_z = get_analytical_solution(parameters_dict, equilibrium_type)
    x_err = abs(round(E_x, digits=ROUND) - round(xₛₒₗ, digits=ROUND))
    y_err = abs(round(E_y, digits=ROUND) - round(yₛₒₗ, digits=ROUND))
    z_err = abs(round(E_z, digits=ROUND) - round(zₛₒₗ, digits=ROUND))
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
    println("")
end;

function analyze(model, D, vars_dict, params_dict, t, tₘₐₓ, eq_type; method="eig")
    trials = []
    for i in 1:ITER
        j = 1
        while true
            params_keys = collect(keys(params_dict))
            params_vals = generate_sutible_parameters(params_dict, eq_type, method)
            params_dict_temp = OrderedDict(zip(params_keys, params_vals))
            sol = solve_model(model, D, vars_dict, params_dict_temp, t, tₘₐₓ)
            if compare_results(sol, params_dict_temp, eq_type)
                xₛₒₗ, yₛₒₗ, zₛₒₗ = last(sol.u)
                if !((xₛₒₗ > MAX_VAL) || (yₛₒₗ > MAX_VAL) || (zₛₒₗ > MAX_VAL))
                    push!(trials, j)
                    j = 1
                    break
                end
            end
            j += 1
        end
    end
    println("Expected number of iterations before true solution is found: $(sum(trials)/length(trials))")
    println("")
end;

end
