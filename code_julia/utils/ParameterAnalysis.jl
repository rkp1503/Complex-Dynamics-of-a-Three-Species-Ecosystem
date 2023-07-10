#=============================================================================
Author: Rayla Kurosaki
GitHub: https://github.com/rkp1503
=============================================================================#

module ParameterAnalysis

import OrderedCollections: OrderedDict

include("../equilibria/axial_z.jl")
import .axial_z
include("../equilibria/boundary_xy.jl")
import .boundary_xy
include("../equilibria/boundary_xz.jl")
import .boundary_xz
include("../equilibria/boundary_yz.jl")
import .boundary_yz
include("../equilibria/interior.jl")
import .interior
include("ComputeData.jl")
import .ComputeData

const MAX_ITER = 1e6
const ROUND = 3
const SOL_MAX = 10
const TOL = 1e-3

function generate_parameters(model, D, variables_dict, parameters_dict, t, tₘₐₓ, equilibrium_type; debug=false)
    for iteration in 1:MAX_ITER
        parameter_vals = generate_parameters_helper(parameters_dict, equilibrium_type; debug)
        parameters_dict_temp = OrderedDict(zip(collect(keys(parameters_dict)), parameter_vals))
        solutions = ComputeData.solve_model(model, D, variables_dict, parameters_dict_temp, t, tₘₐₓ)
        xₛₒₗ, yₛₒₗ, zₛₒₗ = last(solutions.u)
        if !((xₛₒₗ > SOL_MAX) || (yₛₒₗ > SOL_MAX) || (zₛₒₗ > SOL_MAX))
            if compare_results(solutions, parameters_dict_temp, equilibrium_type)
                println("Sutible parameters found at iteration $(iteration):\n$(parameter_vals)")
                print_solutions(solutions, parameters_dict_temp, equilibrium_type)
                return parameter_vals
            end
        end
    end
    println("Failed to find sutible parameters after $(iteration) iterations.")
    return nothing
end;

function generate_parameters_helper(parameters_dict, equilibrium_type; debug=false)
    for _ in 1:MAX_ITER
        parameter_values = [round(2*rand(), digits=ROUND) for _ in 1:length(parameters_dict)]
        if (0 ∉ parameter_values) && (parameter_values[3] < 1)
            if eq_exist(parameter_values, equilibrium_type)
                if eq_stable(parameter_values, equilibrium_type)
                    return parameter_values
                else
                    if debug
                        println("The $(equilibrium_type) equilibrium is not stable for the parameters $(parameter_values)")
                    end
                end
            else
                if debug
                    println("The $(equilibrium_type) equilibrium does not exist for the parameters $(parameter_values)")
                end
            end
        end
    end
end;

function eq_exist(parameter_values, equilibrium_type)
    if equilibrium_type == "z-axial"
        return axial_z.equilibrium_exist(parameter_values)
    elseif equilibrium_type == "xy-boundary"
        return boundary_xy.equilibrium_exist(parameter_values)
    elseif equilibrium_type == "xz-boundary"
        return boundary_xz.equilibrium_exist(parameter_values)
    elseif equilibrium_type == "yz-boundary"
        return boundary_yz.equilibrium_exist(parameter_values)
    else
        return interior.equilibrium_exist(parameter_values)
    end
end;

function eq_stable(parameter_values, equilibrium_type)
    if equilibrium_type == "z-axial"
        return axial_z.equilibrium_stable(parameter_values)
    elseif equilibrium_type == "xy-boundary"
        return boundary_xy.equilibrium_stable(parameter_values)
    elseif equilibrium_type == "xz-boundary"
        return boundary_xz.equilibrium_stable(parameter_values)
    elseif equilibrium_type == "yz-boundary"
        return boundary_yz.equilibrium_stable(parameter_values)
    else
        return interior.equilibrium_stable(parameter_values)
    end
end;

function compare_results(solutions, parameters_dict_temp, equilibrium_type)
    xₛₒₗ, yₛₒₗ, zₛₒₗ = last(solutions.u)
    E_x, E_y, E_z = get_analytical_solution(parameters_dict_temp, equilibrium_type)
    x_err = round(abs(E_x-xₛₒₗ), digits=ROUND)
    y_err = round(abs(E_y-yₛₒₗ), digits=ROUND)
    z_err = round(abs(E_z-zₛₒₗ), digits=ROUND)
    return (x_err <= TOL) && (y_err <= TOL) && (z_err <= TOL)
end;

function validate_parameters(parameter_values, equilibrium_type)
    if !eq_exist(parameter_values, equilibrium_type)
        println("The $(equilibrium_type) equilibrium does not exist for the parameters $(parameter_values)")
        return nothing
    end
    if !eq_stable(parameter_values, equilibrium_type)
        println("The $(equilibrium_type) equilibrium is not stable for the parameters $(parameter_values)")
        return nothing
    end
    println("The $(equilibrium_type) equilibrium exists and is stable for the parameters $(parameter_values)")
    return nothing
end;

function print_solutions(solutions, parameters_dict, equilibrium_type)
    if compare_results(solutions, parameters_dict, equilibrium_type)
        xₛₒₗ, yₛₒₗ, zₛₒₗ = last(solutions.u)
        println("Equilibrium:\n\tX: $(round(xₛₒₗ, digits=ROUND))\n\tY: $(round(yₛₒₗ, digits=ROUND))\n\tZ: $(round(zₛₒₗ, digits=ROUND))")
    end
    return nothing
end;

function get_analytical_solution(parameters_dict, equilibrium_type)
    parameter_vals = collect(values(parameters_dict))
    if equilibrium_type == "z-axial"
        return axial_z.analytical_solution(parameter_vals)
    elseif equilibrium_type == "xy-boundary"
        return boundary_xy.analytical_solution(parameter_vals)
    elseif equilibrium_type == "xz-boundary"
        return boundary_xz.analytical_solution(parameter_vals)
    elseif equilibrium_type == "yz-boundary"
        return boundary_yz.analytical_solution(parameter_vals)
    else
        return interior.analytical_solution(parameter_vals)
    end
end;

function accuracy_test(model, D, variables_dict, parameters_dict, t, tₘₐₓ, equilibrium_type)
    sol_too_large = 0
    numerical_differ_analytical = 0
    for _ in 1:MAX_ITER
        parameter_vals = generate_parameters_helper(parameters_dict, equilibrium_type; debug)
        parameters_dict_temp = OrderedDict(zip(collect(keys(parameters_dict)), parameter_vals))
        solutions = solve_model(model, D, variables_dict, parameters_dict_temp, t, tₘₐₓ)
        xₛₒₗ, yₛₒₗ, zₛₒₗ = last(solutions.u)
        if ((xₛₒₗ > SOL_MAX) || (yₛₒₗ > SOL_MAX) || (zₛₒₗ > SOL_MAX))
            sol_too_large += 1
            continue
        end
        if !compare_results(solutions, parameters_dict_temp, equilibrium_type)
            numerical_differ_analytical += 1
            continue
        end
    end
    println("Equilibrium: $(equilibrium_type)")
    success_rate = 1-(sol_too_large+numerical_differ_analytical)/MAX_ITER
    println("Success Rate: $(success_rate)")
    return nothing
end;
    
end