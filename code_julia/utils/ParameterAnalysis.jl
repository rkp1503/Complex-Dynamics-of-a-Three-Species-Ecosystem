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
        parameter_values = generate_parameters_helper(parameters_dict, equilibrium_type; debug)
        if !isnothing(parameter_values)
            parameters_dict_temp = OrderedDict(zip(collect(keys(parameters_dict)), parameter_values))
            solutions = ComputeData.solve_model(model, D, variables_dict, parameters_dict_temp, t, tₘₐₓ)
            xₛₒₗ, yₛₒₗ, zₛₒₗ = last(solutions.u)
            # if !((xₛₒₗ > SOL_MAX) || (yₛₒₗ > SOL_MAX) || (zₛₒₗ > SOL_MAX))
                if compare_results(solutions, parameter_values, equilibrium_type)
                    println("Sutible parameters found at iteration $(iteration):")
                    println("Parameter Values: $(parameter_values)")
                    println("Equilibrium:")
                    println("\tX: $(round(xₛₒₗ, digits=ROUND))")
                    println("\tY: $(round(yₛₒₗ, digits=ROUND))")
                    println("\tZ: $(round(zₛₒₗ, digits=ROUND))")
                    return parameter_values
                else
                    if debug
                        println("True and approximated solutions are not close enough!\n")
                    end
                end
            # else
            #     if debug
            #         println("Solution is too large!\n")
            #     end
            # end
        end
    end
    println("Failed to find sutible parameters after $(iteration) iterations.")
    return nothing
end;

function generate_parameters_helper(parameters_dict, equilibrium_type; debug=false)
    while true
        parameter_values = [round(2*rand(), digits=ROUND) for _ in 1:length(parameters_dict)]
        if (0 ∉ parameter_values) && (parameter_values[3] < 1)
            equilibria = get_equilibria(parameter_values, equilibrium_type)
            if debug
                println("Number of Equilibria: $(length(equilibria))")
            end
            if !isempty(equilibria)
                get_stable_equilibria!(equilibria, parameter_values, equilibrium_type)
                if debug
                    println("Number of Stable Equilibria: $(length(equilibria))")
                end
                if !isempty(equilibria)
                    return parameter_values
                else
                    if debug
                        println("The $(equilibrium_type) equilibrium is not stable for the parameters $(parameter_values)\n")
                    end
                end
            else
                if debug
                    println("The $(equilibrium_type) equilibrium does not exist for the parameters $(parameter_values)\n")
                end
            end
            return nothing
        end
    end
end;

function get_equilibria(parameter_values, equilibrium_type)
    if equilibrium_type == "z-axial"
        return axial_z.get_equilibria(parameter_values)
    elseif equilibrium_type == "xy-boundary"
        return boundary_xy.get_equilibria(parameter_values)
    elseif equilibrium_type == "xz-boundary"
        return boundary_xz.get_equilibria(parameter_values)
    elseif equilibrium_type == "yz-boundary"
        return boundary_yz.get_equilibria(parameter_values)
    else
        return interior.get_equilibria(parameter_values)
    end
end;

function get_stable_equilibria!(equilibria_lst, parameter_values, equilibrium_type)
    if equilibrium_type == "z-axial"
        return axial_z.get_stable_equilibria!(equilibria_lst, parameter_values)
    elseif equilibrium_type == "xy-boundary"
        return boundary_xy.get_stable_equilibria!(equilibria_lst, parameter_values)
    elseif equilibrium_type == "xz-boundary"
        return boundary_xz.get_stable_equilibria!(equilibria_lst, parameter_values)
    elseif equilibrium_type == "yz-boundary"
        return boundary_yz.get_stable_equilibria!(equilibria_lst, parameter_values)
    else
        return interior.get_stable_equilibria!(equilibria_lst, parameter_values)
    end
end;

function compare_results(solutions, parameter_values, equilibrium_type)
    xₛₒₗ, yₛₒₗ, zₛₒₗ = last(solutions.u)
    equilibria = get_equilibria(parameter_values, equilibrium_type)
    get_stable_equilibria!(equilibria, parameter_values, equilibrium_type)
    E_x, E_y, E_z = equilibria[1]
    x_err = round(abs(E_x-xₛₒₗ), digits=ROUND)
    y_err = round(abs(E_y-yₛₒₗ), digits=ROUND)
    z_err = round(abs(E_z-zₛₒₗ), digits=ROUND)
    return (x_err <= TOL) && (y_err <= TOL) && (z_err <= TOL)
end;

function validate_parameters(parameter_values, equilibrium_type; debug=false)
    equilibria = get_equilibria(parameter_values, equilibrium_type)
    if debug
        println("Number of Equilibria: $(length(equilibria))")
    end
    if isempty(equilibria)
        println("The $(equilibrium_type) equilibrium does not exist for the parameters $(parameter_values)")
    end
    get_stable_equilibria!(equilibria, parameter_values, equilibrium_type)
    if debug
        println("Number of Stable Equilibria: $(length(equilibria))")
    end
    if isempty(equilibria)
        println("The $(equilibrium_type) equilibrium is not stable for the parameters $(parameter_values)")
    end
    println("The $(equilibrium_type) equilibrium exists and is stable for the parameters $(parameter_values)")
    return nothing
end;
    
end