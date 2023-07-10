#=============================================================================
Author: Rayla Kurosaki
GitHub: https://github.com/rkp1503
=============================================================================#

module ComputeData

import DifferentialEquations: ODEProblem, solve, Rodas5P
import ModelingToolkit: @named, ODESystem, structural_simplify
import OrderedCollections: OrderedDict

const ROUND = 3
const SOL_MAX = 10
const TOL_SOLUTION = 1e-12
    
function solve_model(model, D, variables_dict, parameters_dict, t, tₘₐₓ)
    vars_keys = collect(keys(variables_dict))
    vars_vals = collect(values(variables_dict))
    params_keys = collect(keys(parameters_dict))
    params_vals = collect(values(parameters_dict))
    @named sys = ODESystem(model(D, vars_keys, params_keys), t, vars_keys, params_keys)
    sys = structural_simplify(sys)
    prob = ODEProblem(sys, vars_vals, (0.0, tₘₐₓ), params_vals, jac = true)
    return solve(prob, Rodas5P(), reltol=TOL_SOLUTION, abstol=TOL_SOLUTION; verbose=false)
end;

function get_bifurcation_data(model, D, variables_dict, parameters, t, tₘₐₓ, bifurcation_parameter, parameter_bounds, parameter_values; data_points=10_000)
	parameters_dict = OrderedDict(zip(parameters, parameter_values))
	x_arr_t, x_arr_s, x_arr_u = [], [], []
	y_arr_t, y_arr_s, y_arr_u = [], [], []
	z_arr_t, z_arr_s, z_arr_u = [], [], []
	prev_x, prev_y, prev_z = "", "", ""
	curr_x, curr_y, curr_z = "", "", ""
	stability_change_values = []
    bifurcation_values = range(parameter_bounds[1], parameter_bounds[2], length=1+data_points)
	for bifurcation_value in bifurcation_values
		parameters_dict[bifurcation_parameter] = bifurcation_value
		solutions = solve_model(model, D, variables_dict, parameters_dict, t, tₘₐₓ)
		n_start = floor(Int, length(solutions.t)*0.9)
		x_min, x_max = extrema(solutions'[:, 1][n_start:end])
		y_min, y_max = extrema(solutions'[:, 2][n_start:end])
		z_min, z_max = extrema(solutions'[:, 3][n_start:end])
		if (abs(x_min-x_max) < 1/data_points)
			curr_x = "Stable"
			if (0 < x_min < SOL_MAX)
				push!(x_arr_t, [bifurcation_value, x_min])
                push!(x_arr_s, [bifurcation_value, x_min])
			end
			if (0 < x_max < SOL_MAX)
				push!(x_arr_t, [bifurcation_value, x_max])
                push!(x_arr_s, [bifurcation_value, x_max])
			end
		else
			curr_x = "Unstable"
			if (0 < x_min < SOL_MAX)
				push!(x_arr_t, [bifurcation_value, x_min])
                push!(x_arr_u, [bifurcation_value, x_min])
			end
			if (0 < x_max < SOL_MAX)
				push!(x_arr_t, [bifurcation_value, x_max])
                push!(x_arr_u, [bifurcation_value, x_max])
			end
		end
		if (abs(y_min-y_max) < 1/data_points)
			curr_y = "Stable"
			if (0 < y_min < SOL_MAX)
				push!(y_arr_t, [bifurcation_value, y_min])
                push!(y_arr_s, [bifurcation_value, y_min])
			end
			if (0 < y_max < SOL_MAX)
				push!(y_arr_t, [bifurcation_value, y_max])
                push!(y_arr_s, [bifurcation_value, y_max])
			end
		else
			curr_y = "Unstable"
			if (0 < y_min < SOL_MAX)
				push!(y_arr_t, [bifurcation_value, y_min])
                push!(y_arr_u, [bifurcation_value, y_min])
			end
			if (0 < y_max < SOL_MAX)
				push!(y_arr_t, [bifurcation_value, y_max])
                push!(y_arr_u, [bifurcation_value, y_max])
			end
		end
		if (abs(z_min-z_max) < 1/data_points)
			curr_z = "Stable"
			if (0 < z_min < SOL_MAX)
				push!(z_arr_t, [bifurcation_value, z_min])
                push!(z_arr_s, [bifurcation_value, z_min])
			end
			if (0 < z_max < SOL_MAX)
				push!(z_arr_t, [bifurcation_value, z_max])
                push!(z_arr_s, [bifurcation_value, z_max])
			end
		else
			curr_z = "Unstable"
			if (0 < z_min < SOL_MAX)
				push!(z_arr_t, [bifurcation_value, z_min])
                push!(z_arr_u, [bifurcation_value, z_min])
			end
			if (0 < z_max < SOL_MAX)
				push!(z_arr_t, [bifurcation_value, z_max])
                push!(z_arr_u, [bifurcation_value, z_max])
			end
		end
		if (curr_x == curr_y == curr_z) && (curr_x != prev_x)
			if prev_x != ""
				println("$(prev_x) -> $(curr_x): $(bifurcation_value)")
				push!(stability_change_values, bifurcation_value)
			end
			prev_x = curr_x
			prev_y = curr_y
			prev_z = curr_z
		end
	end
    x_arr_t = reformat_data(x_arr_t)
    x_arr_s = reformat_data(x_arr_s)
    x_arr_u = reformat_data(x_arr_u)
	y_arr_t = reformat_data(y_arr_t)
    y_arr_s = reformat_data(y_arr_s)
    y_arr_u = reformat_data(y_arr_u)
	z_arr_t = reformat_data(z_arr_t)
    z_arr_s = reformat_data(z_arr_s)
    z_arr_u = reformat_data(z_arr_u)
	# println(stability_change_values)
    return (x_arr_t, x_arr_s, x_arr_u, y_arr_t, y_arr_s, y_arr_u, z_arr_t, z_arr_s, z_arr_u)
end;

function reformat_data(data)
    if data == []
		return []
	else
		return return permutedims(reshape(hcat(data...), (length(data[1]), length(data))))
	end
end;

function approximate_bounds(model, D, variables_dict, parameters_dict, t, tₘₐₓ, equilibrium_type)
    solutions_x, solutions_y, solutions_z = [], [], []
    for _ in 1:MAX_ITER
        parameter_vals = generate_parameters_helper(parameters_dict, equilibrium_type; debug)
        parameters_dict_temp = OrderedDict(zip(keys(parameters_dict), parameter_vals))
        solutions = solve_model(model, D, variables_dict, parameters_dict_temp, t, tₘₐₓ)
        xₛₒₗ, yₛₒₗ, zₛₒₗ = last(solutions.u)
        if !((xₛₒₗ > SOL_MAX) || (yₛₒₗ > SOL_MAX) || (zₛₒₗ > SOL_MAX))
            if compare_results(solutions, parameters_dict_temp, equilibrium_type)
                push!(solutions_x, xₛₒₗ)
                push!(solutions_y, yₛₒₗ)
                push!(solutions_z, zₛₒₗ)
            end
        end
    end
    println("X ∈ $([round(e, digits=ROUND) for e in extrema(solutions_x)])")
    println("Y ∈ $([round(e, digits=ROUND) for e in extrema(solutions_y)])")
    println("Z ∈ $([round(e, digits=ROUND) for e in extrema(solutions_z)])")
    return nothing
end;

function get_common_bounds(data)
	a = findmax([findmin(data[1][:, 1])[1], findmin(data[4][:, 1])[1], findmin(data[7][:, 1])[1]])[1]
	b = findmin([findmax(data[1][:, 1])[1], findmax(data[4][:, 1])[1], findmax(data[7][:, 1])[1]])[1]
    println("Full Bounds: [$(a), $(b)]")
	a = findmax([findmin(data[2][:, 1])[1], findmin(data[5][:, 1])[1], findmin(data[8][:, 1])[1]])[1]
	b = findmin([findmax(data[2][:, 1])[1], findmax(data[5][:, 1])[1], findmax(data[8][:, 1])[1]])[1]
    println("Stable Bounds: [$(a), $(b)]")
	a = findmax([findmin(data[3][:, 1])[1], findmin(data[6][:, 1])[1], findmin(data[9][:, 1])[1]])[1]
	b = findmin([findmax(data[3][:, 1])[1], findmax(data[6][:, 1])[1], findmax(data[9][:, 1])[1]])[1]
    println("Unstable Bounds: [$(a), $(b)]")
    return nothing
end;


function compute_bifurcation_points(data)
	
end;

end