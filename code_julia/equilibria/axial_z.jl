#=============================================================================
Author: Rayla Kurosaki
GitHub: https://github.com/rkp1503
=============================================================================#

module axial_z

function get_equilibria(parameter_values)
    equilibria = get_analytical_solutions(parameter_values)
    r₂₁, r₃₁, p, ϕ₁₂, ϕ₂₁, ϕ₁₃, u₁, u₂, u₃, u₄ = parameter_values
    for (i, equilibrium) in enumerate(equilibria)
        x, y, z = equilibrium
        # if !(z > 0)
        #     removeat!(equilibria, i)
        # end
        cond = r₃₁ > u₄
        if !all([cond])
            deleteat!(equilibria, i)
        end
    end
    return equilibria
end;

function get_stable_equilibria!(equilibria, parameter_values)
    r₂₁, r₃₁, p, ϕ₁₂, ϕ₂₁, ϕ₁₃, u₁, u₂, u₃, u₄ = parameter_values
    for (i, equilibrium) in enumerate(equilibria)
        x, y, z = equilibrium
        cond1 = (u₄/r₃₁)+(1/ϕ₁₃) < 1
        cond2 = (u₄/r₃₁)+((r₂₁*u₂)/(u₁*(1-p))) < 1
        cond3 = (u₄/r₃₁)<(1/2)
        if !all([cond1, cond2, cond3])
            deleteat!(equilibria, i)
        end
    end
    return nothing
end;

function get_analytical_solutions(parameter_values)
    r₂₁, r₃₁, p, ϕ₁₂, ϕ₂₁, ϕ₁₃, u₁, u₂, u₃, u₄ = parameter_values
    z = 1-(u₄/r₃₁)
    return [(0, 0, z)]
end;
    
end