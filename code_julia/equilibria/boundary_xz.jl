#=============================================================================
Author: Rayla Kurosaki
GitHub: https://github.com/rkp1503
=============================================================================#

module boundary_xz

function get_equilibria(parameter_values)
    equilibria = get_analytical_solutions(parameter_values)
    r₂₁, r₃₁, p, ϕ₁₂, ϕ₂₁, ϕ₁₃, u₁, u₂, u₃, u₄ = parameter_values
    for (i, equilibrium) in enumerate(equilibria)
        x, y, z = equilibrium
        # if !((x > 0) && (z > 0))
        #     removeat!(equilibria, i)
        # end
        cond1 = (u₄/r₃₁)+(1/ϕ₁₃) > 1
        cond2 = r₃₁ > u₄
        if !all([cond1, cond2])
            deleteat!(equilibria, i)
        end
    end
    return equilibria
end;

function get_stable_equilibria!(equilibria, parameter_values)
    r₂₁, r₃₁, p, ϕ₁₂, ϕ₂₁, ϕ₁₃, u₁, u₂, u₃, u₄ = parameter_values
    for (i, equilibrium) in enumerate(equilibria)
        x, y, z = equilibrium
        cond1 = 1-(u₄/r₃₁) > (1-2*(1-ϕ₁₃*(1-(u₄/r₃₁))))/ϕ₁₃
        cond2 = 1-(u₄/r₃₁) > r₂₁*u₂*(1+ϕ₂₁*(1-ϕ₁₃*(1-(u₄/r₃₁)))^2)/(u₁*(1-p))
        cond3 = 1-(u₄/r₃₁) > (r₃₁-u₄)/(2*r₃₁)
        if !all([cond1, cond2, cond3])
            deleteat!(equilibria, i)
        end
    end
    return nothing
end;

function get_analytical_solutions(parameter_values)
    r₂₁, r₃₁, p, ϕ₁₂, ϕ₂₁, ϕ₁₃, u₁, u₂, u₃, u₄ = parameter_values
    x = 1-ϕ₁₃*(1-(u₄/r₃₁))
    z = 1-(u₄/r₃₁)
    return [(x, 0, z)]
end;
    
end