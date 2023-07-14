#=============================================================================
Author: Rayla Kurosaki
GitHub: https://github.com/rkp1503
=============================================================================#

module boundary_yz

import Polynomials: Polynomial, roots

function get_equilibria(parameter_values)
    equilibria = get_analytical_solutions(parameter_values)
    r₂₁, r₃₁, p, ϕ₁₂, ϕ₂₁, ϕ₁₃, u₁, u₂, u₃, u₄ = parameter_values
    for (i, equilibrium) in enumerate(equilibria)
        x, y, z = equilibrium
        # if !(z > 0)
        #     removeat!(equilibria, i)
        # end
        cond1 = y > (u₂*(u₄-r₃₁))/((u₃-u₄+r₃₁)*(1-p))
        cond2 = (r₂₁*r₃₁*u₂)/(u₁*(r₃₁-u₄)*(1-p)) > 1
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
        j₁₁ = 1+ϕ₁₂*y^2-ϕ₁₃*z
        j₂₂ = r₂₁*(1-2*y)-((u₁*u₂*(1-p)*z)/(u₂+(1-p)*y)^2)
        j₂₃ = -((u₁*(1-p)*y)/(u₂+(1-p)*y))
        j₃₂ = ((u₂*u₃*(1-p)*z)/(u₂+(1-p)*y)^2)
        j₃₃ = r₃₁*(1-2*z)-u₄+((u₃*(1-p)*y)/(u₂+(1-p)*y))
        C₂ = -j₁₁-j₂₂-j₃₃
        C₁ = j₁₁*j₂₂+j₁₁*j₃₃+j₂₂*j₃₃-j₂₃*j₃₂
        C₀ = j₁₁*(j₂₃*j₃₂-j₂₂*j₃₃)
        cond1 = C₂ > 0
        cond2 = C₁ > 0
        cond3 = C₀ > 0
        cond4 = C₂*C₁ > C₀
        if !all([cond1, cond2, cond3, cond4])
            deleteat!(equilibria, i)
        end
    end
    return nothing
end;

function get_analytical_solutions(parameter_values)
    sols_real = []
    r₂₁, r₃₁, p, ϕ₁₂, ϕ₂₁, ϕ₁₃, u₁, u₂, u₃, u₄ = parameter_values
    Y₀ = u₂*(r₂₁*r₃₁*u₂+u₁*(u₄-r₃₁)*(1-p))
    Y₁ = r₂₁*r₃₁*u₂*(2*(1-p)-u₂)+u₁*(u₄-u₃-r₃₁)*(1-p)^2
    Y₂ = r₂₁*r₃₁*((1-p)-2*u₂)*(1-p)
    Y₃ = -r₂₁*r₃₁*(1-p)^2
    if Y₀ > 0
        for root in roots(Polynomial([Y₀, Y₁, Y₂, Y₃], :y))
            if (real(root) > 0) && (imag(root) == 0)
                y = real(root)
                z = 1+((((u₃*(1-p)*y)/(u₂+(1-p)*y))-u₄)/r₃₁)
                push!(sols_real, (0, y, z))
            end
        end
    end
    return sols_real
end;

end
