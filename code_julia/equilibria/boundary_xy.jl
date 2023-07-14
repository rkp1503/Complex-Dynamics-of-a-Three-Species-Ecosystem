#=============================================================================
Author: Rayla Kurosaki
GitHub: https://github.com/rkp1503
=============================================================================#

module boundary_xy

import Polynomials: Polynomial, roots

function get_equilibria(parameter_values)
    equilibria = get_analytical_solutions(parameter_values)
    r₂₁, r₃₁, p, ϕ₁₂, ϕ₂₁, ϕ₁₃, u₁, u₂, u₃, u₄ = parameter_values
    β = 11
    for (i, equilibrium) in enumerate(equilibria)
        x, y, z = equilibrium
        # if !(x > 0)
        #     removeat!(equilibria, i)
        # end
        cond = ϕ₂₁ < (β-1)/((ϕ₁₂*β^2+1)^2)
        # cond = ϕ₁₂^2*ϕ₂₁*β^4+2*ϕ₁₂*ϕ₂₁*β^2-β+ϕ₂₁+1 < 0
        if all([cond])
            deleteat!(equilibria, i)
        end
    end
    return equilibria
end;

function get_stable_equilibria!(equilibria, parameter_values)
    r₂₁, r₃₁, p, ϕ₁₂, ϕ₂₁, ϕ₁₃, u₁, u₂, u₃, u₄ = parameter_values
    for (i, equilibrium) in enumerate(equilibria)
        x, y, z = equilibrium
        j₁₁ = 1-2*x+ϕ₁₂*y^2
        j₁₂ = 2*ϕ₁₂*x*y
        j₂₁ = 2*r₂₁*ϕ₂₁*x*y
        j₂₂ = r₂₁*(1-2*y+ϕ₂₁*x^2)
        j₃₃ = r₃₁-u₄+((u₃*(1-p)*y)/(u₂+(1-p)*y))
        C₂ = -j₁₁-j₂₂-j₃₃
        C₁ = j₁₁*j₂₂+j₁₁*j₃₃+j₂₂*j₃₃-j₁₂*j₂₁
        C₀ = j₃₃*(j₁₂*j₂₁-j₁₁*j₂₂)
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
    Y₀ = ϕ₂₁+1
    Y₁ = -1
    Y₂ = 2*ϕ₁₂*ϕ₂₁
    Y₃ = 0
    Y₄ = ϕ₁₂^2*ϕ₂₁
    for root in roots(Polynomial([Y₀, Y₁, Y₂, Y₃, Y₄], :y))
        if (real(root) > 0) && (imag(root) == 0)
            y = real(root)
            x = 1+ϕ₁₂*y^2
            push!(sols_real, (x, y, 0))
        end
    end
    return sols_real
end;

end