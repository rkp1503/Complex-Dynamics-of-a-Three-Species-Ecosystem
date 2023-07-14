#=============================================================================
Author: Rayla Kurosaki
GitHub: https://github.com/rkp1503
=============================================================================#

module interior

import Polynomials: Polynomial, roots

function get_equilibria(parameter_values)
    equilibria = get_analytical_solutions(parameter_values)
    r₂₁, r₃₁, p, ϕ₁₂, ϕ₂₁, ϕ₁₃, u₁, u₂, u₃, u₄ = parameter_values
    for (i, equilibrium) in enumerate(equilibria)
        x, y, z = equilibrium
        # if !((x > 0) && (z > 0))
        #     deleteat!(equilibria, i)
        # end
        cond1 = (1+ϕ₁₂*y^2)/ϕ₁₃ > z
        cond2 = y > (u₂*(u₄-r₃₁))/((u₃-u₄+r₃₁)*(1-p))
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
        j₁₁ = 1-2*x+ϕ₁₂*y^2-ϕ₁₃*z
        j₁₂ = 2*ϕ₁₂*x*y
        j₁₃ = -ϕ₁₃*x
        j₂₁ = 2*ϕ₂₁*r₂₁*x*y
        j₂₂ = r₂₁*(1-2*y+ϕ₂₁*x^2)-((u₁*u₂*(1-p)*z)/(u₂+(1-p)*y)^2)
        j₂₃ = -((u₁*(1-p)*y)/(u₂+(1-p)*y))
        j₃₂ = (u₂*u₃*(1-p)*z)/((u₂+(1-p)*y)^2)
        j₃₃ = r₃₁*(1-2*z)-u₄+((u₃*(1-p)*y)/(u₂+(1-p)*y))
        C₂ = -j₁₁-j₂₂-j₃₃
        C₁ = j₁₁*j₂₂+j₁₁*j₃₃+j₂₂*j₃₃-j₁₂*j₂₁-j₂₃*j₃₂
        C₀ = (j₁₁*j₂₃*j₃₂+j₁₂*j₂₁*j₃₃)-(j₁₁*j₂₂*j₃₃+j₁₃*j₂₁*j₃₂)
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
    Y₆ = r₂₁*r₃₁^2*ϕ₁₂^2*ϕ₂₁*(1-p)^2
    Y₅ = 2*r₂₁*r₃₁^2*u₂*ϕ₁₂^2*ϕ₂₁*(1-p)
    Y₄ = -2*p^2*r₂₁*r₃₁^2*ϕ₁₂*ϕ₁₃*ϕ₂₁+2*p^2*r₂₁*r₃₁^2*ϕ₁₂*ϕ₂₁-2*p^2*r₂₁*r₃₁*u₃*ϕ₁₂*ϕ₁₃*ϕ₂₁+2*p^2*r₂₁*r₃₁*u₄*ϕ₁₂*ϕ₁₃*ϕ₂₁+4*p*r₂₁*r₃₁^2*ϕ₁₂*ϕ₁₃*ϕ₂₁-4*p*r₂₁*r₃₁^2*ϕ₁₂*ϕ₂₁+4*p*r₂₁*r₃₁*u₃*ϕ₁₂*ϕ₁₃*ϕ₂₁-4*p*r₂₁*r₃₁*u₄*ϕ₁₂*ϕ₁₃*ϕ₂₁+r₂₁*r₃₁^2*u₂^2*ϕ₁₂^2*ϕ₂₁-2*r₂₁*r₃₁^2*ϕ₁₂*ϕ₁₃*ϕ₂₁+2*r₂₁*r₃₁^2*ϕ₁₂*ϕ₂₁-2*r₂₁*r₃₁*u₃*ϕ₁₂*ϕ₁₃*ϕ₂₁+2*r₂₁*r₃₁*u₄*ϕ₁₂*ϕ₁₃*ϕ₂₁
    Y₃ = -p^2*r₂₁*r₃₁^2+4*p*r₂₁*r₃₁^2*u₂*ϕ₁₂*ϕ₁₃*ϕ₂₁-4*p*r₂₁*r₃₁^2*u₂*ϕ₁₂*ϕ₂₁+2*p*r₂₁*r₃₁^2+2*p*r₂₁*r₃₁*u₂*u₃*ϕ₁₂*ϕ₁₃*ϕ₂₁-4*p*r₂₁*r₃₁*u₂*u₄*ϕ₁₂*ϕ₁₃*ϕ₂₁-4*r₂₁*r₃₁^2*u₂*ϕ₁₂*ϕ₁₃*ϕ₂₁+4*r₂₁*r₃₁^2*u₂*ϕ₁₂*ϕ₂₁-r₂₁*r₃₁^2-2*r₂₁*r₃₁*u₂*u₃*ϕ₁₂*ϕ₁₃*ϕ₂₁+4*r₂₁*r₃₁*u₂*u₄*ϕ₁₂*ϕ₁₃*ϕ₂₁
    Y₂ = p^2*r₂₁*r₃₁^2*ϕ₁₃^2*ϕ₂₁-2*p^2*r₂₁*r₃₁^2*ϕ₁₃*ϕ₂₁+p^2*r₂₁*r₃₁^2*ϕ₂₁+p^2*r₂₁*r₃₁^2+2*p^2*r₂₁*r₃₁*u₃*ϕ₁₃^2*ϕ₂₁-2*p^2*r₂₁*r₃₁*u₃*ϕ₁₃*ϕ₂₁-2*p^2*r₂₁*r₃₁*u₄*ϕ₁₃^2*ϕ₂₁+2*p^2*r₂₁*r₃₁*u₄*ϕ₁₃*ϕ₂₁+p^2*r₂₁*u₃^2*ϕ₁₃^2*ϕ₂₁-2*p^2*r₂₁*u₃*u₄*ϕ₁₃^2*ϕ₂₁+p^2*r₂₁*u₄^2*ϕ₁₃^2*ϕ₂₁+2*p*r₂₁*r₃₁^2*u₂-2*p*r₂₁*r₃₁^2*ϕ₁₃^2*ϕ₂₁+4*p*r₂₁*r₃₁^2*ϕ₁₃*ϕ₂₁-2*p*r₂₁*r₃₁^2*ϕ₂₁-2*p*r₂₁*r₃₁^2-4*p*r₂₁*r₃₁*u₃*ϕ₁₃^2*ϕ₂₁+4*p*r₂₁*r₃₁*u₃*ϕ₁₃*ϕ₂₁+4*p*r₂₁*r₃₁*u₄*ϕ₁₃^2*ϕ₂₁-4*p*r₂₁*r₃₁*u₄*ϕ₁₃*ϕ₂₁-2*p*r₂₁*u₃^2*ϕ₁₃^2*ϕ₂₁+4*p*r₂₁*u₃*u₄*ϕ₁₃^2*ϕ₂₁-2*p*r₂₁*u₄^2*ϕ₁₃^2*ϕ₂₁-2*r₂₁*r₃₁^2*u₂^2*ϕ₁₂*ϕ₁₃*ϕ₂₁+2*r₂₁*r₃₁^2*u₂^2*ϕ₁₂*ϕ₂₁-2*r₂₁*r₃₁^2*u₂+r₂₁*r₃₁^2*ϕ₁₃^2*ϕ₂₁-2*r₂₁*r₃₁^2*ϕ₁₃*ϕ₂₁+r₂₁*r₃₁^2*ϕ₂₁+r₂₁*r₃₁^2+2*r₂₁*r₃₁*u₂^2*u₄*ϕ₁₂*ϕ₁₃*ϕ₂₁+2*r₂₁*r₃₁*u₃*ϕ₁₃^2*ϕ₂₁-2*r₂₁*r₃₁*u₃*ϕ₁₃*ϕ₂₁-2*r₂₁*r₃₁*u₄*ϕ₁₃^2*ϕ₂₁+2*r₂₁*r₃₁*u₄*ϕ₁₃*ϕ₂₁+r₂₁*u₃^2*ϕ₁₃^2*ϕ₂₁-2*r₂₁*u₃*u₄*ϕ₁₃^2*ϕ₂₁+r₂₁*u₄^2*ϕ₁₃^2*ϕ₂₁
    Y₁ = -p^2*r₃₁^2*u₁-p^2*r₃₁*u₁*u₃+p^2*r₃₁*u₁*u₄-2*p*r₂₁*r₃₁^2*u₂*ϕ₁₃^2*ϕ₂₁+4*p*r₂₁*r₃₁^2*u₂*ϕ₁₃*ϕ₂₁-2*p*r₂₁*r₃₁^2*u₂*ϕ₂₁-2*p*r₂₁*r₃₁^2*u₂-2*p*r₂₁*r₃₁*u₂*u₃*ϕ₁₃^2*ϕ₂₁+2*p*r₂₁*r₃₁*u₂*u₃*ϕ₁₃*ϕ₂₁+4*p*r₂₁*r₃₁*u₂*u₄*ϕ₁₃^2*ϕ₂₁-4*p*r₂₁*r₃₁*u₂*u₄*ϕ₁₃*ϕ₂₁+2*p*r₂₁*u₂*u₃*u₄*ϕ₁₃^2*ϕ₂₁-2*p*r₂₁*u₂*u₄^2*ϕ₁₃^2*ϕ₂₁+2*p*r₃₁^2*u₁+2*p*r₃₁*u₁*u₃-2*p*r₃₁*u₁*u₄-r₂₁*r₃₁^2*u₂^2+2*r₂₁*r₃₁^2*u₂*ϕ₁₃^2*ϕ₂₁-4*r₂₁*r₃₁^2*u₂*ϕ₁₃*ϕ₂₁+2*r₂₁*r₃₁^2*u₂*ϕ₂₁+2*r₂₁*r₃₁^2*u₂+2*r₂₁*r₃₁*u₂*u₃*ϕ₁₃^2*ϕ₂₁-2*r₂₁*r₃₁*u₂*u₃*ϕ₁₃*ϕ₂₁-4*r₂₁*r₃₁*u₂*u₄*ϕ₁₃^2*ϕ₂₁+4*r₂₁*r₃₁*u₂*u₄*ϕ₁₃*ϕ₂₁-2*r₂₁*u₂*u₃*u₄*ϕ₁₃^2*ϕ₂₁+2*r₂₁*u₂*u₄^2*ϕ₁₃^2*ϕ₂₁-r₃₁^2*u₁-r₃₁*u₁*u₃+r₃₁*u₁*u₄
    Y₀ = p*r₃₁^2*u₁*u₂-p*r₃₁*u₁*u₂*u₄+r₂₁*r₃₁^2*u₂^2*ϕ₁₃^2*ϕ₂₁-2*r₂₁*r₃₁^2*u₂^2*ϕ₁₃*ϕ₂₁+r₂₁*r₃₁^2*u₂^2*ϕ₂₁+r₂₁*r₃₁^2*u₂^2-2*r₂₁*r₃₁*u₂^2*u₄*ϕ₁₃^2*ϕ₂₁+2*r₂₁*r₃₁*u₂^2*u₄*ϕ₁₃*ϕ₂₁+r₂₁*u₂^2*u₄^2*ϕ₁₃^2*ϕ₂₁-r₃₁^2*u₁*u₂+r₃₁*u₁*u₂*u₄
    if Y₀ < 0
        for root in roots(Polynomial([Y₀, Y₁, Y₂, Y₃, Y₄, Y₅, Y₆], :y))
            if (real(root) > 0) && (imag(root) == 0)
                y = real(root)
                z = 1+((((u₃*(1-p)*y)/(u₂+(1-p)*y))-u₄)/r₃₁)
                x = 1+ϕ₁₂*y^2-ϕ₁₃*z
                push!(sols_real, (x, y, z))
            end
        end
    end
    return sols_real
end;

function descarte(lst)
    sign_changes = 0
    for i in 1:(length(lst)-1)
        if lst[i]*lst[i+1] < 0
            sign_changes += 1
        end
    end
    return isodd(sign_changes)
end;

end
