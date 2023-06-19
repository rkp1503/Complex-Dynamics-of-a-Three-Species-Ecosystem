#=============================================================================
Author: Ramsey (Rayla) Phuc
Alias: Rayla Kurosaki
GitHub: https://github.com/rkp1503
Co-author: Ephraim Agyingi
=============================================================================#

module boundary_yz

import Polynomials: Polynomial, roots

function equilibrium_exist(param_vals)
    r₂₁, r₃₁, p, ϕ₁₂, ϕ₂₁, ϕ₁₃, u₁, u₂, u₃, u₄ = param_vals
    x, y, z = analytical_solution(param_vals)
    cond1 = y > (u₂*(u₄-r₃₁))/((u₃-u₄+r₃₁)*(1-p))
    cond2 = ((r₁*u₂)/(u₁*(1-p)))+(u₄/r₃₁) < 1
    return all([cond1, cond2])
end;

function equilibrium_stable(param_vals)
    r₁, r₃₁, p, ϕ₁₂, ϕ₂₁, ϕ₁₃, u₁, u₂, u₃, u₄ = param_vals
    x, y, z = analytical_solution(param_vals)
    j₁₁ = 1+ϕ₁₂*y^2-ϕ₁₃*z
    j₂₂ = r₁*(1-2*y)-((u₁*u₂*(1-p)*z)/(u₂+(1-p)*y)^2)
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
    return all([cond1, cond2, cond3, cond4])
end;

function analytical_solution(param_vals)
    r₁, r₃₁, p, ϕ₁₂, ϕ₂₁, ϕ₁₃, u₁, u₂, u₃, u₄ = param_vals
    y = -1
    z = -1
    Y₀ = u₂*(r₁*r₃₁*u₂+u₁*(u₄-r₃₁)*(1-p))
    Y₁ = r₁*r₃₁*u₂*(2*(1-p)-u₂)+u₁*(u₃+u₄-r₃₁)*(1-p)^2
    Y₂ = r₁*r₃₁*((1-p)-2*u₂)*(1-p)
    Y₃ = -r₁*r₃₁*(1-p)^2
    coeff = [Y₀, Y₁, Y₂, Y₃]
    if Y₀ < 0
        poly_roots = roots(Polynomial(coeff, :y))
        for root in poly_roots
            if (real(root) > 0) && (imag(root) == 0)
                y = real(root)
                z = 1+(((u₃*(1-p)*y)/(u₂+(1-p)*y))-u₄)/r₃₁
                break
            end
        end
    end
    return (0, y, z)
end;

function descarte(lst)
    unit_lst = [n/abs(n) for n in lst]
    sign_changes = 0
    for i in 1:(length(unit_lst)-1)
        if unit_lst[i]*unit_lst[i+1] < 0
            sign_changes += 1
        end
    end
    return isodd(sign_changes)
end;

end
