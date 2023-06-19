#=============================================================================
Author: Ramsey (Rayla) Phuc
Alias: Rayla Kurosaki
GitHub: https://github.com/rkp1503
Co-author: Ephraim Agyingi
=============================================================================#

module interior

import Polynomials: Polynomial, roots

function equilibrium_exist(param_vals)
    r₂₁, r₃₁, p, ϕ₁₂, ϕ₂₁, ϕ₁₃, u₁, u₂, u₃, u₄ = param_vals
    x, y, z = analytical_solution(param_vals)
    cond1 = (1+ϕ₁₂*y^2)/ϕ₁₃ > z
    cond2 = y > (u₂*(u₄-r₃₁))/((u₃-u₄+r₃₁)*(1-p))
    cond3 = r₃₁*u₁*(u₄-r₃₁)*(1-p)+r₂₁*u₂*(ϕ₁₃^2*ϕ₂₁*(r₃₁-u₄)^2+2*r₃₁*ϕ₁₃*ϕ₂₁*(u₄-r₃₁)+r₃₁*(1+ϕ₂₁)) < 0
    return all([cond1, cond2, cond3])
end;

function equilibrium_stable(param_vals)
    r₂₁, r₃₁, p, ϕ₁₂, ϕ₂₁, ϕ₁₃, u₁, u₂, u₃, u₄ = param_vals
    x, y, z = analytical_solution(param_vals)
    j₁₁ = 1-2*x+ϕ₁₂*y^2-ϕ₁₃*z
    j₁₂ = 2*ϕ₁₂*x*y
    j₁₃ = -ϕ₁₃*x
    j₂₁ = 2*ϕ₂₁*r₂₁*x*y
    j₂₂ = r₂₁*(1-2*y+ϕ₂₁*x^2)-((u₁*u₂*(1-p)*z)/(u₂+(1-p)*y)^2)
    j₂₃ = -((u₁*(1-p)*y)/(u₂+(1-p)*y))
    j₃₂ = ((u₂*u₃*(1-p)*z)/(u₂+(1-p)*y)^2)
    j₃₃ = r₃₁*(1-2*z)-u₄+((u₃*(1-p)*y)/(u₂+(1-p)*y))

    C₂ = -j₁₁-j₂₂-j₃₃
    C₁ = j₁₁*j₂₂+j₁₁*j₃₃+j₂₂*j₃₃-j₁₂*j₂₁-j₂₃*j₃₂
    C₀ = (j₁₁*j₂₃*j₃₂+j₁₂*j₂₁*j₃₃)-(j₁₁*j₂₂*j₃₃+j₁₃*j₂₁*j₃₂)
    cond1 = C₂ > 0
    cond2 = C₁ > 0
    cond3 = C₀ > 0
    cond4 = C₂*C₁ > C₀
    return all([cond1, cond2, cond3, cond4])
end;

function analytical_solution(param_vals)
    r₂₁, r₃₁, p, ϕ₁₂, ϕ₂₁, ϕ₁₃, u₁, u₂, u₃, u₄ = param_vals
    x = -1
    y = -1
    z = -1
    Y₀ = u₂^2*(r₃₁*u₁*(u₄-r₃₁)*(1-p)+r₂₁*u₂*(ϕ₁₃^2*ϕ₂₁*(r₃₁-u₄)^2+2*r₃₁*ϕ₁₃*ϕ₂₁*(u₄-r₃₁)+r₃₁^2*(1+ϕ₂₁)))
    Y₁ = -u₂*(r₃₁*u₁*(2*(r₃₁-u₄)+u₃)*(1-p)^2+r₂₁*u₂*(ϕ₁₃^2*ϕ₂₁*(3*(r₃₁-u₄)+2*u₃)*(u₄-r₃₁)+2*r₃₁*ϕ₁₃*ϕ₂₁*(3*(r₃₁-u₄)+u₃)-3*r₃₁^2*(1+ϕ₂₁))*(1-p)+r₂₁*r₃₁^2*u₂^2)
    Y₂ = r₃₁*u₁*(u₄-r₃₁-u₃)*(1-p)^3+r₂₁*u₂*(ϕ₁₃^2*ϕ₂₁*(3*(r₃₁-u₄)+u₃)*(r₃₁+u₃-u₄)+2*r₃₁*ϕ₁₃*ϕ₂₁*(3*(r₃₁-u₄)+2*u₃)+3*r₃₁^2*(1+ϕ₂₁))*(1-p)^2-3*r₂₁*r₃₁^2*u₂^2*(1-p)+2*r₂₁*r₃₁*u₂^3*ϕ₁₂*ϕ₂₁*(r₃₁+ϕ₁₃*(u₄-r₃₁))
    Y₃ = r₂₁*(1-p)*(ϕ₂₁*(ϕ₁₃*(r₃₁+u₃-u₄)-r₃₁)^2*(1-p)^2-3*r₃₁^2*u₂*(1-p)-2*ϕ₁₂*ϕ₂₁*r₃₁*u₂^2*(ϕ₁₃*(3*(r₃₁-u₄)+u₃)-r₃₁))
    Y₄ = r₂₁*r₃₁*(ϕ₁₂*ϕ₂₁*u₂*(2*(3*r₃₁*(1-ϕ₁₃)+ϕ₁₃*(3*u₄-2*u₃))*(1-p)^2+r₃₁*ϕ₁₂*u₂^2)-r₃₁*(1-p)^3)
    Y₅ = -r₂₁*r₃₁*ϕ₁₂*ϕ₂₁*(2*(ϕ₁₃*(r₃₁+u₃-u₄)-r₃₁)*(1-p)^2-3*ϕ₁₂*r₃₁*u₂^2)*(1-p)
    Y₆ = 3*r₂₁*r₃₁^2*ϕ₁₂^2*ϕ₂₁*u₂*(1-p)^2
    Y₇ = r₂₁*r₃₁^2*ϕ₁₂^2*ϕ₂₁*(1-p)^3
    coeff = [Y₀, Y₁, Y₂, Y₃, Y₄, Y₅, Y₆, Y₇]
    if Y₀ < 0
        poly_roots = roots(Polynomial(coeff, :y))
        for root in poly_roots
            if (real(root) > 0) && (imag(root) == 0)
                y = real(root)
                z = 1+(((u₃*(1-p)*y)/(u₂+(1-p)*y))-u₄)/r₃₁
                x = 1+ϕ₁₂*y^2-ϕ₁₃*z
                break
            end
        end
    end
    return (x, y, z)
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
