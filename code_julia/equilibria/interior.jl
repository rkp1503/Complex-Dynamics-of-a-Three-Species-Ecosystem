#=============================================================================
Author: Ramsey (Rayla) Phuc
Alias: Rayla Kurosaki
GitHub: https://github.com/rkp1503
Co-author: Ephraim Agyingi
=============================================================================#

module interior

import Polynomials: Polynomial, roots

function equilibrium_exist(param_vals)
    r₁, r₂, p, γ₁₂, γ₂₁, γ₁₃, γ₃₁, v₁, v₂, v₃ = param_vals
    E_x, E_y, E_z = analytical_solution(param_vals)
    cond1 = E_x > 0
    cond2 = E_y > 0
    cond3 = E_z > 0
    return all([cond1, cond2, cond3])
end;

function equilibrium_stable(param_vals, method)
    r₁, r₂, p, γ₁₂, γ₂₁, γ₁₃, γ₃₁, v₁, v₂, v₃ = param_vals
    E_x, E_y, E_z = analytical_solution(param_vals)
    j₁₁ = γ₁₂*E_y^2-γ₁₃*E_z-2*E_x+1
    j₁₂ = 2*γ₁₂*E_x*E_y
    j₁₃ = -γ₁₃*E_x
    j₂₁ = 2*γ₂₁*r₁*E_x*E_y
    j₂₂ = r₁*(1-2*E_y+γ₂₁*E_z^2)-((v₁*(1-p)*E_z)/(v₁+(1-p)*E_y)^2)
    j₂₃ = -((E_y*(1-p))/(v₁+E_y*(1-p)))
    j₃₂ = ((v₁*v₃*(1-p)*E_z)/(v₁+E_y*(1-p))^2)
    j₃₃ = r₂-v₂-2γ₃₁*r₂*E_z+((v₃*E_y*(1-p))/(v₁+E_y*(1-p)))
    if method == "RHC"
        C₂ = -j₁₁-j₂₂-j₃₃
        C₁ = j₁₁*j₂₂+j₁₁*j₃₃+j₂₂*j₃₃-j₁₂*j₂₁-j₂₃*j₃₂
        C₀ = (j₁₁*j₂₃*j₃₂+j₁₂*j₂₁*j₃₃)-(j₁₁*j₂₂*j₃₃+j₁₃*j₂₁*j₃₂)
        cond1 = C₂ > 0
        cond2 = C₁ > 0
        cond3 = C₀ > 0
        cond4 = C₂*C₁ > C₀
        return all([cond1, cond2, cond3, cond4])
    end
end;

function analytical_solution(param_vals)
    r₁, r₂, p, γ₁₂, γ₂₁, γ₁₃, γ₃₁, v₁, v₂, v₃ = param_vals
    E_x = -1
    E_y = -1
    E_z = -1
    Y₀ = v₁*(r₁*v₁*γ₁₃^2*γ₂₁*(r₂-v₂)^2+γ₃₁^2*r₁*r₂^2*v₁*(1+γ₂₁)+γ₃₁*r₂*(2*γ₁₃*γ₂₁*r₁*v₁+(1-p))*(v₂-r₂))
    Y₁ = 2*γ₁₃^2*γ₂₁*r₁*v₁*(r₂-v₂)*(r₂-v₂+v₃)*(1-p)+2*γ₁₃*γ₂₁*γ₃₁*r₁*r₂*v₁*(2*(v₂-r₂)-v₃)*(1-p)+γ₃₁^2*r₁*r₂^2*v₁*(2*(γ₂₁+1)*(1-p)-v₁)+γ₃₁*r₂*(v₂-r₂-v₃)*(1-p)^2
    Y₂ = r₁*(2*γ₁₂*γ₂₁*γ₁₃*γ₃₁*r₂*v₁^2*(v₂-r₂)+γ₁₃^2*γ₂₁*(r₂-v₂+v₃)^2*(1-p)^2+2*γ₂₁*γ₁₃*γ₃₁*r₂*(v₂-r₂-v₃)*(1-p)^2+γ₂₁*γ₃₁^2*r₂^2*(2*γ₁₂*v₁^2+(1-p)^2)+γ₃₁^2*r₂^2*((1-p)-2*v₁)*(1-p))
    Y₃ = γ₃₁*r₁*r₂*(2*γ₁₂*γ₁₃*γ₂₁*v₁*(2*(v₂-r₂)-v₃)+γ₃₁*r₂*(4*γ₁₂*γ₂₁*v₁-(1-p)))*(1-p)
    Y₄ = γ₁₂*γ₂₁*γ₃₁*r₁*r₂*(2*γ₁₃*(v₂-r₂-v₃)*(1-p)^2+γ₃₁*r₂*(γ₁₂*v₁^2+2*(1-p)^2))
    Y₅ = 2*γ₁₂^2*γ₂₁*γ₃₁^2*r₁*r₂^2*v₁*(1-p)
    Y₆ = γ₁₂^2*γ₂₁*γ₃₁^2*r₁*r₂^2*(1-p)^2
    coeff = [Y₀, Y₁, Y₂, Y₃, Y₄, Y₅, Y₆]
    if Y₀ < 0
        poly_roots = roots(Polynomial(coeff, :y))
        for root in poly_roots
            if (real(root) > 0) && (imag(root) == 0)
                E_y = real(root)
                E_z = (1/(r₂*γ₃₁))*(r₂-v₂+((v₃*(1-p)*E_y)/(v₁+(1-p)*E_y)))
                E_x = 1+γ₁₂*E_y^2-γ₁₃*E_z
                break
            end
        end
    end
    return (E_x, E_y, E_z)
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
