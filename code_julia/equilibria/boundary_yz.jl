#=============================================================================
Author: Ramsey (Rayla) Phuc
Alias: Rayla Kurosaki
GitHub: https://github.com/rkp1503
Co-author: Ephraim Agyingi
=============================================================================#

module boundary_yz

import Polynomials: Polynomial, roots

function equilibrium_exist(param_vals)
    r₁, r₂, p, γ₁₂, γ₂₁, γ₁₃, γ₃₁, v₁, v₂, v₃ = param_vals
    E_x, E_y, E_z = analytical_solution(param_vals)
    cond1 = r₂-v₂ > 0
    cond2 = E_y > 0
    return all([cond1, cond2])
end;

function equilibrium_stable(param_vals, method)
    r₁, r₂, p, γ₁₂, γ₂₁, γ₁₃, γ₃₁, v₁, v₂, v₃ = param_vals
    E_x, E_y, E_z = analytical_solution(param_vals)
    j₁₁ = 1+γ₁₂*E_y-γ₁₃*E_z
    j₂₂ = r₁*(1-2*E_y)-((v₁*(1-p)*E_z)/(v₁+(1-p)*E_y)^2)
    j₂₃ = -(((1-p)*E_y)/(v₁+(1-p)*E_y))
    j₃₂ = ((v₁*v₃*(1-p)*E_z)/(v₁+(1-p)*E_y)^2)
    j₃₃ = r₂-v₂-2*r₂*γ₃₁*E_z+((v₃*(1-p)*E_y)/(v₁+(1-p)*E_y))
    if method == "RHC"
        C₂ = -j₁₁-j₂₂-j₃₃
        C₁ = j₁₁*j₂₂+j₁₁*j₃₃+j₂₂*j₃₃-j₂₃*j₃₂
        C₀ = j₁₁*(j₂₃*j₃₂-j₂₂*j₃₃)
        cond1 = C₂ > 0
        cond2 = C₁ > 0
        cond3 = C₀ > 0
        cond4 = C₂*C₁ > C₀
        return all([cond1, cond2, cond3, cond4])
    else
        λ₁ = j₁₁
        λ₂ = (j₂₂+j₃₃+sqrt((j₂₂-j₃₃)^2+4*j₂₃*j₃₂))/2
        λ₃ = (j₂₂+j₃₃-sqrt((j₂₂-j₃₃)^2+4*j₂₃*j₃₂))/2
        cond1 = λ₁ < 0
        cond2 = λ₂ < 0
        cond3 = λ₃ < 0
        return all([cond1, cond2, cond3])
    end
end;

function analytical_solution(param_vals)
    r₁, r₂, p, γ₁₂, γ₂₁, γ₁₃, γ₃₁, v₁, v₂, v₃ = param_vals
    E_x = 0
    E_y = -1
    E_z = -1
    Y₀ = v₁*(r₁*r₂*v₁*γ₃₁+(v₂-r₂)*(1-p))
    Y₁ = r₁*r₂*v₁*γ₃₁*(2*(1-p)-v₁)+(v₂-v₃-r₂)*(1-p)^2
    Y₂ = r₁*r₂*γ₃₁*((1-p)-2*v₁)*(1-p)
    Y₃ = -r₁*r₂*γ₃₁*(1-p)^2
    coeff = [Y₀, Y₁, Y₂, Y₃]
    if Y₀ > 0
        poly_roots = roots(Polynomial(coeff, :y))
        for root in poly_roots
            if (real(root) > 0) && (imag(root) == 0)
                E_y = real(root)
                E_z = (1/(r₂*γ₃₁))*(r₂-v₂+((v₃*(1-p)*E_y)/(v₁+(1-p)*E_y)))
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
