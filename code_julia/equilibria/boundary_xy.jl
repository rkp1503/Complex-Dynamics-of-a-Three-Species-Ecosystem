#=============================================================================
Author: Ramsey (Rayla) Phuc
Alias: Rayla Kurosaki
GitHub: https://github.com/rkp1503
Co-author: Ephraim Agyingi
=============================================================================#

module boundary_xy

import Polynomials: Polynomial, roots

function equilibrium_exist(param_vals)
    r₁, r₂, p, γ₁₂, γ₂₁, γ₁₃, γ₃₁, v₁, v₂, v₃ = param_vals
    β = 9
    # cond = γ₂₁ < (β-1)/((γ₁₂*β^2+1)^2)
    # cond = γ₁₂ < (1/β^2)*(sqrt((β-1)/γ₂₁)-1)
    cond = γ₁₂^2*γ₂₁*β^4+2*γ₁₂*γ₂₁*β^2-β+γ₂₁+1 < 0
    return all([cond])
end;

function equilibrium_stable(param_vals, method)
    r₁, r₂, p, γ₁₂, γ₂₁, γ₁₃, γ₃₁, v₁, v₂, v₃ = param_vals
    E_x, E_y, E_z = analytical_solution(param_vals)
    j₁₁ = 1-2*E_x+γ₁₂*E_y^2
    j₁₂ = 2*γ₁₂*E_x*E_y
    j₂₁ = 2*r₁*γ₂₁*E_x*E_y
    j₂₂ = r₁*(1-2*E_y+γ₂₁*E_x^2)
    j₃₃ = r₂-v₂+((v₃*(1-p)*E_y)/(v₁+(1-p)*E_y))
    if method == "RHC"
        C₂ = -j₁₁-j₂₂-j₃₃
        C₁ = j₁₁*j₂₂+j₁₁*j₃₃+j₂₂*j₃₃-j₁₂*j₂₁
        C₀ = j₃₃*(j₁₂*j₂₁-j₁₁*j₂₂)
        cond1 = C₂ > 0
        cond2 = C₁ > 0
        cond3 = C₀ > 0
        cond4 = C₂*C₁ > C₀
        return all([cond1, cond2, cond3, cond4])
    else
        λ₁ = j₃₃
        λ₂ = (j₁₁+j₂₂+sqrt((j₁₁-j₂₂)^2+4*j₁₂*j₂₁))/2
        λ₃ = (j₁₁+j₂₂-sqrt((j₁₁-j₂₂)^2+4*j₁₂*j₂₁))/2
        cond1 = λ₁ < 0
        cond2 = λ₂ < 0
        cond3 = λ₃ < 0
        return all([cond1, cond2, cond3])
    end
end;

function analytical_solution(param_vals)
    r₁, r₂, p, γ₁₂, γ₂₁, γ₁₃, γ₃₁, v₁, v₂, v₃ = param_vals
    E_x = -1
    E_y = -1
    E_z = 0
    Y₀ = γ₂₁+1
    Y₁ = -1
    Y₂ = 2*γ₁₂*γ₂₁
    Y₃ = 0
    Y₄ = γ₁₂^2*γ₂₁
    poly_roots = roots(Polynomial([Y₀, Y₁, Y₂, Y₃, Y₄], :y))
    for root in poly_roots
        if (real(root) > 0) && (imag(root) == 0)
            E_y = real(root)
            E_x = 1+γ₁₂*E_y^2
            break
        end
    end
    return (E_x, E_y, E_z)
end;

end