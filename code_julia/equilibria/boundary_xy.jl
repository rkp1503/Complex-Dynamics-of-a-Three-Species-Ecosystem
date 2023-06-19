#=============================================================================
Author: Ramsey (Rayla) Phuc
Alias: Rayla Kurosaki
GitHub: https://github.com/rkp1503
Co-author: Ephraim Agyingi
=============================================================================#

module boundary_xy

import Polynomials: Polynomial, roots

function equilibrium_exist(param_vals)
    r₂₁, r₃₁, p, ϕ₁₂, ϕ₂₁, ϕ₁₃, u₁, u₂, u₃, u₄ = param_vals
    β = 11
    # cond = ϕ₂₁ < (β-1)/((ϕ₁₂*β^2+1)^2)
    # cond = ϕ₁₂ < (1/β^2)*(sqrt((β-1)/ϕ₂₁)-1)
    cond = ϕ₁₂^2*ϕ₂₁*β^4+2*ϕ₁₂*ϕ₂₁*β^2-β+ϕ₂₁+1 < 0
    return all([cond])
end;

function equilibrium_stable(param_vals)
    r₂₁, r₃₁, p, ϕ₁₂, ϕ₂₁, ϕ₁₃, u₁, u₂, u₃, u₄ = param_vals
    x, y, z = analytical_solution(param_vals)
    j₁₁ = 1-2*x+ϕ₁₂*y^2
    j₁₂ = 2*ϕ₁₂*x*y
    j₂₁ = 2*r₂₁*ϕ₂₁*x*y
    j₂₂ = r₂₁*(1-2*y+ϕ₂₁*x^2)
    j₃₃ = r₃₁-u₄+((u₃*(1-p)*y)/(u₂+(1-p)*y))
    # if method == "RHC"
    #     C₂ = -j₁₁-j₂₂-j₃₃
    #     C₁ = j₁₁*j₂₂+j₁₁*j₃₃+j₂₂*j₃₃-j₁₂*j₂₁
    #     C₀ = j₃₃*(j₁₂*j₂₁-j₁₁*j₂₂)
    #     cond1 = C₂ > 0
    #     cond2 = C₁ > 0
    #     cond3 = C₀ > 0
    #     cond4 = C₂*C₁ > C₀
    #     return all([cond1, cond2, cond3, cond4])
    # else
    #     λ₁ = j₃₃
    #     λ₂ = (j₁₁+j₂₂+sqrt((j₁₁-j₂₂)^2+4*j₁₂*j₂₁))/2
    #     λ₃ = (j₁₁+j₂₂-sqrt((j₁₁-j₂₂)^2+4*j₁₂*j₂₁))/2
    #     cond1 = λ₁ < 0
    #     cond2 = λ₂ < 0
    #     cond3 = λ₃ < 0
    #     return all([cond1, cond2, cond3])
    # end
    C₂ = -j₁₁-j₂₂-j₃₃
    C₁ = j₁₁*j₂₂+j₁₁*j₃₃+j₂₂*j₃₃-j₁₂*j₂₁
    C₀ = j₃₃*(j₁₂*j₂₁-j₁₁*j₂₂)
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
    Y₀ = ϕ₂₁+1
    Y₁ = -1
    Y₂ = 2*ϕ₁₂*ϕ₂₁
    Y₃ = 0
    Y₄ = ϕ₁₂^2*ϕ₂₁
    poly_roots = roots(Polynomial([Y₀, Y₁, Y₂, Y₃, Y₄], :y))
    for root in poly_roots
        if (real(root) > 0) && (imag(root) == 0)
            y = real(root)
            x = 1+ϕ₁₂*y^2
            break
        end
    end
    return (x, y, 0)
end;

end