#=============================================================================
Author: Ramsey (Rayla) Phuc
Alias: Rayla Kurosaki
GitHub: https://github.com/rkp1503
Co-author: Ephraim Agyingi
=============================================================================#

module boundary_xz

function equilibrium_exist(param_vals)
    r₂₁, r₃₁, p, ϕ₁₂, ϕ₂₁, ϕ₁₃, u₁, u₂, u₃, u₄ = param_vals
    cond1 = 1-ϕ₁₃*(1-(u₄/r₃₁)) > 0
    cond2 = 1-(u₄/r₃₁) > 0
    return all([cond1, cond2])
end;

function equilibrium_stable(param_vals)
    r₂₁, r₃₁, p, ϕ₁₂, ϕ₂₁, ϕ₁₃, u₁, u₂, u₃, u₄ = param_vals
    # j₁₁ = 1-2*x-ϕ₁₃*(1-(u₄/r₃₁))
    # j₂₂ = r₂₁*(1+ϕ₂₁*(1-ϕ₁₃*(1-(u₄/r₃₁)))^2)-((u₁*(1-p)*(1-(u₄/r₃₁)))/u₂)
    # j₃₃ = r₃₁*(1-2*(1-(u₄/r₃₁)))-u₄
    # if method == "RHC"
    #     C₂ = -j₁₁-j₂₂-j₃₃
    #     C₁ = j₁₁*j₂₂+j₁₁*j₃₃+j₂₂*j₃₃
    #     C₀ = -j₁₁*j₂₂*j₃₃
    #     cond1 = C₂ > 0
    #     cond2 = C₁ > 0
    #     cond3 = C₀ > 0
    #     cond4 = C₂*C₁ > C₀
    #     return all([cond1, cond2, cond3, cond4])
    # else
    #     λ₁ = j₁₁
    #     λ₂ = j₂₂
    #     λ₃ = j₃₃
    #     cond1 = λ₁ < 0
    #     cond2 = λ₂ < 0
    #     cond3 = λ₃ < 0
    #     return all([cond1, cond2, cond3])
    # end
    cond1 = 1-(u₄/r₃₁) > (1-2*(1-ϕ₁₃*(1-(u₄/r₃₁))))/ϕ₁₃
    cond2 = 1-(u₄/r₃₁) > r₂₁*u₂*(1+ϕ₂₁*(1-ϕ₁₃*(1-(u₄/r₃₁)))^2)/(u₁*(1-p))
    cond3 = 1-(u₄/r₃₁) > (r₃₁-u₄)/(2*r₃₁)
    return all([cond1, cond2, cond3])
end;

function analytical_solution(param_vals)
    r₂₁, r₃₁, p, ϕ₁₂, ϕ₂₁, ϕ₁₃, u₁, u₂, u₃, u₄ = param_vals
    x = 1-ϕ₁₃*(1-(u₄/r₃₁))
    z = 1-(u₄/r₃₁)
    return (x, 0, z)
end;
    
end