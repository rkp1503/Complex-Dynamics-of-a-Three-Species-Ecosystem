#=============================================================================
Author: Ramsey (Rayla) Phuc
Alias: Rayla Kurosaki
GitHub: https://github.com/rkp1503
Co-author: Ephraim Agyingi
=============================================================================#

module axial_z

function equilibrium_exist(param_vals)
    r₂₁, r₃₁, p, ϕ₁₂, ϕ₂₁, ϕ₁₃, u₁, u₂, u₃, u₄ = param_vals
    cond = 1-(u₄/r₃₁) > 0
    return all([cond])
end;

function equilibrium_stable(param_vals)
    r₂₁, r₃₁, p, ϕ₁₂, ϕ₂₁, ϕ₁₃, u₁, u₂, u₃, u₄ = param_vals
    # j₁₁ = 1-(ϕ₁₃*(1-(u₄/r₃₁)))
    # j₂₂ = r₂₁-(u₁*(1-p)/u₂)*(1-(u₄/r₃₁))
    # j₃₃ = r₃₁*(1-2*(1-(u₄/r₃₁)))
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
    #     cond1 = (u₄/r₃₁)+ϕ₁₃ < 1
    #     cond2 = (u₄/r₃₁)+((r₂₁*u₂)/(u₁*(1-p))) < 1
    #     cond2 = u₄/r₃₁ < 1/2
    #     return all([cond1, cond2, cond3])
    # end
    cond1 = 1-(ϕ₁₃*(1-(u₄/r₃₁))) < 0
    cond2 = r₂₁-(u₁*(1-p)/u₂)*(1-(u₄/r₃₁)) < 0
    cond3 = r₃₁*(1-2*(1-(u₄/r₃₁))) < 0
    return all([cond1, cond2, cond3])
end;

function analytical_solution(param_vals)
    r₂₁, r₃₁, p, ϕ₁₂, ϕ₂₁, ϕ₁₃, u₁, u₂, u₃, u₄ = param_vals
    z = 1-(u₄/r₃₁)
    return (0, 0, z)
end;
    
end