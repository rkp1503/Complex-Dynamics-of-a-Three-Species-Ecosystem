#=============================================================================
Author: Ramsey (Rayla) Phuc
Alias: Rayla Kurosaki
GitHub: https://github.com/rkp1503
Co-author: Ephraim Agyingi
=============================================================================#

module axial_z

function equilibrium_exist(param_vals)
    r₂₁, r₃₁, p, ϕ₁₂, ϕ₂₁, ϕ₁₃, u₁, u₂, u₃, u₄ = param_vals
    cond = r₃₁ > u₄
    return all([cond])
end;

function equilibrium_stable(param_vals)
    r₂₁, r₃₁, p, ϕ₁₂, ϕ₂₁, ϕ₁₃, u₁, u₂, u₃, u₄ = param_vals
    cond1 = (u₄/r₃₁)+(1/ϕ₁₃) < 1
    cond2 = (u₄/r₃₁)+((r₂₁*u₂)/(u₁*(1-p))) < 1
    cond3 = (u₄/r₃₁)<(1/2)
    return all([cond1, cond2, cond3])
end;

function analytical_solution(param_vals)
    r₂₁, r₃₁, p, ϕ₁₂, ϕ₂₁, ϕ₁₃, u₁, u₂, u₃, u₄ = param_vals
    z = 1-(u₄/r₃₁)
    return (0, 0, z)
end;
    
end