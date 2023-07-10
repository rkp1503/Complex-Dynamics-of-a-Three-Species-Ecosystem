#=============================================================================
Author: Rayla Kurosaki
GitHub: https://github.com/rkp1503
=============================================================================#

module boundary_xz

function equilibrium_exist(param_vals)
    r₂₁, r₃₁, p, ϕ₁₂, ϕ₂₁, ϕ₁₃, u₁, u₂, u₃, u₄ = param_vals
    cond1 = (u₄/r₃₁)+(1/ϕ₁₃) > 1
    cond2 = r₃₁ > u₄
    return all([cond1, cond2])
end;

function equilibrium_stable(param_vals)
    r₂₁, r₃₁, p, ϕ₁₂, ϕ₂₁, ϕ₁₃, u₁, u₂, u₃, u₄ = param_vals
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