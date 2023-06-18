#=============================================================================
Author: Ramsey (Rayla) Phuc
Alias: Rayla Kurosaki
GitHub: https://github.com/rkp1503
Co-author: Ephraim Agyingi
=============================================================================#

module boundary_xz

function equilibrium_exist(param_vals)
    r₁, r₂, p, γ₁₂, γ₂₁, γ₁₃, γ₃₁, v₁, v₂, v₃ = param_vals
    cond1 = (r₂-v₂)/(r₂*γ₃₁) > 0
    cond2 = 1-((γ₁₃*(r₂-v₂))/(r₂*γ₃₁)) > 0
    return all([cond1, cond2])
end;

function equilibrium_stable(param_vals, method)
    r₁, r₂, p, γ₁₂, γ₂₁, γ₁₃, γ₃₁, v₁, v₂, v₃ = param_vals
    j₁₁ = ((γ₁₃*(r₂-v₂))/(r₂*γ₃₁))-1
    j₂₂ = r₁*(γ₂₁*j₁₁^2+1)-(((r₂-v₂)*(1-p))/(r₂*v₁*γ₃₁))
    j₃₃ = v₂-r₂
    if method == "RHC"
        C₂ = -j₁₁-j₂₂-j₃₃
        C₁ = j₁₁*j₂₂+j₁₁*j₃₃+j₂₂*j₃₃
        C₀ = -j₁₁*j₂₂*j₃₃
        cond1 = C₂ > 0
        cond2 = C₁ > 0
        cond3 = C₀ > 0
        cond4 = C₂*C₁ > C₀
        return all([cond1, cond2, cond3, cond4])
    else
        λ₁ = j₁₁
        λ₂ = j₂₂
        λ₃ = j₃₃
        cond1 = λ₁ < 0
        cond2 = λ₂ < 0
        cond3 = λ₃ < 0
        return all([cond1, cond2, cond3])
    end
end;

function analytical_solution(param_vals)
    r₁, r₂, p, γ₁₂, γ₂₁, γ₁₃, γ₃₁, v₁, v₂, v₃ = param_vals
    E_x = 1-((γ₁₃*(r₂-v₂))/(r₂*γ₃₁))
    E_y = 0
    E_z = (r₂-v₂)/(r₂*γ₃₁)
    return (E_x, E_y, E_z)
end;
    
end