"""
Author: Ramsey (Rayla) Phuc

Alias: Rayla Kurosaki

GitHub: https://github.com/rkp1503

Co-author: Ephraim Agyingi
"""

import sympy as sym


def xy_boundary(vars_dict: dict[sym.core.Symbol, float], params_dict: dict[sym.core.Symbol, float], eq_vals_temp: list[int | sym.core.Symbol], model_type="Proposed"):
    x, y, z = list(vars_dict.keys())
    params_symb = list(params_dict.keys())
    if model_type == "Original":
        r, p, gamma12, gamma21, gamma, v1, v2, v3 = params_symb[0] if len(params_symb) == 1 else params_symb
        eq_vals_temp[1] = 1-gamma21*x**2
        pass
    else:
        r_yx, r_zx, p, phi_xy, phi_yx, phi_xz, u_1, u_2, u_3, u_4 = params_symb[0] if len(params_symb) == 1 else params_symb
        eq_vals_temp[0] = 1+phi_xy*y**2
        pass
    return dict(zip([x, y, z], [x, y, 0]))
    # return dict(zip([x, y, z], eq_vals_temp))


def yz_boundary(vars_dict, params_dict, eq_vals_temp):
    x, y, z = list(vars_dict.keys())
    params_symb = list(params_dict.keys())
    r_yx, r_zx, p, phi_xy, phi_yx, phi_xz, u_1, u_2, u_3, u_4 = params_symb[0] if len(params_symb) == 1 else params_symb
    eq_vals_temp[2] = 1+((1/r_zx)*(((u_3*(1-p)*y)/(u_2+(1-p)*y))-u_4))
    return dict(zip([x, y, z], [0, y, z]))
    # return dict(zip([x, y, z], eq_vals_temp))


def interior(vars_dict, params_dict, eq_vals_temp, model_type="Proposed"):
    x, y, z = list(vars_dict.keys())
    params_symb = list(params_dict.keys())
    if model_type == "Original":
        r, p, gamma12, gamma21, gamma, v1, v2, v3 = params_symb[0] if len(params_symb) == 1 else params_symb
        eq_vals_temp[1] = v1*v2/((1-p)*(v3-v2))
        eq_vals_temp[2] = ((x-1)*(v3-v2)**2+gamma12*(v1*v2)**2)/(gamma*((1-p)*(v3-v2))**2)
        pass
    else:
        r_yx, r_zx, p, phi_xy, phi_yx, phi_xz, u_1, u_2, u_3, u_4 = params_symb[0] if len(params_symb) == 1 else params_symb
        eq_vals_temp[2] = 1+((1/r_zx)*(((u_3*(1-p)*y)/(u_2+(1-p)*y))-u_4))
        eq_vals_temp[0] = 1+phi_xy*y**2-phi_xz*eq_vals_temp[2]
        pass
    return dict(zip([x, y, z], [x, y, z]))
    # return dict(zip([x, y, z], eq_vals_temp))
