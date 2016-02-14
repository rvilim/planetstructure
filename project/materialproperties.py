import numpy as np
import params
import phases


def get_K_S(r, rho, rho_0, T, K_0, Kp_0, alpha, gamma):
    if T is False:
        return get_K_T(rho, rho_0, T, K_0, Kp_0, gamma) * (1 + alpha * gamma * 300.0)

    return get_K_T(rho, rho_0, T, K_0, Kp_0, gamma) * (1 + alpha * gamma * T)


def get_K_T(rho, rho_0, T, K_0, Kp_0, gamma):
    if T is False:
        return get_K_T_rho300(rho, rho_0, K_0, Kp_0)
    else:
        return get_K_T_rho300(rho, rho_0, K_0, Kp_0) + get_DeltaK_th(rho, T, gamma)


def get_K_T_rho300(rho, rho_0, K_0, Kp_0):
    x = rho / rho_0
    # return 3*K_0*(pow(x,2.0/3)-pow(x,1.0/3))*math.exp(1.5*(Kp_0-1)*(1-pow(x,-1.0/3)))
    return (K_0 / 2.0) * (
        7 * pow(x, 7.0 / 3.0) - 5 * pow(x, 5.0 / 3.0) * (1 + .75 * (Kp_0 - 4) * (pow(x, 2.0 / 3) - 1)) + 1.5 * (
            pow(x, 3.0) - pow(x, 7.0 / 3)) * (Kp_0 - 4.0))


def get_DeltaK_th(rho, T, gamma):
    return 3 * n * R * gamma * rho * (get_f(T, gamma) - get_f(T_0, gamma))


def get_f_integrand(xi, gamma):
    print "Check get_f_integrand"
    return pow(xi, 3) / (math.exp(xi) - 1) + 3 * theta * gamma / (math.exp(theta / T) - 1)


def get_f(T, gamma):
    theta = theta0 * math.exp((gamma0 - gamma) / q)

    return (1 - q - 3 * gamma) * (pow(T, 4) / pow(theta, 3)) * integrate.quad(get_f_integrand, 0, theta / T)


def get_gamma(gamma_0, rho_0, rho, q):
    # return phase['gamma_0']*pow(phase['rho_0']/rho,phase['q'])
    return gamma_0 * pow(rho_0 / rho, q)


def get_eta(eta_0, T=1.0, T_0=1.0, n=1.0):
    return eta_0 * pow(T / T_0, -n)


def get_Ra(layer_solution, layer_number):
    eta = get_eta(10.0 ** 21)
    rho_avg = layer_solution['rho_avg']
    g_avg = layer_solution['g_avg']
    alpha_avg = layer_solution['alpha_avg']
    k_avg = layer_solution['k_avg']
    kappa_avg = layer_solution['kappa_avg']
    q = params.recipe['q_surf_est']

    return rho_avg * g_avg * alpha_avg * q * pow(layer_solution['top_r'] - layer_solution['bottom_r'], 4.0) / (
        kappa_avg * k_avg * eta)


def get_param_avg(layer_solution, layer_number, param):
    param_vals = np.zeros_like(layer_solution['r'])
    i = 0

    if ((params.recipe['layer_' + str(layer_number)]['type']).upper() == 'GEOM'):
        for (r, T, P) in zip(layer_solution['r'], layer_solution['T'], layer_solution['y'][:, 3]):
            param_vals[i] = phases.get_phase(layer_number, None, P, T)[param]
            i = i + 1

    return -np.trapz(param_vals * (layer_solution['r'] ** 2), x=layer_solution['r']) * 3.0 / (
        layer_solution['top_r'] ** 3 - layer_solution['bottom_r'] ** 3)


def get_Rac():
    return 1000.0


def get_boundary_layer(layer_solution):
    return params.a * ((layer_solution['top_r'] - layer_solution['bottom_r']) / 2.0) * pow(
        layer_solution['Ra_avg'] / layer_solution['Ra_c'], -0.25)
