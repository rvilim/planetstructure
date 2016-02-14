from scipy.integrate import odeint

import params


def temperature(structure):
    q_top_bl = solution[layer_number]['q_surf']
    T0 = params.recipe['T_surf']

    for layer_number in xrange(0, params.recipe['layers']):
        x = np.linspace(0.0, structure[layer_number]['top_r'] - structure[layer_number]['bottom_r'], params.layerpoints)

        layer_solution = odeint(temperature_equations, [T0], x, args=(phase))

        dTdr_top_bl =
        x_upper = np.array([0.0, layer['top_bl']])
        x_lower = np.array(
            [layer['top_r'] - layer['bottom_r'] - layer['bottom_bl'], layer['top_r'] - layer['bottom_r']])


def temperature_equations(y, x, layer):
    """
    y=[ T ]
    """

    T = y[0]

    r = layer['top_r'] - x
    phase = layer['phase']

    rho = layer['y_interp']['rho'](r)
    g = layer['y_interp']['g'](r)

    gamma = get_gamma(phase['gamma_0'], phase['rho_0'], rho, phase['q'])
    K_S = get_K_S(r, rho, phase['rho_0'], T, phase['K_0'], phase['Kp_0'], phase['alpha'], gamma)

    if (x < layer['top_bl']):
        # Then we are in the upper boundary layer
        dTdr = -layer['top_q'] / k
    elif (x > layer['top_r'] - layer['bottom_bl']):
        # Then we are in the lower boundary layer
        dTdr = -layer[ / k
        else:
        dTdr = -rho * g * T * gamma / K_S
