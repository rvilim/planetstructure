from scipy.interpolate import UnivariateSpline, interp1d, splrep
import numpy as np


def interp_solution(solution):
    # This function takes in a solution profile and interpolates it for use in the T ode

    layer_num = 0

    for layer in solution:
        # print layer
        # rho=UnivariateSpline(layer['r'], layer['y'][:,0])
        # g=UnivariateSpline(layer['r'], layer['y'][:,1])
        # m=UnivariateSpline(layer['r'], layer['y'][:,2])
        # P=UnivariateSpline(layer['r'], layer['y'][:,3])
        # T=UnivariateSpline(layer['r'], layer['T'])

        rho = interp1d(layer['r'], layer['y'][:, 0], kind='cubic')
        g = interp1d(layer['r'], layer['y'][:, 1], kind='cubic')
        m = interp1d(layer['r'], layer['y'][:, 2], kind='cubic')
        P = interp1d(layer['r'], layer['y'][:, 3], kind='cubic')
        T = interp1d(layer['r'], layer['T'], kind='cubic')

        solution[layer_num]['y_interp'] = {'rho': rho, 'g': g, 'm': m, 'P': P, 'T': T}

        layer_num = layer_num + 1

    return solution
