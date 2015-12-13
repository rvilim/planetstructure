import scipy.integrate as integrate
import numpy as np

import phases
import params
import math

def structure(recipe, R):

	M=recipe['mass']

	get_initial_profile(recipe, R)

def get_initial_profile(recipe, R):
	M=recipe['mass']

	g_surf=params.G*M/pow(R,2.0)
	P_surf=0.1 # This is a bit of a hack. Pressure at the surface should really be zero, update the phase relations
	T_surf=recipe['T_surf']
	rho_surf=phases.get_phase(recipe, 1, None, P_surf, T_surf)['rho_0']

	x=np.linspace(0, R, 100)

	print integrate.odeint(structure_equations, [rho_surf, g_surf, M, P_surf], x, args=(recipe, R, False))

def get_layer(recipe, m):

    layer=0
    layer_masses=recipe['layer_masses'][layer]
    m_above=recipe['mass']-m
    
    while layer<(recipe['layers']-1) and (m_above>layer_masses):
        layer=layer+1
        layer_masses+=recipe['layer_masses'][layer]

    # print "Layer Info", m,layer+1
    return layer+1


def structure_equations(y, x, recipe, R, Temperature):
	"""
	y=[ rho
	    g
	    m
	    P ]
	"""

	r=R-x

	rho=y[0]
	g=y[1]
	m=y[2]
	P=y[3]

	layer=get_layer(recipe, m)

	# print layer, m, recipe['mass'], recipe['mass']-m

	if (Temperature is False):
		phase=phases.get_phase(recipe, layer, rho, P, recipe['T_surf'])

		alpha=phase['alpha']
		gamma=phase['gamma_0']
		rho_0=phase['rho_0']
		K_0=phase['K_0']
		Kp_0=phase['Kp_0']

		phi=K_S(r, rho, rho_0, False, K_0, Kp_0, alpha, gamma)/rho
	else:
		print "Gimme time"
		# phi=K_S(r, rho, rho_0, T, K_0, Kp_0, alpha, gamma)/rho

	print str(x)+", "+phase['name']+", "+str(rho)+", "+str(g)+", "+str(m)+", "+str(P)

	drhodx=-(-rho*g/(phi*pow(10,9)))
	dgdx=-(4*np.pi*params.G*rho-2*params.G*m/pow(r,3))
	dmdx=-(4*np.pi*pow(r,2)*rho)
	dPdx=-(-rho*g)

	return [drhodx, dgdx, dmdx, dPdx]

def K_S(r, rho, rho_0, T, K_0, Kp_0, alpha, gamma):
	if T is False:
		return K_T(rho, rho_0, T, K_0, Kp_0)*(1+alpha*gamma*300.0)

	return K_T(rho, rho_0, T, K_0, Kp_0)*(1+alpha*gamma*T)

def K_T(rho, rho_0, T, K_0, Kp_0):
	if T is False:
		return K_T_rho300(rho, rho_0, K_0, Kp_0)
	else:
		return K_T_rho300(rho, rho_0, K_0, Kp_0)+DeltaK_th(rho, T)

def K_T_rho300(rho, rho_0, K_0, Kp_0):
	x=rho/rho_0
	# return 3*K_0*(pow(x,2.0/3)-pow(x,1.0/3))*math.exp(1.5*(Kp_0-1)*(1-pow(x,-1.0/3)))
	return (K_0/2.0)*(7*pow(x,7.0/3.0)-5*pow(x,5.0/3.0)*(1+.75*(Kp_0-4)*(pow(x,2.0/3)-1))+1.5*(pow(x,3.0)-pow(x,7.0/3))*(Kp_0-4.0))

def DeltaK_th(rho, T):
	return 3*n*R*gamma*rho*(f(T)-f(T_0))

def f_integrand(xi):
	return pow(xi,3)/(math.exp(xi)-1)+3*theta*gamma/(math.exp(theta/T)-1)

def f(T):
	theta=theta0*math.exp((gamma0-gamma)/q)

	return (1-q-3*gamma)*(pow(T,4)/pow(theta,3))*integrate.quad(f_integrand,0,theta/T)

