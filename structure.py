# import scipy.integrate as integrate
from pycse import odelay

import numpy as np

import phases
import params
import math
import matplotlib.pyplot as plt
from interp import interp_solution

def structure(recipe, R):

	M=recipe['mass']

	solution=get_profile(recipe, R)

	# for i in xrange(0,len(solution)):

	# 	for a in xrange(0,len(solution[i]['y'])):
	# 		print solution[i]['y'][a][0]


	# print solution[0]['y']
	# print solution[1]['y']

	r=[]
	rho=[]
	g=[]
	m=[]
	P=[]

	a=[]


	for layer in solution:
		r=r+[a/1000.0 for a in layer['r']]
		rho=rho+[y[0] for y in layer['y']]
		g=g+[y[1] for y in layer['y']]
		m=m+[y[2] for y in layer['y']]
		P=P+[y[3] for y in layer['y']]

	f, ((ax1, ax2),(ax3, ax4)) = plt.subplots(2, 2)
	ax1.plot(r,rho)
	ax1.set_title('Rho')
	ax2.plot(r,g)
	ax2.set_title('g')
	ax3.plot(r,m)
	ax3.set_title('m (below)')
	ax4.plot(r,P)
	ax4.set_title('P')

	plt.show()

def get_profile(recipe, R):
	M=recipe['mass']

	g=params.G*M/pow(R,2.0)
	P=recipe['P_surf'] 
	T=None
	phase=phases.get_phase(recipe, 0, None, P, recipe['T_surf'])
	rho=phase['rho_0']

	layersbelow_mass=recipe['mass']
	radius=R

	solution=[None]*2

	for layer_number in xrange(0,recipe['layers']):
		x=np.linspace(0, radius, 100)

		layersbelow_mass=layersbelow_mass-recipe['layer_masses'][layer_number]
		define_mass_event(layersbelow_mass)

		layer_soln=odelay(structure_equations, [rho,g,M,P], x, events=[mass_event], args=(recipe, phase, radius, False))

		layer_end_x=layer_soln[2][0]
		layer_end_conditions=layer_soln[3][0]

		solution[layer_number]={}
		solution[layer_number]['X']=layer_soln[0]
		solution[layer_number]['y']=layer_soln[1]
		solution[layer_number]['r']=radius-layer_soln[0]

		if(layer_number<recipe['layers']-1):
			if(T is None):
				phase=phases.get_phase(recipe, layer_number+1, None, P, 300)

				radius=radius-layer_end_x
				P=layer_end_conditions[3]
				rho=phases.vinet_density(P, phase)
				g=layer_end_conditions[1]
				M=layer_end_conditions[2]
			else:
				print "Gimme time"

	# Interp all my fields for use in the temperature ODE
	solution=interp_solution(solution)

	return solution

def del_mass_event():
	try:
		del mass_event
	except:
		pass

def define_mass_event(m):
	
	del_mass_event()

	mass_template = """
def mass_event(Y,x):
	rho, g, m, P = Y
	value= m-%f
	isterminal = True
	direction  = 0
	return value, isterminal, direction
	"""

	exec(mass_template % (m)) in globals()

def define_radius_event(x):

	radius_template = """
		def radius_event(Y,x):
			rho, g, m, P = Y
			value= x-%f
		    isterminal = True
		    direction  = 0
		    return value, isterminal, direction
	"""

	exec(radius_template % (x))

def get_layer(recipe, m):

    layer=0
    layer_masses=recipe['layer_masses'][layer]
    m_above=recipe['mass']-m
    
    while layer<(recipe['layers']-1) and (m_above>layer_masses):
        layer=layer+1
        layer_masses+=recipe['layer_masses'][layer]

    return layer+1


def structure_equations(y, x, recipe, phase, R, Temperature):
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

	alpha=phase['alpha']
	gamma=phase['gamma_0']
	rho_0=phase['rho_0']
	K_0=phase['K_0']
	Kp_0=phase['Kp_0']

	if (Temperature is False):
		# phase=phases.get_phase(recipe, layer_number, None, recipe['P_surf'], recipe['T_surf'])
		phi=K_S(r, rho, rho_0, False, K_0, Kp_0, alpha, gamma)/rho
	else:
		print "Gimme time"
		# phi=K_S(r, rho, rho_0, T, K_0, Kp_0, alpha, gamma)/rho

	# print str(x)+", "+phase['name']+", "+str(rho)+", "+str(g)+", "+str(m)+", "+str(P)

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

