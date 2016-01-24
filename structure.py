# import scipy.integrate as integrate
from pycse import odelay

import numpy as np

import phases
import params
import math
import matplotlib.pyplot as plt
from interp import interp_solution
from materialproperties import get_K_S, get_gamma, get_Ra, get_param_avg, get_eta

from scipy.interpolate import interp1d

def structure(R):

	M=params.recipe['mass']

	solution=get_profile(R)

	# print solution
	plot_solution(solution)

def plot_solution(solution):
	r=[]
	rho=[]
	g=[]
	m=[]
	P=[]

	a=[]

	f, ((ax1, ax2),(ax3, ax4)) = plt.subplots(2, 2)
	# ax1.hold(True)	

	ax1.set_title('Rho')
	ax2.set_title('g')
	ax3.set_title('m (below)')
	ax4.set_title('P')

	# for layer in solution:
	# 	r=r+[a/1000.0 for a in layer['r']]
	# 	rho=rho+[y[0] for y in layer['y']]
	# 	g=g+[y[1] for y in layer['y']]
	# 	m=m+[y[2] for y in layer['y']]
	# 	P=P+[y[3] for y in layer['y']]

	# ax1.plot(r,rho)
	# ax2.plot(r,g)
	# ax3.plot(r,m)
	# ax4.plot(r,P)

	for layer in solution:

		ax1.plot(layer['r'],layer['y_interp']['rho'](layer['r']))
		ax2.plot(layer['r'],layer['y_interp']['g'](layer['r']))
		ax3.plot(layer['r'],layer['y_interp']['m'](layer['r']))
		ax4.plot(layer['r'],layer['y_interp']['P'](layer['r']))

	plt.show()

def get_profile(R):
	M=params.recipe['mass']

	g=params.G*M/pow(R,2.0)
	P=params.recipe['P_surf'] 
	T=None
	phase=phases.get_phase(0, None, P, params.recipe['T_surf'])
	rho=phase['rho_0']

	layersbelow_mass=params.recipe['mass']
	radius=R

	solution=[None]*2

	for layer_number in xrange(0,params.recipe['layers']):

		x=np.linspace(0, radius, params.layerpoints)

		layersbelow_mass=layersbelow_mass-params.recipe['layer_masses'][layer_number]
		define_mass_event(layersbelow_mass)

		layer_solution=odelay(structure_equations, [rho,g,M,P], x, events=[mass_event], args=(phase, radius, False))

		if(T is None):
			# If this is our first iteration we don't have an initial temperature profile yet
			# create a constant profile of temperature T_surf. 
			# This means phase relations won't go insane in subsaquent steps

			T=np.zeros_like(layer_solution[0])+params.recipe['T_surf']

		solution[layer_number]=get_layer_properties(layer_solution, T, layer_number, radius)

		if(layer_number<params.recipe['layers']-1):
			phase=phases.get_phase(layer_number+1, None, P, 300)

			radius=solution[layer_number]['bottom_r']
			P=solution[layer_number]['bottom_P']
			rho=phases.vinet_density(P, phase)
			g=solution[layer_number]['bottom_g']
			M=solution[layer_number]['bottom_m']

	# Interp all my fields for use in the temperature ODE
	solution=interp_solution(solution)

	return solution

def get_layer_properties(layer_solution, T, layer_number, radius):
	layer={}
	layer['T']=T
	layer['X']=layer_solution[0]
	layer['y']=layer_solution[1]
	layer['r']=radius-layer_solution[0]
	layer['top_r']=radius
	layer['bottom_r']=radius-layer_solution[2][0]
	layer['bottom_rho']=layer_solution[3][0][0]
	layer['bottom_g']=layer_solution[3][0][1]
	layer['bottom_m']=layer_solution[3][0][2]
	layer['bottom_P']=layer_solution[3][0][3]

	layer['g_avg']=-np.trapz(layer['y'][:,1]*(layer['r']**2), x=layer['r'])*3.0/(layer['top_r']**3-layer['bottom_r']**3)
	layer['rho_avg']=-np.trapz(layer['y'][:,0]*(layer['r']**2), x=layer['r'])*3.0/(layer['top_r']**3-layer['bottom_r']**3)


	layer['alpha_avg']=get_param_avg(layer, layer_number, 'alpha')
	layer['k_avg']=get_param_avg(layer, layer_number, 'k')
	layer['kappa_avg']=get_param_avg(layer, layer_number, 'kappa')
	layer['eta_avg']=get_param_avg(layer, layer_number, 'eta_0')


	layer['Ra_avg']=get_Ra(layer, layer_number)
	print 'Ra', layer['Ra_avg']

	# print layer['g_avg'], layer['rho_avg'], layer['alpha_avg'], layer['k_avg'], layer['kappa_avg'], layer['eta_avg']

	return layer

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

def get_layer(m):

    layer=0
    layer_masses=params.recipe['layer_masses'][layer]
    m_above=params.recipe['mass']-m
    
    while layer<(params.recipe['layers']-1) and (m_above>layer_masses):
        layer=layer+1
        layer_masses+=params.recipe['layer_masses'][layer]

    return layer+1


def structure_equations(y, x, phase, R, Temperature):
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
	rho_0=phase['rho_0']
	K_0=phase['K_0']
	Kp_0=phase['Kp_0']
	gamma=get_gamma(phase['gamma_0'], phase['rho_0'], rho, phase['q'])

	if (Temperature is False):
		phi=get_K_S(r, rho, rho_0, False, K_0, Kp_0, alpha, gamma )/rho
	else:
		print "Gimme time"
		# phi=K_S(r, rho, rho_0, T, K_0, Kp_0, alpha, gamma)/rho


	drhodx=-(-rho*g/(phi*pow(10,9)))
	dgdx=-(4*np.pi*params.G*rho-2*params.G*m/pow(r,3))
	dmdx=-(4*np.pi*pow(r,2)*rho)
	dPdx=-(-rho*g)

	return [drhodx, dgdx, dmdx, dPdx]


