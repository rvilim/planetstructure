from materialproperties
import K_S, gamma
import params

def temperature(structure):
	for layer in structure:
		x=np.linspace(0, layer['top_r']-layer['bottom_r']), params.layerpoints)


def temperature_equations(y,x,layer):
	"""
	y=[ T ]
	"""

	T=y[0]

	r=layer['top_r']-x
	phase=layer['phase']

	rho=layer['y_interp']['rho'](r)
	g=layer['y_interp']['g'](r)

	gamma=get_gamma(phase['gamma_0'], phase['rho_0'], rho, phase['q'])
	K_S=get_K_S(r, rho, phase['rho_0'], T, phase['K_0'], phase['Kp_0'], phase['alpha'], gamma)
	
	layer['upper_bl_end']
	layer['lower_bl_start']

	if( not in boundary layer):
		dTdr=-rho*g*T*gamma/K_S
	else:
		dTdr=-q/k