import numpy as np
import math

def ironsolidus(rho, config):
	Tm0=float(config.get('ironsolidus', 'Tm0'))
	rho0=float(config.get('ironsolidus', 'rho0'))
	gamma0=float(config.get('ironsolidus', 'gamma0'))
	q=float(config.get('ironsolidus', 'q'))
	
	solidus=Tm0*math.pow(rho0/rho, 2.0/3.0)
	solidus*=math.exp((2*gamma0/q)*(1-math.pow(rho0/rho,q)))

	return solidus
