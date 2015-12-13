from shapely.wkt import loads
import shapely
import math

def load_phases(recipe):
	# print recipe['phases']
	for layer_number in xrange(1,recipe['layers']+1):

		layer_name='layer_'+str(layer_number)

		for i, phase in enumerate(recipe[layer_name]['phases']):
			if('geom' in phase):
				recipe[layer_name]['phases'][i]['geom']=loads(phase['geom'])

	return recipe

def get_phase(recipe, layer, rho, P, T):
	# So this function is a bit complicated. It tries to return the correct phase for the
	# current thermodynamic equations. It's complicated because we can either specify Lindeman
	# phase relations (for calculating melting/solid) or a geom type (for calculating different
	# phases within a solid). 

	# If the layer is set as linedmans then we need a densit. Otherwise we don't. The whole
	# rho is not None thing comes from the fact that at the surface which I am assuming to be
	# solid, we need an initial density. We get that from our GEOM phase relations.

	# This could probably be improved.

	point=shapely.geometry.Point(P/pow(10,9.0), T)

	layer=recipe['layer_'+str(layer)]

	if ((layer['type'].upper()=="LINDEMAN") and (rho is not None)):
		phaseindex={}
		for i, phase in enumerate(layer['phases']):
			phaseindex[phase['phase'].upper()]=i

			if(phase['phase'].upper()=='SOLID'):
				gamma_0=phase['gamma_0']
				rho_0=phase['rho_0']
				q=phase['q']
				Tm_0=phase['Tm_0']
				T_m=Tm_0*pow(rho_0/rho,2.0/3)*math.exp((2*gamma_0/q)*(1-pow(rho_0/rho,q)))

		print "Melting Tempereature", rho, P, T, T_m
		if(T>=T_m):
			return layer['phases'][phaseindex['LIQUID']]
		else:
			return layer['phases'][phaseindex['SOLID']]

	elif (layer['type'].upper()=="GEOM"):
		for phase in layer['phases']:
			if phase['geom'].contains(point):
				return phase

	print "Unknown phase or GEOM type is not found for layer "+str(layer)

