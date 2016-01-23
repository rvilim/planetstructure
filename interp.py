from scipy.interpolate import interp1d

def interp_solution(solution):
	# This function takes in a solution profile and interpolates it for use in the T ode

	layer_num=0

	for layer in solution:
		rho=interp1d(layer['r'], layer['y'][:,0], kind='cubic')
		g=interp1d(layer['r'], layer['y'][:,1], kind='cubic')
		m=interp1d(layer['r'], layer['y'][:,2], kind='cubic')
		P=interp1d(layer['r'], layer['y'][:,3], kind='cubic')

		solution[layer_num]['y_interp']=[rho,g,m,P]

		layer_num=layer_num+1

	return solution