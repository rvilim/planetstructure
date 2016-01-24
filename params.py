import ConfigParser
import os
import json
import phases

G=0.0
layerpoints=100
recipe=None

def get_params(params_file, recipe_file):
	global G
	global recipe
	global layerpoints
	
	parser=ConfigParser.ConfigParser()

	params_filename="params.ini"
	params_filepath=os.path.join(os.path.dirname(__file__),"params",params_filename)
	parser.read(params_filepath)

	G=float(parser.get('Constants', 'G'))
	layerpoints=int(parser.get('Model', 'layerpoints'))

	recipe=get_recipe(recipe_file)

	# print recipe
	
	recipe=phases.load_phases(recipe)


def get_recipe(recipe_file):
	with open(recipe_file) as data_file:    
		recipe = json.load(data_file)

	return calc_recipe(recipe)

def calc_recipe(recipe):
	# In here goes things that should be calculated from the recipe file
	# and are completely determined by the recipe file

	# layer masses (calculatable from planetary mass and layer fracions)

	recipe['layer_masses']=[recipe['mass']*fraction for fraction in recipe['layer_fractions']]

	return recipe