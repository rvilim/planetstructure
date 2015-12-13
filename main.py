import json


import phases
import params
import structure

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

def main():
	recipe_file="recipes/earth.json"

	recipe=get_recipe(recipe_file)
	recipe=phases.load_phases(recipe)

	params.get_params('params/params.ini')
	
	structure.structure(recipe, 6371000.0)

if __name__=="__main__":
	main()