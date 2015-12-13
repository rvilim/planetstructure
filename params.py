import ConfigParser
import os

G=0.0

def get_params(params_file):
	global G

	parser=ConfigParser.ConfigParser()

	params_filename="params.ini"
	params_filepath=os.path.join(os.path.dirname(__file__),"params",params_filename)
	parser.read(params_filepath)

	G=float(parser.get('Constants', 'G'))