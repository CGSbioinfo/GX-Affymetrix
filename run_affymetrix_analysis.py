import os
import argparse
import re
def read_parameters_file(params_file):
	true = re.compile("true", flags = re.IGNORECASE)
	false = re.compile("false", flags = re.IGNORECASE)
	params_file = params_file
	args = {}


	with open(params_file, 'r') as f:
		for line in f:
			entry = line.strip().split("=")
			if entry[0]:
				args[entry[0].strip(" ")] = entry[1].strip()
	global samples
	samples = args['sample file']
	global comparisons
	comparisons = args['comparisons file']
	global qccFile
	qccFile = args['qcc file']
	global rawDataFolder
	rawDataFolder = args['raw data folder']


	return (args)

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description = 'Run the Affymetrix GX pipeline')
	parser.add_argument('--arguments_file', help='Argument file. Default=arguments.txt', default='arguments.txt')
	args=parser.parse_args()
	fileArgs = read_parameters_file(args.arguments_file)
	#os.system('Rscript affy_pipeline_v1.2.R -s '+samples+' -c '+comparisons+' -q '+qccFile+' -r '+rawDataFolder)
	os.system('Rscript generate_report.R')