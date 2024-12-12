#! /usr/bin/env python

import os
import glob
import shutil
import argparse

######################################
import sys
sys.path.append('path/to/smth')
from analysis import PropertyAnalyser
######################################

def ParseArguments():

	parser = argparse.ArgumentParser(
		description="MemAssembler. A tool to setup and analyse coarse-grained membrane self-assembly simulations with GROMACS."
	)

	parser.add_argument(
		"-f",
		"--file",
		help="""Path to PDB file (e.g. 'path/to/protein1.pdb')""",
		type=str,
		# required=True,
	)

	# parser_setup.add_argument(
	# 	"-f",
	# 	"--file",
	# 	help="""Multiple files for processing?""",
	# 	type=str,
	# 	required=True,
	# )

	args = parser.parse_args()

	return args
	
def main(args):
	"""
	args - arguments parsed by ParseArguments()
	"""

	# do analysis and plotting and everything

	return 

if __name__ == "__main__":
    
	# launch as a command-line tool

	# parse arguments
	arguments = ParseArguments()

	# run analysis and plotting
	main(arguments)
