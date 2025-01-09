#! /usr/bin/env python

import argparse
from MDPropTrack.analysis import PropertyAnalyser
from MDPropTrack.analysis import LipidPropertyCalculator, ProteinPropertyCalculator
import matplotlib.pyplot as plt

def ParseArguments():
	"""
	Parse arguments for the MDPropTrack run
	"""

	parser = argparse.ArgumentParser(
		description="MDPropTrack. Extract and plot properties from GROMACS simulation."
	)

	### KEY INPUT PARAMETERS ###
	parser.add_argument(
		"-e",
		"--edr",
		help = """Path to edr file(s) comma-separated
		'path/to/step1.edr,path/to/step2.edr'""",
		type = str
	)

	parser.add_argument(
		"-t",
		"--trj",
		help = """Path to trajectory file(s) comma-separated
		'path/to/step1.xtc,path/to/step2.xtc'""",
		type = str
	)

	parser.add_argument(
		"-p",
		"--topol",
		help = """Path to topology file, usually 
		neccessary for trajectory analysis""",
		type = str
	)

	parser.add_argument(
		"-o",
		"--output",
		help = """Otput figure. E.g. 'save/as/multiplot.jpg'""",
		type = str
	)

	### TRAJECTORY ANALYSIS PARAMETERS ###

	parser.add_argument(
		"--center_group",
		help = """Atom selection to center trajectory""",
		default = 'protein',
		type = str
	)

	parser.add_argument(
		"--fit_group",
		help = """Atom selection for rot-trans fit""",
		default = 'protein',
		type = str
	)

	parser.add_argument(
		"--tu",
		help = """Time units, ns or ps. Default ns""",
		default = 'ns',
		choices = ['ns', 'ps'],
		type = str
	)

	parser.add_argument(
		"--step",
		help = """Step for trajectory analysis. Default 1""",
		default = 1,
		type = int
	)

	parser.add_argument(
		"--sequential",
		help = """Treat files as sequential steps. Default True""",
		default = True,
		type = bool
	)

	parser.add_argument(
		"-v", 
		"--verbose",
		help = """Default False""",
		action='store_true'
	)

	### PROPERTIES FOR CALCULATION ###
	parser.add_argument(
		"--apl",
		help = """Average area per lipid. Provide lipid selection for leaflet identification
		and Voronoi tesselation. E.g. 'name PO4' for MARTINI phospholipids.""",
		type = str
	)

	parser.add_argument(
		"--thickness",
		help = """Average bilayer thickness. Provide lipid selection
		for leaflet identification E.g. 'name PO4' for MARTINI phospholipids.""",
		type = str
	)

	parser.add_argument(
		"--scc",
		help = """Average order parameter. Provide one or two lipid tail selections
		E.g. 'name ??A,name ??B' for sn1 and sn2 lipid tails in MARTINI.""",
		type = str
	)

	parser.add_argument(
		"--rg",
		help = """Radius of gyration. Provide atom selection for Rg calculation.
		E.g. 'protein'""",
		type = str
	)

	parser.add_argument(
		"--rmsd",
		help = """RMSD for an atom group. Provide two comma-separated selections
		for a leas-square fit and for rmsd calculation.
		E.g. 'backbone,resid 12-40'""",
		type = str
	)

	### PROPERTIES TO PLOT ###

	parser.add_argument(
		"--plot",
		help = """Comma-separated list of properties to plot. Will override the default list.
		Default: 'Potential,Temperature,Pressure,Volume' + requested trj properties.""",
		default = 'Potential,Temperature,Pressure,Volume',
		type = str
	)

	parser.add_argument(
		"--plot_convergence",
		help = """Plot autocorrelation time vs simulation time""",
		action='store_true'
	)

	args = parser.parse_args()

	return args
	
def DefineFuncs(args):
	"""
	Define list of functions and respective names
	for PropertyAnalyser
	"""

	funcs, func_names = [], []

	if args.apl is not None:
		funcs.append(
			LipidPropertyCalculator(lipid_sel=args.apl).CalcAreaPerLipid
		)
		func_names.append('APL')

	if args.thickness is not None:
		funcs.append(
			LipidPropertyCalculator(lipid_sel=args.thickness).CalcBilayerThickness
		)
		func_names.append('Bilayer Thickness')

	if args.scc is not None:
		tail_sel = args.scc.split(',') if ',' in args.scc \
				   else args.scc
		funcs.append(
			LipidPropertyCalculator(tail_sel=tail_sel).CalcOrderParameter
		)
		func_names.append('Lipid tail order parameter')

	if args.rg is not None:
		funcs.append(
			ProteinPropertyCalculator(protein_sel=args.rg).CalcGyrationRadius
		)
		func_names.append('Rg')

	if args.rmsd is not None:
		sel1, sel2 = args.rmsd.split(',')
		funcs.append(
			ProteinPropertyCalculator(protein_sel=sel1, fit_sel=sel2).CalcAreaPerLipid
		)
		func_names.append('RMSD')

	return funcs, func_names

def main():
	"""
	MDPropTrack run ./	
	args - arguments parsed by ParseArguments()
	"""

	# parse arguments
	args = ParseArguments()

	# define functions for trajectory analysis
	funcs, func_names = DefineFuncs(args)

	# initiate PropertyAnalyser
	pa = PropertyAnalyser(
		
		edr = args.edr.split(','),
		trj = args.trj.split(','),
		topol = args.topol,
		
		funcs = funcs,
		func_names = func_names,

		center_group = args.center_group,
		rot_trans_group = args.fit_group
	)

	# run analyses
	pa.extract_properties(
		tu = args.tu,
		step = args.step,
		sequential = args.sequential,
		verbose = args.verbose
	)

	# plot
	properties_to_plot = args.plot.split(',')
	if func_names is not None:
		properties_to_plot.extend(func_names)

	pa.plot(
		properties_to_plot = properties_to_plot,
		plot_convergence=args.plot_convergence,
		x_lab = f"Time, {args.tu}"
	)

	# save
	output = 'properties_timeplot.jpg' \
			 if args.output is None else args.output

	plt.savefig(
		output,
		dpi=300,
		bbox_inches='tight'
	)
	
	plt.close()

	return 

if __name__ == "__main__":
	
	# launch as a command-line tool
	main()
