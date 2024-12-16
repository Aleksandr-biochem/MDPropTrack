import panedr
import numpy as np
import pandas as pd
import seaborn as sns
from tqdm import tqdm
import MDAnalysis as mda
import lipyphilic as lpp
import matplotlib.pyplot as plt

class PropertyAnalyser:
	"""
	Extract properties of a simulation system
	from edr files and/or trajectories
	"""

	def __init__(self, edr=None, trj=None, topol=None, funcs=None, func_names=None):
		"""
		edr - str or list(str), path or list of paths
		to edr files in a sequential order
	
		trj - str or list(str), path or list of paths
		to trajectory files in a sequential order

		topol - str, path to a topology file 
		to facilitate MDAnalysis Universe loading
		optional but may be required for your custom functions

		funcs - list(functions),
		list of functions to apply along the trajectory
		Functions are expected to adhere to a specific input-output structure
		See examples for details
		
		func_names - list(str), names for the properties
		computed by funcs
		"""

		# energy and trajectory files for analysis
		self.edr = self._check_type(edr)
		self.trj = self._check_type(trj)
		self.topol = self._check_type(topol)

		# functions to be applied along the trajectory
		self.funcs = funcs
		self.func_names = func_names

		####################
		# transformations for trajectories
		# self.transformations
		####################
		
		# pandas DataFrame with all the data 
		self.data = None
		
	def _check_type(self, var):
		"""
		Check that var is None, str or list(str)
		If str, convert to list(str)
	
		Returns
		list(srt)
		"""
	
		# is None
		if var is None:
			checked_var = var
	
		# is str
		elif isinstance(var, str):
			checked_var = [var]
	
		# if it's list are all elements str 
		elif isinstance(var, list):
			check_list = [isinstance(v, str) for v in var]
			if False not in check_list:
				checked_var = var
			else: 
				raise Exception(f"input should be str or list(str)")
				
		# if neither then raise exception
		else:
			raise Exception(f"input should be str or list(str)")
	
		return checked_var

	def _read_edrs(self, tu, sequential):
		"""
		Read data from edr files and 
		append it to self.data

		tu - str, time units option, ns or ps
		
		sequential - bool, if True then supplied files
		are considered sequential step and the Time 
		is adjusted accordingly
		deffault, True

		Returns self
		"""

		# to make continuous timeline
		# from several edr files
		end_time = 0

		for edr in self.edr:

			# read edr data into pd.DataFrame
			df = panedr.edr_to_df(edr)

			# convert time if needed
			if tu == 'ns':
				df['Time'] = df['Time'] / 1000
				
			# update time for sequential steps
			if sequential:
				step_duration = df.Time.iloc[-1]
				df['Time'] += end_time
				end_time += step_duration
			
			# add step_name column
			df['Step_name'] = edr.split('/')[-1][:-4]
			
			# append to self.data
			self.data = pd.concat([self.data, df]).reset_index(drop=True)
			
		return self

	def _apply_funcs_trj(self, system, tu, step, verbose):
		"""
		Apply all functions to the trajectory
		Returns calculated properties as pd.DataFrame
		
		system - MDAnalysis Universe, trajectory for analysis
		tu - str, time units to use, 'ns' or 'ps'
		step - int, step for trajectory analysis
		verbose - bool, report trajectory analysis progress

		Returns
		trj_dat - pd.DataFrame with Time and Properties
		"""

		# data from trajectories 
		trj_dat = {
			'Time': [ts.time for ts in system.trajectory[::step]]
		}

		# apply each function to the trajectory
		for n, func in enumerate(self.funcs):

			# save property values
			trj_dat[self.func_names[n]] = func(
				system=system,
				step=step,
				verbose=verbose
			)

		# conver to pd.DataFrame
		trj_dat = pd.DataFrame.from_dict(trj_dat)
		
		return trj_dat
		
	def _analyse_trjs(self, tu, step, verbose, sequential):
		"""
		Calculate properties along the trajectories
		and appends them to self.data
		
		tu - str, time units to use, 'ns' or 'ps'
		
		step - int, step for trajectory analysis

		verbose - bool, report trajectory analysis progress

		sequential - bool, if True then supplied files
		are considered sequential step and the Time 
		is adjusted accordingly
		deffault, True
		
		Returns self	
		"""

		# define time units conversion
		t_norm = 1000 if tu == 'ns' else 1
		end_time = 0
		
		# list of DataFames with data
		# from trajectories
		trj_dat_combined = []
		
		# iterate over trajectories 
		for trj in self.trj:

			# load trajectory
			system = mda.Universe(trj) if self.topol is None \
					 else mda.Universe(self.topol[0], trj)

			#################
			# apply transformations
			################# 

			# calculate properties over the trajectory
			trj_dat = self._apply_funcs_trj(
				system = system,
				tu = tu,
				step = step,
				verbose = verbose
			)
			
			# adjust time for sequential steps
			if sequential:
				trj_dat['Time'] = trj_dat['Time'] / t_norm
				system.trajectory[-1]
				step_duration = system.trajectory.time / t_norm
				trj_dat['Time'] += end_time
				end_time += step_duration
			
			# add step_name column
			trj_dat['Step_name'] = trj.split('/')[-1][:-4]

			# add to list
			trj_dat_combined.append(trj_dat)
		
		# concatenate data from all trajectories
		trj_dat_combined = pd.concat(trj_dat_combined).reset_index(drop=True)
		
		# bind trj_data do self.data on time
		if self.data is not None:
			self.data = pd.merge(
				self.data,
				trj_dat_combined,
				on = ['Time', 'Step_name'],
				how = 'outer'
			)
		# of just save trj data as self.data
		else:
			self.data = trj_dat_combined
		
		return self
	
	def extract_properties(self, tu='ns', step=1, sequential=True, verbose=False):
		"""
		Extract data from edr and/or trj files
		
		tu - str, time units option, ns or ps
		default ns

		For trajectory analysis only:

		step - int, step for trajectory analysis
		default 1
		
		sequential - bool, if True then supplied files
		are considered sequential step and the Time 
		is adjusted accordingly
		deffault, True

		verbose - bool, show traj analysis progress
		default False
		
		Returns self
		"""

		# check time units
		if tu not in ['ns', 'ps']:
			raise Exception("Unrecognised tu input")
		
		# check that one of the inputs is there
		if (self.edr is None) and (self.trj is None):
			raise Exception("Neither edr or trj file were supplied")
		
		# extract properties from edr
		if self.edr is not None:
			self._read_edrs(
				tu = tu,
				sequential = sequential
			)

		# extract properties from trajectories
		if self.trj is not None:

			if self.funcs is None:
				raise Exception(f"Trajectiries supplied without functions for anaysis. Provide `func` argument")

			self._analyse_trjs(
				tu = tu,
				step = step,
				sequential = sequential, 
				verbose = verbose
			)
		
		return self

	def _construct_multiplot(self, n_prop, figure_kwargs):
		"""
		Construct subplot

		n_prop - int, number of properties to plot
		figure_kwargs - dict, matplotlib figure kwargs

		Return
		fig, axs - mplt Figure and list(Axes)
		"""

		# custom parameters for multiplot
		if figure_kwargs is not None:
			fig, axs = plt.subplots(**figure_kwargs)
		
		# define from n_prop
		else: 
			
			if n_prop > 1:
				ncols = 2
				nrows = (n_prop // 2) + (n_prop % 2)
			else:
				ncols = 1
				nrows = 1

			fig, axs = plt.subplots(
				ncols = ncols,
				nrows = nrows,
				figsize = (10 * ncols, 5 * nrows)
			)

		# transform to 1D array of Axes
		if isinstance(axs, np.ndarray):
			axs = axs.flatten()
		else:
			axs = np.array([axs])
					
		return fig, axs
		
	def plot(self,
			 properties_to_plot = ['Potential', 'Temperature', 'Pressure', 'Volume'],
			 x_lab='Time, ns',
			 cmap='Set1',
			 figure_kwargs=None, 
			 style_kwargs={"style": "darkgrid", "rc": {"grid.color": ".6", "grid.linestyle": ":"}},
			 sns_kwargs={'alpha': 0.7}):
		"""
		Plot properties from self.data using the 'Time' column as X-axis
		
		properties_to_plot - list(str) list of features to plot
		Only works with column names from self.data
		default list: 'Total Energy', 'Temperature', 'Pressure', 'Volume'
		
		x_lab - str, x axis label, defalut 'Time, ns'
		
		figure_kwargs - dict, matplotlib figure kwargs
		These can be used to adjust output subplot
		
		style_kwargs - dict, seaborn style kwargs
		sns_kwargs - dict, seaborn lineplot kwargs
		
		Returns:
		matplotlib figure object
		"""
		
		# do we have data to plot
		if self.data is None:
			raise Exception("self.data is None")

		# check whether requested properties
		# are in self.data
		prop_list = []
		for prop in properties_to_plot:
			if prop not in self.data.columns:
				print(f"Skipping {prop}, not in self.data.columns")
			else:
				prop_list.append(prop)
				
		# apply seaborn style to the plot
		with sns.axes_style(**style_kwargs):
	
			fig, axs = self._construct_multiplot(
				n_prop = len(prop_list),
				figure_kwargs = figure_kwargs
			)
				
			# plot each property on a different subplot
			for i, prop in enumerate(prop_list):
			
				sns.lineplot(
					data = self.data,
					x = 'Time',
					y = prop,
					hue = 'Step_name',
					palette = cmap,
					ax = axs[i],
					**sns_kwargs
				)
			
				# set title and axes labels 
				axs[i].set_title(prop, fontweight='bold', fontsize=18, pad=10)
				axs[i].set_xlabel(x_lab, fontsize=15, labelpad=10)
				axs[i].set_ylabel(prop,  fontsize=15, labelpad=10)
	
		plt.tight_layout()
		
		return fig, axs


class LipidPropertyCalculator:
	"""
	A class of with methods to calculate
	key lipid properties from the trajectory
	"""

	def __init__(self, lipid_sel, leaflet=0):
		"""
		Class atributes hold parameters to be used in methods

		lipid_sel - str, MDAnalysis selection for lipid group analysis
		
		leaflet - int, leaflets to use for property averaging 
		-1 - lower
		1  - upper
		0  - both
		"""
		self.lipid_sel = lipid_sel
		self.leaflet = leaflet

	def CalcAreaPerLipid(self, system, step=1, verbose=False):
		"""
		Calculate average area per lipid

		system - MDAnalysis Universe, trajectory for analysis
		step - int, step for trajectory analysis
		verbose - bool, report trajectory analysis progress

		Returns
		list(floats)
		"""

		# binning 
		n_bins = int(system.dimensions[0] // 10)
		
		# assign leaflets
		leaflets = lpp.leaflets.assign_leaflets.AssignLeaflets(
		  universe = system,
		  lipid_sel = self.lipid_sel,
		  n_bins = n_bins 
		)
		leaflets.run(
			step = step,
			verbose = verbose,
		)
		
		# compute apl
		apl = lpp.analysis.area_per_lipid.AreaPerLipid(
			universe = system,
			lipid_sel = self.lipid_sel,
			leaflets = leaflets.leaflets
		)
		apl.run(
			step = step,
			verbose = verbose,
		)

		# choose the leaflet(s) to compute metrics for
		if self.leaflet == 0:
			leaflet_vals = [-1, 1]
		else:
			leaflet_vals = [self.leaflet]

		# compute mean in selected group for each frame
		mask = np.isin(leaflets.leaflets, leaflet_vals)
		apl_by_frame = [
			np.nanmean(apl.areas[mask[:, i], i]) for i in range(apl.areas.shape[1])
		]
		
		return apl_by_frame
	
