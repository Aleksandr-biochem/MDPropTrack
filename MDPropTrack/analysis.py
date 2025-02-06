import panedr
import numpy as np
import pandas as pd
import seaborn as sns
from tqdm import tqdm
import MDAnalysis as mda
import lipyphilic as lpp
import matplotlib.pyplot as plt
from MDAnalysis import transformations

class PropertyAnalyser:
	"""
	Class formalising extraction and analysis of properties of a simulation system
	using edr files and/or trajectories as inputs
	
	Attributes
	----------

	simulations: str - path to .edr file
				 tuple(str, str) - paths to topology and trajectory files
				 mda.Universe - loaded and transformed Universe
				 list(str, tuple(str, str), mda.Universe)
		simulation data for analysis

	tags: list(str), dict(str: list(str))
		tags describing each simulation step

	funcs: list(functions)
		list of functions to apply along the trajectory
		Functions are expected to adhere to a specific input-output structure
		
	func_names: list(str)
		names for the properties computed by funcs

	transformation arguments:
		center_group: str
			atom selection to center in the box
		rot_trans_group: str
			atom selection to fit alonmg the trj
	"""

	def __init__(
		self,
		simulations=None,
		tags=None,
		funcs=None,
		func_names=None
	):

		# energy files, trajectory files and Universe instances for analysis
		self.simulations = self._check_input(simulations)
		
		# autofill names for simulation steps
		if tags is None:
			self.tags = [
				{'name': f"Simulation_{i}"} for i in range(1, len(self.simulations) + 1)
			]
		# tags from input
		else:
			self.tags = self._transform_tags(tags)

		# functions to be applied along the trajectory/Universe
		self.funcs = funcs
		self.func_names = func_names

		# pandas DataFrame with extracted data
		self.data = None
	
	@classmethod
	def _check_var_type(self, var):
		"""
		Check that var is one of the supported types:
		- None
		- str
		- tuple(str, str)
		- mda.Universe
	
		Returns
		----------
		check: bool
			True if correct type
		"""
	
		# is None
		if var is None:
			check = True
	
		# is str or mda. Universe
		elif isinstance(var, str) \
		  or isinstance(var, mda.Universe):
			check = True

		# tuple(str, str)
		elif isinstance(var, tuple) and \
		     isinstance(var[0], str) and isinstance(var[1], str):
		    check = True
				
		# if neither then raise exception
		else:
			check = False
	
		return check

	@classmethod
	def _check_input(self, var):
		"""
		Format input to list of the supported types:
		- None
		- str
		- tuple(str, str)
		- mda.Universe
	
		Returns
		----------
		list(str, tuple(str, str), mda.Universe)
		"""

		# if list check taht every element is of supported type
		if isinstance(var, list):
			for v in var:
				if self._check_var_type(v):
					continue
				else: 
					raise Exception("Couldn't process `simulations` arg. Check types")
		
		# if one element, check and convert to list
		elif self._check_var_type(var):
			checked_var = [var]

		# if neither then raise exception
		else:
			raise Exception("Couldn't process `simulations` arg. Check types")
	
		return checked_var

	@classmethod
	def _transform_tags(self, tags):
		"""
		Convert the format of tags so that every simulation
		has a corresponding dict of tags
		
		Returns
		----------
		list(dict(str: str))
		"""

		# if list of names
		if isinstance(tags, list):
			tags = [ {'name': str(tag)} ]

		# if dictionary, transform it
		elif isinstance(tags, dict):
			tags = [
				dict(zip(tags, i)) for i in zip(*tags.values())
			]

		else:
			raise Exception("Couldn't process `tags`. Check format")

		return tags


	def _get_transformations(self, system):
		"""
		Define transformations for the trajectory
		
		system - MDAnalysis Universe
		"""

		# unwrap atoms by default
		ag = system.atoms
		workflow = [transformations.unwrap(ag)]

		# check additional transformation groups
		if self.center_group is not None:
			
			center_group = system.select_atoms(self.center_group)

			# add centering
			if len(center_group.atoms) > 0:
				workflow.extend([
					transformations.center_in_box(center_group),
					transformations.wrap(ag, compound='residues')
				])	

		if self.rot_trans_group is not None:
			
			fit_group = system.select_atoms(self.rot_trans_group)
			system_ref = system.copy()
			fit_group_ref = system_ref.select_atoms(self.rot_trans_group)

			# add rotation and translation
			if len(fit_group.atoms) > 0:
				workflow.extend([
					transformations.fit_rot_trans(fit_group, fit_group_ref)
				])

		return workflow

	def _read_edr(self, file, tags, tu, sequential):
		"""
		Read data from edr file and append it to self.data

		Parameters
		----------
	
		file: str
			path to edr file

		tu: str
			time units option, ns or ps
		
		sequential:
			bool, if True then `Time` is adjusted when appending to self.data
			to indicate sequential simulations

		Returns
		----------
		self
		"""

		# read edr data into pd.DataFrame
		df = panedr.edr_to_df(edr)

		# convert time to ns if needed
		if tu == 'ns':
			df['Time'] = df['Time'] / 1000
			
		# shift Time for sequential steps
		if sequential and (self.data is not None):
			df['Time'] += self.data.Time.iloc[-1]
		
		# add tag columns
		for tag in tags:
			df[tag] = tags[tag]
		
		# append to self.data
		self.data = pd.concat([self.data, df]).reset_index(drop=True)
			
		return self

	def _apply_funcs_trj(self, system, step, verbose):
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

		# initiate data with 'Time' column
		
		# time step
		dt = system.trajectory.dt

		# start and end time
		system.trajectory[0]
		time1 = system.trajectory.time
		system.trajectory[-1]
		time2 = system.trajectory.time

		# generate the times
		trj_dat = np.arange(time1, time2 + dt, dt * step).reshape((-1, 1))

		# apply each function to the trajectory
		for func in self.funcs:

			# save property values
			vals = func(
				system=system,
				step=step,
				verbose=verbose
			)

			# check shape
			if len(vals.shape) == 1:
				vals = vals.reshape((-1, 1))

			trj_dat = np.concatenate((trj_dat, vals), axis=1)

		# return default property names if none supplied
		# or if the number of names is incorrect
		if (self.func_names is None) or \
		   (len(self.func_names) != (trj_dat.shape[1] - 1)):
			self.func_names = [
				f"Prop{i}" for i in range(1, trj_dat.shape[1])
			]

		# conver to pd.DataFrame
		trj_dat = pd.DataFrame(
			data = trj_dat,
			columns = ['Time'] + self.func_names
		)
		
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

			if verbose:
				print(f"Loading and transforming '{trj}'...")

			# load trajectory
			system = mda.Universe(trj) if self.topol is None \
					 else mda.Universe(self.topol[0], trj)
		
			# system.trajectory.add_transformations(*workflow)
			system.trajectory.add_transformations(
				*self._get_transformations(system)
			)

			# calculate properties over the trajectory
			trj_dat = self._apply_funcs_trj(
				system = system,
				step = step,
				verbose = verbose
			)
			
			# shift time for sequential steps
			if sequential:
				system.trajectory[-1]
				trj_dat['Time'] += end_time
				end_time += system.trajectory.time
			
			# add step_name column
			trj_dat['Step_name'] = trj.split('/')[-1][:-4]

			# add to list
			trj_dat_combined.append(trj_dat)
		
		# concatenate data from all trajectories
		trj_dat_combined = pd.concat(trj_dat_combined).reset_index(drop=True)
		
		# adjust time units
		trj_dat_combined['Time'] = trj_dat_combined['Time'] / t_norm

		# bind trj_data do self.data on time
		if self.data is not None:
			self.data = pd.merge(
				self.data,
				trj_dat_combined,
				on = ['Time', 'Step_name'],
				how = 'outer'
			)

		# or just save trj data as self.data
		else:
			self.data = trj_dat_combined
		
		return self

	def extract_properties(self, tu='ns', step=1, sequential=False, verbose=False):
		"""
		Extract data from edr files and trajectories

		Parameters
		----------
		
		tu - str
			convert to these time units; ns or ps, default ns

		For trajectory analysis only:

		step: int
			step for trajectory analysis, default 1
		
		sequential: bool
			if True then supplied `simulations` are considered sequential steps
			and `Time` column in .data is adjusted accordingly
			deffault, False

		verbose: bool
			verbose traj analysis progress, default False
		
		Returns
		----------
		self
		"""

		# check time units
		if tu not in ['ns', 'ps']:
			raise Exception("Unrecognised tu option")
		
		# check that one of the inputs is there
		if self.simulations is None:
			raise Exception("No simulation data provided")
		
		# analyse each simulation input
		for sim, tags in zip(self.simulations, self.tags):

			# single file is expected to be .edr
			if isinstance(sim, str):

				self._read_edr(
					file = sim,
					tags = tags,
					tu = tu,
					sequential = sequential
				)

			# analyse trajectory
			elif:

				# transform trj file into MDAnalysis Universe
				# if input_type == 'trj':

				# 	sim = # load universe

				self._analyse_trjs(
					sim
					tu = tu,
					step = step,
					sequential = sequential, 
					verbose = verbose
				)
		
		return self

	################################
	# The following group of methods for convergence assessment
	# is implemented with the code adapted from
	# https://emcee.readthedocs.io/en/stable/tutorials/autocorr/
	# Copyright 2012-2021, Dan Foreman-Mackey & contributors

	def _get_next_pow_two(self, n):
		"""
		Return right nearest 
		int that is power of 2
		"""
		i = 1
		while i < n:
			i = i << 1
		return i

	def _get_autocorr_func_1d(self, y, norm=True):
		"""
		Estimate autocorrelation function from
		time series data using one-dimensional
		discrete Fourier Transform
		
		y - np.array(floats), time series data
		norm - bool, normalise acf, default True

		Returns
		np.array(floats) empirical acf values
		""" 

		n = self._get_next_pow_two(len(y))

		# Compute the FFT and then
		# the auto-correlation function
		f = np.fft.fft(y - np.mean(y), n = (2 * n))
		acf = np.fft.ifft(f * np.conjugate(f))[: len(y)].real
		acf /= 4 * n

		# Normilise if requested
		if norm:
			acf /= acf[0]

		return acf
	
	def _auto_window(self, taus, c):
		"""
		Automated windowing procedure following Sokal (1989)
		
		taus - np.array of tau estimates
		c - float, coefficient in tau estimation

		Return
		int, window number
		"""
		m = np.arange(len(taus)) < c * taus
		if np.any(m):
			window = np.argmin(m)
		else:
			window = len(taus) - 1
		return window

	def _estimate_autocorr_tau(self, y, c=5.0):
		"""
		Estimate autocorrelation time tau
		from time-series
		
		y - np.array(floats), time series data
		c - float, coefficient in tau estimation
		default 0.5

		Returns
		"""

		# get acf estimate
		f = self._get_autocorr_func_1d(y)
		
		# compute and return tau estimate
		taus = 2.0 * np.cumsum(f) - 1.0
		window = self._auto_window(taus, c)
		
		return taus[window]

	def _estimate_convergence(self, prop):
		"""
		Estimate autocorrelation time (tau) vs simulation length
		to assess convergence

		prop - str, property name from self.data 

		Returns 
		two np.array(float)
		array of time points and taus 
		"""

		# get the time-series values
		dat = self.data[prop].values

		# generate step points and
		# time points for tau estimates
		N = np.exp(
			np.linspace(np.log(100), np.log(dat.shape[0]), 10)
		).astype(int)
		ts = self.data['Time'].values[N - 1]

		# estimate tau from trj slices
		tau_data = np.array([
			self._estimate_autocorr_tau(dat[:n]) for n in N
		])
		
		return ts, tau_data

	##################################
	
	def _construct_multiplot(self, n_prop, figure_kwargs):
		"""
		Construct subplot grid

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
			 properties_to_plot=['Potential', 'Temperature', 'Pressure', 'Volume'],
			 plot_convergence=False,
			 labels=None,
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
		
		labels - list(str), list of cunstom names for plotted steps

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
				
				# for convergence
				if plot_convergence:
					ts, tau_data = self._estimate_convergence(prop)

					sns.lineplot(
						x = ts,
						y = tau_data,
						ax = axs[i],
						marker='o',
						**sns_kwargs
					)

					axs[i].set_ylabel(
						r"Autocorrelation time $\tau$", 
						fontsize=15,
						labelpad=10
					)

				# for regular time plots
				else:
					sns.lineplot(
						data = self.data,
						x = 'Time',
						y = prop,
						hue = 'Step_name',
						palette = cmap,
						ax = axs[i],
						**sns_kwargs
					)
					
					# change labels if requested
					if labels is not None:
						handles, previous_labels = axs[i].get_legend_handles_labels()
						axs[i].legend(
							handles = handles,
							labels = labels
						)

					axs[i].set_ylabel(prop,  fontsize=15, labelpad=10)

				# set title and axes labels 
				axs[i].set_title(prop, fontweight='bold', fontsize=18, pad=10)
				axs[i].set_xlabel(x_lab, fontsize=15, labelpad=10)
	
		plt.tight_layout()
		
		return fig, axs


class LipidPropertyCalculator:
	"""
	A class of with methods to calculate
	key lipid properties from the trajectory
	"""

	def __init__(
			self,
			lipid_sel=None,
			apl_sel=None,
			tail_sel=None,
			calculate=['apl', 'thickness', 'order_param'],
			filter_lipid=['all'],
			leaflet_to_average=0,
			bin_len_leaflets=10,
			bin_len_thickness=20
		):
		"""
		Class attributes hold parameters to be used in methods

		lipid_sel - str, MDAnalysis selection of a lipid group for leaflet identification
		and membrane thickness calculation (usually phosphate)
		
		apl_sel - atr, atom selection for the group to perform Voronoi Tesselation on

		tail_sel - str, MDAnalysis selection for lipid tails
		that will be used for order parameter calculation
		
		calculate - str or list(str), keywords of properties to calculate
		choices: ['apl', 'thickness', 'order_param'], default all 3

		filter_lipid - str or list(str), one or multiple atom selections
		to filter lipids in APL and order parameter calculation

		leaflet_to_average - int, leaflets to use for area per lipid averaging
		-1 - lower
		1  - upper
		0  - both

		bin_len_leaflets - float, bin width for leaflet identification, default 10
		bin_len_thickness - float, bin width for membrane thickness calculation, default 20
		"""
		self.lipid_sel = lipid_sel
		self.apl_sel=apl_sel
		self.tail_sel = tail_sel
		self.calculate = [calculate] if isinstance(calculate, str) \
						 else calculate
		self.filter_lipid = [filter_lipid] if isinstance(filter_lipid, str) \
						 	else filter_lipid
		self.leaflet_to_average = leaflet_to_average
		self.bin_len_leaflets = bin_len_leaflets
		self.bin_len_thickness = bin_len_thickness
		self.leaflets = None

	def _assign_leaflets(self, system, step=1, verbose=False):
		"""
		Run LiPyPhilic leaflet assignment over the trajectory
		Labels stored in self.leaflets
		"""

		# binning 
		n_bins_leaflets  = int(system.dimensions[0] // self.bin_len_leaflets)
		
		leaflets = lpp.leaflets.assign_leaflets.AssignLeaflets(
		  universe = system,
		  lipid_sel = self.lipid_sel,
		  n_bins = n_bins_leaflets 
		)

		if verbose:
			print('Assigning leaflets...')

		leaflets.run(
			step = step,
			verbose = verbose
		)

		self.leaflets = leaflets.leaflets

		return self

	def _filter_lipids(self, system, main_sel, filter_sel='all'):
		"""
		Generate a mask to filter lipid species from bilayer
		
		system - MDAnalysis Universe, trajectory for analysis
		main_sel - str, main atom selection for filtering
		filter_sel - str, atom selection to combine with main selection

		returns
		np.array(bool)
		"""

		# get ids of all lipid particles from the main selection
		lipid_ids = system.select_atoms(main_sel).residues.resids

		# modify selection with additional filter
		filtered_ids = system.select_atoms(
			main_sel + ' and ' + filter_sel
		).residues.resids
		
		# return bool mask
		return np.isin(lipid_ids, filtered_ids)

	def CalcProps(self, system, step=1, verbose=False):
		"""
		Calculate averaged properties from list:
		- area per lipid
		- bilayer thickness
		- orientational order parameter of lipid tails

		system - MDAnalysis Universe, trajectory for analysis
		step - int, step for trajectory analysis
		verbose - bool, report trajectory analysis progress

		Returns
		np.array(floats)
		"""

		props = []

		# Area Per Lipid
		if 'apl' in self.calculate:
			props.append(
				self.CalcAreaPerLipid(
					system=system,
					step=step,
					verbose=verbose
				)
			)

		# Bilayer Thickness
		if 'thickness' in self.calculate:
			props.append(
				self.CalcBilayerThickness(
					system=system,
					step=step,
					verbose=verbose
				)
			)

		# Order Parameter
		if 'order_param' in self.calculate:
			props.append(
				self.CalcOrderParameter(
					system=system,
					step=step,
					verbose=verbose
				)
			)

		return np.hstack(props)

	def CalcAreaPerLipid(self, system, step=1, verbose=False):
		"""
		Calculate average area per lipid over trajectory 

		system - MDAnalysis Universe, trajectory for analysis
		step - int, step for trajectory analysis
		verbose - bool, report trajectory analysis progress
		
		Also requires:
		self.lipid_sel - atom selection for lipids in the bilayer.
		These atoms will also be used to perform the Voronoi tessellation
		
		leaflet_to_average - int, leaflets to use for area per lipid averaging
		-1 - lower
		1  - upper
		0  - both

		Returns
		np.array(floats)
		"""

		# check leaflet assignment
		if self.leaflets is None:
			self._assign_leaflets(
				system=system,
				step=step,
				verbose=verbose
			)

		# configure apl calculation
		apl = lpp.analysis.area_per_lipid.AreaPerLipid(
			universe = system,
			lipid_sel = self.apl_sel,
			leaflets = self.leaflets
		)

		if verbose:
			print('Calculating area per lipid...')

		apl.run(
			step = step,
			verbose = verbose
		)

		# choose the leaflet(s) to compute metrics for
		if self.leaflet_to_average == 0:
			leaflet_vals = [-1, 1]
		else:
			leaflet_vals = [self.leaflet_to_average]
		mask_leaflet = np.isin(self.leaflets, leaflet_vals)

		# compute mean by frame in lipid groups
		apl_by_frame = []
		for filter_sel in self.filter_lipid:
			
			# define lipid mask
			mask_lipid = self._filter_lipids(
				system=system,
				main_sel=self.lipid_sel,
				filter_sel=filter_sel
			)

			# average apl in lipid group
			apl_by_frame.append(
				[np.nanmean(apl.areas[mask_leaflet[:, i] * mask_lipid, i]) \
				 for i in range(apl.areas.shape[1])]
			)
		
		return np.array(apl_by_frame).T
	
	def CalcBilayerThickness(self, system, step=1, verbose=False):
		"""
		Calculate bilayer_thickness over trajectory 

		system - MDAnalysis Universe, trajectory for analysis
		step - int, step for trajectory analysis
		verbose - bool, report trajectory analysis progress
		
		Also requires:
		self.lipid_sel - atom selection for lipids in the bilayer.
		Atoms used to identify leaflets.
		These atoms will also be used to define thickness.

		Returns
		np.array(floats)
		"""

		# check leaflet assignment
		if self.leaflets is None:
			self._assign_leaflets(
				system=system,
				step=step,
				verbose=verbose
			)

		# binning 
		n_bins_thickness = int(system.dimensions[0] // self.bin_len_thickness)
		
		# compute thickness
		memb_thickness = lpp.analysis.MembThickness(
			universe = system,
		 	leaflets = self.leaflets,
			lipid_sel = self.lipid_sel,
			n_bins = n_bins_thickness
		)

		if verbose:
			print('Calculating membrane thickness...')

		memb_thickness.run(
			step = step,
			verbose = verbose
		)

		return memb_thickness.memb_thickness.reshape((-1, 1))

	def CalcOrderParameter(self, system, step=1, verbose=False):
		"""
		Calculate average orientational order parameter
		of lipid tails over trajectory 

		system - MDAnalysis Universe, trajectory for analysis
		step - int, step for trajectory analysis
		verbose - bool, report trajectory analysis progress

		Also requires:
		self.tail_sel - atom selection(s) for lipid tails for 
		order parameter calculation

		Returns
		np.array(floats)
		"""

		# check leaflet assignment
		if self.leaflets is None:
			self._assign_leaflets(
				system=system,
				step=step,
				verbose=verbose
			)

		# check number of selections to work with
		if isinstance(self.tail_sel, list):
			sn1_sel, sn2_sel = self.tail_sel
		else:
			sn1_sel = self.tail_sel
			sn2_sel = None

		# run analysis for SN1 tail
		scc_sn1 = lpp.analysis.order_parameter.SCC(
			universe = system,
			tail_sel = sn1_sel
		)

		if verbose:
			print('Calculating order parameter for SN1 tail...')

		scc_sn1.run(
			step = step,
			verbose = verbose
		)

		# run for SN2 tail if any 
		# and return weighted average
		if sn2_sel is not None:

			scc_sn2 = lpp.analysis.order_parameter.SCC(
				universe = system,
				tail_sel = sn2_sel
			)

			if verbose:
				print('Calculating order parameter for SN2 tail...')

			scc_sn2.run(
				step = step,
				verbose = verbose
			)

			scc_av = lpp.analysis.order_parameter.SCC.weighted_average(scc_sn1, scc_sn2)

		# just one tail
		else:
			scc_av = scc_sn1

		# compute mean by frame for lipid groups
		scc_by_frame = []
		for filter_sel in self.filter_lipid:
			
			# define lipid mask
			mask_lipid = self._filter_lipids(
				system=system,
				main_sel=sn1_sel,
				filter_sel=filter_sel
			)

			# mask and average by frame
			scc_by_frame.append(
				np.nanmean(scc_av.SCC[mask_lipid, :], axis=0)
			)

		return np.array(scc_by_frame).T

class ProteinPropertyCalculator:
	"""
	A class of with methods to calculate
	key protein properties from the trajectory
	"""

	def __init__(self, protein_sel=None, fit_sel=None):
		"""
		Class atributes hold parameters to be used in methods

		protein_sel - str, MDAnalysis atom selection 
		for protein group for analysis
		
		fit_sel - str, atom selection for superimposition in CalcRMSD
		"""
		self.protein_sel = [protein_sel] if isinstance(protein_sel, str) \
						   else protein_sel
		self.fit_sel = fit_sel

	def CalcGyrationRadius(self, system, step=1, verbose=False):
		"""
		Calculate radius of gyration along the trajectory

		system - MDAnalysis Universe, trajectory for analysis
		step - int, step for trajectory analysis
		verbose - bool, report trajectory analysis progress

		Also requires:
		self.protein_sel - one or multiple atom selections for Rg calculation

		Returns
		np.array(floats)
		"""

		# make selection(s)
		selections = [
			system.select_atoms(sel) for sel in self.protein_sel
		]

		# iterate over trj and calculate R
		Rg_by_frame = []
		iterator = system.trajectory[::step]	
		if verbose:
			print('Calculating gyration radius...')
		for ts in (tqdm(iterator) if verbose else iterator):
			Rg_by_frame.append(
				[sel.radius_of_gyration() for sel in selections]
			)

		return np.array(Rg_by_frame)

	def CalcRMSD(self, system, step=1, verbose=False):
		"""
		Calculate RMSD for a group 

		system - MDAnalysis Universe, trajectory for analysis
		step - int, step for trajectory analysis
		verbose - bool, report trajectory analysis progress

		Also requires:
		self.protein_sel - one or multiple atom selections for RMSD calculation
		self.fit_sel - atom selection for superimposition

		Returns
		np.array(floats)
		"""

		ref = system.copy()

		rms = mda.analysis.rms.RMSD(
			system,
			ref,
			select = self.fit_sel,
			center = True,
			groupselections = self.protein_sel
		) 

		if verbose:
			print('Calculating RMSD...')
		
		rms.run(
			step = step,
			verbose = True
		)

		return rms.results.rmsd[:, -len(self.protein_sel):]