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
	Extract properties of a simulation system
	from edr files and/or trajectories
	and analyse convergence
	"""

	def __init__(
		self,
		edr=None,
		trj=None,
		topol=None,
		funcs=None,
		func_names=None,
		center_group='protein',
		rot_trans_group='protein'
	):
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

		transformation arguments
		center_group - str, atom selection to center in the box
		rot_trans_group - str, atom selection to fit alonmg the trj
		"""

		# energy and trajectory files for analysis
		self.edr = self._check_type(edr)
		self.trj = self._check_type(trj)
		self.topol = self._check_type(topol)

		# functions to be applied along the trajectory
		self.funcs = funcs
		self.func_names = self._check_type(func_names)

		# transformation arguments
		self.center_group = center_group
		self.rot_trans_group = rot_trans_group
		
		# pandas DataFrame with extracted data
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
				raise Exception("input should be str or list(str)")
				
		# if neither then raise exception
		else:
			raise Exception("input should be str or list(str)")
	
		return checked_var

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

		if self.rot_trans_group is not None:
			fit_group = system.select_atoms(self.rot_trans_group)
		
		# add centering
		if len(center_group.atoms) > 0:
			workflow.extend([
				transformations.center_in_box(center_group),
				transformations.wrap(ag, compound='residues')
			])	

		# add rotation and translation
		if len(fit_group.atoms) > 0:
			workflow.extend([
				transformations.fit_rot_trans(fit_group, fit_group)
			])

		return workflow

	def _read_edrs(self, tu, sequential):
		"""
		Read data from edr files and 
		append it to self.data

		tu - str, time units option, ns or ps
		
		sequential - bool, if True then supplied files
		are considered sequential steps and `Time` 
		is adjusted accordingly, deffault True

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

		# data from trajectories 
		trj_dat = {
			'Time': [ts.time for ts in system.trajectory[::step]]
		}

		# return default property names if none supplied
		if self.func_names is None:
			self.func_names = [
				f"Prop{i}" for i in range(1, len(self.funcs) + 1)
			]

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

			if verbose:
				print(f"Loading and transforming {trj}...")

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
				raise Exception("`func` argument is required for trajectory analysis")

			self._analyse_trjs(
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

	def __init__(self, lipid_sel=None, tail_sel=None, leaflet=0):
		"""
		Class atributes hold parameters to be used in methods

		lipid_sel - str, MDAnalysis selection for lipid group for analysis
		
		tail_sel - str, MDAnalysis selection for lipid tails
		that will be used by CalcOrderParameter

		leaflet - int, leaflets to use for property averaging 
		Valied for CalcAreaPerLipid
		-1 - lower
		1  - upper
		0  - both
		"""
		self.lipid_sel = lipid_sel
		self.tail_sel = tail_sel
		self.leaflet = leaflet

	def CalcAreaPerLipid(self, system, step=1, verbose=False):
		"""
		Calculate average area per lipid

		system - MDAnalysis Universe, trajectory for analysis
		step - int, step for trajectory analysis
		verbose - bool, report trajectory analysis progress

		Also requires:
		self.lipid_sel - atom selection for lipids in the bilayer.
		Atoms used to identify leaflets.
		These atoms will also be used to perform the Voronoi tessellation.
		self.leaflet - leaflet(s) to perform averaging for 

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
			verbose = verbose
		)
		
		# compute apl
		apl = lpp.analysis.area_per_lipid.AreaPerLipid(
			universe = system,
			lipid_sel = self.lipid_sel,
			leaflets = leaflets.leaflets
		)
		apl.run(
			step = step,
			verbose = verbose
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
	
	def CalcBilayerThickness(self, system, step=1, verbose=False):
		"""
		Calculate average bilayer thickness over trajectory 

		system - MDAnalysis Universe, trajectory for analysis
		step - int, step for trajectory analysis
		verbose - bool, report trajectory analysis progress
		
		Also requires:
		self.lipid_sel - atom selection for lipids in the bilayer.
		Atoms used to identify leaflets.
		These atoms will also be used to define thickness.

		Returns
		list(floats)
		"""

		# binning 
		n_bins_leaflets  = int(system.dimensions[0] // 10)
		n_bins_thickness = int(system.dimensions[0] // 20)
		
		# assign leaflets
		leaflets = lpp.leaflets.assign_leaflets.AssignLeaflets(
		  universe = system,
		  lipid_sel = self.lipid_sel,
		  n_bins = n_bins_leaflets 
		)
		
		leaflets.run(
			step = step,
			verbose = verbose
		)

		# compute thickness
		memb_thickness = lpp.analysis.MembThickness(
			universe = system,
		 	leaflets = leaflets.leaflets,
			lipid_sel = self.lipid_sel,
			n_bins = n_bins_thickness
		)

		memb_thickness.run(
			step = step,
			verbose = verbose
		)

		return memb_thickness.memb_thickness

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
		list(floats)
		"""

		# number of selections to work with?
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

			scc_sn2.run(
				step = step,
				verbose = verbose
			)

			scc_av = lpp.analysis.order_parameter.SCC.weighted_average(scc_sn1, scc_sn2)
			mean_order_parameters = np.nanmean(scc_av.SCC, axis=0)
		
		else:
			mean_order_parameters = np.nanmean(scc_sn1.SCC, axis=0)

		return mean_order_parameters

class ProteinPropertyCalculator:
	"""
	A class of with methods to calculate
	key protein properties from the trajectory
	"""

	def __init__(self, protein_sel, fit_sel):
		"""
		Class atributes hold parameters to be used in methods

		protein_sel - str, MDAnalysis atom selection 
		for protein group for analysis
		
		fit_sel - str, atom selection for superimposition in CalcRMSD
		"""
		self.protein_sel = protein_sel
		self.fit_sel = fit_sel

	def CalcGyrationRadius(self, system, step=1, verbose=False):
		"""
		Calculate radius of gyration along the trajectory

		system - MDAnalysis Universe, trajectory for analysis
		step - int, step for trajectory analysis
		verbose - bool, report trajectory analysis progress

		Also requires:
		self.protein_sel - atom selection for Rg calculation

		Returns
		list(floats)
		"""

		# make selection 
		prot_selection = system.select_atoms(self.protein_sel)

		# iterate over trj and calculate R
		Rg_by_frame = []
		iterator = system.trajectory[::step]	
		for ts in (tqdm(iterator) if verbose else iterator):
			Rg_by_frame.append(
				prot_selection.radius_of_gyration()			
			)

		return Rg_by_frame

	def CalcRMSD(self, system, step=1, verbose=False):
		"""
		Calculate RMSD for a group 

		system - MDAnalysis Universe, trajectory for analysis
		step - int, step for trajectory analysis
		verbose - bool, report trajectory analysis progress

		Also requires:
		self.protein_sel - atom selection for RMSD calculation
		self.fit_sel - atom selection for superimposition

		Returns
		list(floats)
		"""

		ref = system.copy()

		rms = mda.analysis.rms.RMSD(
			system,
			ref,
			select = self.fit_sel,
			center = True,
			groupselections = [self.protein_sel]
		) 

		rms.run(
			step = step,
			verbose = True
		)

		return rms.rmsd[:, -1]