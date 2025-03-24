import panedr
import numpy as np
import pandas as pd
import seaborn as sns
import MDAnalysis as mda
import matplotlib.pyplot as plt

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

	tags: list(str), dict(str: list(str)
		tags describing each simulation
		list(str) is assumend to be a list of simulation names

	funcs: list(functions)
		list of functions to apply along the trajectory
		Functions are expected to adhere to a specific input-output structure
		
	func_names: list(str)
		names for the properties computed by funcs

	data: pd.DataFRame
		data extracted from simulations
	
	tag_names: list(str)
		list of tags in self.data

	_tag_combinations: list(list(str))
		list of tag value combinations in self.data

	properties: list(str)
		list of properties in self.data

	tau_data: pd.DataFrame
		autocorrelation time (tau) vs simulation length for each property

	"""

	def __init__(
		self,
		simulations=None,
		tags=None,
		funcs=None,
		func_names=None
	):

		# energy files, trajectory files and Universe instances for analysis
		self.simulations = [] if (simulations is None) else self._check_input(simulations)
		
		# autofill names for simulation steps
		if tags is None:
			self.tags = [
				{'name': f"Simulation_{i}"} for i in range(1, len(self.simulations) + 1)
			]
			self.tag_names  = list(self.tags[0].keys()) if len(self.tags) > 0 else []
		# tags from input
		else:
			self.tags = self._transform_tags(tags)
			self.tag_names  = list(self.tags[0].keys())

		# functions to be applied along the trajectory/Universe
		self.funcs = funcs
		self.func_names = func_names

		# pandas DataFrame with extracted data
		self.data = None
		self.properties = None
		self._tag_combinations = []

		# pandas dataframe with convergence data
		self.tau_data = None

		# internal colour-blind friendly cmap
		self._custom_palette = [
			'#648fff', '#dc267f', 
			'#785ef0', '#fe6100',
			'#ffb000', '#000000'
		]
	
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
		list(str, tuple(str, str), mda.Universe) or None
		"""

		# if Noen return None
		if var is None:
			var = var

		# if list check taht every element is of supported type
		elif isinstance(var, list):
			for v in var:
				if self._check_var_type(v):
					continue
				else: 
					raise Exception("Couldn't process `simulations` arg. Check types")
		
		# if one element, check and convert to list
		elif self._check_var_type(var):
			var = [var]

		# if neither then raise exception
		else:
			raise Exception("Couldn't process `simulations` arg. Check types")
	
		return var

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
			tags = [ {'name': str(t)} for t in tags]

		# if dictionary, transform it
		# make sure that all tag values are str
		elif isinstance(tags, dict):
			tags = [
				 {t: str(v) for t, v in zip(tags, val_comb)} \
				 for val_comb in zip(*tags.values())
			]

		else:
			raise Exception("Couldn't process `tags`. Check format")

		return tags

	def load_data(self, data, tag_columns=['name'], read_csv_kwargs={'index_col': None}):
		"""
		Load DataFrame from file into self.data
		Populate self.tag_names, self.properties and self._tag_combinations

		Parameters
		----------
		
		data: str or pd.DataFrame
			path to data file or or pd.DataFrame to use as data

		read_csv_kwargs: dict
			kwargs for pd.read_csv; default {'index_col': None}

		tag_columns: list(str)
			name of the columns that are tags,
			Default ['name']

		Returns
		----------
		self
		"""

		# read data file
		if isinstance(data, str):
			self.data = pd.read_csv(data, **read_csv_kwargs)
		elif isinstance(data, pd.DataFrame):
			self.data = data
		else:
			raise Exception("data should be file path or pd.DataFrame")


		# separate tags and properties
		self.tag_names = []
		self.properties = []
		for col in self.data.columns:
			if col != 'Time':
				if col in tag_columns:
					self.tag_names.append(col)
				else:
					self.properties.append(col)
			else:
				continue

		# make sure that tag columns are string
		for tag in self.tag_names:
			self.data[tag] = self.data[tag].astype('string')

		# get tag combinations
		self._tag_combinations = list(
			self.data.groupby(self.tag_names).count().index
		)

		return self

	def _append_data(self, df):
		"""
		Append DataFrame to self.data with proper tag merger
		
		Parameters
		----------
		
		trj: pd.DataFrame
			dta to append

		Returns
		----------
		self
		"""

		# tag combination in this df
		tag_combination = set(np.unique(
			df[self.tag_names].values.astype(str),
			axis=0
		)[0])

		# if self.data is empty
		if self.data is None:
			self.data = df

		else:
			
			# merge if there is a tag match
			if tag_combination in self._tag_combinations:
				self.data = pd.merge(
					self.data,
					df,
					on = ['Time'] + self.tag_names,
					how = 'outer'
				)

			# concat if there is no tag match
			else:
				self.data = pd.concat([self.data, df]) \
							.reset_index(drop=True)

		# add tag combination of the new DataFrame piece
		self._tag_combinations.append(tag_combination)

		return self

	def _read_edr(self, edr, tags, tu, sequential):
		"""
		Read data from edr file and append it to self.data

		Parameters
		----------
	
		file: str
			path to edr file

		tags: dict(str: str)
			tag comuns to be added to extarcted data

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
			df[tag] = df[tag].astype('string')
		
		# append to self.data
		self._append_data(df)
			
		return self

	def _apply_funcs_trj(self, system, step, verbose):
		"""
		Apply all functions to the trajectory
		Returns calculated properties as pd.DataFrame
		
		Parameters
		----------

		system: mda.Universe
			trajectory for analysis
		
		step: int
			step for analysis

		verbose: bool
			verbose the trj analysis process

		Returns
		----------
		trj_dat - pd.DataFrame with Time and calculated properties
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
		
	def _analyse_trj(self, trj, tags, tu, step, verbose, sequential):
		"""
		Calculate properties along the trajectory
		and append data to self.data
		
		Parameters
		----------
		
		trj: mda.Universe
			system fopr analysis

		tags: dict(str: str)
			tag comuns to be added to extarcted data

		tu: str
			time units option, ns or ps

		step: int
			step for analysis

		verbose: bool
			verbose the trj analysis process
		
		sequential:
			bool, if True then `Time` is adjusted when appending to self.data
			to indicate sequential simulations

		Returns
		----------
		self	
		"""

		if verbose:
			trj_name = ', '.join(
				[f'{t}: {tags[t]}' for t in tags]
			)
			print(f"Analysing trj {trj_name}...")

		# calculate properties over the trajectory
		trj_dat = self._apply_funcs_trj(
			system = trj,
			step = step,
			verbose = verbose
		)
		
		# adjust time units
		if tu == 'ns':
			trj_dat['Time'] = trj_dat['Time'] / 1000

		# shift Time for sequential steps
		if sequential and (self.data is not None):
			trj_dat['Time'] += self.data.Time.iloc[-1]
		
		# add tag columns
		for tag in tags:
			trj_dat[tag] = tags[tag]
			trj_dat[tag] = trj_dat[tag].astype('string')

		# append to self.data
		self._append_data(trj_dat)

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
		
		# reset self.data
		self.data = None

		# analyse each simulation input
		for sim, tags in zip(self.simulations, self.tags):

			# single file is expected to be .edr
			if isinstance(sim, str):

				self._read_edr(
					edr = sim,
					tags = tags,
					tu = tu,
					sequential = sequential
				)

			# analyse trajectory
			else:

				# transform trj file into MDAnalysis Universe
				if isinstance(sim, tuple):
					trj = mda.Universe(sim[0], sim[1])
				else:
					trj = sim

				# analyse trajectory
				self._analyse_trj(
					trj=trj,
					tags=tags,
					tu = tu,
					step = step,
					sequential = sequential, 
					verbose = verbose
				)
		
		# collect property names 
		self.properties = []
		for col in self.data.columns:
			if (col != 'Time') and (col not in self.tag_names):
				self.properties.append(col)

		return self

	################################
	# The following group of methods for convergence assessment
	# is implemented with the code adapted from
	# https://emcee.readthedocs.io/en/stable/tutorials/autocorr/
	# Copyright 2012-2021, Dan Foreman-Mackey & contributors

	@classmethod
	def _get_next_pow_two(self, n):
		"""
		Return right nearest 
		int that is power of 2

		Parameters
		----------
		
		n: float
			number to find nearest power of 2 int for

		Returns
		----------
		int - nearest int that is power of 2
		"""
		i = 1
		while i < n:
			i = i << 1
		return i

	@classmethod
	def _get_autocorr_func_1d(self, y, norm=True):
		"""
		Estimate autocorrelation function from
		time series data using one-dimensional
		discrete Fourier Transform

		Parameters
		----------
		
		y: np.array(floats)
			time series data
		norm: bool
			normalise acf, default True

		Returns
		----------
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
	
	@classmethod
	def _auto_window(self, taus, c):
		"""
		Automated windowing procedure following Sokal (1989)
		
		Parameters
		----------
		taus: np.array(float)
			tau estimates

		c: float
			coefficient in tau estimation

		Return
		int, window number
		"""

		m = np.arange(len(taus)) < c * taus
		
		if np.any(m):
			window = np.argmin(m)
		else:
			window = len(taus) - 1
		
		return window

	@classmethod
	def _estimate_autocorr_tau(self, y, c=5.0):
		"""
		Estimate autocorrelation time tau
		from time-series

		Parameters
		----------
		
		y: np.array(floats)
			time series data

		c: float
			coefficient in tau estimation, default 0.5

		Returns
		----------
		tau - float
		"""

		# get acf estimate
		f = self._get_autocorr_func_1d(y)
		
		# compute and return tau estimate
		taus = 2.0 * np.cumsum(f) - 1.0
		window = self._auto_window(taus, c)
		
		return taus[window]

	def estimate_convergence(self):
		"""
		Estimate autocorrelation time (tau) vs simulation length
		to assess convergence

		Assigns self.tau_data based on self.data

		Returns
		----------
		self 
		"""

		if self.data is None:
			raise Exception('.data is empty')
		
		# reset self.tau_data
		self.tau_data = []

		# analyse convergence for every property
		for prop in self.properties:
			
			# and for for every tag combination 
			for tag_comb in self._tag_combinations:

				# subset data
				tag_filter = ' & '.join(
					[
						f"{tag} == '{tag_val}'" for tag, tag_val \
						in zip(self.tag_names, tag_comb)
					]
				)
				dat = self.data.query(tag_filter)[['Time', prop]].values

				# generate step points and
				# time points for tau estimates
				N = np.exp(
					np.linspace(np.log(100), np.log(dat.shape[0]), 10)
				).astype(int)
				ts = dat[N - 1, 0]

				# estimate tau from trj slices
				tau_data = np.array([
					self._estimate_autocorr_tau(dat[:n, 1]) for n in N
				])

				# generate DataFrame
				tau_data = pd.DataFrame.from_dict(
					data = {
						'Time': ts,
						'tau': tau_data
					}
				)
				tau_data['Property'] = prop
				for tag, tag_val in zip(self.tag_names, tag_comb):
					tau_data[tag] = tag_val

				# add to self.tau_data
				self.tau_data.append(tau_data)

		# merge all dataframes in self.tau_data
		self.tau_data = pd.concat(self.tau_data) \
						.reset_index(drop=True)

		return self

	##################################
	
	# from matplotlib.colors import LinearSegmentedColormap, to_rgba_array
	# @staticmethod
	# def hex_to_cmap(self, hex_colours):
	# 	"""
	# 	Produce a colormap from a list of discrete colors without interpolation
		
	# 	Parameters
	# 	----------

	# 	hex_colours: list(str)
	# 		list of hex colour codes

	# 	Returns
	# 	----------
	# 	colourmap
	# 	"""

	# 	# covert to rgb and reshape
	# 	clrs = to_rgba_array(hex_colours)
	# 	clrs = np.vstack([clrs[0], clrs, clrs[-1]])

	# 	colour_dict = {
	# 		prime_color : [
	# 			(i / (len(clrs) - 2.), clrs[i, j], clrs[i + 1, j]) for i in range(len(clrs) - 1)
	# 		] for j, prime_color in enumerate(['red','green','blue'])
	# 	}
		
	# 	return LinearSegmentedColormap('Custom_cmap', colour_dict)

	def _construct_multiplot(self, n_subplots, figure_kwargs):
		"""
		Construct subplot grid

		Parameters
		----------

		n_subplots: int
			number of subplots on the grid
		
		figure_kwargs: dict
			matplotlib figure kwargs

		Return
		----------
		fig, axs - mplt Figure and list(Axes)
		"""
	
		if n_subplots > 1:
			ncols = 2
			nrows = (n_subplots // 2) + (n_subplots % 2)
		else:
			ncols = 1
			nrows = 1

		internal_figure_kwargs = {
			'ncols': ncols,
			'nrows': nrows,
			'figsize': (10 * ncols, 5 * nrows)
		}

		# overwrite from figure_kwargs
		if figure_kwargs is not None:
			for key in figure_kwargs:
				internal_figure_kwargs[key] = figure_kwargs[key]

		fig, axs = plt.subplots(
			**internal_figure_kwargs
		)

		# transform to 1D array of Axes
		if isinstance(axs, np.ndarray):
			axs = axs.flatten()
		else:
			axs = np.array([axs])

		return fig, axs
	
	def _get_subplot_specs(self, properties_to_plot, subplot_by, query, plot_convergence):
		"""
		Define specifications for each subplot:
		- property to plot
		- data subset to work on
		
		Parameters
		----------
		
		properties_to_plot: str or list(str)
			columns from self.data to plot, default None

		subplot_by: str or list(str),
			self.data column(s) to use for subplot separation
			Default None, will subplot only by properties_to_plot

		query: str
			The query string to evaluate for pd.query().
			Used to subset a part of self.data for plotting

		plot_convergence: bool
			Plot convergence instread of time series

		Returns
		----------
		subplot_specs: dict(subplot_name: {'query': str, 'prop': str})
		"""

		# make sure that we have list(str)
		prop_list  = self._check_input(properties_to_plot)
		subplot_by = self._check_input(subplot_by)

		# for convergence plot we use properties for subplotting
		if plot_convergence:
			subplot_by = ['Property'] if subplot_by is None \
						 else ['Property'] + subplot_by

		# define subplot grouping from 'subplot_by' arg
		if subplot_by is None:
			subplot_queries = [(None, query)]
		else:
			
			subplot_queries = []

			# get unique combinations of tags values to query
			if plot_convergence:
				# we also need to filter Property column first
				tag_val_combs = list(
					self.tau_data[self.tau_data.Property.isin(prop_list)] \
					.groupby(subplot_by).count().index
				) 
			else:
				tag_val_combs = list(
					self.data.groupby(subplot_by).count().index
				)

			# this should be list of tuples
			if isinstance(tag_val_combs[0], str):
				tag_val_combs = [(v, ) for v in tag_val_combs]

			# define all queries
			for tag_val_comb in tag_val_combs:
				
				# define query name as tag value combination
				query_name = ' '.join(tag_val_comb)

				# construct query expression
				query_expr = ' & '.join(
					[f"{t} == '{v}'" for t, v in zip(subplot_by, tag_val_comb)]
				)

				# add additional query if requested
				if query is not None:
					query_expr= f"({query_expr}) & ({query})"

				subplot_queries.append(
					(query_name, query_expr)
				)

		# define specifications for each subplot
		subplot_specs = {}
		for prop in (['tau'] if plot_convergence else prop_list):
			for query_name, query_expr in subplot_queries:
				
				# define subplot title
				subplot_name = prop if query_name is None \
							   else f'{prop}, {query_name}'
				
				# define property toplot and query for self.data
				subplot_specs[subplot_name] = {
					'prop': prop,
					'query': query_expr
				}

		return subplot_specs

	def plot(
			self,
			properties_to_plot=None,
			plot_convergence=False,
			subplot_by=None,
			query=None,
			hue='name',
			x='Time',
			x_lab='Time, ns',
			palette=None,
			figure_kwargs=None, 
			style_kwargs={
				"style": "darkgrid",
				"rc": {"grid.color": ".6", "grid.linestyle": ":"}
			},
			sns_kwargs={'alpha': 0.7}
		):
		"""
		Plot properties from self.data using the 'Time' column as X-axis
		
		Parameters
		----------

		properties_to_plot: str or list(str)
			columns from self.data to plot, default None

		plot_convergence: bool
			Plot convergence instread of time series; default False
		
		subplot_by: str or list(str),
			self.data column(s) to use for subplot separation
			Default None, will subplot only by properties_to_plot
		
		query: str
			The query string to evaluate for pd.query().
			Used to subset a part of self.data for plotting. Default None

		hue: str or list(str),
			self.data column(s) to use for hue, default 'name'
		
		x: str,
			Column name to use as x-axis; default: 'Time'

		x_lab: str,
			x axis label, defalut 'Time, ns'

		palette: str, list(str), mplt cmap
			seaborn colourmap name, a list of colours or mplt cmap
		
		figure_kwargs: dict
			matplotlib figure kwargs
		
		style_kwargs: dict
			seaborn style kwargs

		sns_kwargs: dict
			seaborn lineplot kwargs
		
		Returns:
		----------
		fig, axs - mplt Figure and list(Axes)
		"""
		
		# do we have data to plot
		if self.data is None:
			raise Exception("self.data is None")
		# and for convergence plots
		if plot_convergence and (self.tau_data is None):
			self.estimate_convergence()

		# define subplot_specs
		subplot_specs = self._get_subplot_specs(
			properties_to_plot=properties_to_plot,
			subplot_by=subplot_by,
			query=query,
			plot_convergence=plot_convergence
		)

		# check if we need a new column for hue
		if isinstance(hue, list):
			hue_col = ' '.join(hue)
			self.data[hue_col] = self.data[hue].apply(
				lambda row: ' '.join([str(el) for el in row]), axis=1
			)
		else:
			hue_col = hue
		
		# define palette
		if palette is not None:
			sns_kwargs['palette'] = palette
		else:
			sns_kwargs['palette'] = self._custom_palette

		# use clear marker for convergence plot
		if plot_convergence and ('marker' not in sns_kwargs.keys()):
			sns_kwargs['marker'] = 'o'

		# use seaborn style for this multiplot
		with sns.axes_style(**style_kwargs):
			
			# define subplot grid
			fig, axs = self._construct_multiplot(
				n_subplots = len(subplot_specs),
				figure_kwargs = figure_kwargs
			)

			# populate each subplot
			for i, subplot_name in enumerate(subplot_specs):

				dat = 'tau_data' if plot_convergence else 'data'

				sns.lineplot(
					data = getattr(self, dat) if subplot_specs[subplot_name]['query'] is None \
						   else getattr(self, dat).query(subplot_specs[subplot_name]['query']),
					x = x,
					y = subplot_specs[subplot_name]['prop'],
					hue = hue_col,
					ax = axs[i],
					**sns_kwargs
				)

				# set title
				axs[i].set_title(
					subplot_name,
					fontweight='bold',
					fontsize=18,
					pad=10
				)

				# set y-label 
				axs[i].set_ylabel(
					ylabel= r"Autocorrelation time $\tau$" if plot_convergence \
							else subplot_specs[subplot_name]['prop'], 
					fontsize=15,
					labelpad=10
				)

				# set x-axis label 
				axs[i].set_xlabel(
					x_lab,
					fontsize=15,
					labelpad=10
				)
		
		# remove special hue column from data if it was constructed
		if hue_col != hue:
			self.data.drop(columns=[hue_col])

		plt.tight_layout()
		
		return fig, axs