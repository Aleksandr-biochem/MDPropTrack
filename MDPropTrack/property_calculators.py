import numpy as np
from tqdm import tqdm
import MDAnalysis as mda
import lipyphilic as lpp

class BaseCalculator:
	"""
	Base class to construct property calculators
	that are applied to trajectories

	Class atributes hold parameters to be used in _single_frame method

	Attributes
	----------
	"""

	def __init__(self):
		pass

	def _single_frame(self, system):
		"""
		This method should perform some calculations on one frame
		and return a property value

		This method can use class attributes as parameters

		Returns
		----------
		float, list(floats) or np.array(floats)
		"""
		pass

	def Calc(self, system, step=1, verbose=False):
		"""
		Run calculation by applying self._single_frame()
		along the trajectory

		Parameters
		----------
		
		system: MDAnalysis Universe
			universe for analysis

		step: int
			trajectory analysis step

		verbose: bool
			verbose progress

		Returns
		----------
		np.array(floats)
		"""

		# to store results
		results = []

		if verbose:
			print('Running property calculator...')

		# iterate over trj and calculate property
		iterator = system.trajectory[::step]
		for ts in (tqdm(iterator) if verbose else iterator):
			results.append(
				self._single_frame(system)
			)

		return np.array(results)

class GyrationRadiusCalculator(BaseCalculator):
	"""
	Class describing calculation of Rg along the trajectory

	Class atributes hold parameters to be used in _single_frame method

	Attributes
	----------
	protein_sel: str or list(str)
		one of several selections for analysis

	"""
	
	def __init__(self, protein_sel=None):

		# protein_sel is a paramether for our calcvulator defining
		# the atoms groups for Rg calculation
		self.protein_sel = [protein_sel] if isinstance(protein_sel, str) \
						   else protein_sel

		# we can also define additional attributes to assist calculation
		self.at_groups = None
	
	# we can create additional methods to assist calculation
	# this method allows to make atom selections once
	# to use them in every frame instead of calling select_atoms() every frame
	def _make_selections(self, system):
		"""
		Make atom selection for Rg calculation
		They are stored in self.at_groups

		Parameters
		----------
		
		system: MDAnalysis Universe
			universe for analysis

		Returns
		----------
		self
		"""

		self.at_groups = [
			system.select_atoms(sel) for sel in self.protein_sel
		]

		return self

	def _single_frame(self, system):
		"""
		This method should perform some calculations on one frame
		and return a property value

		This method can use class attributes as parameters

		Returns
		----------
		float, list(floats) or np.array(floats)
		"""

		# make atom selections if None
		if self.at_groups is None:
			self._make_selections(system)

		# calculate Rg for every atom group
		property_val = [
			sel.radius_of_gyration() for sel in self.at_groups
		]

		return property_val

class LipidPropertyCalculator:
	"""
	A calculator class of with methods to compute
	key lipid properties from the trajectory
	
	Class attributes hold parameters to be used in methods

	Attributes
	----------
	lipid_sel: str
		lipid group selection for leaflet identification
		and membrane thickness calculation (usually phosphate)
	
	apl_sel: str,
		atom selection to perform Voronoi Tesselation on

	tail_sel: str
		lipid tail selection that will be used for order parameter calculation
	
	calculate: str or list(str)
		keywords of properties to calculate
		choices: ['apl', 'thickness', 'order_param'], default all 3

	filter_lipid: str or list(str)
		one or multiple atom selections to filter lipids in APL and order parameter calculation

	leaflet_to_average: int
		leaflets to use for area per lipid averaging
		-1 - lower
		1  - upper
		0  - both

	bin_len_leaflets: float
		bin width for leaflet identification, default 15
	
	bin_len_thickness: float
		bin width for membrane thickness calculation, default 20

	"""

	def __init__(
			self,
			lipid_sel=None,
			apl_sel=None,
			tail_sel=None,
			calculate=['apl', 'thickness', 'order_param'],
			filter_lipid=['all'],
			leaflet_to_average=0,
			bin_len_leaflets=15,
			bin_len_thickness=20
		):
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
		Assigns self.leaflets
		
		Parameters
		----------
		
		system: MDAnalysis Universe
			universe for analysis

		step: int
			trajectory analysis step

		verbose: bool
			verbose progress

		Returns
		----------
		self
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
		
		Parameters
		----------
		
		system: MDAnalysis Universe
			universe for analysis

		main_sel: str
			main atom selection for filtering
			
		filter_sel: str
			atom selection to combine with main selection

		Returns
		----------
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

		Parameters
		----------
		
		system: MDAnalysis Universe
			universe for analysis

		step: int
			trajectory analysis step

		verbose: bool
			verbose progress

		Returns
		----------
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

		Parameters
		----------
		
		system: MDAnalysis Universe
			universe for analysis

		step: int
			trajectory analysis step

		verbose: bool
			verbose progress

		Also requires:
		self.lipid_sel: str
			atom selection for lipids in the bilayer.
			These atoms will also be used to perform the Voronoi tessellation
		
		self.leaflet_to_average: int
			leaflets to use for area per lipid averaging
			-1 - lower
			1  - upper
			0  - both

		Returns
		----------
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

		Parameters
		----------
		
		system: MDAnalysis Universe
			universe for analysis

		step: int
			trajectory analysis step

		verbose: bool
			verbose progress
		
		Also requires:
		self.lipid_sel: str
			atom selection for lipids in the bilayer.
			Atoms used to identify leaflets.
			These atoms will also be used to define thickness.

		Returns
		----------
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

		Parameters
		----------
		
		system: MDAnalysis Universe
			universe for analysis

		step: int
			trajectory analysis step

		verbose: bool
			verbose progress

		Also requires:
		self.tail_sel: str, list(str)
			atom selection(s) for lipid tails for order parameter calculation

		Returns
		----------
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

class RMSDCalculator:
	"""
	Class to calculare RMSD

	Class atributes hold parameters to be used in methods

	Attributes
	----------
	protein_sel: str or list(str)
		one of several selections for analysis
		
	fit_sel: str
		atom selection for least squeare fit
	"""

	def __init__(
		self,
		protein_sel=None,
		fit_sel=None
	):
		self.protein_sel = [protein_sel] if isinstance(protein_sel, str) \
						   else protein_sel
		self.fit_sel = fit_sel

	def Calc(self, system, step=1, verbose=False):
		"""
		calculate RMSD along the trajectory 

		Parameters
		----------
		
		system: MDAnalysis Universe
			universe for analysis

		step: int
			trajectory analysis step

		verbose: bool
			verbose progress

		Also requires:
		self.protein_sel: str
			one or multiple atom selections for RMSD calculation
		self.fit_sel: str
			atom selection for least squeare fit

		Returns
		----------
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
			verbose = verbose
		)

		return rms.results.rmsd[:, -len(self.protein_sel):]

	