# Support color bars (and maybe more later on).

from matplotlib.collections import PathCollection
from matplotlib.lines import Line2D

class Handle:
	"""
	This class represents a handle of a constituent of
	a Geoplot.

	It can be used in the colorbar routine.
	"""
	def __init__(self, routine, args, kwargs):
		self._routine = routine
		self._done = False
		self._args = args
		self._kwargs = kwargs
		self._h = None


	def routine(self):
		return self._routine


	def register(self, h):
		"""
		Register a matplotlib handle.
		"""
		self._h = h


	def __call__(self):
		"""
		Return the matplotlib handle.
		"""
		return self._h


	def has_cmap(self):
		"""
		Return whether this handle has an associated color map,
		i.e. can be used to plot a
		"""
		return self._type == 'scatter' or self._type == 'imshow'
