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
		self._cbar_label = None


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
		return self._routine == 'scatter' or self._routine == 'imshow'


	def cbar_label(self):
		"""
		Return the label of the color map data dimension.
		"""
		return self._cbar_label


	def cbar_handle(self):
		"""
		Return the handle useable for color bar plotting.
		"""
		if not self.has_cmap():
			raise RuntimeError("Handle does not contain color bar information!")

		if self._routine == 'imshow':
			return self._h[1]

		return self._h




class JointHandle(Handle):
	"""
	This class represents a handle of a constituent of
	a Geoplot made up of different constituents.

	It can be used in the colorbar routine.
	"""
	def __init__(self, routine, handles):
		self._routine = routine
		self._done = False
		self._args = None
		self._kwargs = None
		self._h = None
		self._subhandles = handles


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
		return [h() for h in self._subhandles]


	def subhandles(self):
		return self._subhandles


	def has_cmap(self):
		"""
		Return whether this handle has an associated color map,
		i.e. can be used to plot a
		"""
		return any(h.has_cmap() for h in self._subhandles)


	def cbar_label(self):
		"""
		Return the label of the color map data dimension.
		"""
		for h in self._subhandles:
			if h.has_cmap():
				return h.cbar_label()

		raise RuntimeError("Handle does not contain color bar information!")


	def cbar_handle(self):
		"""
		Return the handle useable for color bar plotting.
		"""
		for h in self._subhandles:
			if h.has_cmap():
				return h.cbar_handle()

		raise RuntimeError("Handle does not contain color bar information!")
