# -*- coding: utf-8 -*-
#
# Plot tools projection submodule module init file.
# Copyright (C) 2018 Malte Ziebarth
# 
# This software is distributed under the MIT license.
# See the LICENSE file in this repository.

from .backend import _generate_ticks, _maximum_geographic_extents, _unit_vector

class Projection():
	"""
	Projection base class, i.e. interface definition.
	"""

	def __init__(self):
		pass


	def project(self, lon, lat):
		"""
		Project longitude/latitude coordinates to
		two-dimensional Euclidean space.
		"""
		return self._project(lon, lat)


	def inverse(self, x, y):
		"""
		Inverse projection: Project x,y coordinates
		back to longitude/latitude.
		"""
		return self._inverse(x,y)


	def unit_vector_east(self, lon=None, lat=None, x=None, y=None):
		"""
		Calculate the unit vector in longitude direction
		"""
		if (lon is None) != (lat is None) or (x is None != y is None) \
		   or (lon is None == x is None):
			raise ValueError("Either 'lon' and 'lat' or 'x' and 'y' have to be "
			                 "given!")

		if lon is None:
			lon,lat = self._projection.inverse(x,y)

		if hasattr(self, "_unit_vector_east"):
			return self._unit_vector_east(lon=lon, lat=lat, x=x, y=y)

		return _unit_vector(lon, lat, self, "east")


	def unit_vector_north(self, lon=None, lat=None, x=None, y=None):
		"""
		Calculate the unit vector in latitude direction
		"""
		if (lon is None) != (lat is None) or (x is None != y is None) \
		   or (lon is None == x is None):
			raise ValueError("Either 'lon' and 'lat' or 'x' and 'y' have to be "
			                 "given!")

		if lon is None:
			lon,lat = self._projection.inverse(x,y)

		if hasattr(self, "_unit_vector_north"):
			return self._unit_vector_north(lon, lat)

		return _unit_vector(lon, lat, self, "north")


	def is_global(self):
		"""
		Return true if the projection is global, i.e.
		always spans the whole globe.
		"""
		return self._is_global()

	def maximum_geographic_extents(self, xlim, ylim):
		"""
		Calculate the geographic extents of given x/y
		limits:
		"""
		if hasattr(self, '_max_geographic_extents'):
			return self._max_geographic_extents(xlim,ylim)

		# If not, use default backend that uses numerical optimization:
		return _maximum_geographic_extents(self, xlim, ylim)

	def identifier(self):
		"""
		A hashable identifier, a list L consisting of:
		L[0] : name of the projection (string)
		L[1] : *args to recreate an object of that class.
		"""
		return self._identifier()

	def generate_ticks(self, xlim, ylim, tick_delta_degree):
		"""
		Generate ticks. Can be adjuste to special projections'
		needs by supplying the _generate_ticks method.

		Required arguments:
		   xlim, ylim :       2d-array-like objects that contain
		                      the minimum coordinate at index 0
		                      and maximum coordinate at index 1.
		   tick_delta_degree: The desired spacing of ticks in
		                      degrees.

		Returns:
		   A dictionary of tick positions. The dictionary
		   contains entries "top", "bot", "left", and "right",
		   each of which contains a tuple of lists.
		   The first list of each tuple contains the longitude
		   ticks of that border, the second tuple the latitude
		   ticks.
		"""
		if hasattr(self, '_generate_ticks'):
			# May become useful for global projections that
			# have funny (read: non-rectangular) outlines.
			return self._generate_ticks(xlim,ylim)
		else:
			# Call the default backend (it uses optimization
			# to find bounds for a general projection)
			return _generate_ticks(self, xlim, ylim, tick_delta_degree)
