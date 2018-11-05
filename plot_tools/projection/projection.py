# -*- coding: utf-8 -*-
#
# Plot tools projection submodule module init file.
# Copyright (C) 2018 Malte Ziebarth
# 
# This software is distributed under the MIT license.
# See the LICENSE file in this repository.

from .backend import _generate_ticks

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


	def unit_vector_east(self, lon,lat):
		"""
		Calculate the unit vector in longitude direction
		"""
		return self._unit_vector_east(lon,lat)


	def unit_vector_north(self, lon,lat):
		"""
		Calculate the unit vector in latitude direction
		"""
		return self._unit_vector_north(lon,lat)


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
			self._generate_ticks(xlim,ylim)
		else:
			# Call the default backend (it uses optimization
			# to find bounds for a general projection)
			
			# TODO!
			_generate_ticks(self, xlim, ylim, tick_delta_degree)
