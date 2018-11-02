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


	def inverse(x, y):
		"""
		Inverse projection: Project x,y coordinates
		back to longitude/latitude.
		"""
		return self._inverse(x,y)


	def generate_ticks(xlim, ylim):
		"""
		Generate ticks. Can be adjuste to special projections'
		needs by 
		"""
		if hasattr(self, _generate_ticks):
			# May become useful for global projections that
			# have funny (read: non-rectangular) outlines.
			self._generate_ticks(xlim,ylim)
		else:
			# Call the default backend (it uses optimization
			# to find bounds for a general projection)
			
			# TODO!
			_generate_ticks(self, xlim, ylim)
