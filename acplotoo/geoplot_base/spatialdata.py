# -*- coding: utf-8 -*-
#
# Acplotoo SpatialData class file.
# Copyright (C) 2018 Malte Ziebarth
# 
# This software is distributed under the MIT license.
# See the LICENSE file in this repository.

import numpy as np


class SpatialData:
	"""
	A class representing a spatial data set. It is used to
	ensure that data is given in a format suitable for plotting.
	"""

	def __init__(self, lon=None, lat=None, data=None, x=None, y=None):
		"""
		Initialize a SpatialData instance.
		"""

		# 1) Exactly one of the pairs (lon,lat) and (x,y) has to be given.
		#    Which one determines the coordinate type.
		if (lon is None) == (x is None) or (lon is None) != (lat is None) \
		or (x is None) != (y is None):
			raise ValueError("Exactly one of the pairs (lon,lat) or (x,y) have "
			                 "to be given!")

		self._coordinate_type = "geographic" if lon is not None else "projected"
		if self._coordinate_type == "geographic":
			x0 = lon
			x1 = lat
		else:
			x0 = x
			x1 = y

		# 2) Ensure we can use all the numpy methods:
		x0 = np.atleast_1d(x0)
		x1 = np.atleast_1d(x1)
		data = np.atleast_1d(data)

		# 3) Classify input coordinate type.
		#    This is where most of the magic happens cleaning the input data
		#    and which allows us to easily use the data later on, hurray!
		#    Be sure to do a good job here.
		self._indexing = None
		self._ndim_x0 = x0.squeeze().ndim
		self._ndim_x1 = x1.squeeze().ndim
		self._ndim_data = data.squeeze().ndim
		if ndim_data < max(ndim_x0,ndim_x1):
			raise TypeError("Data dimensions must be at least that of the coordinate "
			                "dimensions!")
		if self._ndim_x0 > 2 or self._ndim_x1 > 2 or self._ndim_data > 3:
			raise TypeError("Support only ndim <= 2 for coordinate arrays and ndim <= 3 "
			                "for data arrays!")
		if self._ndim_data == 3:
			self._data_point_size = data.shape[2]
		else:
			self._data_point_size = 1

		if x0.size != x1.size:
			# Both coordinate arrays are of different size.
			# This can only be accomplished by x0 and x1
			# describing the coordinates of a grid.
			if ndim_x0 != 1 or ndim_x1 != 1 or (ndim_data != 2 and ndim_data != 3):
				raise TypeError("If sizes of coordinate arrays differ, they must describe "
				                "one-dimensional axes of data grid!")
			if x0.size == data.shape[0]:
				self._indexing = 'ij'
			elif x1.size == data.shape[0]:
				self._indexing = 'xy'
			else:
				raise TypeError("Coordinates/data shape not compliant!")

			self._src_type = 'llg'
		else:
			# Both coordinate arrays are of same size.
			# This can either be accomplished by x0 and x1 describing
			# the coordinate axes of a grid, or the coordinates of each
			# data points. This is decided by comparing to the data shape.
			if self._data_point_size * x0.size == data.size:
				# Each data point is associated with one coordinate point.
				if ndim_x0 == 1:
					self._src_type = 'lll'
				else:
					self._src_type == 'ggg'
			else:
				# Grid mode again.
				self._src_type = 'llg'
				if ndim_x0 != 1 or ndim_x1 != 1 or len(x0.shape) != len(x1.shape):
					raise TypeError("Coordinate shape not compliant!")
				if len(x0.shape) == 1:
					self._indexing = 'ij'
				else:
					if x0.shape[0] == 1:
						self._indexing = 'xy'
					else:
						self._indexing = 'ij'

		# Now we have read all the information we might need from the shape.
		# Squeeze the arrays and see if everything is really consistent with
		# what we expect:
		x0 = x0.squeeze()
		x1 = x1.squeeze()
		data = data.squeeze()
		if self._src_type == 'llg':
			if self._ndim_data not in [2,3] or (self._ndim_data != 3 and
			                                    self._data_point_size != 1):
				raise TypeError("Data shape not compliant (dim=" + str(self._ndim_data)
				                + ", shape=" + str(data.shape)
				                + ", point_size=" + str(self._data_point_size)
				                + ", x0.size=" + str(x0.size)
				                + ", x1.size=" + str(x1.size) + ")" )
			if self._indexing == 'ij':
				if x0.size != data.shape[0] or x1.size != data.shape[1]:
					raise TypeError("Data shape not compliant (dim=" + str(self._ndim_data)
					                + ", shape=" + str(data.shape)
					                + ", point_size=" + str(self._data_point_size)
					                + ", x0.size=" + str(x0.size)
					                + ", x1.size=" + str(x1.size) + ")" )
			elif self._indexing == 'xy':
				if x0.size != data.shape[1] or x1.size != data.shape[0]:
					raise TypeError("Data shape not compliant (dim=" + str(self._ndim_data)
					                + ", shape=" + str(data.shape)
					                + ", point_size=" + str(self._data_point_size)
					                + ", x0.size=" + str(x0.size)
					                + ", x1.size=" + str(x1.size) + ")" )
		elif self._src_type in ['lll','ggg']:
			if not np.array_equal(x0.shape, data.shape[:self._ndim_x0]) \
			or not np.array_equal(x1.shape, x0.shape):
				raise TypeError("When indexing each data point, the "
				                "coordinates have to be of same shape "
				                "as the data!")

		# For grid-grid-grid source type, determine the indexing order
		# and ascertain that the coordinates are according to grid definition:
		if self._src_type == 'ggg':
			if x0[0,0] != x0[1,0]:
				self._indexing = 'ij'
				# Test:
				assert np.all(np.diff(x0[:,0] > 0))
				assert np.all(np.diff(x0,axis=1) == 0)
				assert np.all(np.diff(x1[0,:] > 0))
				assert np.all(np.diff(x1,axis=0) == 0)
			else:
				self._indexing = 'xy'
				# Test:
				assert np.all(np.diff(x0[0,:] > 0))
				assert np.all(np.diff(x0,axis=0) == 0)
				assert np.all(np.diff(x1[:,0] > 0))
				assert np.all(np.diff(x1,axis=1) == 0)
