# -*- coding: utf-8 -*-
#
# Plot tools Geoplot submodule module init file.
# Copyright (C) 2018 Malte Ziebarth
# 
# This software is distributed under the MIT license.
# See the LICENSE file in this repository.

from .geoplot.backend import _generate_axes_boxes

class Geoplot:
	
	
	def __init__(self, ax, projection):
		"""
		Init method.
		
		Required arguments:
		   ax         :
		   projection :
		"""
		# TODO : Checks!
		self._ax = ax
		self._projection = projection
		
		
		# Setup internal data:
		self._data_xlim = None
		self._data_ylim = None
		self._xlim = None
		self._ylim = None
		self._ticks = None
		
		# Setup configuration:
		self._box_axes = True
	
	
	def _adjust_axes(self):
		"""
		This method is called whenever the data situation
		has changed and we need to update the axes.
		"""
		
		# 1) See whether we have to adjust ticks:
		if self._data_xlim is not None:
			if self._data_xlim[0] != self._xlim[0] or
			   self._data_xlim[1] != self._xlim[1] or
			   self._data_ylim[0] != self._ylim[0] or
			   self._data_ylim[1] != self._ylim[1]:
			
			self._xlim = self._data_xlim
			self._ylim = self._data_ylim
			
			# Compute ticks:
			ticks = self._projection.generate_ticks(self._xlim, self._ylim)

			# Plot ticks:
			self._plot_axes(tick_dict)

	def _plot_axes(self, tick_dict):
		"""
		This method draws the axes artists.
		"""
		if self._box_axes:
			# Draw box axes as polygons. In the backend, we generate
			# the poly collection in a coordinate system x:[0,1]
			# and y:[0,1].
			axes_boxes = _generate_axes_boxes(tick_dict, self._xlim, self._ylim)
