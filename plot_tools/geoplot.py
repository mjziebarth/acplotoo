# -*- coding: utf-8 -*-
#
# Plot tools Geoplot submodule module init file.
# Copyright (C) 2018 Malte Ziebarth
# 
# This software is distributed under the MIT license.
# See the LICENSE file in this repository.



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
	
	
	def _adjust_axes():
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
			
	
	def _generate_axes_ticks():
		# Generate the ticks on the axes.
