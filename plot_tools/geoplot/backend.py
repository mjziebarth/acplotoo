# -*- coding: utf-8 -*-
#
# Plot tools geoplot backend file.
# Copyright (C) 2018 Malte Ziebarth
# 
# This software is distributed under the MIT license.
# See the LICENSE file in this repository.

import numpy as np

def _generate_axes_boxes(tick_dict, xlim, ylim):
	"""
	Generates axes boxes frame in a coordinate system
	[0,1]x[0,1].
	"""
	# First, for each axis generate a sorted list of ticks:
	axes = ["bot","top","left","right"]
	ticks_sorted = []
	for i in range(4):
		ticks_unsorted = tick_dict[axes[i]]
		n0 = ticks_unsorted[0].size
		n1 = ticks_unsorted[1].size
		n_ticks = n0 + n1
		# In tick_array, we save all ticks and mark them by
		# coordinate:
		tick_array = np.zeros((n_ticks,2))
		tick_array[0:n0,0] = ticks_unsorted[0]
		tick_array[n0:,0] = tick_unsorted[1]
		tick_array[n0:,1] = 1
		# Sort:
		order = np.argsort(tick_array[:,0])
		tick_array = tick_array[order,:]
		# Save:
		ticks_sorted[i] = tick_array
	
	# Now generate the boxes:
	for i in range(4):
		# Convert the tick positions to interval [0,1]:
		lim = xlim if i < 2 else ylim
		normalized = (ticks_sorted[i] - lim[0]) / (lim[1]-lim[0])
		x = np.concantenate(np.zeros(1),normalized,np.ones(1))
	
	
	
	var_index = 
