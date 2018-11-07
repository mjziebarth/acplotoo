# -*- coding: utf-8 -*-
#
# Plot tools geoplot backend file.
# Copyright (C) 2018 Malte Ziebarth
# 
# This software is distributed under the MIT license.
# See the LICENSE file in this repository.

import numpy as np

from .rect import Rect
from matplotlib.patches import Rectangle

def _create_tick_arrays(tick_dict, which_ticks):
	"""
	Generate arrays of ticks to be used with _generate_axes stuff
	"""
	axes = ["bot","top","left","right"]
	tick_arrays = []
	for i in range(4):
		# First, for each axis generate a sorted list of ticks:
		ticks_unsorted = tick_dict[axes[i]]
		n0 = ticks_unsorted[0].size
		n1 = ticks_unsorted[1].size

		# Depending on the tick style, we have to continue differently:
		if which_ticks == 'both':
			J = [0,1]
			n_ticks = n0 + n1
			tick_array = np.zeros((n_ticks,2))
			tick_array[0:n0,0] = ticks_unsorted[0]
			tick_array[n0:,0] = ticks_unsorted[1]
			tick_array[n0:,1] = 1
		else:
			# Choose one:
			if which_ticks == 'lonlat':
				j = int(i/2)
			elif which_ticks == 'latlon':
				j = (int(i/2)+1) % 2
			elif which_ticks == 'significant':
				j = int(np.argmax([n0,n1]))
			J = [j]
			tick_array = np.zeros((ticks_unsorted[j].shape[0],2))
			tick_array[:,0] = ticks_unsorted[j]
			tick_array[:,1] = j

		# In tick_array, we save all ticks and mark them by
		# coordinate:
		# Sort:
		order = np.argsort(tick_array[:,0])
		tick_arrays += [tick_array[order,:]]

	return tick_arrays

def _generate_axes_boxes(tick_arrays, xlim, ylim, width, canvas, linewidth):
	"""
	Generates axes boxes frame in a coordinate system
	[0,1]x[0,1].
	
	   canvas : class `Rect`
	"""

	# Sanity checks:
	if not isinstance(canvas,Rect):
		raise RuntimeError("_generate_axes_boxes: 'canvas' has to be Rect!")

	# Preparations:
	ticks_sorted = []
	LW = linewidth
	boxes = []
	colors = []
	for i in range(4):
		# Convert the tick positions to interval [w_rel,1-w_rel]:
		lim = xlim if i < 2 else ylim
		normalized = (tick_arrays[i][:,0] - lim[0]) / (lim[1]-lim[0])
		if i < 2:
			x = np.concatenate([np.zeros((1,)),normalized,np.ones((1,))]) \
			    * (canvas.width()-2*width-LW)
		else:
			x = np.concatenate([np.zeros((1,)),normalized,np.ones((1,))]) \
			    * (canvas.height()-2*width-LW)

		# Create polygons:
		if i == 0:
			# Bottom axis:
			xy = [(canvas.x0+x[j]+width+0.5*LW, canvas.y0+0.5*LW) for j in range(len(x)-1)]
			boxes += [Rectangle(xy[j],
			                    x[j+1]-x[j], width)
			          for j in range(len(x)-1)]
		elif i == 1:
			# Top axis:
			xy = [(canvas.x0+x[j]+width+0.5*LW, canvas.y1-width-0.5*LW) for j in range(len(x)-1)]
			boxes += [Rectangle(xy[j],
			                    x[j+1]-x[j], width,
			                    facecolor = 'k' if j % 2 == 0 else 'w')
			          for j in range(len(x)-1)]
		elif i == 2:
			# Left axis:
			xy = [(canvas.x0+0.5*LW, canvas.y1-width-x[j+1]-0.5*LW) for j in range(len(x)-1)]
			boxes += [Rectangle(xy[j],
			                    width, x[j+1]-x[j],
			                    facecolor = 'k' if j % 2 == 0 else 'w')
			          for j in range(len(x)-1)]
		elif i == 3:
			# Right axis:
			xy = [(canvas.x1-width-0.5*LW,canvas.y1-width-x[j+1]-0.5*LW) for j in range(len(x)-1)]
			boxes += [Rectangle(xy[j],
			                    width, x[j+1]-x[j],
			                    facecolor = 'k' if j % 2 == 0 else 'w')
			          for j in range(len(x)-1)]
		
		print("i =",i,", xy =",xy)
		print("x:",x)
		
		# Save colors:
		colors += [('k' if j % 2 == 0 else 'w') for j in range(len(x)-1)]
	
	# Calculate the remainder of the canvas:
	canvas_remainder = canvas.strip_margin(width)
	
	# Return the boxes:
	return boxes, colors, canvas_remainder
