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

def _generate_axes_boxes(tick_dict, xlim, ylim, width, canvas, linewidth):
	"""
	Generates axes boxes frame in a coordinate system
	[0,1]x[0,1].
	
	   canvas : class `Rect`
	"""
	# TODO remove!
	print("tick_dict:",tick_dict)
	print("canvas:   ",canvas)

	# Sanity checks:
	if not isinstance(canvas,Rect):
		raise RuntimeError("_generate_axes_boxes: 'canvas' has to be Rect!")
	
	# First, for each axis generate a sorted list of ticks:
	axes = ["bot","top","left","right"]
	ticks_sorted = []
	LW = linewidth
	for i in range(4):
		ticks_unsorted = tick_dict[axes[i]]
		n0 = ticks_unsorted[0].size
		n1 = ticks_unsorted[1].size
		n_ticks = n0 + n1
		# In tick_array, we save all ticks and mark them by
		# coordinate:
		tick_array = np.zeros((n_ticks,2))
		tick_array[0:n0,0] = ticks_unsorted[0]
		tick_array[n0:,0] = ticks_unsorted[1]
		tick_array[n0:,1] = 1
		# Sort:
		order = np.argsort(tick_array[:,0])
		tick_array = tick_array[order,:]
		# Save:
		ticks_sorted += [tick_array]
	
	# Now generate the boxes:
	boxes = []
	colors = []
	for i in range(4):
		# Convert the tick positions to interval [w_rel,1-w_rel]:
		lim = xlim if i < 2 else ylim
		normalized = (ticks_sorted[i][:,0] - lim[0]) / (lim[1]-lim[0])
		if i < 2:
			x = np.concatenate([np.zeros((1,)),normalized,np.ones((1,))]) * (canvas.width()-2*width-LW)
		else:
			x = np.concatenate([np.zeros((1,)),normalized,np.ones((1,))]) * (canvas.height()-2*width-LW)
		
		
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
