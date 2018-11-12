# -*- coding: utf-8 -*-
#
# Plot tools geoplot backend file.
# Copyright (C) 2018 Malte Ziebarth
# 
# This software is distributed under the MIT license.
# See the LICENSE file in this repository.

import numpy as np
import struct

from .rect import Rect
from ..cache import has_joblib, plot_tools_cache
from matplotlib.patches import Rectangle
from shapely.geometry.polygon import Polygon as SPoly
from shapely.geometry import LineString, MultiPolygon, box, MultiLineString, MultiPolygon
from shapely.ops import polygonize_full, unary_union
from shapely.prepared import prep
from shapely import speedups

from datetime import datetime

# Enable speedups:
if speedups.available:
	speedups.enable()


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
	LW = linewidth / 72.
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
			xy = [(canvas.x0+x[j]+width+0.5*LW, canvas.y0+0.5*LW)
			      for j in range(len(x)-1)]
			boxes += [Rectangle(xy[j],
			                    x[j+1]-x[j], width)
			          for j in range(len(x)-1)]
		elif i == 1:
			# Top axis:
			xy = [(canvas.x0+x[j]+width+0.5*LW, canvas.y1-width-0.5*LW)
			      for j in range(len(x)-1)]
			boxes += [Rectangle(xy[j],
			                    x[j+1]-x[j], width,
			                    facecolor = 'k' if j % 2 == 0 else 'w')
			          for j in range(len(x)-1)]
		elif i == 2:
			# Left axis:
			xy = [(canvas.x0+0.5*LW, canvas.y0+x[j]+width+0.5*LW)
			      for j in range(len(x)-1)]
			boxes += [Rectangle(xy[j],
			                    width, x[j+1]-x[j],
			                    facecolor = 'k' if j % 2 == 0 else 'w')
			          for j in range(len(x)-1)]
		elif i == 3:
			# Right axis:
			xy = [(canvas.x1-width-0.5*LW,canvas.y0+x[j]+width+0.5*LW)
			      for j in range(len(x)-1)]
			boxes += [Rectangle(xy[j],
			                    width, x[j+1]-x[j],
			                    facecolor = 'k' if j % 2 == 0 else 'w')
			          for j in range(len(x)-1)]
		
		# Save colors:
		if i < 2:
			colors += [('k' if j % 2 == 0 else 'w') for j in range(len(x)-1)]
		else:
			colors += [('k' if j % 2 == 0 else 'w') for j in range(len(x)-1)][::-1]
	
	# Calculate the remainder of the canvas:
	canvas_remainder = canvas.strip_margin(width+0.5*LW)
	
	# Return the boxes:
	return boxes, colors, canvas_remainder


def _generate_axes_ticks(tick_arrays, grid_lons, grid_lats, ticks_between,
                         xlim, ylim, canvas, projection, box_axes_width, linewidth):
	# Generate axes tick lines!
	XY = []
	LW = linewidth
	grid_ticks = [grid_lons, grid_lats]

	for i in range(4):
		# See whether we have ticks on this axis:
		tick_array = tick_arrays[i]
		x = tick_array[:,0]
		if tick_array.shape[0] == 0:
			continue

		# Now see if the ticks are part of the grid ticks:
		if i < 2:
			lons, lats = projection.inverse(x, ylim[i]*np.ones_like(x))
		else:
			lons, lats = projection.inverse(xlim[i-2]*np.ones_like(x),x)

		is_grid_tick = [np.any(np.isclose([lons,lats][int(tick_array[i,1])][i],
		                                  grid_ticks[int(tick_array[i,1])]))
		                for i in range(tick_array.shape[0])]
		tick_array = tick_array[is_grid_tick,:]
		x = tick_array[:,0]

		# Convert the tick positions to interval [w_rel,1-w_rel]:
		# This should be mirrored from the axes boxes code!
		lim = xlim if i < 2 else ylim
		normalized = (x - lim[0]) / (lim[1]-lim[0])

		# Calculate window coordinates w:
		if i < 2:
			w = normalized * (canvas.width()-2*box_axes_width-LW)
		else:
			w = normalized * (canvas.height()-2*box_axes_width-LW)

		# Determine lon/lat coordinates:
		if i < 2:
			# Bottom & top:
			lon, lat = projection.inverse(x, lim[i]*np.ones_like(x))
		else:
			# Left & right:
			lon, lat = projection.inverse(lim[i-2]*np.ones_like(x), x)

		# The tick array (created in _create_tick_array) contains an index at
		# positions [:,1]. This index is 0 for longitude ticks and 1 for latitude
		# ticks. In the former case, we need the eastwards unit vector, in the
		# latter the northward.
		east = (tick_array[:,1] == 1)
		v = np.zeros((w.size, 2))
		v[east, :] = projection.unit_vector_east(lon=lon[east],lat=lat[east])
		v[~east,:] = projection.unit_vector_north(lon=lon[~east],lat=lat[~east])

		sign_v = np.sign(v)

		# Create lines:
		xy = np.zeros((w.size,ticks_between+2,2))
		if i == 0:
			# Bottom axis:
			x0 = canvas.x0 + w + box_axes_width + 0.5*LW
			x1 = canvas.x0 + w + box_axes_width + 0.5*LW \
			     - v[:,0] * box_axes_width * sign_v[:,1]
			y0 = (canvas.y0 + 0.5*LW + box_axes_width) * np.ones_like(w)
			y1 = canvas.y0 + 0.5*LW \
			     + (1.0 - np.abs(v[:,1])) * box_axes_width
		elif i == 1:
			# Top axis:
			x0 = canvas.x0 + w + box_axes_width + 0.5*LW
			x1 = canvas.x0 + w + box_axes_width + 0.5*LW \
			     + v[:,0] * box_axes_width * sign_v[:,1]
			y0 = (canvas.y1 - box_axes_width - 0.5*LW) * np.ones_like(w)
			y1 = canvas.y1 - 0.5*LW \
			     - (1.0 - np.abs(v[:,1])) * box_axes_width
		elif i == 2:
			# Left axis:
			x0 = (canvas.x0 + 0.5*LW + box_axes_width) * np.ones_like(w)
			x1 = canvas.x0 + 0.5*LW \
			     + (1.0 - np.abs(v[:,0])) * box_axes_width
			y0 = canvas.y0 + w + box_axes_width + 0.5*LW
			y1 = canvas.y0 + w + box_axes_width + 0.5*LW \
			                   - v[:,1] * box_axes_width * sign_v[:,0]
		elif i == 3:
			# Right axis:
			x0 = (canvas.x1 - 0.5*LW - box_axes_width) * np.ones_like(w)
			x1 = canvas.x1 - 0.5*LW \
			     - (1.0 - np.abs(v[:,0])) * box_axes_width
			y0 = canvas.y0 + w + box_axes_width + 0.5*LW
			y1 = canvas.y0 + w + box_axes_width + 0.5*LW \
			     + v[:,1] * box_axes_width * sign_v[:,0]

		# Add ticks between:
		for k in range(w.size):
			xy[k,:,0] = np.linspace(x0[k],x1[k],ticks_between+2)
			xy[k,:,1] = np.linspace(y0[k],y1[k],ticks_between+2)

		XY += [xy]

	# Connect all XY:
	if len(XY) > 0:
		xy = np.concatenate(XY,axis=0)
	else:
		xy = None

	# Calculate the remainder of the canvas:
	# TODO if XY[i] is empty, it is possible not to strip a margin off that
	# side!
	canvas_remainder = canvas.strip_margin(box_axes_width+0.5*LW)

	return xy, canvas_remainder


def identify_jumps(c, lim):
	# Identify points in a line where a jump between left
	# and right or top and bottom happens outside the plot
	# canvas:
	c1 = np.roll(c,-1)
	cmin = np.minimum(c,c1)
	cmax = np.maximum(c,c1)
	jump = np.logical_and(cmin < lim[0], cmax > lim[1])

	# Return boolean array that is True at the index right
	# before the jump:
	return jump


# GSHHG:
def gshhg_read_header(bytes,i0):
	# Read a GSHHG header:
	i, n, flag, west, east, south, north, area, \
	area_full,container,ancestor \
	    = struct.unpack('>IIIiiiiIIII', bytes[i0:i0+44])

	# Return:
	return i0+44, (i, n, flag, west, east, south, north, area, area_full,
	               container, ancestor)


@plot_tools_cache.cache(ignore=["projection", "verbose"])
def read_coastlines(gshhg_path, projection, projection_identifier, geographic_extents,
                    xlim, ylim, verbose):
	with open(gshhg_path, 'rb') as f:
		bytes = f.read()

	if verbose > 0:
		print("Reading coastlines...")
		t0 = datetime.now()

	N = len(bytes)
	coasts = []
	raw_data = []
	i=0
	while i < N:
		# Read header:
		i, res = gshhg_read_header(bytes,i)

		# Obtain level and number of points:
		n = res[1]
		level = res[2] % 256

		# Read points:
		xy = np.frombuffer(bytes[i:i+8*n],dtype=np.dtype('>i4'),count=2*n)
		i += 8*n

		# Create lat/lon array:
		latlon = 1e-6*xy.reshape((n,2))

		# Ensure correct coordinate range:
		latlon[:,0] = (latlon[:,0] + 180.0) % 360.0 - 180.0
		latlon[:,1] = (latlon[:,1] + 90.0) % 180.0 - 90.0

		raw_data += [(level, latlon)]

	if verbose > 1:
		t1 = datetime.now()
		print("took:",t1-t0)
		print("Filtering coastline coordinates...")

	filtered_data = []
	for d in raw_data:
		outside = np.logical_or(
		             np.logical_or(d[1][:,0] < geographic_extents[0,0],
		                           d[1][:,0] > geographic_extents[0,1]),
		             np.logical_or(d[1][:,1] < geographic_extents[1,0],
		                           d[1][:,1] > geographic_extents[1,1])
		                        )
		if not np.all(outside):
			filtered_data += [d]

	if verbose > 1:
		t2 = datetime.now()
		print("   took:",t2-t1)
		print("   remaining:",len(filtered_data))
		print("Converting coastline coordinates...")

	for dat in filtered_data:
		# Convert to xy:
		x, y = projection.project(dat[1][:,0], dat[1][:,1])

		# Do a second filtering, now in xy coordinates:
		if np.any(np.logical_and(
		              np.logical_and(x >= xlim[0], x <= xlim[1]),
		              np.logical_and(y >= ylim[0], y <= ylim[1]))
		          ):
			# At least one point is inside canvas:
			coasts += [(dat[0], x, y)]


	if verbose > 1:
		t3 = datetime.now()
		print("Coastlines read! Took:",t3-t0)

	return coasts


@plot_tools_cache.cache
def coast_line_patches_and_path(coords, cnvs_x, cnvs_y, xclip, yclip):
	"""
	coords:          In canvas coordinates.
	cnvs_x, cnvs_y : Visible canvas coordinates in canvas coordinates
	                 (the part behind the axes is invisible)
	xclip, yclip :   In data coordinates.
	"""
	# Optimize the paths:
	patches_xy = []
	coast_path = []
	for i in range(len(coords)):
		c = coords[i]
		outside = np.logical_or(np.logical_or(c[0] < cnvs_x[0], c[0] > cnvs_x[1]),
		                        np.logical_or(c[1] < cnvs_y[0], c[1] > cnvs_y[1]))
		if np.any(outside):
			# Split polygons at points where a coordinate jump over the plot
			# are is performed. We assume this can only happen by wrapping around
			# back of sphere. This assumption should be good if the polygons are
			# 'reasonably' sampled.

			both_outside = np.logical_and(outside, np.roll(outside,1))
			jump_x = identify_jumps(c[0], cnvs_x)
			jump_y = identify_jumps(c[1],cnvs_y)

			split = np.logical_and(both_outside, np.logical_or(jump_x,jump_y))

			# Split polygons at those points:
			XY = np.concatenate([c[0][:,np.newaxis],
			                     c[1][:,np.newaxis]],axis=1)
			XY = np.split(XY, np.argwhere(split).flatten()+1, axis=0)

			# Add sub polygons:
			for xy in XY:
				patches_xy += [xy]

				# Add polygon for sea shore to global polygon:
				if c[2] == 1 and xy.shape[0] > 2:
					# This is a level 1 polygon, AKA sea border!
					# Polygonize, thanks to Mike T on stackoverflow!
					closed = np.concatenate([xy,xy[0:2]],axis=0)
					ls = LineString(closed)
					mls = unary_union(ls)
					coast_path += polygonize_full(mls)[0]
		else:
			# If all inside, all is well!

			# Create patch:
			xy = np.concatenate([c[0][:,np.newaxis],c[1][:,np.newaxis]], axis=1)
			patches_xy += [xy]

			# Add polygon for sea shore to global polygon:
			if c[2] == 1 and xy.shape[0] > 2:
				spoly = SPoly(xy)
				coast_path += [spoly]


	# Finalize coast path:
	clip_box = box(xclip[0],yclip[0],xclip[1],yclip[1])
	coast_path = MultiPolygon([clip_box.intersection(poly)
	                           for poly in coast_path])
	coast_path = unary_union(coast_path)

	return patches_xy, coast_path
