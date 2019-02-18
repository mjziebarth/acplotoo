# -*- coding: utf-8 -*-
#
# Acplotoo geoplot backend file.
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

def _ensure_coordinate_format(x0, x1, data, target_type, target_indexing='ij'):
	"""
	Required Parameters:
		x0,x1 : First and second coordinates.
		data  : Data array.
		target_type : One of:
		              - 'lll' a coordinate pair for each data point)
		              - 'llg' data given in grid, coordinates describe one of the axes.
		                      If data grid axes are of different size, coordinates can
		                      be matched to the same-sized axes. Otherwise, x0 is
		                      assumed to iterate over index i and x1 over index j of
		                      data[i,j], i.e. data[i,j] corresponds to (x0[i],x1[j])
		              - 'ggg' data and coordinates given as equal-sized grids each.
	"""
	# 1) Ensure we can use all the numpy methods:
	x0 = np.atleast_1d(x0)
	x1 = np.atleast_1d(x1)
	data = np.atleast_1d(data)

	# 2) Classify input coordinate type:
	ndim_x0 = x0.squeeze().ndim
	ndim_x1 = x1.squeeze().ndim
	ndim_data = data.squeeze().ndim
	if ndim_data < max(ndim_x0,ndim_x1):
		raise TypeError("Data dimensions must be at least that of the coordinate "
		                "dimensions!")
	if ndim_x0 > 2 or ndim_x1 > 2 or ndim_data > 3:
		raise TypeError("Support only ndim <= 2 for coordinate arrays and ndim <= 3 "
		                "for data arrays!")
	if ndim_data == 3:
		data_point_size = data.shape[2]
	else:
		data_point_size = 1

	if x0.size != x1.size:
		# Both coordinate arrays are of different size.
		# This can only be accomplished by x0 and x1
		# describing the coordinates of a grid.
		if ndim_x0 != 1 or ndim_x1 != 1 or (ndim_data != 2 and ndim_data != 3):
			raise TypeError("If sizes of coordinate arrays differ, they must describe "
			                "one-dimensional axes of data grid!")
		if x0.size == data.shape[0]:
			indexing = 'ij'
		elif x1.size == data.shape[0]:
			indexing = 'xy'
		else:
			raise TypeError("Coordinates/data shape not compliant!")

		src_type = 'llg'
	else:
		# Both coordinate arrays are of same size.
		# This can either be accomplished by x0 and x1 describing
		# the coordinate axes of a grid, or the coordinates of each
		# data points. This is decided by comparing to the data shape.
		if data_point_size * x0.size == data.size:
			# Each data point is associated with one coordinate point.
			if ndim_x0 == 1:
				src_type = 'lll'
			else:
				src_type == 'ggg'
		else:
			# Grid mode again.
			src_type = 'llg'
			if ndim_x0 != 1 or ndim_x1 != 1 or len(x0.shape) != len(x1.shape):
				raise TypeError("Coordinate shape not compliant!")
			if len(x0.shape) == 1:
				indexing = 'ij'
			else:
				if x0.shape[0] == 1:
					indexing = 'xy'
				else:
					indexing = 'ij'

	x0 = x0.squeeze()
	x1 = x1.squeeze()
	data = data.squeeze()
	if src_type == 'llg':
		if ndim_data not in [2,3] or (ndim_data != 3 and data_point_size != 1):
			raise TypeError("Data shape not compliant (dim=" + str(ndim_data)
				            + ", shape=" + str(data.shape)
				            + ", point_size=" + str(data_point_size)
				            + ", x0.size=" + str(x0.size)
				            + ", x1.size=" + str(x1.size) + ")" )
		if indexing == 'ij':
			if x0.size != data.shape[0] or x1.size != data.shape[1]:
				raise TypeError("Data shape not compliant (dim=" + str(ndim_data)
					            + ", shape=" + str(data.shape)
					            + ", point_size=" + str(data_point_size)
					            + ", x0.size=" + str(x0.size)
					            + ", x1.size=" + str(x1.size) + ")" )
		elif indexing == 'xy':
			if x0.size != data.shape[1] or x1.size != data.shape[0]:
				raise TypeError("Data shape not compliant (dim=" + str(ndim_data)
					            + ", shape=" + str(data.shape)
					            + ", point_size=" + str(data_point_size)
					            + ", x0.size=" + str(x0.size)
					            + ", x1.size=" + str(x1.size) + ")" )
	elif src_type in ['lll','ggg']:
		if not np.array_equal(x0.shape, data.shape[:ndim_x0]) \
		or not np.array_equal(x1.shape, x0.shape):
			raise TypeError("When indexing each data point, the "
			                "coordinates have to be of same shape "
			                "as the data!")

	# For grid-grid-grid source type, determine the indexing order
	# and ascertain that the coordinates are according to grid definition:
	if src_type == 'ggg':
		if x0[0,0] != x0[1,0]:
			indexing = 'ij'
			# Test:
			assert np.all(np.diff(x0[:,0] > 0))
			assert np.all(np.diff(x0,axis=1) == 0)
			assert np.all(np.diff(x1[0,:] > 0))
			assert np.all(np.diff(x1,axis=0) == 0)
		else:
			indexing = 'xy'
			# Test:
			assert np.all(np.diff(x0[0,:] > 0))
			assert np.all(np.diff(x0,axis=0) == 0)
			assert np.all(np.diff(x1[:,0] > 0))
			assert np.all(np.diff(x1,axis=1) == 0)

	# 3) Now convert to desired output coordinate:
	if src_type == 'lll':
		if target_type != 'lll':
			raise NotImplementedError("Linear coordinates to anything else not yet "
			                          "implemented.")
		else:
			# Nothing to do!
			pass
	elif src_type == 'llg':
		if target_type in ['lll','ggg']:
			x0,x1 = np.meshgrid(x0,x1,indexing=indexing)

			if target_type == 'lll':
				x0 = x0.flatten()
				x1 = x1.flatten()
				data = data.flatten()
		else:
			# Nothing to do!
			pass
	elif src_type == 'ggg':
		if target_type == 'lll':
			x0 = x0.flatten()
			x1 = x1.flatten()
			data = data.flatten()
		elif target_type == 'llg':
			if indexing == 'ij':
				x0 = x0[:,0].flatten()
				x1 = x1[0,:].flatten()
			elif indexing == 'xy':
				x0 = x0[0,:].flatten()
				x1 = x1[:,0].flatten()
		else:
			# Nothing to do:
			pass

	# See if we have to transpose:
	if indexing != target_indexing:
		data = data.T
		if target_type == 'ggg':
			x0 = x0.T
			x1 = x1.T

	# Return the data!
	return x0, x1, data



def _create_tick_arrays(tick_dict, which_ticks):
	"""
	Generate arrays of ticks to be used with _generate_axes stuff
	"""
	axes = ["bot","top","left","right"]
	ticks_unsorted = [tick_dict[a] for a in axes]
	tick_arrays = []
	tick_arrays = [[],[],[],[]]
	# Determine number of ticks by type and axis:
	Nlon = np.zeros(4,dtype=int)
	Nlat = np.zeros(4,dtype=int)
	for i in range(4):
		Nlon[i] = tick_dict[axes[i]][0].size
		Nlat[i] = tick_dict[axes[i]][1].size

	# Here determine which ticks j (0=lon, 1=lat) we want to display
	# at axis i.
	if which_ticks == 'both':
		J = [(0,1) for i in range(4)]
	elif which_ticks == 'lonlat':
		J = [(i // 2,) for i in range(4)]
	elif which_ticks == 'latlon':
		J = [((i // 2 + 1) % 2,) for i in range(4)]
	elif which_ticks == 'significant':
		# Optimize the number of ticks under the constraint of having
		# two axes for each tick.
		I = set(range(4))
		IJ = [(i,j) for i in I for j in [l for l in I if l > i]]
		KL = [tuple(I-set(ij)) for ij in IJ]

		cost = [(Nlon[IJ[m][0]] + Nlon[IJ[m][1]]) * (Nlat[KL[m][0]] + Nlat[KL[m][1]])
		        for m in range(len(IJ))]
		m_max = cost.index(max(cost))
		J = list(range(4))
		J[IJ[m_max][0]] = (0,)
		J[IJ[m_max][1]] = (0,)
		J[KL[m_max][0]] = (1,)
		J[KL[m_max][1]] = (1,)


	# Now use that info to save the corresponding ticks in the tick arrays
	# in a linearly ordered fashion:
	for i in range(4):
		# Create an array holding the linear coordinates tick_array[:,0]
		# and the tick type tick_array[:,1] (0=lon, 1=lat).
		n_ticks = 0
		array = []
		use_lon = 0 in J[i]
		use_lat = 1 in J[i]
		if use_lon:
			n_ticks += Nlon[i]
			nlon = Nlon[i]
		else:
			nlon = 0
		if use_lat:
			n_ticks += Nlat[i]

		tick_array = np.zeros((n_ticks,2))

		if use_lon:
			tick_array[:nlon,0] = ticks_unsorted[i][0]
		if use_lat:
			tick_array[nlon:,0] = ticks_unsorted[i][1]
			tick_array[nlon:,1] = 1

		# Sort the ticks by linear coordinate at axis:
		order = np.argsort(tick_array[:,0])
		tick_arrays[i] = tick_array[order,:]

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


def _generate_axes_ticks(tick_arrays, grid_lons, grid_lats,
                         xlim, ylim, canvas, projection, box_axes_width, linewidth,
                         axes_margin, gplt):
	# Generate axes tick lines!
	XY = []
	LW = linewidth / 72.
	AM = axes_margin / 72.
	grid_ticks = [grid_lons, grid_lats]
	tick_masks = []
	ticks = []
	has_ticks = []
	label_texts = []
	labels = []
	required_label_space = []

	# Step 1: Pre-render the axes labels to determine the margin we need.
	for i in range(4):
		tick_array = tick_arrays[i]
		has_ticks += [tick_array.shape[0] != 0]
		if not has_ticks[i]:
			continue
		
		x = tick_array[:,0]

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

		# Save the grid tick mask:
		tick_masks += [is_grid_tick]

		# Determine lon/lat coordinates:
		lim = xlim if i < 2 else ylim
		if i < 2:
			# Bottom & top:
			lon, lat = projection.inverse(x, lim[i]*np.ones_like(x))
		else:
			# Left & right:
			lon, lat = projection.inverse(lim[i-2]*np.ones_like(x), x)

		# Save ticks:
		ticks += [(lon, lat)]

		# The tick array (created in _create_tick_array) contains an index at
		# positions [:,1]. This index is 0 for longitude ticks and 1 for latitude
		# ticks. In the former case, we need the northwards unit vector, in the
		# latter the eastward.
		is_lon = tick_array[:,1] == 0

		# Create and prerender the tick labels:
		lonlat_tick_array = [(int(lon[k]) if tick_array[k,1] == 0 else int(lat[k]),
		                          tick_array[k,1])
		                     for k in range(tick_array.shape[0])]
		label_texts += [gplt._determine_tick_labels(lonlat_tick_array)]
		labels += [[gplt._ax.text(0, 0, txt) for txt in label_texts[i]]]

		
		# Determine the width / height we have to reserve for each axis:
		renderer = gplt._ax.figure.canvas.get_renderer()
		if i < 2:
			required_label_space += [max(txt.get_window_extent(renderer).height
			                             for txt in labels[i]) / 72.0]
		else:
			required_label_space += [max(txt.get_window_extent(renderer).width
			                             for txt in labels[i]) / 72.0]


	# Step 2: Place the axes label and the ticks.
	for i in range(4):
		# See whether we have ticks on this axis:
		tick_array = tick_arrays[i]
		x = tick_array[:,0]
		if not has_ticks[i]:
			XY += [np.zeros((0,2,2))]
			tick_masks += [np.zeros(0)]
			continue

		# Convert the tick positions to interval [w_rel,1-w_rel]:
		# This should be mirrored from the axes boxes code!
		lim = xlim if i < 2 else ylim
		normalized = (x - lim[0]) / (lim[1]-lim[0])

		# Calculate window coordinates w:
		if i < 2:
			w = normalized * (canvas.width()-2*(box_axes_width+AM)-LW
			                  - required_label_space[2] - required_label_space[3])
		else:
			w = normalized * (canvas.height()-2*(box_axes_width+AM)-LW
			                  -required_label_space[0] - required_label_space[1])

		# Determine lon/lat coordinates:
		lon, lat = ticks[i]

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
		xy = np.zeros((w.size,2,2))
		if i == 0:
			# Bottom axis:
			xy[:,0,0] = canvas.x0 + w + box_axes_width + AM +0.5*LW \
			            + required_label_space[2]
			xy[:,1,0] = canvas.x0 + w + box_axes_width + AM + 0.5*LW \
			            - v[:,0] * box_axes_width * sign_v[:,1] + required_label_space[2]
			xy[:,0,1] = canvas.y0 + 0.5*LW + box_axes_width + AM \
			            + required_label_space[0]
			xy[:,1,1] = canvas.y0 + 0.5*LW + required_label_space[0] + AM \
			            + (1.0 - np.abs(v[:,1])) * box_axes_width
		elif i == 1:
			# Top axis:
			xy[:,0,0] = canvas.x0 + w + box_axes_width + 0.5*LW + AM \
			            + required_label_space[2]
			xy[:,1,0] = canvas.x0 + w + box_axes_width + 0.5*LW + AM \
			            + v[:,0] * box_axes_width * sign_v[:,1] + required_label_space[2]
			xy[:,0,1] = canvas.y1 - box_axes_width - 0.5*LW - AM - required_label_space[1]
			xy[:,1,1] = canvas.y1 - 0.5*LW - AM - required_label_space[i]\
			     - (1.0 - np.abs(v[:,1])) * box_axes_width
		elif i == 2:
			# Left axis:
			xy[:,0,0] = canvas.x0 + 0.5*LW + box_axes_width + AM \
			            + required_label_space[2]
			xy[:,1,0] = canvas.x0 + 0.5*LW + required_label_space[2] + AM \
			            + (1.0 - np.abs(v[:,0])) * box_axes_width
			xy[:,0,1] = canvas.y0 + w + box_axes_width + 0.5*LW + AM \
			            + required_label_space[0]
			xy[:,1,1] = canvas.y0 + w + box_axes_width + 0.5*LW + AM \
			            + required_label_space[0] - v[:,1] * box_axes_width * sign_v[:,0]
		elif i == 3:
			# Right axis:
			xy[:,0,0] = canvas.x1 - 0.5*LW - box_axes_width - AM \
			             - required_label_space[3]
			xy[:,1,0] = canvas.x1 - 0.5*LW - required_label_space[3] - AM \
			            - (1.0 - np.abs(v[:,0])) * box_axes_width
			xy[:,0,1] = canvas.y0 + w + box_axes_width + 0.5*LW + AM \
			            + required_label_space[0]
			xy[:,1,1] = canvas.y0 + w + box_axes_width + 0.5*LW + AM \
			            + v[:,1] * box_axes_width * sign_v[:,0] \
			            + required_label_space[0]

		# Move the labels:
		if i == 0:
			for k,lbl in enumerate(labels[i]):
				lbl.set_position((xy[k,1,0]-0.5*lbl.get_window_extent(renderer).width / 72.0,
				                 canvas.y0))
		elif i == 1:
			for k,lbl in enumerate(labels[i]):
				lbl.set_position((xy[k,1,0]-0.5*lbl.get_window_extent(renderer).width / 72.0,
				                 canvas.y1 - lbl.get_window_extent(renderer).height / 72.0))
		elif i == 2:
			for k,lbl in enumerate(labels[i]):
				lbl.set_position((canvas.x0,
				                  xy[k,1,1] - 0.5*lbl.get_window_extent(renderer).height / 72.0))
		elif i == 3:
			for k,lbl in enumerate(labels[i]):
				lbl.set_position((canvas.x1 - lbl.get_window_extent(renderer).width / 72.0,
				                  xy[k,1,1] - 0.5*lbl.get_window_extent(renderer).height / 72.0))

		XY += [xy]

	# Calculate the remainder of the canvas:
	# TODO if XY[i] is empty, it is possible not to strip a margin off that
	# side!
	if all([len(xy)==0 for xy in XY]):
		XY = None

	margin = np.array(required_label_space) + box_axes_width + 0.5*LW + AM
	canvas_remainder = canvas.strip_margin((margin[2], margin[0], margin[3], margin[1]))


	return XY, canvas_remainder, tick_masks


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
					polys, dangles, cutedges, invalidrings = polygonize_full(mls)
					coast_path += polys
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
	coast_path = [clip_box.intersection(poly) for poly in coast_path]
	coast_path = [obj.geoms if isinstance(obj,MultiPolygon) else [obj]
	              for obj in coast_path]
	coast_path = MultiPolygon([c for obj in coast_path for c in obj])
	coast_path = unary_union(coast_path)

	return patches_xy, coast_path
