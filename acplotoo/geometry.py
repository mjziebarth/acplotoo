# -*- coding: utf-8 -*-
#
# Geometry methods for the acplotoo Python module.
# Copyright (C) 2018 Malte Ziebarth
# 
# This software is distributed under the MIT license.
# See the LICENSE file in this repository.
import numpy as np
from .sphere import azimuth as _azimuth, great_circle_distance,\
                    to_euclidean_3d
from .euclidean import rotation_matrix, rotate_vectors

# TODO These methods are probably better to include in sphere.py and
# euclidean.py

def _line_coords(lon1, lat1, lon2, lat2, seg_len, view_center):
	"""
	Calculate coordinates of points on the great circle segment
	between two points.
	"""
	# Start by creating a great circle through point 1
	# and poles:
	d = great_circle_distance(lon1, lat1, lon2, lat2)
	n = int(np.ceil(d / seg_len))
	lons = np.ones(n)*lon1
	lats = np.linspace(90, 90-d, n)
	x,y,z = to_euclidean_3d(lons, lats, view_center)
	axis = np.array(to_euclidean_3d(lon1+90, 0, view_center)).reshape((3,))
	x,y,z = rotate_vectors(x,y,z, axis,90-d-lat1)
	
	# Rotate the great circle around point one to point 2:
	axis = np.array(to_euclidean_3d(lon1, lat1, view_center)).reshape((3,))
	angle = _azimuth(lon1, lat1, lon2, lat2)
	x,y,z = rotate_vectors(x,y,z, axis, -angle)
	
	return x[::-1],y[::-1],z[::-1]


def _z_zero_crossing_xy(x0, y0, z0, x1, y1, z1):
	"""
	Calculate x and y coordinate of a z zero crossing
	using linear interpolation and norming xy vector
	to 1.
	"""
	# We expect z1 > 0, z0 < 0:
	if z0 > 0 or z1 < 0:
		raise RuntimeError("Zero crossing wrong order!")
	
	# Linear interpolation between the two xy coordinate
	# pairs with z as decision parameter:
	dz = z1 - z0
	x = (z1 * x1 - z0 * x0) / dz
	y = (z1 * y1 - z0 * y0) / dz
	
	# Make sure that the coordinates are scaled correctly
	# to one:
	scale = 1.0/np.sqrt(x*x + y*y)
	x *= scale
	y *= scale
	
	# Return results:
	return x,y
	
	

def _spherically_clip_polygon(x, y, z, seg_len):
	"""
	Calculate the path of an object clipped at the z=0 great
	circle.
	"""
	# Identify visible nodes:
	visible = z >= 0
	
	# If all visible, return arrays trivially:
	if np.all(visible):
		return x,y
	
	# Convenience: Invisible nodes:
	invisible = ~visible
	
	# If no visible nodes, return empty:
	if np.all(invisible):
		return None,None
	
	# Identify crossings:
	crs_pos = np.logical_and(invisible[:-1], visible[1:])
	crs_neg = np.logical_and(invisible[1:],visible[:-1])
	
	# Also for the last element the periodic wrapping:
	if invisible[-1] and visible[0]:
		crs_pos[-1] = True
	elif visible[-1] and invisible[0]:
		crs_neg[-1] = True
	
	# Determine the indices:
	crs_pos_id = np.argwhere(crs_pos)
	crs_neg_id = np.argwhere(crs_neg)
	
	# For all crossings, identify the z zero crossing points:
	N = len(z)
	crs_pos_xy = [_z_zero_crossing_xy(x[i], y[i], z[i], x[(i+1) % N],
	                                  y[(i+1) % N], z[(i+1) % N])
	              for i in crs_pos_id]
	crs_neg_xy = [_z_zero_crossing_xy(x[(i+1) % N], y[(i+1) % N], z[(i+1) % N],
	                                  x[i], y[i], z[i]) for i in crs_neg_id]
	
	# Now we can connect between different crossings starting with the first
	# negative crossing:
	output = [] # Ordered list of output coordinates.
	M = len(crs_pos_id) # Number of zero crossings (div 2)
	j_pos = (0 if crs_pos_id[0] > crs_neg_id[0] else 1) % M
	j_neg = 0
	for i in range(M):
		# The two crossing points:
		pt0 = crs_neg_xy[j_neg]
		pt1 = crs_pos_xy[j_pos]
		
		# Determine xy angles:
		theta0 = np.arctan2(pt0[1],pt0[0])
		theta1 = np.arctan2(pt1[1],pt1[0])
		
		# Determine arc len and start/end:
		if theta0 > theta1:
			theta1 += 2*np.pi
		arc_len = np.rad2deg(theta1-theta0)
		n_arc = int(np.ceil(arc_len / seg_len))
		
		# Create arc:
		theta = np.linspace(theta0,theta1,n_arc) % (2*np.pi)
		x_arc = np.cos(theta)
		y_arc = np.sin(theta)
		
		# Append first the zero crossing to output
		# coordinates:
		output += [np.concatenate([x_arc[:,np.newaxis], y_arc[:,np.newaxis]],axis=1)]
		
		# Append the following coordinates from the original
		# shape:
		id0 = crs_pos_id[j_pos]+1
		j_neg = (j_neg+1) % M
		id1 = crs_neg_id[j_neg]
		if id1 > id0:
			# Include both id0 and id1 as well as all between:
			ids = np.arange(id0,id1+1)
		else:
			# If id0 > id1, we have reached the periodic wrapping
			# point of the polygon. We thus need to use the ranges
			# {id0, ..., M-1} and {0, ..., id1}.
			ids = np.concatenate([np.arange(id0,M),np.arange(0,id1+1)], axis=0)
		output += [np.concatenate([x[ids][:,np.newaxis], y[ids][:,np.newaxis]], axis=1)]
		
		# Next crossing:
		j_pos = (j_pos+1) % M
	
	# Since we always include both the hidden and the visible part of each
	# wrap and take care of the periodic wrapping, we are now finished!
	# Combine the output arrays and split into x and y:
	output = np.concatenate(output, axis=0)
	x_out = output[:,0]
	y_out = output[:,1]
	return x_out, y_out
