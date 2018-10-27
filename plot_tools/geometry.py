# -*- coding: utf-8 -*-
#
# Geometry methods for the sphereplot Python module.
# Copyright (C) 2018 Malte Ziebarth
# 
# This software is distributed under the MIT license.
# See the LICENSE file in this repository.
import numpy as np
from .sphere import azimuth as _azimuth, great_circle_distance

def rotation_matrix(axis, theta):
	"""
	Return the rotation matrix associated with counterclockwise
	rotation about the given axis.
	
	Required arguments:
	   axis  : A three-dimensional Euclidean vector, the
	           rotational axis. Need not be normed.
	   theta : The rotational angle in degrees.
	
	Returns:
	   A 3x3 rotation matrix. Can be dot-multiplied with a
	   vector to yield 
	"""
	theta = theta*np.pi/180.0
	axis = np.asarray(axis)
	axis = axis/np.sqrt(np.dot(axis, axis))
	a = np.cos(theta/2.0)
	b, c, d = -axis*np.sin(theta/2.0)
	aa, bb, cc, dd = a*a, b*b, c*c, d*d
	bc, ad, ac, ab, bd, cd = b*c, a*d, a*c, a*b, b*d, c*d
	return np.array([[aa+bb-cc-dd, 2*(bc+ad), 2*(bd-ac)],
	                 [2*(bc-ad), aa+cc-bb-dd, 2*(cd+ab)],
	                 [2*(bd+ac), 2*(cd-ab), aa+dd-bb-cc]])


def rotate_vectors(x,y,z,axis,theta):
	"""
	Rotate vectors around around a given rotational axis
	by the rotational angle theta.
	
	Required arguments:
	   x     : Numpy array of x coordinates.
	   y     : Numpy array of y coordinates.
	   z     : Numpy array of z coordinates.
	   axis  : A three-dimensional Euclidean vector, the
	           rotational axis. Need not be normed.
	   theta : The rotational angle in degrees.
	
	Returns:
	   x,y,z coordinate arrays of same input shape rotated
	   accordingly.
	"""
	shape = x.shape
	M = rotation_matrix(axis,theta)
	V = np.dot(M, np.array([np.array(x.flat),np.array(y.flat),np.array(z.flat)]))
	return V[0].reshape(shape), V[1].reshape(shape), V[2].reshape(shape)


def convert_coordinates_3d(lon, lat, view_center):
	"""
	Convert a set of lon-lat-coordinates to three dimensional
	Euclidean space and rotate that Euclidean space according
	to a view.
	
	Required arguments:
	   lon         : Numpy array of longitude coordinates in
	                 degrees.
	   lat         : Numpy array of latitude coordinates in
	                 degrees.
	   view_center : A pair (lon,lat) giving the position of
	                 the view vector in the same coordinate
	                 system as the coordinates to be converted.
	
	Returns:
	   x,y,z coordinate arrays rotated accordingly.
	"""
	D2R = np.pi/180.0
	# First rotation: Longitude.
	if view_center[0] != 0:
		lon = lon - view_center[0]
	
	# Convert to 3d coordinates:
	x = np.array((np.sin(D2R*lon)*np.cos(D2R*lat)).flat)
	z = np.array((np.cos(D2R*lon)*np.cos(D2R*lat)).flat)
	y = np.array(np.sin(D2R*lat).flat)
	
	# Rotate according to latitude of view_center:
	if view_center[1] != 0:
		axis = (1.0, 0, 0)
		#axis = (np.cos(D2R*view_center[0]), 0, np.sin(D2R*view_center[0]))
		x,y,z = rotate_vectors(x,y,z, axis, view_center[1])
	
	return x,y,z


def convert_coordinates(lon, lat, view_center):
	"""
	Convert a set of lon-lat-coordinates to three dimensional
	Euclidean space and rotate that Euclidean space according
	to a view.
	
	Required arguments:
	   lon         : Numpy array of longitude coordinates in
	                 degrees.
	   lat         : Numpy array of latitude coordinates in
	                 degrees
	   view_center : A pair (lon,lat) giving the position of
	                 the view vector in the same coordinate
	                 system as the coordinates to be converted.
	
	Returns:
	   [x,y] : coordinate arrays rotated accordingly. Only
	           coordinates at points where z>0 are returned.
	"""
	x,y,z = convert_coordinates_3d(lon, lat, view_center)
	
	# With the choice of coordinates, y>=0 is visible:
	mask = z >= 0.0
	
	return [x[mask],y[mask]]


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
	x,y,z = convert_coordinates_3d(lons, lats, view_center)
	axis = np.array(convert_coordinates_3d(lon1+90, 0, view_center)).reshape((3,))
	x,y,z = rotate_vectors(x,y,z, axis,90-d-lat1)
	
	# Rotate the great circle around point one to point 2:
	axis = np.array(convert_coordinates_3d(lon1, lat1, view_center)).reshape((3,))
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
