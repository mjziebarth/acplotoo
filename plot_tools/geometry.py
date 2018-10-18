# -*- coding: utf-8 -*-
#
# Geometry methods for the sphereplot Python module.
# Copyright (C) 2018 Malte Ziebarth
# 
# This software is distributed under the MIT license.
# See the LICENSE file in this repository.
import numpy as np

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


def great_circle_distance(lon1, lat1, lon2, lat2):
	"""
	Return the pairwise great circle distance between
	two sets of points.
	
	Required arguments:
	   lon1 : Longitude coordinates of first set in
	          degrees.
	   lat1 : Latitude coordinates of first set in
	          degrees.
	   lon2 : Longitude coordinates of second set in
	          degrees.
	   lat2 : Latitude coordinates of second set in
	          degrees.
	
	Returns:
	   Set of great circle distances. Numpy array.
	
	The dimensions of the two sets have to be either
	equal or one of the sets may only be a single
	point.
	
	Equation is combining eqns. (14)-(16) for spherical geometry (f=0) from:
	   T. Vincenty, Direct and Inverse Solutions of Geodesics on the Ellipsoid
	   with Application of Nested Equations, Survey Review 23 (176),
	   (Directorate of Overseas Survey, Kingston Road, Tolworth, Surrey 1975)
	"""
	ctj = np.cos(np.deg2rad(lat2))
	stj = np.sin(np.deg2rad(lat2))
	cti = np.cos(np.deg2rad(lat1))
	sti = np.sin(np.deg2rad(lat1))
	slon = np.sin(np.deg2rad(lon2-lon1))
	clon = np.cos(np.deg2rad(lon2-lon1))
	return np.rad2deg(np.arctan2(np.sqrt((ctj*slon)**2+(cti*stj - sti*ctj*clon)**2),
	                             sti*stj + cti*ctj*clon))



def _gc_rotation_angle(lon1, lat1, lon2, lat2, tolerance=1e-8):
	# Return angle alpha between isolongitude great circle through
	# point 1 and great circle between point 1 and 2.
	# cot(alpha)* tan(delta lon * cos(lat2)) = sin(delta lat)
	if np.abs(lat2-90.0) < tolerance:
		return 0
	elif np.abs(lat2+90.0) < tolerance:
		return 180.0
	D2R = np.pi/180.0
	R2D = 1.0/D2R
	delta_lat = lat2-lat1
	d12 = great_circle_distance(lon1, lat1, lon2, lat2)
	
	# Sides and angles of the spherical triangle in which
	# we want to determine the angle B:
	a = d12
	b = 90.0-lat2
	c = 90.0-lat1
	A = lon2-lon1
	
	arg = ((np.cos(D2R*b)*np.sin(D2R*c) - np.sin(D2R*b)*np.cos(D2R*c)*np.cos(D2R*A)) /
	       np.sin(D2R*a))
	if arg > 1:
		arg = 1
	elif arg < -1:
		arg = -1
	angle = R2D*np.arccos(arg)
	
	# Depending on which longitude comes "before" the other,
	# set the angle to positive or negative:
	if (lon2-lon1) % 360 > 180.0:
		angle = -angle
	
	return -angle



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
	angle = _gc_rotation_angle(lon1, lat1, lon2, lat2)
	x,y,z = rotate_vectors(x,y,z, axis, angle)
	
	return x[::-1],y[::-1],z[::-1]
