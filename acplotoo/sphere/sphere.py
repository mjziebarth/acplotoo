# -*- coding: utf-8 -*-
#
# Spherical geometry methods for the acplotoo Python module.
# Copyright (C) 2018 Malte Ziebarth
# 
# This software is distributed under the MIT license.
# See the LICENSE file in this repository.
import numpy as np

from ..euclidean import rotate_vectors

def to_euclidean_3d(lon, lat, view_center):
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

def to_euclidean_2d(lon, lat, view_center):
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
	x,y,z = to_euclidean_3d(lon, lat, view_center)
	
	# With the choice of coordinates, y>=0 is visible:
	mask = z >= 0.0
	
	return [x[mask],y[mask]]

def from_euclidean(x, y, z, view_center):
	"""
	Convert a set of three-dimensional Euclidean coordinates
	to longitude/latitude coordinates according to the view
	the coordinates were created in.

	Required arguments:
	   x, y, z     : Euclidean coordinates, normed to 1.
	   view_center : View center (lon, lat).
	                 Note: Currently only (0,0) is implemented.
	                 Other view_center will raise error.

	Returns:
	   lon, lat    : Longitude / latitude coordinates.
	"""
	if view_center[0] != 0 or view_center[1] != 0:
		raise RuntimeError("from_euclidean: view_center is not implemented!")

	lon = np.rad2deg(np.arctan2(x,z))
	lat = np.rad2deg(np.arcsin(np.maximum(np.minimum(y,1.0),-1.0)))
	
	return lon, lat
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


def displace(lon, lat, azimuth, distance):
	"""
	Displace one or more points a certain distance along
	a specified azimuth.
	
	Required arguments:
	   lon      : Longitude coordinate(s) of point(s) to
	              displace in degrees.
	   lat      : Latitude coordinate(s) of point(s) to
	              displace in degrees.
	   azimuth  : Azimuth of displacement direction(s) in
	              degrees. It is evaluate at each point's
	              position, where it describes the azimuth
	              of the great circle on which the point
	              is displaced.
	   distance : Distance (in degree) to displace the
	              points by.
	
	Returns:
	   lon, lat : Displaced coordinate(s).
	
	Because of the spherical geometry, displacements of
	multiple points are in general not parallel.
	"""
	
	# Spherical trigonometry.
	# alpha := |delta lon|
	# beta  := azimuth
	# a = distance
	# b = 90 - lat1
	# c = 90 - lat0
	sign_beta = np.sign(180.0 - (azimuth % 360.0))
	cosa = np.cos(np.deg2rad(distance))
	sina = np.sin(np.deg2rad(distance))
	cosc = np.cos(np.deg2rad(90.0-lat))
	sinc = np.sin(np.deg2rad(90.0-lat))
	cosb = cosc*cosa + sinc*sina*np.cos(np.deg2rad(azimuth))
	sinb = sign_beta * np.sqrt(1.0-cosb**2)
	alpha = np.rad2deg(np.arccos((cosa - cosb*cosc) / sinb*sinc))
	
	# Results:
	lat1 = 90.0 - np.rad2deg(np.arccos(cosb))
	dlon = sign_beta * alpha
	lon1 = lon + dlon
	
	return lon1, lat1


def azimuth(lon1, lat1, lon2, lat2, tolerance=1e-8):
	"""
	Return at point 1 the azimuth of the great circle
	between point 1 and 2.
	
	Required arguments:
	   lon1      : Longitude coordinate of point 1 in degrees.
	   lat1      : Latitude coordinate of point 1 in degrees.
	   lon2      : Longitude coordinate of point 2 in degrees.
	   lat2      : Latitude coordinate of point 2 in degrees.
	
	Optional arguments:
	   tolerance : 
	
	Returns:
	   azimuth
	"""
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
	
	return angle


def small_circle_points(lon, lat, r, seg_len, min_points=16):
	"""
	Return a segmentation of a small (or great) circle on
	a sphere, i.e. a set of uniform density on the circle.

	Required arguments:
	   lon, lat : Coordinates of the circle's center in degrees.
	   r        : Radius of the circle in degrees.
	   seg_len  : Spacing of the points in degrees (approximately).

	Optional arguments:
	   min_points : Minimum number of segments.
	                (Default: 16)

	Returns:
	   lon, lat : Coordinates of the circle segmentation.
	"""
	
	# Determine number of points:
	N = max(int(np.ceil(np.sin(np.deg2rad(r)) * 360.0/seg_len)), min_points)
	
	print("N:",N)
	
	# Create the circle segmentation:
	lons, lats = np.linspace(0,360,N), np.ones(N)*(90.0-r)
	x,y,z = to_euclidean_3d(lons, lats, (0.0,0.0))

	# Rotate to target latitude:
	axis = np.array(to_euclidean_3d(lon-90.0, 0.0, (0.0, 0.0)))\
	       .reshape((3,))
	x,y,z = rotate_vectors(x,y,z, axis, -(90.0-lat))

	# Obtain longitude and latitude coordinates:
	return from_euclidean(x,y,z, (0.0, 0.0))
