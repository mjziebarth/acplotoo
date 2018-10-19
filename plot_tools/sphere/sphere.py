# -*- coding: utf-8 -*-
#
# Spherical geometry methods for the plot_tools Python module.
# Copyright (C) 2018 Malte Ziebarth
# 
# This software is distributed under the MIT license.
# See the LICENSE file in this repository.
import numpy as np


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
