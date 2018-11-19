# -*- coding: utf-8 -*-
#
# Euclidean geometry methods for the acplotoo Python module.
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
