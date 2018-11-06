# -*- coding: utf-8 -*-
#
# Plot tools Hotine oblique Mercator projection submodule
# init file.
# Copyright (C) 2018 Malte Ziebarth
# 
# This software is distributed under the MIT license.
# See the LICENSE file in this repository.

import numpy as np

from .projection import Projection


# Some helper functions:
def _tfun(phi0,e):
	# Method to calculate t or t0 from Snyder (1987):
	return np.tan(0.25*np.pi - 0.5*phi0) / \
		   np.power((1.0 - e*np.sin(phi0)) / (1.0 + e*np.sin(phi0)),0.5*e)

def _arcsin(x):
	# An arcsin fixing numerical instability under the assumption
	# that |x|>1 can only happen by numerical instability, not
	# due to wrong math.
	if abs(x) <= 1.0:
		return np.arcsin(x)
	elif x > 1.0:
		return 0.5*np.pi
	else:
		return -0.5*np.pi


class HotineObliqueMercator(Projection):
	"""
	Hotine Oblique Mercator projection.
	"""
	
	def __init__(self, lon0, lat0, azimuth, k0=1.0, ellipsoid='WGS84',
	             tolerance=1e-7):
		"""
		Optional parameters:
		   k0        : Scale factor along latitude-parallel.
		               Default: 1.0
		   ellipsoid : Determines the reference ellipsoid to use for
		               projection. Determines a (great half axis in m)
		               and f (flattening).
		               Possible values:
		               'WGS84' :
		                    a = 6378.137e3
		                    f = 1/298.257223563
		               (Default: 'WGS84')
		   tolerance : Numerical tolerance (in Degrees). Is used to
		               determine whether a value is zero.
		               (Default: 1e-7)
		"""
		# Set parameters:
		self._lon0 = float(lon0)
		self._lat0 = float(lat0)
		self._azimuth = float(azimuth)
		self._k0 = float(k0)
		self._tolerance = float(tolerance)
		
		if ellipsoid == 'WGS84':
			self._a = 6378.137e3
			self._f = 1/298.257223563
		else:
			raise RuntimeError("Hotine Oblique Mercator: Unknown ellipsoid!")
		
		# Calculate additional parameters we need for calculations:
		self._check_and_set_constants()


	def project_to_uvk(self, lon, lat, return_k=False):
		"""
		Project lon/lat to x/y using a hotine oblique mercator
		projection. Follows alternte B from Snyder (1987), which is
		also used in proj4.
		Projection center lies at (0,0) in projected coordinates.
	
		Required parameters:
		   lon          : Array of longitude coordinates to convert
		   lat          : Array of latitude coordinates to convert
	
		Optional parameters:
		   return_k : Whether to calculate and return the scale factor
			          k. Default: False
	
		Returns:
		   If return_k == False:
			  x,y 
		   otherwise:
			  x,y,k
		"""
		# Unravel constants:
		e,phi0,lambda_c,alpha_c,A,B,t0,D,E,F,G,gamma0,lambda0 = self._constants
		tol = self._tolerance

		# Degree to radians:
		phi = np.deg2rad(lat)
		lbda = np.deg2rad(lon)

		# (9-25) to (9-30), respecting note:
		mask = np.logical_and(lat < 90.0-tol,lat > -90.0+tol)
		t = _tfun(phi[mask],e)
		Q = E / np.power(t,B)
		S = 0.5*(Q - 1.0/Q)
		T = 0.5*(Q + 1.0/Q)

		# Note about handling longitude wrapping:
		dlambda = lbda-lambda0
		dlambda[dlambda < -np.pi] += 2*np.pi
		dlambda[dlambda > np.pi] -= 2*np.pi
		V = np.sin(B*dlambda)
		U = (-V * np.cos(gamma0) + S*np.sin(gamma0)) / T

		# Calculate v
		v = np.zeros_like(phi)
		v[mask] = A * np.log((1.0-U) / (1.0+U)) / (2.0*B)
		v[~mask] = (A/B) * np.log(np.tan(0.25*np.pi + 0.5 * np.sign(phi[~mask])
			                                          * gamma0))

		# Calculate u:
		u = np.zeros_like(phi)
		W = np.cos(B*dlambda)
		mask2 = np.abs(W) > tol
		mask3 = mask2[mask]
		u[np.logical_and(mask,mask2)] = A/B * np.arctan2(S[mask3]*np.cos(gamma0)
			                                             + V[mask3]*np.sin(gamma0), W[mask3])
		u[~mask2] = A*B*dlambda[~mask2]
		u[~mask] = A*phi[~mask]/B

		# Correct for offset to have map centered on u=0:
		uc = np.sign(phi0) * A/B * np.arctan2(np.sqrt(D**2-1), np.cos(alpha_c))

		# Return coordinates:
		if not return_k:
			u -= uc
			return u,v

		# Calculate k, the scale factor:
		k = A * np.cos(B*u/A) * np.sqrt(1-e**2 * np.sin(phi)**2) / \
			(self._a * np.cos(phi) * np.cos(B*dlambda))
		u -= uc

		# Circumvent numerical problems:
		k = np.maximum(k,k0)

		# Now return coordinates:
		return u, v, k


	def inverse_from_uv(self, u, v):
		# Unravel the constants:
		e,phi0,lambda_c,alpha_c,A,B,t0,D,E,F,G,gamma0,lambda0 = self._constants
		tol = self._tolerance

		# Correct u:
		uc = np.sign(phi0) * A/B * np.arctan2(np.sqrt(D**2-1), np.cos(alpha_c))
		u += uc

		# Now calculate (9-42) to (9-47):
		Q = np.exp(-B*v/A)
		S = 0.5*(Q - 1./Q)
		T = 0.5*(Q + 1./Q)
		V = np.sin(B*u/A)
		U = (V*np.cos(gamma0) + S*np.sin(gamma0)) / T
		t = np.power(E / np.sqrt((1. + U)/(1. - U)), 1./B)

		# Consider U= +/-1 case:
		mask = np.logical_or(np.abs(U-1.0) < tol, np.abs(U+1.0) < tol)
		phi = np.zeros_like(v)
		lbda = np.zeros_like(u)

		# The shortcut:
		phi[mask] = np.sign(u[mask]) * 0.25 * np.pi
		lbda[mask] = lambda0

		# The rest iteratively:
		phi_ = np.atleast_1d(0.5*np.pi - 2*np.arctan(t))
		phi_fun = lambda x : 0.5*np.pi - 2*np.arctan(t * np.power((1.-e*np.sin(x))/
			                                        (1.+e*np.sin(x)), 0.5*e))
		phi_new = phi_fun(phi_)
		itermask = np.abs(phi_new - phi_).max() > tol
		while np.any(itermask):
			phi_[itermask] = phi_new[itermask]
			phi_new = phi_fun(phi_[itermask])
			itermask = np.abs(phi_new - phi_[itermask]).max() > tol

		phi[~mask] = phi_

		# Calculate lambda usin (9-48):
		lbda = np.atleast_1d(np.arctan2(S*np.cos(gamma0) - V*np.sin(gamma0),
		                                np.cos(B*u/A)))
		lbda[lbda < -np.pi] += 2*np.pi
		lbda[lbda > np.pi] -= 2*np.pi
		lbda = lambda0 - (lbda / B)

		# Return coordinates:
		lat = np.rad2deg(phi)
		lon = np.rad2deg(lbda)
		lon[lon > 180.0] -= 360.0
		lon[lon < -180.0] += 360.0

		return lon, lat

	### The interface implementation: ###

	def _project(self, lon, lat):
		"""
		Class implementation of _project method.
		"""
		return self.project_to_uvk(lon, lat)


	def _inverse(self, x, y):
		"""
		Class implementation of _inverse method.
		"""
		return self.inverse_from_uv(x, y)
	
	
	### Helper methods: ###

	def _check_and_set_constants(self):
		"""
		Perform conversions from degree to radians and do some
		consistency checks (from Snyder (1987))!
		"""
		# For readability: Attributes to local variables.
		f = self._f
		lat0 = self._lat0
		lon0 = self._lon0
		azimuth = self._azimuth
		tol = self._tolerance

		# Convert flattening to eccentricity:
		e = np.sqrt(2*f - f**2)

		# Convert to radians.
		# Follow notation from Snyder (1987)
		phi0 = np.deg2rad(lat0)
		lambda_c = np.deg2rad(lon0)
		alpha_c  = np.deg2rad(azimuth)

		# Works for -90 < azimuth_c < 90, so transform to that range:
		if azimuth > 90.0:
			alpha_c -= np.pi
		if azimuth < -90.0:
			alpha_c += np.pi

		# Check some limitations from Snyder (1987):
		if np.abs(phi0) < tol:
			raise ValueError("hotine_oblique_mercator: lat_c cannot be zero!")
		if np.abs(phi0-0.5*np.pi) < tol or np.abs(phi0+0.5*np.pi) < tol:
			raise ValueError("hotine_oblique_mercator: lat_c cannot be +/-90°!")

		# Calculate the constants A,B,D,E,F,G,gamma0, and lambda0.

		# Calculate intermediate values from Snyder (1987):
		B = np.sqrt(1.0 + e**2 * np.cos(phi0)**4 / (1.0 - e**2))
		A = self._a * B * self._k0 * np.sqrt(1.0 - e**2) / (1.0 - e**2 * np.sin(phi0)**2)

		# Lambda to calculate t and t0:
		t0 = _tfun(phi0,e)
		D = B * np.sqrt(1.0 - e**2) / (np.cos(phi0) * np.sqrt(1.0 - e**2 * np.sin(phi0)**2))

		# Ensure D >= 1:
		if abs(D) < 1.0:
			D = np.sign(D)

		F = D + np.sign(phi0) * np.sqrt(D**2 - 1.0)
		E = F * np.power(t0,B)
		G = 0.5 * (F - 1.0/F)
		gamma0 = _arcsin(np.sin(alpha_c) / D)
		lambda0 = lambda_c - _arcsin(G * np.tan(gamma0)) / B
		
		# Save constants:
		self._constants = (e, phi0, lambda_c, alpha_c, A, B, t0, D, E, F, G, gamma0, lambda0)
