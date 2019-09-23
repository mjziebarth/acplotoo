# -*- coding: utf-8 -*-
#
# Acplotoo Hotine oblique Mercator projection submodule
# init file.
# Copyright (C) 2018 Malte Ziebarth
# 
# This software is distributed under the MIT license.
# See the LICENSE file in this repository.

import numpy as np

from .projection import Projection

# Try to use proj4:
try:
	from pyproj import Proj
	_has_proj = True
except:
	_has_proj = False

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

def _preprocess_coordinates(x0,x1, cname0, cname1):
	# Make sure we can apply all the numpy routines and ndarray methods:
	x0 = np.atleast_1d(x0)
	x1 = np.atleast_1d(x1)

	# Handle incompatible shapes:
	if not np.array_equal(x0.shape,x1.shape):
		if x0.size != 1 and x1.size != 1:
			raise RuntimeError("Incompatible shapes of %s (%s) and %s (%s)!"
			                   % (cname0, str(x0.shape), cname1, str(x1.shape)))
		elif x0.size == 1:
			x0 = x0 * np.ones_like(x1)
		else:
			x1 = x1 * np.ones_like(x0)

	# Work on one-dimensional arrays:
	shape = list(x0.shape)
	if x0.squeeze().ndim != 1:
		x0 = x0.flatten()
		x1 = x1.flatten()

	return x0, x1, shape


class HotineObliqueMercator(Projection):
	"""
	Hotine Oblique Mercator projection.
	"""

	def __init__(self, lon0, lat0, azimuth, k0=1.0, ellipsoid='WGS84',
	             tolerance=1e-7, invert_v=True, invert_u=False, no_rot=True,
	             use_proj=True):
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
		   invert_v  : Whether to invert the v coordinate. If so, the
		               uv-axes for azimuth=90° equal a typical 
		               xy-coordinate system with x in longitude and
		               y in latitude direction. Otherwise, the v axis
		               will be oriented in negative y direction.
		               (Default: True)
		"""
		# Set parameters:
		self._lon0 = (float(lon0) + 180.0) % 360.0 - 180.0
		self._lat0 = float(lat0)
		self._azimuth = float(azimuth)
		self._k0 = float(k0)
		self._tolerance = float(tolerance)
		self._invert_v = bool(invert_v)
		self._invert_u = bool(invert_u)
		self._ellipsoid = ellipsoid
		self._no_rot = bool(no_rot)
		self._use_proj = bool(use_proj)
		
		if ellipsoid == 'WGS84':
			self._a = 6378.137e3
			self._f = 1/298.257223563
		else:
			raise RuntimeError("Hotine Oblique Mercator: Unknown ellipsoid!")
		
		# Calculate additional parameters we need for calculations:
		self._check_and_set_constants()


	def __hash__(self):
		return hash(self._identifier()[1:])



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
		lon, lat, shape = _preprocess_coordinates(lon, lat, "lon", "lat")

		# Unravel constants:
		e,phi0,lambda_c,alpha_c,A,B,t0,D,E,F,G,gamma0,lambda0,uc = self._constants
		tol = self._tolerance

		# Shortcut if Proj4 exists:
		if self._proj is not None:
			u_,v_ = self._proj(lon, lat)
			u_ -= uc
			if self._invert_v:
				v_ = -v_
			if self._invert_u:
				u_ = -u_
			if not return_k:
				return u_.reshape(shape), v_.reshape(shape)
		else:
			u_,v_ = None, None

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
		V = np.sin(B*dlambda[mask])
		U = (-V * np.cos(gamma0) + S*np.sin(gamma0)) / T

		# Calculate v
		v = np.zeros_like(phi)
		v[mask] = A * np.log((1.0-U) / (1.0+U)) / (2.0*B)
		v[~mask] = (A/B) * np.log(np.tan(0.25*np.pi + 0.5 * np.sign(phi[~mask])
			                                          * gamma0))

		# Invert:
		if self._invert_v:
			v = -v

		# Calculate u:
		u = np.zeros_like(phi)
		W = np.cos(B*dlambda)
		mask2 = np.abs(W) > tol
		mask3 = mask2[mask]
		mask4 = np.logical_and(mask,mask2)
		u[mask4] = A/B * np.arctan2(S[mask3]*np.cos(gamma0) + V[mask3]*np.sin(gamma0), 
			                        W[mask4])
		u[~mask2] = A*B*dlambda[~mask2]
		u[~mask] = A*phi[~mask]/B

		# Return coordinates:
		if not return_k:
			u -= uc
			if self._invert_u:
				return -u.reshape(shape), v.reshape(shape)
			return u.reshape(shape), v.reshape(shape)

		# Calculate k, the scale factor:
		k = A * np.cos(B*u/A) * np.sqrt(1-e**2 * np.sin(phi)**2) / \
			(self._a * np.cos(phi) * np.cos(B*dlambda))
		u -= uc

		# Circumvent numerical problems:
		k = np.maximum(k,self._k0)

		# Take proj4 if found:
		if u_ is not None:
			return u_.reshape(shape), v_.reshape(shape), k.reshape(shape)

		# Now return coordinates:
		if self._invert_u:
			return -u.reshape(shape), v.reshape(shape), k.reshape(shape)
		return u.reshape(shape), v.reshape(shape), k.reshape(shape)


	def inverse_from_uv(self, u, v):
		"""
		Converts projected coordinates (u,v) back to (lon,lat).
		"""
		# Preprocess coordinates:
		u, v, shape = _preprocess_coordinates(u, v, "u", "v")

		# Unravel the constants:
		e,phi0,lambda_c,alpha_c,A,B,t0,D,E,F,G,gamma0,lambda0,uc = self._constants
		tol = self._tolerance

		# Shortcut if Proj4 exists:
		if self._proj is not None:
			if self._invert_u:
				u = -u
			if self._invert_v:
				v = -v
			u = u + uc
			lon, lat = self._proj(u,v,inverse=True)
			return lon.reshape(shape), lat.reshape(shape)


		# Correct u:
		if self._invert_u:
			u = -u
		u = u + uc

		# Correct v:
		if self._invert_v:
			v = -v

		# Now calculate (9-42) to (9-47):
		Q = np.exp(-B*v/A)
		S = 0.5*(Q - 1./Q)
		T = 0.5*(Q + 1./Q)
		V = np.sin(B*u/A)
		U = (V*np.cos(gamma0) + S*np.sin(gamma0)) / T
		t = np.power(E / np.sqrt((1. + U)/(1. - U)), 1./B)

		# Consider U= +/-1 case:
		mask = np.logical_or(np.abs(U-1.0) < tol, np.abs(U+1.0) < tol)
		phi = np.zeros_like(v,dtype=float)
		lbda = np.zeros_like(u,dtype=float)

		# The shortcut:
		phi[mask] = np.sign(u[mask]) * 0.25 * np.pi
		lbda[mask] = lambda0

		# The rest iteratively:
		phi_ = np.atleast_1d(0.5*np.pi - 2*np.arctan(t))
		phi_fun = lambda sinx : 0.5*np.pi - 2*np.arctan(t * np.power((1.-e*sinx)/
			                                        (1.+e*sinx), 0.5*e))
		phi_new = phi_fun(np.sin(phi_))
		itermask = np.abs(phi_new - phi_).max() > tol
		while np.any(itermask):
			phi_[itermask] = phi_new[itermask]
			phi_new = phi_fun(np.sin(phi_[itermask]))
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

		return lon.reshape(shape), lat.reshape(shape)

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

	def _scale_k(self, lon, lat):
		""""
		Class implementation of _scale_k method.
		"""
		return self.project_to_uvk(lon, lat, return_k=True)[2]


	def _is_global(self):
		return False

	def _identifier(self):
		return ("HotineObliqueMercator", (self._lon0, self._lat0, self._azimuth,
		        self._k0, self._ellipsoid, self._tolerance, self._invert_v,
		        self._invert_u))


	def __eq__(self, other):
		if isinstance(other, HotineObliqueMercator):
			return self._identifier() == other._identifier()

		return False


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

		# Compute offset of center in u:
		uc = np.sign(phi0) * A/B * np.arctan2(np.sqrt(D**2-1), np.cos(alpha_c))

		# Save constants:
		self._constants = (e, phi0, lambda_c, alpha_c, A, B, t0, D, E, F, G, gamma0, lambda0,
		                   uc)

		# Initialize Proj object, if possible:
		if _has_proj and self._use_proj:
			if azimuth == 0:
				self._proj = None
			else:
				projstr = "+proj=omerc " \
				          "+lat_0=" + str(lat0) + " " \
				          "+lonc=" + str(lon0) + " " \
				          "+alpha=" + str(azimuth) + " " \
				          "+k_0=" + str(self._k0) + " " \
				          "+ellps=" + self._ellipsoid + " +no_rot +no_off"
				self._proj = Proj(projstr)
		else:
			self._proj = None
