# -*- coding: utf-8 -*-
#
# Plot tools projection submodule backend file.
# Copyright (C) 2018 Malte Ziebarth
# 
# This software is distributed under the MIT license.
# See the LICENSE file in this repository.
import numpy as np
from scipy.optimize import minimize, brentq

def _generate_ticks(projection, xlim, ylim, tick_delta_degree):
	# Standard method for tick generation: Sample lon/lat
	# coordinates at border and detect increments in tick_delta_degree
	# steps.
	
	# Create sampling coordinates for axes:
	N = 100
	x_sample = [np.linspace(xlim[0],xlim[1],N), np.linspace(xlim[0],xlim[1],N),
	            xlim[0]*np.ones(N), xlim[1]*np.ones(N)]
	y_sample = [ylim[0]*np.ones(N), ylim[1]*np.ones(N),
	            np.linspace(ylim[0],ylim[1],N), np.linspace(ylim[0],ylim[1],N)]
	axes = ["bot","top","left","right"]

	
	# Now iterate over axes. Setup the variable "locations" as a two-dimensional
	# matrix of lists: locations[i][j]; i in range(4); j in range(2)
	# Index i iterates over the four axes, index j iterates over lon/lat
	locations = [[None,None] for i in range(4)]
	ret = dict()
	for i in range(4):
		# Do sampling:
		inv = projection.inverse(x_sample[i],y_sample[i])
		# Continue both for lon and lat:
		for j in range(2):
			sample = inv[j]
			# Now calculate which tick each of those sampling points
			# belongs to:
			ticks = np.floor(sample / tick_delta_degree)

			# Identify the transition points. Note: Here we assume that
			# there are no steps where more than one transition occurs.
			# If that was the case, we should sample more densely to
			# begin with - or adjust the tick_delta_degree (if we have
			# more than 100 ticks, for example, the ticks may be hard
			# to read).
			delta = ticks[1:] - ticks[:-1]
			transition = np.argwhere(delta != 0).flatten()

			# Check the assumption:
			if np.any(np.abs(delta) > 1):
				# Maybe raising an exception is too strict.
				raise RuntimeError("Tick resolution too fine!")

			# For each transition point, find the exact location, i.e. the
			# location where sample[k+1] is reached:
			if i < 2:
				# Top & bottom axes: Vary x.
				fun = lambda x,y : projection.inverse(x,ylim[i])[j] / tick_delta_degree \
				                 - y

				locations[i][j] = np.array([brentq(fun, x_sample[i][k], x_sample[i][k+1],
				                                   args=(ticks[k:k+2].max(),))
				                            for k in transition])
			else:
				# Left & right axes: Vary y.
				fun = lambda y,x : projection.inverse(xlim[i-2],y)[j] / \
				                   tick_delta_degree - x
				locations[i][j] = np.array([brentq(fun, y_sample[i][k], y_sample[i][k+1],
				                                   args=(ticks[k:k+2].max(),))
			                               for k in transition])

			print("locations[",i,"][",j,"]:",locations[i][j])

		# Save in return dictionary:
		ret[axes[i]] = locations[i]

	return ret


def _maximum_geographic_extents(projection, xlim, ylim):
	# Standard method for calculating the geographic extends
	# of limits in projection coordinates.

	# Determine numerically the maximum latitude / longitude coordinates at the
	# boundary:
	top = lambda x,j,sign   : float(sign*projection.inverse(np.atleast_1d(x),
	                                                        np.atleast_1d(ylim[1]))[j])
	bot = lambda x,j,sign   : float(sign*projection.inverse(np.atleast_1d(x),
	                                                        np.atleast_1d(ylim[0]))[j])
	left = lambda y,j,sign  : float(sign*projection.inverse(np.atleast_1d(xlim[0]),
	                                                        np.atleast_1d(y))[j])
	right = lambda y,j,sign : float(sign*projection.inverse(np.atleast_1d(xlim[1]),
	                                                        np.atleast_1d(y))[j])

	funs = [bot, top, left, right]
	x0   = [np.mean(xlim), np.mean(xlim), np.mean(ylim), np.mean(ylim)]
	bnd  = [xlim, xlim, ylim, ylim]
	lim_lon = [180.0, -180.0]
	lim_lat = [90.0, -90.0]
	method = 'TNC' # TNC seems to work!
	for i in range(4):
		# Minimize and maximize lon/lat in current axis:
		x_lonmin = minimize(funs[i], (x0[i],), args=(0,1.0,), bounds=(bnd[i],),
		                    method=method).x
		x_lonmax = minimize(funs[i], (x0[i],), args=(0,-1.0,), bounds=(bnd[i],),
		                    method=method).x
		x_latmin = minimize(funs[i], (x0[i],), args=(1,1.0,), bounds=(bnd[i],),
		                    method=method).x
		x_latmax = minimize(funs[i], (x0[i],), args=(1,-1.0,), bounds=(bnd[i],),
		                    method=method).x

		# Check if we have new extrema:
		if funs[i](x_lonmin, 0, 1.0) < lim_lon[0]:
			lim_lon[0] = funs[i](x_lonmin, 0, 1.0)

		if funs[i](x_lonmax, 0, 1.0) > lim_lon[1]:
			lim_lon[1] = funs[i](x_lonmax, 0, 1.0)

		if funs[i](x_latmin, 1, 1.0) < lim_lat[0]:
			lim_lat[0] = funs[i](x_latmin, 1, 1.0)

		if funs[i](x_latmax, 1, 1.0) > lim_lat[1]:
			lim_lat[1] = funs[i](x_latmax, 1, 1.0)

	print("_maximum_geographic_extents: {",lim_lon,",",lim_lat,"}")

	# Return extents:
	return lim_lon, lim_lat
