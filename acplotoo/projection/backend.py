# -*- coding: utf-8 -*-
#
# Acplotoo projection submodule backend file.
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

	# TODO : This method could be improved in terms of sophisticatedness.

	axes = ["bot","top","left","right"]
	has_results = [False,False,False,False]
	locations = [[None,None] for i in range(4)]
	N = [100,100,100,100]
	ret = dict()
	while not all(has_results):
		# Create sampling coordinates for axes:
		x_sample = [np.linspace(xlim[0],xlim[1],N[0]),
		            np.linspace(xlim[0],xlim[1],N[1]),
		            xlim[0]*np.ones(N[2]), xlim[1]*np.ones(N[3])]
		y_sample = [ylim[0]*np.ones(N[0]), ylim[1]*np.ones(N[1]),
		            np.linspace(ylim[0],ylim[1],N[2]),
		            np.linspace(ylim[0],ylim[1],N[3])]


		# Now iterate over axes. Setup the variable "locations" as a two-dimensional
		# matrix of lists: locations[i][j]; i in range(4); j in range(2)
		# Index i iterates over the four axes, index j iterates over lon/lat
		for i in range(4):
			if has_results[i]:
				continue

			# Do sampling:
			inv = projection.inverse(x_sample[i],y_sample[i])
			# Continue both for lon and lat:
			results_i = 0
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
					ids = np.argwhere(np.abs(delta) > 1)
					if np.any(np.sign(ticks[ids]) == np.sign(ticks[ids+1])):
						# Assumption failed. Increase sampling points.
						# (this is an inefficient hotfix)
						N[i] *= 10
						if N[i] >= 1e5:
							# Maybe raising an exception is too strict.
							raise RuntimeError("Tick resolution too fine!")

						break

				# For each transition point, find the exact location, i.e. the
				# location where sample[k+1] is reached:
				if i < 2:
					# Top & bottom axes: Vary x.
					fun = lambda x,y : projection.inverse(x,ylim[i])[j] / \
					                   tick_delta_degree - y

					locations[i][j] = np.array([brentq(fun, x_sample[i][k],
					                                   x_sample[i][k+1],
					                                   args=(ticks[k:k+2].max(),))
					                            for k in transition])
				else:
					# Left & right axes: Vary y.
					fun = lambda y,x : projection.inverse(xlim[i-2],y)[j] / \
					                   tick_delta_degree - x
					locations[i][j] = np.array([brentq(fun, y_sample[i][k],
					                                   y_sample[i][k+1],
					                                   args=(ticks[k:k+2].max(),))
				                               for k in transition])

				results_i += 1

			# Save in return dictionary:
			ret[axes[i]] = locations[i]
			has_results[i] = results_i == 2

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

	# Return extents:
	return lim_lon, lim_lat


def _unit_vector(lon, lat, projection, cardinal_direction,
                 delta=1e-8):
	# Determine the gradient of x/y in lon/lat directions
	# and of those construct the unit vector!
	if cardinal_direction == "north":
		grad = (projection.project(lon,lat), projection.project(lon,lat+delta))
	elif cardinal_direction == "east":
		grad = (projection.project(lon,lat), projection.project(lon+delta,lat))
	else:
		raise RuntimeError("Only cardinal directions 'north' and 'east' supported!")

	# Calculate deltas:
	dx = grad[1][0] - grad[0][0]
	dy = grad[1][1] - grad[0][1]

	# Unnormalized vector:
	vec = np.concatenate((dx[...,np.newaxis], dy[...,np.newaxis]),-1) / delta

	# Calculate and return gradients:
	return vec / np.linalg.norm(vec,axis=-1)[...,np.newaxis]
