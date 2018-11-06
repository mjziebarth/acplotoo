# -*- coding: utf-8 -*-
#
# Plot tools projection submodule backend file.
# Copyright (C) 2018 Malte Ziebarth
# 
# This software is distributed under the MIT license.
# See the LICENSE file in this repository.
import numpy as np
from scipy.optimize import minimize, brentq


#def _determine_ticks_optimize(fun, b0, b1):
#	"""
#	Determine maximum and minimum along a parameter.
#	"""
#	# Minimum longitude:
#	f = lambda x : fun(x)[0]
#	lon0 = minimize(f, 0.5*(b0+b1), method='SLSQP', bounds=[b0, b1]).x
#	# Maximum longitude:
#	f = lambda x : -fun(x)[0]
#	lon1 = minimize(f, 0.5*(b0+b1), method='SLSQP', bounds=[b0, b1]).x
#	# Minimum latitutde:
#	f = lambda x : fun(x)[1]
#	lat0 = minimize(f, 0.5*(b0+b1), method='SLSQP', bounds=[b0, b1]).x
#	# Maximum latitude:
#	f = lambda x : -fun(x)[1]
#	lat1 = minimize(f, 0.5*(b0+b1), method='SLSQP', bounds=[b0, b1]).x
#	
#	# Bounds:
#	return (lon0, lon1), (lat0,lat1)


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
				# sample = projection.inverse(x_sample[i],y_sample[i])[j]
				for k in transition:
					print("k         =",k)
					print("   fun(1) =",fun(x_sample[i][k],ticks[k+1]))
					print("   fun(2) =",fun(x_sample[i][k+1],ticks[k+1]))
					print("ticks[k]  =",ticks[k])
					print("ticks[k+1]=",ticks[k+1])
				
				locations[i][j] = [brentq(fun, x_sample[i][k], x_sample[i][k+1],
				                          args=(ticks[k+1]),)
				                          )[0]
			                       for k in transition]
			else:
				# Left & right axes: Vary y.
				fun = lambda y,x : projection.inverse(xlim[i-2],y)[j] / \
				                   tick_delta_degree - x
				locations[i][j] = [brentq(fun, y_sample[i][k], y_sample[i][k+1],
				                          args=(np.floor(sample[k+1]/tick_delta_degree),)
				                          )[0]
			                       for k in transition]

		# Save in return dictionary:
		ret[axes[i]] = locations[i]
	
	return ret
