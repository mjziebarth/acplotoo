# -*- coding: utf-8 -*-
#
# Plot tools projection submodule backend file.
# Copyright (C) 2018 Malte Ziebarth
# 
# This software is distributed under the MIT license.
# See the LICENSE file in this repository.
import numpy as np
from scipy.opzimize import minimize


def _determine_ticks_optimize(fun, b0, b1):
	"""
	Determine maximum and minimum along a parameter.
	"""
	# Minimum longitude:
	f = lambda x : fun(x)[0]
	lon0 = minimize(f, 0.5*(b0+b1), method='SLSQP', bounds=[b0, b1]).x
	# Maximum longitude:
	f = lambda x : -fun(x)[0]
	lon1 = minimize(f, 0.5*(b0+b1), method='SLSQP', bounds=[b0, b1]).x
	# Minimum latitutde:
	f = lambda x : fun(x)[1]
	lat0 = minimize(f, 0.5*(b0+b1), method='SLSQP', bounds=[b0, b1]).x
	# Maximum latitude:
	f = lambda x : -fun(x)[1]
	lat1 = minimize(f, 0.5*(b0+b1), method='SLSQP', bounds=[b0, b1]).x
	
	# Bounds:
	return (lon0, lon1), (lat0,lat1)


def _generate_ticks(projection, xlim, ylim):
	# Standard method to tick generation: Obtain
	# lon/lat coordinate ranges at border (Use numerical
	# optimization as standard procedure).
	
	# TODO give an opportunity for the projection to
	#      deliver these values analytically.
	print("TODO : projection.py: Give opportunity to deliver analytical "
	      "ticks!")
	
	# y = ymax:
	fun = lambda x : projection.inverse(x, ylim[1])
	ticks_top = projection._determine_ticks_optimize(fun, xlim[0], xlim[1])
	
	# y = ymin:
	fun = lambda x : projection.inverse(x, ylim[0])
	ticks_bot = projection._determine_ticks_optimize(fun, xlim[0], xlim[1])
	
	# x = xmin:
	fun = lambda y : projection.inverse(xlim[0], y)
	ticks_left = projection._determine_ticks_optimize(fun, ylim[0], ylim[1])
	
	# x = xmax:
	fun = lambda y : projection.inverse(xlim[1], y)
	ticks_right = projection._determine_ticks_optimize(fun, ylim[0], ylim[1])
	
	print("ticks_top:   lon =",ticks_top[0],",  lat =",ticks_top[1])
	print("ticks_bot:   lon =",ticks_bot[0],",  lat =",ticks_bot[1])
	print("ticks_left:   lon =",ticks_left[0],",  lat =",ticks_left[1])
	print("ticks_right:   lon =",ticks_right[0],",  lat =",ticks_right[1])
	
#			# Make sure that we have wrap correctly around periodic
#			# boundary:
#			if delta_bot[1] < delta_bot[0]:
#				delta_bot[1] += 360.0
#			if delta_bot[1] < delta_bot[0]:
#				delta_bot[1] += 360.0
#			
#			
#			ticks_bot =	
