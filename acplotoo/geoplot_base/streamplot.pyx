# -*- coding: utf-8 -*-
# distutils : language = c++
#
# Acplotoo streamplot cython implementation.
# Copyright (C) 2018 Malte Ziebarth
# 
# This software is distributed under the MIT license.
# See the LICENSE file in this repository.

# Interface to python world:
import numpy as np

# Cython imports:
import cython
cimport numpy as np
from libcpp cimport bool
from libcpp.vector cimport vector
from libcpp.list cimport list as clist
from libcpp.utility cimport pair
from libcpp.memory cimport unique_ptr, make_unique
from cython.operator cimport dereference, preincrement
np.import_array()

####################################################################################
#                                    C++ types                                     #
####################################################################################

# The classes. Define only what we need:
cdef extern from "streamplot.hpp" namespace "acplotoo":
	ctypedef pair[size_t,size_t] index_t

	cdef cppclass Grid:
		Grid(double x0, double x1, double y0, double y1, size_t nx, size_t ny) except +

	cdef cppclass Vectorfield:
		Vectorfield(size_t M, size_t N) except +
		pair[double,double]& operator[](const index_t& index)


# The method:
	cdef clist[vector[pair[double,double]]]\
	         streamplot_polygons(
		                    const vector[pair[double,double]]& start,
		                    const Vectorfield& velocity, const Grid& velocity_grid,
		                    const pair[size_t,size_t]& mask_size,
		                    double min_len, double max_len, size_t max_steps,
		                    bool forward, bool backward, double tol)



def streamplot_calculate_polygons(np.ndarray[double, ndim=2] x, 
        np.ndarray[double, ndim=2] y, np.ndarray[double, ndim=2] vx,
        np.ndarray[double, ndim=2] vy,
        double min_len, double max_len,
        np.ndarray[double, ndim=1] start_x, np.ndarray[double, ndim=1] start_y,
        bool forward=True, bool backward=True, size_t max_steps=1000, 
        double tolerance = 1e-10):
	# Ensure that shapes are consistent:
	assert x.shape[0] == y.shape[0] and x.shape[1] == y.shape[1]
	assert x.shape[0] == vx.shape[0] and x.shape[1] == vx.shape[1]
	assert x.shape[0] == vy.shape[0] and x.shape[1] == vy.shape[1]
	assert start_x.size() ==  start_y.size()

	# Preparation:
	cdef size_t M = x.shape[0]
	cdef size_t N = x.shape[1]
	cdef double xmin = x.min()
	cdef double xmax = x.max()
	cdef double ymin = y.min()
	cdef double ymax = y.max()

	# Create grids and vector fields:
	cdef unique_ptr[Vectorfield] velocity
	cdef unique_ptr[Grid] velocity_grid
	velocity.reset(new Vectorfield(M, N))
	velocity_grid.reset(new Grid(xmin, xmax, ymin, ymax, M, N))

	if not velocity or not velocity_grid:
		raise RuntimeError("Error: Could not localize Vectorfield or Grid!")

	# Fill the velocity field:
	cdef size_t i
	cdef pair[size_t,size_t] index
	if x[0,0] != x[1,0]:
		for i in range(M):
			for j in range(N):
				index.first = i
				index.second = j
				dereference(velocity)[index].first = vx[i,j]
				dereference(velocity)[index].second = vy[i,j]
	elif x[0,0] != x[0,1]:
		for j in range(N):
			for i in range(M):
				index.first = i
				index.second = j
				dereference(velocity)[index].first = vx[j,i]
				dereference(velocity)[index].second = vy[j,i]
	else:
		raise RuntimeError("Could not determine grid setup!")

	# Mask size:
	cdef pair[size_t,size_t] mask_size
	mask_size.first = 30
	mask_size.second = 30

	# Fill start values:
	cdef vector[pair[double,double]] start
	start.resize(start_x.size())
	for i in range(start_x.size()):
		start[i].first = start_x[i]
		start[i].second = start_y[i]

	# Call algorithm:
	cdef clist[vector[pair[double,double]]] res = \
	    streamplot_polygons(start, dereference(velocity), dereference(velocity_grid), 
	                        mask_size, min_len, max_len, max_steps, forward, backward,
	                        tolerance,)
