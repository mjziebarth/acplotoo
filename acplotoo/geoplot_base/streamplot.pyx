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

# Woraround courtesy of user985030
cdef extern from "<array>" namespace "std" nogil:
	cdef cppclass ArrayDbl4 "std::array<double,4>":
		ArrayDbl4() except+
		double& operator[](size_t)

# The classes. Define only what we need:
cdef extern from "streamplot.hpp" namespace "acplotoo":
	ctypedef pair[size_t,size_t] index_t

	cdef cppclass Grid:
		Grid(double x0, double x1, double y0, double y1, size_t nx, size_t ny) except +
		pair[double,double]& operator[](const index_t& index)

	cdef cppclass Vectorfield:
		Vectorfield(size_t M, size_t N) except +
		pair[double,double]& operator[](const index_t& index)

	cdef cppclass Scalarfield:
		Scalarfield(size_t M, size_t N) except +
		double& operator[](const index_t& index)

# The method:
	cdef pair[clist[vector[pair[double,double]]],
	          vector[ArrayDbl4]]\
	         streamplot_polygons(
		                    const vector[pair[double,double]]& start,
		                    const Vectorfield& velocity, const Grid& velocity_grid,
		                    const Scalarfield& width_grid,
		                    double min_len, double max_len, double step_len_min,
		                    double arrow_head_step, double collision_radius,
		                    size_t max_steps, bool forward,
		                    bool backward, double tol, size_t tile_history_size) nogil



def _streamplot_calculate_polygons(np.ndarray[double, ndim=2] x, 
         np.ndarray[double, ndim=2] y, np.ndarray[double, ndim=2] vx,
         np.ndarray[double, ndim=2] vy, np.ndarray[double, ndim=2] widths,
         double min_len, double max_len, double step_len_min,
         double arrow_head_step, double collision_radius,
         np.ndarray[double, ndim=1] start_x, np.ndarray[double, ndim=1] start_y,
         bool forward=True, bool backward=True, size_t max_steps=10000, 
         double tolerance = 1e-5, size_t tile_history_size=10):
	# Ensure that shapes are consistent:
	assert x.shape[0] == y.shape[0] and x.shape[1] == y.shape[1]
	assert x.shape[0] == vx.shape[0] and x.shape[1] == vx.shape[1]
	assert x.shape[0] == vy.shape[0] and x.shape[1] == vy.shape[1]
	assert x.shape[0] == widths.shape[0] and x.shape[1] == widths.shape[1]
	assert start_x.size ==  start_y.size

	# Preparation:
	cdef size_t M
	cdef size_t N
	cdef double xmin = x.min()
	cdef double xmax = x.max()
	cdef double ymin = y.min()
	cdef double ymax = y.max()

	if x[0,0] != x[1,0]:
		M = x.shape[0]
		N = x.shape[1]
	else:
		M = x.shape[1]
		N = x.shape[0]

	# Create grids and vector fields:
	cdef unique_ptr[Vectorfield] velocity
	cdef unique_ptr[Grid] velocity_grid
	cdef unique_ptr[Scalarfield] width_field
	velocity.reset(new Vectorfield(M, N))
	velocity_grid.reset(new Grid(xmin, xmax, ymin, ymax, M, N))
	width_field.reset(new Scalarfield(M, N))

	if not velocity or not velocity_grid:
		raise RuntimeError("Error: Could not localize Vectorfield or Grid!")

	# Fill the velocity field:
	cdef size_t i, j
	cdef double vx_,vy_
	cdef pair[size_t,size_t] index
	if x[0,0] != x[1,0]:
		for i in range(M):
			for j in range(N):
				index.first = i
				index.second = j
				dereference(velocity)[index].first  = vx[i,j]
				dereference(velocity)[index].second = vy[i,j]
				dereference(width_field)[index] = widths[i,j]
	elif x[0,0] != x[0,1]:
		for j in range(N):
			for i in range(M):
				index.first = i
				index.second = j
				dereference(velocity)[index].first  = vx[j,i]
				dereference(velocity)[index].second = vy[j,i]
				dereference(width_field)[index] = widths[j,i]
	else:
		raise RuntimeError("Could not determine grid setup!")

	# Fill start values:
	cdef vector[pair[double,double]] start
	cdef size_t n_start = start_x.size
	start.resize(n_start)
	for i in range(n_start):
		start[i].first = start_x[i]
		start[i].second = start_y[i]

	# Call algorithm:
	cdef pair[clist[vector[pair[double,double]]],
	          vector[ArrayDbl4]] res
	with nogil:
		res = \
		    streamplot_polygons(start, dereference(velocity), dereference(velocity_grid),
		                        dereference(width_field),
		                        min_len, max_len, step_len_min, arrow_head_step,
		                        collision_radius,
		                        max_steps, forward, backward, tolerance,
		                        tile_history_size)

	# Return a python list:
	polygons = list()
	cdef np.ndarray[ndim=2,dtype=double] tmp
	cdef clist[vector[pair[double,double]]].iterator it = res.first.begin()
	for i in range(res.first.size()):
		tmp = np.zeros([dereference(it).size(),2])
		for j in range(dereference(it).size()):
			tmp[j,0] = dereference(it)[j].first
			tmp[j,1] = dereference(it)[j].second
		polygons += [tmp]
		preincrement(it)

	cdef vector[ArrayDbl4].iterator it2 = res.second.begin()
	cdef np.ndarray[ndim=2,dtype=double] arrows = np.zeros((res.second.size(),4))
	for i in range(res.second.size()):
		for j in range(4):
			arrows[i,j] = dereference(it2)[j]
		preincrement(it2)

	return polygons, arrows
