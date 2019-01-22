# -*- coding: utf-8 -*-
#
# Acplotoo Geoplot class file.
# Copyright (C) 2018 Malte Ziebarth
# 
# This software is distributed under the MIT license.
# See the LICENSE file in this repository.

from .geoplot_base.rect import Rect
from .geoplot_base.base import GeoplotBase
from .geoplot_base.backend import _ensure_coordinate_format
from .projection.projection import Projection

from warnings import warn

import numpy as np
from scipy.interpolate import interp2d
from scipy.spatial import cKDTree




# GEOPLOT:

class Geoplot(GeoplotBase):


	def __init__(self, ax, projection, limits_xy=None, gshhg_path=None,
	             which_ticks='significant', water_color='lightblue',
	             land_color='white', coast_color='black', verbose=0,
	             use_joblib=False, resize_figure=False):
		"""
		Init method.

		Required arguments:
		   ax         :
		   projection :
		   limits_xy  : [xlim, ylim]

		Optional arguments:
		   which_ticks : Determines which ticks to display at which
		                 axis. One of:
		                 'both' - draw both lon and lat ticks on all
		                     axes.
		                 'significant' - draw either lon or lat
		                     what has more ticks.
		                 'lonlat' - draw lon ticks at x- and lat ticks
		                     at y-axis
		                 'latlon' - reverse 'lonlat'
		   use_joblib  : Whether to use joblib to cache some intermediat
		                 results (e.g. coastlines). Can be useful if a
		                 lot of plots are created for the same projection.
		"""

		if not isinstance(projection,Projection):
			raise RuntimeError("The projection has to be of class 'Projection'!")

		super().__init__(ax, projection, gshhg_path, which_ticks,
		                 water_color, land_color, coast_color, verbose, use_joblib,
		                 resize_figure)

		self._gshhg_path = gshhg_path


		self._canvas = None
		self._plot_canvas = self._canvas

		# Setup configuration:
		self._box_axes = False
		self._box_axes_width = 0.1

		# If limits are given, set them:
		if limits_xy is not None:
			self._user_xlim = limits_xy[0]
			self._user_ylim = limits_xy[1]
			self._xlim = self._user_xlim
			self._ylim = self._user_ylim
			self._schedule_callback()


	def set_xlim(self, xlim):
		# TODO Sanity checks.
		self._user_xlim = xlim
		self._schedule_callback()

	def set_ylim(self, ylim):
		self._user_ylim = ylim
		self._schedule_callback()

	def coastline(self, level, water_color=None, land_color=None, coast_color=None,
	              zorder=0, **kwargs):
		"""
		Plot the coast line.
		"""
		self._coast_level=level

		if self._gshhg_path is None:
			raise RuntimeError("GSHHG not loaded!")

		if water_color is not None:
			self._water_color = water_color

		if land_color is not None:
			self._land_color = land_color

		if coast_color is not None:
			self._coast_color = coast_color

		# Remove previous coastline:
		for i in range(len(self._scheduled)):
			if self._scheduled[i] == 'coastline':
				del self._scheduled[i]
				break

		# Schedule coastline:
		self._scheduled += [['coastline', False, (level,zorder,kwargs)]]
		self._schedule_callback()

	def grid(self, on=True, grid_constant=1.0, anchor_lon=0.0, anchor_lat=0.0, **kwargs):
		"""
		Set grid on or off.
		"""
		# Save configuration:
		self._grid_on = on
		self._grid_constant = grid_constant
		self._grid_kwargs = {**self._grid_kwargs_base, **kwargs}
		self._grid_anchor = (anchor_lon, anchor_lat)

		if not "linewidth" in kwargs:
			self._grid_kwargs["linewidth"] = 0.5

		# Schedule callback:
		self._update_grid = True
		self._schedule_callback()


	def scatter(self, lon, lat, **kwargs):
		"""
		Scatter plot.
		"""
		# Schedule marker plot:
		self._scheduled += [['scatter', False, (lon, lat, kwargs)]]
		self._schedule_callback()

	def quiver(self, lon, lat, u, v, c=None, **kwargs):
		"""
		Quiver (arrow) plot.

		Required arguments:
		   lon, lat : Geodetic coordinates of arrow origins.
		   u        : Vector components in longitude direction.
		   v        : Vector components in latitude direction.

		Optional arguments:
		   c        : Vector of colors.
		              (Default: None)
		   kwargs   : Passed to matplotlib quiver.
		"""
		# Schedule quiver:
		self._scheduled += [['quiver', False, (lon, lat, u, v, c, kwargs)]]
		self._schedule_callback()

	def streamplot_projected(self, x, y, u, v, backend='matplotlib', **kwargs):
		"""
		Streamplot.

		Required arguments:
		   x, y : 1d-arrays defining the grid in projected coordinates
		   u    : 2d-grid of vector components in longitude direction.
		   v    : 2d-grid of vector components in latitude direction.

		Optional arguments:
		   backend : The backend to use. Either 'matplotlib' or 'custom'.
		             (Default: 'matplotlib')
		   kwargs  : Passed to matplotlib streamplot or PolyCollection,
		             depending on backend.

		The 'custom' backend uses a custom C++ implementation of the
		streamplot method. The result of the RKF45 integration of
		trajectories on the vector field is converted to polygons
		representing the envelope of the trajectories. The width of
		the envelope is controlled by a width-field passed as the
		'' keyword parameter.
		"""
		# Checks:
		if backend not in ['matplotlib','custom']:
			raise ValueError("Backend must be either 'matplotlib' or 'custom'")

		# Schedule streamplot:
		self._scheduled += [['streamplot', False, (x, y, u, v, backend, kwargs)]]
		self._schedule_callback()

	def streamplot(self, lon, lat, u, v, backend='matplotlib', **kwargs):
		# TODO Interpolate data to grid and do streamplot on grid!
		# TODO : Also convert start point!
		raise NotImplementedError("Geoplot.streamplot() not implemented yet.")

	def tensorfield_symmetric_2d(self, lon=None, lat=None, t1=None, t2=None, angle=None,
	                             x=None, y=None, linewidth=1.0, streamcolor='white',
	                             coastcolor='lightgray', watercolor="black",
	                             landcolor="none", cmap='inferno', coastmask=True,
	                             resample=True, n_resample=400, 
	                             resample_method='nearest',
	                             tensor=None,
	                             **kwargs):
		"""
		Plot a two-dimensional field of a symmetric two-dimensional tensor
		using streamplot. The tensor's principal axis direction is visualized
		using the streamplot line direction. The difference between first
		and second principal component is visualized using the line widths.
		The first principal component's amplitude is visualized using a
		color map applied to the lines.

		Call signatures:
		================

		1) Using numpy arrays.
		   ...

		2) Using unephy SymmetricTensorField:

		Required keyword argument:
		   tensor: A SymmetricTensorField
		"""
		if tensor is not None:
			# See if unephy is installed:
			try:
				from unephy import SymmetricTensorField, CoordinateSystem,\
				                   MapProjectionSystem
			except ImportError:
				raise RuntimeError("If 'tensor' is set, it has to be an unephy "
				                   "SymmetricTensorField instance. Could not "
				                   "import unephy!")

			# Sanity checks:
			if not (x is None and y is None and t1 is None and t2 is None and \
			        angle is None):
				raise RuntimeError("If tensor is given, all of x, y, t1, t2, "
				                   "and angle have to be None!")

			if not isinstance(tensor,SymmetricTensorField):
				raise RuntimeError("'tensor' has to be a SymmetricTensorField "
				                   "instance!")

			system = CoordinateSystem.current()
			if not isinstance(system,MapProjectionSystem) or \
			   system._projection._projection != self._projection:
				print("Coordinate system projection:",system._projection)
				print("Own projection:              ",self._projection)
				raise RuntimeError("We need to be in a projection environment fitting to "
				                   "this Geoplot's projection!")

			# Now obtain coordinates and data from the tensor field:
			with system:
				coordinates = tensor.coordinates()
				if not coordinates.is_grid():
					raise RuntimeError("Tensor needs to be a grid in plot projection "
					                   "system.")
				xy = coordinates.raw(system.default_unit())
				t1 = tensor.principal_component("first").raw(tensor.unit())
				t2 = tensor.principal_component("second").raw(tensor.unit())
				angle = tensor.principal_azimuth().raw("arcdegree")
				if tensor.usability_mask() is not None:
					mask = tensor.usability_mask()
					# Determine resulting array shape:
					ids = np.argwhere(mask)
					i0 = ids[...,0].min()
					i1 = ids[...,0].max()
					j0 = ids[...,1].min()
					j1 = ids[...,1].max()
					shape = (i1-i0+1, j1-j0+1)

					xy = xy[mask,:]
					xy = xy.reshape((*shape, xy.shape[-1]))
					if not np.all(xy[:,0,0].reshape(-1,1) == xy[:,:,0]) or \
					   not np.all(xy[0,:,1].reshape(1,-1) == xy[:,:,1]):
						raise RuntimeError("Do not have a grid after applying mask!")
					t1 = t1[mask].reshape(shape)
					t2 = t2[mask].reshape(shape)
					angle = angle[mask].reshape(shape)


				x = xy[:,0,0]
				y = xy[0,:,1]

		# Exactly one of the pairs (lon,lat) and (x,y) has to be given:
		if (lon is None) == (x is None) or (lon is None) != (lat is None) \
		or (x is None) != (y is None):
			raise ValueError("Exactly one of the pairs (lon,lat) or (x,y) have "
			                 "to be given!")

		# Determine coordinate type given:
		coordinate_type = 'projected' if x is not None else 'geographic'

		# Tensor has to be given:
		if t1 is None or t2 is None or angle is None:
			raise ValueError("Tensor has to be given!")

		# Bring coordinates into correct form:
		data = np.concatenate([t1[...,np.newaxis],t2[...,np.newaxis],
		                       angle[...,np.newaxis]], axis=-1)
		if coordinate_type == 'geographic':
			lon,lat,data = _ensure_coordinate_format(lon, lat, data, 'llg',
			                                         target_indexing='ij')
		else:
			x,y,data = _ensure_coordinate_format(x, y, data, 'llg',
			                                     target_indexing='ij')
		t1 = data[...,0]
		t2 = data[...,1]
		angle = data[...,2]

		# Resample if needed:
		if coordinate_type == 'geographic' or resample:
			# In case lon/lat is given, we need to project the coordinates.
			# Also flatten everything:
			if coordinate_type == 'geographic':
				x,y = Projection.project(self._projection, lon.flatten(), lat.flatten())
			elif x.ndim != 1:
				x = x.flatten()
				y = y.flatten()
			t1 = t1.flatten()
			t2 = t2.flatten()
			angle = angle.flatten()

			# Create source grid:
			xg,yg = np.meshgrid(x, y, indexing='ij')

			# Create target grid:
			xvals = np.linspace(x.min(),x.max(),n_resample)
			yvals = np.linspace(y.min(),y.max(),n_resample)
			xg_target,yg_target = np.meshgrid(xvals, yvals, indexing='ij')

			if resample_method == 'nearest':
				tree = cKDTree(np.concatenate([xg.flatten()[:,np.newaxis],
				                               yg.flatten()[:,np.newaxis]],
				                              axis=-1))
				ids = tree.query(np.concatenate([xg_target[...,np.newaxis],
				                                 yg_target[...,np.newaxis]],
				                                axis=-1))[1]
				t1 = t1[ids]
				t2 = t2[ids]
				angle = angle[ids]
			elif resample_method == 'spline':
				splines_angle = [interp2d(x, y, np.cos(np.deg2rad(angle))),
					             interp2d(x, y, np.sin(np.deg2rad(angle)))]
				spline_t1 = interp2d(x, y, t1)
				spline_t2 = interp2d(x, y, t2)
				t1 = spline_t1(xvals,yvals).T
				t2 = spline_t2(xvals,yvals).T
				angle = np.rad2deg(np.arctan2(splines_angle[1](xvals,yvals),
					                          splines_angle[0](xvals,yvals))).T
			else:
				raise ValueError("resample_method must be one of 'nearest' or 'spline'!")

		else:
			xvals = x
			yvals = y

		# Calculate relevant properties from tensor:
		color = t1
		width = t1-t2
		u = np.sin(np.deg2rad(angle))
		v = np.cos(np.deg2rad(angle))

		# Save keys in addition to old ones:
		kwdict = dict(kwargs)
		kwdict["linewidth"] = linewidth * (width - width.min())/(width.max() - width.min())
		kwdict["facecolor"] = streamcolor
		kwdict["quiver"] = False

		# Obtain zorder:
		zorder = kwdict.pop("zorder",1)

		# If coastline is set, we show the coastline in the same color as
		# the streamplot:
		if coastcolor is not None or watercolor is not None or landcolor is not None:
			self.coastline(self._coast_level, water_color=watercolor,
			               land_color=landcolor, coast_color=coastcolor,
			               zorder=zorder+1)

		# Call imshow:
		if coordinate_type == 'geographic':
			raise NotImplementedError("imshow not yet implemented!")
		else:
			self.imshow_projected(t1.T, [x.min(),x.max()], [y.min(),y.max()], cmap=cmap,
			                      origin='lower', zorder=zorder, coastmask=coastmask)

		# Call streamplot:
		self.streamplot_projected(xvals, yvals, u, v,
		                          backend='custom', zorder=zorder+2, **kwdict)


	def imshow_projected(self, z, xlim, ylim, **kwargs):
		"""
		Plot a field (in projected coordinates) using imshow.
		"""
		# Convert to numpy:
		if not isinstance(xlim,np.ndarray):
			xlim = np.array(xlim)
		if not isinstance(ylim,np.ndarray):
			ylim = np.array(ylim)

		# Check data limits:
		if self._data_xlim is None:
			self._data_xlim = xlim
		else:
			if xlim[0] < self._data_xlim[0]:
				self._data_xlim[0] = xlim[0]
			if xlim[1] > self._data_xlim[1]:
				self._data_xlim[1] = xlim[1]
		if self._data_ylim is None:
			self._data_ylim = ylim
		else:
			if ylim[0] < self._data_ylim[0]:
				self._data_ylim[0] = ylim[0]
			if ylim[1] > self._data_ylim[1]:
				self._data_ylim[1] = ylim[1]

		# Schedule plot:
		self._scheduled += [['imshow', False, (z, xlim,ylim,kwargs)]]
		self._schedule_callback()
