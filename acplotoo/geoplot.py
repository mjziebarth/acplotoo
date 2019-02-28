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
	             use_joblib=False, resize_figure=False, axes_margin_pt=5.0):
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
		   axes_margin_pt : (Default: 5.0)
		"""

		if not isinstance(projection,Projection):
			raise RuntimeError("The projection has to be of class 'Projection'!")

		super().__init__(ax, projection, gshhg_path, which_ticks,
		                 water_color, land_color, coast_color, verbose, use_joblib,
		                 resize_figure, axes_margin_pt)

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
		# Some sanity checks:
		if isinstance(on,bool):
			pass
		elif on in ('on','off'):
			on = (on == 'on')
		else:
			raise ValueError("First argument or 'on' keyword parameter has to be "
			                 "boolean or one of ('on','off')!")

		if not isinstance(grid_constant, float) or isinstance(grid_constant,int) \
		   or grid_constant <= 0:
			raise TypeError("grid_constant has to be a positive number!")

		if not isinstance(anchor_lon, float) or isinstance(anchor_lon,int):
			raise TypeError("anchor_lon has to be a number!")

		if not isinstance(anchor_lat, float) or isinstance(anchor_lat,int):
			raise TypeError("anchor_lat has to be a number!")

		# Save configuration:
		self._grid_on = on
		self._grid_constant = grid_constant
		self._grid_kwargs = {**self._grid_kwargs_base, **kwargs}
		self._grid_anchor = (((anchor_lon+180.0) % 360.0) - 180.0,
		                     ((anchor_lat+90.0) % 180.0) - 90.0)

		if not "linewidth" in kwargs:
			self._grid_kwargs["linewidth"] = 0.5

		# Schedule callback:
		self._update_grid = True
		self._schedule_callback()


	def freeze_axes(self):
		"""
		Freeze the axes limits.
		"""
		if self._xlim is not None:
			self._user_xlim = self._xlim
		else:
			self._user_xlim = self._data_xlim

		if self._ylim is not None:
			self._user_ylim = self._ylim
		else:
			self._user_ylim = self._ylim


	def filter_ticks(self, all_axes=None, left=None, right=None, top=None,
	                 bot=None):
		"""
		Refine the auto-calculated ticks. Dictionaries can be
		passed that contain all ticks allowed. Intended to be used
		for fine-tuning the axes ticks.

		Optional keyword parameters:
		   all_axes              : Dictionary of ticks that should
		                           be selected on all axes.
		                           (Default: None)
		   left, right, top, bot : Dictionaries of ticks that should
		                           be selected on the respective
		                           axis.
		                           (Default: None)

		Each dictionary may contain longitude and latitude ticks,
		under the respective keys "lon" and "lat". The ticks shall
		be given in arcdegrees in ranges [-180,180] and [-90,90],
		respectively. Only the ticks given in the dictionary are
		accepted from the automatically generated ticks.
		The specific sets of keys given by left, right, top, and bot
		are each combined with the general set all_axes.

		Example dictionary:
		   left={"lon" : (-113,), "lat" : (50, 48)}
		"""
		if left is None:
			left = dict()
		if right is None:
			right = dict()
		if top is None:
			top = dict()
		if bot is None:
			bot = dict()
		if all_axes is None:
			all_axes = dict()
		if not isinstance(left, dict) or not isinstance(right,dict) \
		    or not isinstance(top,dict) or not isinstance(bot,dict) \
		    or not isinstance(all_axes,dict):
			raise TypeError("All arguments to filter_ticks have to be "
			                "dictionaries!")
		tick_dicts = ({**all_axes, **bot}, {**all_axes, **top},
		              {**all_axes, **left}, {**all_axes, **right})
		if any(any(key not in ("lon","lat") or not (isinstance(D[key],tuple)
		                                            or isinstance(D[key],list))
		                                    or not all(isinstance(i,int) for i in D[key])
		           for key in D) for D in tick_dicts):
			raise RuntimeError("Dictionaries have to consist of keys 'lon' or 'lat' "
			                   "and sequences of int!")

		# Save tick dicts:
		if tick_dicts != self._tick_filters:
			self._tick_filters = tick_dicts
			self._update_axes = True
			self._schedule_callback()


	def scatter(self, *args, **kwargs):
		"""
		Scatter plot. Has two call signatures:

		scatter(lon, lat):
		    Call scatter plot with two arrays of longitude
		    and latitude coordinates (in arcdegrees).

		scatter(point_set):
		    Call with a single argument, a unephy PointSet.

		Optional keyword arguments:
		    **kwargs : Matplotlib scatter kwargs.
		"""
		# Schedule marker plot:
		if len(args) == 2:
			lon = args[0]
			lat = args[1]
			if not isinstance(lon,np.ndarray) or not isinstance(lat, np.ndarray):
				lon = np.array(lon)
				lat = np.array(lat)
		elif len(args) == 1:
			# Make sure unephy is installed:
			try:
				from unephy import PointSet, GeographicSystem
			except ImportError:
				raise RuntimeError("If 'tensor' is set, it has to be an unephy "
				                   "SymmetricTensorField instance. Could not "
				                   "import unephy!")

			pointset = args[0]
			if not isinstance(pointset, PointSet):
				raise TypeError("point_set has to be a unephy PointSet!")
			with GeographicSystem():
				coords = pointset.coordinates().raw("arcdegree")
				lon = coords[...,0].reshape(-1)
				lat = coords[...,1].reshape(-1)

		else:
			raise RuntimeError("Invalid call signature!")

		self._add_data(lon=lon, lat=lat)
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


	def orientations(self, *args, error=None, symmetric=False,
	                 c=None, markerfill=False, **kwargs):
		"""
		Indicate a set of located orientations given by (lon, lat, azimuth)
		tuples.

		Call signature 1 arguments:

		   lon, lat  : Geodetic coordinates of sample points.
		   azimuth   : Azimuth values at sample points.

		Call signature 2 arguments:

		   azimuth   : unephy SpatialDataSet of azimuth values.

		Optional arguments:
		   error      : Measure of uncertainty of azimuth values.
		                (Default: None)
		   symmetric  : Whether to draw an arrow head (False) or not (True).
		                (Default: False)
		   c          : Vector of arrow / bar colors.
		                (Default: None)
		   markerfill : Whether to fill the markers.
		                (Default: False)
		   kwargs     : Passed to matplotlib quiver.
		"""
		if len(args) == 3:
			# Call signature 1.
			lon = args[0]
			lat = args[1]
			azimuth = args[2]
			if not isinstance(lon,np.ndarray):
				lon = np.array(lon)
			if not isinstance(lat,np.ndarray):
				lat = np.array(lat)
			if not isinstance(azimuth, np.ndarray):
				azimuth = np.array(azimuth)

		elif len(args) == 1:
			# Call signature 2.
			# Sanity checks:
			from unephy import SpatialDataSet, CoordinateSystem,\
			                   MapProjectionSystem, GeographicSystem
			data = args[0]
			assert isinstance(data,SpatialDataSet)
			assert data.data_shape() == (1,)

			# Obtain raw data:
			with GeographicSystem():
				coords = data.coordinates().raw("arcdegree")
			lon = coords[...,0].reshape(-1)
			lat = coords[...,1].reshape(-1)
			azimuth = data.raw("arcdegree").reshape(-1)

			system = CoordinateSystem.current()
			if system._projection._projection != self._projection:
				raise RuntimeError("We need to be in a projection environment fitting to "
				                   "this Geoplot's projection!")

		else:
			raise RuntimeError("Invalid call signature!")

		# TODO some more sanity checks on keyword parameters?

		# Adjust and fill the keyword arguments:
		kw = {"units" : "dots", "pivot" : "mid", "zorder" : 1, "color" : "k",
		      "width" : 1.0}
		if symmetric:
			kw["headwidth"] = 0.0
			kw["headlength"] = 0.0
			kw["headaxislength"] = 0.0
		if "linewidth" in kwargs:
			kw["width"] = kwargs["linewidth"]

		kw = {**kw, **kwargs}

		# Create the quiver plot:
		az_rad = np.deg2rad(azimuth)
		self.quiver(lon, lat, np.sin(az_rad), np.cos(az_rad), c=c, **kw)

		# Uncertainties:
		if error is not None:
			raise NotImplementedError("Error in orientations plot not yet implemented!")

		else:
			# Do scatter plot to mark origin of quivers:
			# Empirical equation for s:
			s = 0.8 * np.sqrt(kw["width"] / 0.25)
			self.scatter(lon, lat, c=kw["color"] if markerfill else "none",
			             edgecolor=kw["color"], zorder=kw["zorder"], s=s,
			             linewidth=kw["width"])


	def streamplot_projected(self, x, y, u, v, backend='matplotlib',
	                         show_border=False, **kwargs):
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
		if show_border:
			if backend != 'custom':
				raise RuntimeError("Border only possible in 'custom' mode!")

		# Schedule streamplot:
		self._scheduled += [['streamplot', False, (x, y, u, v, backend, show_border,
		                                           kwargs)]]
		self._schedule_callback()


	def streamplot(self, lon, lat, u, v, backend='matplotlib', **kwargs):
		# TODO Interpolate data to grid and do streamplot on grid!
		# TODO : Also convert start point!
		raise NotImplementedError("Geoplot.streamplot() not implemented yet.")


	def scalar_field(self, lon=None, lat=None, scalar=None, x=None, y=None,
	                 coastcolor='lightgray', watercolor="black", landcolor="none",
	                 cmap='default', coastmask=True, resample=True, n_resample=400, 
	                 resample_method='nearest', show_border=True,
	                 colorbar='horizontal', cax=None, cbar_label=None, **kwargs):
		"""
		Plot a two-dimensional scalar field using imshow.

		Call signatures:
		================

		1) Using numpy arrays.
		   ...

		2) Using unephy ScalarField:

		Required keyword argument:
		   scalar: A ScalarField

		Optional keyword arguments:
		   coastcolor     :
		                    (Default: 'lightgray')
		   watercolor     :
		                    (Default: 'black')
		   landcolor      :
		                    (Default: 'none')
		   cmap           :
		                    (Default: 'default')
		   coastmask      : Whether to clip the image at the land borders.
		                    (Default: True)
		   resample       : Whether to resample the field.
		                    (Default: True)
		   n_resample     :
		                    (Default: 400)
		   resample_method: One of 'nearest' and 'spline'
		                    (Default: 'nearest')
		   show_border    :
		                    (Default: True)
		   colorbar       :
		                    (Default: None)
		"""
		# Try to import unephy:
		try:
			from unephy import ScalarField, CoordinateSystem,\
			                   MapProjectionSystem
			has_unephy = True
		except:
			has_unephy = False


		if has_unephy and isinstance(scalar,ScalarField):
			# Sanity checks:
			if not (x is None and y is None and lon is None and lat is None):
				raise RuntimeError("If scalar is given, all of x, y, lon, and lat "
				                   " have to be None!")

			system = CoordinateSystem.current()
			if not isinstance(system,MapProjectionSystem) or \
			   system._projection._projection != self._projection:
				raise RuntimeError("We need to be in a projection environment fitting to "
				                   "this Geoplot's projection!")

			# Now obtain coordinates and data from the scalar field:
			with system:
				coordinates = scalar.coordinates()
				if not coordinates.is_grid():
					raise RuntimeError("Scalar needs to be a grid in plot projection "
					                   "system.")
				xy = coordinates.raw(system.default_unit())
				scalar_ = scalar.raw(scalar.unit())
				if scalar.usability_mask() is not None:
					mask = scalar.usability_mask()
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
					scalar_ = scalar_[mask].reshape(shape)


				x = xy[:,0,0]
				y = xy[0,:,1]

		else:
			assert isinstance(scalar_, np.ndarray)
			scalar_ = scalar

		# Exactly one of the pairs (lon,lat) and (x,y) has to be given:
		if (lon is None) == (x is None) or (lon is None) != (lat is None) \
		or (x is None) != (y is None):
			raise ValueError("Exactly one of the pairs (lon,lat) or (x,y) have "
			                 "to be given!")

		# Determine coordinate type given:
		coordinate_type = 'projected' if x is not None else 'geographic'

		# Scalar has to be given:
		if scalar_ is None:
			raise ValueError("Scalar has to be given!")

		# Bring coordinates into correct form:
		if coordinate_type == 'geographic':
			lon,lat,scalar_ = _ensure_coordinate_format(lon, lat, scalar_, 'llg',
			                                            target_indexing='ij')
		else:
			x,y,scalar_ = _ensure_coordinate_format(x, y, scalar_, 'llg',
			                                        target_indexing='ij')

		# Resample if needed:
		if coordinate_type == 'geographic' or resample:
			# In case lon/lat is given, we need to project the coordinates.
			# Also flatten everything:
			if coordinate_type == 'geographic':
				x,y = Projection.project(self._projection, lon.flatten(), lat.flatten())
			elif x.ndim != 1:
				x = x.flatten()
				y = y.flatten()
			scalar_ = scalar_.flatten()

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
				scalar_ = scalar_[ids]
			elif resample_method == 'spline':
				spline = interp2d(x, y, scalar_)
				scalar_ = spline(xvals,yvals).T
			else:
				raise ValueError("resample_method must be one of 'nearest' or 'spline'!")

		else:
			xvals = x
			yvals = y

		# Default color map:
		if cmap == 'default':
			cmap = 'inferno'

		# Save keys in addition to old ones:
		kwdict = dict(kwargs)

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
			self.imshow_projected(scalar_.T, [x.min(),x.max()], [y.min(),y.max()],
			                      cmap=cmap, origin='lower', zorder=zorder,
			                      coastmask=coastmask, colorbar=colorbar,
			                      cax=cax, cbar_label=cbar_label,
			                      **kwdict)


	def tensorfield_symmetric_2d(self, lon=None, lat=None, t1=None, t2=None, angle=None,
	                             x=None, y=None, linewidth=1.0, streamcolor='white',
	                             coastcolor='lightgray', watercolor="black",
	                             landcolor="none", cmap='default', coastmask=True,
	                             resample=True, n_resample=400, 
	                             resample_method='nearest',
	                             tensor=None, colormode='max',
	                             direction='max', colorbar='horizontal',
	                             cax=None, cbar_label=None,
	                             thickness='difference', show_border=True,
	                             **kwargs):
		"""
		Plot a two-dimensional field of a symmetric two-dimensional tensor
		using streamplot. The tensor's principal axis direction is visualized
		using the streamplot line direction. The difference between first
		and second principal component is visualized using the line widths.
		The first principal component's amplitude is visualized using a
		background color map.

		Call signatures:
		================

		1) Using numpy arrays.
		   ...

		2) Using unephy SymmetricTensorField:

		Required keyword argument:
		   tensor: A SymmetricTensorField

		Optional keyword arguments:
		   colormode : Controls how the background color is determined. One of
		               'max'    : Amplitude of the maximum principal component.
		               'maxabs' : Amplitude of the maximum absolute principal
		                          component.
		               'sum'    : First invariant: Sum of the two principal
		                          components.
		   direction : Determines by which criterion the direction of the
		               integrated vector field is chosen. One of
		               'max'    : The direction of the biggest principal component
		                          is chosen.
		               'min'    : The direction of the smallest principal
		                          component is chosen.
		               'maxabs' : The direction of biggest or smallest principal
		                          component is chosen by maximum absolute
		                          magnitude.
		   thickness : Determines by which criterion the stream line thickness
		               is chosen. One of:
		               'difference' : The amplitude of the difference between
		                              biggest and smallest principal component
		                              determines the thickness.
		               'abs'        : The absolute value of the component
		                              determining the direction is chosen.
		   colorbar  : One of 'horizontal' or 'vertical'.
		   cax       : None or an axis on which to draw the color bar.
		               (Default: None)
		"""
		if not colormode in ['max','min','maxabs','sum','angle','second_moment']:
			raise ValueError("colormode must be one of 'max', 'min', 'maxabs', "
			                 "'sum', 'second_moment', or 'angle'.")
		if not direction in ['max','min','maxabs']:
			raise ValueError("direction must be one of 'max', 'min', or 'maxabs'.")
		if not thickness in ['difference','abs']:
			raise ValueError("thickness must be one of 'difference', or 'abs'.")

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
			        angle is None and lon is None and lat is None):
				raise RuntimeError("If tensor is given, all of x, y, lon, lat, "
				                   "t1, t2, and angle have to be None!")

			if not isinstance(tensor,SymmetricTensorField):
				raise RuntimeError("'tensor' has to be a SymmetricTensorField "
				                   "instance!")

			system = CoordinateSystem.current()
			if not isinstance(system,MapProjectionSystem) or \
			   system._projection._projection != self._projection:
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
		assert np.all(t1 >= t2)
		if colormode == 'max':
			color = t1
		elif colormode == 'min':
			color = t2
		elif colormode == 'maxabs':
			# At every point, choose the principal component with bigget
			# absolute magnitude:
			color = t1.copy()
			mask = np.abs(t2) > np.abs(t1)
			color[mask] = t2[mask]
		elif colormode == 'sum':
			color = t1+t2
		elif colormode == 'second_moment':
			color = t1*t2
		elif colormode == 'angle':
			color = angle.view()
			if cmap == 'default':
				# Use a cyclic color map.
				# cmocean cm_phase seems to be a more pleasing color
				# map than the only matplotlib one, hsv.
				try:
					from cmocean.cm import phase as cm_phase
					cmap = cm_phase
				except:
					cmap = 'hsv'

		if direction == 'max':
			pass
		elif direction == 'min':
			# Since t1 is the biggest principal component, we have to
			# rotate the direction by 90°:
			angle += 90.0
		elif direction == 'maxabs':
			# Identify all ids where |t2| > |t1|:
			mask = np.abs(t2) > np.abs(t1)
			angle[mask] += 90.0

		if thickness == 'difference':
			width = t1-t2
		elif thickness == 'abs':
			width = t1

		u = np.sin(np.deg2rad(angle))
		v = np.cos(np.deg2rad(angle))

		# Default color map:
		if cmap == 'default':
			if np.all(color >= 0.0) or np.all(color <= 0.0):
				if np.all(color >= 0.0):
					cmap = 'inferno'
					vmin = 0.0
					vmax = color.max()
				else:
					cmap = 'inferno_r'
					vmax = 0.0
					vmin = color.min()
			else:
				cmap = 'seismic'
				vmax = np.abs(color).max()
				vmin = -vmax

			# Produce default color limits:
			if 'vmax' not in kwargs:
				kwargs['vmax'] = vmax
			if 'vmin' not in kwargs:
				kwargs['vmin'] = vmin

		# Keys for imshow:
		imshow_kwargs = {}
		if 'vmin' in kwargs:
			imshow_kwargs["vmin"] = kwargs.pop("vmin",0)
		if 'vmax' in kwargs:
			imshow_kwargs["vmax"] = kwargs.pop("vmax",0)

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
			self.imshow_projected(color.T, [x.min(),x.max()], [y.min(),y.max()],
			                      cmap=cmap, origin='lower', zorder=zorder,
			                      colorbar=colorbar, cax=cax, cbar_label=cbar_label,
			                      coastmask=coastmask, **imshow_kwargs)

		# Call streamplot:
		self.streamplot_projected(xvals, yvals, u, v, backend='custom',
		                          show_border=show_border, zorder=zorder+2,
		                          **kwdict)


	def imshow_projected(self, z, xlim, ylim, colorbar='horizontal', cax=None,
	                     cbar_label=None, **kwargs):
		"""
		Plot a field (in projected coordinates) using imshow.
		"""
		# Convert to numpy:
		if not isinstance(xlim,np.ndarray):
			xlim = np.array(xlim)
		if not isinstance(ylim,np.ndarray):
			ylim = np.array(ylim)

		# Check data limits:
		self._add_data(x=xlim, y=ylim)

		# Schedule plot:
		self._scheduled += [['imshow', False, (z, xlim, ylim, colorbar, cax, cbar_label,
		                                       kwargs)]]
		self._schedule_callback()


	def polygon(self, x=None, y=None, lon=None, lat=None, **kwargs):
		"""
		Plot a polygon.
		"""
		# Sanity checks:
		if (x is None) == (lon is None):
			if x is None:
				raise ValueError("Coordinates have to be given!")
			else:
				raise ValueError("Coordinates have to be given in projected "
				                 "xor geographic coordinates!")

		if (x is None) != (y is None):
			raise ValueError("Both of x and y have to be given!")
		if (lon is None) != (lat is None):
			raise ValueError("Both longitudes and latitudes have to be given!")

		# Convert to numpy:
		if x is not None and not isinstance(x,np.ndarray):
			x = np.array(x)
		if y is not None and not isinstance(y,np.ndarray):
			y = np.array(y)
		if lon is not None and not isinstance(lon,np.ndarray):
			lon = np.array(lon)
		if lat is not None and not isinstance(lat,np.ndarray):
			lat = np.array(lat)

		# Check data limits:
		self._add_data(x=x, y=y, lon=lon, lat=lat)

		# Schedule plot:
		self._scheduled +=[['polygon', False, (x, y, lon, lat, kwargs)]]
		self._schedule_callback()
