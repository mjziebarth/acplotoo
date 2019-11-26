# -*- coding: utf-8 -*-
#
# Acplotoo Geoplot class file.
# Copyright (C) 2018 Malte Ziebarth
# 
# This software is distributed under the MIT license.
# See the LICENSE file in this repository.

from .geoplot_base.rect import Rect
from .geoplot_base.base import GeoplotBase
from .geoplot_base.backend import _ensure_coordinate_format, has_unephy, \
                                  is_land_backend
from .geoplot_base.handle import Handle, JointHandle
from .projection.projection import Projection
from .projection import Rotation

from warnings import warn

import numpy as np
from scipy.interpolate import interp2d
from scipy.spatial import cKDTree, ConvexHull

# Monkey patch pyfftw:
try:
	import pyfftw
	pyfftw.interfaces.cache.enable()
	def fftgrid(shape):
		return pyfftw.empty_aligned(shape, dtype='float64')
	fftjobs = lambda n_jobs : {"threads" : n_jobs}
	from pyfftw.interfaces.scipy_fftpack import fft2, fftfreq, ifft2, fftshift

except:
	fftgrid = np.zeros
	fftjobs = lambda n_jobs : dict()
	warn("Could not import pyFFTW: Monkey patching failed.")
	from scipy.fftpack import fft2, fftfreq, ifft2, fftshift

from math import ceil

from matplotlib.axes import Axes


# GEOPLOT:

class Geoplot(GeoplotBase):


	def __init__(self, ax, projection, limits_xy=None, gshhg_path=None,
	             which_ticks='significant', water_color='lightblue',
	             land_color='white', coast_color='black',verbose=0,
	             use_joblib=False, axes_margin_pt=5.0,
	             rotation=0, label_sign='label', tick_spacing=1.0,
	             secondary_tick_spacing_km=100,
	             secondary_tick_color='lightgray',
	             # Debugging:
	             _ax_background=None):
		"""
		Init method.

		Required arguments:
		   ax         :
		   projection :
		   limits_xy  : [xlim, ylim]

		Optional arguments:
		   rotation    : Select a rotation angle from 90, 180, or 270/-90.
		                 The canvas will be rotated accordingly.
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
		   tick_spacing :
		   secondary_tick_spacing_km : Number or None.
		   secondary_tick_color :
		"""

		if not isinstance(projection,Projection):
			raise RuntimeError("The projection has to be of class 'Projection'!")

		self._rotation = rotation = int(rotation)
		if rotation != 0:
			if not rotation in (90,180,270,-90):
				raise ValueError("Only rectangular rotations supported: 0, "
				                 "90, 180, 270/-90.")

			# Wrap projection by a rotation projection:
			projection = Rotation(projection, rotation)

		assert tick_spacing is None or \
		       (isinstance(tick_spacing,float) and tick_spacing > 0.0)

		secondary_tick_spacing_km = float(secondary_tick_spacing_km)


		super().__init__(ax, projection, gshhg_path, which_ticks, tick_spacing,
		                 secondary_tick_spacing_km, water_color, land_color,
		                 coast_color, secondary_tick_color, verbose, use_joblib,
		                 axes_margin_pt, label_sign, _ax_background)

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


	def set_ticks_off(self):
		"""
		Set ticks off.
		"""
		self._tick_spacing = None
		self._update_axes = True
		self._schedule_callback()


	def encompass(self, lon, lat):
		"""
		Set the x and y limits in the current projection so that the points given by
		lon and lat are contained.
		"""
		x,y = self._projection.project(lon, lat)
		self._user_xlim = (np.min(x), np.max(x))
		self._user_ylim = (np.min(y), np.max(y))
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

		Call signature:

		coastline(level, water_color=None, land_color=None,
		          coast_color=None, zorder=0)
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
			if self._scheduled[i].routine() == 'coastline':
				del self._scheduled[i]
				break

		# Schedule coastline:
		h = Handle('coastline', (level, zorder), kwargs)
		self._scheduled += [h]
		self._schedule_callback()

		return h


	def compass(self, lon=None, lat=None, x=None, y=None, color='black',
	            size=1):
		"""
		Plot a map compass at a specific location.

		Keyword parameters:
		   lon, lat : Longitude and latitude of North arrow origin.
		   x,y      : X- and y-coordinates of North arrow origin in
		              projected coordinates. Exclusive with lon/lat.
		   color    : Color of compass.
		              (Default: 'black')
		   size     : Size of compass cross in inches.
		              (Default: 1)
		"""
		if lon is None:
			if x is None or y is None or lat is not None:
				raise RuntimeError("Either lon,lat or x,y have to be given, "
				                   "exclusively.")
			lon, lat = self._projection.inverse(x,y)
		else:
			if not isinstance(lon,float) and not isinstance(lon,int) or \
			   not isinstance(lat,float) and not isinstance(lat,int):
				raise TypeError("'lon' and 'lat' have to be numbers when "
				                "given.")

		if not isinstance(size,float) and not isinstance(size,int):
			raise TypeError("'size' has to be a number.")

		h = Handle('compass', ((lon,lat), color, size), dict())
		self._scheduled += [h]
		self._schedule_callback()

		return h


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


	def scatter(self, *args, sort_by=None, **kwargs):
		"""
		Scatter plot. Has two call signatures:

		scatter(lon, lat):
		    Call scatter plot with two arrays of longitude
		    and latitude coordinates (in arcdegrees).

		scatter(point_set):
		    Call with a single argument, a unephy PointSet.

		Optional keyword arguments:
		    sort_by  : An array of parameters whereby to sort the plot
		               parameters.
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

		if sort_by is not None:
			order = np.argsort(sort_by)
			N = len(order)
			lon = lon[order]
			lat = lat[order]
			if "s" in kwargs:
				if isinstance(kwargs["s"], np.ndarray):
					kwargs["s"] = kwargs["s"][order]
			if "c" in kwargs:
				if isinstance(kwargs["c"], np.ndarray) and kwargs["c"].shape[0] == N:
					kwargs["c"] = kwargs["c"][order,...]

		self._add_data(lon=lon, lat=lat)
		h = Handle('scatter', (lon, lat), kwargs)
		self._scheduled += [h]
		self._schedule_callback()

		return h


	def plot(self, *args, **kwargs):
		"""
		Plot. Has two call signatures:

		plot(lon, lat):
		    Call plot with two arrays of longitude
		    and latitude coordinates (in arcdegrees).

		Plot(point_set):
		    Call with a single argument, a unephy PointSet.

		Optional keyword arguments:
		    **kwargs : Matplotlib plot kwargs.
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
		h = Handle('plot', (lon, lat), kwargs)
		self._scheduled += [h]
		self._schedule_callback()

		return h


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

		self._add_data(lon=lon, lat=lat)

		# Schedule quiver:
		h = Handle('quiver', (lon, lat, u, v, c), kwargs)
		self._scheduled += [h]
		self._schedule_callback()

		return h


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
		# Test for unephy:
		_has_unephy, ScalarField, CoordinateSystem,\
		MapProjectionSystem, GeographicSystem, SpatialDataSet,\
		SymmetricTensorField \
		   = has_unephy()

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
			if not _has_unephy:
				raise RuntimeError("Call signature with a spatial "
				                   "dataset requires unephy!")
			data = args[0]
			assert isinstance(data,SpatialDataSet)
			assert data.data_shape() == (1,)

			# Obtain raw data:
			with GeographicSystem():
				coords = data.coordinates().raw("arcdegree")
			lon = coords[...,0].reshape(-1)
			lat = coords[...,1].reshape(-1)
			azimuth = data.raw("arcdegree").reshape(-1)

			# Make sure that the tensor is actually a grid in the current
			# coordinate system:
			system = CoordinateSystem.current()
			plot_projection = self._projection._projection \
			                  if isinstance(self._projection, Rotation) else \
			                  self._projection
			if not isinstance(system,MapProjectionSystem) or \
			   system._projection._projection != plot_projection:
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
		handles = [self.quiver(lon, lat, np.sin(az_rad), np.cos(az_rad), c=c, **kw)]

		# Uncertainties:
		if error is not None:
			raise NotImplementedError("Error in orientations plot not yet implemented!")

		else:
			# Do scatter plot to mark origin of quivers:
			# Empirical equation for s:
			s = 0.8 * np.sqrt(kw["width"] / 0.25)
			handles += [self.scatter(lon, lat, edgecolor=kw["color"],
			                         c=kw["color"] if markerfill else "none",
			                         zorder=kw["zorder"], s=s,
			                         linewidth=kw["width"])]

			return JointHandle("orientations", handles)


	def streamplot_projected(self, x, y, u, v, data_mask=None,
	                         backend='matplotlib', show_border=False,
	                         **kwargs):
		"""
		Streamplot.

		Required arguments:
		   x, y : 1d-arrays defining the grid in projected coordinates
		   u    : 2d-grid of vector components in longitude direction.
		   v    : 2d-grid of vector components in latitude direction.

		Optional arguments:
		   data_mask : Masks parts of the vector field if set to False.
		               (Default: None)
		   backend   : The backend to use. Either 'matplotlib' or 'custom'.
		               (Default: 'matplotlib')
		   kwargs    : Passed to matplotlib streamplot or PolyCollection,
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
		h = Handle('streamplot', (x, y, u, v, backend, show_border,
		                          data_mask),
		           kwargs)
		self._scheduled += [h]
		self._schedule_callback()

		return h


	def streamplot(self, lon, lat, u, v, backend='matplotlib', **kwargs):
		# TODO Interpolate data to grid and do streamplot on grid!
		# TODO : Also convert start point!
		raise NotImplementedError("Geoplot.streamplot() not implemented yet.")


	def scalar_field(self, lon=None, lat=None, scalar=None, x=None, y=None,
	                 coastcolor='lightgray', watercolor="black", landcolor="none",
	                 cmap='default', coastmask=True, resample=False, n_resample=400,
	                 resample_method='nearest', show_border=True,
	                 data_mask=None, fadeout_distance=None,
	                 background='white', hillshade=None,
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
		   coastcolor       :
		                      (Default: 'lightgray')
		   watercolor       :
		                      (Default: 'black')
		   landcolor        :
		                      (Default: 'none')
		   cmap             :
		                      (Default: 'default')
		   coastmask        : Whether to clip the image at the land borders.
		                      (Default: True)
		   resample         : Whether to resample the field.
		                      (Default: False)
		   n_resample       :
		                      (Default: 400)
		   resample_method  : One of 'nearest' and 'spline'
		                      (Default: 'nearest')
		   show_border      :
		                      (Default: True)
		   colorbar         :
		                      (Default: None)
		   data_mask        :
		                      (Default: None)
		   fadeout_distance :
		                      (Default: None)
		   background       : Which background to fade out into for
		                      masked data.
		                      (Default: 'white')
		   hillshade        :
		                      (Default: None)
		"""
		# Test for unephy:
		_has_unephy, ScalarField, CoordinateSystem,\
		MapProjectionSystem, GeographicSystem, SpatialDataSet,\
		SymmetricTensorField \
		   = has_unephy()

		# Preprocess the hillshade:
		if _has_unephy and isinstance(hillshade,SpatialDataSet):
			hillshade = hillshade.raw('1').reshape(hillshade.shape())

		if _has_unephy and isinstance(scalar,ScalarField):
			# Sanity checks:
			if not (x is None and y is None and lon is None and lat is None):
				raise RuntimeError("If scalar is given, all of x, y, lon, and lat "
				                   " have to be None!")

			# Make sure that the scalar is actually a grid in the current
			# coordinate system:
			system = CoordinateSystem.current()
			plot_projection = self._projection._projection \
			                  if isinstance(self._projection, Rotation) else \
			                  self._projection
			if not isinstance(system,MapProjectionSystem) or \
			   system._projection._projection != plot_projection:
				raise RuntimeError("We need to be in a projection environment "
				                   "fitting to this Geoplot's projection!")
			with system:
				if not scalar.coordinates().is_grid():
					raise RuntimeError("Tensor is not a grid in plot coordinates!")

			# Now obtain coordinates and data from the scalar field:
			with GeographicSystem():
				coordinates = scalar.coordinates()
				lonlat = coordinates.raw("arcdegree")
				x,y = self._projection.project(lonlat[...,0], lonlat[...,1])
				scalar_ = scalar.raw(scalar.unit())

				# Assert the grid setup:
				xerr_xfirst = np.abs(x[:,0].reshape(-1,1) - x).max()
				xerr_yfirst = np.abs(x[0,:].reshape(1,-1) - x).max()
				yerr_xfirst = np.abs(y[0,:].reshape(1,-1) - y).max()
				yerr_yfirst = np.abs(y[:,0].reshape(-1,1) - y).max()

				if xerr_xfirst > 1e4*xerr_yfirst:
					if yerr_xfirst > 1e4*yerr_yfirst:
						# Have y first!
						grid_order = 'yx'
					else:
						raise RuntimeError("Could not identify grid order!")
				elif xerr_yfirst > 1e4*xerr_xfirst:
					if yerr_yfirst > 1e4*yerr_xfirst:
						# Have x first!
						grid_order = 'xy'
					else:
						raise RuntimeError("Could not identify grid order!")
				else:
					raise RuntimeError("Could not identify grid order!")

				if scalar.usability_mask() is not None:
					mask = scalar.usability_mask()
					# Determine resulting array shape:
					ids = np.argwhere(mask)
					i0 = ids[...,0].min()
					i1 = ids[...,0].max()
					j0 = ids[...,1].min()
					j1 = ids[...,1].max()
					shape = (i1-i0+1, j1-j0+1)

					x = x[mask].reshape(shape)
					y = y[mask].reshape(shape)
					scalar_ = scalar_[mask].reshape(shape)

					if hillshade is not None:
						hillshade = hillshade[mask].reshape(shape)
				else:
					mask = Ellipsis

				# Apply the data mask:
				if data_mask is None:
					data_mask = scalar.data_mask()[mask]
					if np.any(~data_mask):
						if mask is not Ellipsis:
							data_mask = data_mask.reshape(shape)
					else:
						data_mask = None

				# Now rotate the data if y first:
				if grid_order == 'yx':
					x = x.mean(axis=0)
					y = y.mean(axis=1)
					scalar_ = scalar_.T
					if hillshade is not None:
						hillshade = hillshade.T
				else:
					x = x.mean(axis=1)
					y = y.mean(axis=0)

				# Handle rotation in this Geoplot's projection:
				if self._rotation == -90:
					scalar_ = scalar_[:,::-1]
					if hillshade is not None:
						hillshade = hillshade[:,::-1]
				elif self._rotation == 180:
					pass
				else:
					raise NotImplementedError('Image may not be rotated '
					                          'correctly.')

		else:
			assert isinstance(scalar, np.ndarray)
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
				raise ValueError("resample_method must be one of 'nearest' or "
				                 "'spline'!")

			if data_mask is not None:
				raise NotImplementedError("Resampling of data_mask "
				         "not yet implemented.")

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
		handles = []
		if coastcolor is not None or watercolor is not None or landcolor is not None:
			handles += [self.coastline(self._coast_level, water_color=watercolor,
			                           land_color=landcolor, coast_color=coastcolor,
			                           zorder=zorder+1)]

		# Make sure that hillshade is properly oriented:
		if hillshade is not None:
			hillshade = hillshade[::-1, ::-1]


		# Call imshow:
		if coordinate_type == 'geographic':
			raise NotImplementedError("imshow not yet implemented!")
		else:
			handles += [self.imshow_projected(scalar_, [x.min(),x.max()],
			                                 [y.min(),y.max()], cmap=cmap,
			                                 origin='lower', zorder=zorder,
			                                 data_mask=data_mask,
			                                 hillshade=hillshade,
			                                 fadeout_distance=fadeout_distance,
			                                 coastmask=coastmask, **kwdict,
			                                 cax=cax,cbar_label=cbar_label)]

		return JointHandle("scalar_field", handles)



	def tensorfield_symmetric_2d(self, lon=None, lat=None, t1=None, t2=None, angle=None,
	                             x=None, y=None, linewidth=1.0, streamcolor='white',
	                             coastcolor='lightgray', watercolor="black",
	                             landcolor="none", cmap='default', coastmask=True,
	                             resample=False, n_resample=400,
	                             resample_method='nearest',
	                             tensor=None, colormode='max',
	                             direction='max', colorbar='horizontal',
	                             cax=None, colorscale='lin',
	                             cbar_label=None, thickness='difference',
	                             data_mask=None, fadeout_distance=None,
	                             background='white', hillshade=None,
	                             show_border=True, **kwargs):
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
		   colorscale       : One of 'lin', 'log', and 'sqrt'.
		                      (Default: 'lin')
		   data_mask        :
		                      (Default: None)
		   fadeout_distance :
		                      (Default: None)
		   background       : Which background to fade out into for
		                      masked data.
		                      (Default: 'white')
		   hillshade        :
		                      (Default: None)
		"""
		if not colormode in ['max','min','maxabs','sum','angle','second_moment']:
			raise ValueError("colormode must be one of 'max', 'min', 'maxabs', "
			                 "'sum', 'second_moment', or 'angle'.")
		if not direction in ['max','min','maxabs']:
			raise ValueError("direction must be one of 'max', 'min', or 'maxabs'.")
		if not thickness in ['difference','abs']:
			raise ValueError("thickness must be one of 'difference', or 'abs'.")
		if not colorscale in ('lin','sqrt','log'):
			raise ValueError("'colorscale' must be one of 'lin', 'sqrt', "
			                 "or 'log'.")

		# Import unephy:
		_has_unephy, ScalarField, CoordinateSystem,\
		MapProjectionSystem, GeographicSystem, SpatialDataSet,\
		SymmetricTensorField \
		   = has_unephy()


		if tensor is not None:
			# See if unephy is installed:
			if not _has_unephy:
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

			# Make sure that the tensor is actually a grid in the current
			# coordinate system:
			system = CoordinateSystem.current()
			plot_projection = self._projection._projection \
			                  if isinstance(self._projection, Rotation) else \
			                  self._projection
			if not isinstance(system,MapProjectionSystem) or \
			   system._projection._projection != plot_projection:
				raise RuntimeError("We need to be in a projection environment fitting to "
				                   "this Geoplot's projection!")
			with system:
				if not tensor.coordinates().is_grid():
					raise RuntimeError("Tensor is not a grid in plot coordinates!")

			# Now obtain coordinates and data from the tensor field:
			with GeographicSystem():
				coordinates = tensor.coordinates()
				lonlat = coordinates.raw("arcdegree")
#				xy = coordinates.raw(system.default_unit())
				t1 = tensor.principal_component("first").raw(tensor.unit())
				t2 = tensor.principal_component("second").raw(tensor.unit())
				angle = tensor.principal_azimuth().raw("arcdegree")

				# Project using this plot's projection:
				x,y = self._projection.project(lonlat[...,0], lonlat[...,1])

				# Assert the grid setup:
				xerr_xfirst = np.abs(x[:,0].reshape(-1,1) - x).max()
				xerr_yfirst = np.abs(x[0,:].reshape(1,-1) - x).max()
				yerr_xfirst = np.abs(y[0,:].reshape(1,-1) - y).max()
				yerr_yfirst = np.abs(y[:,0].reshape(-1,1) - y).max()


				if xerr_xfirst > 1e4*xerr_yfirst:
					if yerr_xfirst > 1e4*yerr_yfirst:
						# Have y first!
						grid_order = 'yx'
					else:
						raise RuntimeError("Could not identify grid order!")
				elif xerr_yfirst > 1e4*xerr_xfirst:
					if yerr_yfirst > 1e4*yerr_xfirst:
						# Have x first!
						grid_order = 'xy'
					else:
						raise RuntimeError("Could not identify grid order!")
				else:
					raise RuntimeError("Could not identify grid order!")

				if tensor.usability_mask() is not None:
					mask = tensor.usability_mask()
					# Determine resulting array shape:
					ids = np.argwhere(mask)
					i0 = ids[...,0].min()
					i1 = ids[...,0].max()
					j0 = ids[...,1].min()
					j1 = ids[...,1].max()
					shape = (i1-i0+1, j1-j0+1)

					x = x[mask].reshape(shape)
					y = y[mask].reshape(shape)
					t1 = t1[mask].reshape(shape)
					t2 = t2[mask].reshape(shape)
					angle = angle[mask].reshape(shape)

					if hillshade is not None:
						hillshade = hillshade[mask].reshape(shape)
				else:
					mask = Ellipsis

				# Apply the data mask:
				if data_mask is None:
					data_mask = tensor.data_mask()[mask]
					if np.any(~data_mask):
						if mask is not Ellipsis:
							data_mask = data_mask.reshape(shape)
					else:
						data_mask = None


				# Now Rotate the data if y first:
				if grid_order == 'yx':
					x = x.mean(axis=0)
					y = y.mean(axis=1)
					t1 = t1.T
					t2 = t2.T
					angle = angle.T
				else:
					x = x.mean(axis=1)
					y = y.mean(axis=0)


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
			xg, yg = np.meshgrid(x, y, indexing='ij')

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

		if colorscale == 'log':
			color = np.log10(np.abs(color))
		elif colorscale == 'sqrt':
			color = np.sqrt(np.abs(color))

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

		# Handle the width of the streamlines.
		# Set them to zero in masked areas.
		lw_streamplot = linewidth * (width[::-1,::-1] - width.min()) \
		                /(width.max() - width.min())

		# Save keys in addition to old ones:
		kwdict = dict(kwargs)
		kwdict["linewidth"] = lw_streamplot
		kwdict["facecolor"] = streamcolor
		kwdict["quiver"] = False

		# Obtain zorder:
		zorder = kwdict.pop("zorder",1)

		# If coastline is set, we show the coastline in the same color as
		# the streamplot:
		handles = []
		if coastcolor is not None or watercolor is not None or landcolor is not None:
			handles += [self.coastline(self._coast_level, water_color=watercolor,
			                           land_color=landcolor, coast_color=coastcolor,
			                           zorder=zorder+1)]

		# Call imshow:
		if coordinate_type == 'geographic':
			raise NotImplementedError("imshow not yet implemented!")
		else:
			handles += [self.imshow_projected(color, [x.min(),x.max()],
			                                  [y.min(),y.max()], cmap=cmap,
			                                  origin='lower', zorder=zorder,
			                                  coastmask=coastmask,
			                                  data_mask=data_mask,
			                                  fadeout_distance=fadeout_distance,
			                                  cax=cax,cbar_label=cbar_label,
			                                  hillshade=hillshade,
			                                  **imshow_kwargs)]

		# Call streamplot:
		handles += [self.streamplot_projected(xvals, yvals, u, v, backend='custom',
		                                      show_border=show_border,
		                                      data_mask=data_mask[::-1,::-1],
		                                      zorder=zorder+2, **kwdict)]

		return JointHandle("tensorfield_symmetric_2d",handles)



	def convex_hull(self, point_set=None, x=None, y=None, lon=None, lat=None,
	                usability_mask=True, **kwargs):
		"""
		Plot the convex hull of a point set.

		Returns a handle.
		"""
		if point_set is not None:
			# See if unephy is installed:
			try:
				from unephy import PointSet, CoordinateSystem,\
				                   MapProjectionSystem, GeographicSystem,\
				                   Field
			except ImportError:
				raise RuntimeError("If 'point_set' is set, it has to be an unephy "
				                   "PointSet instance. Could not "
				                   "import unephy!")

			# Sanity checks:
			if not (x is None and y is None and lon is None and lat is None):
				raise RuntimeError("If point_set is given, all of x, y, lon, and "
				                   "lat have to be None!")

			if not isinstance(point_set, PointSet):
				raise RuntimeError("'point_set' has to be a PointSet "
				                   "instance!")

			# Make sure that the tensor is actually a grid in the current
			# coordinate system:
			system = CoordinateSystem.current()
			plot_projection = self._projection._projection \
			                  if isinstance(self._projection, Rotation) else \
			                  self._projection
			if not isinstance(system,MapProjectionSystem) or \
			   system._projection._projection != plot_projection:
				raise RuntimeError("We need to be in a projection environment "
				                   "fitting to this Geoplot's projection!")

			# Now obtain coordinates and data from the tensor field:
			with GeographicSystem():
				coordinates = point_set.coordinates()
				lonlat = coordinates.raw("arcdegree")

			# For fields, handle usability mask:
			if usability_mask and isinstance(point_set, Field):
				if point_set.usability_mask() is not None:
					mask = point_set.usability_mask()
					lonlat = lonlat[mask]

			lonlat = lonlat.reshape((-1,2))

			# Project coordinates:
			x,y = self._projection.project(lonlat[:,0], lonlat[:,1])

			# Compute convex hull:
			xy = np.stack((x,y), axis=-1)
			hull = ConvexHull(xy)
			x = xy[hull.vertices, 0]
			y = xy[hull.vertices, 1]

		return self.polygon(x=x, y=y, lon=lon, lat=lat, **kwargs)


	def distortion(self, xlim=None, ylim=None,
	               cmap='inferno', cax=None,
	               contours='percent', labels=None,
	               min_samples=100, ):
		"""
		Plot the current projection's distortion.

		Keyword parameters:
		   xlim, ylim  : Specify custom rectangular limits for plotting the
		                 distortion. (Default: None)
		   contours    : Number of equidistortion contours to plot.
		                 (Default: 5)
		   cmap        : If given, plot distortion as a color image.
		                 (Default: None)
		   cax         : Axis where to plot the color bar.
		                 (Default: None)
		   min_samples : Minimum number of distortion samples on any axis.
		                 (Default: 100)
		"""
		# First determine the projection area:
		if xlim is None:
			xlim = self.xlim()
		else:
			if len(xlim) != 2:
				raise RuntimeError("'xlim' parameter has to consist of two floats!")
			xlim = (float(xlim[0]),float(xlim[1]))
		if ylim is None:
			ylim = self.ylim()
		else:
			if len(ylim) != 2:
				raise RuntimeError("'ylim' parameter has to consist of two floats!")
			ylim = (float(ylim[0]),float(ylim[1]))

		# Then evaluate the distortion on a grid:
		dx = xlim[1] - xlim[0]
		dy = ylim[1] - ylim[0]
		if dx > dy:
			gc = dy / min_samples
		else:
			gc = dx / min_samples
		nx = max(round(dx/gc),min_samples)
		ny = max(round(dy/gc),min_samples)
		gx, gy = np.meshgrid(np.linspace(xlim[0],xlim[1],nx),
		                     np.linspace(ylim[0],ylim[1],ny),
		                     indexing='ij')
		glon,glat = self._projection.inverse(gx,gy)
		k = self._projection.scale_k(glon, glat) - 1.0

		# Now continue with these distortion data.
		# First color image representation:
		handles = []
		if cmap is not None:
			handles += [self.imshow_projected(k[:,::-1], xlim, ylim, cmap=cmap)]

		# Contours:
		if contours in ('percent','%'):
			scale=50.
			contours= []
			while len(contours) < 5:
				scale *= 2
				cmin = np.floor(k.min() * scale)
				cmax = np.ceil(k.max() * scale)
				if cmin == cmax:
					contours = np.array([cmin])
				else:
					contours = np.arange(cmin,cmax)
				contours /= scale

		# Then contour:
		if isinstance(contours,float) or isinstance(contours,int):
			if contours > 0:
				handles += [self.contour(k, contours, x=gx, y=gy, labels=labels,
				                         colors='k', fmt='%1.2f')]
			else:
				pass
		elif isinstance(contours,list) or isinstance(contours,np.ndarray):
			if scale == 100.0:
				fmt='%1.2f'
			else:
				fmt='%1.2e'
			handles += [self.contour(k, contours, x=gx, y=gy, labels=labels,
			                         colors='k', fmt=fmt)]
		elif contours is None:
			pass
		else:
			raise NotImplementedError()

		return JointHandle("distortion",handles)



	def imshow_projected(self, z, xlim, ylim, cax=None,
	                     cbar_label=None, data_mask=None,
			             fadeout_distance=None, background_color='white',
			             hillshade=None, hillshade_strength=0.3,
			             hsmin=0.0, hsmax=1.0,
			             n_jobs=10, **kwargs):
		"""
		Plot a field (in projected coordinates) using imshow.
		"""
		# Convert to numpy:
		if not isinstance(xlim,np.ndarray):
			xlim = np.array(xlim)
		if not isinstance(ylim,np.ndarray):
			ylim = np.array(ylim)

		if data_mask is not None:
			if not isinstance(data_mask,np.ndarray) \
			   or not data_mask.dtype == bool:
				raise RuntimeError("'data_mask' has to be a boolean numpy "
				                   "array!")
			if data_mask.shape != z.shape:
				print("data_mask.shape:",data_mask.shape)
				print("z.shape:        ",z.shape)
				raise RuntimeError("'data_mask' has to be of the same shape as "
				                   "the z-data!")
			if data_mask is not None and fadeout_distance is None:
				raise RuntimeError("If 'data_mask' is (implicitly) given, "
				                   "'fadeout_distance' has to be given as well!")

		# Convert fadeout distance:
		try:
			from unephy import CoordinateSystem, Datum
			_has_unephy = True
			if isinstance(fadeout_distance, Datum):
				unit = CoordinateSystem.current().default_unit()
				fadeout_distance = fadeout_distance.raw(unit)
		except ImportError:
			_has_unephy = False

		# Handle data mask:
		if data_mask is not None:
			# First create a margin. Create a margined-mask
			# that is set to True whenever a grid cell of
			# the data_mask is set to True and within
			# 0.5*fadeout_distance to the cell:
			x = np.linspace(xlim[0],xlim[1],z.shape[0])
			y = np.linspace(ylim[0],ylim[1],z.shape[1])
			xy = np.stack(np.meshgrid(x, y, indexing='ij'),axis=2)
			tree = cKDTree(xy[data_mask])
			margined_mask = data_mask.copy()
			query = tree.query_ball_point(xy[~data_mask],
			                              fadeout_distance)
			margined_mask[~data_mask] = [len(q) > 0 for q in query]

			# Grid constant:
			height = ylim[1]-ylim[0]
			width = xlim[1] - xlim[0]
			shape = z.shape
			dx = width / (shape[0]-1.)
			dy = height / (shape[1]-1.)

			# Obtain padded grid:
			L_pad = 1.5 * fadeout_distance
			Nx_pad = ceil(L_pad / width * shape[0])
			Ny_pad = ceil(L_pad / height * shape[1])
			padded_shape = (2*Nx_pad+shape[0], 2*Ny_pad+shape[1])
			padded_grid = fftgrid(padded_shape)
			padded_grid[...] = 0.0
			grid_index = (slice(Nx_pad,Nx_pad+shape[0]),
			              slice(Ny_pad,Ny_pad+shape[1]))

			# Fourier transform the mask:
			fx = fftfreq(padded_shape[0], d=dx)
			fy = fftfreq(padded_shape[1], d=dy)
			fx,fy = np.meshgrid(fx, fy, indexing='ij')

			# Transfer to padded grid:
			padded_grid[grid_index] = margined_mask

			# Fourier transform:
			fft = fft2(padded_grid, **fftjobs(n_jobs))

			# Window filter:
			def hann_fourier(L,f):
				window = (np.sinc(2*L*f) + 0.5*np.sinc(2*L*f - 1)
				            + 0.5*np.sinc(2*L*f + 1))
				return window

			window = hann_fourier(fadeout_distance, fx) \
			       * hann_fourier(fadeout_distance, fy)
			fft *= window

			# Retransform:
			ifft = ifft2(fft, **fftjobs(n_jobs)).reshape(padded_shape).astype(float)
			transition = ifft[grid_index]

#			transition = margined_mask

			# Some cosmetics:
			transition[data_mask] = 1.0
			transition[transition < 0] = 0.0
			transition[transition > 1] = 1.0

		else:
			transition = data_mask

		# Hillshade:
		if hillshade is not None:
			if not isinstance(hillshade, np.ndarray):
				raise TypeError("'hillshade' has to be a numpy array!")
			if hillshade.shape != z.shape:
				print("hillshade.shape:",hillshade.shape)
				print("z.shape:        ",z.shape)
				raise ValueError("'hillshade' has to be of same shape as "
				                 "image data!")
			hillshade = hillshade - hillshade.min()
			hillshade /= hillshade.max()
			hillshade = hsmin + (hsmax-hsmin) * hillshade
			hillshade = hillshade.T

		# Check data limits:
		self._add_data(x=xlim, y=ylim)

		# Schedule plot:
		h = Handle('imshow', (z, xlim, ylim, cbar_label, transition,
		                      background_color, hillshade,
		                      hillshade_strength),
		           kwargs)
		self._scheduled += [h]
		self._schedule_callback()

		# Colorbar:
		if cbar_label is not None or cax is not None:
			self.colorbar(h, label=cbar_label, cax=cax)

		return h


	def contour(self, z, levels=None, x=None, y=None, lon=None, lat=None,
	            labels=None, colors=None, **kwargs):
		"""
		Plot contours of a scalar field.

		Returns a handle.
		"""

		# Unephy:
		try:
			from unephy import SpatialDataSet, GeographicSystem
			if isinstance(z, SpatialDataSet):
				with GeographicSystem():
					ll = z.coordinates().raw('arcdegree')
					lon = ll[...,0]
					lat = ll[...,1]
				z = z.raw(z.unit())
		except:
			pass

		# Check coordinate consistency:
		x,y = self._process_coordinates('xy', lon, lat, x, y)
		self._add_data(x=x, y=y)

		# Check data consistency:
		if not isinstance(z, np.ndarray):
			z = np.array(z)
		if x.ndim == 1:
			if y.ndim != 1:
				raise RuntimeError("Both coordinate arrays need to be of same dimension!")
			if z.ndim != 2:
				raise RuntimeError("Value array has to be two-dimensional!")
			if not x.size == z.shape[0] or not y.size == z.shape[1]:
				raise RuntimeError("Value array must be size (M,N) where M is "
				                   "x coordinate size and N y coordinate size.")
		elif x.ndim == 2:
			if y.ndim != 2:
				raise RuntimeError("Both coordinates need to have same "
				                   "two-dimensional shape!")
			if not np.array_equal(z.shape, x.shape) or \
			   not np.array_equal(x.shape, y.shape):
				raise RuntimeError("Two-dimensional array shapes need to be equal!")
		else:
			raise RuntimeError("Coordinate shape not compatible!")

		# Check labels:
		if labels == True or isinstance(labels,list):
			if labels == True:
				labels = [str(l) for l in levels]
			elif not all(isinstance(l,str) for l in labels):
				raise RuntimeError("Labels have to be strings!")
		else:
			labels = None

		if colors is not None:
			kwargs["colors"] = colors

		# Schedule contour:
		h = Handle('contour', (x, y, z, levels, labels), kwargs)
		self._scheduled += [h]
		self._schedule_callback()

		return h


	def _process_coordinates(self, destination, lon=None, lat=None, x=None, y=None):
		"""
		Check coordinate consistency.

		Parameters:
		   destination : One of 'keep', 'lonlat', 'xy'
		"""
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

		# Project if required:
		if destination == 'xy':
			if x is not None:
				return x, y
			else:
				return self._projection.project(lon, lat)
		elif destination == 'lonlat':
			if lon is not None:
				return lon, lat
			else:
				return self._projection.inverse(x,y)
		elif destination == 'keep':
			return lon, lat, x, y
		else:
			raise RuntimeError()


	def polygon(self, x=None, y=None, lon=None, lat=None, **kwargs):
		"""
		Plot a polygon.
		"""
		# Sanity checks:
		lon, lat, x, y = self._process_coordinates('keep', lon, lat, x, y)

		# Check data limits:
		self._add_data(x=x, y=y, lon=lon, lat=lat)

		# Schedule plot:
		h = Handle('polygon', (x, y, lon, lat), kwargs)
		self._scheduled += [h]
		self._schedule_callback()

		return h


	def is_land(self, lon=None, lat=None, x=None, y=None):
		"""
		Query whether a series of coordinates are within the
		land mass.
		"""
		# Sanity checks:
		x,y = self._process_coordinates('xy', lon, lat, x, y)

		x,y = self._plot_canvas.obtain_coordinates(x, y, self._xlim, self._ylim)

		xy = np.stack((x, y), axis=-1)

		# Now check if coast path contains these coordinates:
		mask = is_land_backend(xy, self._projection.identifier(),
		                       self._xlim, self._ylim,
		                       coastpath=self._coast_path_prep)

		return mask



	def colorbar(self, handle=None, label=None, cax=None,
	             location='upper center',
	             orientation='horizontal', **kwargs):
		"""
		Plot a color bar.

		Call signature:

		colorbar(handle=None, label=None, cax=None)
		   handle : Handle with color map information.
		   label  : Axis label.
		   cax    : Matplotlib axis to plot the colorbar on.

		Optional keyword arguments:
		   orientation : One of 'horizontal' and 'vertical'
		                 (Default: 'horizontal')
		   location    : One of 'upper center', ... . Not used
		                 when cax is given.
		"""
		# Sanity checks:
		if handle is not None:
			if not isinstance(handle, Handle):
				raise TypeError("'handle' has to be a Handle instance!")
			handle = handle.cbar_handle()

		if label is not None and not isinstance(label,str):
			try:
				label = str(label)
			except:
				raise RuntimeError("Could not convert 'label' to string!")

		if cax is not None and not isinstance(cax,Axes):
			raise TypeError("'cax' has to be a matplotlib.axes.Axes "
			                "instance!")

		if orientation not in ('horizontal','vertical'):
			raise RuntimeError("Orientation must be one of 'horizontal' "
			                   "or 'vertical'!")

		if location not in ('upper center',):
			raise ValueError("Location must be one of 'upper center'.")

		# Obtain an axis for color bar if none given:
		if cax is None:
			pos = self._ax.get_position()
			cax = self._ax.figure.add_axes(pos)
			self._ax.set_position(pos)

			if location == 'upper center':
				l,b,w,h = pos.bounds
				# Hotfix:
				cax.set_position((l+0.33*w, b+0.9*h, 0.34*w, 0.03*h))

		# Schedule plot:
		h = Handle('colorbar', (handle, label, cax, orientation), kwargs)
		self._scheduled += [h]
		self._schedule_callback()

		return h
