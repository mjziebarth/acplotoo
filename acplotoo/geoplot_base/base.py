# -*- coding: utf-8 -*-
#
# Acplotoo GeoplotBase class file.
# Copyright (C) 2018 Malte Ziebarth
# 
# This software is distributed under the MIT license.
# See the LICENSE file in this repository.

from .backend import _generate_axes_boxes, _create_tick_arrays, has_joblib,\
                     identify_jumps, read_coastlines, coast_line_patches_and_path,\
                     _generate_axes_ticks
from .rect import Rect
from .streamplot import _streamplot_calculate_polygons


from matplotlib.collections import PatchCollection, LineCollection, PolyCollection
from matplotlib.patches import Polygon, Rectangle, PathPatch
from matplotlib.path import Path
from matplotlib import rcParams
from shapely.geometry.polygon import Polygon as SPoly
from shapely.geometry import LineString, MultiPolygon, box, MultiLineString, MultiPolygon
from shapely.ops import polygonize_full, unary_union, split
from shapely.prepared import prep
from shapely import speedups
import numpy as np
from scipy.interpolate import RectBivariateSpline

from warnings import warn

from datetime import datetime

# TODO fix!
from matplotlib.transforms import Bbox

# Enable speedups:
if speedups.available:
	speedups.enable()


# Utility:
def _adjust_axes_position_callback(ax0, ax1):
	pass
	if ax0.get_position() != ax1.get_position():
		ax1.set_position(ax0.get_position())


class GeoplotBase:
	def __init__(self, ax, projection, gshhg_path, which_ticks,
	             water_color, land_color, coast_color, verbose, use_joblib,
	             resize_figure):

		if resize_figure:
			self._fig = ax.get_figure()
			if len(self._fig.get_axes()) != 1 or ax is not self._fig.get_axes()[0]:
				raise RuntimeError("Cannot automatically resize figure if figure "
				                   "has more than one axis!")
			self._figsize_original = self._fig.get_size_inches()

		# TODO : Checks!
		self._ax = ax
		pos = self._ax.get_position()
		self._ax_pos = [pos.x0, pos.y0, pos.width, pos.height]

		# Obtain the original axes size in inches:
		size = self._ax.get_window_extent().transformed(self._ax.get_figure()
		                                                .dpi_scale_trans.inverted())
		self._ax_dx = size.x1 - size.x0
		self._ax_dy = size.y1 - size.y0

		# Because of set_aspect restriction, we cannot use twinx/twiny. We thus
		# have to emulate that behaviour of a twin axis. Only the labels of the
		# ax2 are relevant, so we can set it beneath the other axis:
		self._ax2 = ax.get_figure().add_axes(self._ax_pos, zorder=self._ax.zorder-1)
		self._ax.add_callback(lambda : _adjust_axes_position_callback(self._ax,self._ax2))
		self._projection = projection

		self._which_ticks = which_ticks
		self._water_color = water_color
		self._land_color = land_color
		self._coast_color = coast_color

		if use_joblib and not has_joblib:
			raise ImportError("Could not import joblib but was "
			                  "requested to do so!")
		self._use_joblib = use_joblib


		# Initialize axes:
		self._ax2.spines['top'].set_visible(False)
		self._ax2.spines['bottom'].set_visible(False)
		self._ax2.spines['left'].set_visible(False)
		self._ax2.spines['right'].set_visible(False)
		self._ax2.xaxis.tick_top()
		self._ax2.yaxis.tick_right()
		self._ax2.tick_params(top=False, bottom=False, left=False, right=False)

		# See if we can load GSHHG:
		if gshhg_path is None:
			self._gshhg_path = None
		else:
			self._gshhg_path = gshhg_path

		# Setup internal data:
		self._resize_figure = resize_figure
		self._use_latex = rcParams["text.usetex"]
		self._data_xlim = None
		self._data_ylim = None
		self._xlim = None
		self._ylim = None
		self._user_xlim = None
		self._user_ylim = None
		self._aspect = 1.0
		self._ticks = None
		self._scheduled = []
		self._coasts = None
		self._coast_level=4
		self._coast_kwargs_base = {'linewidth' : 0.5}
		self._clip_rect = Rectangle([0.0, 0.0], 1.0, 1.0)
		self._grid_on = True
		self._grid_constant = 1.0
		self._grid_kwargs_base = {'color' : 'gray', 'linewidth' : 0.5}
		self._grid_kwargs = {**self._grid_kwargs_base}
		self._grid_anchor = (0.0,0.0)
		self._grid_ticks_between = 8
		self._verbose = verbose
		self._adjusted = False
		self._grid_handles = []
		self._tick_dict = None
		self._update_axes = False
		self._update_grid = False
		self._box_axes_linewidth = 0.5
		self._zorder_axes_0 = 19
		self._streamplot_config = {'interpolation_points_per_axis' : 1000,
		                           'start_points_per_axis' : 10,
		                           'forward' : True, 'backward' : True,
		                           'minlength' : 0.01, 'maxlength' : 1.0,
		                           'step_len_min' : 1e-4,
		                           'collision_radius' : 1e-2,
		                           'arrow_step_len' : 0.05, 'arrow_len' : 0.025,
		                           'max_steps' : 2e4, 'tolerance' : 1e-3,
		                           'linewidth_base' : 1./72.,
		                           'kwargs' : {'edgecolor' : 'none'},
		                           'quiver' : True}

	def _has_initialized_axes(self):
		"""
		Returns true if we can plot. This basically means that we
		have initialized the projection space.
		That initialization is inherent to global projections.
		"""
		return self._projection.is_global() \
		    or (self._data_xlim is not None and self._data_ylim is not None) \
		    or (self._user_xlim is not None and self._user_ylim is not None)

	def _read_gshhg(self):
		if self._gshhg_path is None:
			raise RuntimeError("GSHHG not loaded!")

		# Call the backend which provides a joblib caching:
		self._coasts = read_coastlines(self._gshhg_path, self._projection,
		                   self._projection.identifier(),
		                   np.array(self._geographic_extents),
		                   np.array(self._xlim), np.array(self._ylim),
		                   self._verbose)

	def _imshow(self, z, xlim, ylim, **kwargs):
		# Actually we would have to add half the gridsize for
		# pixel-registered grid, but for small grids it's
		# okay for now.
		coastmask = kwargs.pop("coastmask",None)
		xlim, ylim = self._plot_canvas.obtain_coordinates(xlim, ylim,
		                      self._xlim, self._ylim)
		h = self._ax.imshow(z, extent=(xlim[0],xlim[1],ylim[0],ylim[1]),
		                    **kwargs)
		clip_path = self._clip_rect
		if coastmask is not None and coastmask and self._coast_patch is not None:
			h.set_clip_path(self._coast_patch())
		else:
			h.set_clip_path(clip_path)

		return h

	def _quiver(self, lon, lat, u, v, c, **kwargs):
		# Quiver plot!

		# Checks:
		if isinstance(lon,list):
			lon = np.array(lon)
		if isinstance(lat,list):
			lat = np.array(lat)

		# Convert coordinates:
		x,y = self._projection.project(lon,lat)

		# Obtain plot coordinates:
		x,y = self._plot_canvas.obtain_coordinates(x, y, self._xlim, self._ylim)

		# Obtain unit vectors at locations:
		east = self._projection.unit_vector_east(lon=lon, lat=lat)
		north = self._projection.unit_vector_north(lon=lon, lat=lat)

		# Obtain vector components in current projection:
		vx = east[:,0] * u + north[:,0] * v
		vy = east[:,1] * u + north[:,1] * v

		# Quiver:
		if c is None:
			h = self._ax.quiver(x, y, vx, vy, **kwargs)
		else:
			h = self._ax.quiver(x, y, vx, vy, c, **kwargs)
		h.set_clip_path(self._clip_rect)

		return h


	def _streamplot(self, x, y, u, v, backend, show_border, **kwargs):
		# Streamplot!

		# Checks:
		if isinstance(x,list):
			x = np.array(x)
		if isinstance(y,list):
			y = np.array(y)

		# Obtain plot coordinates:
		x,y = self._plot_canvas.obtain_coordinates(x, y, self._xlim, self._ylim)

		# Convert starting points:
		if "start_points" in kwargs:
			sp = kwargs["start_points"]
			sp = self._plot_canvas.obtain_coordinates(sp[:,0], sp[:,1], self._xlim,
			                                          self._ylim)
			kwargs["start_points"] = np.concatenate([sp[0][:,np.newaxis],
			                                         sp[1][:,np.newaxis]],axis=1)

		# Obtain unit vectors at locations:
		east = self._projection.unit_vector_east(x=x, y=y)
		north = self._projection.unit_vector_north(x=x, y=y)

		# Obtain vector components in current projection:
		vx = east[:,0] * u + north[:,0] * v
		vy = east[:,1] * u + north[:,1] * v

		# Streamplot:
		if backend == 'matplotlib':
			h = self._ax.streamplot(x, y, vx, vy, **kwargs)
			h.lines.set_clip_path(self._clip_rect)
			h.lines.set_capstyle('round')
			h.arrows.set_clip_path(self._clip_rect)
		elif backend == 'custom':
			# Local copy of config, override with kwargs:
			conf = dict(self._streamplot_config)
			kwargs_poly = conf["kwargs"]
			for k in conf.keys():
				if k in kwargs.keys():
					conf[k] = kwargs[k]
					kwargs.pop(k,None)
			for k in kwargs.keys():
				kwargs_poly[k] = kwargs[k]

			# Interpolate:
			nr = conf['interpolation_points_per_axis']
			xr,yr = np.linspace(x.min(), x.max(), nr), np.linspace(y.min(), y.max(), nr)
			bispline_x = RectBivariateSpline(x,y,vx)
			bispline_y = RectBivariateSpline(x,y,vy)
			vxg = bispline_x(xr, yr)
			vyg = bispline_y(xr, yr)
			xg,yg = np.meshgrid(xr,yr,indexing='ij')

			# Obtain linewidth:
			if "linewidth" in kwargs:
				linewidth = kwargs.pop("linewidth",1.0)
				if isinstance(linewidth,np.ndarray):
					bispline_lw = RectBivariateSpline(x,y,linewidth)
					lw = bispline_lw(xr,yr)
					linewidth = conf["linewidth_base"] * linewidth.max() * lw / lw.max()
			else:
				linewidth = np.ones(xg.shape)

			# Remove some unused kwargs that do or may exist:
			unused_keywords = ['cmap','linewidth']
			for kw in unused_keywords:
				if kw in kwargs_poly.keys():
					warn("The keyword '" + kw + "' does not have any effect in "
						 "streamplot.")
			for kw in ["linestyle","density","start_points"] + unused_keywords:
				kwargs_poly.pop(kw,None)

			# Determine kwargs of arrowheads:
			kwargs_quiver = {"pivot" : "mid", "width" : .07, "minshaft" : 1.0,
			                 "minlength" : 0, "units" : "xy", "headaxislength" : 5}
			if "zorder" in kwargs_poly:
				kwargs_quiver["zorder"] = kwargs_poly["zorder"]
			if "facecolor" in kwargs_poly:
				kwargs_quiver["color"] = kwargs_poly["facecolor"]

			# Obtain start points:
			if "start_points" in kwargs:
				start_xy = kwargs.pop("start_points",None)
				start_x = start_xy[:,0]
				start_y = start_xy[:,1]
			else:
				ns = conf["start_points_per_axis"]
				start_x, start_y = np.meshgrid(np.linspace(x.min(), x.max(), ns),
				                               np.linspace(y.min(), y.max(), ns))

			# Scale the relative parameters of trajectory lengths to the
			# current data length scale:
			xlim,ylim = self._plot_canvas.obtain_coordinates(np.array(self._xlim), 
			                                                 np.array(self._ylim),
			                                                 self._xlim, self._ylim)
			lenscale = np.sqrt((xlim[1]-xlim[0])**2
			                   +(ylim[1]-ylim[0])**2)
			maxlength = conf["maxlength"] * lenscale
			minlength = conf["minlength"] * lenscale
			step_len_min = conf["step_len_min"] * lenscale
			arrow_step_len = conf["arrow_step_len"] * lenscale
			collision_radius = conf["collision_radius"] * lenscale

			# Calculate polygons:
			polygons, arrows = _streamplot_calculate_polygons(xg, yg, vxg, vyg,
			                        linewidth, minlength, maxlength,
			                        step_len_min, arrow_step_len,
			                        collision_radius,
			                        start_x.flatten(), start_y.flatten(),
			                        conf["forward"], conf["backward"],
			                        conf["max_steps"], conf["tolerance"])

			# Scale arrows:
			arrows[:,2:] *= conf["arrow_len"] * lenscale

			# Add poly collection:
			h = []

			# If we are to draw the polygon boundaries, draw polygons twice:
			# Once with edgecolor and once, overlaid, without. This prevents
			# the background color drowning the foreground color.
			if show_border:
				kwargs_poly2 = dict(**kwargs_poly)
				kwargs_poly2["edgecolor"] = 'black'
				h += [self._ax.add_collection(PolyCollection(polygons,
				                                  **kwargs_poly2))]
				h[-1].set_clip_path(self._clip_rect)

			h += [self._ax.add_collection(PolyCollection(polygons, **kwargs_poly))]
			h[-1].set_clip_path(self._clip_rect)
			if conf["quiver"]:
				h += [self._ax.quiver(arrows[:,0],arrows[:,1],arrows[:,2],arrows[:,3],
				                      **kwargs_quiver)]
				h[1].set_clip_path(self._clip_rect)
		return h


	def _coastline(self, level, zorder, **kwargs):
		# Create keyword dictionary:
		kws = dict(self._coast_kwargs_base)
		kws.update(kwargs)

		# See whether we have to read GSHHS:
		if self._coasts is None \
		  or self._xlim[0] < self._coast_xlim[0] \
		  or self._xlim[1] > self._coast_xlim[1] \
		  or self._ylim[0] < self._coast_ylim[0] \
		  or self._ylim[1] > self._coast_ylim[1]:
			# Read GSHHS:
			self._read_gshhg()
			self._coast_xlim = self._xlim
			self._coast_ylim = self._ylim

		if self._verbose > 0:
			print("_coastline ...")
		t0 = datetime.now()

		# Of loaded polygons, plot all polygons with level <= level:
		coords = []
		for poly in self._coasts:
			if poly[0] <= level:
				# Obtain coordinates in canvas coordinates:
				x,y = self._plot_canvas.obtain_coordinates(poly[1],poly[2], self._xlim,
				                                           self._ylim)
				coords += [(x,y,poly[0])]

		if self._verbose > 1:
			print("self._coasts.size:",len(self._coasts))

		cnvs_x = np.array([self._plot_canvas.x0, self._plot_canvas.x1])
		cnvs_y = np.array([self._plot_canvas.y0, self._plot_canvas.y1])

		# Call backend to generate patches and path clipped to canvas.
		# The backend employs joblib if wished and possible to load for
		# repeated calls.
		patches_xy, coast_path =\
		    coast_line_patches_and_path(coords, cnvs_x, cnvs_y, self._xclip,
		                                self._yclip)

		# Add also base rectangle showing water:
		water_rect = Rectangle([self._clip_rect.get_x(),self._clip_rect.get_y()],
		                                self._clip_rect.get_width(),
		                                self._clip_rect.get_height(),
		                                edgecolor='none',
			                            facecolor=self._water_color,
			                            zorder=-1)
		self._water_artist = self._ax.add_patch(water_rect)

		# Prepare coast path geometry and create patches:
		self._coast_path = coast_path
		self._coast_path_prep = prep(coast_path)
		patches = [Polygon(xy) for xy in patches_xy]

		# Plot all polygons:
		if len(patches) > 0:
			colors = [[self._water_color, self._land_color][c[2] % 2] for c in coords]
			h = self._ax.add_collection(PatchCollection(patches,facecolors=colors,
				                                        edgecolors=self._coast_color,
				                                        zorder=zorder, **kws))
			h.set_clip_path(self._clip_rect)
		else:
			h = None

		if self._verbose > 0:
			print("   _coastlines done. Took:",datetime.now()-t0)

		return h

	def _coast_patch(self):
		# Lazily create coast patch only if needed.
		if self._coast_patch_ is None:
			if self._coast_path is None:
				return None
			polys = [xy.exterior.coords for xy in 
			         self._coast_path.geoms]
			self._coast_patch = PathPatch(Path.make_compound_path(*PolyCollection(polys)
			                                                      .get_paths()), 
			                              edgecolor='none', facecolor='none',
			                              zorder=-2)
			#self._ax.add_collection(pc)
			self._ax.add_patch(self._coast_patch)

		return self._coast_patch

	def _scatter(self, lon, lat, **kwargs):
		"""
		Plot markers on map.
		"""

		# Checks:
		if isinstance(lon,list):
			lon = np.array(lon)
		if isinstance(lat,list):
			lat = np.array(lat)

		# Convert coordinates:
		x,y = self._projection.project(lon,lat)

		# Obtain plot coordinates:
		x,y = self._plot_canvas.obtain_coordinates(x, y, self._xlim, self._ylim)

		# Plot marker:
		h = self._ax.scatter(x, y, clip_on=True, **kwargs)
		h.set_clip_path(self._clip_rect)

		return h


	def _schedule_callback(self):
		"""
		The callback!
		"""

		# 0) If axes are not initialized, we need not
		#    continue:
		if not self._has_initialized_axes():
			return
		
		# 1) Check if we need to adjust axes:
		need_readjust = self._canvas_change()
		need_readjust = need_readjust or self._update_axes or self._update_grid
		if need_readjust:
			# Also delete existing coast line path:
			self._coast_path = None
			self._coast_path_prep = None
			self._coast_patch_ = None

		# 2) Determine all grid ticks:
		if need_readjust:
			self._determine_grid_ticks()

		# 3) Plot ticks:
		if need_readjust:
			self._plot_axes()

		# 4) Plot grid:
		if need_readjust:
			self._plot_grid()

		# 5) Iterate over everything and plot:
		for job in self._scheduled:
			if need_readjust or not job[1]:
				self._plot(job[0],job[2])
				job[1] = True


	def _canvas_change(self):
		"""
		This method is called whenever the data situation
		has changed and we need to update the axes.
		"""

		# See whether we have to adjust ticks:
		if self._data_xlim is None and self._user_xlim is None:
			# Nothing to be done!
			return False

		need_readjust = not self._adjusted

		if self._user_xlim is None:
			need_readjust = need_readjust or not np.array_equal(self._data_xlim,
			                                                    self._xlim)
			self._xlim = self._data_xlim
		else:
			need_readjust = need_readjust or not np.array_equal(self._user_xlim,
			                                                    self._xlim)
			self._xlim = self._user_xlim

		if self._user_ylim is None:
			need_readjust = need_readjust or not np.array_equal(self._data_ylim,
			                                                    self._ylim)
			self._ylim = self._data_ylim
		else:
			need_readjust = need_readjust or not np.array_equal(self._user_ylim,
			                                                    self._ylim)
			self._ylim = self._user_ylim

		if not need_readjust:
			# Nothing to be done!
			return False

		# Set axes aspect:
		self._calculate_canvas_size()
		if self._axes_aspect > 0:
			self._ax2.xaxis.tick_top()
			pos = self._ax_pos
			if self._aspect > 1.0:
				pos = (pos[0], pos[1], pos[3] / self._axes_aspect, pos[3])
			elif self._aspect < 1.0:
				pos = (pos[0], pos[1], pos[2], self._axes_aspect * pos[2])
			self._ax.set_position(pos)
			self._ax2.set_position(pos)
			# Resize figure if wished:
			if self._resize_figure:
				if self._axes_aspect > 1.0:
					self._fig.set_size_inches((self._figsize_original[1]/self._axes_aspect,
					                           self._figsize_original[1]))
				elif self._axes_aspect < 1.0:
					self._fig.set_size_inches((self._figsize_original[0],
					                           self._axes_aspect * 
					                           self._figsize_original[1]))
				else:
					self._fig.set_size_inches(self._figsize_original)

		# Compute ticks:
		self._tick_dict = self._projection.generate_ticks(self._xlim, self._ylim, 1.0)

		# Determine geographic extents:
		self._geographic_extents = \
		    self._projection.maximum_geographic_extents(self._xlim, self._ylim)

		# Reschedule coast lines:
		self._coasts = None

		# Save state:
		self._adjusted = True

		# Re-schedule everything:
		return True

	def _determine_grid_ticks(self):
		# Determine all grid ticks visible in the canvas and save them
		# in self._grid_lons and self._grid_lats.
		# TODO : This fails if the longitudes wrap around 360°.
		#        This could be solved by creating a GeographicExtents class.
		dlon = self._geographic_extents[0][1] - self._geographic_extents[0][0]
		dlat = self._geographic_extents[1][1] - self._geographic_extents[1][0]
		n_lons = int(np.ceil(dlon/self._grid_constant))+2
		n_lats = int(np.ceil(dlat/self._grid_constant))+2
		i0_lon = int(np.floor(self._geographic_extents[0][0] / self._grid_constant))
		i0_lat = int(np.floor(self._geographic_extents[1][0] / self._grid_constant))
		lons = (self._grid_constant * (np.arange(n_lons) + i0_lon)
		        + self._grid_anchor[0] + 180.0) % 360.0 - 180.0
		lats = (self._grid_constant * (np.arange(n_lats) + i0_lat)
		        + self._grid_anchor[1] + 90.0) % 180.0 - 90.0
		self._grid_lons = lons
		self._grid_lats = lats


	def _determine_tick_labels(self, tick_vals):
		# TODO Increase sophisticatedness!
		self._label_sign = 'label'
		# Calculate degrees, minutes, seconds:
		values = np.array([val[0] for val in tick_vals])
		sign = np.sign(values)
		values *= sign
		latlon_switch = np.array([val[1] for val in tick_vals], dtype=int)
		degrees = np.round(values).astype(int)
		values -= degrees
		minutes = np.round(values / 60.0).astype(int)
		values -= 60 * minutes
		seconds = np.round(values / 3600.0).astype(int)

		# Base labels:
		labels = []
		prefix = "$" if self._use_latex else ""
		degree_symbol = "^\\circ" if self._use_latex else "°"
		minute_symbol = "'"
		second_symbol = "\""
		space = "\\," if self._use_latex else " "
		suffix = "$" if self._use_latex else ""
		for i in range(len(values)):
			if seconds[i] == 0 and minutes[i] == 0:
				labels += [str(degrees[i]) + degree_symbol]
			elif seconds[i] == 0:
				labels += [str(degrees[i]) + degree_symbol + space + str(minutes[i])
					       + minute_symbol]
			else:
				labels += [str(degrees[i]) + degree_symbol + space + str(minutes[i])
					       + minute_symbol + space + str(seconds[i]) + second_symbol]

		if self._label_sign == 'label':
			for i in range(len(values)):
				if tick_vals[i][1] == 0:
					if sign[i] < 0:
						labels[i] = prefix + labels[i] + space + \
						            ("\\mathrm{W}" if self._use_latex else "W") + suffix
					else:
						labels[i] = prefix + labels[i] + space + \
						            ("\\mathrm{E}" if self._use_latex else "E") + suffix
				else:
					if sign[i] < 0:
						labels[i] = prefix + labels[i] + space + \
						            ("\\mathrm{S}" if self._use_latex else "S") + suffix
					else:
						labels[i] = prefix + labels[i] + space + \
						            ("\\mathrm{N}" if self._use_latex else "N") + suffix
		elif self._label_sign == 'sign':
			for i in range(len(values)):
				if sign[i] == -1:
					labels[i] = prefix + "-" + labels[i] + suffix

		return labels


	def _plot_axes(self):
		"""
		This method draws the axes artists.
		"""
		
		# TODO linewidth!
		linewidth = self._box_axes_linewidth
		
		# First, clear all:
		self._ax.clear()
		# Hide some axis stuff used for usual plotting:
		self._ax.spines['top'].set_visible(False)
		self._ax.spines['bottom'].set_visible(False)
		self._ax.spines['left'].set_visible(False)
		self._ax.spines['right'].set_visible(False)
		self._ax.tick_params(top=False, bottom=False, left=False, right=False)
		self._ax.set_xlim([self._canvas.x0,self._canvas.x1])
		self._ax.set_ylim([self._canvas.y0,self._canvas.y1])
		self._ax2.set_xlim([self._canvas.x0,self._canvas.x1])
		self._ax2.set_ylim([self._canvas.y0,self._canvas.y1])
		self._ax2.tick_params(top=False, bottom=False, left=False, right=False)
		canvas = self._canvas

		# Now determine how much space we need:
		# TODO

		# Generate tick arrays:
		tick_arrays = _create_tick_arrays(self._tick_dict,self._which_ticks)

		# Plot box axes if wished:
		if self._box_axes:
			# Draw box axes as polygons.
			axes_boxes, colors, canvas = \
			    _generate_axes_boxes(tick_arrays, self._xlim, self._ylim,
			                         self._box_axes_width, canvas, linewidth)
			self._ax.add_collection(PatchCollection(axes_boxes, facecolors=colors,
			                                        edgecolors='k',
			                                        zorder=self._zorder_axes_0,
			                                        linewidth=linewidth))
			tick_mask = None
		else:
			# Otherwise plot ticks:
			axes_ticks_xy, canvas, tick_mask = \
			    _generate_axes_ticks(tick_arrays, self._grid_lons, self._grid_lats,
			                         self._xlim, self._ylim,
			                         canvas, self._projection, self._box_axes_width,
			                         linewidth)
			if axes_ticks_xy is not None:
				self._ax.add_collection(LineCollection(np.concatenate(axes_ticks_xy,
				                                                      axis=0),
				                                       zorder=self._zorder_axes_0,
				                                       linewidth=linewidth,
				                                       color='k'))

		# Plot ticks:
		if axes_ticks_xy is not None:
			for i in range(4):
				# Obtain relevant ticks:
				ta = tick_arrays[i]
				if len(ta) == 0:
					continue
				if tick_mask is not None:
					ta = ta[tick_mask[i],:]

				# Obtain tick values in geographic coordinates:
				if i < 2:
					inv = self._projection.inverse(ta[:,0],
					                               self._ylim[i]*np.ones(ta.shape[0]))
				else:
					inv = self._projection.inverse(self._xlim[i-2]*np.ones(ta.shape[0]),
					                               ta[:,0])

				tick_vals = [(inv[int(ta[k,1])][k], int(ta[k,1]))
				             for k in range(ta.shape[0])]

				# Obtain tick positions in canvas coordinates:
				x = axes_ticks_xy[i][:,1,0 if i < 2 else 1]

				# Set axes ticks:
				labels = self._determine_tick_labels(tick_vals)
				if i == 0:
					self._ax.set_xticks(x)
					self._ax.set_xticklabels(labels)
				elif i == 1:
					self._ax2.set_xticks(x)
					self._ax2.set_xticklabels(labels)
				elif i == 2:
					self._ax.set_yticks(x)
					self._ax.set_yticklabels(labels)
				elif i == 3:
					self._ax2.set_yticks(x)
					self._ax2.set_yticklabels(labels)

		# The remaining canvas can be plotted on:
		self._plot_canvas = canvas

		# Obtain clipping rectangle:
		xclip, yclip = self._plot_canvas.obtain_coordinates(
		                        np.array([self._xlim[0], self._xlim[1]]),
		                        np.array([self._ylim[0], self._ylim[1]]),
		                        self._xlim, self._ylim)
		self._xclip = xclip
		self._yclip = yclip

		# Add border rect to plot. This is also to make sure that the
		# transformation is initialized:
		self._clip_rect = Rectangle([xclip[0],yclip[0]], xclip[1]-xclip[0],
			                        yclip[1]-yclip[0],edgecolor='k',
			                        facecolor='none',
			                        zorder=self._zorder_axes_0+1)
		self._clip_box = box(xclip[0],yclip[0],xclip[1],yclip[1])
		self._ax.add_artist(self._clip_rect)


	def _plot_grid(self):
		# Delete old grid:
		for h in self._grid_handles:
			h.remove()

		if self._verbose > 1:
			print("plot grid ...")
			print("   extents:",self._geographic_extents)
		# Plot grid:
		lons = self._grid_lons
		lats = self._grid_lats
		n_lons = lons.size
		n_lats = lats.size
		gridlines = []

		# Create ticks between:
		tb = self._grid_ticks_between
		lats_with_between = np.zeros((n_lats-1) * tb + n_lats)
		for i in range(n_lats-1):
			lats_with_between[i*(tb+1):(i+1)*(tb+1)+1] = np.linspace(lats[i],lats[i+1],
			                                                       tb+2)
		lons_with_between = np.zeros((n_lons-1) * tb + n_lons)
		for i in range(n_lons-1):
			lons_with_between[i*(tb+1):(i+1)*(tb+1)+1] = np.linspace(lons[i],lons[i+1],
			                                                       tb+2)

		# Create all grid coordinates:
		ones = np.ones_like(lats_with_between)
		for i in range(n_lons):
			# Project and transform to canvas coordinates
			x,y = self._projection.project(lons[i]*ones, lats_with_between)
			x,y = self._plot_canvas.obtain_coordinates(x, y, self._xlim, self._ylim)

			# Ignore outside points:
			if np.all(x < self._xclip[0]) or np.all(x > self._xclip[1]) \
			or np.all(y < self._yclip[0]) or np.all(y > self._yclip[1]):
				continue

			# Split lines at backside of sphere if needed:
			jump_x = identify_jumps(x, self._xclip)
			jump_y = identify_jumps(y, self._yclip)
			jump = np.logical_or(jump_x,jump_y)
			jump[-1] = False
			if np.any(jump):
				# Split:
				ids = np.argwhere(jump).flatten()+1
				x = np.split(x,ids)
				y = np.split(y,ids)
				for j in range(len(x)):
					if x[j].size <= 1:
						continue

					# Save all in one:
					half_meridian = np.zeros((x[j].size,2))
					half_meridian[:,0] = x[j]
					half_meridian[:,1] = y[j]
					gridlines += [half_meridian]
			else:
				# Save all in one:
				half_meridian = np.zeros((x.size,2))
				half_meridian[:,0] = x
				half_meridian[:,1] = y
				gridlines += [half_meridian]

		ones = np.ones_like(lons_with_between)
		for i in range(n_lats):
			# Project and transform to canvas coordinates
			x,y = self._projection.project(lons_with_between,lats[i]*ones)
			x,y = self._plot_canvas.obtain_coordinates(x, y, self._xlim, self._ylim)

			# Ignore outside points:
			if np.all(x < self._xclip[0]) or np.all(x > self._xclip[1]) \
			or np.all(y < self._yclip[0]) or np.all(y > self._yclip[1]):
				continue

			# Split lines at backside of sphere if needed:
			jump_x = identify_jumps(x, self._xclip)
			jump_y = identify_jumps(y, self._yclip)
			jump = np.logical_or(jump_x,jump_y)
			jump[-1] = False
			if np.any(jump):
				# Split:
				ids = np.argwhere(jump).flatten()+1
				x = np.split(x,ids)
				y = np.split(y,ids)
				for j in range(len(x)):
					if x[j].size <= 1:
						continue

					# Save all in one:
					col = np.zeros((x[j].size,2))
					col[:,0] = x[j]
					col[:,1] = y[j]
					gridlines += [col]
			else:
				# Save all in one:
				col = np.zeros((x.size,2))
				col[:,0] = x
				col[:,1] = y
				gridlines += [col]

		gridlines = split(MultiLineString(gridlines), self._clip_box)
		filter_box = prep(self._clip_box)
		filter_ = filter(filter_box.contains, gridlines)

		gridlines = LineCollection([np.array(g.xy).T for g in filter_], 
		                           clip_path=self._clip_rect, clip_on=True,
		                           **self._grid_kwargs)
		h = self._ax.add_collection(gridlines)
		h.set_clip_path(self._clip_rect)
		self._grid_handles = [h]

		if self._verbose > 1:
			print(" ... done!")

	def _get_axes_space(self, ax):
		"""
		Return the space an axis ('bot', 'top', 'left', or 'right')
		occupies.
		"""
		return self._box_axes_width + self._box_axes_linewidth / 72.


	def _calculate_canvas_size(self):
		"""
		Calculate the axes aspect.
		"""
		dy = self._ylim[1]-self._ylim[0]
		dx = self._xlim[1]-self._xlim[0]

		# For both, determine the ratio of projection to display
		# coordinates:
		rx = (self._ax_dx - self._get_axes_space("left") - self._get_axes_space("right")) / dx
		ry = (self._ax_dy - self._get_axes_space("top") - self._get_axes_space("bot")) / dy

		# The minimum ratio determines the size of the canvas:
		if rx < ry:
			ax_dy_target = rx * dy + self._get_axes_space("top") + \
			               self._get_axes_space("bot")
			ax_dx_target = self._ax_dx
		else:
			ax_dx_target = ry * dx + self._get_axes_space("left") + \
			               self._get_axes_space("right")
			ax_dy_target = self._ax_dy

		# Obtain target axes aspect:
		self._axes_aspect = ax_dy_target / ax_dx_target

		# Obtain canvas size:
		self._canvas = Rect(0, 0, ax_dx_target, ax_dy_target)
		self._plot_canvas = self._canvas


	def _plot(self, cmd, args):
		"""
		Main plotting code is here!
		"""
		if cmd == "imshow":
			self._imshow(*args[0:3],**args[3])
		elif cmd == "coastline":
			self._coastline(*args[0:2],**args[2])
		elif cmd == "scatter":
			self._scatter(args[0],args[1],**args[2])
		elif cmd == "quiver":
			self._quiver(*args[0:5], **args[5])
		elif cmd == "streamplot":
			self._streamplot(*args[0:6], **args[6])

		# Reset aspect:
		self._ax.set_aspect(self._aspect)
