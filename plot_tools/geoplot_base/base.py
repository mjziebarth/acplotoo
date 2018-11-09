# -*- coding: utf-8 -*-
#
# Plot tools GeoplotBase class file.
# Copyright (C) 2018 Malte Ziebarth
# 
# This software is distributed under the MIT license.
# See the LICENSE file in this repository.

from .backend import _generate_axes_boxes, _create_tick_arrays, has_joblib,\
                     identify_jumps, read_coastlines, coast_line_patches_and_path
from .rect import Rect


from matplotlib.collections import PatchCollection, LineCollection
from matplotlib.patches import Polygon, Rectangle
from shapely.geometry.polygon import Polygon as SPoly
from shapely.geometry import LineString, MultiPolygon, box, MultiLineString, MultiPolygon
from shapely.ops import polygonize_full, unary_union
from shapely.prepared import prep
from shapely import speedups
import numpy as np


from datetime import datetime

# TODO fix!
from matplotlib.transforms import Bbox

# Enable speedups:
if speedups.available:
	speedups.enable()


# Utility:



class GeoplotBase:
	def __init__(self, ax, projection, gshhg_path, which_ticks,
	             water_color, land_color, verbose, use_joblib):

		# TODO : Checks!
		self._ax = ax
		self._projection = projection

		self._which_ticks = which_ticks
		self._water_color = water_color
		self._land_color = land_color

		if use_joblib and not has_joblib:
			raise ImportError("Could not import joblib but was "
			                  "requested to do so!")
		self._use_joblib = use_joblib


		# Initialize axes:
		ax.set_axis_off()

		# See if we can load GSHHG:
		if gshhg_path is None:
			self._gshhg_path = None
		else:
			self._gshhg_path = gshhg_path

		# Setup internal data:
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
		self._clip_rect = Rectangle([0.0, 0.0], 1.0, 1.0)
		self._grid_on = True
		self._grid_constant = 1.0
		self._grid_kwargs_base = {'color' : 'gray', 'linewidth' : 0.5}
		self._grid_kwargs = {**self._grid_kwargs_base}
		self._grid_anchor = (0.0,0.0)
		self._verbose = verbose
		self._adjusted = False
		self._grid_handles = []
		self._tick_dict = None
		self._update_axes = False
		self._update_grid = False
		self._box_axes_linewidth = 0.5

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
		xlim, ylim = self._plot_canvas.obtain_coordinates(xlim, ylim,
		                      self._xlim, self._ylim)
		h = self._ax.imshow(z, extent=(xlim[0],xlim[1],ylim[0],ylim[1]),
		                    **kwargs)
		h.set_clip_path(self._clip_rect)

		return h


	def _coastline(self, level, zorder, **kwargs):
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

		print("self._coasts.size:",len(self._coasts))

		cnvs_x = np.array([self._plot_canvas.x0, self._plot_canvas.x1])
		cnvs_y = np.array([self._plot_canvas.y0, self._plot_canvas.y1])

		# Call backend to generate patches and path clipped to canvas.
		# The backend employs joblib if wished and possible to load for
		# repeated calls.
		patches_xy, coast_path =\
		    coast_line_patches_and_path(coords, cnvs_x, cnvs_y, self._xclip,
		                                self._yclip)

		# Prepare coast path geometry and create patches:
		self._coast_path = prep(coast_path)
		patches = [Polygon(xy, **kwargs) for xy in patches_xy]


		# Plot all polygons:
		if len(patches) > 0:
			colors = [[self._water_color, self._land_color][c[2] % 2] for c in coords]
			h = self._ax.add_collection(PatchCollection(patches,facecolors=colors,
				                                        edgecolors='k', zorder=zorder))
			h.set_clip_path(self._clip_rect)
		else:
			h = None

		if self._verbose > 0:
			print("   _coastlines done. Took:",datetime.now()-t0)

		return h

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
		need_readjust = need_readjust or self._update_axes
		if need_readjust:
			# Also delete existing coast line path:
			self._coast_path = None

		# 2) Plot ticks:
		if need_readjust:
			self._plot_axes()

		# 3) Plot grid:
		if need_readjust or self._update_grid:
			self._plot_grid()

		# 4) Iterate over everything and plot:
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
		#self._aspect = (self._ylim[1]-self._ylim[0]) / (self._xlim[1]-self._xlim[0])
		self._ax.set_aspect(self._axes_aspect)

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


	def _plot_axes(self):
		"""
		This method draws the axes artists.
		"""
		
		# TODO linewidth!
		linewidth = self._box_axes_linewidth
		
		# First, clear all:
		self._ax.clear()
		self._ax.set_axis_off()
		self._ax.set_xlim([self._canvas.x0,self._canvas.x1])
		self._ax.set_ylim([self._canvas.y0,self._canvas.y1])
		canvas = self._canvas
		
		print("\n\ncanvas:",canvas)

		# Now determine how much space we need:
		# TODO

		# Generate tick arrays:
		tick_arrays = _create_tick_arrays(self._tick_dict,self._which_ticks)

		# Plot ticks:
		# TODO

		# Plot box axes if wished:
		if self._box_axes:
			# Draw box axes as polygons.
			axes_boxes, colors, canvas = \
			    _generate_axes_boxes(tick_arrays, self._xlim, self._ylim,
			                         self._box_axes_width, canvas, linewidth)
			self._ax.add_collection(PatchCollection(axes_boxes, facecolors=colors,
			                                        edgecolors='k',zorder=20))

		# The remaining canvas can be plotted on:
		self._plot_canvas = canvas
		print("_plot_canvas:",canvas)

		# Obtain clipping rectangle:
		xclip, yclip = self._plot_canvas.obtain_coordinates(
		                        np.array([self._xlim[0], self._xlim[1]]),
		                        np.array([self._ylim[0], self._ylim[1]]),
		                        self._xlim, self._ylim)
		self._xclip = xclip
		self._yclip = yclip

		print("xclip:",xclip)
		print("yclip:",yclip)
		print("\n\n")

		# Invisible rect added to plot. This is just to make sure that the
		# transformation is initialized:
		# TODO : Make sure this is actually invisible!
		self._clip_rect = Rectangle([xclip[0],yclip[0]], xclip[1]-xclip[0],
			                        yclip[1]-yclip[0],edgecolor='none',
			                        facecolor=self._water_color,
			                        zorder=-1)
		self._clip_box = box(xclip[0],yclip[0],xclip[1],yclip[1])
		self._ax.add_artist(self._clip_rect)


	def _plot_grid(self):
		# Delete old grid:
		for h in self._grid_handles:
			h.remove()

		if self._verbose > 0:
			print("plot grid ...")
			print("   extents:",self._geographic_extents)
		# Plot grid:
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
		gridlines = []
		ones = np.ones_like(lats)
		for i in range(n_lons):
			# Project and transform to canvas coordinates
			x,y = self._projection.project(lons[i]*ones,lats)
			x,y = self._plot_canvas.obtain_coordinates(x, y, self._xlim, self._ylim)

			# Ignore outside points:
			outside = np.logical_or(np.logical_or(x < self._xclip[0],
			                                      x > self._xclip[1]),
			                        np.logical_or(y < self._yclip[0],
			                                      y > self._yclip[1]))
			if np.all(outside):
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
				half_meridian = np.zeros((n_lats,2))
				half_meridian[:,0] = x
				half_meridian[:,1] = y
				gridlines += [half_meridian]

		ones = np.ones_like(lons)
		for i in range(n_lats):
			# Project and transform to canvas coordinates
			x,y = self._projection.project(lons,lats[i]*ones)
			x,y = self._plot_canvas.obtain_coordinates(x, y, self._xlim, self._ylim)

			# Ignore outside points:
			outside = np.logical_or(np.logical_or(x < self._xclip[0],
			                                      x > self._xclip[1]),
			                        np.logical_or(y < self._yclip[0],
			                                      x > self._yclip[1]))
			if np.all(outside):
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
				col = np.zeros((n_lons,2))
				col[:,0] = x
				col[:,1] = y
				gridlines += [col]

		gridlines = MultiLineString(gridlines)
		gridlines = self._clip_box.intersection(gridlines)

		gridlines = LineCollection([np.array(g.xy).T for g in gridlines.geoms], 
		                           clip_path=self._clip_rect, clip_on=True,
		                           **self._grid_kwargs)
		h = self._ax.add_collection(gridlines)
		h.set_clip_path(self._clip_rect)
		self._grid_handles = [h]

		if self._verbose > 0:
			print(" ... done!")

	def _get_axes_space(self, ax):
		"""
		Return the space an axis ('bot', 'top', 'left', or 'right')
		occupies.
		"""
		return self._box_axes_width + self._box_axes_linewidth


	def _calculate_canvas_size(self):
		"""
		Calculate the axes aspect.
		"""
		dy = self._ylim[1]-self._ylim[0]
		dx = self._xlim[1]-self._xlim[0]
		aspect = dy / dx
		
		# Obtain the axes size in inches:
		size = self._ax.get_window_extent().transformed(self._ax.get_figure()
		                                                .dpi_scale_trans.inverted())
		ax_dx = size.x1 - size.x0
		ax_dy = size.y1 - size.y0
		
		# For both, determine the ratio of projection to display
		# coordinates:
		rx = (ax_dx - self._get_axes_space("left") - self._get_axes_space("right")) / dx
		ry = (ax_dy - self._get_axes_space("top") - self._get_axes_space("bot")) / dy
		
		# The minimum ratio determines the size of the canvas:
		if rx < ry:
			ax_dy_target = rx * dy + self._get_axes_space("top") + self._get_axes_space("bot")
			ax_dx_target = ax_dx
		else:
			ax_dx_target = ry * dx + self._get_axes_space("left") + self._get_axes_space("right")
			ax_dy_target = ax_dy
		
		# Obtain target axes aspect:
		self._axes_aspect = ax_dy_target / ax_dx_target
		
		# Obtain canvas size:
		print("size:",size)
		
		self._canvas = Rect(0, 0, ax_dx_target, ax_dy_target)
		#self._canvas = Rect(0,0,1,1)
		# TODO
		self._plot_canvas = self._canvas


	def _plot(self, cmd, args):
		"""
		Main plotting code is here!
		"""
		print("args:",args)
		if cmd == "imshow":
			self._imshow(*args[0:3],**args[3])
		elif cmd == "coastline":
			self._coastline(*args[0:2],**args[2])
		elif cmd == "scatter":
			self._scatter(args[0],args[1],**args[2])
		
		# Reset aspect:
		self._ax.set_aspect(self._aspect)
