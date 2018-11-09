# -*- coding: utf-8 -*-
#
# Plot tools GeoplotBase class file.
# Copyright (C) 2018 Malte Ziebarth
# 
# This software is distributed under the MIT license.
# See the LICENSE file in this repository.

from .backend import _generate_axes_boxes, _create_tick_arrays
from .rect import Rect


from matplotlib.collections import PatchCollection, LineCollection
from matplotlib.patches import Polygon, Rectangle
from shapely.geometry.polygon import Polygon as SPoly
from shapely.geometry import LineString, MultiPolygon, box, MultiLineString, MultiPolygon
from shapely.ops import polygonize_full, unary_union
from shapely.prepared import prep
from shapely import speedups
import struct
import numpy as np


from datetime import datetime

# TODO fix!
from matplotlib.transforms import Bbox

# Enable speedups:
if speedups.available:
	speedups.enable()


# Utility:


def gshhg_read_header(bytes,i0):
	# Read a GSHHG header:
	i, n, flag, west, east, south, north, area, \
	area_full,container,ancestor \
	    = struct.unpack('>IIIiiiiIIII', bytes[i0:i0+44])
	
	# Return:
	return i0+44, (i, n, flag, west, east, south, north, area, area_full,
	               container, ancestor)

def identify_jumps(c, lim):
	# Identify points in a line where a jump between left
	# and right or top and bottom happens outside the plot
	# canvas:
	c1 = np.roll(c,-1)
	cmin = np.minimum(c,c1)
	cmax = np.maximum(c,c1)
	jump = np.logical_and(cmin < lim[0], cmax > lim[1])

	# Return boolean array that is True at the index right
	# before the jump:
	return jump



class GeoplotBase:
	def __init__(self, ax, projection, gshhg_path, which_ticks,
	             water_color, land_color, verbose):

		# TODO : Checks!
		self._ax = ax
		self._projection = projection

		self._which_ticks = which_ticks
		self._water_color = water_color
		self._land_color = land_color

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
		self._ticks = None
		self._scheduled = []
		self._coasts = None
		self._all_coasts = None
		self._clip_rect = Rectangle([0.0, 0.0], 1.0, 1.0)
		self._grid_on = True
		self._grid_constant = 1.0
		self._grid_kwargs = {'color' : 'gray'}
		self._grid_anchor = (0.0,0.0)
		self._verbose = verbose

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

		with open(self._gshhg_path, 'rb') as f:
			bytes = f.read()

		if self._verbose > 0:
			print("Reading coastlines...")
			t0 = datetime.now()

		N = len(bytes)
		self._all_coasts = []
		raw_data = []
		i=0
		while i < N:
			# Read header:
			i, res = gshhg_read_header(bytes,i)

			# Obtain level and number of points:
			n = res[1]
			level = res[2] % 256

			# Read points:
			xy = np.frombuffer(bytes[i:i+8*n],dtype=np.dtype('>i4'),count=2*n)
			i += 8*n

			# Create lat/lon array:
			latlon = 1e-6*xy.reshape((n,2))

			raw_data += [(level, latlon)]

		if self._verbose > 1:
			t1 = datetime.now()
			print("took:",t1-t0)
			print("Filtering coastline coordinates...")

		filtered_data = []
		for d in raw_data:
			outside = np.logical_or(
			             np.logical_or(d[1][:,0] < self._geographic_extents[0][0],
			                           d[1][:,0] > self._geographic_extents[0][1]),
			             np.logical_or(d[1][:,1] < self._geographic_extents[1][0],
			                           d[1][:,1] > self._geographic_extents[1][1])
			                        )
			if not np.all(outside):
				filtered_data += [d]

		if self._verbose > 1:
			t2 = datetime.now()
			print("   took:",t2-t1)
			print("   remaining:",len(filtered_data))
			print("Converting coastline coordinates...")

		for dat in filtered_data:
			# Convert to xy and save:
			x, y = self._projection.project(dat[1][:,0], dat[1][:,1])

			self._all_coasts += [(dat[0], x, y)]

		if self._verbose > 1:
			t3 = datetime.now()
			print("Coastlines read! Took:",t3-t2)


	def _imshow(self, z, xlim, ylim, **kwargs):
		# Actually we would have to add half the gridsize for
		# pixel-registered grid, but for small grids it's
		# okay for now.
		xlim, ylim = self._plot_canvas.obtain_coordinates(xlim, ylim,
		                      self._xlim, self._ylim)
		self._ax.imshow(z, extent=(xlim[0],xlim[1],ylim[0],ylim[1]),
		                **kwargs)


	def _coastline(self, level, **kwargs):
		if self._verbose > 0:
			print("_coastline ...")

		if self._all_coasts is None \
		  or self._xlim[0] < self._coast_xlim[0] \
		  or self._xlim[1] > self._coast_xlim[1] \
		  or self._ylim[0] < self._coast_ylim[0] \
		  or self._ylim[1] > self._coast_ylim[1]:
			# Read GSHHS:
			self._read_gshhg()
			self._coast_xlim = self._xlim
			self._coast_ylim = self._ylim

		t0 = datetime.now()
		# See if a cached version exists:
		if self._coasts is None:
			# We have to recalculate.
			self._coasts = []
			for poly in self._all_coasts:
				if np.any(np.logical_and(
				              np.logical_and(poly[1] >= self._xlim[0],
				                             poly[1] <= self._xlim[1]),
				              np.logical_and(poly[2] >= self._ylim[0],
				                             poly[2] <= self._ylim[1]))
				          ):
					# At least one point is inside canvas:
					self._coasts += [poly]

		# Of those, plot all polygons with level <= level:
		coords = []
		for poly in self._coasts:
			if poly[0] <= level:
				# Obtain coordinates in canvas coordinates:
				x,y = self._plot_canvas.obtain_coordinates(poly[1],poly[2], self._xlim,
				                                           self._ylim)
				coords += [(x,y,poly[0])]

		# Optimize the paths:
		patches = []
		coast_path = []
		for i in range(len(coords)):
			c = coords[i]
			cnvs_x = [self._plot_canvas.x0, self._plot_canvas.x1]
			cnvs_y = [self._plot_canvas.y0, self._plot_canvas.y1]
			outside = np.logical_or(np.logical_or(c[0] < cnvs_x[0], c[0] > cnvs_x[1]),
			                        np.logical_or(c[1] < cnvs_y[0], c[1] > cnvs_y[1]))
			if np.any(outside):
				# Split polygons at points where a coordinate jump over the plot
				# are is performed. We assume this can only happen by wrapping around
				# back of sphere. This assumption should be good if the polygons are
				# 'reasonably' sampled.
				
				both_outside = np.logical_and(outside, np.roll(outside,1))
				jump_x = identify_jumps(c[0], cnvs_x)
				jump_y = identify_jumps(c[1],cnvs_y)
				
				split = np.logical_and(both_outside, np.logical_or(jump_x,jump_y))
				
				# Split polygons at those points:
				XY = np.concatenate([c[0][:,np.newaxis],
				                     c[1][:,np.newaxis]],axis=1)
				XY = np.split(XY, np.argwhere(split).flatten()+1, axis=0)
				
				# Add sub polygons:
				for xy in XY:
					patches += [Polygon(xy, **kwargs)]

					# Add polygon for sea shore to global polygon:
					if c[2] == 1 and xy.shape[0] > 2:
						# This is a level 1 polygon, AKA sea border!
						# Polygonize, thanks to Mike T on stackoverflow!
						closed = np.concatenate([xy,xy[0:2]],axis=0)
						ls = LineString(closed)
						mls = unary_union(ls)
						coast_path += polygonize_full(mls)[0]
			else:
				# If all inside, all is well!

				# Create patch:
				xy = np.concatenate([c[0][:,np.newaxis],c[1][:,np.newaxis]], axis=1)
				patches += [Polygon(xy, **kwargs)]

				# Add polygon for sea shore to global polygon:
				if c[2] == 1 and xy.shape[0] > 2:
					spoly = SPoly(xy)
					coast_path += [spoly]

		# Plot all polygons:
		if len(patches) > 0:
			colors = [[self._water_color, self._land_color][c[2] % 2] for c in coords]
			h = self._ax.add_collection(PatchCollection(patches,facecolors=colors,
			                                            edgecolors='k',
			                                            clip_path=self._clip_rect,
			                                            clip_on=True))
			h.set_clip_path(self._clip_rect)

		# Finalize coast path:
		coast_path = MultiPolygon([self._clip_box.intersection(poly)
		                           for poly in coast_path])
		coast_path = unary_union(coast_path)

		self._coast_path = prep(coast_path)

		if self._verbose > 0:
			print("   _coastlines done. Took:",datetime.now()-t0)


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
		self._ax.scatter(x,y, clip_path=self._clip_rect, clip_on=True, **kwargs)


	def _schedule_callback(self):
		"""
		The callback!
		"""
		
		# 1) Check if we need to adjust axes:
		if self._adjust_axes():
			# Axes have been adjusted, replot everything!
			for job in self._scheduled:
				job[1] = False
			
			# Also delete existing coast line path:
			self._coast_path = None

		# 2) If axes are not initialized, we need not
		#    continue:
		if not self._has_initialized_axes():
			return
		
		# 2) Iterate over everything and plot:
		for job in self._scheduled:
			if not job[1]:
				self._plot(job[0],job[2])
				job[1] = True


	def _adjust_axes(self):
		"""
		This method is called whenever the data situation
		has changed and we need to update the axes.
		"""

		# See whether we have to adjust ticks:
		if self._data_xlim is None and self._user_xlim is None:
			# Nothing to be done!
			return False
		
		need_readjust = False

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

		# Compute ticks:
		tick_dict = self._projection.generate_ticks(self._xlim, self._ylim, 1.0)

		# Plot ticks:
		self._plot_axes(tick_dict)

		# Determine geographic extents:
		self._geographic_extents = \
		    self._projection.maximum_geographic_extents(self._xlim, self._ylim)

		# Plot grid:
		self._plot_grid()

		# Reschedule coast lines:
		self._coasts = None

		# Re-schedule everything:
		return True


	def _plot_axes(self, tick_dict):
		"""
		This method draws the axes artists.
		"""
		
		# TODO linewidth!
		linewidth = 0.005
		
		# First, clear all:
		self._ax.clear()
		self._ax.set_axis_off()
		self._ax.set_xlim([0.0,1.0])
		self._ax.set_ylim([0.0,1.0])
		canvas = self._canvas

		# Now determine how much space we need:
		# TODO

		# Generate tick arrays:
		tick_arrays = _create_tick_arrays(tick_dict,self._which_ticks)

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

		# Obtain clipping rectangle:
		xclip, yclip = self._plot_canvas.obtain_coordinates(
		                        np.array([self._xlim[0], self._xlim[1]]),
		                        np.array([self._ylim[0], self._ylim[1]]),
		                        self._xlim, self._ylim)
		self._xclip = xclip
		self._yclip = yclip

		# Invisible rect added to plot. This is just to make sure that the
		# transformation is initialized:
		# TODO : Make sure this is actually invisible!
		self._clip_rect = Rectangle([xclip[0],yclip[0]], xclip[1]-xclip[0],
			                        yclip[1]-yclip[0],edgecolor='none',
			                        facecolor=self._water_color,
			                        zorder=0)
		self._clip_box = box(xclip[0],yclip[0],xclip[1],yclip[1])
		self._ax.add_artist(self._clip_rect)


	def _plot_grid(self):
		if self._verbose > 0:
			print("plot grid ...")
		# Plot grid:
		n_lons = int(np.floor(360.0/self._grid_constant))
		n_lats = int(np.floor(180.0/self._grid_constant))+1
		lons = (self._grid_constant * np.arange(n_lons) + self._grid_anchor[0]) % 360.0
		lats = (self._grid_constant * np.arange(n_lats) + self._grid_anchor[1]
		        + 90.0) % 180.0 - 90.0
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

		if self._verbose > 0:
			print(" ... done!")
			print("\n\n")


	def _plot(self, cmd, args):
		"""
		Main plotting code is here!
		"""
		if cmd == "imshow":
			self._imshow(*args)
		elif cmd == "coastline":
			self._coastline(*args)
		elif cmd == "scatter":
			print("args:",args)
			self._scatter(args[0],args[1],**args[2])
