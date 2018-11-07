# -*- coding: utf-8 -*-
#
# Plot tools GeoplotBase class file.
# Copyright (C) 2018 Malte Ziebarth
# 
# This software is distributed under the MIT license.
# See the LICENSE file in this repository.

from .backend import _generate_axes_boxes
from .rect import Rect


from matplotlib.collections import PatchCollection
from matplotlib.patches import Polygon, Rectangle
import struct
import numpy as np

# TODO fix!
from matplotlib.transforms import Bbox

# Utility:


def gshhg_read_header(bytes,i0):
	# Read a GSHHG header:
	i, n, flag, west, east, south, north, area, \
	area_full,container,ancestor \
	    = struct.unpack('>IIIiiiiIIII', bytes[i0:i0+44])
	
	# Return:
	return i0+44, (i, n, flag, west, east, south, north, area, area_full,
	               container, ancestor)

class GeoplotBase:
	def __init__(self, ax, projection, gshhg_path):
		
		# TODO : Checks!
		self._ax = ax
		self._projection = projection
		
		# Initialize axes:
		ax.set_axis_off()
		
		# See if we can load GSHHG:
		if gshhg_path is None:
			self._gshhg_path = None
		else:
			self._gshhg_path = gshhg_path
			self._read_gshhg()
		
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


	def _read_gshhg(self):
		if self._gshhg_path is None:
			raise RuntimeError("GSHHG not loaded!")
		
		with open(self._gshhg_path, 'rb') as f:
			bytes = f.read()
		
		print("Reading coastlines...")
		
		N = len(bytes)
		self._all_coasts = []
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
			
			# Convert to xy and save:
			x, y = self._projection.project(latlon[:,0],latlon[:,1])
			
			self._all_coasts += [(level, x, y)]

		print("done!")


	def _imshow(self, z, xlim, ylim, **kwargs):
		# Actually we would have to add half the gridsize for
		# pixel-registered grid, but for small grids it's
		# okay for now.
		xlim, ylim = self._plot_canvas.obtain_coordinates(xlim, ylim,
		                      self._xlim, self._ylim)
		self._ax.imshow(z, extent=(xlim[0],xlim[1],ylim[0],ylim[1]),
		                **kwargs)

	def _coastline(self, level, **kwargs):
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
		
		# Obtain clipping rectangle:
		print("\n\nOBTAIN CLIPPING RECTANGLE")
		xclip, yclip = self._plot_canvas.obtain_coordinates(
		                        np.array([self._xlim[0], self._xlim[1]]),
		                        np.array([self._ylim[0], self._ylim[1]]),
		                        self._xlim, self._ylim)
		clip_bbox = Bbox(np.array([[xclip[i],yclip[i]] for i in [0,1]]))
		# TODO Make this rect invisible!
		clip_rect = Rectangle([xclip[0],yclip[0]], xclip[1]-xclip[0],
			                   yclip[1]-yclip[0],edgecolor='none',alpha=0.8,zorder=10)
		print("\n\n")
		
		# Of those, plot all polygons with level <= level:
		patches = []
		for poly in self._coasts:
			if poly[0] <= level:
				# Obtain coordinates in canvas coordinates:
				print("poly[1]:",[poly[1].min(),poly[1].max()])
				print("poly[2]:",[poly[2].min(),poly[2].max()])
				x,y = self._plot_canvas.obtain_coordinates(poly[1],poly[2], self._xlim,
				                                           self._ylim)

				print("have Polygon in bounds x (rel):[",x.min(),",",x.max(),
				      "], y (rel):[",y.min(),",",y.max(),"]")


				# Create polygon:
				patches += [Polygon(np.concatenate([x[:,np.newaxis],y[:,np.newaxis]],
				                                   axis=1),
				#                                  clip_box=clip_bbox,
				#                                  clip_on=True,
				                                   **kwargs)]
				print("last polygon:",patches[-1].xy)
				print(patches[-1].get_clip_box())
				print(patches[-1].get_clip_on())
				print("clip_box:    ",clip_bbox)

		# Plot all polygons:
		if len(patches) > 0:
			print(" ---- ADDING PATCH COLLECTION ----")
			#print("patches:",patches)
			self._ax.add_artist(clip_rect)
			h = self._ax.add_collection(PatchCollection(patches,facecolors='none',
			                                            edgecolors='k',
			                                            clip_path=clip_rect,
			                                            clip_on=True))
			h.set_clip_path(clip_rect)
			print("h.cp:",h.get_clip_path())
			print("xlim:",self._ax.get_xlim())
			print("ylim:",self._ax.get_ylim())
			print("bbox: ",clip_bbox)
			print("rect: ",clip_rect)
			print("xclip:",xclip)
			print("yclip:",yclip)
		else:
			print("xclip:",xclip)
			print("yclip:",yclip)
			print("bbox: ",clip_bbox)
			print("len(self._coasts):",len(self._coasts))


	def _schedule_callback(self):
		"""
		The callback!
		"""
		
		# 1) Check if we need to adjust axes:
		if self._adjust_axes():
			# Axes have been adjusted, replot everything!
			for job in self._scheduled:
				job[1] = False
		
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

		print("self._user_xlim:",self._user_xlim)
		print("self._user_ylim:",self._user_ylim)

		if not need_readjust:
			# Nothing to be done!
			return False

		# Compute ticks:
		tick_dict = self._projection.generate_ticks(self._xlim, self._ylim, 1.0)

		# Plot ticks:
		self._plot_axes(tick_dict)
		
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
		print("cleared axis!")
		canvas = self._canvas
		print("canvas:",canvas)

		# Now determine how much space we need:
		# TODO

		# Plot ticks:
		# TODO

		# Plot box axes if wished:
		if self._box_axes:
			# Draw box axes as polygons.
			axes_boxes, colors, canvas = \
			    _generate_axes_boxes(tick_dict, self._xlim, self._ylim,
			                         self._box_axes_width, canvas, linewidth)
			self._ax.add_collection(PatchCollection(axes_boxes, facecolors=colors,
			                                        edgecolors='k'))
		print("plot_canvas:",canvas)
		
		# The remaining canvas can be plotted on:
		self._plot_canvas = canvas


	def _plot(self, cmd, args):
		"""
		Main plotting code is here!
		"""
		if cmd == "imshow":
			self._imshow(*args)
		elif cmd == "coastline":
			self._coastline(*args)