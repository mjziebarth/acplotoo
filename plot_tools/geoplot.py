# -*- coding: utf-8 -*-
#
# Plot tools Geoplot class file.
# Copyright (C) 2018 Malte Ziebarth
# 
# This software is distributed under the MIT license.
# See the LICENSE file in this repository.

from .geoplot_base.rect import Rect
from .geoplot_base.base import GeoplotBase





# GEOPLOT:

class Geoplot(GeoplotBase):


	def __init__(self, ax, projection, limits_xy=None, gshhg_path=None,
	             which_ticks='significant', water_color='lightblue',
	             land_color='white', verbose=0, use_joblib=False):
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

		super().__init__(ax, projection, gshhg_path, which_ticks,
		                 water_color, land_color, verbose, use_joblib)

		self._gshhg_path = gshhg_path


		self._canvas = None
		self._plot_canvas = self._canvas

		# Setup configuration:
		self._box_axes = True
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

	def coastline(self, level, water_color=None, land_color=None,
	              zorder=0, **kwargs):
		"""
		Plot the coast line.
		"""
		if self._gshhg_path is None:
			raise RuntimeError("GSHHG not loaded!")

		if water_color is not None:
			self._water_color = water_color

		if land_color is not None:
			self._land_color = land_color


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

	def imshow_projected(self, z, xlim, ylim, **kwargs):
		"""
		Plot a field (in projected coordinates) using imshow.
		"""
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
