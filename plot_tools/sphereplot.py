# -*- coding: utf-8 -*-
#
# Main class of the sphereplot Python module.
# Copyright (C) 2018 Malte Ziebarth
# 
# This software is distributed under the MIT license.
# See the LICENSE file in this repository.
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.collections import PolyCollection, LineCollection
from matplotlib.patches import Polygon
from matplotlib.axes import Axes

# Library imports:
from .geometry import rotation_matrix, rotate_vectors,\
    convert_coordinates_3d, convert_coordinates,\
    great_circle_distance, _gc_rotation_angle, _line_coords\

from .helpers import connect_masked_sequence


# The sphereplot class

class Sphereplot:
	#TODO docstring!
	"""
	A convenience class.
	
	"""
	
	def __init__(self, ax, view_center=(0.0,0.0)):
		#TODO docstring!
		"""
		Init method.
		"""
		
		# Make sure ax is a matplotlib Axes object:
		if not isinstance(ax,Axes):
			raise TypeError("Sphereplot has to be initialized with a "
			                "matplotlib.axes.Axes object!")
		
		# Make sure that view center is of right type:
		try:
			self.view_center = (float(view_center[0]),
			                    float(view_center[1]))
		except:
			raise TypeError("Sphereplot has to be initialized with "
			                "a lon/lat view_center tuple!")
		
		# Save ax:
		self.ax = ax
		
		# Init fake 3d view:
		ax.set_axis_off()
		ax.set_xlim([-1.01,1.01])
		ax.set_ylim([-1.01,1.01])
		ax.set_aspect('equal')
	
	
	def great_circle(self, lon1, lat1, lon2, lat2,
	                 tolerance=1e-8, **kwargs):
		"""
		Plot a great circle through the points (lon1, lat1)
		and (lon2, lat2)
		"""
		#TODO docstring!
		#TODO add segment_len keyword argument!
	
		if np.abs(lat1-90.0) < tolerance or np.abs(lat1+90.0) < tolerance:
			raise ValueError("Unhandled degeneracy: Point one is a pole!")
	
		# Start by creating a great circle through point 1
		# and poles:
		lons = np.concatenate([np.ones(50)*lon1,np.ones(50)*lon1+180.0])
		lats = np.concatenate([np.linspace(-90,90,50), np.linspace(90,-90,50)])
		x,y,z = convert_coordinates_3d(lons, lats, self.view_center)
	
		# Rotate the great circle around point one to point 2:
		axis = np.array(convert_coordinates_3d(lon1, lat1, self.view_center)).reshape((3,))
		angle = _gc_rotation_angle(lon1, lat1, lon2, lat2)
		x,y,z = rotate_vectors(x,y,z, axis, angle)
	
		# Mask to hide points on backside of sphere:
		mask = z >= 0.0
		
		# We know that half of the values are invisible. If both first and last
		# index are visible, the order is mixed up, since array is structured like
		# (v=visible, o=obscured):
		# [v, ..., v,o, ..., o, v, ..., v]
		x = connect_masked_sequence(x,mask)
		y = connect_masked_sequence(y,mask)
		lines = LineCollection([np.array([x,y]).T], **kwargs)
		return self.ax.add_collection(lines)
	
	
	def scatter(self, lon, lat, **kwargs):
		"""
		Scatter points on the sphere.
		"""
		x,y = convert_coordinates(np.array(lon), np.array(lat),
		                          self.view_center)
		return self.ax.scatter(x, y, **kwargs)
	
	
	def line(self, lon1, lat1, lon2, lat2, seg_len=5.0, 
	          **kwargs):
		"""
		Plot a great circle segment between the points (lon1, lat1)
		and (lon2, lat2)
		"""
		# Create the line coordinates:
		x,y,z = _line_coords(lon1, lat1, lon2, lat2, seg_len,
		                     self.view_center)
	
		# Mask:
		mask = z >= 0.0
		lines = LineCollection([np.array([x[mask],y[mask]]).T], **kwargs)
		return self.ax.add_collection(lines)
	
	
	def triangle(self, lon0, lat0, lon1, lat1, lon2, lat2,
	             seg_len=2.5, **kwargs):
		"""
		Plots a spherical triangle polygon.
		"""
		# Get coordinates of the three lines:
		x0,y0,z0 = _line_coords(lon0, lat0, lon1, lat1, seg_len,
		                        self.view_center)
		x1,y1,z1 = _line_coords(lon1, lat1, lon2, lat2, seg_len,
		                        self.view_center)
		x2,y2,z2 = _line_coords(lon2, lat2, lon0, lat0, seg_len,
		                        self.view_center)
	
		# Connect lines:
		x = np.concatenate([x0,x1,x2])
		y = np.concatenate([y0,y1,y2])
		z = np.concatenate([z0,z1,z2])
	
		# Mask:
		mask = z >= 0.0
		if not np.all(mask):
			x = x[mask]
			y = y[mask]
	
		# Plot:
		handles = []
		poly = Polygon(np.array([x,y]).T, **kwargs)
		handles += [self.ax.add_patch(poly)]
	
		return handles
	
	
	def disk(self, lon, lat, r, seg_len=2.5, radius_angle=None, **kwargs):
		"""
		Plot a disk centered on point (lon, lat) with radius r.
		"""
		D2R = np.pi/180.0
		
		# Create radius line:
		if radius_angle is not None:
			n = int(np.ceil(r/seg_len))
			xr,yr,zr = convert_coordinates_3d(np.zeros(n),
			               np.linspace(90.0, 90.0-r,n),
			               self.view_center)
			axis = np.array(convert_coordinates_3d(0.0, 90.0,
			                self.view_center)).reshape((3,))
			xr,yr,zr = rotate_vectors(xr,yr,zr, axis, radius_angle)
		
		# Create circle around north pole:
		n = int(np.ceil(360*np.sin(D2R*r)/seg_len))
		lons, lats = np.linspace(0,360,n), np.ones(n)*(90.0-r)
		x,y,z = convert_coordinates_3d(lons, lats, self.view_center)
		
		# Rotate circle to latitude:
		axis = np.array(convert_coordinates_3d(90.0, 0.0, self.view_center))\
		       .reshape((3,))
		x,y,z = rotate_vectors(x,y,z, axis, 90.0-lat)
		if radius_angle is not None:
			xr,yr,zr = rotate_vectors(xr,yr,zr, axis, 90.0-lat)
		
		# Rotate circle to longitude:
		axis = np.array(convert_coordinates_3d(0.0, 90.0, self.view_center))\
		       .reshape((3,))
		x,y,z = rotate_vectors(x,y,z, axis, lon)
		if radius_angle is not None:
			xr,yr,zr = rotate_vectors(xr,yr,zr, axis, lon)
		
		# Mask:
		mask = z < 0.0
		if np.any(mask):
			raise Exception("Error: Cannot handle obscured points in disk plotting "
			                "yet!")
		
		# Plot:
		handles = []
		poly = Polygon(np.array([x,y]).T,**kwargs)
		handles += [self.ax.add_patch(poly)]
		
		# Plot radius:
		if radius_angle is None:
			handles += [None]
		else:
			kwargs_radius = dict()
			if 'linewidth' in kwargs:
				kwargs_radius["linewidth"] = kwargs["linewidth"]
			lines = LineCollection([np.array([xr,yr]).T], color='k', **kwargs_radius)
			handles += [self.ax.add_collection(lines)]
		
		
		return handles
	
	
	def disk_sector(self, lon, lat, r, azi0, azi1, seg_len=2.5,
	                mode='sector', **kwargs):
		"""
		Plots a circular sector or segment on a sphere.
		
		
		Optional parameters:
		   mode : Choose the shape to plot. One of
		          'sector'  : arc + wedge
		          'segment' : connect ends of arc by line
		          Default: 'sector'
		
		"""
		# Create the circle section so that it only has to be rotated
		# in latitude direction:
		D2R = np.pi/180.0
		delta_lon = (azi1-azi0) % 360.0
		lon0 = (lon + 180.0 - azi0) % 360.0
		lon1 = lon0 - delta_lon
		# Number of segments:
		Nc = int(np.ceil(delta_lon / seg_len * np.sin(D2R*r)))
		# The coordinates of the circle section:
		circ_lon = np.mod(np.linspace(lon0, lon1, Nc), 360.0)
		circ_lat = (90.0-r) * np.ones(Nc)
	
		if mode == 'segment':
			# Convert to Euclidean coordinates:
			x,y,z = convert_coordinates_3d(circ_lon, circ_lat,
					                       self.view_center)
		
			# Connect the end points of the arc directly without going to
			# the center of the circle:
			x2,y2,z2 = _line_coords(lon1, 90.0-r, lon0, 90.0-r, seg_len,
			                        self.view_center)
		
			# Combined coordinates of the polygon:
			x = np.concatenate([x, x2])
			y = np.concatenate([y, y2])
			z = np.concatenate([z, z2])
		elif mode == 'sector':
			# The coordinates of the lines connecting the segments to the center:
			Nl = int(np.ceil(r / seg_len))
			l1_lon = lon0*np.ones(Nl)
			l1_lat = np.linspace(90.-r,90.,Nl)[::-1]
			l2_lon = lon1*np.ones(Nl)
			l2_lat = np.linspace(90.-r,90.,Nl)
	
			# Combined coordinates of the polygon:
			lon_poly = np.concatenate([l1_lon, circ_lon, l2_lon])
			lat_poly = np.concatenate([l1_lat, circ_lat, l2_lat])
	
			# Convert to Euclidean coordinates:
			x,y,z = convert_coordinates_3d(lon_poly, lat_poly,
					                       self.view_center)
	
		# Rotate to target latitude:
		axis = np.array(convert_coordinates_3d(lon-90.0, 0.0, self.view_center))\
		       .reshape((3,))
		x,y,z = rotate_vectors(x,y,z, axis, -(90.0-lat))
	
		# Mask:
		mask = z < 0.0
		if np.any(mask):
			raise Exception("Error: Cannot handle obscured points in disk section plotting "
				            "yet!")
	
		# Plot:
		handles = []
		poly = Polygon(np.array([x,y]).T,**kwargs)
		handles += [self.ax.add_patch(poly)]
	
		return handles
	
	
	def disk_intersection(self, lon1, lat1, lon2, lat2, r,
	                      seg_len=2.5, delta_r=0.5, hatchcolor=None, **kwargs):
		"""
		Plot the intersection of two disks.
		"""
		# Check whether an intersection exists:
		if great_circle_distance(lon1, lat1, lon2, lat2) > r:
			return
		
		D2R = np.pi/180.0
	
		# Create circle around north pole:
		n = int(np.ceil(360*np.sin(D2R*r)/seg_len))
		lons, lats = np.linspace(0,360,n), np.ones(n)*(90.0-r)
		x1,y1,z1 = convert_coordinates_3d(lons.copy(), lats.copy(),
		                                  self.view_center)
		x2,y2,z2 = convert_coordinates_3d(lons.copy(), lats.copy(),
		                                  self.view_center)
		xc1,yc1,zc1 = convert_coordinates_3d(lon1, lat1,
		                                     self.view_center)
		xc2,yc2,zc2 = convert_coordinates_3d(lon2, lat2,
		                                     self.view_center)
	
		# Rotate circles to latitude:
		axis = np.array(convert_coordinates_3d(90.0, 0.0, self.view_center))\
		       .reshape((3,))
		x1,y1,z1 = rotate_vectors(x1,y1,z1, axis, 90.0-lat1)
		x2,y2,z2 = rotate_vectors(x2,y2,z2, axis, 90.0-lat2)
	
		# Rotate circles to longitude:
		axis = np.array(convert_coordinates_3d(0.0, 90.0, self.view_center))\
		       .reshape((3,))
		x1,y1,z1 = rotate_vectors(x1,y1,z1, axis, lon1)
		x2,y2,z2 = rotate_vectors(x2,y2,z2, axis, lon2)
	
		# Combine circles:
		x,y,z = np.concatenate([x1,x2]),np.concatenate([y1,y2]),np.concatenate([z1,z2])
	
		# Filter out all points close enough to both centers:
		thresh_dotp = np.cos((r+delta_r) * D2R)
		dotp1 = x*xc1 + y*yc1 + z*zc1
		dotp2 = x*xc2 + y*yc2 + z*zc2
		mask = np.logical_and(dotp1 >= thresh_dotp, dotp2 >= thresh_dotp)
		x = connect_masked_sequence(x,mask)
		y = connect_masked_sequence(y,mask)
		z = connect_masked_sequence(z,mask)
		# Mask:
		mask = z < 0.0
		if np.any(mask):
			raise Exception("Error: Cannot handle obscured points in disk plotting "
				            "yet!")
	
		# Plot:
		handles = []
		vertices = np.array([x,y]).T
		hlw = kwargs.pop("hatchlinewidth",None)
		if hatchcolor is None or not "hatch" in kwargs or kwargs["hatch"] is None:
			poly = Polygon(vertices,**kwargs)
			handles += [self.ax.add_patch(poly)]
		else:
			kwargs1 = kwargs.copy()
			kwargs2 = kwargs.copy()
			kwargs2["edgecolor"] = hatchcolor
			if hlw is not None:
				mpl.rcParams['hatch.linewidth'] = hlw
			poly2 = Polygon(vertices,**kwargs,edgecolor=hatchcolor)
			handles += [self.ax.add_patch(poly2)]
			kwargs1["hatch"]=None
			poly1 = Polygon(vertices,**kwargs1)
			handles += [self.ax.add_patch(poly1)]
		
		return handles
	
	
	def arc_segment(self, lon, lat, r, azi0, azi1, seg_len=2.5,
	                **kwargs):
		"""
		Plots an arc segment on a sphere.
		"""
		# Create the circle section so that it only has to be rotated
		# in latitude direction:
		D2R = np.pi/180.0
		delta_lon = (azi1-azi0) % 360.0
		lon0 = (lon + 180.0 + azi0) % 360.0
		lon1 = lon0 + delta_lon
		# Number of segments:
		Nc = int(np.ceil(delta_lon / seg_len * np.sin(D2R*r)))
		# The coordinates of the circle section:
		circ_lon = np.mod(np.linspace(lon0, lon1, Nc), 360.0)
		circ_lat = (90.0-r) * np.ones(Nc)
	
		# Convert to Euclidean coordinates:
		x,y,z = convert_coordinates_3d(circ_lon, circ_lat,
			                           self.view_center)
	
		# Rotate to target latitude:
		axis = np.array(convert_coordinates_3d(lon-90.0, 0.0, self.view_center)).reshape((3,))
		x,y,z = rotate_vectors(x,y,z, axis, -(90.0-lat))
	
		# Mask:
		mask = z < 0.0
		if np.any(mask):
			raise Exception("Error: Cannot handle obscured points in disk section plotting "
				            "yet!")
	
		# Plot:
		return self.ax.plot(x, y, **kwargs)
	
	
	def wireframe(self, lon_ticks=18, lat_ticks=11, ticks_between=10,
	              vc_override=None, **kwargs):
		lon = np.linspace(0, 360, lon_ticks)
		lat = np.linspace(-90, 90, lat_ticks)
		lon = np.concatenate([lon, lon[-1]+ [lon[1]-lon[0]]])
		line_segments = []
		
		# View center:
		if vc_override is not None:
			vc = vc_override
		else:
			vc = self.view_center
		
		for i in range(lon_ticks):
			for j in range(lat_ticks):
				if j != 0 and j != lat_ticks-1:
					lon_loc = np.linspace(lon[i], lon[i+1], ticks_between)
					lat_loc = lat[j]*np.ones(ticks_between)
					x,y = convert_coordinates(lon_loc, lat_loc, vc)
					if len(x) > 1:
						line_segments += [np.array([x,y]).T]
				if j != lat_ticks-1:
					lat_loc = np.linspace(lat[j], lat[j+1], ticks_between)
					lon_loc = lon[i]*np.ones(ticks_between)
					x,y = convert_coordinates(lon_loc, lat_loc, vc)
					if len(x) > 1:
						line_segments += [np.array([x,y]).T]
		
		handles = []
		lines = LineCollection(line_segments, **kwargs)
		handles += [self.ax.add_collection(lines)]
	
		# Plot boundary:
		n_b = lon_ticks*ticks_between
		x,y = np.sin(np.linspace(0,2*np.pi,n_b)), np.cos(np.linspace(0,2*np.pi,n_b))
		handles += [self.ax.add_collection(LineCollection([np.array([x,y]).T], **kwargs))]
		
		return handles
	
	
	def project(self, lon, lat, three_d=False):
		"""
		Return projection coordinates of some lon-lat points.
		"""
		if three_d:
			return convert_coordinates_3d(lon, lat, self.view_center)
		
		return convert_coordinates(lon, lat, self.view_center)
