# -*- coding: utf-8 -*-
#
# Plot tools geoplot rect file.
# Copyright (C) 2018 Malte Ziebarth
# 
# This software is distributed under the MIT license.
# See the LICENSE file in this repository.

import numpy as np

class Rect():
	
	def __init__(self,x0,y0,x1,y1):
		self.x0 = x0
		self.y0 = y0
		self.x1 = x1
		self.y1 = y1
		self.dx = x1-x0
		self.dy = y1-y0
	
	def __str__(self):
		return "((" + str(self.x0) + "," + str(self.y0) + "), ("\
		       + str(self.x1) + "," + str(self.y1) + "))"
	
	def from_bbox(Bbox):
		return Rect(Bbox.bounds[0],Bbox.bounds[1],
		            Bbox.bounds[0]+Bbox.bounds[2],
		            Bbox.bounds[0]+Bbox.bounds[3])
	
	
	def strip_margin(self, margin):
		if isinstance(margin,float) or isinstance(margin,np.float64):
			return Rect(self.x0-margin, self.y0-margin, self.x1-margin,
			            self.y1-margin)
		else:
			if len(margin) != 4:
				raise RuntimeError("Rect.strip_margin(): 'margin' has to be "
				                   "float or a 4-element iterable!")
			return Rect(self.x0-margin[0], self.y0-margin[1], self.x1-margin[2],
			            self.y1-margin[3])
	
	def obtain_coordinates(self, x, y, xlim, ylim):
		print("xlim:",xlim)
		print("x:",type(x))
		x_ = self.x0 + self.dx * (x - xlim[0]) / (xlim[1]-xlim[0])
		y_ = self.y0 + self.dy * (y - ylim[0]) / (ylim[1]-ylim[0])
		return x_, y_
	
	def width(self):
		return self.dx
	
	def height(self):
		return self.dy
