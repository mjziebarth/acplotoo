# -*- coding: utf-8 -*-
#
# Acplotoo geoplot helpers.
# Copyright (C) 2018 Malte Ziebarth
# 
# This software is distributed under the MIT license.
# See the LICENSE file in this repository.


from matplotlib.transforms import Affine2D
from matplotlib.textpath import TextPath
from matplotlib.patches import PathPatch
from matplotlib import rcParams
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable


def rotated_text(text, pos, rotation, anchor, color, zorder, size=None,
                 margin_rel=0.15, edgecolor='none'):
	"""
	Create a rotated text a certain position.

	Keyword arguments:
	   anchor : Anchor position, one of 'll', 'ul',
	            'lr', and 'ur', refering to 'lower left'
	            etc.
	"""
	if size is None:
		size = rcParams['font.size']

	# Step 1: Create the text and measure it:
	path = TextPath((0,0), text, size=size/72.0, usetex=rcParams['text.usetex'])
	xmin = path.vertices[:,0].min()
	xmax = path.vertices[:,0].max()
	ymin = path.vertices[:,1].min()
	ymax = path.vertices[:,1].max()
	width = xmax - xmin
	height = ymax - ymin

	margin = margin_rel * height

	# Step 2: Move the text such that the anchor is at (0,0):
	
	if anchor == 'll':
		offset = ( -xmin+margin, -ymin+margin)
	elif anchor == 'ul':
		offset = (-xmin+margin, -ymax-margin)
	elif anchor == 'lr':
		offset = (-xmax-margin, -ymin+margin)
	elif anchor == 'ur':
		offset = (-xmax-margin, -ymax-margin)
	transform = Affine2D().translate(offset[0],offset[1])

	# Now rotate the text and afterwards move to
	# position:
	transform.rotate_deg(-rotation).translate(pos[0],pos[1])

	return PathPatch(path.transformed(transform), facecolor=color,
	                 edgecolor=edgecolor, zorder=zorder)



def field2color(field, cmap, **kwargs):
	"""
	Convert a scalar field to RGBA values using
	a color map.
	"""
	kwargs = dict(kwargs)
	# Extract possibly relevant keyword arguments:
	vmin = kwargs.pop("vmin",None)
	vmax = kwargs.pop("vmax",None)

	# Normalization:
	if vmin is None:
		vmin = field.min()
	if vmax is None:
		vmax = field.max()
	norm = Normalize(vmin, vmax)
	c = norm(field)

	# Create scalar mappable and obtain the color values
	# from color map:
	mappable = ScalarMappable(norm, cmap)
	mappable.set_array(field)

	return mappable.to_rgba(field), mappable
