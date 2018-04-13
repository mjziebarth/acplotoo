# -*- coding: utf-8 -*-
#
# Plot tools module file.
# Copyright (C) 2018 Malte Ziebarth
# 
# This software is distributed under the MIT license.
# See the LICENSE file in this repository.

import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cmx


def get_cm_colors(colormap, n_colors):
	""" 
	Obtain a number of colors from a colormap.

	Parameters:
	
	colormap: The name of the colormap.
	
	n_colors: The number of colors uniformly evenly
	          spaced along the colormap.
	"""
	cm = plt.get_cmap(colormap)
	cNorm  = colors.Normalize(vmin=0, vmax=n_colors-1)
	scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=cm)
	return scalarMap.to_rgba(range(n_colors))
	
