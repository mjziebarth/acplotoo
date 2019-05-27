import numpy as np
import matplotlib.pyplot as plt

from matplotlib.patches import Polygon, PathPatch
from matplotlib.textpath import TextPath
from matplotlib.path import Path
from matplotlib.transforms import Affine2D
from matplotlib import rcParams

# Points of the fourarrow:
outer = ((-6.0,0.0), (0.0,6.0), (6.0,0.0), (0.0,-6.0))
inner = ((-1.0,1.0), (1.0,1.0), (1.0,-1.0), (-1.0,-1.0))
direction = ['W','N','E','S']

def plot_north_arrow(ax, pos0, size, angle, flipped, color, zorder=4):
	# Create polygon data. First the triangles:
	triangles = [((0,0), outer[i], inner[i]) for i in range(4)]

	# Then the boundary polygon:
	boundary = [None]*8
	boundary[::2] = outer
	boundary[1::2] = inner

	# Norm to area (-1,1)x(-1,1)
	boundary = np.array(boundary) / 6.0
	triangles = np.array(triangles) / 6.0

	# Create the general transform:
	rotation = Affine2D()
	rotation.rotate_deg(-angle)
	rotation.translate(*pos0)

	# Plot the polygons:
	triangles = 0.5 * size * triangles
	boundary = 0.5 * size * boundary

	for poly in triangles:
		path = Polygon(poly, closed=True).get_path()
		h = ax.add_patch(PathPatch(path.transformed(rotation), facecolor=color,
		                           edgecolor='none', zorder=zorder))

	path = Polygon(boundary, closed=True).get_path()
	ax.add_patch(PathPatch(path.transformed(rotation), linestyle='-', facecolor='none',
	                     edgecolor=color, zorder=zorder))

	# Labels of directions:
	fig = ax.get_figure()
	renderer = ax.figure.canvas.get_renderer()
	for i in range(4):
		_dir = direction[i]
		if flipped and i in (0,2):
			if _dir == 'W':
				_dir = 'E'
			elif _dir == 'E':
				_dir = 'W'
		pos = (size*outer[i][0]/12.0, size*outer[i][1]/12.0)
		path = TextPath(pos, _dir, size=rcParams['font.size']/72.0,
		                usetex=rcParams['text.usetex'])
		width = path.vertices[:,0].max() - path.vertices[:,0].min()
		height = path.vertices[:,1].max() - path.vertices[:,1].min()
		transform = Affine2D()
		if i == 0:
			# West:
			transform.translate(-1.3*width, -0.5*height)
			if flipped:
				# West offset does not work for East label.
				transform.translate(-0.3*width, 0)
		elif i ==3:
			# South:
			transform.translate(-0.5*width, -1.4*height)
		elif i == 1:
			# North:
			transform.translate(-0.5*width, 0.4*height)
		elif i == 2:
			# East:
			transform.translate(0.3*width, -0.5*height)
		transform.rotate_deg(-angle)
		transform.translate(*pos0)
		h = ax.add_patch(PathPatch(path.transformed(transform),
		                           facecolor=color, edgecolor='none', zorder=zorder))
