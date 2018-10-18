#!/bin/python
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

from plot_tools import get_cm_colors, Sphereplot
from plot_tools.sphere import displace


# Use some nice LaTeX:
mpl.rcParams["text.usetex"] = True
mpl.rc('font',family='serif')
mpl.rcParams['text.latex.preamble'].append(r'\usepackage{amsfonts}')
mpl.rcParams.update({'font.size' : 9})


# Configuration:
vc = (55.0,25)
LW = 0.8


# Create plot:
fig = plt.figure(figsize=(3.37, 3.37))
ax = fig.add_subplot(111)
ax.set_position([0.05,0.05, 0.9, 0.9])
plot = Sphereplot(ax, view_center=vc)



# Get some colors from a colormap:
n_colors = 25
colors2 = get_cm_colors('magma',n_colors)

# Use some of them for qualitative plotting:
colors = colors2[9::5,:]


# Coordinates and geometry:
lon0 = 35.0
lat0 = 13.0
azimuth = 37.0
R = 40.0
d01 = 60.0
lon1,lat1 = displace(lon0, lat0, azimuth, d01)


# Plot:
plot.wireframe(linewidth=0.5*LW, color='gray', lon_ticks=10, lat_ticks=9, zorder=1,
               ticks_between=20)

# Plot a disk:
plot.disk(lon0, lat0, R, seg_len=1.5, facecolor=colors[3], zorder=0)

# Overlay the disk with a series of sectors:
delta_angle = 45.0 / n_colors
for i in range(n_colors):
	angle = 45.0 - i*delta_angle
	plot.disk_sector(lon0, lat0, R, azimuth-angle, azimuth+angle, seg_len=2.5,
		             facecolor=colors2[-(i+1)], zorder=0)


# Plot scatter marks at the centers:
plot.scatter([lon0,lon1],[lat0,lat1],marker='.',c='k', zorder=3)

# Plot a great circle through the points:
plot.great_circle(lon0, lat0, lon1, lat1, color='gray')





# Print plot:
fig.savefig("readme_image.svg")

