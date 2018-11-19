#!/bin/python
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

from acplotoo import get_cm_colors, Sphereplot
from acplotoo.sphere import displace, azimuth


# Reproducibility:
np.random.seed(777)

# Use some nice LaTeX:
mpl.rcParams["text.usetex"] = True
mpl.rc('font',family='serif')
mpl.rcParams['text.latex.preamble'].append(r'\usepackage{amsfonts}')
mpl.rcParams.update({'font.size' : 9})


# Configuration:
vc = (55.0,25)
LW = 0.8


### Plot 1: ###
# A triangle as the intersection of three great circles.

# Create plot 1:
fig = plt.figure(figsize=(6.74, 3.37))
ax1 = fig.add_subplot(121)
ax1.set_position([0.025,0.05, 0.45, 0.9])
plot1 = Sphereplot(ax1, view_center=vc)



# Get some colors from a colormap:
n_colors = 5
colors = get_cm_colors('YlGnBu_r',n_colors)


# Coordinates and geometry:
lon0 = 35.0
lat0 = 13.0

lon1 = 65.0
lat1 = 68.0

lon2 = 83.0
lat2 = -7.0

R=7.0


# Plot:
plot1.wireframe(linewidth=0.5*LW, color='gray', lon_ticks=10, lat_ticks=9, zorder=1,
                ticks_between=20)

# Plot disks around the points:
plot1.disk(lon0, lat0, R, seg_len=0.5, facecolor=colors[2], zorder=0)
plot1.disk(lon1, lat1, R, seg_len=0.5, facecolor=colors[2], zorder=0)
plot1.disk(lon2, lat2, R, seg_len=0.5, facecolor=colors[2], zorder=0)


# Plot a triangle:
plot1.triangle(lon0, lat0, lon1, lat1, lon2, lat2, color=colors[1], zorder=0)

# Plot three great circles:
plot1.great_circle(lon0, lat0, lon1, lat1, color=colors[0], linewidth=LW, zorder=3)
plot1.great_circle(lon1, lat1, lon2, lat2, color=colors[0], linewidth=LW, zorder=3)
plot1.great_circle(lon2, lat2, lon0, lat0, color=colors[0], linewidth=LW, zorder=3)

# Plot sectors that mark the angles of the triangle.
# To do so, first determine those angles in the global
# reference frame:
azi01 = azimuth(lon0, lat0, lon1, lat1)
azi02 = azimuth(lon0, lat0, lon2, lat2)
azi10 = azimuth(lon1, lat1, lon0, lat0)
azi12 = azimuth(lon1, lat1, lon2, lat2)
azi20 = azimuth(lon2, lat2, lon0, lat0)
azi21 = azimuth(lon2, lat2, lon1, lat1)

# Plot the sectors:
plot1.disk_sector(lon0, lat0, R, azi01, azi02, seg_len=0.5,
                  facecolor=colors[3], zorder=0)
plot1.disk_sector(lon1, lat1, R, azi12, azi10, seg_len=0.5,
                  facecolor=colors[3], zorder=0)
plot1.disk_sector(lon2, lat2, R, azi20, azi21, seg_len=0.5,
                  facecolor=colors[3], zorder=0)



### Plot 2: ###
# A subset of the sphere in lon/lat bounds and all points of a randomly
# distributed set that are in it.
ax2 = fig.add_subplot(122)
ax2.set_position([0.525,0.05, 0.45, 0.9])
plot2 = Sphereplot(ax2, view_center=vc)


# Configuration and generation of points:
lon_bounds = [37.0, 72.0]
lat_bounds = [3.0, 64.0]

N = 400
lons = 360.0*np.random.random(N)
lats = np.rad2deg(np.arcsin(2.0*np.random.random(N)-1.0))


# Plot:
plot2.wireframe(linewidth=0.5*LW, color='gray', lon_ticks=10, lat_ticks=9, zorder=1,
                ticks_between=20)

# Plot the bounds:
plot2.bounds(lon_bounds[0], lon_bounds[1], lat_bounds[0], lat_bounds[1],
             facecolor=colors[3], edgecolor='none', zorder=0)

# Plot the points. Points inside the bounds are highlighted
# by using a color of the color map, the others are indicated
# in gray:
c = 0.75*np.ones((N,3))
mask = np.logical_and(np.logical_and(lons >= lon_bounds[0], lons <= lon_bounds[1]),
                      np.logical_and(lats >= lat_bounds[0], lats <= lat_bounds[1]))
c[mask,:] = colors[2][:3]
plot2.scatter(lons, lats, c=c, s=3)

# Plot bounds again to clip points:
plot2.bounds(lon_bounds[0], lon_bounds[1], lat_bounds[0], lat_bounds[1],
             facecolor='none', edgecolor=colors[0], zorder=2)


# Print plot:
fig.savefig("readme_image.svg")
