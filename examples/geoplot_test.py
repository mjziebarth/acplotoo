#!/bin/python
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from acplotoo import Geoplot
from acplotoo.projection import HotineObliqueMercator
from acplotoo.sphere import azimuth
import numpy as np
from matplotlib import rcParams
from matplotlib import rc

from datetime import datetime

# Measure script execution time:
t0 = datetime.now()

rcParams['text.usetex'] = True
rc('font',family='serif')

# Create projection:
lon = 11.1539
lat = 54.4565
lat = 57.4565
projection = HotineObliqueMercator(lon, lat, 60)

# Create a path that we can visualize using quiver:
path_lon = np.array([5.0,   6.0,  6.0,  7.0,  8.0,  9.0, 10.0])
path_lat = np.array([54.5, 55.0, 56.7, 57.0, 57.5, 58.0, 58.5])
angles = [azimuth(path_lon[i], path_lat[i], path_lon[i+1], path_lat[i+1])
          for i in range(len(path_lon)-1)]
U = np.sin(np.deg2rad(angles))
V = np.cos(np.deg2rad(angles))


# Create plot:
fig = plt.figure(figsize=(3.37, 3.37))
fig = plt.figure(figsize=(7.0, 7.0))
ax = fig.add_subplot(111)
gplt = Geoplot(ax, projection, limits_xy=[[-1000,1000],[-100,100]],
               gshhg_path='gshhg-bin-2.3.7/gshhs_i.b', verbose=1, use_joblib=True,
               resize_figure=True)


gplt.set_xlim([-500e3,500e3])
gplt.set_ylim([-250e3,250e3])
gplt.grid(grid_constant=2.0)

gplt.scatter([11.1539, 11.1539, 12., lon, 9.], [54.4565, 55.,54.4565, lat, 52.], marker='*', zorder=2)

# Plot coastline:
gplt.coastline(4)

gplt.scatter([10.0],[56.0], marker='s', zorder=21, color='r')

# Plot the path using quiver:
gplt.quiver(path_lon[:-1], path_lat[:-1], U, V)



fig.savefig('geoplot_test.svg')

print("took:",datetime.now()-t0)
