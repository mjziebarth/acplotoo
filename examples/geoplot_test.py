#!/bin/python
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from plot_tools import Geoplot
from plot_tools.projection import HotineObliqueMercator
import numpy as np



# Create projection:
lon = 11.1539
lat = 54.4565
lat = 57.4565
projection = HotineObliqueMercator(lon, lat, 60)

# Create plot:
fig = plt.figure(figsize=(3.37, 3.37))
fig = plt.figure(figsize=(7.0, 7.0))
ax = fig.add_subplot(111)
ax.set_position([0.01,0.01,0.98,0.98])
gplt = Geoplot(ax, projection, limits_xy=[[-1000,1000],[-100,100]], gshhg_path='gshhg-bin-2.3.7/gshhs_i.b')




# Plot coastline:
gplt.set_xlim([-500e3,500e3])
gplt.set_ylim([-500e3,500e3])

gplt.scatter([11.1539, 11.1539, 12., lon, 9.], [54.4565, 55.,54.4565, lat, 52.], marker='*', zorder=2)

gplt.coastline(4)




fig.savefig('geoplot_test.svg')
