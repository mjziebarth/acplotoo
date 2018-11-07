#!/bin/python
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from plot_tools import Geoplot
from plot_tools.projection import HotineObliqueMercator
import numpy as np



# Create projection:
lon = 11.1539
lat = 54.4565
projection = HotineObliqueMercator(lon, lat, 90)

# Create plot:
fig = plt.figure(figsize=(3.37, 3.37))
ax = fig.add_subplot(111)
gplt = Geoplot(ax, projection, limits_xy=[[-1000,1000],[-100,100]], gshhg_path='gshhg-bin-2.3.7/gshhs_l.b')


# Plot coastline:
gplt.set_xlim([-500e3,500e3])
gplt.set_ylim([-500e3,500e3])

gplt.coastline(1)




fig.savefig('geoplot_test.svg')
