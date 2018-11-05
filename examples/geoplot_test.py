#!/bin/python
import matplotlib.pyplot as plt
from plot_tools import Geoplot
from plot_tools.projection import HotineObliqueMercator



# Create projection:
projection = HotineObliqueMercator(120, 10, -10)

# Create plot:
fig = plt.figure(figsize=(3.37, 3.37))
ax = fig.add_subplot(111)
gplt = Geoplot(ax, projection, limits_xy=[[-1000,1000],[-100,100]], gshhg_path='gshhg-bin-2.3.7/gshhs_c.b')

# Plot coastline:
gplt.set_xlim([-1000,1000])
gplt.set_ylim([-100,100])

gplt.coastline(1)


fig.savefig('geoplot_test.svg')
