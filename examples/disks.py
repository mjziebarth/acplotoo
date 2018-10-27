#!/bin/python
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

from plot_tools import get_cm_colors, Sphereplot
from plot_tools.sphere import displace, azimuth


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


### Plot some random disks on a sphere. ###

# Create plot:
fig = plt.figure(figsize=(3.37, 3.37))
ax = fig.add_subplot(111)
ax.set_position([0.025,0.025, 0.95, 0.95])
plot = Sphereplot(ax, view_center=vc)


# Get some colors from a colormap:
n_colors = 20
colors = get_cm_colors('YlGnBu_r',n_colors)


# Coordinates and geometry:
N = 100
np.random.seed(2000)
lon = 360.0 * np.random.random(N)
lat = np.rad2deg(np.arcsin(2*np.random.random(N)-1.0))

# Use varying radii and colors:
R=7.0 * (1.5*np.random.random(N) + 1.)


# Plot:
plot.wireframe(linewidth=0.5*LW, color='gray', lon_ticks=10, lat_ticks=9, zorder=1,
                ticks_between=20)

# Plot disks around the points:
for i in range(N):
	plot.disk(lon[i], lat[i], R[i], seg_len=0.5, facecolor=colors[i % n_colors],
	          zorder=(1.0-R[i]/R.max()))

# Print plot:
fig.savefig("disks.svg")
