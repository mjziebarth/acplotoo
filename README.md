![Spherical sketch](examples/readme_image.svg)

# plot_tools
A sporadically growing collection of code snippets for plotting and creating sketches with matplotlib.

## Current features:
* Spherical sketches
* Obtaining a small set of colors from a color map


## Supplies:
The **plot_tools** module which exposes the following methods:

### Top level functions
```python
    get_cm_colors(colormap, n_colors):
        # Obtain a number 'n_colors' of colors from a colormap
        # named 'colormap'.

```

### Sphereplot class
A class to plot and draw on a sphere. Designed to visualize
spherical geometry.
```python
   plot_tools.Sphereplot(ax, view_center)
      # Init a spherical plot on axis 'ax' with a top-down view
      # centered at 'view_center' (in lon/lat coordinates)
   
      great_circle(lon1, lat1, lon2, lat2, **kwargs)
         # Plot a great circle fixed at two points.
      
      scatter(lon, lat, **kwargs)
         # Scatter plot on a sphere
      
      line(lon1, lat1, lon2, lat2, **kwargs)
         # Plot a geodesic between two points.
      
      triangle(lon0, lat0, lon1, lat1, lon2, lat2, **kwargs)
         # Plot a triangle between three points.
      
      disk(lon, lat, r, radius_angle=None, **kwargs)
         # Plot a disk of radius 'r' around a point.
         # If 'radius_angle' is given, plot a radius-indicating
         # line from the center to the disk border.
      
      disk_sector(lon, lat, r, azi0, azi1, mode='sector', **kwargs)
         # Plot a circular sector (mode='sector') or circular
         # segment (mode='segment') on the sphere. 'azi0' and
         # 'azi1' denote the starting and ending azimuths
         # respectively.
      
      disk_intersection(lon1, lat1, lon2, lat2, r)
         # Plot the intersection area of two disks on the
         # sphere.
      
      arc_segment(lon, lat, r, azi0, azi1, **kwargs)
         # Plot an arc segment. Parameters refer to same
         # properties as disk_sector.
      
      bounds(lonmin, lonmax, latmin, latmax, **kwargs)
         # Visualize a longitude/latitude coordinate bound
         # interval on the sphere.
      
      wireframe(lon_ticks=18, lat_ticks=11, ticks_between=10,
                vc_override=None, **kwargs)
         # Plot a wireframe to visualize the sphere.
         # The wireframe can be rotated relative to other
         # plots using an override visual center.
      
      project(lon, lat, three_d=False)
         # Project lon/lat coordinates to axis coordinates.
```


## Requires:
* matplotlib
* numpy

## Installation:
Install with pip inside the root folder:
```bash
    pip install .
```

## Example code
This features the [example plot](examples/readme_image.py) at
the top and shows how to use the Sphereplot class and how to
obtain color map themed colors using the ```YlGnBu``` color map.

First, do the required imports from matplotlib, numpy, and
plot_tools:
```python
   import numpy as np
   import matplotlib.pyplot as plt
   import matplotlib as mpl

   from plot_tools import get_cm_colors, Sphereplot
   from plot_tools.sphere import displace, azimuth
```
Secondly, note that the **Sphereplot** class is designed as a
wrapper on the matplotlib.axes.Axes class. A typical plot
script will thus begin by setting up a matplotlib figure and
axis and creating a Sphereplot instance on the axis:
```python
   fig = plt.figure(figsize=(6.74, 3.37))
   ax1 = fig.add_subplot(121)
   ax1.set_position([0.025,0.05, 0.45, 0.9])
   plot1 = Sphereplot(ax1, view_center=vc)
```
Using ```get_cm_colors```, obtain some ```YlGnBu``` themed
colors:
```python
   n_colors = 5
   colors = get_cm_colors('YlGnBu_r',n_colors)
```
Then, set up your spherical geometry. Some helpful methods
are included:
