import geopandas as gpd
from shapely.geometry import Point, LineString, Polygon, MultiPoint
from os.path import join
import numpy as np
import matplotlib.pyplot as plt
import sys
from numpy import asarray

def distance_matrix(x0, y0, x1, y1):
    obs = np.vstack((x0, y0)).T
    interp = np.vstack((x1, y1)).T

    # Make a distance matrix between pairwise observations
    # Note: from <http://stackoverflow.com/questions/1871536>
    # (Yay for ufuncs!)
    d0 = np.subtract.outer(obs[:,0], interp[:,0])
    d1 = np.subtract.outer(obs[:,1], interp[:,1])

    return np.hypot(d0, d1)

# Setup: Generate data...
def simple_idw(x, y, z, xi, yi, p):
    dist = distance_matrix(x,y, xi,yi)

    # In IDW, weights are 1 / distance
    weights = 1.0 / dist ** p

    # Make weights sum to one
    weights /= weights.sum(axis=0)

    # Multiply the weights for each interpolated point by all observed Z-values
    zi = np.dot(weights.T, z)
    return zi
def plot(x,y,z,grid):
    plt.figure()
    plt.imshow(grid, extent=(x.min(), x.max(), y.max(), y.min()))
    # plt.hold(True)
    # plt.scatter(x,y,c=z)
    plt.colorbar()

gdf = gpd.read_file(join('shapefiles', 'nitrate_wgs84.shp'))

gdf['x'] = gdf.geometry.x
gdf['y'] = gdf.geometry.y
# xy = gdf[['x', 'y']].to_numpy()
# z = gdf.nitr_ran.to_numpy()

xmin, ymin, xmax, ymax = gdf.total_bounds

resolution = 50
# nx = ny = resolution
xi = np.linspace(xmin, xmax, resolution)
yi = np.linspace(ymin, ymax, resolution)
points = []
xi, yi = np.meshgrid(xi, yi)
positions = np.vstack([xi.ravel(), yi.ravel()])
for i in range(50):
   for j in range(50):
       point = yi[i, j], xi[i, j]
    #    point = Point(point)
       points.append(point)
    #    points.append((xi[i,j],yi[i,j]))
    #    Z[i,j] = f(X[i,j],Y[i,j])
# xi, yi = xi.flatten(), yi.flatten()

# ----------------------------------------------------------------------
# n = 10
# nx, ny = 50, 50
# x, y, z = map(np.random.random, [n, n, n])
x = gdf.geometry.x
y = gdf.geometry.y
z = gdf.nitr_ran

resolution = 100
# nx = ny = resolution
xi = np.linspace(xmin, xmax, resolution)
yi = np.linspace(ymin, ymax, resolution)
xi, yi = np.meshgrid(xi, yi)
xi, yi = xi.flatten(), yi.flatten()

# Calculate IDW
p = 2
grid1 = simple_idw(x,y,z,xi,yi, p)
grid1_reshape = grid1.reshape((resolution, resolution))
    
# Comparisons...
plot(x,y,z,grid1_reshape)
plt.title('Homemade IDW')

