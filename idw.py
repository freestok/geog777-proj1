""" invdisttree.py: inverse-distance-weighted interpolation using KDTree
    fast, solid, local
"""
from __future__ import division
import numpy as np
from scipy.spatial import cKDTree as KDTree
    # http://docs.scipy.org/doc/scipy/reference/spatial.html
from os.path import join
import geopandas as gpd
import matplotlib.pyplot as plt
import gdal
import osr
from osgeo import gdal_array
__date__ = "2010-11-09 Nov"  # weights, doc

#...............................................................................
class Invdisttree:
    """ inverse-distance-weighted interpolation using KDTree:
invdisttree = Invdisttree( X, z )  -- data points, values
interpol = invdisttree( q, nnear=3, eps=0, p=1, weights=None, stat=0 )
    interpolates z from the 3 points nearest each query point q;
    For example, interpol[ a query point q ]
    finds the 3 data points nearest q, at distances d1 d2 d3
    and returns the IDW average of the values z1 z2 z3
        (z1/d1 + z2/d2 + z3/d3)
        / (1/d1 + 1/d2 + 1/d3)
        = .55 z1 + .27 z2 + .18 z3  for distances 1 2 3

    q may be one point, or a batch of points.
    eps: approximate nearest, dist <= (1 + eps) * true nearest
    p: use 1 / distance**p
    weights: optional multipliers for 1 / distance**p, of the same shape as q
    stat: accumulate wsum, wn for average weights

How many nearest neighbors should one take ?
a) start with 8 11 14 .. 28 in 2d 3d 4d .. 10d; see Wendel's formula
b) make 3 runs with nnear= e.g. 6 8 10, and look at the results --
    |interpol 6 - interpol 8| etc., or |f - interpol*| if you have f(q).
    I find that runtimes don't increase much at all with nnear -- ymmv.

p=1, p=2 ?
    p=2 weights nearer points more, farther points less.
    In 2d, the circles around query points have areas ~ distance**2,
    so p=2 is inverse-area weighting. For example,
        (z1/area1 + z2/area2 + z3/area3)
        / (1/area1 + 1/area2 + 1/area3)
        = .74 z1 + .18 z2 + .08 z3  for distances 1 2 3
    Similarly, in 3d, p=3 is inverse-volume weighting.

Scaling:
    if different X coordinates measure different things, Euclidean distance
    can be way off.  For example, if X0 is in the range 0 to 1
    but X1 0 to 1000, the X1 distances will swamp X0;
    rescale the data, i.e. make X0.std() ~= X1.std() .

A nice property of IDW is that it's scale-free around query points:
if I have values z1 z2 z3 from 3 points at distances d1 d2 d3,
the IDW average
    (z1/d1 + z2/d2 + z3/d3)
    / (1/d1 + 1/d2 + 1/d3)
is the same for distances 1 2 3, or 10 20 30 -- only the ratios matter.
In contrast, the commonly-used Gaussian kernel exp( - (distance/h)**2 )
is exceedingly sensitive to distance and to h.

    """
# anykernel( dj / av dj ) is also scale-free
# error analysis, |f(x) - idw(x)| ? todo: regular grid, nnear ndim+1, 2*ndim

    def __init__( self, X, z, leafsize=10, stat=0 ):
        assert len(X) == len(z), "len(X) %d != len(z) %d" % (len(X), len(z))
        self.tree = KDTree( X, leafsize=leafsize )  # build the tree
        self.z = z
        self.stat = stat
        self.wn = 0
        self.wsum = None

    def __call__( self, q, nnear=6, eps=0, p=1, weights=None ):
            # nnear nearest neighbours of each query point --
        q = np.asarray(q)
        qdim = q.ndim
        if qdim == 1:
            q = np.array([q])
        if self.wsum is None:
            self.wsum = np.zeros(nnear)

        self.distances, self.ix = self.tree.query( q, k=nnear, eps=eps )
        interpol = np.zeros( (len(self.distances),) + np.shape(self.z[0]) )
        jinterpol = 0
        for dist, ix in zip( self.distances, self.ix ):
            if nnear == 1:
                wz = self.z[ix]
            elif dist[0] < 1e-10:
                wz = self.z[ix[0]]
            else:  # weight z s by 1/dist --
                w = 1 / dist**p
                if weights is not None:
                    w *= weights[ix]  # >= 0
                w /= np.sum(w)
                wz = np.dot( w, self.z[ix] )
                if self.stat:
                    self.wn += 1
                    self.wsum += w
            interpol[jinterpol] = wz
            jinterpol += 1
        return interpol if qdim > 1  else interpol[0]

def generate_query_points(gdf: gpd.GeoDataFrame, resolution: int):
    xmin, ymin, xmax, ymax = gdf.total_bounds
    # nx = ny = resolution
    xi = np.linspace(xmin, xmax, resolution)
    yi = np.linspace(ymin, ymax, resolution)
    # points = []
    # xi, yi = np.meshgrid(xi, yi)
    # xi, yi = np.meshgrid(xi,yi)
    g = np.meshgrid(xi, yi)
    flattened = zip(*(i.flat for i in g))
    return np.asarray(list(flattened))
    # points = set(t)
    # points = []
    # for i in range(resolution):
    #     row = []
    #     for j in range(resolution):
    #         point = yi[i, j], xi[i, j]
    #         #    point = Point(point)
    #         row.append(point)
    #     points.append(row)
    # return np.asarray(points)

def plot(x,y,z,grid):
    plt.figure()
    plt.imshow(grid, extent=(x.min(), x.max(), y.max(), y.min()))
    # plt.hold(True)
    # plt.scatter(x,y,c=z)
    plt.colorbar()
#...............................................................................
if __name__ == "__main__":
    import sys
    gdf = gpd.read_file(join('shapefiles', 'well_nitrate.shp'))

    N = gdf.shape[0] # number of rows
    Ndim = 2 # for mimicing points
    Nask = N  # N Nask 1e5: 24 sec 2d, 27 sec 3d on mac g4 ppc
    Nnear = 11  # 8 2d, 11 3d => 5 % chance one-sided -- Wendel, mathoverflow.com
    leafsize = 10
    eps = .1  # approximate nearest, dist <= (1 + eps) * true nearest
    p = 2  # weights ~ 1 / distance**p
    # cycle = .25
    seed = 1
    resolution = 500

    # exec "\n".join( sys.argv[1:] )  # python this.py N= ...
    np.random.seed(seed )
    np.set_printoptions( 3, threshold=100, suppress=True )  # .3f

    # print("\nInvdisttree:  N %d  Ndim %d  Nask %d  Nnear %d  leafsize %d  eps %.2g  p %.2g" % (
    #     N, Ndim, Nask, Nnear, leafsize, eps, p))

    # def terrain(x):
    #     """ ~ rolling hills """
    #     return np.sin( (2*np.pi / cycle) * np.mean( x, axis=-1 ))

    gdf['x'] = gdf.geometry.x
    gdf['y'] = gdf.geometry.y
    xy = gdf[['x', 'y']].to_numpy()
    z1 = gdf.nitr_ran.to_numpy()
    ask = generate_query_points(gdf, resolution)
    ask_reshape = ask.reshape((resolution, resolution, 2))
    # known = np.random.uniform( size=(N,Ndim) ) ** .5  # 1/(p+1): density x^p
    # z = terrain( known )

    ## random points to get?
    # xmin, ymin, xmax, ymax = gdf.total_bounds
    # ask = np.random.uniform( size=(Nask,Ndim) )

#...............................................................................
    invdisttree = Invdisttree( xy, z1, leafsize=leafsize, stat=1 )
    interpol = invdisttree(ask, nnear=Nnear, eps=eps, p=p )
    final_grid = interpol.reshape((resolution, resolution))

    # final_grid = np.transpose(final_grid)
    final_grid = np.flip(final_grid, 0)
    plot(gdf.x, gdf.y, 0, final_grid)


    
    # xmin,ymin,xmax,ymax = [lon.min(),lat.min(),lon.max(),lat.max()]
    nrows,ncols = np.shape(final_grid)
    xmin, ymin, xmax, ymax = gdf.total_bounds

    xres = (xmax-xmin)/float(ncols)
    yres = (ymax-ymin)/float(nrows)
    geotransform=(xmin,xres,0,ymax,0, -yres)   
    # That's (top left x, w-e pixel resolution, rotation (0 if North is up), 
    #         top left y, rotation (0 if North is up), n-s pixel resolution)
    # I don't know why rotation is in twice???

    file_name = 'myraster_flip.tif'
    output_raster = gdal.GetDriverByName('GTiff').Create(file_name,ncols, nrows, 1 ,gdal.GDT_Float32)  # Open the file
    output_raster.SetGeoTransform(geotransform)  # Specify its coordinates
    srs = osr.SpatialReference()                 # Establish its coordinate encoding
    srs.ImportFromEPSG(4326)                     # This one specifies WGS84 lat long.
                                                # Anyone know how to specify the 
                                                # IAU2000:49900 Mars encoding?
    output_raster.SetProjection( srs.ExportToWkt() )   # Exports the coordinate system 
                                                    # to the file
    output_raster.GetRasterBand(1).WriteArray(final_grid)   # Writes my array to the raster

    output_raster.FlushCache()

    dataset = gdal.Open(file_name)
    band1 = dataset.GetRasterBand(1) # Red channel

    # b1 = np.transpose(band1.ReadAsArray())
    # b1 = np.flip(b1, 0)
    b1_og = band1.ReadAsArray()
    # print(b1)
    # img = np.dstack((b1, b2, b3))
    # f = plt.figure()
    plt.imshow(b1_og)
    dataset=None
    # print( "average distances to nearest points: %s" % \
    #     np.mean( invdisttree.distances, axis=0 ))
    # print("average weights: %s" % (invdisttree.wsum / invdisttree.wn))
        # see Wikipedia Zipf's law
    # err = np.abs( terrain(ask) - interpol )
    # print("average |terrain() - interpolated|: %.2g" % np.mean(err))

    # print "interpolate a single point: %.2g" % \
    #     invdisttree( known[0], nnear=Nnear, eps=eps )
