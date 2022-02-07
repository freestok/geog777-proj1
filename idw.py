from os.path import join

import gdal
import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
import osr
import pandas as pd
import rasterio
from rasterstats import zonal_stats
from scipy.spatial import KDTree

#...............................................................................
class IDWGenerator:
    """ 
    https://stackoverflow.com/a/3119544/8947008
    inverse-distance-weighted interpolation using KDTree:
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
    def __init__( self, X, z, bounds, leafsize=10, stat=1 ):
        assert len(X) == len(z), "len(X) %d != len(z) %d" % (len(X), len(z))
        self.tree = KDTree( X, leafsize=leafsize )  # build the tree
        self.z = z
        self.stat = stat
        self.wn = 0
        self.wsum = None
        self.bounds = bounds

    def interpolate( self, q, nnear=6, eps=0, p=1, weights=None ):
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
        if nnear == 1:
            for dist, ix in zip( self.distances, self.ix ):
                wz = self.z[ix]
            interpol[jinterpol] = wz
            jinterpol += 1
        else:
            for dist, ix in zip( self.distances, self.ix ):
                if dist[0] < 1e-10:
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
    
    def plot(self, tif_file):
        # plot
        dataset = gdal.Open(tif_file)
        band1 = dataset.GetRasterBand(1)
        plt.imshow(band1.ReadAsArray())
        dataset=None

    def generate_query_points(self, resolution: int):
        # assumes extent is a perfect square, but that isn't necessarily true
        xmin, ymin, xmax, ymax = self.bounds
        # nx = ny = resolution
        xi = np.linspace(xmin, xmax, resolution)
        yi = np.linspace(ymin, ymax, resolution)
        mesh = np.meshgrid(xi, yi)
        flattened = zip(*(i.flat for i in mesh))
        return np.asarray(list(flattened))

    def export_to_tif(self, grid, output_tif):
        # export to tif
        # https://gis.stackexchange.com/a/37431/78614
        nrows,ncols = np.shape(grid)
        xmin, ymin, xmax, ymax = self.bounds

        xres = (xmax-xmin)/float(ncols)
        yres = (ymax-ymin)/float(nrows)
        geotransform=(xmin,xres,0,ymax,0, -yres)

        output_raster = gdal.GetDriverByName('GTiff').Create(output_tif,ncols, nrows, 1 ,gdal.GDT_Float32)  # Open the file
        output_raster.SetGeoTransform(geotransform)  # Specify its coordinates
        srs = osr.SpatialReference()                 # Establish its coordinate encoding
        srs.ImportFromEPSG(4326)                     # This one specifies WGS84 lat long.
        output_raster.SetProjection( srs.ExportToWkt() )   # Exports the coordinate system 
        output_raster.GetRasterBand(1).WriteArray(grid)   # Writes my array to the raster
        output_raster.FlushCache()
        output_raster = None

    def calculate_zonal_stats(self, vector_file, raster_file, output_file):
        tracts = gpd.read_file(vector_file)
        stats = zonal_stats(vector_file, raster_file)
        df = pd.DataFrame(stats)
        final_df = pd.merge(tracts, df, left_index=True, right_index=True)
        final_df.to_file(output_file)