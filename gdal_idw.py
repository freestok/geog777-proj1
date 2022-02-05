import pandas as pd
import json
import math
import os
from osgeo import gdal
from shapely.geometry import Point
import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
from rasterstats import zonal_stats
from rasterio.plot import show
import rasterio
from zonal_stats import zonal_stats
shp = os.path.join('shapefiles', 'nitrate_wgs84.shp')
gdf = gpd.read_file(shp)
# gdf.geometry = gdf.apply(lambda r: Point(r.geometry.x, r.geometry.y, r.nitr_ran), axis=1)
# new_shp = os.path.join('shapefiles', 'well_nitrate2.shp')
# gdf.to_file(new_shp)
# points = ogr.Open(os.path.join('shapefiles', 'well_nitrate.shp'), 0)
# layer = points.GetLayer()

xmin, ymin, xmax, ymax = gdf.total_bounds
dem = gdal.Open('tifs//points_rasterized_bigger.tif')
gt = dem.GetGeoTransform()
ulx = gt[0]
uly = gt[3]
res = gt[1]
xsize = dem.RasterXSize
ysize = dem.RasterYSize
lrx = ulx + xsize * res
lry = uly - ysize * res
dem = None

# grid settings
# pixel = 100

w1 = xsize
h1 = ysize
w2 = math.floor(w1 * 1)
h2 = math.floor(h1 * 1)
p = 2
radius1 = 2000
radius2 = 2000
min_points = 0
smoothing = 0.05
file_name = f'tifs//idw-s{smoothing}-p{p}.tif'
ds = gdal.Grid(file_name, shp, format='GTiff', zfield='nitr_ran',
            outputBounds=[ulx, uly, lrx, lry], width=w2, height=h2,
            algorithm=f'invdist:power={p}:smoothing={smoothing}:radius1={radius1}:radius2={radius2}:min_points={min_points}')
ds = None

# ---------------------------- just plotting stuff -----------------------------
dataset = gdal.Open(file_name)
band1 = dataset.GetRasterBand(1) # Red channel

# b1 = np.transpose(band1.ReadAsArray())
# b1 = np.flip(b1, 0)
b1_og = band1.ReadAsArray()
# print(b1)
# img = np.dstack((b1, b2, b3))
f = plt.figure()
plt.imshow(b1_og)
# plt.savefig(f'idw{int(smoothing)*100}.png', transparent=False)
plt.show()
# gdf.plot(column='nitr_ran')
dataset = band1 = None

# x = gdf.geometry.x
# y = gdf.geometry.y
# z = gdf.nitr_ran

# plt.scatter(x, y, c=z)
# ----------------------------- plot end ---------------------------------------

## zonal stats
# tracts = gpd.read_file('shapefiles//cancer_tracts_wgs84.shp')
# stats = zonal_stats('shapefiles//cancer_tracts_wgs84.shp', file_name)
# df = pd.DataFrame(stats)
# final_df = pd.merge(tracts, df, left_index=True, right_index=True)
# final_df.to_file('shapefiles//aggregateNitrate.shp')

print('complete!')