from os.path import join

import geopandas as gpd
import numpy as np

from idw import IDWGenerator

#...............................................................................
NNEAR = 8  # number of nearest neighbors to look for
EPS = .1  # approximate nearest, dist <= (1 + eps) * true nearest
P = 2  # weights ~ 1 / distance**p
RESOLUTION = 500

gdf = gpd.read_file(join('shapefiles', 'well_nitrate.shp'))
gdf['x'] = gdf.geometry.x
gdf['y'] = gdf.geometry.y
xy = gdf[['x', 'y']].to_numpy()
z = gdf.nitr_ran.to_numpy()
# ask = generate_query_points(gdf, RESOLUTION)


idw_gen = IDWGenerator(xy, z, gdf.total_bounds)
ask = idw_gen.generate_query_points(resolution=RESOLUTION)
interpol = idw_gen.interpolate(ask, nnear=NNEAR, eps=EPS, p=P )

## get it into a grid, flip it (for some reason it is upside-down)
final_grid = interpol.reshape((RESOLUTION, RESOLUTION))
# final_grid = np.transpose(final_grid)
final_grid = np.flip(final_grid, 0)

# export to tif
raster_file = 'myraster_flip.tif'
idw_gen.export_to_tif(final_grid, raster_file)

# calc zonal stats by reading in tif
# possible it could just read in an array instead of first exporting to 
# tif first
vector_file = join('shapefiles', 'cancer_tracts_wgs84.shp')
output_file = join('shapefiles', 'aggregateNitrate.shp')
idw_gen.calculate_zonal_stats(vector_file, raster_file, output_file)