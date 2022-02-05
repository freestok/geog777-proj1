import geopandas as gpd
from os.path import join
import idw
import matplotlib.pyplot as plt

def plot(x,y,z,grid):
    plt.figure()
    plt.imshow(grid, extent=(x.min(), x.max(), y.max(), y.min()))
    # plt.hold(True)
    # plt.scatter(x,y,c=z)
    plt.colorbar()

print('read file...')
gdf = gpd.read_file(join('shapefiles', 'well_nitrate.shp'))
x = gdf.geometry.x.to_numpy()
y = gdf.geometry.y.to_numpy()
z = gdf.nitr_ran.to_numpy()

print('Start idw.main')
idw_results = idw.main(x, y, z, p=1.5, n=50)
grid1 = idw_results['z_idw']
print('done!')

# for plotting
plot(x, y, z, grid1)