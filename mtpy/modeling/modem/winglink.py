import lib
import interp
import nc
import numpy as np

source_crs = lib.epsg_to_crs(32753)

geo_grid = lib.Grid3D.from_extent_and_resolution(left=(131.5, -31., -4500.),
                                                 right=(132.5, -30., -125.),
                                                 resolution=(0.005, 0.005, 100.),
                                                 crs=lib.wgs84_crs)

utm_points = geo_grid.flatten().to_crs(source_crs)

utm_grid = utm_points.to_gridded_array(geo_grid.shape)

arr = np.loadtxt('04_LP_3dmod_it100.xyzv')

source_points = lib.Points3D(arr[:, 0:3], source_crs).fix_false_origin()

source_val = arr[:, 3]

source_in_geo = source_points.to_crs(lib.wgs84_crs)

mask = source_in_geo.clipping_mask(geo_grid)

source_points.points = source_points.points[mask]
source_val = source_val[mask]

# resistivity = interp.nearest(source_points, source_val, utm_points)
resistivity = interp.linear(source_points, source_val, utm_grid)
# resistivity = resistivity.reshape(geo_grid.shape[:3])
print(resistivity.shape)

nc.write_resistivity_grid('nearest.nc', lib.WGS84_crs, geo_grid.y, geo_grid.x, geo_grid.z, resistivity.transpose([2, 0, 1]))
