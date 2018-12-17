#!/usr/bin/env python

"""
Created by U.Zannat from GA, Nov, 2018
+++++++++++++++++
This code works on MTpy winGlink format output and generate a netCDF file.
In WinGlink format the data are in points and in utm (meter). An input
file has been added named 'winGlink_3dmod_input.xyzv' where the col1 is lat,
col2 is lon, col3 is depth and col4 is resistivity. Here we interpolated the
irregular points into a grid in global crs using the inverse distance weighing as the interpolator function.
+++++++++++++++++
"""
from mtpy.utils import interp, nc, lib
import numpy as np

def main():
    source_epsg = 32753
    wgs84_epsg = 4326

    geo_grid = lib.Grid3D.from_extent_and_resolution(left=(131.5, -31., -4500.),
                                                     right=(132.5, -30., -125.),
                                                     resolution=(0.005, 0.005, 100.),
                                                     epsg=wgs84_epsg)

    utm_points = geo_grid.flatten().to_epsg(source_epsg)

    arr = np.loadtxt('winGlink_3dmod_input.xyzv')

    false_easting = 500000.0
    false_northing = 10000000.0

    source_points = lib.Points3D(np.array([arr[:, 0] - false_easting, arr[:, 1] - false_northing, arr[:, 2]]).T, source_epsg)
    source_val = arr[:, 3]

    source_in_geo = source_points.to_epsg(wgs84_epsg)

    mask = source_in_geo.clipping_mask(geo_grid)

    source_points.points = source_points.points[mask]
    source_val = source_val[mask]

    # resistivity = interp.nearest(source_points, source_val, utm_points)
    resistivity = interp.IDW(source_points, source_val, utm_points)
    resistivity = resistivity.reshape(geo_grid.shape[:3])

    nc.write_resistivity_grid('IDW.nc', wgs84_epsg, geo_grid.y, geo_grid.x, geo_grid.z, resistivity.transpose([2, 0, 1]))

if __name__ == "__main__":
    main()
