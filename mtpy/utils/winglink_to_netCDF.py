#!/usr/bin/env python

"""
This code works on MTpy winGlink format output and generate a netCDF file.
In WinGlink format the data are in points and in utm (meter). An input
file has extension '.xyzv' where the col1 is lat,
col2 is lon, col3 is depth and col4 is resistivity. Here we interpolated the
irregular points into a grid in global CRS using the inverse distance weighing as the interpolator function.

@author Umma Zannat
"""

import argparse

import numpy as np
from pyproj import Proj

from mtpy.utils import nc


def main(input_file, output_file):
    source_proj = Proj(init='epsg:32753')
    wgs84_proj = Proj(init='epsg:4326')

    (geo_x, geo_y, geo_z, geo_grid) = nc.grid_from_extent_and_resolution(left=(131.5, -31., -4500.),
                                                                         right=(132.5, -30., -125.),
                                                                         resolution=(0.005, 0.005, 100.))

    utm_points = nc.transform_3d(wgs84_proj, source_proj, nc.flatten_grid(geo_grid))

    arr = np.loadtxt(input_file)

    false_easting = 500000.0
    false_northing = 10000000.0

    source_points = np.array([arr[:, 0] - false_easting, arr[:, 1] - false_northing, arr[:, 2]]).T
    source_val = arr[:, 3]

    source_in_geo = nc.transform_3d(source_proj, wgs84_proj, source_points)

    mask = nc.clipping_mask(source_in_geo, geo_grid)

    source_points = source_points[mask]
    source_val = source_val[mask]

    resistivity = nc.IDW(source_points, source_val, utm_points)
    resistivity = resistivity.reshape(geo_grid.shape[:3])

    nc.write_resistivity_grid(output_file, wgs84_proj, geo_y, geo_x, geo_z, resistivity.transpose([2, 0, 1]))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('input_file', help='input .xyzv file')
    parser.add_argument('output_file', help='output .nc file')
    args = parser.parse_args()
    main(args.input_file, args.output_file)
