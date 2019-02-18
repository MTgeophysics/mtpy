#!/usr/bin/env python

"""
This code operates on a MTpy WinGlink format file and generates a NetCDF file.
In the WinGlink format the first three columns are x, y, z coordinates in UTM (meter),
and the fourth is the resistivity. Here we interpolated the irregularly-spaced data
points into a grid in global WGS84 CRS using the Inverse Distance Weighing as the
interpolation method.

@author Umma Zannat, GA, 2019
"""

import argparse

import numpy as np
from pyproj import Proj

from . import nc


def read_winglink_xyzv(input_file, false_easting=0.0, false_northing=0.0):
    """
    Read points and values in WinGlink format.
    Optionally fix false easting and northing.
    """
    arr = np.loadtxt(input_file)
    source_points = np.array([arr[:, 0] - false_easting, arr[:, 1] - false_northing, arr[:, 2]]).T
    source_values = arr[:, 3]
    return source_points, source_values


def mask_interpolate_and_write(output_file,
                               source_points, source_values, source_proj,
                               grid_spec, grid_proj):
    """
    Interpolate point data to grid and write to netCDF file.
    """
    grid_x, grid_y, grid_z, grid_points = grid_spec

    # flatten grid because IDW takes a flattened array as input
    utm_points = nc.transform_3d(grid_proj, source_proj, nc.flatten_grid(grid_points))

    # mask source for grid
    # this is needed because the high resistivity of the air layers contaminate
    # the IDW interpolated data underground
    source_points_in_grid_proj = nc.transform_3d(source_proj, grid_proj, source_points)
    mask = nc.clipping_mask(source_points_in_grid_proj, grid_points)
    source_points = source_points[mask]
    source_values = source_values[mask]

    # interpolate data at grid points and reshape back to the grid
    resistivity = nc.IDW(source_points, source_values, utm_points)
    resistivity = resistivity.reshape(grid_points.shape[:3])

    # write to output NetCDF file
    nc.write_resistivity_grid(output_file, grid_proj, grid_y, grid_x, grid_z, resistivity.transpose([2, 0, 1]))


def main(input_file, source_proj, false_easting, false_northing, output_file, grid_proj, left, right, resolution):

    source_points, source_values = read_winglink_xyzv(input_file, false_easting=false_easting, false_northing=false_northing)

    grid_spec = nc.grid_from_extent_and_resolution(left=left,
                                                   right=right,
                                                   resolution=resolution)

    mask_interpolate_and_write(output_file, source_points, source_values, source_proj, grid_spec, grid_proj)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('input_file', help='input .xyzv file')
    parser.add_argument('output_file', help='output .nc file')
    parser.add_argument('--source-proj', help='EPSG projection code for input file', type=int, default=32753)
    parser.add_argument('--false-easting', help="false easting to correct", type=float, default=500000.0)
    parser.add_argument('--false-northing', help="false northing to correct", type=float, default=10000000.0)
    parser.add_argument('--grid-proj', help='EPSG projection code for output grid', type=int, default=4326)
    parser.add_argument('--grid-corner-1', help="one corner of the output grid in x,y,z format", type=str, default="131.5,-31.0,-4500.0")
    parser.add_argument('--grid-corner-2', help="the other corner of the output grid in x,y,z format", type=str, default="132.5,-30.0,-125.0")
    parser.add_argument('--grid-resolution', help='resolution of the output grid in x,y,z format', type=str, default="0.005,0.005,100.0")
    args = parser.parse_args()
    left = tuple(float(x) for x in args.grid_corner_1.split(','))
    right = tuple(float(x) for x in args.grid_corner_2.split(','))
    resolution = tuple(float(x) for x in args.grid_resolution.split(','))
    source_proj = Proj(init='epsg:' + str(args.source_proj))
    grid_proj = Proj(init='epsg:' + str(args.grid_proj))
    main(args.input_file, source_proj, args.false_easting, args.false_northing, args.output_file, grid_proj, left, right, resolution)
