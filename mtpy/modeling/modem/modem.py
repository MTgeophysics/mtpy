#!/usr/bin/env python

import os
from os.path import join, abspath

import numpy as np
from osgeo import osr
from scipy.interpolate import interp2d
from pyproj import Proj, transform

from mtpy.modeling.modem import Model, Data, write_resistivity_grid, wgs84_crs, get_utm_zone


def mid_point(arr):
    """
    Mid point of a one-dimensional array `arr`.
    Shape of `arr` should be just a tuple with one value in it.
    """
    [shape] = arr.shape
    mid = int(shape / 2)
    return arr[mid]


def uniform_interior_grid(arr, spacing, mid):
    """
    Make a regular grid in lat/lon space.
    """
    end_points = sorted((arr[0], arr[-1]))
    units = int((end_points[0] - mid) / spacing), int((end_points[1] - mid) / spacing) + 1

    candidates = [mid + i * spacing for i in range(units[0], units[1])]
    return np.array([x for x in candidates if end_points[0] < x and x < end_points[1]])


def median_spacing(arr):
    """
    The spacing that occurs the maximum time.
    """
    return np.median(arr[1:] - arr[:-1])


def lon_lat_grid_spacing(center, width, height, to_wgs84):
    """
    Returns center longitude and latitude, and spacing for longitude and latitude.
    """
    center_x, center_y = center.east.item(), center.north.item()
    center_lon, center_lat = to_wgs84(center_x, center_y)
    shifted_lon, _ = to_wgs84(center_x + width, center_y)
    _, shifted_lat = to_wgs84(center_x, center_y + height)

    return center_lon, center_lat, shifted_lon - center_lon, shifted_lat - center_lat


def interpolated_layer(x, y, layer):
    """
    Create the interpolation function for each layer.
    """
    return interp2d(x, y, layer)  # bounds_error=True


# TODO can this be done with lib.Points3D?
def converter(in_spatial_ref, out_spatial_ref):
    """
    Transfrom coordinates from one spatial ref to another.
    """
    in_proj = Proj(in_spatial_ref.ExportToProj4())
    out_proj = Proj(out_spatial_ref.ExportToProj4())

    def result(x, y):
        return transform(in_proj, out_proj, x, y)

    return result


def main():
    # Define Data and Model Paths
    MT_PATH = abspath(join(__file__, '../../../..'))

    data = Data()
    data.read_data_file(data_fn=join(MT_PATH, 'examples/model_files/ModEM/ModEM_Data.dat'))

    # create a model object using the data object and read in model data
    model = Model(data_obj=data)
    model.read_model_file(model_fn=join(MT_PATH, 'examples/model_files/ModEM/ModEM_Model_File.rho'))

    center = data.center_point

    resistivity_data = {
        'x': center.east.item() + (model.grid_east[1:] + model.grid_east[:-1])/2,
        'y': center.north.item() + (model.grid_north[1:] + model.grid_north[:-1])/2,
        'z': (model.grid_z[1:] + model.grid_z[:-1])/2,
        'resistivity': np.transpose(model.res_model, axes=(2, 0, 1))
    }

    spatial_ref = get_utm_zone(center.lat.item(), center.lon.item())

    to_wgs84 = converter(spatial_ref, wgs84_crs)
    from_wgs84 = converter(wgs84_crs, spatial_ref)

    center_lon, center_lat, width, height = lon_lat_grid_spacing(center, median_spacing(model.grid_east),
                                                                 median_spacing(model.grid_north), to_wgs84)

    # TODO use lib.Grid3D?
    lon_list = [to_wgs84(x, y)[0]
                for x in resistivity_data['x']
                for y in resistivity_data['y']]

    lat_list = [to_wgs84(x, y)[1]
                for x in resistivity_data['x']
                for y in resistivity_data['y']]

    interpolation_funcs = [interpolated_layer(resistivity_data['x'],
                                              resistivity_data['y'],
                                              resistivity_data['resistivity'][z_index, :, :])
                           for z_index in range(resistivity_data['z'].shape[0])]

    result = {
        'longitude': uniform_interior_grid(sorted(lon_list), width, center_lon),
        'latitude': uniform_interior_grid(sorted(lat_list), height, center_lat),
        'depth': resistivity_data['z']
    }

    result['resistivity'] = np.zeros(tuple(result[key].shape[0]
                                           for key in ['depth', 'latitude', 'longitude']))

    def uniform_layer(interp_func, latitudes, longitudes):
        """
        Calculate the interpolated values for the layer.
        """
        lats, lons = latitudes.shape[0], longitudes.shape[0]

        result = np.zeros((lats, lons))
        for j in range(lons):
            for i in range(lats):
                lon, lat = longitudes[j], latitudes[i]
                x, y = from_wgs84(lon, lat)

                result[i, j] = interp_func(x, y)

        return result

    for z_index in range(result['depth'].shape[0]):
        print 'layer #', z_index + 1
        result['resistivity'][z_index, :, :] = uniform_layer(interpolation_funcs[z_index],
                                                             result['latitude'], result['longitude'])

    write_resistivity_grid('wgs84.nc', wgs84_crs,
                           result['latitude'], result['longitude'], result['depth'],
                           result['resistivity'], z_label='depth')


if __name__ == '__main__':
    main()
