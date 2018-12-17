#!/usr/bin/env python
"""
# created by U.Zannat from GA, Nov, 2018
+++++++++++++++++
This code works on MTpy modEM format output and generate a netCDF file.
The modEM format generates an irregular grid in local unit(east, north, depth).
Here we converted that grid in global crs (in lon, lat, depth) and in regular
grid. The regular spacing we used here is the mode of our input data. The
interpolation is done in local grid to have the advantage of having the same
unit for our three variables (east, north, depth in km).
+++++++++++++++++
"""
import os
from os.path import join, abspath

import numpy as np
from scipy import stats
from scipy.interpolate import interp2d
from pyproj import Proj, transform

from mtpy.modeling.modem import Model, Data
from mtpy.utils import lib
from mtpy.utils import nc
from mtpy.utils import gis_tools


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
    return stats.mode(arr[1:] - arr[:-1]).mode.item()

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


def converter(in_proj, out_proj):
    """
    Transfrom coordinates from one epsg to another.
    """
    def result(x, y):
        return transform(in_proj, out_proj, x, y)

    return result


def main():
    # Define Data and Model Paths
    MT_PATH = abspath(join(__file__, '../../../'))

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

    zone_number, is_northern, utm_zone = gis_tools.get_utm_zone(center.lat.item(), center.lon.item())

    wgs84_proj = Proj(init='epsg:4326')
    source_proj = Proj('+proj=utm +zone=%d +%s +datum=%s' % (zone_number, 'north' if is_northern else 'south', 'WGS84'))
    to_wgs84 = converter(source_proj, wgs84_proj)
    from_wgs84 = converter(wgs84_proj, source_proj)

    center_lon, center_lat, width, height = lon_lat_grid_spacing(center, median_spacing(model.grid_east),
                                                                 median_spacing(model.grid_north), to_wgs84)

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

    nc.write_resistivity_grid('wgs84.nc', 4326,
                              result['latitude'], result['longitude'], result['depth'],
                              result['resistivity'], z_label='depth')


if __name__ == '__main__':
    main()
