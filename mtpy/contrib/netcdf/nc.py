"""
Utilities to write data to NetCDF files.

@author Umma Zannat
"""
import os

from netCDF4 import Dataset
import numpy as np
from scipy import spatial
#Not work import scipy.spatial.cKDTree as KDTree
from pyproj import transform, Proj

from mtpy.utils import gis_tools


def IDW(source_points, source_values, query_points, k=6, p=5):
    """
    Inverse Distance Weighting interpolation method with power `p`
    and `k` nearest neighbors.
    Returns the interpolated values at `query_points`.
    """
    tree = spatial.cKDTree(source_points, k)
    distances, indices = tree.query(query_points, k=k)
    inv_dist = 1. / np.power(distances, p)
    weights = inv_dist / inv_dist.sum(axis=1)[:, np.newaxis]
    return (weights * source_values[indices]).sum(axis=1)


def create_dataset(filename, overwrite=True):
    """ Create a NetCDF dataset. """
    if os.path.exists(filename):
        if overwrite:
            os.remove(filename)
        else:
            raise ValueError("file {} already exists".format(filename))

    return Dataset(filename, 'w', format='NETCDF4')


def proj_to_epsg(proj):
    for word in proj.srs.split():
        key, value = word.split('=')
        if key == "+init":
            return value

    raise ValueError('EPSG code not found for projection')


def set_grid_mapping_attrs_geographic(crs_var, proj):
    """
    Attach CRS information to a NetCDF variable.
    """
    crs_var.spatial_reference = proj_to_epsg(proj)
    crs_var.grid_mapping_name = 'latitude_longitude'
    crs_var.units = 'degree'


def write_resistivity_grid(output_file, epsg_code,
                           latitude, longitude, elevation, resistivity_data,
                           **kwargs):
    """ Resistivity_data in (elevation, latitude, longitude) grid. """

    with create_dataset(output_file) as dataset:
        dataset.description = 'Resistivity Model'

        z_label = kwargs.get('z_label', 'elevation')

        # dimensions
        dataset.createDimension('latitude', latitude.shape[0])
        dataset.createDimension('longitude', longitude.shape[0])
        dataset.createDimension(z_label, elevation.shape[0])

        # variables
        x = dataset.createVariable('latitude', 'f4', ('latitude',))
        y = dataset.createVariable('longitude', 'f4', ('longitude',))
        z = dataset.createVariable(z_label, 'f4', (z_label,))

        for var in [x, y]:
            var.units = 'degree'
        z.units = 'm'

        resistivity = dataset.createVariable('resistivity', 'f4', (z_label, 'latitude', 'longitude'))
        resistivity.grid_mapping = 'crs'
        resistivity.long_name = 'resistivity'
        resistivity.units = "ohm-m"

        # populate variables
        x[:] = latitude
        y[:] = longitude
        z[:] = elevation

        resistivity[:, :, :] = resistivity_data

        # attach crs info
        crs_var = dataset.createVariable('crs', 'i4', ())
        set_grid_mapping_attrs_geographic(crs_var, epsg_code)


class Interval:
    """
    An interval of real numbers.
    """
    def __init__(self, left, right):
        assert left <= right

        self.left = left
        self.right = right

    @property
    def width(self):
        return self.right - self.left

    @property
    def mid(self):
        return (self.left + self.right) / 2.

    def __contains__(self, other):
        if isinstance(other, Interval):
            return other.left in self and other.right in self

        return self.left < other and other < self.right

    def clipping_mask(self, array):
        return np.logical_and(self.left < array, array < self.right)

    def __repr__(self):
        return '({}, {})'.format(self.left, self.right)


def bounds(arr):
    """ The bounds (min, max) of an array as an Interval. """
    return Interval(left=np.min(arr), right=np.max(arr))


def clipping_mask(arr, grid):
    """
    Create mask for an array of points where the points fall inside the grid.
    """
    xbounds = bounds(grid[:, :, :, 0])
    ybounds = bounds(grid[:, :, :, 1])
    zbounds = bounds(grid[:, :, :, 2])

    and_ = np.logical_and

    return and_(xbounds.clipping_mask(arr[:, 0]),
                and_(ybounds.clipping_mask(arr[:, 1]),
                     zbounds.clipping_mask(arr[:, 2])))


def transform_3d(proj_from, proj_to, arr):
    """
    Transform the x, y coordinates of the points but keep the z coordinate unchanged.
    """
    x = arr[:, 0]
    y = arr[:, 1]
    z = arr[:, 2]
    new_x, new_y = transform(proj_from, proj_to, x, y)
    return np.array([new_x, new_y, z]).T

def grid_from_extent_and_resolution(left, right, resolution):
    """
    Create a grid from the coordinates of two opposite corners, and also resolution.
    The grid coordinate arrays are also returned.
    """
    # coords are for cell centers
    x, y, z = [np.arange(left[dim] + resolution[dim] / 2., right[dim], resolution[dim])
               for dim in range(3)]
    return x, y, z, np.stack(np.meshgrid(x, y, z), axis=-1)


def flatten_grid(arr):
    """
    Flatten the array of point positions, so the resultant array
    is two-dimensional, where the first dimension is the index of the point
    and the second is the vector component (0, 1, 2 for x, y, z). Useful
    for IDW input.
    """
    xdim, ydim, zdim, vec_dim = arr.shape
    return arr.reshape((xdim * ydim * zdim, vec_dim))
