"""
Some utility functions for winGliknk format.
"""
import numpy as np

from osgeo import osr

def epsg_to_crs(epsg):
    """
    A CRS from its EPSG code.
    """
    source_ref = osr.SpatialReference()
    source_ref.ImportFromEPSG(epsg)
    return source_ref

def get_utm_zone(latitude, longitude):
    """
    Returns zone number and hemisphere if lat, lon is provided.
    """
    result = osr.SpatialReference()
    result.SetWellKnownGeogCS('WGS84')
    is_northern = int(latitude >= 0)
    zone_number = int(1 + (longitude + 180.0) / 6.0)
    result.SetUTM(zone_number, is_northern)
    return result

def crs_to_str(crs):
    """
::    Retrieve the crs as string.
    """
    return crs.GetAttrValue('GEOGCS') or crs.GetAttrValue('PROJCS')

wgs84_crs = epsg_to_crs(4326)

def transform_coords(coordinates, from_crs, to_crs):
    """
    A GDAL functionality to transform coordinates from one spatial ref to another.
    """
    if from_crs.IsSame(to_crs):
        return coordinates

    return osr.CoordinateTransformation(from_crs, to_crs).TransformPoints(coordinates)

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
    return Interval(left=np.min(arr), right=np.max(arr))



class Points3D:
    """
    A collection of 3D points in (x, y, z) coordinates.
    """
    def __init__(self, points, crs):
        assert len(points.shape) == 2
        self.points = points
        self.crs = crs
        assert self.dim == 3

    @property
    def count(self):
        return self.points.shape[0]

    @property
    def dim(self):
        return self.points.shape[1]

    @property
    def x(self):
        return self.points[:, 0]

    @property
    def y(self):
        return self.points[:, 1]

    @property
    def z(self):
        return self.points[:, 2]

    def fix_false_origin(self):
        false_easting = self.crs.GetProjParm(osr.SRS_PP_FALSE_EASTING)
        false_northing = self.crs.GetProjParm(osr.SRS_PP_FALSE_NORTHING)
        return Points3D(np.array([self.x - false_easting, self.y - false_northing, self.z]).T, self.crs)

    def clipping_mask(self, grid):
        assert self.crs.IsSame(grid.crs)

        xbounds, ybounds, zbounds = grid.bounds

        and_ = np.logical_and

        return and_(xbounds.clipping_mask(self.x),
                    and_(ybounds.clipping_mask(self.y),
                         zbounds.clipping_mask(self.z)))

    def __repr__(self):
        return "<collection of {} points in {}>".format(self.count, crs_to_str(self.crs))

    def to_crs(self, crs):
        xy = np.array([self.x, self.y]).T
        new_xy = np.array(transform_coords(xy, self.crs, crs))[:, 0:2]
        new_xyz = np.array([new_xy[:, 0], new_xy[:, 1], self.z]).T
        return Points3D(new_xyz, crs)

    def to_gridded_array(self, shape):
        return self.points.reshape(shape)


class Grid3D:
    def __init__(self, x, y, z, crs):
        assert len(x.shape) == 1
        assert len(y.shape) == 1
        assert len(z.shape) == 1

        self.crs = crs
        self.x = x
        self.y = y
        self.z = z
        self.grid = np.stack(np.meshgrid(x, y, z), axis=-1)

    @staticmethod
    def from_extent_and_resolution(left, right, resolution, crs):
        # coords are for cell centers
        x, y, z = [np.arange(left[dim] + resolution[dim] / 2., right[dim], resolution[dim])
                   for dim in range(3)]
        return Grid3D(x, y, z, crs)

    def __repr__(self):
        return "<grid of shape {} in {}>".format(self.grid.shape, crs_to_str(self.crs))

    @property
    def shape(self):
        return self.grid.shape

    @property
    def bounds(self):
        return (bounds(self.x), bounds(self.y), bounds(self.z))

    def flatten(self):
        xdim, ydim, zdim, vec_dim = self.grid.shape
        return Points3D(self.grid.reshape((xdim * ydim * zdim, vec_dim)), self.crs)
