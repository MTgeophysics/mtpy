"""
Some utility functions for winGlink format.
"""
import numpy as np
from pyproj import transform, Proj
from mtpy.utils import gis_tools


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
    def __init__(self, points, epsg):
        assert len(points.shape) == 2
        self.points = points
        self.epsg = epsg
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

    def clipping_mask(self, grid):
        assert self.epsg == grid.epsg

        xbounds, ybounds, zbounds = grid.bounds

        and_ = np.logical_and

        return and_(xbounds.clipping_mask(self.x),
                    and_(ybounds.clipping_mask(self.y),
                         zbounds.clipping_mask(self.z)))

    def __repr__(self):
        return "<collection of {} points in epsg:{}>".format(self.count, self.epsg)

    def to_epsg(self, epsg):
        proj_from = Proj(gis_tools.EPSG_DICT[self.epsg])
        proj_to = Proj(gis_tools.EPSG_DICT[epsg])
        new_x, new_y = transform(proj_from, proj_to, self.x, self.y)
        new_xyz = np.array([new_x, new_y, self.z]).T
        return Points3D(new_xyz, epsg)

    def to_gridded_array(self, shape):
        return self.points.reshape(shape)


class Grid3D:
    def __init__(self, x, y, z, epsg):
        assert len(x.shape) == 1
        assert len(y.shape) == 1
        assert len(z.shape) == 1

        self.epsg = epsg
        self.x = x
        self.y = y
        self.z = z
        self.grid = np.stack(np.meshgrid(x, y, z), axis=-1)

    @staticmethod
    def from_extent_and_resolution(left, right, resolution, epsg):
        # coords are for cell centers
        x, y, z = [np.arange(left[dim] + resolution[dim] / 2., right[dim], resolution[dim])
                   for dim in range(3)]
        return Grid3D(x, y, z, epsg)

    def __repr__(self):
        return "<grid of shape {} in epsg:{}>".format(self.grid.shape, self.epsg)

    @property
    def shape(self):
        return self.grid.shape

    @property
    def bounds(self):
        return (bounds(self.x), bounds(self.y), bounds(self.z))

    def flatten(self):
        xdim, ydim, zdim, vec_dim = self.grid.shape
        return Points3D(self.grid.reshape((xdim * ydim * zdim, vec_dim)), self.epsg)
