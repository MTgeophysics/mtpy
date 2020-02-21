"""
Tests for resistivity model to geotiff.
"""
import pytest
import numpy as np

from mtpy.utils import convert_modem_data_to_geogrid as conv

SIRSAM_RF = 'sirsam_Na_randomforest'


@pytest.fixture
def model_grid():
    grid = [-60, -30, -20, -15, -12, -10., -9., -8., -7., -6., -5., -4., -3., -2., -1., 0.,
            1., 2., 3., 4., 5., 6., 7., 8., 9., 10., 12, 15, 20, 30, 60]
    pad = 5
    cell_size = 1.
    center = 0., 0.
    return grid, pad, cell_size, center


@pytest.fixture
def grid_centres():
    centres = [-45., -25., -17.5, -13.5, -11., -9.5, -8.5, -7.5, -6.5, -5.5, -4.5, -3.5, -2.5,
               -1.5, -0.5, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 11., 13.5,
               17.5, 25., 45.]
    stripped = [-9.5, -8.5, -7.5, -6.5, -5.5, -4.5, -3.5, -2.5, -1.5, -0.5, 0.5, 1.5, 2.5, 3.5,
                4.5, 5.5, 6.5, 7.5, 8.5, 9.5]
    stripped_keep_start = [-45., -25., -17.5, -13.5, -11., -9.5, -8.5, -7.5, -6.5, -5.5, -4.5,
                           -3.5, -2.5, -1.5, -0.5, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5,
                           9.5]
    return centres, stripped, stripped_keep_start


def test_get_centre_points(model_grid, grid_centres):
    grid, _, _, _ = model_grid
    ce = conv._get_centres(grid)

    # Should have one less centre point then there are nodes as the
    #  last node is the terminating boundary and has no centre.
    np.testing.assert_array_equal(ce, grid_centres[0])


def test_strip_padding(model_grid, grid_centres):
    _, pad, _, _ = model_grid
    centres, stripped, stripped_keep_start = grid_centres

    s = conv._strip_padding(centres, pad)
    # Stripped array should have padding removed from either side
    np.testing.assert_array_equal(s, stripped)

    s = conv._strip_padding(centres, pad, keep_start=True)
    # Stripped array while keeping start (for stripping Z padding)
    #  should only have padding removed from end
    np.testing.assert_array_equal(s, stripped_keep_start)


def test_gdal_origin(model_grid, grid_centres):
    _, _, cell_size, center_point = model_grid
    _, stripped, _ = grid_centres
    # Origin should be upper-left, so western-most point shifted 1/2
    #  cell west and northern-most point shifted 1/2 cell north
    origin = conv._get_gdal_origin(stripped, cell_size, center_point[0],
                                   stripped, cell_size, center_point[1])
    assert origin == (-10.0, 10.0)

