"""
Tests for resistivity model to geotiff.
"""
import pytest
import numpy as np

from mtpy.utils import convert_modem_data_to_geogrid as conv

EAST_GRID_PARAMS = [
    # (cell size, padding, start, end, direction)
    # start and end do not include padding cells
    (1., 5, -10., 10.),
    (10., 0, -100., 100.),
]

NORTH_GRID_PARAMS = [
    (1., 5, -10., 10.),
]

CENTER_POINTS = [
    # (x, y)
    (0., 0.),
    (128312.123, -3289321.3)
]


def make_grid_axis(cell_size, padding, start, end):
    # Build a test grid
    grid = np.arange(start, end + cell_size, cell_size)
    # Add padding. Increases exponentially. Rough approximation of
    #  the mesh_tools.get_padding_cells method
    for i in range(padding):
        pcell = (cell_size * 1.2) ** i
        start_pad = grid[0] - pcell
        end_pad = grid[-1] + pcell
        grid = np.insert(grid, 0, start_pad)
        grid = np.append(grid, end_pad)
    return grid, cell_size, padding


@pytest.fixture(params=EAST_GRID_PARAMS)
def model_grid_e(request):
    cs, p, s, e = request.param
    return make_grid_axis(cs, p, s, e)


@pytest.fixture(params=NORTH_GRID_PARAMS)
def model_grid_n(request):
    cs, p, s, e = request.param
    return make_grid_axis(cs, p, s, e)


@pytest.fixture(params=CENTER_POINTS)
def grid_center(request):
    return request.param


def test_get_center_points(model_grid_e):
    grid, cs, _ = model_grid_e
    ce = conv._get_centers(grid)
    centers = np.mean([grid[1:], grid[:-1]], axis=0)
    np.testing.assert_array_equal(ce, centers)


def test_strip_padding(model_grid_e):
    grid, _, padding = model_grid_e
    stripped = conv._strip_padding(grid, padding)
    if padding == 0:
        test_stripped = grid
    else:
        test_stripped = grid[padding:-padding]
    np.testing.assert_array_equal(stripped, test_stripped)


def test_strip_padding_keep_start(model_grid_e):
    grid, _, padding = model_grid_e
    stripped = conv._strip_padding(grid, padding, keep_start=True)
    if padding == 0:
        test_stripped = grid
    else:
        test_stripped = grid[:-padding]
    np.testing.assert_array_equal(stripped, test_stripped)


def test_gdal_origin(model_grid_e, model_grid_n, grid_center):
    grid_e, cell_size_e, _ = model_grid_e
    grid_n, cell_size_n, _ = model_grid_n
    # Origin should be upper-left, so western-most point shifted 1/2
    #  cell west and northern-most point shifted 1/2 cell north
    origin = conv._get_gdal_origin(grid_e, cell_size_e, grid_center[0],
                                   grid_n, cell_size_n, grid_center[1])
    test_origin = (grid_e[0] + grid_center[0] - cell_size_e / 2,
                   grid_n[-1] + grid_center[1] + cell_size_n / 2)
    print(origin, test_origin)
    assert origin == test_origin


#@pytest.fixture
#def grid_centres():
#    centres = [-45., -25., -17.5, -13.5, -11., -9.5, -8.5, -7.5, -6.5, -5.5, -4.5, -3.5, -2.5,
#               -1.5, -0.5, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 11., 13.5,
#               17.5, 25., 45.]
#    stripped = [-9.5, -8.5, -7.5, -6.5, -5.5, -4.5, -3.5, -2.5, -1.5, -0.5, 0.5, 1.5, 2.5, 3.5,
#                4.5, 5.5, 6.5, 7.5, 8.5, 9.5]
#    stripped_keep_start = [-45., -25., -17.5, -13.5, -11., -9.5, -8.5, -7.5, -6.5, -5.5, -4.5,
#                           -3.5, -2.5, -1.5, -0.5, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5,
#                           9.5]
#    return centres, stripped, stripped_keep_start
#
#
#def test_get_centre_points(model_grid, grid_centres):
#    grid, _, _, _ = model_grid
#    ce = conv._get_centres(grid)
#
#    # Should have one less centre point then there are nodes as the
#    #  last node is the terminating boundary and has no centre.
#    np.testing.assert_array_equal(ce, grid_centres[0])
#
#
#def test_strip_padding(model_grid, grid_centres):
#    _, pad, _, _ = model_grid
#    centres, stripped, stripped_keep_start = grid_centres
#
#    s = conv._strip_padding(centres, pad)
#    # Stripped array should have padding removed from either side
#    np.testing.assert_array_equal(s, stripped)
#
#    s = conv._strip_padding(centres, pad, keep_start=True)
#    # Stripped array while keeping start (for stripping Z padding)
#    #  should only have padding removed from end
#    np.testing.assert_array_equal(s, stripped_keep_start)
#
#
#def test_gdal_origin(model_grid, grid_centres):
#    _, _, cell_size, center_point = model_grid
#    _, stripped, _ = grid_centres
#    # Origin should be upper-left, so western-most point shifted 1/2
#    #  cell west and northern-most point shifted 1/2 cell north
#    origin = conv._get_gdal_origin(stripped, cell_size, center_point[0],
#                                   stripped, cell_size, center_point[1])
#    assert origin == (-10.0, 10.0)

