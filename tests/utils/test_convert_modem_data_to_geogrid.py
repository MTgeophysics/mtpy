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
    (1., 0, -10., 10.),
    (1., 2, -7., 7.)
]

NORTH_GRID_PARAMS = [
    (1., 5, 15., 15.),
    (1., 0, -15., 15.),
    (1., 4, -8., 8.)
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

    # Should have one less than center point than elements in grid.
    #  last element of the grid is the terminating boundary so it has
    #  no center pont.
    assert len(ce) == len(grid) - 1
    np.testing.assert_array_equal(ce, centers)


def test_strip_padding(model_grid_e):
    grid, _, padding = model_grid_e
    stripped = conv._strip_padding(grid, padding)
    if padding == 0:
        test_stripped = grid
    else:
        test_stripped = grid[padding:-padding]

    # Padding is removed from start and end.
    assert len(stripped) == len(grid) - padding * 2
    np.testing.assert_array_equal(stripped, test_stripped)


def test_strip_padding_keep_start(model_grid_e):
    grid, _, padding = model_grid_e
    stripped = conv._strip_padding(grid, padding, keep_start=True)
    if padding == 0:
        test_stripped = grid
    else:
        test_stripped = grid[:-padding]

    # Padding is only removed from end.
    assert len(stripped) == len(grid) - padding
    np.testing.assert_array_equal(stripped, test_stripped)


def test_gdal_origin(model_grid_e, model_grid_n, grid_center):
    grid_e, cell_size_e, _ = model_grid_e
    grid_n, cell_size_n, _ = model_grid_n
    origin = conv._get_gdal_origin(grid_e, cell_size_e, grid_center[0],
                                   grid_n, cell_size_n, grid_center[1])

    # Origin should be upper-left, so western-most point shifted 1/2
    #  cell west and northern-most point shifted 1/2 cell north
    test_origin = (grid_e[0] + grid_center[0] - cell_size_e / 2,
                   grid_n[-1] + grid_center[1] + cell_size_n / 2)
    np.testing.assert_array_equal(origin, test_origin)


def test_target_grid_generation(model_grid_e, model_grid_n):
    ge, cse, _ = model_grid_e
    gn, csn, _ = model_grid_n
    target_grid_x, target_grid_y = conv._build_target_grid(ge, cse, gn, csn)
    test_tgx, test_tgy = np.meshgrid(np.arange(ge[0], ge[-1], cse),
                                     np.arange(gn[0], gn[-1], csn))

    # Testing to make sure the axes are in the correct order.
    np.testing.assert_array_equal(target_grid_x, test_tgx)
    np.testing.assert_array_equal(target_grid_y, test_tgy)


def test_strip_resgrid(model_grid_e, model_grid_n):
    grid_n, _, pad_n = model_grid_n
    grid_e, _, pad_e = model_grid_e
    grid_z, _, pad_z = model_grid_e
    resgrid = np.zeros((grid_n.size, grid_e.size, grid_z.size))
    resgrid_stripped = conv._strip_resgrid(resgrid, pad_n, pad_e, pad_z)

    # Resgrid is [[n], [e], [z]]. Need to make sure the correct padding is
    #  bring stripped from the correct axis (this caused bugs in the
    #  past).
    assert resgrid_stripped.shape == (len(grid_n) - pad_n * 2,
                                      len(grid_e) - pad_e * 2,
                                      len(grid_z) - pad_z)
    pad_n = slice(None) if pad_n == 0 else slice(pad_n, -pad_n)
    pad_e = slice(None) if pad_e == 0 else slice(pad_e, -pad_e)
    pad_z = slice(None) if pad_z == 0 else slice(None, -pad_z)
    np.testing.assert_array_equal(resgrid_stripped, resgrid[pad_n, pad_e, pad_z])


def test_get_depth_indicies():
    z_cells = np.asarray([0, 100, 500, 1000, 10000])

    assert conv._get_depth_indicies(z_cells, [49]) == set([0])
    assert conv._get_depth_indicies(z_cells, [500]) == set([2])

    # Should round up to nearest depth
    assert conv._get_depth_indicies(z_cells, [50]) == set([1])

    # If two provided depths are closest to same index, only return
    #  that index (by returning a set)
    assert conv._get_depth_indicies(z_cells, [50, 51]) == set([1])

    # If no depths are provided, return every index in a list
    assert conv._get_depth_indicies(z_cells, []) == [0, 1, 2, 3, 4]


def test_interpolate_depth_slice():
    ce = np.array([-2., 0., 2.])
    cse = 2.
    cn = np.array([-3., -1.5, 0., 1.5, 3.])
    csn = 1.5
    cz = np.array([0., 10.])
    # Dummy resisitvity model
    resgrid = np.array([
        [[100., 100., 100., 100., 100.], [150., 150., 150., 150., 150.], [100., 100., 100, 100., 100.]],
        [[150., 150., 150., 150., 150.], [100., 100., 100., 100., 100.], [150., 150., 150., 150., 150.]]
    ])
    # Transpose it to get [[Y], [X], [Z]]
    resgrid = resgrid.T
    tgx, tgy = np.meshgrid(np.arange(ce[0], ce[-1], cse),
                           np.arange(cn[0], cn[-1], csn))
    
    res_slice = conv._interpolate_slice(ce, cn, resgrid, 0, tgx, tgy, log_scale=True)
    expected = np.array([
        [2., 2.17609126],
        [2., 2.17609126],
        [2., 2.17609126],
        [2., 2.17609126],
    ])
    assert np.allclose(res_slice, expected)

    res_slice = conv._interpolate_slice(ce, cn, resgrid, 0, tgx, tgy, log_scale=False)
    expected = np.array(
        [[1024., 2381.063868],
         [1024., 2381.063868],
         [1024., 2381.063868],
         [1024., 2381.063868]]
    )
    assert np.allclose(res_slice, expected)

