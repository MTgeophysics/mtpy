"""
Tests for resistivity model to geotiff.
"""
import pytest
import numpy as np

from mtpy.utils import convert_modem_data_to_geogrid as conv


def test_strip_padding():
    grid_axis = [-2.0, -1.0, 0.0, 1.0, 2.0]
    test = conv._strip_padding(grid_axis, 1)
    expected = [-1.0, 0.0, 1.0]
    np.testing.assert_array_equal(test, expected)

    # For 0 padding case, the array should come back untouched.
    test = conv._strip_padding(grid_axis, 0)
    expected = grid_axis
    np.testing.assert_array_equal(test, expected)

    # When keep_start=True, padding should only be removed from the
    #  end. This is used for stripping padding from Z-axis (only has
    #  padding cells on the bottom).
    test = conv._strip_padding(test, 1, keep_start=True)
    expected = [-2.0, -1.0, 0.0, 1.0]
    np.testing.assert_array_equal(test, expected)


def test_gdal_origin():
    grid_east = [-2.0, -1.0, 0.0, 1.0, 2.0]
    grid_north = [-1.5, 0.0, 1.5]
    # Center of the survey site, puts the grid into reference
    center_point = 20.0, 150.0
    test = conv._get_gdal_origin(
        grid_east, 1.0, center_point[1], grid_north, 1.5, center_point[0]
    )
    # Origin should be upper-left, so western-most point shifted 1/2
    #  cell west and northern-most point shifted 1/2 cell north.
    expected = [147.5, 22.25]

    np.testing.assert_array_equal(test, expected)


def test_target_grid_generation():
    grid_east = [-2.0, -1.0, 0.0, 1.0, 2.0]
    grid_north = [-1.5, 0.0, 1.5]
    test_grid_x, test_grid_y = conv._build_target_grid(grid_east, 1.0, grid_north, 1.5)
    expected_x = np.array(
        [
            [-2.0, -1.0, 0.0, 1.0, 2],
            [-2.0, -1.0, 0.0, 1.0, 2],
            [-2.0, -1.0, 0.0, 1.0, 2],
        ]
    )
    expected_y = np.array(
        [
            [-1.5, -1.5, -1.5, -1.5, -1.5],
            [0.0, 0.0, 0.0, 0.0, 0.0],
            [1.5, 1.5, 1.5, 1.5, 1.5],
        ]
    )

    # Testing to make sure the axes are in the correct order.
    np.testing.assert_array_equal(test_grid_x, expected_x)
    np.testing.assert_array_equal(test_grid_y, expected_y)


def test_strip_resgrid():
    resgrid = np.array(
        [
            [
                [100.0, 100.0, 100.0, 100.0, 100.0],
                [150.0, 150.0, 150.0, 150.0, 150.0],
                [100.0, 100.0, 100, 100.0, 100.0],
            ],
            [
                [150.0, 150.0, 150.0, 150.0, 150.0],
                [100.0, 100.0, 100.0, 100.0, 100.0],
                [150.0, 150.0, 150.0, 150.0, 150.0],
            ],
        ]
    )
    # Transpose to get [[Y], [X], [Z]]
    resgrid = resgrid.T
    test = conv._strip_resgrid(resgrid, 1, 1, 1)
    # We're padding one cell off either side of X and Y, and one cell
    #  off the end of Z. So we're left with the [150.] array in middle
    #  of the top row, with two elements from either side removed.
    expected = np.array([[[150.0, 150.0, 150.0]]]).T  # Again, note the transpose.

    np.testing.assert_array_equal(test, expected)


def test_get_depth_indicies():
    z_cells = np.asarray([0, 100, 500, 1000, 10000])

    assert conv._get_depth_indicies(z_cells, [49]) == set([0])
    assert conv._get_depth_indicies(z_cells, [500]) == set([2])

    # Should round up to nearest depth
    assert conv._get_depth_indicies(z_cells, [50]) == set([1])

    # If two provided depths are closest to same index, only return
    #  that index (by returning a set)
    assert conv._get_depth_indicies(z_cells, [50, 51]) == set([1])

    # If no depths are provided, return every index in a set
    assert conv._get_depth_indicies(z_cells, []) == {0, 1, 2, 3, 4}


def test_interpolate_depth_slice():
    ce = np.array([-2.0, 0.0, 2.0])
    cse = 2.0
    cn = np.array([-3.0, -1.5, 0.0, 1.5, 3.0])
    csn = 1.5
    # Dummy resisitvity model
    resgrid = np.array(
        [
            [
                [100.0, 100.0, 100.0, 100.0, 100.0],
                [150.0, 150.0, 150.0, 150.0, 150.0],
                [100.0, 100.0, 100, 100.0, 100.0],
            ],
            [
                [150.0, 150.0, 150.0, 150.0, 150.0],
                [100.0, 100.0, 100.0, 100.0, 100.0],
                [150.0, 150.0, 150.0, 150.0, 150.0],
            ],
        ]
    )
    # Transpose it to get [[Y], [X], [Z]]
    resgrid = resgrid.T
    tgx, tgy = np.meshgrid(np.arange(ce[0], ce[-1], cse), np.arange(cn[0], cn[-1], csn))

    res_slice = conv._interpolate_slice(ce, cn, resgrid, 0, tgx, tgy, log_scale=True)
    expected = np.array(
        [[2.0, 2.17609126], [2.0, 2.17609126], [2.0, 2.17609126], [2.0, 2.17609126],]
    )
    assert np.allclose(res_slice, expected)

    res_slice = conv._interpolate_slice(ce, cn, resgrid, 0, tgx, tgy, log_scale=False)
    expected = np.array(
        [
            [1024.0, 2381.063868],
            [1024.0, 2381.063868],
            [1024.0, 2381.063868],
            [1024.0, 2381.063868],
        ]
    )
    assert np.allclose(res_slice, expected)


def test_rotate_geotransform():
    pixel_width, pixel_height = 5.0, 5.0
    origin_x, origin_y = 50.0, 100.0
    angle = 90.0
    # GDAL transform:
    #  upperleft X, pixel width, row rotation, upperleft Y, column rotation, pixel height
    gt = [origin_x, pixel_width, 0, origin_y, 0, pixel_height]

    # Rotate about the upper-left
    test = conv._rotate_transform(gt, angle, origin_x, origin_y)
    expected = [gt[0], 3.061616997868383e-16, -5.0, gt[3], 5.0, 3.061616997868383e-16]

    assert test == expected

    # Rotate about a center point of (0., 0.)
    test = conv._rotate_transform(gt, angle, 0.0, 0.0)
    expected = [
        100.0,
        3.061616997868383e-16,
        -5.0,
        -49.99999999999999,
        5.0,
        3.061616997868383e-16,
    ]
