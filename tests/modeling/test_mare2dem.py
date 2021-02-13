"""
Tests for EDI to Mare2D conversion.
"""
import os
import filecmp
import tempfile
import shutil
import difflib
from copy import deepcopy

import pytest
import numpy as np

from mtpy.utils import convert_modem_data_to_geogrid as conv
from mtpy.modeling import occam2d as o2d
from mtpy.modeling import mare2dem as m2d
from tests import M2D_DIR, EDI_DATA_DIR2, AUS_TOPO_FILE


@pytest.fixture()
def ref_output():
    return os.path.join(M2D_DIR, "Mare2Ddata.dat")


@pytest.fixture()
def test_output():
    tmpdir = tempfile.mkdtemp()

    # Full path to save Occam2D data file
    o2d_path = os.path.join(tmpdir, "o2d_data.dat")
    rot_o2d_path = os.path.join(tmpdir, "rot_o2d_data.dat")

    # Full path to save Mare2D data file
    m2d_path = os.path.join(tmpdir, "mare2dem_test.txt")

    # Generate an Occam2D data object from EDI data
    o2d_data = o2d.Data(
        edi_path=EDI_DATA_DIR2,
        model_mode="1",
        optimize_line=True,
        interpolate_freq=False,
        res_te_err=20.0,
        phase_te_err=10.0,
        res_tm_err=10.0,
        phase_tm_err=5.0,
    )

    # Save the data file
    # We need a data file with the non-rotated profile and the rotated
    # profile. This is because the elevation will be interpolated over
    # the non-projected profile and stations.
    rot_o2d_data = deepcopy(
        o2d_data
    )  # Make a copy because 'write_data_file' will populate data
    o2d_data._rotate_to_strike = False
    o2d_data.save_path = o2d_path
    rot_o2d_data.save_path = rot_o2d_path
    o2d_data.write_data_file(data_fn=o2d_path)
    rot_o2d_data.write_data_file(data_fn=rot_o2d_path)

    gstrike = o2d_data.geoelectric_strike

    # Convert the Occam2D profile to Mare2D
    (
        mare_origin,
        utm_zone,
        site_locations,
        site_elevations,
        site_names,
        m2d_profile,
        profile_elevation,
    ) = m2d.occam2d_to_mare2dem(
        o2d_data, rot_o2d_data, AUS_TOPO_FILE, elevation_sample_n=300
    )

    m2d.write_mare2dem_data(
        o2d_path,
        site_locations,
        site_elevations,
        site_names,
        mare_origin,
        utm_zone,
        gstrike,
        solve_statics=False,
        savepath=m2d_path,
    )

    yield m2d_path

    if os.path.exists(tmpdir):
        shutil.rmtree(tmpdir)


def test_mare2dem_data(ref_output, test_output):
    files_are_same = filecmp.cmp(ref_output, test_output)
    if not files_are_same:
        print("File comparison failed, testing within tolerance")
        with open(ref_output) as r, open(test_output) as t:
            diff = list(difflib.unified_diff(r.readlines(), t.readlines()))
            # Test X, Y, Z are within tolerance (2 decimal places)
            for i, line in enumerate(diff):
                if line.startswith("-"):
                    if diff[i + 1].startswith("+"):
                        a = line.split()
                        b = diff[i + 1].split()
                        ax, ay, az = float(a[0]), float(a[1]), float(a[2])
                        bx, by, bz = float(b[0]), float(b[1]), float(b[2])
                        files_are_same = (
                            np.testing.assert_almost_equal(ax, bx, decimal=2)
                            and np.testing.assert_almost_equal(ay, by, decimal=2)
                            and np.testing.assert_almost_equal(az, bz, decimal=2)
                        )
            if not files_are_same:
                print(
                    "File comparison failed and values out of tolerance, printing diff"
                )
                diff = difflib.unified_diff(r.readlines(), t.readlines())
                for line in diff:
                    print(line)
    assert files_are_same


@pytest.fixture()
def coords():
    x = [1, 2, 3, 4]
    y = [2, 5, 7, 8]
    return x, y


def test_line_length(coords):
    x, y = coords
    test = m2d.line_length(x[0], y[0], x[-1], y[-1])
    expected = 6.708
    np.testing.assert_almost_equal(test, expected, decimal=3)


def test_points_o2d_to_m2d(coords):
    x, y = coords
    test = m2d.points_o2d_to_m2d(x, y)
    expected = [-3.354, -0.192, 2.031, 3.354]
    np.testing.assert_almost_equal(test, expected, decimal=3)
