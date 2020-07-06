"""
Tests for EDI to Mare2D conversion.
"""
import os
import filecmp
import tempfile
import shutil

import pytest
import numpy as np

from mtpy.utils import convert_modem_data_to_geogrid as conv
from mtpy.modeling import occam2d as o2d
from mtpy.modeling import mare2dem as m2d
from tests import M2D_DIR, EDI_DATA_DIR2, AUS_TOPO_FILE


@pytest.fixture()
def ref_output():
    return os.path.join(M2D_DIR, 'Mare2D_data.txt')


@pytest.fixture()
def test_output():
    tmpdir = tempfile.mkdtemp()

    # Full path to save Occam2D data file
    o2d_path = os.path.join(tmpdir, 'o2d_data.dat')

    # Full path to save Mare2D data file
    m2d_path = os.path.join(tmpdir, 'mare2dem_test.txt')

    # Generate an Occam2D data object from EDI data
    gstrike = -72
    station_list = m2d.station_list(EDI_DATA_DIR2)
    o2d_data = o2d.Data(edi_path=EDI_DATA_DIR2, model_mode='1', station_list=station_list,
                        interpolate_freq=False, geoelectric_strike=gstrike, res_te_err=20.,
                        phase_te_err=10., res_tm_err=10., phase_tm_err=5.)

    # Save the data file
    o2d_data.save_path = o2d_path
    o2d_data.write_data_file(data_fn=o2d_path)

    # Convert the Occam2D profile to Mare2D
    mare_origin_x, mare_origin_y, utm_zone, site_locations, site_elevations, m2d_profile, profile_elevation = \
        m2d.occam2d_to_mare2dem(o2d_data, AUS_TOPO_FILE, elevation_sample_n=300)

    m2d.write_mare2dem_data(o2d_path, site_locations, site_elevations,
                            (mare_origin_x, mare_origin_y, utm_zone),
                            gstrike, solve_statics=False, savepath=m2d_path)

    yield m2d_path

    if os.path.exists(tmpdir):
        shutil.rmtree(tmpdir)


def test_mare2dem_data(ref_output, test_output):
    assert (filecmp.cmp(ref_output, test_output))
