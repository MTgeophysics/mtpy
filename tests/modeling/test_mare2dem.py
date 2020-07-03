"""
Tests for resistivity model to geotiff.
"""
import os
import filecmp

import pytest
import numpy as np

from mtpy.utils import convert_modem_data_to_geogrid as conv
from tests import M2D_DIR


@pytest.fixture(scope='module')
def ref_output():
    return os.path.join(M2D_DIR, 'Mare2D_data.txt')


@pytest.fixture(scope='module')
def test_output():
    return '/tmp/mare2dem_test.txt'


def test_mare2dem_data(ref_output, test_output):
    assert (filecmp.cmp(ref_output, test_output))


