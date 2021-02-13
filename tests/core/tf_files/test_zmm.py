# -*- coding: utf-8 -*-
"""
Created on Mon Feb  1 15:54:38 2021

:copyright: 
    Jared Peacock (jpeacock@usgs.gov)

:license: MIT

"""
# =============================================================================
# Imports
# =============================================================================
import unittest
from mtpy.core import mt
from mtpy.core.io.zmm import ZMM
from tests import EDI_DATA_DIR_BB

# =============================================================================
# test j file
# =============================================================================
class TestZMM(unittest.TestCase):
    """
    test j file
    """

    def setUp(self):
        self.zmm_fn = list(EDI_DATA_DIR_BB.glob("*.zmm"))[0]
        self.zmm_obj = ZMM(self.zmm_fn)

    def test_location(self):
        self.assertAlmostEqual(34.727, self.zmm_obj.latitude, 3)
        self.assertAlmostEqual(-115.735, self.zmm_obj.longitude, 3)
        self.assertAlmostEqual(13.10, self.zmm_obj.declination, 3)
        self.assertEqual(0.0, self.zmm_obj.elevation)

    def test_station(self):
        self.assertEqual("300", self.zmm_obj.station)

    def test_n_periods(self):
        # file has 2 null periods -999 so size is 12 not 14 like in the file.
        self.assertEqual(38, self.zmm_obj.periods.size)

    def test_has_z(self):
        self.assertTrue(self.zmm_obj.Z.z.all() != 0)

    def test_has_tipper(self):
        self.assertTrue(self.zmm_obj.Tipper.tipper.all() != 0)

    def test_fail_input_fn(self):
        def set_fn(fn):
            self.zmm_obj.fn = fn

        self.assertRaises(ValueError, set_fn, r"/home/test.edi")


class TestReadZMM(unittest.TestCase):
    """
    test reading j file into MT object
    """

    def setUp(self):
        self.zmm_fn = list(EDI_DATA_DIR_BB.glob("*.zmm"))[0]
        self.zmm_obj = ZMM(self.zmm_fn)
        self.mt_obj = mt.MT(self.zmm_fn)

    def test_location(self):
        self.assertEqual(
            self.zmm_obj.station_metadata.location.latitude, self.mt_obj.latitude
        )
        self.assertEqual(
            self.zmm_obj.station_metadata.location.longitude, self.mt_obj.longitude
        )
        self.assertEqual(
            self.zmm_obj.station_metadata.location.elevation, self.mt_obj.elevation
        )

    def test_station(self):
        self.assertEqual(self.zmm_obj.station, self.mt_obj.station)

    def test_z(self):
        self.assertEqual(self.zmm_obj.Z, self.mt_obj.Z)

    def test_tipper(self):
        self.assertEqual(self.zmm_obj.Tipper, self.mt_obj.Tipper)

    # def test_birrp_parameters(self):
    #     bp_list = [f"{k} = {v}" for k, v in self.zmm_obj.header_dict.items()]
    #     self.assertEqual(bp_list,
    #                      self.mt_obj.station_metadata.transfer_function.processing_parameters)


# =============================================================================
# Run
# =============================================================================
if __name__ == "__main__":
    unittest.main()
