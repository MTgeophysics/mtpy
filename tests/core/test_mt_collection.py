# -*- coding: utf-8 -*-
"""
Collection of MT stations

Created on Mon Jan 11 15:36:38 2021

:copyright: 
    Jared Peacock (jpeacock@usgs.gov)

:license: MIT

"""
# =============================================================================
# Imports
# =============================================================================
import unittest

import pandas as pd

from mtpy import MT, MTCollection

import mt_metadata
from mth5.helpers import validate_name


# =============================================================================
#
# =============================================================================
class TestMTCollection(unittest.TestCase):
    def setUp(self):
        self.fn_list = [
            value for key, value in mt_metadata.__dict__.items() if key.startswith("TF")
        ]

        self.mc = MTCollection()
        self.mc.initialize_collection("test_collection")

        self.mc.add_tf(self.fn_list)

        self.true_dataframe = pd.DataFrame(
            {
                "station": {
                    0: "14-IEB0537A",
                    1: "NMX20",
                    2: "24",
                    3: "300",
                    4: "BP05",
                    5: "YSW212abcdefghijkl",
                    6: "GEO858",
                    7: "TEST01",
                    8: "Geoscience_Australia",
                    9: "s08",
                    10: "SAGE_2005",
                },
                "survey": {
                    0: "BOULIA",
                    1: "CONUS_South",
                    2: "unknown_survey",
                    3: "unknown_survey",
                    4: "unknown_survey",
                    5: "unknown_survey",
                    6: "unknown_survey_001",
                    7: "unknown_survey_002",
                    8: "unknown_survey_003",
                    9: "unknown_survey_004",
                    10: "unknown_survey_005",
                },
                "latitude": {
                    0: -22.823722222222223,
                    1: 34.470528,
                    2: 32.83331167,
                    3: 34.727,
                    4: 0.0,
                    5: 44.631,
                    6: 22.691378333333333,
                    7: -30.930285,
                    8: -23.051133333333333,
                    9: -34.646,
                    10: 35.55,
                },
                "longitude": {
                    0: 139.29469444444445,
                    1: -108.712288,
                    2: -107.08305667,
                    3: -115.735,
                    4: 0.0,
                    5: -110.44,
                    6: 139.70504,
                    7: 127.22923,
                    8: 139.46753333333334,
                    9: 137.006,
                    10: -106.28333333333333,
                },
                "elevation": {
                    0: 158.0,
                    1: 1940.05,
                    2: 0.0,
                    3: 0.0,
                    4: 0.0,
                    5: 0.0,
                    6: 181.0,
                    7: 175.27,
                    8: 122.0,
                    9: 0.0,
                    10: 0.0,
                },
                "tf_id": {
                    0: "14-IEB0537A",
                    1: "NMX20",
                    2: "24",
                    3: "300",
                    4: "BP05",
                    5: "ysw212abcdefghijkl",
                    6: "GEO858",
                    7: "TEST01",
                    8: "Geoscience_Australia",
                    9: "s08",
                    10: "SAGE_2005",
                },
                "units": {
                    0: "none",
                    1: "none",
                    2: "none",
                    3: "none",
                    4: "none",
                    5: "none",
                    6: "none",
                    7: "none",
                    8: "none",
                    9: "none",
                    10: "none",
                },
                "has_impedance": {
                    0: True,
                    1: True,
                    2: True,
                    3: True,
                    4: True,
                    5: False,
                    6: True,
                    7: True,
                    8: True,
                    9: True,
                    10: True,
                },
                "has_tipper": {
                    0: True,
                    1: True,
                    2: False,
                    3: True,
                    4: False,
                    5: True,
                    6: True,
                    7: True,
                    8: True,
                    9: False,
                    10: True,
                },
                "has_covariance": {
                    0: True,
                    1: True,
                    2: False,
                    3: True,
                    4: False,
                    5: True,
                    6: False,
                    7: False,
                    8: True,
                    9: False,
                    10: True,
                },
                "period_min": {
                    0: 0.003125,
                    1: 4.65455,
                    2: 0.0009765625,
                    3: 1.16364,
                    4: 1.333333,
                    5: 0.01818,
                    6: 0.005154639175257732,
                    7: 0.0012115271966653925,
                    8: 0.00010061273153504844,
                    9: 0.007939999015440123,
                    10: 0.00419639110365086,
                },
                "period_max": {
                    0: 2941.176470588235,
                    1: 29127.11,
                    2: 42.6657564638621,
                    3: 10922.66699,
                    4: 64.55,
                    5: 4096.0,
                    6: 1449.2753623188407,
                    7: 1211.5274902250933,
                    8: 1.0240026214467108,
                    9: 2730.8332372990308,
                    10: 209.73154362416108,
                },
                "hdf5_reference": {
                    0: None,
                    1: None,
                    2: None,
                    3: None,
                    4: None,
                    5: None,
                    6: None,
                    7: None,
                    8: None,
                    9: None,
                    10: None,
                },
                "station_hdf5_reference": {
                    0: None,
                    1: None,
                    2: None,
                    3: None,
                    4: None,
                    5: None,
                    6: None,
                    7: None,
                    8: None,
                    9: None,
                    10: None,
                },
            }
        )

    def test_filename(self):
        self.assertEqual(
            self.mc.working_directory.joinpath("test_collection.h5"),
            self.mc.mth5_filename,
        )

    def test_dataframe(self):
        h5_df = self.mc.dataframe[
            self.mc.dataframe.columns[
                ~self.mc.dataframe.columns.isin(
                    ["hdf5_reference", "station_hdf5_reference"]
                )
            ]
        ]

        true_df = self.true_dataframe[
            self.true_dataframe.columns[
                ~self.true_dataframe.columns.isin(
                    ["hdf5_reference", "station_hdf5_reference"]
                )
            ]
        ]

        self.assertTrue((h5_df == true_df).all().all())

    def test_get_tf(self):
        for tf_fn in self.fn_list:
            original = MT(tf_fn)

            h5_tf = self.mc.get_tf(validate_name(original.tf_id))

            original.survey_metadata.id = h5_tf.survey_metadata.id
            original.survey_metadata.hdf5_reference = (
                h5_tf.survey_metadata.hdf5_reference
            )
            original.survey_metadata.mth5_type = h5_tf.survey_metadata.mth5_type
            original.station_metadata.acquired_by.author = (
                h5_tf.station_metadata.acquired_by.author
            )

            if tf_fn.stem in ["spectra_in", "spectra_out"]:
                self.assertTrue((original.dataset == h5_tf.dataset).all())
                continue
            with self.subTest(original.tf_id):
                self.mc.logger.info(f"testing: {original.tf_id} from {tf_fn.name}")
                self.assertEqual(h5_tf, original)

    def tearDown(self):
        self.mc.mth5_collection.close_mth5()
        self.mc.mth5_filename.unlink()


# =============================================================================
# run
# =============================================================================
if __name__ == "__main__":
    unittest.main()
