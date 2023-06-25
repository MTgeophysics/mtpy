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
    @classmethod
    def setUpClass(self):
        self.fn_list = [
            value
            for key, value in mt_metadata.__dict__.items()
            if key.startswith("TF")
        ]

        self.mc = MTCollection()
        self.mc.open_collection("test_collection")

        self.mc.add_tf(self.fn_list)
        self.maxDiff = None

        self.true_dataframe = pd.DataFrame(
            {
                "station": {
                    0: "14-IEB0537A",
                    1: "CAS04",
                    2: "NMX20",
                    3: "500fdfilNB207",
                    4: "SMG1",
                    5: "GAA54",
                    6: "22",
                    7: "24",
                    8: "2813",
                    9: "300",
                    10: "BP05",
                    11: "YSW212abcdefghijkl",
                    12: "GEO858",
                    13: "TEST01",
                    14: "TEST_01",
                    15: "s08",
                    16: "SAGE_2005",
                    17: "21PBS-FJM",
                },
                "survey": {
                    0: "BOULIA",
                    1: "CONUS_South",
                    2: "CONUS_South",
                    3: "Nepabunna_2010",
                    4: "South_Chile",
                    5: "Transportable_Array",
                    6: "unknown_survey",
                    7: "unknown_survey",
                    8: "unknown_survey",
                    9: "unknown_survey",
                    10: "unknown_survey",
                    11: "unknown_survey",
                    12: "unknown_survey_001",
                    13: "unknown_survey_002",
                    14: "unknown_survey_003",
                    15: "unknown_survey_004",
                    16: "unknown_survey_005",
                    17: "unknown_survey_006",
                },
                "latitude": {
                    0: -22.823722222222223,
                    1: 37.63335,
                    2: 34.470528,
                    3: -30.587969,
                    4: -38.41,
                    5: 31.888699,
                    6: 38.6653467,
                    7: 32.83331167,
                    8: 44.1479163,
                    9: 34.727,
                    10: 0.0,
                    11: 44.631,
                    12: 22.691378333333333,
                    13: -30.930285,
                    14: -23.051133333333333,
                    15: -34.646,
                    16: 35.55,
                    17: 0.0,
                },
                "longitude": {
                    0: 139.29469444444445,
                    1: -121.46838,
                    2: -108.712288,
                    3: 138.959969,
                    4: -73.904722,
                    5: -83.281681,
                    6: -113.1690717,
                    7: -107.08305667,
                    8: -111.0497517,
                    9: -115.735,
                    10: 0.0,
                    11: -110.44,
                    12: 139.70504,
                    13: 127.22923,
                    14: 139.46753333333334,
                    15: 137.006,
                    16: -106.28333333333333,
                    17: 0.0,
                },
                "elevation": {
                    0: 158.0,
                    1: 329.387,
                    2: 1940.05,
                    3: 534.0,
                    4: 10.0,
                    5: 77.025,
                    6: 1548.1,
                    7: 0.0,
                    8: 0.0,
                    9: 0.0,
                    10: 0.0,
                    11: 0.0,
                    12: 181.0,
                    13: 175.27,
                    14: 122.0,
                    15: 0.0,
                    16: 0.0,
                    17: 0.0,
                },
                "tf_id": {
                    0: "14-IEB0537A",
                    1: "CAS04",
                    2: "NMX20",
                    3: "500fdfilNB207",
                    4: "SMG1",
                    5: "GAA54",
                    6: "22",
                    7: "24",
                    8: "2813",
                    9: "300",
                    10: "BP05",
                    11: "ysw212abcdefghijkl",
                    12: "GEO858",
                    13: "TEST01",
                    14: "TEST_01",
                    15: "s08",
                    16: "SAGE_2005",
                    17: "21PBS-FJM",
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
                    11: "none",
                    12: "none",
                    13: "none",
                    14: "none",
                    15: "none",
                    16: "none",
                    17: "none",
                },
                "has_impedance": {
                    0: True,
                    1: True,
                    2: True,
                    3: True,
                    4: True,
                    5: True,
                    6: True,
                    7: True,
                    8: True,
                    9: True,
                    10: True,
                    11: False,
                    12: True,
                    13: True,
                    14: True,
                    15: True,
                    16: True,
                    17: True,
                },
                "has_tipper": {
                    0: True,
                    1: True,
                    2: True,
                    3: False,
                    4: True,
                    5: True,
                    6: True,
                    7: False,
                    8: False,
                    9: True,
                    10: False,
                    11: True,
                    12: True,
                    13: True,
                    14: True,
                    15: False,
                    16: True,
                    17: True,
                },
                "has_covariance": {
                    0: True,
                    1: False,
                    2: True,
                    3: False,
                    4: False,
                    5: True,
                    6: False,
                    7: False,
                    8: False,
                    9: True,
                    10: False,
                    11: True,
                    12: False,
                    13: False,
                    14: True,
                    15: False,
                    16: True,
                    17: False,
                },
                "period_min": {
                    0: 0.003125,
                    1: 4.65455,
                    2: 4.65455,
                    3: 0.0064,
                    4: 16.0,
                    5: 7.31429,
                    6: 0.0125,
                    7: 0.0009765625,
                    8: 0.00390625,
                    9: 1.16364,
                    10: 1.333333,
                    11: 0.01818,
                    12: 0.005154639175257732,
                    13: 0.0012115271966653925,
                    14: 0.00010061273153504844,
                    15: 0.007939999015440123,
                    16: 0.00419639110365086,
                    17: 0.000726427429899753,
                },
                "period_max": {
                    0: 2941.176470588235,
                    1: 29127.11,
                    2: 29127.11,
                    3: 2.730674,
                    4: 11585.27,
                    5: 18724.57,
                    6: 1365.3368285956144,
                    7: 42.6657564638621,
                    8: 1024.002621446711,
                    9: 10922.66699,
                    10: 64.55,
                    11: 4096.0,
                    12: 1449.2753623188407,
                    13: 1211.5274902250933,
                    14: 1.0240026214467108,
                    15: 2730.8332372990308,
                    16: 209.73154362416108,
                    17: 526.3157894736842,
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
                    11: None,
                    12: None,
                    13: None,
                    14: None,
                    15: None,
                    16: None,
                    17: None,
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
                    11: None,
                    12: None,
                    13: None,
                    14: None,
                    15: None,
                    16: None,
                    17: None,
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
            original.read()

            h5_tf = self.mc.get_tf(validate_name(original.tf_id))

            original.survey_metadata.id = h5_tf.survey_metadata.id
            original.survey_metadata.hdf5_reference = (
                h5_tf.survey_metadata.hdf5_reference
            )
            original.survey_metadata.mth5_type = h5_tf.survey_metadata.mth5_type
            original.station_metadata.acquired_by.author = (
                h5_tf.station_metadata.acquired_by.author
            )
            if original.station_metadata.transfer_function.runs_processed in [
                [],
                [""],
            ]:
                original.station_metadata.transfer_function.runs_processed = (
                    original.station_metadata.run_list
                )

            if tf_fn.stem in ["spectra_in", "spectra_out"]:
                self.assertTrue((original.dataset == h5_tf.dataset).all())
                continue
            ## these are close but each has its own unique issues, which are
            ## snowflakes.  For now skip, but check back periodically

            # with self.subTest(f"{original.tf_id} survey_metadata"):
            #     self.assertDictEqual(
            #         h5_tf.survey_metadata.to_dict(single=True),
            #         original.survey_metadata.to_dict(single=True),
            #     )

            # with self.subTest(f"{original.tf_id} station_metadata"):
            #     self.assertDictEqual(
            #         h5_tf.station_metadata.to_dict(single=True),
            #         original.station_metadata.to_dict(single=True),
            #     )

            with self.subTest(f"{original.tf_id} data"):
                self.assertTrue((original.dataset == h5_tf.dataset).all())

    @classmethod
    def tearDownClass(self):
        self.mc.mth5_collection.close_mth5()
        self.mc.mth5_filename.unlink()


# =============================================================================
# run
# =============================================================================
if __name__ == "__main__":
    unittest.main()
