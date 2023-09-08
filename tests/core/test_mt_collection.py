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

from mtpy import MT, MTCollection, MTData

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
                    0: "21PBS-FJM",
                    1: "22",
                    2: "24",
                    3: "2813",
                    4: "300",
                    5: "701_merged_wrcal",
                    6: "BP05",
                    7: "GEO858",
                    8: "SAGE_2005",
                    9: "TEST01",
                    10: "TEST_01",
                    11: "YSW212abcdefghijkl",
                    12: "s08",
                    13: "14-IEB0537A",
                    14: "CAS04",
                    15: "NMX20",
                    16: "KAK",
                    17: "500fdfilNB207",
                    18: "SMG1",
                    19: "GAA54",
                },
                "survey": {
                    0: "0",
                    1: "0",
                    2: "0",
                    3: "0",
                    4: "0",
                    5: "0",
                    6: "0",
                    7: "0",
                    8: "0",
                    9: "0",
                    10: "0",
                    11: "0",
                    12: "0",
                    13: "BOULIA",
                    14: "CONUS_South",
                    15: "CONUS_South",
                    16: "JMA",
                    17: "Nepabunna_2010",
                    18: "South_Chile",
                    19: "Transportable_Array",
                },
                "latitude": {
                    0: 0.0,
                    1: 38.6653467,
                    2: 32.83331167,
                    3: 44.1479163,
                    4: 34.727,
                    5: 40.64811111111111,
                    6: 0.0,
                    7: 22.691378333333333,
                    8: 35.55,
                    9: -30.930285,
                    10: -23.051133333333333,
                    11: 44.631,
                    12: -34.646,
                    13: -22.823722222222223,
                    14: 37.63335,
                    15: 34.470528,
                    16: 36.232,
                    17: -30.587969,
                    18: -38.41,
                    19: 31.888699,
                },
                "longitude": {
                    0: 0.0,
                    1: -113.1690717,
                    2: -107.08305667,
                    3: -111.0497517,
                    4: -115.735,
                    5: -106.21241666666667,
                    6: 0.0,
                    7: 139.70504,
                    8: -106.28333333333333,
                    9: 127.22923,
                    10: 139.46753333333334,
                    11: -110.44,
                    12: 137.006,
                    13: 139.29469444444445,
                    14: -121.46838,
                    15: -108.712288,
                    16: 140.186,
                    17: 138.959969,
                    18: -73.904722,
                    19: -83.281681,
                },
                "elevation": {
                    0: 0.0,
                    1: 1548.1,
                    2: 1414.37487793,
                    3: 1954.861816406,
                    4: 948.158935547,
                    5: 2489.0,
                    6: 0.0,
                    7: 181.0,
                    8: 1674.992797852,
                    9: 175.27,
                    10: 122.0,
                    11: 2352.3984375,
                    12: 0.0,
                    13: 158.0,
                    14: 329.387,
                    15: 1940.05,
                    16: 36.0,
                    17: 534.0,
                    18: 10.0,
                    19: 77.025,
                },
                "tf_id": {
                    0: "21PBS-FJM",
                    1: "22",
                    2: "24",
                    3: "2813",
                    4: "300",
                    5: "701_merged_wrcal",
                    6: "BP05",
                    7: "GEO858",
                    8: "SAGE_2005",
                    9: "TEST01",
                    10: "TEST_01",
                    11: "ysw212abcdefghijkl",
                    12: "s08",
                    13: "14-IEB0537A",
                    14: "CAS04",
                    15: "NMX20",
                    16: "KAK",
                    17: "500fdfilNB207",
                    18: "SMG1",
                    19: "GAA54",
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
                    18: "none",
                    19: "none",
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
                    18: True,
                    19: True,
                },
                "has_tipper": {
                    0: True,
                    1: True,
                    2: False,
                    3: False,
                    4: True,
                    5: True,
                    6: False,
                    7: True,
                    8: True,
                    9: True,
                    10: True,
                    11: True,
                    12: False,
                    13: True,
                    14: True,
                    15: True,
                    16: False,
                    17: False,
                    18: True,
                    19: True,
                },
                "has_covariance": {
                    0: False,
                    1: False,
                    2: False,
                    3: False,
                    4: True,
                    5: False,
                    6: False,
                    7: False,
                    8: True,
                    9: False,
                    10: True,
                    11: True,
                    12: False,
                    13: True,
                    14: False,
                    15: True,
                    16: False,
                    17: False,
                    18: False,
                    19: True,
                },
                "period_min": {
                    0: 0.000726427429899753,
                    1: 0.0125,
                    2: 0.0009765625,
                    3: 0.00390625,
                    4: 1.16364,
                    5: 0.0001,
                    6: 1.333333,
                    7: 0.005154639175257732,
                    8: 0.00419639110365086,
                    9: 0.0012115271966653925,
                    10: 0.00010061273153504844,
                    11: 0.01818,
                    12: 0.007939999015440123,
                    13: 0.003125,
                    14: 4.65455,
                    15: 4.65455,
                    16: 6.4,
                    17: 0.0064,
                    18: 16.0,
                    19: 7.31429,
                },
                "period_max": {
                    0: 526.3157894736842,
                    1: 1365.3368285956144,
                    2: 42.6657564638621,
                    3: 1024.002621446711,
                    4: 10922.66699,
                    5: 2912.710720057042,
                    6: 64.55,
                    7: 1449.2753623188407,
                    8: 209.73154362416108,
                    9: 1211.5274902250933,
                    10: 1.0240026214467108,
                    11: 4096.0,
                    12: 2730.8332372990308,
                    13: 2941.176470588235,
                    14: 29127.11,
                    15: 29127.11,
                    16: 614400.0,
                    17: 2.730674,
                    18: 11585.27,
                    19: 18724.57,
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
                    18: None,
                    19: None,
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
                    18: None,
                    19: None,
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

    def test_to_mt_data(self):
        mt_data_01 = self.mc.to_mt_data()

        mt_data_02 = MTData()
        for tf_fn in self.fn_list:
            original = MT(tf_fn)
            original.read()
            mt_data_02.add_station(original, compute_relative_location=False)
        mt_data_02.compute_relative_locations()

        self.assertEqual(mt_data_01, mt_data_02)

    @classmethod
    def tearDownClass(self):
        self.mc.mth5_collection.close_mth5()
        self.mc.mth5_filename.unlink()


# =============================================================================
# run
# =============================================================================
if __name__ == "__main__":
    unittest.main()
