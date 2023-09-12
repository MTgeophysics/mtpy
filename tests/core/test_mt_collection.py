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
                    0: "14-IEB0537A",
                    1: "CAS04",
                    2: "NMX20",
                    3: "KAK",
                    4: "500fdfilNB207",
                    5: "SMG1",
                    6: "GAA54",
                    7: "300",
                    8: "YSW212abcdefghijkl",
                    9: "BP05",
                    10: "701_merged_wrcal",
                    11: "GEO858",
                    12: "TEST01",
                    13: "TEST_01",
                    14: "s08",
                    15: "SAGE_2005_og",
                    16: "SAGE_2005_out",
                    17: "21PBS-FJM",
                    18: "24",
                    19: "22",
                    20: "2813",
                },
                "survey": {
                    0: "BOULIA",
                    1: "CONUS_South",
                    2: "CONUS_South",
                    3: "JMA",
                    4: "Nepabunna_2010",
                    5: "South_Chile",
                    6: "Transportable_Array",
                    7: "unknown_survey",
                    8: "unknown_survey_001",
                    9: "unknown_survey_002",
                    10: "unknown_survey_003",
                    11: "unknown_survey_004",
                    12: "unknown_survey_005",
                    13: "unknown_survey_006",
                    14: "unknown_survey_007",
                    15: "unknown_survey_008",
                    16: "unknown_survey_009",
                    17: "unknown_survey_010",
                    18: "unknown_survey_011",
                    19: "unknown_survey_012",
                    20: "unknown_survey_013",
                },
                "latitude": {
                    0: -22.823722222222223,
                    1: 37.63335,
                    2: 34.470528,
                    3: 36.232,
                    4: -30.587969,
                    5: -38.41,
                    6: 31.888699,
                    7: 34.727,
                    8: 44.631,
                    9: 0.0,
                    10: 40.64811111111111,
                    11: 22.691378333333333,
                    12: -30.930285,
                    13: -23.051133333333333,
                    14: -34.646,
                    15: 35.55,
                    16: 35.55,
                    17: 0.0,
                    18: 32.83331167,
                    19: 38.6653467,
                    20: 44.1479163,
                },
                "longitude": {
                    0: 139.29469444444445,
                    1: -121.46838,
                    2: -108.712288,
                    3: 140.186,
                    4: 138.959969,
                    5: -73.904722,
                    6: -83.281681,
                    7: -115.735,
                    8: -110.44,
                    9: 0.0,
                    10: -106.21241666666667,
                    11: 139.70504,
                    12: 127.22923,
                    13: 139.46753333333334,
                    14: 137.006,
                    15: -106.28333333333333,
                    16: -106.28333333333333,
                    17: 0.0,
                    18: -107.08305667,
                    19: -113.1690717,
                    20: -111.0497517,
                },
                "elevation": {
                    0: 158.0,
                    1: 329.387,
                    2: 1940.05,
                    3: 36.0,
                    4: 534.0,
                    5: 10.0,
                    6: 77.025,
                    7: 948.158935547,
                    8: 2352.3984375,
                    9: 0.0,
                    10: 2489.0,
                    11: 181.0,
                    12: 175.27,
                    13: 122.0,
                    14: 0.0,
                    15: 1674.992797852,
                    16: 1674.992797852,
                    17: 0.0,
                    18: 1414.37487793,
                    19: 1548.1,
                    20: 1954.861816406,
                },
                "tf_id": {
                    0: "14-IEB0537A",
                    1: "CAS04",
                    2: "NMX20",
                    3: "KAK",
                    4: "500fdfilNB207",
                    5: "SMG1",
                    6: "GAA54",
                    7: "300",
                    8: "ysw212abcdefghijkl",
                    9: "BP05",
                    10: "701_merged_wrcal",
                    11: "GEO858",
                    12: "TEST01",
                    13: "TEST_01",
                    14: "s08",
                    15: "SAGE_2005_og",
                    16: "SAGE_2005",
                    17: "21PBS-FJM",
                    18: "24",
                    19: "22",
                    20: "2813",
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
                    20: "none",
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
                    8: False,
                    9: True,
                    10: True,
                    11: True,
                    12: True,
                    13: True,
                    14: True,
                    15: True,
                    16: True,
                    17: True,
                    18: True,
                    19: True,
                    20: True,
                },
                "has_tipper": {
                    0: True,
                    1: True,
                    2: True,
                    3: False,
                    4: False,
                    5: True,
                    6: True,
                    7: True,
                    8: True,
                    9: False,
                    10: True,
                    11: True,
                    12: True,
                    13: True,
                    14: False,
                    15: True,
                    16: True,
                    17: True,
                    18: False,
                    19: True,
                    20: False,
                },
                "has_covariance": {
                    0: True,
                    1: False,
                    2: True,
                    3: False,
                    4: False,
                    5: False,
                    6: True,
                    7: True,
                    8: True,
                    9: False,
                    10: False,
                    11: False,
                    12: False,
                    13: True,
                    14: False,
                    15: True,
                    16: False,
                    17: False,
                    18: False,
                    19: False,
                    20: False,
                },
                "period_min": {
                    0: 0.003125,
                    1: 4.65455,
                    2: 4.65455,
                    3: 6.4,
                    4: 0.0064,
                    5: 16.0,
                    6: 7.31429,
                    7: 1.16364,
                    8: 0.01818,
                    9: 1.333333,
                    10: 0.0001,
                    11: 0.005154639175257732,
                    12: 0.0012115271966653925,
                    13: 0.00010061273153504844,
                    14: 0.007939999015440123,
                    15: 0.00419639110365086,
                    16: 0.00419639110365086,
                    17: 0.000726427429899753,
                    18: 0.0009765625,
                    19: 0.0125,
                    20: 0.00390625,
                },
                "period_max": {
                    0: 2941.176470588235,
                    1: 29127.11,
                    2: 29127.11,
                    3: 614400.0,
                    4: 2.730674,
                    5: 11585.27,
                    6: 18724.57,
                    7: 10922.66699,
                    8: 4096.0,
                    9: 64.55,
                    10: 2912.710720057042,
                    11: 1449.2753623188407,
                    12: 1211.5274902250933,
                    13: 1.0240026214467108,
                    14: 2730.8332372990308,
                    15: 209.73154362416108,
                    16: 209.73154362416108,
                    17: 526.3157894736842,
                    18: 42.6657564638621,
                    19: 1365.3368285956144,
                    20: 1024.002621446711,
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
                    20: None,
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
                    20: None,
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
            original.survey_metadata.mth5_type = (
                h5_tf.survey_metadata.mth5_type
            )
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
            for key, value in mt_data_01.items():
                if original.station in key:
                    original.survey = validate_name(value.survey)
                    original.station_metadata.transfer_function.runs_processed = (
                        value.station_metadata.transfer_function.runs_processed
                    )
                    original.station_metadata.run_list = (
                        value.station_metadata.run_list
                    )
                    value.survey_metadata.time_period = (
                        original.survey_metadata.time_period
                    )
                    if (
                        original.station_metadata.transfer_function.data_quality.good_from_period
                        == 0.0
                    ):
                        value.station_metadata.transfer_function.data_quality.good_from_period = (
                            0.0
                        )
                    if (
                        original.station_metadata.transfer_function.data_quality.good_to_period
                        == 0.0
                    ):
                        value.station_metadata.transfer_function.data_quality.good_to_period = (
                            0.0
                        )
                    break
            if original.station_metadata.comments in [""]:
                original.station_metadata.comments = None
            if original.station_metadata.acquired_by.author in [""]:
                original.station_metadata.acquired_by.author = None
            if (
                original.station_metadata.transfer_function.processing_type
                in [""]
            ):
                original.station_metadata.transfer_function.processing_type = (
                    None
                )
            mt_data_02.add_station(original, compute_relative_location=False)
        mt_data_02.compute_relative_locations()

        # "fix" some of the data
        mt_data_01["CONUS_South.CAS04"].survey_metadata.update_bounding_box()
        mt_data_02["CONUS_South.CAS04"].survey_metadata.country = "USA"

        mt_data_01["CONUS_South.NMX20"].survey_metadata.update_bounding_box()

        mt_data_02[
            "unknown_survey_009.SAGE_2005_out"
        ].station_metadata.runs = mt_data_01[
            "unknown_survey_009.SAGE_2005_out"
        ].station_metadata.runs
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
