import glob
import os

import matplotlib.pyplot as plt
import pytest

from mtpy.imaging.plotnresponses import PlotMultipleResponses
from tests.imaging import ImageTestCase, ImageCompare


def _expected_compare_fail():
    pytest.xfail("expected to be different, check the result image manually")


test_params = [
    (
        "data/edifiles",
        {
            "style_all": {
                "fig_size": (16, 6),
                "savefig_kwargs": {"dpi": 200},
                "on_compare_fail": _expected_compare_fail,
            },
            "style_compare": {
                "fig_size": (8, 6),
                "savefig_kwargs": {"dpi": 200},
                "on_compare_fail": _expected_compare_fail,
            },
        },
    ),
    (
        "examples/data/edi2",
        {
            "style_all": {
                "fig_size": (16, 6),
                "savefig_kwargs": {"dpi": 200},
                "on_compare_fail": _expected_compare_fail,
            },
            "style_compare": {
                "fig_size": (8, 6),
                "savefig_kwargs": {"dpi": 200},
                "on_compare_fail": _expected_compare_fail,
            },
        },
    ),
    (
        "examples/data/edi_files",
        {
            "style_all": {
                "fig_size": (16, 6),
                "savefig_kwargs": {"dpi": 200},
                "on_compare_fail": _expected_compare_fail,
            },
            "style_compare": {
                "fig_size": (8, 6),
                "savefig_kwargs": {"dpi": 200},
                "on_compare_fail": _expected_compare_fail,
            },
        },
    ),
    (
        "../MT_Datasets/3D_MT_data_edited_fromDuanJM",
        {
            "style_all": {"fig_size": (16, 6), "savefig_kwargs": {"dpi": 200}},
            "style_compare": {"fig_size": (8, 6), "savefig_kwargs": {"dpi": 200}},
        },
    ),
    (
        "../MT_Datasets/GA_UA_edited_10s-10000s",
        {
            "style_all": {
                "fig_size": (32, 6),
                "tolerance": 10,
                "savefig_kwargs": {"dpi": 200},
            },
            "style_compare": {"fig_size": (8, 6), "savefig_kwargs": {"dpi": 200}},
        },
    ),
]


class TestPlotMultipleResponses(ImageTestCase):
    pass


def _test_gen(edi_path):
    def style_all(self):
        edi_file_list = glob.glob(os.path.join(edi_path, "*.edi"))
        pt_obj = PlotMultipleResponses(
            fn_list=edi_file_list,
            plot_num=1,
            plot_tipper="yr",
            plot_style="all",
            plot_yn="n",
            fig_size=(8, 6),
            fig_dpi=100,
        )
        pt_obj.plot()
        plt.pause(0.5)

    def style_compare(self):
        edi_file_list = glob.glob(os.path.join(edi_path, "*.edi"))
        # compare
        pt_obj = PlotMultipleResponses(
            fn_list=edi_file_list,
            plot_num=1,
            plot_tipper="yr",
            plot_style="compare",
            plot_yn="n",
            fig_size=(8, 6),
            fig_dpi=100,
        )
        pt_obj.plot()
        plt.pause(0.5)

    return style_all, style_compare


# generate tests
for edi_path, img_kwargs in test_params:
    if os.path.isdir(edi_path):
        test_name = os.path.basename(edi_path)
        for _test_func in _test_gen(edi_path):
            plot_name = _test_func.__name__
            _test_func.__name__ = "test_{test_name}_{plot_name}".format(
                test_name=test_name, plot_name=plot_name
            )
            setattr(
                TestPlotMultipleResponses,
                _test_func.__name__,
                ImageCompare(**img_kwargs[plot_name]).__call__(_test_func),
            )
