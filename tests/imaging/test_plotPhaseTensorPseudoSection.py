import glob
import os

import matplotlib.pyplot as plt
import pytest

from mtpy.imaging.phase_tensor_pseudosection import PlotPhaseTensorPseudoSection
from tests.imaging import ImageTestCase, ImageCompare


def _expected_compare_fail():
    pytest.xfail(
        "expected the image to be different on different platform, please check the image manually."
    )


test_params = [
    (
        "data/edifiles",
        {
            "fig_size": (8, 8),
            "savefig_kwargs": {"dpi": 100},
            "on_compare_fail": _expected_compare_fail,
        },
    ),
    (
        "examples/data/edi2",
        {
            "fig_size": (5, 8),
            "savefig_kwargs": {"dpi": 100},
            "on_compare_fail": _expected_compare_fail,
        },
    ),
    ("examples/data/edi_files", {"fig_size": (8, 6), "savefig_kwargs": {"dpi": 100}}),
    (
        "../MT_Datasets/3D_MT_data_edited_fromDuanJM",
        {"fig_size": (8, 6), "savefig_kwargs": {"dpi": 100}},
    ),
    (
        "../MT_Datasets/GA_UA_edited_10s-10000s",
        {"fig_size": (8, 6), "savefig_kwargs": {"dpi": 100}},
    ),
]


class TestPlotPhaseTensorPseudoSection(ImageTestCase):
    pass


def _test_gen(edi_path):
    def default(self):
        edi_file_list = glob.glob(os.path.join(edi_path, "*.edi"))
        freq = 1
        ptpObj = PlotPhaseTensorPseudoSection(
            fn_list=edi_file_list,
            tscale="period",
            # ylim=(1e-1, 1e3),  # orig period
            # range to plot
            # period range to plot
            ylim=(0, 10000),
            # xlim = (0,10000),
            stretch=(2000, 40),
            # determines (x,y) aspect ratio of plot
            station_id=(0, 10),  # indices for showing station names
            ellipse_dict={"size": 6},
            plot_tipper="yri",
            arrow_dict={"size": 5, "head_length": 0.2, "head_width": 0.1, "lw": 0.5},
            # arrow parameters, adjust as
            # necessary. lw = linewidth
            font_size=4,
            fig_size=(8, 6),
            fig_dpi=100,
        )
        ptpObj.plot()
        plt.pause(1)
        save_figure_name = "{}.png".format(default.__name__)
        save_figure_path = os.path.join(self._temp_dir, save_figure_name)
        ptpObj.save_figure2(save_fn=save_figure_path, close_plot="n")
        assert os.path.isfile(save_figure_path)

    return (default,)


# generate tests
for edi_path, img_kwargs in test_params:
    if os.path.isdir(edi_path):
        test_name = os.path.basename(edi_path)
        for _test_func in _test_gen(edi_path):
            _test_func.__name__ = "test_{test_name}_{plot_name}".format(
                test_name=test_name, plot_name=_test_func.__name__
            )
            setattr(
                TestPlotPhaseTensorPseudoSection,
                _test_func.__name__,
                ImageCompare(**img_kwargs).__call__(_test_func),
            )
