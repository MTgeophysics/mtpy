import glob
import os

# configure matplotlib for testing
import matplotlib.pyplot as plt

from mtpy.imaging.plotstrike import PlotStrike
from tests.imaging import ImageTestCase, ImageCompare

edi_paths = [
    "data/edifiles",
    "examples/data/edi2",
    "examples/data/edi_files",
    "../MT_Datasets/3D_MT_data_edited_fromDuanJM",
    "../MT_Datasets/GA_UA_edited_10s-10000s",
    "data/edifiles2",
]


class TestPlotStrike(ImageTestCase):
    pass


def _test_gen(edi_path):
    def default(self):
        edi_file_list = glob.glob(os.path.join(edi_path, "*.edi"))
        pt_obj = PlotStrike(
            fn_list=edi_file_list,
            plot_yn="n",
            save_figure_path=self._temp_dir,
            fig_size=(8, 6),
            fig_dpi=100,
        )
        pt_obj.plot()
        plt.pause(1)
        save_figure_name = "{}.png".format(default.__name__)
        save_figure_path = os.path.join(self._temp_dir, save_figure_name)
        pt_obj.save_plot(save_figure_path, file_format="png", close_plot="n")
        assert os.path.isfile(save_figure_path)

    def rotation(self):
        edi_file_list = glob.glob(os.path.join(edi_path, "*.edi"))
        # change rotation
        pt_obj = PlotStrike(
            fn_list=edi_file_list, plot_yn="n", rot_z=90, fig_size=(8, 6), fig_dpi=100
        )
        pt_obj.plot()
        plt.pause(1)

    def type(self):
        edi_file_list = glob.glob(os.path.join(edi_path, "*.edi"))
        # plot type
        pt_obj = PlotStrike(
            fn_list=edi_file_list,
            plot_yn="n",
            plot_type=1,
            fig_size=(8, 6),
            fig_dpi=100,
        )
        pt_obj.plot()
        plt.pause(1)

    def tipper(self):
        edi_file_list = glob.glob(os.path.join(edi_path, "*.edi"))
        # plot_tipper
        pt_obj = PlotStrike(
            fn_list=edi_file_list,
            plot_yn="n",
            plot_tipper="y",
            fig_size=(8, 6),
            fig_dpi=100,
        )
        pt_obj.plot()
        plt.pause(1)

    def fold(self):
        edi_file_list = glob.glob(os.path.join(edi_path, "*.edi"))
        # fold
        pt_obj = PlotStrike(
            fn_list=edi_file_list, plot_yn="n", fold=False, fig_size=(8, 6), fig_dpi=100
        )
        pt_obj.plot()
        plt.pause(1)

    return default, rotation, type, tipper, fold


# generate tests
for edi_path in edi_paths:
    if os.path.isdir(edi_path):
        test_name = os.path.basename(edi_path)
        for _test_func in _test_gen(edi_path):
            _test_func.__name__ = "test_{test_name}_{plot_name}".format(
                test_name=test_name, plot_name=_test_func.__name__
            )
            setattr(
                TestPlotStrike,
                _test_func.__name__,
                ImageCompare(fig_size=(8, 6), savefig_kwargs={"dpi": 100}).__call__(
                    _test_func
                ),
            )
