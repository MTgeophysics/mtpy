import glob
import os

import matplotlib.pyplot as plt

from mtpy.core import mt
from mtpy.imaging.plot_mt_response import PlotMTResponse
from tests.imaging import ImageTestCase

edi_paths = [
    "data/edifiles",
    "examples/data/edi2",
    "examples/data/edi_files",
    "../MT_Datasets/3D_MT_data_edited_fromDuanJM",
    "../MT_Datasets/GA_UA_edited_10s-10000s",
    "data/edifiles2",
]


class TestPlotMTResponse(ImageTestCase):
    pass


def _test_gen(edi_path):
    def default(self):
        edi_file_list = glob.glob(os.path.join(edi_path, "*.edi"))
        for edi_file in edi_file_list:
            plt.clf()
            mt_obj = mt.MT(edi_file)
            pt_obj = mt_obj.plot_mt_response(plot_yn="n")
            pt_obj.plot()
            plt.pause(0.5)
            save_figure_name = "{}.png".format(default.__name__)
            save_figure_path = os.path.join(self._temp_dir, save_figure_name)
            pt_obj.save_plot(save_figure_path)
            assert os.path.isfile(save_figure_path)

    return default


for edi_path in edi_paths:
    if os.path.isdir(edi_path):
        test_name = os.path.basename(edi_path)
        _test_func = _test_gen(edi_path)
        _test_func.__name__ = "test_{test_name}_{plot_name}".format(
            test_name=test_name, plot_name=_test_func.__name__
        )
        setattr(TestPlotMTResponse, _test_func.__name__, _test_func)
