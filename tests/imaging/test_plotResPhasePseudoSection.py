import glob
import os

# configure matplotlib for testing
import matplotlib.pyplot as plt

from mtpy.imaging.plotpseudosection import PlotResPhasePseudoSection
from mtpy.utils.decorator import ImageCompare
from tests.imaging import ImageTestCase

edi_paths = [
    "tests/data/edifiles",
    "examples/data/edi2",
    "examples/data/edi_files",
    "../MT_Datasets/3D_MT_data_edited_fromDuanJM",
    "../MT_Datasets/GA_UA_edited_10s-10000s",
    "../MT_Datasets/GA_UA_edited_10s-10000s",
    "tests/data/edifiles2"
]


class TestPlotResPhasePseudoSection(ImageTestCase):
    pass


def test_gen(edi_path):
    def imshow(self):
        edi_file_list = glob.glob(os.path.join(edi_path, "*.edi"))
        pt_obj = PlotResPhasePseudoSection(fn_list=edi_file_list, plot_yn='n', plot_style='imshow')
        pt_obj.plot()
        plt.pause(1)
        save_figure_name = "{}.png".format(imshow.__name__)
        save_figure_path = os.path.join(self._temp_dir, save_figure_name)
        pt_obj.save_plot(save_figure_path, close_plot='n')
        assert (os.path.isfile(save_figure_path))

    def pcolormesh(self):
        edi_file_list = glob.glob(os.path.join(edi_path, "*.edi"))
        pt_obj = PlotResPhasePseudoSection(fn_list=edi_file_list, plot_yn='n', plot_style='pcolormesh')
        pt_obj.plot()
        plt.pause(1)
        save_figure_name = "{}.png".format(imshow.__name__)
        save_figure_path = os.path.join(self._temp_dir, save_figure_name)
        pt_obj.save_plot(save_figure_path, close_plot='n')
        assert (os.path.isfile(save_figure_path))

    return imshow, pcolormesh


# generate tests
for edi_path in edi_paths:
    if os.path.isdir(edi_path):
        test_name = os.path.basename(edi_path)
        for test_func in test_gen(edi_path):
            test_func.__name__ = "test_{test_name}_{plot_name}".format(
                test_name=test_name, plot_name=test_func.__name__)
            setattr(
                TestPlotResPhasePseudoSection,
                test_func.__name__,
                ImageCompare(fig_size=(8, 6)).__call__(test_func))

if 'test_gen' in globals():
    del globals()['test_gen']
if 'test_func' in globals():
    del globals()['test_func']
