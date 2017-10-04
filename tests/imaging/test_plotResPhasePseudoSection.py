import glob
import inspect
import os
import unittest

# configure matplotlib for testing
import matplotlib.pyplot as plt

from mtpy.imaging.plotpseudosection import PlotResPhasePseudoSection
from tests.imaging import ImageTestCase

edi_paths = [
    "",
    "tests/data/edifiles",
    "examples/data/edi2",
    "examples/data/edi_files",
    "../MT_Datasets/3D_MT_data_edited_fromDuanJM/",
    "../MT_Datasets/GA_UA_edited_10s-10000s/",
    "tests/data/edifiles2"
]


class TestPlotResPhasePseudoSection(ImageTestCase):
    # def setUp(self):
    #     plt.clf()

    def tearDown(self):
        plt.pause(1)
        plt.close()
        plt.clf()

    @unittest.skipUnless(os.path.isdir(edi_paths[1]), "data file not found")
    def test_plot_01_imshow(self):
        edi_path = edi_paths[1]
        self._plot_imshow(edi_path, "%s.png" % inspect.currentframe().f_code.co_name)

    @unittest.skipUnless(os.path.isdir(edi_paths[1]), "data file not found")
    def test_plot_01_pcolormesh(self):
        edi_path = edi_paths[1]
        self._plot_pcolormesh(edi_path, "%s.png" % inspect.currentframe().f_code.co_name)

    @unittest.skipUnless(os.path.isdir(edi_paths[2]), "data file not found")
    def test_plot_02_imshow(self):
        edi_path = edi_paths[2]
        self._plot_imshow(edi_path,
                          "%s.png" % inspect.currentframe().f_code.co_name)

    @unittest.skipUnless(os.path.isdir(edi_paths[2]), "data file not found")
    def test_plot_02_pcolormesh(self):
        edi_path = edi_paths[2]
        self._plot_pcolormesh(edi_path,
                              "%s.png" % inspect.currentframe().f_code.co_name)

    @unittest.skipUnless(os.path.isdir(edi_paths[3]), "data file not found")
    def test_plot_03_imshow(self):
        edi_path = edi_paths[3]
        self._plot_imshow(edi_path,
                          "%s.png" % inspect.currentframe().f_code.co_name)

    @unittest.skipUnless(os.path.isdir(edi_paths[3]), "data file not found")
    def test_plot_03_pcolormesh(self):
        edi_path = edi_paths[3]
        self._plot_pcolormesh(edi_path,
                              "%s.png" % inspect.currentframe().f_code.co_name)

    @unittest.skipUnless(os.path.isdir(edi_paths[4]), "data file not found")
    def test_plot_04_imshow(self):
        edi_path = edi_paths[4]
        self._plot_imshow(edi_path,
                          "%s.png" % inspect.currentframe().f_code.co_name)

    @unittest.skipUnless(os.path.isdir(edi_paths[4]), "data file not found")
    def test_plot_04_pcolormesh(self):
        edi_path = edi_paths[4]
        self._plot_pcolormesh(edi_path,
                          "%s.png" % inspect.currentframe().f_code.co_name)

    @unittest.skipUnless(os.path.isdir(edi_paths[5]), "data file not found")
    def test_plot_05_imshow(self):
        edi_path = edi_paths[5]
        self._plot_imshow(edi_path,
                          "%s.png" % inspect.currentframe().f_code.co_name)

    @unittest.skipUnless(os.path.isdir(edi_paths[5]), "data file not found")
    def test_plot_05_pcolormesh(self):
        edi_path = edi_paths[5]
        self._plot_pcolormesh(edi_path,
                          "%s.png" % inspect.currentframe().f_code.co_name)

    def _plot_imshow(self, edi_path, save_figure_path):
        edi_file_list = glob.glob(os.path.join(edi_path, "*.edi"))
        save_figure_path = os.path.join(self._temp_dir, save_figure_path)
        pt_obj = PlotResPhasePseudoSection(fn_list=edi_file_list, plot_yn='n', plot_style='imshow')
        pt_obj.plot()
        pt_obj.save_plot(save_figure_path, close_plot='n')
        assert(os.path.isfile(save_figure_path))

    def _plot_pcolormesh(self, edi_path, save_figure_path):
        edi_file_list = glob.glob(os.path.join(edi_path, "*.edi"))
        save_figure_path = os.path.join(self._temp_dir, save_figure_path)
        pt_obj = PlotResPhasePseudoSection(fn_list=edi_file_list, plot_yn='n', plot_style='pcolormesh')
        pt_obj.plot()
        pt_obj.save_plot(save_figure_path, close_plot='n')
        assert(os.path.isfile(save_figure_path))
