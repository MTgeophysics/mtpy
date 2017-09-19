import glob
import inspect
import os
import unittest
from unittest import TestCase

import matplotlib.pyplot as plt

plt.ion()

from mtpy.imaging.phase_tensor_pseudosection import PlotPhaseTensorPseudoSection

edi_paths = [
    "",
    "tests/data/edifiles",
    "examples/data/edi2",
    "examples/data/edi_files",
    "../MT_Datasets/3D_MT_data_edited_fromDuanJM/",
    "../MT_Datasets/GA_UA_edited_10s-10000s/",
    "tests/data/edifiles2"
]

class TestPlotPhaseTensorPseudoSection(TestCase):
    @classmethod
    def setUpClass(cls):
        cls._temp_dir = "tests/temp"
        if not os.path.isdir(cls._temp_dir):
            os.mkdir(cls._temp_dir)

    @classmethod
    def tearDownClass(cls):
        plt.close('all')

    def setUp(self):
        plt.clf()

    def test_plot_01(self):
        edi_path = edi_paths[1]
        self._plot(edi_path,
                   "%s.png" % inspect.currentframe().f_code.co_name)

    @unittest.skipUnless(os.path.isdir(edi_paths[2]), "data file not found")
    def test_plot_02(self):
        edi_path = edi_paths[2]
        self._plot(edi_path,
                       "%s.png" % inspect.currentframe().f_code.co_name)

    @unittest.skipUnless(os.path.isdir(edi_paths[3]), "data file not found")
    def test_plot_03(self):
        edi_path = edi_paths[3]
        self._plot(edi_path,
                   "%s.png" % inspect.currentframe().f_code.co_name)

    @unittest.skipUnless(os.path.isdir(edi_paths[4]), "data file not found")
    def test_plot_04(self):
        edi_path = edi_paths[4]
        self._plot(edi_path,
                   "%s.png" % inspect.currentframe().f_code.co_name)

    @unittest.skipUnless(os.path.isdir(edi_paths[5]), "data file not found")
    def test_plot_05(self):
        edi_path = edi_paths[5]
        self._plot(edi_path,
                   "%s.png" % inspect.currentframe().f_code.co_name)

    def _plot(self, edi_path, save_figure_path):
        edi_file_list = glob.glob(os.path.join(edi_path, "*.edi"))
        freq = 1
        save_figure_path = os.path.join(self._temp_dir, save_figure_path)
        ptpObj = PlotPhaseTensorPseudoSection(fn_list=edi_file_list,
                                              tscale='period',
                                              # ylim=(1e-1, 1e3),  # orig period
                                              # range to plot
                                              # period range to plot
                                              ylim=(0, 10000),
                                              # xlim = (0,10000),
                                              stretch=(2000, 40),
                                              # determines (x,y) aspect ratio of plot
                                              station_id=(
                                                  0, 10),  # indices for showing station names
                                              ellipse_dict={'size': 6},
                                              plot_tipper='yri',
                                              arrow_dict={'size': 5, 'head_length': 0.2,
                                                          'head_width': 0.1, 'lw': 0.5},
                                              # arrow parameters, adjust as
                                              # necessary. lw = linewidth
                                              font_size=4,
                                              dpi=300)
        ptpObj.plot()
        plt.pause(1)
        ptpObj.save_figure2(save_fn=save_figure_path)
