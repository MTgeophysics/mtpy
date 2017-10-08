import glob
import inspect
import os
import unittest
from unittest import TestCase

# configure matplotlib for testing
import matplotlib.pyplot as plt
import pytest

from mtpy.imaging.plotstrike import PlotStrike
from mtpy.utils.decorator import ImageCompare
from tests.imaging import ImageTestCase

plt.ion()

edi_paths = [
    "tests/data/edifiles",
    "examples/data/edi2",
    "examples/data/edi_files",
    "../MT_Datasets/3D_MT_data_edited_fromDuanJM",
    "../MT_Datasets/GA_UA_edited_10s-10000s",
    "tests/data/edifiles2"
]


class TestPlotStrike(ImageTestCase):
    def setUp(self):
        plt.clf()


def test_gen(edi_path):
    def default(self):
        edi_file_list = glob.glob(os.path.join(edi_path, "*.edi"))
        pt_obj = PlotStrike(fn_list=edi_file_list, plot_yn='n', save_figure_path=self._temp_dir)
        pt_obj.plot()
        plt.pause(1)
        save_figure_name = "{}.png".format(default.__name__)
        save_figure_path = os.path.join(self._temp_dir, save_figure_name)
        pt_obj.save_plot(save_figure_path, file_format='png', close_plot='n')
        assert (os.path.isfile(save_figure_path))

    def rotation(self):
        edi_file_list = glob.glob(os.path.join(edi_path, "*.edi"))
        # change rotation
        pt_obj = PlotStrike(fn_list=edi_file_list, plot_yn='n', rot_z=90)
        pt_obj.plot()
        plt.pause(1)

    def type(self):
        edi_file_list = glob.glob(os.path.join(edi_path, "*.edi"))
        # plot type
        pt_obj = PlotStrike(fn_list=edi_file_list, plot_yn='n', plot_type=1)
        pt_obj.plot()
        plt.pause(1)

    def tipper(self):
        edi_file_list = glob.glob(os.path.join(edi_path, "*.edi"))
        # plot_tipper
        pt_obj = PlotStrike(fn_list=edi_file_list, plot_yn='n', plot_tipper='y')
        pt_obj.plot()
        plt.pause(1)

    def fold(self):
        edi_file_list = glob.glob(os.path.join(edi_path, "*.edi"))
        # fold
        pt_obj = PlotStrike(fn_list=edi_file_list, plot_yn='n', fold=False)
        pt_obj.plot()
        plt.pause(1)

    return default, rotation, type, tipper, fold


# generate tests
for edi_path in edi_paths:
    if os.path.isdir(edi_path):
        test_name = os.path.basename(edi_path)
        for test_func in test_gen(edi_path):
            test_func.__name__ = "test_{test_name}_{plot_name}".format(
                test_name=test_name, plot_name=test_func.__name__)
            setattr(
                TestPlotStrike,
                test_func.__name__,
                ImageCompare(fig_size=(8, 6), savefig_kwargs={"dpi": 100}).__call__(test_func))

if 'test_gen' in globals():
    del globals()['test_gen']
if 'test_func' in globals():
    del globals()['test_func']
