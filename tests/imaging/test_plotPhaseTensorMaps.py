import glob
import inspect
import os.path
import unittest
from unittest import TestCase

# configure matplotlib for testing
import matplotlib.pyplot as plt

plt.ion()
from mtpy.imaging.penetration import load_edi_files
from mtpy.imaging.phase_tensor_maps import PlotPhaseTensorMaps

# use non-interactive backend 'Agg', so that you do not have to see figure pop-out.

edi_paths = [
    "",
    "tests\\data\\edifiles",
    "../MT_Datasets/3D_MT_data_edited_fromDuanJM/",
    "../MT_Datasets/GA_UA_edited_10s-10000s/",
    "tests\\data\\edifiles2"
]


class TestPlotPhaseTensorMaps(TestCase):
    @classmethod
    def setUpClass(cls):
        cls._temp_dir = "tests/temp"
        if not os.path.isdir(cls._temp_dir):
            os.mkdir(cls._temp_dir)

        # 1) Define plots params
        # parameters describing ellipses, differ for different map scales: deg, m, km
        # Try different size to find a suitable value for your case. as a
        # guidance: 1 degree=100KM
        cls.ellipse_dict = {
            'size': 0.2,
            'colorby': 'phimin',
            'range': (
                0,
                90,
                1),
            'cmap': 'mt_bl2gr2rd'}

        # adjust to suitable size: parameters describing the induction vector arrows
        cls.arrow_dict = {'size': 0.5,
                          'lw': 0.2,
                          'head_width': 0.04,
                          'head_length': 0.04,
                          'threshold': 0.8,
                          'direction': 0}

        # parameters describing the arrow legend (not necessarily used)
        # self.arrow_legend_dict = {'position': 'upper right',
        #                      'fontpad': 0.0025,
        #                      'xborderpad': 0.07,
        #                      'yborderpad': 0.015}

    @classmethod
    def tearDownClass(cls):
        plt.close('all')

    def setUp(self):
        plt.clf()

    def test_plot_01(self):
        edi_path = edi_paths[1]
        freq = 1
        self._plot_gen(edi_path, freq,
                       "%s.png" % inspect.currentframe().f_code.co_name,
                       "params_%s" % inspect.currentframe().f_code.co_name)

    @unittest.skipUnless(os.path.isdir(edi_paths[2]), "data file not found")
    def test_plot_02(self):
        edi_path = edi_paths[2]
        freq = 10
        self._plot_gen(edi_path, freq,
                       "%s.png" % inspect.currentframe().f_code.co_name,
                       "params_%s" % inspect.currentframe().f_code.co_name)

    @unittest.skipUnless(os.path.isdir(edi_paths[3]), "data file not found")
    def test_plot_03_01(self):
        edi_path = edi_paths[3]
        freq = 0.025
        self._plot_gen(edi_path, freq,
                       "%s.png" % inspect.currentframe().f_code.co_name,
                       "params_%s" % inspect.currentframe().f_code.co_name)

    @unittest.skipUnless(os.path.isdir(edi_paths[3]), "data file not found")
    def test_plot_03_02(self):
        edi_path = edi_paths[3]
        freq = 0.01
        self._plot_gen(edi_path, freq,
                       "%s.png" % inspect.currentframe().f_code.co_name,
                       "params_%s" % inspect.currentframe().f_code.co_name)

    @unittest.skipUnless(os.path.isdir(edi_paths[3]), "data file not found")
    def test_plot_03_03(self):
        edi_path = edi_paths[3]
        freq = 0.0625
        self._plot_gen(edi_path, freq,
                       "%s.png" % inspect.currentframe().f_code.co_name,
                       "params_%s" % inspect.currentframe().f_code.co_name)

    @unittest.skipUnless(os.path.isdir(edi_paths[3]), "data file not found")
    def test_plot_03_04(self):
        edi_path = edi_paths[3]
        freq = 0.0005
        self._plot_gen(edi_path, freq,
                       "%s.png" % inspect.currentframe().f_code.co_name,
                       "params_%s" % inspect.currentframe().f_code.co_name)

    @unittest.skipUnless(os.path.isdir(edi_paths[3]), "data file not found")
    def test_plot_04_01(self):
        edi_path = edi_paths[4]
        freq = 1
        self._plot_gen(edi_path, freq,
                       "%s.png" % inspect.currentframe().f_code.co_name,
                       "params_%s" % inspect.currentframe().f_code.co_name)

    @unittest.skipUnless(os.path.isdir(edi_paths[3]), "data file not found")
    def test_plot_04_02(self):
        """
        testing to use Z and tipper objects as input

        this fails because the constructor of PlotPhaseTensorMaps only initialize the Mplot object properly when reading from files
        :return:
        """
        edi_path = edi_paths[4]
        freq = 1
        mt_objs = load_edi_files(edi_path)
        z_objs = [mt.Z for mt in mt_objs]
        tipper = [mt.Tipper for mt in mt_objs]
        save_figure_path = os.path.join(self._temp_dir, "%s.png" % inspect.currentframe().f_code.co_name)
        save_param_path = os.path.join(self._temp_dir, "params_%s" % inspect.currentframe().f_code.co_name)
        pt_obj = PlotPhaseTensorMaps(z_object_list=z_objs,
                                     tipper_object_list=tipper,
                                     plot_freq=freq,
                                     ftol=0.10,  # freq tolerance,which will decide how many data points included
                                     mapscale='deg',  # deg or m, or km
                                     xpad=0.4,  # plot margin; change according to lat-lon in edifiles
                                     ypad=0.4,  # ~ 2* ellipse size
                                     ellipse_dict=self.ellipse_dict,
                                     plot_tipper='yr',
                                     arrow_dict=self.arrow_dict,
                                     # arrow_legend_dict=arrow_legend_dict,
                                     # fig_spython examples/plot_phase_tensor_map.py tests/data/edifiles/ 10 /e/MTPY2_Outputs/ptmap3deg.pngize=(6, 5),
                                     # fig_dpi=300, the default is OK. Higher dpi
                                     # may distort figure
                                     save_fn=save_figure_path)
        path2figure = pt_obj.plot()
        plt.pause(1)
        pt_obj.save_figure(save_figure_path)
        pt_obj.export_params_to_file(save_path=save_param_path)

    def _plot_gen(self, edi_path, freq, save_figure_path=None, save_param_path=None):
        edi_file_list = glob.glob(os.path.join(edi_path, "*.edi"))
        save_figure_path = os.path.join(self._temp_dir, save_figure_path)
        save_param_path = os.path.join(self._temp_dir, save_param_path)
        pt_obj = PlotPhaseTensorMaps(fn_list=edi_file_list,
                                     plot_freq=freq,
                                     ftol=0.10,  # freq tolerance,which will decide how many data points included
                                     mapscale='deg',  # deg or m, or km
                                     xpad=0.4,  # plot margin; change according to lat-lon in edifiles
                                     ypad=0.4,  # ~ 2* ellipse size
                                     ellipse_dict=self.ellipse_dict,
                                     plot_tipper='yr',
                                     arrow_dict=self.arrow_dict,
                                     # arrow_legend_dict=arrow_legend_dict,
                                     # fig_spython examples/plot_phase_tensor_map.py tests/data/edifiles/ 10 /e/MTPY2_Outputs/ptmap3deg.pngize=(6, 5),
                                     # fig_dpi=300, the default is OK. Higher dpi
                                     # may distort figure
                                     save_fn=save_figure_path)
        # 3) do the plot and save figure - if the param save_path provided
        path2figure = pt_obj.plot(show=True)
        plt.pause(1)
        pt_obj.save_figure(save_figure_path)
        pt_obj.export_params_to_file(save_path=save_param_path)
