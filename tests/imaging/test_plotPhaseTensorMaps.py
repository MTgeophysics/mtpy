import glob
import inspect
import os.path
import unittest

import pytest

from mtpy.imaging.penetration import load_edi_files
from mtpy.imaging.phase_tensor_maps import PlotPhaseTensorMaps
from tests.imaging import ImageTestCase, ImageCompare


# configure matplotlib for testing
def _expected_compare_fail():
    pytest.xfail(
        "expected the image to be different on different platform, please check the image manually."
    )


test_params = [
    (
        "data/edifiles",
        1,
        {
            "fig_size": (7, 8),
            "savefig_kwargs": {"dpi": 100},
            "on_compare_fail": _expected_compare_fail,
        },
    ),
    (
        "../MT_Datasets/3D_MT_data_edited_fromDuanJM",
        10,
        {"fig_size": (7, 8), "savefig_kwargs": {"dpi": 100}},
    ),
    (
        "../MT_Datasets/GA_UA_edited_10s-10000s",
        0.025,
        {"fig_size": (8, 5), "savefig_kwargs": {"dpi": 150}},
    ),
    (
        "../MT_Datasets/GA_UA_edited_10s-10000s",
        0.01,
        {"fig_size": (8, 7), "savefig_kwargs": {"dpi": 100}},
    ),
    (
        "../MT_Datasets/GA_UA_edited_10s-10000s",
        0.0625,
        {"fig_size": (8, 5), "savefig_kwargs": {"dpi": 150}},
    ),
    (
        "../MT_Datasets/GA_UA_edited_10s-10000s",
        0.0005,
        {"fig_size": (8, 5), "savefig_kwargs": {"dpi": 150}},
    ),
    (
        "data/edifiles2",
        1,
        {
            "fig_size": (7, 8),
            "savefig_kwargs": {"dpi": 100},
            "on_compare_fail": _expected_compare_fail,
        },
    ),
]


class TestPlotPhaseTensorMaps(ImageTestCase):
    @classmethod
    def setUpClass(cls):
        super(TestPlotPhaseTensorMaps, cls).setUpClass()
        # 1) Define plots params
        # parameters describing ellipses, differ for different map scales: deg, m, km
        # Try different size to find a suitable value for your case. as a
        # guidance: 1 degree=100KM
        cls.ellipse_dict = {
            "size": 0.2,
            "colorby": "phimin",
            "range": (0, 90, 1),
            "cmap": "mt_bl2gr2rd",
        }

        # adjust to suitable size: parameters describing the induction vector arrows
        cls.arrow_dict = {
            "size": 0.5,
            "lw": 0.2,
            "head_width": 0.04,
            "head_length": 0.04,
            "threshold": 0.8,
            "direction": 0,
        }

        # parameters describing the arrow legend (not necessarily used)
        # self.arrow_legend_dict = {'position': 'upper right',
        #                      'fontpad': 0.0025,
        #                      'xborderpad': 0.07,
        #                      'yborderpad': 0.015}

    @unittest.skipUnless(os.path.isdir("data/edifiles2"), "data file not found")
    # @unittest.expectedFailure
    @pytest.mark.skip(reason="no way of currently testing this")
    def test_edifiles2_input(self):
        """
        testing to use Z and tipper objects as input

        this fails because the constructor of PlotPhaseTensorMaps only initialize the Mplot object properly when reading from files
        :return:
        """
        edi_path = test_params[4]
        freq = 1
        mt_objs = load_edi_files(edi_path)
        z_objs = [mt.Z for mt in mt_objs]
        tipper = [mt.Tipper for mt in mt_objs]
        save_figure_path = os.path.join(
            self._temp_dir, "%s.png" % inspect.currentframe().f_code.co_name
        )
        save_param_path = os.path.join(
            self._temp_dir, "params_%s" % inspect.currentframe().f_code.co_name
        )
        pt_obj = PlotPhaseTensorMaps(
            z_object_list=z_objs,
            tipper_object_list=tipper,
            plot_freq=freq,
            ftol=0.10,  # freq tolerance,which will decide how many data points included
            mapscale="deg",  # deg or m, or km
            xpad=0.4,  # plot margin; change according to lat-lon in edifiles
            ypad=0.4,  # ~ 2* ellipse size
            # ellipse_dict=self.ellipse_dict, # not implemented
            ellipse_size=0.2,
            ellipse_colorby="phimin",
            ellipse_range=(0, 90, 1),
            ellipse_cmap="mt_bl2gr2rd",
            plot_tipper="yr",
            # arrow_dict=self.arrow_dict, # not implemented
            arrow_size=0.5,
            arrow_lw=0.2,
            arrow_head_width=0.04,
            arrow_head_length=0.04,
            arrow_direction=0,
            arrow_threshold=0.8,
            # arrow_legend_dict=arrow_legend_dict,
            # fig_spython examples/plot_phase_tensor_map.py data/edifiles/ 10 /e/MTPY2_Outputs/ptmap3deg.pngize=(6, 5),
            # fig_dpi=300, the default is OK. Higher dpi
            # may distort figure
            save_fn=save_figure_path,
            fig_size=(8, 6),
            fig_dpi=100,
        )
        path2figure = pt_obj.plot()
        pt_obj.save_figure(save_figure_path)
        assert os.path.isfile(save_figure_path)
        pt_obj.export_params_to_file(save_path=save_param_path)
        assert os.path.isdir(save_param_path)


def _test_gen(edi_path, freq):
    def default(self):
        save_figure_path = os.path.join(self._temp_dir, "%s.png" % default.__name__)
        save_param_path = os.path.join(self._temp_dir, "params_%s" % default.__name__)
        edi_file_list = glob.glob(os.path.join(edi_path, "*.edi"))
        pt_obj = PlotPhaseTensorMaps(
            fn_list=edi_file_list,
            plot_freq=freq,
            ftol=0.10,  # freq tolerance,which will decide how many data points included
            mapscale="deg",  # deg or m, or km
            xpad=0.4,  # plot margin; change according to lat-lon in edifiles
            ypad=0.4,  # ~ 2* ellipse size
            # ellipse_dict=self.ellipse_dict, # Not implemented
            ellipse_size=0.2,
            ellipse_colorby="phimin",
            ellipse_range=(0, 90, 1),
            ellipse_cmap="mt_bl2gr2rd",
            plot_tipper="yr",
            # arrow_dict=self.arrow_dict, # Not implemented
            arrow_size=0.5,
            arrow_lw=0.2,
            arrow_head_width=0.04,
            arrow_head_length=0.04,
            arrow_direction=0,
            arrow_threshold=0.8,
            # arrow_legend_dict=arrow_legend_dict,
            # fig_spython examples/plot_phase_tensor_map.py data/edifiles/ 10 /e/MTPY2_Outputs/ptmap3deg.pngize=(6, 5),
            # fig_dpi=300, the default is OK. Higher dpi
            # may distort figure
            save_fn=save_figure_path,
        )
        # 3) do the plot and save figure - if the param save_path provided
        path2figure = pt_obj.plot(show=True)
        pt_obj.save_figure(save_figure_path, close_plot="n")
        assert os.path.isfile(save_figure_path)
        pt_obj.export_params_to_file(save_path=save_param_path)
        assert os.path.isdir(save_param_path)

    return (default,)


# generate tests
for edi_path, freq, img_kwargs in test_params:
    if os.path.isdir(edi_path):
        test_name = os.path.basename(edi_path)
        for _test_func in _test_gen(edi_path, freq):
            plot_name = _test_func.__name__
            _test_func.__name__ = "test_{test_name}_{freq}_{plot_name}".format(
                test_name=test_name,
                freq=str(freq).replace(".", "_"),
                plot_name=plot_name,
            )
            setattr(
                TestPlotPhaseTensorMaps,
                _test_func.__name__,
                ImageCompare(**img_kwargs).__call__(_test_func),
            )
