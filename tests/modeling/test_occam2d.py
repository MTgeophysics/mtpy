# -*- coding: utf-8 -*-
# ! /usr/bin/env python
"""
Description:
    Test create input files required to run occam2d inversion, which include:
    ('Occam2DMesh',    'Occam2DModel',    'Occam2DStartup')
    These input files are created from standard edi data file set.

References:
    examples/tests/occam2d_buildinputfiles.py

CreationDate:   31/10/2017
Developer:      fei.zhang@ga.gov.au

Revision History:
    LastUpdate:  31/10/2017   FZ
"""

import os
from unittest import TestCase

import mtpy.modeling.occam2d as occam2d
from tests import EDI_DATA_DIR, SAMPLE_DIR, make_temp_dir

# import matplotlib.pyplot as plt
# plt.ion() # make figure disappear automatically:
# plt.ioff()  # make figure show normally and need to click to close the figure to continue the proc
from tests.imaging import reset_matplotlib, plt_wait, plt_close
from tests.modeling import diff_files


class TestOccam2D(TestCase):
    @classmethod
    def setUpClass(cls):
        reset_matplotlib()
        cls._temp_dir = make_temp_dir(cls.__name__)

    def setUp(self):

        # set the dir to the output from the previously correct run
        self._expected_output_dir = os.path.join(SAMPLE_DIR, "Occam2d")

        if not os.path.isdir(self._expected_output_dir):
            self._expected_output_dir = None

        # directory to save created input files
        self._output_dir = make_temp_dir("Occam2d", self._temp_dir)

    def _main_func(self, edipath):
        """
        test function should be successful with a edipath
        :return:
        """
        # path to save to
        savepath = self._output_dir

        # list of stations
        slst = [edi[0:-4] for edi in os.listdir(edipath) if edi.find(".edi") > 0]

        # create an occam data object
        ocd = occam2d.Data(
            edi_path=edipath,
            station_list=slst,
            #                  interpolate_freq=True,
            #                  freq=np.logspace(-3,1,30)
        )

        ocd.save_path = savepath

        # choose frequency range to invert

        # ocd.freq_num = 50
        ocd.freq_min = 1
        ocd.freq_max = 10000
        # ocd.freq_num = 50 # number of frequencies to invert for
        ###########make data file

        # error floors
        ocd.res_te_err = 10
        ocd.res_tm_err = 10
        ocd.phase_te_err = 5
        ocd.phase_tm_err = 5
        # ocd.model_mode= 4
        ocd.write_data_file()

        # make model and mesh files
        ocr = occam2d.Regularization(ocd.station_locations)
        # number of layers
        ocr.n_layers = 60
        # cell width to aim for, note this is the mesh size (2 mesh blocks per model block)
        ocr.cell_width = 200
        # controls number and size of padding
        ocr.num_x_pad_cells = 9
        ocr.x_pad_multiplier = 1.9
        # controls aspect ratio of blocks
        ocr.trigger = 0.25

        # z1 layer and target depth in metres
        ocr.z1_layer = 20
        ocr.z_target_depth = 10000
        ocr.num_z_pad_cells = 10
        ocr.z_bottom = 100000
        ocr.save_path = ocd.save_path
        ocr.build_mesh()
        ocr.build_regularization()
        ocr.write_mesh_file()
        ocr.write_regularization_file()

        ocr.plot_mesh()
        plt_wait(1)
        plt_close()

        # make startup file
        ocs = occam2d.Startup()
        ocs.iterations_to_run = 40
        ocs.data_fn = os.path.join(ocd.save_path, "OccamDataFile.dat")
        ocs.resistivity_start = 2.0
        ocr.get_num_free_params()
        ocs.param_count = ocr.num_free_param
        ocs.save_path = ocd.save_path
        ocs.model_fn = ocr.reg_fn
        ocs.write_startup_file()

        return savepath

    def test_fun(self):
        """

        :return:
        """

        outdir = self._main_func(edipath=EDI_DATA_DIR)

        for afile in ("Occam2DMesh", "Occam2DModel", "Occam2DStartup"):
            output_data_file = os.path.join(outdir, afile)
            self.assertTrue(
                os.path.isfile(output_data_file), "output data file not found"
            )

            expected_data_file = os.path.join(self._expected_output_dir, afile)

            self.assertTrue(
                os.path.isfile(expected_data_file),
                "Ref output data file does not exist, nothing to compare with",
            )

            print(("Comparing", output_data_file, "and", expected_data_file))

            is_identical, msg = diff_files(
                output_data_file, expected_data_file, ignores=["Date/Time:"]
            )
            print(msg)
            self.assertTrue(
                is_identical, "The output file is not the same with the baseline file."
            )

    def test_plot_model_and_responses(self):
        """
            test function
            :return: T/F
            """

        # path to directory containing inversion files
        idir = os.path.join(SAMPLE_DIR, "Occam2d")

        # save path, to save plots to
        savepath = self._temp_dir
        offset = 0

        # go to model results directory and find the latest iteration file
        iterfile = "ITER12.iter"
        respfile = "RESP12.resp"

        datafn = "OccamDataFile.dat"
        # get the iteration number
        iterno = iterfile[-7:-5]
        outfilename = iterfile[:-5]

        plotmodel = True  # set to True to plot the resistivity model
        plotresponses = True  # set to True to plot the responses
        save = True

        # horizontal padding on the edges for plotting, in km
        xpad = 1

        # plot the model
        if plotmodel:
            plotm = occam2d.PlotModel(
                iter_fn=os.path.join(idir, iterfile),
                data_fn=os.path.join(idir, datafn),
                station_font_pad=0.5,
                station_font_size=6,
                station_font_rotation=75,
                climits=(0.0, 2.5),  # colour scale limits
                xpad=xpad,
                dpi=300,  # resolution of figure
                fig_aspect=0.5,  # aspect ratio between horizontal and vertical scale
                ylimits=(0, 10),  # depth limits
                stationid=(-1, 3),  # index of station name to plot
                plot_yn="n",
            )
            plotm.plot()
            if save:
                plotm.save_figure(
                    os.path.join(savepath, outfilename + "_resmodel.png"), close_fig="n"
                )  # this will produce 1 figure .png
            plt_wait(1)

        # plot the responses
        if plotresponses:
            plotresponse = occam2d.PlotResponse(
                os.path.join(idir, datafn),
                resp_fn=os.path.join(idir, respfile),
                plot_type=["pb35", "pb40"],
            )
            if save:
                plotresponse.save_figures(
                    savepath, close_fig="n"
                )  # this will produce 2 .pdf file
            plt_wait(1)

        plt_close()
