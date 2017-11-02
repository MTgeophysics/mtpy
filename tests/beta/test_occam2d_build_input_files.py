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
import mtpy.modeling.occam2d as occam2d

import tests.util_functions as ufun
from unittest import TestCase

import matplotlib.pyplot as plt
plt.ion() # make figure disappear automatically:
#plt.ioff()  # make figure show normally and need to click to close the figure to continue the proc

class TestOccam2D(TestCase):
    def setUp(self):

        # set the dir to the output from the previously correct run
        self._expected_output_dir = r'E:/Githubz/mtpy/examples/model_files/Occam2d'

        if not os.path.isdir(self._expected_output_dir):
            self._expected_output_dir = None

        # directory to save created input files
        self._output_dir = r'E:/Githubz/mtpy/tests/beta/Occam2d'
        # ufun.clean_recreate(self._output_dir) # this may remove other test functions' output
        if not os.path.exists(self._output_dir):
            os.mkdir(self._output_dir)

    def _main_func(self, edipath):
        """
        test function should be successful with a edipath
        :return:
        """
        # path to save to
        savepath = self._output_dir

        # list of stations
        slst=[edi[0:-4] for edi in os.listdir(edipath) if edi.find('.edi')>0]


        # create an occam data object
        ocd = occam2d.Data(edi_path=edipath,
                           station_list=slst,
        #                  interpolate_freq=True,
        #                  freq=np.logspace(-3,1,30)
                           )

        ocd.save_path = savepath

        # choose frequency range to invert

        #ocd.freq_num = 50
        ocd.freq_min = 1
        ocd.freq_max = 10000
        #ocd.freq_num = 50 # number of frequencies to invert for
        ###########make data file

        # error floors
        ocd.res_te_err = 10
        ocd.res_tm_err = 10
        ocd.phase_te_err = 5
        ocd.phase_tm_err = 5
        #ocd.model_mode= 4
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
        ocr.trigger= 0.25

        # z1 layer and target depth in metres
        ocr.z1_layer = 20
        ocr.z_target_depth = 10000
        ocr.num_z_pad_cells = 10
        ocr.z_bottom = 100000
        ocr.save_path=ocd.save_path
        ocr.build_mesh()
        ocr.build_regularization()
        ocr.write_mesh_file()
        ocr.write_regularization_file()

        ocr.plot_mesh()

        #make startup file
        ocs=occam2d.Startup()
        ocs.iterations_to_run=40
        ocs.data_fn=os.path.join(ocd.save_path,'OccamDataFile.dat')
        ocs.resistivity_start=2.0
        ocr.get_num_free_params()
        ocs.param_count=ocr.num_free_param
        ocs.save_path=ocd.save_path
        ocs.model_fn=ocr.reg_fn
        ocs.write_startup_file()

        return savepath

    def test_fun(self):
        """

        :return:
        """
        outdir = self._main_func(edipath = r"E:/Githubz/mtpy/examples/data/edi_files")

        for afile in ('Occam2DMesh',    'Occam2DModel',    'Occam2DStartup'):

            output_data_file =  os.path.join(outdir, afile)
            self.assertTrue(os.path.isfile(output_data_file), "output data file not found")

            expected_data_file = os.path.join(self._expected_output_dir, afile)

            self.assertTrue(os.path.isfile(expected_data_file),
                            "Ref output data file does not exist, nothing to compare with"
                            )

            print ("Comparing", output_data_file, "and", expected_data_file)

            count = ufun.diffiles(output_data_file, expected_data_file)
            if afile == "Occam2DStartup":
                self.assertTrue(count == 1, "Only-1 different line in for this file %s" % afile)
            else:
                self.assertTrue(count == 0, "The output files different in %s lines" % count)
