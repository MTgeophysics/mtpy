#! /usr/bin/env python
"""
Description:

    This script is to create input files required to run occam1d inversion, which include:

    OccamStartup1D
    Model1D
    Occam1d_DataFile_DET.dat.

    These input files are created from a standard edi data file.
    
References: 
    examples/tests/occam1d_buildinputfiles.py

CreationDate:   31/10/2017
Developer:      fei.zhang@ga.gov.au

Revision History:
    LastUpdate:     31/10/2017   FZ
    LastUpdate:     dd/mm/yyyy  Who     Optional description
"""

# import section

import os
import sys
import difflib
import mtpy.modeling.occam1d as mtoc1d  # Wrapper class to interact with Occam1D
import shutil
from unittest import TestCase

class TestOccam1D(TestCase):

    def setUp(self):
        # for each test, setup a different output dir
        self._output_dir = r'E:/Githubz/mtpy/tests/Occam1d'
        if os.path.exists(self._output_dir):
        # clear dir if it already exist
            shutil.rmtree(self._output_dir)

        os.mkdir(self._output_dir)

        # set the dir to the output from the previously correct run
        self._expected_output_dir = r'E:/Githubz/mtpy/examples/model_files/Occam1d'

        if not os.path.isdir(self._expected_output_dir):
            self._expected_output_dir = None


    def test_all(self, path2edifile = r'E:/Githubz/mtpy/examples/data/edi_files/pb23c.edi',
                 savepath = r'E:/Githubz/mtpy/tests/Occam1d'):
        """
        test function
        :return:
        """

        # directory to save created input files
        # savepath = r'C:\Git\mtpy\examples\model_files\Occam1d'
        # savepath = r'E:/Githubz/mtpy/examples/model_files/Occam1d'

        if not os.path.exists(savepath):
            os.mkdir(savepath)

        # create data file
        ocd = mtoc1d.Data()  # create an object and assign values to arguments

        ocd.write_data_file(edi_file=path2edifile,
                            mode='det',
                            # mode, can be te, tm, det (for res/phase) or tez, tmz, zdet for real/imag impedance tensor values
                            save_path=savepath,
                            res_errorfloor=5,  # percent error floor
                            phase_errorfloor=1,  # error floor in degrees
                            z_errorfloor=2.5,
                            remove_outofquadrant=True)

        # create model file
        ocm = mtoc1d.Model(n_layers=100,  # number of layers
                           target_depth=10000,  # target depth in metres, before padding
                           z1_layer=10  # first layer thickness in metres
                           )
        ocm.write_model_file(save_path=savepath)

        # create startup file
        ocs = mtoc1d.Startup(data_fn=ocd.data_fn,  # basename of data file *default* is Occam1DDataFile
                             model_fn=ocm.model_fn,  # basename for model file *default* is Model1D
                             max_iter=200,  # maximum number of iterations to run
                             target_rms=0.0)

        ocs.write_startup_file()

    # Wrote Data File to : E:/Githubz/mtpy/tests/Occam1d\Occam1d_DataFile_DET.dat
    # Wrote Model file: E:/Githubz/mtpy/tests/Occam1d\Model1D
    # Wrote Input File: E:/Githubz/mtpy/tests/Occam1d\OccamStartup1D

        output_data_file = os.path.normpath(os.path.join(savepath, "Occam1d_DataFile_DET.dat"))
        self.assertTrue(os.path.isfile(output_data_file), "output data file does not exist")
        expected_data_file = os.path.normpath(os.path.join(self._expected_output_dir,
                                                           "Occam1d_DataFile_DET.dat"))
        self.assertTrue(
            os.path.isfile(expected_data_file),
            "expected output data file does not exist, nothing to compare"
        )

        print ("comparing", output_data_file, "and", expected_data_file)

        # with open(output_data_file, 'r') as output:
        #     with open(expected_data_file, 'r') as expected:
        #         diff = difflib.unified_diff(
        #             expected.readlines(),
        #             output.readlines(),
        #             fromfile='expected',
        #             tofile='output'
        #         )
        # count = 0
        # for line in diff:
        #     sys.stdout.write(line)
        #     count += 1
        # self.assertTrue(count == 0, "output different!")


