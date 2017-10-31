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

        # directory to save created input files
        self._output_dir = r'E:/Githubz/mtpy/tests/beta/Occam1d'
        if os.path.exists(self._output_dir):
        # clear dir if it already exist
            shutil.rmtree(self._output_dir)

        os.mkdir(self._output_dir)

        # set the dir to the output from the previously correct run
        self._expected_output_dir = r'E:/Githubz/mtpy/examples/model_files/Occam1d'

        if not os.path.isdir(self._expected_output_dir):
            self._expected_output_dir = None


    def test_all(self, path2edifile = r'E:/Githubz/mtpy/examples/data/edi_files/pb23c.edi'):
        """
        test function
        :return:
        """

        # create data file
        ocd = mtoc1d.Data()  # create an object and assign values to arguments

        ocd.write_data_file(edi_file=path2edifile,
                            mode='det',
                            # mode, can be te, tm, det (for res/phase) or tez, tmz, zdet for real/imag impedance tensor values
                            save_path=self._output_dir,
                            res_errorfloor=5,  # percent error floor
                            phase_errorfloor=1,  # error floor in degrees
                            z_errorfloor=2.5,
                            remove_outofquadrant=True)

        # create model file
        ocm = mtoc1d.Model(n_layers=100,  # number of layers
                           target_depth=10000,  # target depth in metres, before padding
                           z1_layer=10  # first layer thickness in metres
                           )
        ocm.write_model_file(save_path= self._output_dir)

        # create startup file
        ocs = mtoc1d.Startup(data_fn=ocd.data_fn,  # basename of data file *default* is Occam1DDataFile
                             model_fn=ocm.model_fn,  # basename for model file *default* is Model1D
                             max_iter=200,  # maximum number of iterations to run
                             target_rms=0.0)

        ocs.write_startup_file()


        for afile in ("Model1D", "Occam1d_DataFile_DET.dat",  "OccamStartup1D"):

            output_data_file = os.path.normpath(os.path.join(self._output_dir, afile))
            self.assertTrue(os.path.isfile(output_data_file), "output data file not found")

            expected_data_file = os.path.normpath(os.path.join(self._expected_output_dir, afile))

            self.assertTrue(os.path.isfile(expected_data_file),
                "Ref output data file does not exist, nothing to compare with"
            )

            print ("Comparing", output_data_file, "and", expected_data_file)

            count = diffiles(output_data_file,expected_data_file)
            if afile == "OccamStartup1D":
                self.assertTrue(count==1, "The output files different in %s lines"%count)
            else:
                self.assertTrue(count==0, "The output files different in %s lines"%count)



def diffiles(f1, f2):
    """
    compare two files
    :param f1:
    :param f2:
    :return: the number count of different lines
    """
    test_lines = open(f1).readlines()
    correct_lines = open(f2).readlines()

    count=0
    for test, correct in zip(test_lines, correct_lines):
        if test != correct:
            print ("Diffiles() Failure@: Expected %r; BUT Got %r." % (correct, test))
            count= count+1
        else:
            pass

    return count



