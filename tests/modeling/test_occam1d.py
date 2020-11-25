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
    LastUpdate:     07/10/2017   YG
                    31/10/2017   FZ



Current issue:
    test func2  failed when comparing output files (RhoZxy value different)
"""

# import section

import os

import numpy as np

from unittest import TestCase

import mtpy.modeling.occam1d as mtoc1d  # Wrapper class to interact with Occam1D
import tests.imaging
import tests.modeling
from tests import EDI_DATA_DIR, SAMPLE_DIR, make_temp_dir
from tests.imaging import reset_matplotlib


class TestOccam1D(TestCase):
    @classmethod
    def setUpClass(cls):
        reset_matplotlib()
        cls._temp_dir = make_temp_dir(cls.__name__)

    def setUp(self):
        # set the dir to the output from the previously correct run
        self._expected_output_dir = os.path.join(SAMPLE_DIR, "Occam1d")

        if not os.path.isdir(self._expected_output_dir):
            self._expected_output_dir = None

        # directory to save created input files
        self._output_dir = make_temp_dir(self._testMethodName, base_dir=self._temp_dir)

    def _main_func(self, path2edifile):
        """
        test function should be successful with a default path2edifile
        :return:
        """
        edifile_name = os.path.basename(path2edifile)
        tmpdir = edifile_name[:-4] + "_dir"  # remove the trailing .edi
        tmp_save_path = os.path.join(self._output_dir, tmpdir)
        tests.modeling._clean_recreate(tmp_save_path)

        # create data file
        ocd = mtoc1d.Data()  # create an object and assign values to arguments

        ocd.write_data_file(
            edi_file=path2edifile,
            mode="det",
            # mode, can be te, tm, det (for res/phase) or tez, tmz, zdet for real/imag impedance tensor values
            save_path=tmp_save_path,
            res_errorfloor=5,  # percent error floor
            phase_errorfloor=1,  # error floor in degrees
            z_errorfloor=2.5,
            remove_outofquadrant=True,
        )

        # create model file
        ocm = mtoc1d.Model(
            n_layers=100,  # number of layers
            target_depth=10000,  # target depth in metres, before padding
            z1_layer=10,  # first layer thickness in metres
        )
        ocm.write_model_file(save_path=tmp_save_path)

        # create startup file
        ocs = mtoc1d.Startup(
            data_fn=ocd.data_fn,  # basename of data file *default* is Occam1DDataFile
            model_fn=ocm.model_fn,  # basename for model file *default* is Model1D
            max_iter=200,  # maximum number of iterations to run
            target_rms=0.0,
        )

        ocs.write_startup_file()

        return tmp_save_path

    def test_fun1(self):
        """ use the same pb23c.edi to reproduce previous run results"""

        outdir = self._main_func(os.path.join(EDI_DATA_DIR, "pb23c.edi"))

        for afile in ("Model1D", "Occam1d_DataFile_DET.dat", "OccamStartup1D"):
            output_data_file = os.path.join(outdir, afile)
            self.assertTrue(
                os.path.isfile(output_data_file), "output data file not found"
            )

            expected_data_file = os.path.join(self._expected_output_dir, afile)

            self.assertTrue(
                os.path.isfile(expected_data_file),
                "Ref output data file does not exist, nothing to compare with",
            )
            is_identical, msg = tests.modeling.diff_files(
                output_data_file, expected_data_file, ignores=["Date/Time"]
            )

            print(msg)
            self.assertTrue(
                is_identical, "The output file is not the same with the baseline file."
            )

    def test_fun2(self):
        """ another test edi case: The output files should be different !!!"""

        # outdir = self._main_func(r'E:/Githubz/mtpy/examples/data/edi_files/pb25c.edi')
        outdir = self._main_func(os.path.join(EDI_DATA_DIR, "pb25c.edi"))

        # for afile in ("Model1D", "Occam1d_DataFile_DET.dat", "OccamStartup1D"):
        for afile in [
            "Occam1d_DataFile_DET.dat"
        ]:  # only one file is different, the other 2 files same?

            output_data_file = os.path.join(outdir, afile)
            self.assertTrue(
                os.path.isfile(output_data_file), "output data file not found"
            )

            expected_data_file = os.path.join(self._expected_output_dir, afile)

            self.assertTrue(
                os.path.isfile(expected_data_file),
                "Ref output data file does not exist, nothing to compare with",
            )

            is_identical, msg = tests.modeling.diff_files(
                output_data_file, expected_data_file
            )
            print(msg)
            self.assertFalse(
                is_identical, "The output file is the same with the baseline file."
            )

    def test_view_outputs(self):
        # FZ's workdir
        savepath = self._temp_dir

        # model and data file names
        modelfn = os.path.join(self._expected_output_dir, "Model1D")
        datafn = os.path.join(self._expected_output_dir, "Occam1d_DataFile_DET.dat")

        # list of all the 'smooth' files, exclude the log file
        fn_list = np.array([ff for ff in os.listdir(self._expected_output_dir)])

        # go to model results directory and find the latest iteration file
        # iterfile = 'ITER_135.iter'
        # respfile = 'ITER_135.resp'
        iterfile = "ITER_97.iter"
        respfile = "ITER_97.resp"

        # get maximum iteration
        iterfn = os.path.join(self._expected_output_dir, iterfile)
        respfn = os.path.join(self._expected_output_dir, respfile)

        # read in the model, don't need these lines to view the output but useful if you want to analyse the data
        oc1m = mtoc1d.Model(model_fn=modelfn)
        oc1m.read_iter_file(iterfn)

        # read in the data file
        oc1d = mtoc1d.Data(data_fn=datafn)
        oc1d.read_data_file(data_fn=datafn)
        oc1d.read_resp_file(resp_fn=respfn, data_fn=datafn)

        # plot the model output
        pr = mtoc1d.Plot1DResponse(
            data_te_fn=datafn,
            data_tm_fn=datafn,
            model_fn=modelfn,
            resp_te_fn=respfn,
            iter_te_fn=iterfn,
            depth_limits=(0, 1),
        )
        pr.axm.set_xlim(1e-1, 1e3)
        pr.axr.set_ylim(1, 100)

        p2file = os.path.join(savepath, "occam1dplot.png")
        pr.save_figure(p2file, close_plot="n")

        tests.imaging.plt_wait(1)
        tests.imaging.plt_close()

        assert os.path.exists(p2file)
