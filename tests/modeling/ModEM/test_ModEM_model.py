from unittest import TestCase
import numpy as np
import pytest
import os

from mtpy.modeling.modem import Model, Data
from tests import make_temp_dir, SAMPLE_DIR
from tests.modeling import diff_files

class TestModEM_Model(TestCase):
    def setUp(self):
        self.model_epsg = 28355
        self._temp_dir = make_temp_dir(self.__name__)
        # directory to save created input files
        self._output_dir = make_temp_dir(self._testMethodName, base_dir=self._temp_dir)

        
        
    def test_read_gocad_sgrid_file(self):
        
        # set the dir to the output from the previously correct run
        expected_output_dir = os.path.join(SAMPLE_DIR, 'ModEM')        
        
        self._sgrid_file = os.path.join(SAMPLE_DIR,'gocad','ModEM_Model_File.sg')
        
        if not os.path.isdir(expected_output_dir):
            expected_output_dir = None
            
        output_fn = 'ModEM_Model_File.rho'

        # read data file to get centre position
        data_fn = os.path.join(SAMPLE_DIR,'ModEM','ModEM_Data.dat')
        dObj = Data()
        dObj.read_data_file(data_fn=data_fn)
        
        mObj = Model(data_obj = dObj,save_path=self._output_dir)
        mObj.read_gocad_sgrid_file(self._sgrid_file)
        mObj.write_model_file()
            
        output_data_file = os.path.normpath(os.path.join(self._output_dir, output_fn))

        self.assertTrue(os.path.isfile(output_data_file), "output data file not found")

        expected_data_file = os.path.normpath(os.path.join(expected_output_dir, 'ModEM_Model_File.rho'))

        self.assertTrue(os.path.isfile(expected_data_file),
                        "Ref output data file does not exist, nothing to compare with"
                        )

        # print ("Comparing", output_data_file, "and", expected_data_file)

        is_identical, msg = diff_files(output_data_file, expected_data_file)
        print msg
        self.assertTrue(is_identical, "The output file is not the same with the baseline file.")
        
    
    def test_write_gocad_sgrid_file(self):
        return