"""
Test the classes of the module modem.py
"""
import glob
import os
import shutil
import tempfile

import unittest
from unittest import TestCase
from mtpy.modeling.modem import Data
from mtpy.modeling.modem import Model
from mtpy.modeling.modem import Covariance

from mtpy.core.edi_collection import EdiCollection

import numpy as np

def select_periods(edifiles_list):
    """
    FZ: Use edi_collection to analyse the whole set of EDI files
    :param edifiles:
    :return:
    """

    edis_obj = EdiCollection(edifiles_list)

    uniq_period_list = edis_obj.all_unique_periods  # filtered list of periods ?
    print("Unique periods",len(uniq_period_list))

    #1 ASK user to input a Pmin and Pmax

    #2 percetage stats
    # select commonly occured frequencies from all stations.
    # This could miss some slightly varied frquencies in the middle range.
    select_period_list = np.array(edis_obj.get_periods_by_stats(percentage=10.0))
    print("Selected periods ", len(select_period_list))

    return select_period_list


class TestModem(TestCase):

    def setUp(self):

        print ("Calling setUp")

        epsg_code = 3112

        self.inputdir = '../examples/data/edi2'
        self.inputdir = 'E:/Data/MT_Datasets/3D_MT_data_edited_fromDuanJM'
        self.inputdir = 'E:/Data/MT_Datasets/GA_UA_edited_10s-10000s'

        self.topofile = 'E:/Data/MT_Datasets/aussie_etopo1_bedrock.asc'
        # self.topofile = '../examples/etopo1.asc'

        tempdir='E:/temp'  # tempfile.gettempdir()
        self.outputdir = os.path.join(tempdir,'test_out') # '/tmp/test_out/'

        #clean-up the output dir then creat a new one (empty)
        if os.path.exists(self.outputdir):
            shutil.rmtree(self.outputdir)
        os.mkdir(self.outputdir)

        edi_list = glob.glob(self.inputdir + '/*.edi')

        assert (edi_list)>1

        # period list (can take periods from one of the edi files, or just specify
        # periods directly using the logspace function (commented out))

        # eo = mtedi.Edi(edi_list[0])  # this may miss some periods?
        # period_list = 1. / eo.Z.freq # period_list = np.logspace(-3,3)

        period_list = select_periods(edi_list)

        self.datob = Data(edi_list=edi_list,
                  inv_mode='1',
                  period_list=period_list,
                  epsg=epsg_code,
                  error_type='floor',
                  error_floor=10)
        # period_buffer=0.000001)

        self.datob.write_data_file(save_path=self.outputdir)

        # create model file
        self.model = Model(Data=self.datob,
               # cell_size_east=2000, cell_size_north=2000,
               cell_size_east=10000, cell_size_north=10000,
               pad_north=10,  # number of padding cells in each of the north and south directions
               pad_east=10,  # number of east and west padding cells
               pad_z=10,  # number of vertical padding cells
               pad_stretch_v=3,
               # factor to increase by in padding cells (vertical)
               pad_stretch_h=3,
               # factor to increase by in padding cells (vertical)
               n_airlayers=10,  # number of air layers
               res_model=100,  # halfspace resistivity value for reference model
               n_layers=40,  # total number of z layers, including air
               z1_layer=1000,  # first layer thickness
               epsg=epsg_code,  # epsg
               z_target_depth=200000)

        self.model.make_mesh()

        #model.plot_mesh()

        # write a model file to initialise a resistivity model
        self.model.write_model_file(save_path=self.outputdir)

        # add topography to res model
        self.model.add_topography(self.topofile, interp_method='nearest')

        # make covariance file
        self.cov = Covariance(mask_arr=self.model.covariance_mask,
              save_path=self.outputdir,
              smoothing_east=0.3,
              smoothing_north=0.3,
              smoothing_z=0.3)

        self.cov.write_covariance_file(model_fn=self.model.model_fn)

    def tearDown(self):
        print ("Calling tearDown -- to clean up the results/objects of this test run")
        del self.model
        del self.datob
        del self.cov

        print ("please check %s" % self.outputdir)


    # def test_something(self):
    #     # self.assertEqual(True, False)
    #     self.assertEqual(True, 1) # 0=False 1=True

    def test_topo_ascii(self):
        """

        :return:
        """
        print ("testing topo file reading...")

        # check results

        return True


if __name__ == '__main__':

    unittest.main()
