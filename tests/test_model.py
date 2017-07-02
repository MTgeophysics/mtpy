import glob
import os
import shutil
from unittest import TestCase

from mtpy.core.edi import Edi
from mtpy.modeling.modem_data import Data
from mtpy.modeling.modem_model import Model
from mtpy.modeling.modem_covariance import Covariance


class TestModemModel(TestCase):
    def setUp(self):
        """
        Set up test environment.
        :return:
        """

        print ("Calling setUp")

        self.epsg_code = 3112

        self.inputdir = 'tests/data/edifiles'  # '../examples/data/edi2'

        # self.inputdir = 'E:/Data/MT_Datasets/3D_MT_data_edited_fromDuanJM'
        # self.inputdir = '/e/Data/MT_Datasets/GA_UA_edited_10s-10000s'

        self.topofile = '/e/Data/MT_Datasets/concurry_topo/aussie_etopo1_bedrock.asc'
        # self.topofile = '../examples/etopo1.asc'

        tempdir = 'temp'  # tempfile.gettempdir()
        self.outputdir = os.path.join(tempdir, 'test_out')  # '/tmp/test_out/'

        # clean-up the output dir then creat a new one (empty)
        if os.path.exists(self.outputdir):
            shutil.rmtree(self.outputdir)

        os.mkdir(self.outputdir)

        print ("please check %s" % self.outputdir)

        # 1) make a list of all .edi files that will be inverted for
        edi_list = glob.glob(self.inputdir + '/*.edi')

        # 2) create data file
        eobj = Edi(edi_list[0])  # this may miss some periods?
        period_list = 1. / eobj.Z.freq # period_list = np.logspace(-3,3)
        md = Data(edi_list=edi_list,
                  inv_mode='1',
                  period_list = period_list,
                  epsg=self.epsg_code,
                  error_type='floor',
                  error_floor=10)
        # period_buffer=0.000001)

        md.write_data_file(save_path=self.outputdir)

        # 3) make a grid from the stations themselves with 200m cell spacing
        self.mmodel = Model(Data=md, cell_size_east=200, cell_size_north=200)
        self.mmodel.make_mesh()
        # all is good write the mesh file
        self.mmodel.write_model_file(save_path=self.outputdir)

        print ("add topography to res model")
        self.mmodel.add_topography(self.topofile, interp_method='nearest')

        print("make covariance file")
        cov = Covariance(mask_arr=self.mmodel.covariance_mask,
                         save_path=self.outputdir,
                         smoothing_east=0.3,
                         smoothing_north=0.3,
                         smoothing_z=0.3)

        print ("to write cov file:", self.mmodel.model_fn)
        cov.write_covariance_file(model_fn=self.mmodel.model_fn)

    def test_it(self):
        """
        make the setUp to run at least once
        :return:
        """
        assert(True)


    # def test_make_mesh(self):
    #     """
    #     make mesh and output the mesh to a file
    #     :return:
    #     """
    #
    #     self.mmodel.make_mesh()
    #
    #     # all is good write the mesh file
    #     self.mmodel.write_model_file(save_path=self.outputdir)


        # Rotate Mesh
        #
        # self.mmodel.mesh_rotation_angle = 60
        # self.mmodel.make_mesh()


    # def test_plot_mesh(self):
    #     """
    #     check to see if the mesh is what you think it should be
    #     """
    #     #self.fail()
    #
    #     self.mmodel.plot_mesh()
    #
    # def test_add_topo_cov(self):
    #
    #     print ("add topography to res model")
    #     self.mmodel.add_topography(self.topofile, interp_method='nearest')
    #
    #     print("make covariance file")
    #     cov = Covariance(mask_arr=self.mmodel.covariance_mask,
    #                           save_path=self.outputdir,
    #                           smoothing_east=0.3,
    #                           smoothing_north=0.3,
    #                           smoothing_z=0.3)
    #
    #     print ("to write cov file:", self.mmodel.model_fn)
    #     cov.write_covariance_file(model_fn=self.mmodel.model_fn)

        # return True

# def test_add_topography(self):
#     self.fail()
#
#
# def test_project_surface(self):
#     self.fail()
#
#
# def test_assign_resistivity_from_surfacedata(self):
#     self.fail()
#
#
# def test_project_stations_on_topography(self):
#     self.fail()
#
#
#
# def test_write_model_file(self):
#     self.fail()
#
#
# def test_read_model_file(self):
#     self.fail()
#
#
# def test_read_ws_model_file(self):
#     self.fail()
#
#
# def test_write_vtk_file(self):
#     self.fail()
#
#
# def test_write_gocad_sgrid_file(self):
#     self.fail()
#
#
# def test_read_gocad_sgrid_file(self):
#     self.fail()
#
#
# def test_write_xyres(self):
#     self.fail()
#
#
# def test_add_topography_to_model(self):
#     self.fail()
#
#
# def test_change_data_elevation(self):
#     self.fail()


if __name__ == '__main__':
    """
    nosetests tests/test_model.py
    python tests/test_model.py
    """
    import unittest
    unittest.main()
