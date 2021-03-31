"""
Testing Modem.Model
=====================

Make sure that the output is the same as what was previously made.

Should add tests on using the Model class

updated - 2021/02/11 (JP) to use Path
"""
import os
from pathlib import Path
from unittest import TestCase

from mtpy.core.edi_collection import EdiCollection
from mtpy.modeling.modem import Data, Model
from tests import make_temp_dir, EDI_DATA_LIST
from tests.imaging import plt_close

# =============================================================================
#
# =============================================================================


class TestModel(TestCase):
    """
    this test suite only validates the functionality of Model objects but does not verify the output files
    """

    @classmethod
    def setUpClass(cls):
        # setup temp dir
        cls._temp_dir = make_temp_dir(cls.__name__)

    def setUp(self):
        # for each test, setup a different output dir
        self._output_dir = make_temp_dir(self._testMethodName, base_dir=self._temp_dir)

        # set the dir to the output from the previously correct run
        self._expected_output_dir = Path(
            self._temp_dir, "expected_model_output", self._testMethodName
        )

        if not self._expected_output_dir.is_dir():
            self._expected_output_dir = None

    def tearDown(self):
        plt_close("all")


def _test_gen(edi_path):
    def _test_func(self):
        if not Path(edi_path).is_dir():
            # input file does not exist, skip test after remove the output dir
            os.rmdir(self._output_dir)
            self.skipTest(f"edi path does not exist: {edi_path}")

        # epsg to project to. Google epsg 'your projection'
        epsg_code = 3112

        # generate data
        edi_list = list(edi_path.glob("*.edi"))
        period_list = EdiCollection(edi_list).select_periods()

        datob = Data(
            edi_list=edi_list,
            inv_mode="1",
            period_list=period_list,
            epsg=epsg_code,
            error_type_tipper="abs",
            error_type_z="egbert",
            comp_error_type=None,
            error_floor=10,
        )
        datob.write_data_file(save_path=self._output_dir)

        # create mesh grid model object
        model = Model(
            stations_object=datob.station_locations,
            Data=datob,
            epsg=epsg_code,
            cell_size_east=10000,
            cell_size_north=10000,  # GA_VIC
            pad_north=8,  # number of padding cells in each of the north and south directions
            pad_east=8,  # number of east and west padding cells
            pad_z=8,  # number of vertical padding cells
            # factor to increase by in padding cells (vertical)
            pad_stretch_v=1.5,
            # factor to increase by in padding cells (horizontal)
            pad_stretch_h=1.5,
            n_air_layers=0,  # number of air layers 0, 10, 20, depend on topo elev height
            res_model=100,  # halfspace resistivity value for initial reference model
            n_layers=50,  # total number of z layers, including air and pad_z
            z1_layer=50,  # first layer thickness metres, depend
            z_target_depth=500000,
        )
        # the data file will be re-write in this method. No topo elev file used yet
        model.make_mesh()
        model.plot_mesh()
        model.plot_mesh_xy()
        model.plot_mesh_xz()

        # write a model file and initialise a resistivity model
        model.write_model_file(save_path=self._output_dir)

    return _test_func


# generate tests
for edi_path in EDI_DATA_LIST:
    _func = _test_gen(edi_path)
    _func.__name__ = "test_{}".format(os.path.basename(edi_path))
    setattr(TestModel, _func.__name__, _func)
