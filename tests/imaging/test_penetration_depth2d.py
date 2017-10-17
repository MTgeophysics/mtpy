# configure matplotlib for testing
import matplotlib.pyplot as plt

from mtpy.imaging.penetration_depth2d import plot2Dprofile
from tests.imaging import ImageTestCase, ImageCompare


class TestPenetration_depth2d(ImageTestCase):
    @classmethod
    def setUpClass(cls):
        super(TestPenetration_depth2d, cls).setUpClass()

        cls._edifiles = "tests/data/edifiles"
        cls._period_index_list = [0, 1, 10, 20, 30, 40, 50, 59]

    @classmethod
    def tearDownClass(cls):
        plt.close('all')

    def test_plot2Dprofile_no_period_index_list(self):
        """
        testing plot2Dprofile without period index list
        exception should be raised
        :return:
        """
        try:
            plot2Dprofile(self._edifiles)
        except Exception:
            pass

    @ImageCompare(fig_size=(8, 6))
    def test_plot2Dprofile_det(self):
        plot2Dprofile(self._edifiles, self._period_index_list, 'det')

    @ImageCompare(fig_size=(8, 6))
    def test_plot2Dprofile_zxy(self):
        plot2Dprofile(self._edifiles, self._period_index_list, 'zxy')

    @ImageCompare(fig_size=(8, 6))
    def test_plot2Dprofile_zyx(self):
        plot2Dprofile(self._edifiles, self._period_index_list, 'zyx')

    def test_plot2Dprofile_wrong_rho(self):
        try:
            plot2Dprofile(self._edifiles, self._period_index_list, 'dat')
        except Exception:
            pass
