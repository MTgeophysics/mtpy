# configure matplotlib for testing
import matplotlib.pyplot as plt

from mtpy.imaging.penetration_depth3d import plot_latlon_depth_profile
from mtpy.utils.decorator import ImageCompare
from tests.imaging import ImageTestCase


class TestPenetration_depth3d(ImageTestCase):
    @classmethod
    def setUpClass(cls):
        cls._edifiles_small = "tests/data/edifiles"

    @ImageCompare(fig_size=(8, 6))
    def test_plot_latlon_depth_profile_period_index(self):
        plot_latlon_depth_profile(self._edifiles_small, 10, 'det', savefig=False)

    def test_plot_latlon_depth_profile_no_period(self):
        try:
            plot_latlon_depth_profile(self._edifiles_small, savefig=False)
            assert(False)  # if this statement reached, it is wrong
        except Exception as ex:
            print (ex)
            assert(True)

    @ImageCompare(fig_size=(8, 6))
    def test_plot_latlon_depth_profile_period(self):
        plot_latlon_depth_profile(self._edifiles_small, 2.857, savefig=False)
