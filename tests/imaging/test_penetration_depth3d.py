import pytest
from mtpy.imaging.penetration_depth3d import plot_latlon_depth_profile
from mtpy.imaging.penetration_depth3d import plot_many_periods
from tests.imaging import ImageTestCase, ImageCompare


class TestPenetration_depth3d(ImageTestCase):
    @classmethod
    def setUpClass(cls):
        super(TestPenetration_depth3d, cls).setUpClass()
        cls._edifiles_small = "data/edifiles"

    @ImageCompare(fig_size=(8, 6))
    def test_plot_latlon_depth_profile_period_index(self):
        plot_latlon_depth_profile(
            self._edifiles_small, 10, "det", showfig=False, savefig=False
        )

    @ImageCompare(fig_size=(8, 6))
    def test_plot_latlon_depth_profile_period(self):
        plot_latlon_depth_profile(self._edifiles_small, 2.857, savefig=False)

    def test_plot_latlon_depth_profile_no_period(self):
        with self.assertRaises(Exception):
            plot_latlon_depth_profile(
                self._edifiles_small, showfig=False, savefig=False
            )

    def test_plot_many_periods(self):
        plot_many_periods(self._edifiles_small, n_periods=3)
