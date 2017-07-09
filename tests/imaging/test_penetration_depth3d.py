from unittest import TestCase
import os

# configure matplotlib for testing
import matplotlib
# use non-interactive backend 'Agg', so that you do not have to see figure pop-out.
matplotlib.use('Agg')  # comment out this line if you want to see the plots 1-by-1 on screen.

from mtpy.imaging.penetration_depth3d import plot_many_periods
from mtpy.imaging.penetration_depth3d import plot_latlon_depth_profile

class TestPenetration_depth3d(TestCase):
    def setUp(self):
        self._temp_dir = "tests/temp"
        if not os.path.isdir(self._temp_dir):
            os.mkdir(self._temp_dir)
        self._edifiles_small = "tests/data/edifiles"

    # def test_plot_many_periods(self):
    #     plot_many_periods(self._edifiles_small)

    def test_plot_latlon_depth_profile_period_index(self):
        plot_latlon_depth_profile(self._edifiles_small, 10, 'det', savefig=False)

    def test_plot_latlon_depth_profile_no_period(self):
        try:
            plot_latlon_depth_profile(self._edifiles_small, savefig=False)
            assert(False)  # if this statement reached, it is wrong
        except Exception, ex:
            print (ex)
            assert(True)

    def test_plot_latlon_depth_profile_period(self):
        plot_latlon_depth_profile(self._edifiles_small, 2.857, savefig=False)
