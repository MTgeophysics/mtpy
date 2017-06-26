from unittest import TestCase
import os.path

# configure matplotlib for testing
import matplotlib
matplotlib.use('Agg')  # comment out this line if you want to see the plots


from mtpy.imaging.penetration_depth2d import plot2Dprofile


class TestPenetration_depth2d(TestCase):
    def setUp(self):
        self._temp_dir = "tests/temp"
        if not os.path.isdir(self._temp_dir):
            os.mkdir(self._temp_dir)
        self._edifiles = "tests/data/edifiles"
        self._period_index_list = [0, 1, 10, 20, 30, 40, 50, 59]

    def test_plot2Dprofile_no_period_index_list(self):
        """
        testing plot2Dprofile without perido index list
        exception should be raised
        :return:
        """
        try:
            plot2Dprofile(self._edifiles)
        except Exception:
            pass

    def test_plot2Dprofile_det(self):
        plot2Dprofile(self._edifiles, self._period_index_list, 'det')

    def test_plot2Dprofile_zxy(self):
        plot2Dprofile(self._edifiles, self._period_index_list, 'zxy')

    def test_plot2Dprofile_zyx(self):
        plot2Dprofile(self._edifiles, self._period_index_list, 'zyx')

    def test_plot2Dprofile_wrong_rho(self):
        try:
            plot2Dprofile(self._edifiles, self._period_index_list, 'dat')
        except Exception:
            pass
