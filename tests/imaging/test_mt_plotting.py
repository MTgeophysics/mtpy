# -*- coding: utf-8 -*-
"""
Created on Sun Jun 25 14:45:49 2023

@author: jpeacock
"""

# =============================================================================
# Imports
# =============================================================================
import unittest
from unittest import mock

import matplotlib.pyplot as plt

from mtpy import MT
from mtpy.imaging import PlotMTResponse
from mt_metadata import TF_EDI_CGG

# =============================================================================


@mock.patch("mtpy.imaging.PlotMTResponse.plt")
class TestMTPlotResponse(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        self.m = MT(TF_EDI_CGG)
        self.m.read()

        self.p = self.m.plot_mt_response(show_plot=False)

    def test_z(self):
        self.assertEqual(self.p.Z, self.m.Z)

    def test_tipper(self):
        self.assertEqual(self.p.Tipper, self.m.Tipper)

    def test_pt(self):
        self.assertEqual(self.p.pt, self.m.pt)

    def test_period(self):
        self.assertTrue((self.p.period == self.m.period).all())

    def test_has_tipper(self):
        self.p._has_tipper()
        self.assertEqual(self.p.plot_tipper, "yri")

    def test_has_pt(self):
        self.p._has_pt()
        self.assertTrue(self.p.plot_pt)

    def test_figure_called(self, mock_plot):
        assert mock_plot.figure.called


# @mock.patch("mtpy.imaging.PlotMTResponse.plt")
# def test_module(mock_plt):
#     m = MT(TF_EDI_CGG)
#     m.read()
#     p = m.plot_mt_response(show_plot=False)

#     # Assert plt.figure got called
#     assert mock_plt.figure.not_called

#     assert mock_plt.figure.clf.called


# =============================================================================
# run
# =============================================================================
if __name__ == "__main__":
    unittest.main()
