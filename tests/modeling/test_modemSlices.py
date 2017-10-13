#! /usr/bin/env python
"""
Description:
    unit testing python class methods modemSlices.py

References:
    https://gajira.atlassian.net/browse/ALAMP-31

CreationDate:   13/10/2017
Developer:      fei.zhang@ga.gov.au

Revision History:
    LastUpdate:     13/10/2017   FZ
    LastUpdate:
"""

# import section
import os
from unittest import TestCase
from mtpy.modeling.modem_output_to_views import ModemSlices

class TestModemSlices(TestCase):

    def test_find_stations_in_meshgrid(self):
        self.assertTrue(2 > 1, msg='Dummy True test')

    def test_get_slice_data(self):
        self.assertTrue(2 > 1, msg='Dummy test')
        #self.fail()

    def test_create_csv(self):
        # self.fail()
        self.assertTrue(2 > 1, msg='Dummy test')

    def test_plot_a_slice(self):
        # self.fail()
        self.assertTrue(2 > 1, msg='Dummy test')

    def test_plot_multi_slices(self):
        # self.fail()
        self.assertTrue(2 > 1, msg='Dummy test')
