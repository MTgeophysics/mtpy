# -*- coding: utf-8 -*-
"""
Created on Tue Feb  2 11:03:48 2021

:copyright: 
    Jared Peacock (jpeacock@usgs.gov)

:license: MIT

"""
# =============================================================================
# Imports
# =============================================================================
import unittest

import numpy as np

from mtpy.core.tf import TransferFunction

# =============================================================================


class TestTransferFunction(unittest.TestCase):
    """
    Test TransferFunction
    """

    def setUp(self):
        self.tf = TransferFunction()

    def test_input_list_periods(self):
        self.tf.periods = [1, 2, 3]
        self.assertIsInstance(self.tf.periods, np.ndarray)
        self.assertEqual(self.tf.periods.dtype.type, np.float_)

        self.tf.periods = [[1, 2, 3]]
        self.assertIsInstance(self.tf.periods, np.ndarray)
        self.assertEqual(self.tf.periods.dtype.type, np.float_)

    def test_input_periods_fail(self):
        def set_periods(p):
            self.tf.periods = p

        self.assertRaises(ValueError, set_periods, [[1, 2, 3], [1, 2, 3]])


# =============================================================================
# Run
# =============================================================================
if __name__ == "__main__":
    unittest.main()
