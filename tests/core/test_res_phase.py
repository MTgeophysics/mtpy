# -*- coding: utf-8 -*-
"""
Created on Fri Oct  7 11:53:55 2022

@author: jpeacock
"""

# =============================================================================
# Imports
# =============================================================================
import unittest

import numpy as np
from mtpy.core.res_phase import ResPhase

# =============================================================================


class TestResPhase(unittest.TestCase):
    def setUp(self):
        self.rp = ResPhase()

    def test_initialize(self):
        for name in [
            "resistivity",
            "resistivity_err",
            "resistivity_model_err",
            "phase",
            "phase_err",
            "phase_model_err",
        ]:
            with self.subTest(name):
                self.assertEqual(None, getattr(self.rp, name))


# =============================================================================
#
# =============================================================================
if __name__ == "__main__":
    unittest.main()
