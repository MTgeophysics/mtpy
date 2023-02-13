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
        self.rp.resistivity = np.array(
            [
                [np.logspace(-1, 2, 10), np.logspace(1, 3, 10)],
                [np.logspace(1, 3, 10), np.logspace(-1, 2, 10)],
            ]
        ).T

        self.rp.phase = np.array(
            [
                [np.linspace(-110, 110, 10), np.linspace(10, 90, 10)],
                [np.linspace(-170, -110, 10), np.linspace(-110, 110, 10)],
            ]
        ).T

        self.rp.frequency = np.logspace(-3, 3, 10)

    def test_initialize(self):
        rp = ResPhase()
        for name in [
            "resistivity",
            "resistivity_err",
            "resistivity_model_err",
            "phase",
            "phase_err",
            "phase_model_err",
        ]:
            with self.subTest(name):
                self.assertEqual(None, getattr(rp, name))

    def test_z(self):
        with self.subTest("element 0"):
            self.assertTrue(
                np.isclose(
                    self.rp._z[0],
                    np.array(
                        [
                            [
                                -0.0076478 - 0.02101217j,
                                -0.22020971 - 0.03882891j,
                            ],
                            [
                                0.22020971 + 0.03882891j,
                                -0.0076478 - 0.02101217j,
                            ],
                        ]
                    ),
                ).all()
            )

        with self.subTest("element 10"):
            self.assertTrue(
                np.isclose(
                    self.rp._z[-1],
                    np.array(
                        [
                            [
                                -2.41844763e02 + 664.46302439j,
                                -7.64780290e02 - 2101.21657803j,
                            ],
                            [
                                1.36919675e-13 + 2236.0679775j,
                                -2.41844763e02 + 664.46302439j,
                            ],
                        ]
                    ),
                ).all()
            )


# =============================================================================
#
# =============================================================================
if __name__ == "__main__":
    unittest.main()
