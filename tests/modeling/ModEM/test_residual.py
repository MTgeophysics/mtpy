# -*- coding: utf-8 -*-
"""
Created on Thu Jun 18 14:50:12 2020

@author: u64125
"""
import glob
import os
from unittest import TestCase
import numpy as np
from mtpy.modeling.modem import Residual
from tests import make_temp_dir, SAMPLE_DIR


class TestResidual(TestCase):
    def setUp(self):
        self._model_dir = os.path.join(SAMPLE_DIR, "ModEM_2")
        self._residual_fn = os.path.join(self._model_dir, "Modular_MPI_NLCG_004.res")
        self.residual_object = Residual(residual_fn=self._residual_fn)
        self.residual_object.read_residual_file()
        self.test_station = "Synth10"
        self.sidx = np.where(
            self.residual_object.residual_array["station"] == self.test_station
        )[0][0]

    def test_read_residual(self):

        expected_residual_z = np.array(
            [
                [
                    [-1.02604300 - 1.02169400e00j, 0.57546600 - 2.86067700e01j],
                    [-1.66748500 + 2.77643100e01j, 1.10572800 + 7.91478700e-01j],
                ],
                [
                    [-0.65690610 - 7.81721300e-01j, 3.19986900 - 1.32354200e01j],
                    [-3.99973300 + 1.25588600e01j, 0.70137690 + 7.18065500e-01j],
                ],
                [
                    [-0.37482920 - 5.49270000e-01j, 2.99884200 - 5.47220900e00j],
                    [-3.48449700 + 4.94168400e00j, 0.38315510 + 5.34895800e-01j],
                ],
                [
                    [-0.22538380 - 3.74392300e-01j, 2.25689600 - 1.52787900e00j],
                    [-2.48044900 + 1.12663400e00j, 0.24971600 + 3.52619400e-01j],
                ],
                [
                    [-0.10237530 - 2.46537400e-01j, 1.22627900 + 2.40738800e-01j],
                    [-1.30554300 - 4.96037400e-01j, 0.17602110 + 1.96752500e-01j],
                ],
                [
                    [-0.06368253 - 1.50662800e-01j, 0.22056860 + 3.82570600e-01j],
                    [-0.15430870 - 6.20512300e-01j, 0.08838979 + 1.58171100e-01j],
                ],
                [
                    [-0.02782632 - 9.68861500e-02j, 0.01294457 + 1.44029600e-01j],
                    [0.04851138 - 2.18578600e-01j, 0.06537638 + 6.12342500e-02j],
                ],
                [
                    [-0.04536762 - 9.97313200e-02j, 0.23600760 + 6.91851800e-02j],
                    [-0.19451440 - 9.89826000e-02j, 0.04738959 + 4.03351500e-03j],
                ],
                [
                    [0.02097188 - 1.58254400e-02j, 0.37978140 + 2.43604500e-01j],
                    [-0.25993120 - 1.23627200e-01j, 0.06108656 - 4.74416500e-02j],
                ],
                [
                    [0.00851522 + 1.43499600e-02j, 0.03704846 + 2.03981800e-01j],
                    [-0.05785917 - 1.39049200e-01j, 0.12167730 - 7.85872700e-02j],
                ],
                [
                    [0.01532396 + 4.53273900e-02j, -0.18341090 + 1.94764500e-01j],
                    [0.17545040 - 1.09104100e-01j, 0.12086500 + 4.31819300e-02j],
                ],
                [
                    [-0.06653160 + 4.91978100e-03j, -0.25800590 + 5.63935700e-02j],
                    [0.23867870 + 2.23033200e-02j, 0.14759170 + 1.04700000e-01j],
                ],
                [
                    [-0.02225152 - 2.51810700e-02j, -0.24345510 - 8.70721200e-02j],
                    [0.18595070 + 2.62889600e-02j, 0.01455407 + 1.57621400e-01j],
                ],
                [
                    [0.01056665 + 1.59123900e-02j, -0.12207970 - 1.49491900e-01j],
                    [0.14559130 - 8.58512800e-03j, -0.07530988 + 1.22342400e-01j],
                ],
                [
                    [-0.02195550 + 7.05037300e-02j, -0.02480397 - 1.14487400e-01j],
                    [0.16408400 - 3.97682300e-02j, -0.12009660 + 7.95397500e-02j],
                ],
                [
                    [-0.10702040 + 6.75635400e-02j, 0.04605616 - 6.53660100e-02j],
                    [0.23191780 - 2.30734300e-02j, -0.13914140 + 4.37319200e-02j],
                ],
                [
                    [-0.11429350 - 5.96242700e-03j, 0.04939716 + 6.49396100e-03j],
                    [0.23719750 + 2.03324600e-02j, -0.16310670 - 1.33806800e-02j],
                ],
            ]
        )

        expected_residual_tip = np.array(
            [
                [[1.99677300e-02 + 7.51403000e-03j, -9.26884200e-03 + 7.61811900e-04j]],
                [[1.56102400e-02 + 9.40827100e-03j, -9.94980800e-03 + 2.77858500e-05j]],
                [[1.00342300e-02 + 9.31690900e-03j, -9.23218900e-03 - 1.33653600e-04j]],
                [[5.68685000e-03 + 8.56526900e-03j, -8.54973100e-03 + 2.46999200e-04j]],
                [[1.07812200e-03 + 8.77348500e-03j, -8.40960800e-03 + 8.43687200e-04j]],
                [[1.18444500e-04 + 5.94769200e-03j, -1.05262400e-02 + 1.20061900e-04j]],
                [
                    [
                        -6.58875300e-03 + 6.67077300e-03j,
                        -1.49532000e-02 + 7.55132400e-03j,
                    ]
                ],
                [
                    [
                        -1.38596400e-02 + 9.23400700e-03j,
                        -1.75996500e-02 + 7.20104400e-03j,
                    ]
                ],
                [
                    [
                        -2.93723000e-02 - 1.98601600e-02j,
                        -1.66622000e-02 + 1.51911600e-02j,
                    ]
                ],
                [
                    [
                        -3.63411900e-02 - 2.81472900e-02j,
                        -4.21293500e-02 + 1.69332600e-03j,
                    ]
                ],
                [[4.33842300e-03 - 9.09171200e-02j, -5.64092100e-02 - 6.71133000e-02j]],
                [[2.31374400e-02 - 1.23249000e-01j, -4.60595000e-02 - 1.08313400e-01j]],
                [[1.46526100e-01 - 1.77551800e-01j, 2.11499700e-02 - 1.97421400e-01j]],
                [[2.86257700e-01 - 1.26995500e-01j, 1.63245800e-01 - 2.62023200e-01j]],
                [[3.43578700e-01 + 2.36713800e-02j, 3.56271700e-01 - 2.17845400e-01j]],
                [[2.36714300e-01 + 1.50994700e-01j, 4.83549200e-01 - 9.88437900e-02j]],
                [[9.51891900e-02 + 1.36287700e-01j, 4.89243800e-01 + 8.02669000e-02j]],
            ]
        )

        expected_perlist = np.array(
            [
                1.00000000e-02,
                2.05353000e-02,
                4.21697000e-02,
                8.65964000e-02,
                1.77828000e-01,
                3.65174000e-01,
                7.49894000e-01,
                1.53993000e00,
                3.16228000e00,
                6.49382000e00,
                1.33352000e01,
                2.73842000e01,
                5.62341000e01,
                1.15478000e02,
                2.37137000e02,
                4.86968000e02,
                1.00000000e03,
            ]
        )

        assert np.all(
            np.abs(
                self.residual_object.residual_array["z"][self.sidx]
                - expected_residual_z
            )
            / expected_residual_z
            < 1e-6
        )
        assert np.all(
            np.abs(
                self.residual_object.residual_array["tip"][self.sidx]
                - expected_residual_tip
            )
            / expected_residual_tip
            < 1e-6
        )
        assert np.all(
            np.abs(self.residual_object.period_list - expected_perlist)
            / expected_perlist
            < 1e-6
        )

    def test_get_rms(self):
        self.residual_object.get_rms()

        expected_rms = 4.598318
        expected_rms_z = 5.069801
        expected_rms_tip = 3.3062455

        assert np.abs(self.residual_object.rms - expected_rms < 1e-6)
        assert np.abs(self.residual_object.rms_z - expected_rms_z < 1e-6)
        assert np.abs(self.residual_object.rms_tip - expected_rms_tip < 1e-6)

        # expected rms by component for station
        expected_rms_by_component_z = np.array(
            [[1.54969747, 3.45869927], [4.34009684, 2.31467718]]
        )
        expected_rms_by_component_tip = np.array([[3.63839712, 5.18765567]])

        assert np.all(
            np.abs(
                self.residual_object.rms_array["rms_z_component"][self.sidx]
                - expected_rms_by_component_z
            )
            < 1e-6
        )
        assert np.all(
            np.abs(
                self.residual_object.rms_array["rms_tip_component"][self.sidx]
                - expected_rms_by_component_tip
            )
            < 1e-6
        )

        expected_rms_by_period = np.array(
            [
                5.0808857,
                3.47096626,
                2.34709635,
                1.52675356,
                1.09765601,
                0.73261742,
                0.43272026,
                0.62780779,
                1.13149442,
                0.95879571,
                1.71206363,
                2.41720885,
                3.54772784,
                4.66569959,
                5.69325904,
                6.60449661,
                7.3224976,
            ]
        )
        expected_rms_by_period_z = np.array(
            [
                6.21674098,
                4.24399835,
                2.86799775,
                1.86322907,
                1.33659998,
                0.88588157,
                0.47925494,
                0.70885089,
                1.29429634,
                0.91572845,
                1.47599353,
                2.15781454,
                2.45842532,
                2.40733405,
                2.81537801,
                4.54401732,
                6.51544423,
            ]
        )
        expected_rms_by_period_tip = np.array(
            [
                0.38789413,
                0.3460871,
                0.27524824,
                0.22289946,
                0.20383115,
                0.20152558,
                0.31995293,
                0.42129405,
                0.7003091,
                1.03959147,
                2.10626965,
                2.86642089,
                5.066696,
                7.3291025,
                9.02146822,
                9.46371701,
                8.71521005,
            ]
        )

        assert np.all(
            np.abs(
                self.residual_object.rms_array["rms_z_period"][self.sidx]
                - expected_rms_by_period_z
            )
            < 1e-6
        )
        assert np.all(
            np.abs(
                self.residual_object.rms_array["rms_tip_period"][self.sidx]
                - expected_rms_by_period_tip
            )
            < 1e-6
        )
        assert np.all(
            np.abs(
                self.residual_object.rms_array["rms_period"][self.sidx]
                - expected_rms_by_period
            )
            < 1e-6
        )
