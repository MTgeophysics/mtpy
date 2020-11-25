#!/bin/env python
"""
Description:
    Tests VelocityModel class
References:
 
CreationDate:   12/22/17
Developer:      rakib.hassan@ga.gov.au
 
Revision History:
    LastUpdate:     12/22/17   RH
    LastUpdate:     dd/mm/yyyy  Who     Optional description
"""

import os
from unittest import TestCase
from mtpy.imaging.seismic import VelocityModel
import numpy


class Test_VelocityModel(TestCase):
    def test_max_depth_reached(self):
        utils = os.path.dirname(__file__)
        mtpy = os.path.dirname(utils)
        base = os.path.dirname(mtpy)
        examples = os.path.join(base, "examples")
        data = os.path.join(examples, "data")
        ModEM_files = os.path.join(data, "seismic")
        v_fn = os.path.join(ModEM_files, "stacking_velocities.txt")

        v = VelocityModel(v_fn, ni=20)

        print("\n Testing maximum depths reached\n")
        ts = numpy.linspace(0, 20, 100)  # up to 20 s
        maxDepths = []
        maxDepths_5nn = []
        maxDepths_10nn = []
        for cdp in v._cdps:
            maxDepths.append(numpy.max(v.getDepth(cdp, ts)))
            maxDepths_5nn.append(numpy.max(v.getDepth(cdp, ts, nn=5)))
            maxDepths_10nn.append(numpy.max(v.getDepth(cdp, ts, nn=10)))
        # end for
        maxDepths_all = numpy.max(
            v.getDepth(None, ts)
        )  # using mean velocity across the profile

        maxDepths = numpy.array(maxDepths)
        maxDepths_5nn = numpy.array(maxDepths_5nn)
        maxDepths_10nn = numpy.array(maxDepths_10nn)

        mean, std = numpy.mean(maxDepths), numpy.std(maxDepths)
        print(
            "Mean of maximum depths reached, computed for each cdp: %f; std: %f"
            % (mean, std)
        )

        self.assertGreaterEqual(
            mean, 60e3, "Mean of maximum depth reached must >= 60e3 m"
        )
        self.assertGreaterEqual(std, 2e3, "Std. of maximum depth reached must >= 2e3 m")

        mean_5nn, std_5nn = numpy.mean(maxDepths_5nn), numpy.std(maxDepths_5nn)
        print(
            "Mean of maximum depths reached, computed for mean depth profile "
            "of 5 neighbouring cdps, at each cdp location: %f; std: %f"
            % (mean_5nn, std_5nn)
        )

        self.assertGreaterEqual(
            mean_5nn, 60e3, "Mean of maximum depth reached must >= 60e3 m"
        )
        self.assertGreaterEqual(
            std_5nn, 1.5e3, "Std. of maximum depth reached must >= 1.5e3 m"
        )
        self.assertGreater(
            mean_5nn, mean, "Mean should be larger because of lateral averaging"
        )
        self.assertLess(
            std_5nn, std, "Std. should be smaller because of lateral averaging"
        )

        mean_10nn, std_10nn = numpy.mean(maxDepths_10nn), numpy.std(maxDepths_10nn)
        print(
            "Mean of maximum depths reached, computed for mean depth profile "
            "of 10 neighbouring cdps, at each cdp location: %f; std: %f"
            % (mean_10nn, std_10nn)
        )

        self.assertGreaterEqual(
            mean_10nn, 60e3, "Mean of maximum depth reached must >= 60e3 m"
        )
        self.assertGreaterEqual(
            std_10nn, 1.0e3, "Std. of maximum depth reached must >= 1.0e3 m"
        )
        self.assertGreater(
            mean_10nn, mean_5nn, "Mean should be larger because of lateral averaging"
        )
        self.assertLess(
            std_10nn, std_5nn, "Std. should be smaller because of lateral averaging"
        )

        self.assertTrue(
            numpy.allclose(maxDepths_all, mean),
            "Maximum depth computed using mean depth profile of all CDPs must equal"
            " mean of maximum depths computed for depth-profiles at each CDP",
        )

    # end func

    def test_mean_cdp_velocity(self):
        utils = os.path.dirname(__file__)
        mtpy = os.path.dirname(utils)
        base = os.path.dirname(mtpy)
        examples = os.path.join(base, "examples")
        data = os.path.join(examples, "data")
        ModEM_files = os.path.join(data, "seismic")
        v_fn = os.path.join(ModEM_files, "stacking_velocities.txt")

        v = VelocityModel(v_fn, ni=20)

        print("\n Testing mean CDP velocity\n")

        mean, std = (
            numpy.mean(v._cdp_mean_interval_velocity),
            numpy.std(v._cdp_mean_interval_velocity),
        )

        print("Mean interval velocity for all cdps: %f; std: %f" % (mean, std))
        self.assertLessEqual(mean, 6e3, "Mean velocity exceeds 6000 m/s")

    # end func


# end class
