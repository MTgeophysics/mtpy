# -*- coding: utf-8 -*-
"""
Created on Tue Mar  7 18:24:58 2023

@author: jpeacock
"""
# =============================================================================
# Imports
# =============================================================================
from pathlib import Path
import time

import numpy as np

# =============================================================================
class Startup(object):
    """
    Reads and writes the startup file for Occam2D.

    .. note:: Be sure to look at the Occam 2D documentation for description
              of all parameters

    ========================= =================================================
    Key Words/Attributes      Description
    ========================= =================================================
    data_fn                   full path to data file
    date_time                 date and time the startup file was written
    debug_level               [ 0 | 1 | 2 ] see occam documentation
                              *default* is 1
    description               brief description of inversion run
                              *default* is 'startup created by mtpy'
    diagonal_penalties        penalties on diagonal terms
                              *default* is 0
    format                    Occam file format
                              *default* is 'OCCAMITER_FLEX'
    iteration                 current iteration number
                              *default* is 0
    iterations_to_run         maximum number of iterations to run
                              *default* is 20
    lagrange_value            starting lagrange value
                              *default* is 5
    misfit_reached            [ 0 | 1 ] 0 if misfit has been reached, 1 if it
                              has.  *default* is 0
    misfit_value              current misfit value.  *default* is 1000
    model_fn                  full path to model file
    model_limits              limits on model resistivity values
                              *default* is None
    model_value_steps         limits on the step size of model values
                              *default* is None
    model_values              np.ndarray(num_free_params) of model values
    param_count               number of free parameters in model
    resistivity_start         starting resistivity value.  If model_values is
                              not given, then all values with in model_values
                              array will be set to resistivity_start
    roughness_type            [ 0 | 1 | 2 ] type of roughness
                              *default* is 1
    roughness_value           current roughness value.
                              *default* is 1E10
    save_path                 directory path to save startup file to
                              *default* is current working directory
    startup_basename          basename of startup file name.
                              *default* is Occam2DStartup
    startup_fn                full path to startup file.
                              *default* is save_path/startup_basename
    stepsize_count            max number of iterations per step
                              *default* is 8
    target_misfit             target misfit value.
                              *default* is 1.
    ========================= =================================================

    :Example: ::

        >>> startup = occam2d.Startup()
        >>> startup.data_fn = ocd.data_fn
        >>> startup.model_fn = profile.reg_fn
        >>> startup.param_count = profile.num_free_params
        >>> startup.save_path = r"/home/occam2d/Line1/Inv1"
    """

    def __init__(self, **kwargs):
        self.save_path = Path()
        self.startup_basename = "Occam2DStartup"
        self.startup_fn = None
        self.model_fn = None
        self.data_fn = None
        self.format = "OCCAMITER_FLEX"
        self.date_time = time.ctime()
        self.description = "startup created by mtpy"
        self.iterations_to_run = 20
        self.roughness_type = 1
        self.target_misfit = 1.0
        self.diagonal_penalties = 0
        self.stepsize_count = 8
        self.model_limits = None
        self.model_value_steps = None
        self.debug_level = 1
        self.iteration = 0
        self.lagrange_value = 5.0
        self.roughness_value = 1e10
        self.misfit_value = 1000
        self.misfit_reached = 0
        self.param_count = None
        self.resistivity_start = 2
        self.model_values = None

        for key, value in kwargs.items():
            setattr(self, key, value)

    def write_startup_file(
        self, startup_fn=None, save_path=None, startup_basename=None
    ):
        """
        Write a startup file based on the parameters of startup class.
        Default file name is save_path/startup_basename

        Arguments:
        -----------
            **startup_fn** : string
                             full path to startup file. *default* is None

            **save_path** : string
                            directory to save startup file. *default* is None

            **startup_basename** : string
                                   basename of starup file. *default* is None

        """
        if save_path is not None:
            self.save_path = Path(save_path)

        if self.save_path is None:
            self.save_path = self.data_fn.parent
        if startup_basename is not None:
            self.startup_basename = startup_basename

        if startup_fn is None:
            self.startup_fn = self.save_path.joinpath(self.startup_basename)

        # --> check to make sure all the important input are given
        if self.data_fn is None:
            raise ValueError("Need to input data file name")

        if self.model_fn is None:
            raise ValueError("Need to input model/regularization file name")

        if self.param_count is None:
            raise ValueError("Need to input number of model parameters")

        slines = []
        slines.append(f"{'Format:':<20}{self.format}")
        slines.append(f"{'Description:':<20}{self.description}")

        if self.model_fn.parent == self.save_path:
            slines.append(f"{'Model File:':<20}{self.model_fn.name}")
        else:
            slines.append(f"{'Model File:':<20}{self.model_fn}")

        if self.data_fn.parent == self.save_path:
            slines.append(f"{'Data File:':<20}{self.data_fn.name}")

        else:
            slines.append(f"{'Data File:':<20}{self.data_fn}")

        slines.append(f"{'Date/Time:':<20}{self.date_time}")
        slines.append(f"{'Iterations to run:':<20}{self.iterations_to_run}")
        slines.append("{'Target Misfit:':<20}{self.target_mistfit}")
        slines.append(f"{'Roughness Type:':<20}{self.roughness_type}")
        slines.append(f"{'Diagonal Penalties:':<20}{self.diagonal_penalties}")
        slines.append(f"{'Stepsize Cut Count:':<20}{self.stepsize_count}")
        if self.model_limits is None:
            slines.append(f"{'!Model Limits:':<20}{'none'}")
        else:
            slines.append(
                f"{'Model Limits:':<20}{self.model_limits[0]},{self.model_limits[1]}"
            )
        if self.model_value_steps is None:
            slines.append(f"{'!Model Value Steps:':<20}{'none'}")
        else:
            slines.append(
                f"{'Model Value Steps:':<20}{self.model_value_steps}"
            )
        slines.append(f"{'Debug Level:':<20}{self.debug_level}")
        slines.append(f"{'Iteration:':<20}{self.iteration}")
        slines.append(f"{'Lagrange Value:':<20}{self.lagrange_value}")
        slines.append(f"{'Roughness Value:':<20}{self.roughness_value}")
        slines.append(f"{'Misfit Value:':<20}{self.misfit_value}")
        slines.append(f"{'Misfit Reached:':<20}{self.misfit_reached}")
        slines.append(f"{'Param Count:':<20}{self.param_count}")

        # make an array of starting values if not are given
        if self.model_values is None:
            self.model_values = np.zeros(self.param_count)
            self.model_values[:] = self.resistivity_start

        if self.model_values.shape[0] != self.param_count:
            raise ValueError(
                "length of model vaues array is not equal "
                "to param count {0} != {1}".format(
                    self.model_values.shape[0], self.param_count
                )
            )

        # write out starting resistivity values
        sline = []
        for ii, mv in enumerate(self.model_values):
            sline.append(f"{mv:^10.4f}")
            if np.remainder(ii + 1, 4) == 0:
                sline.append("\n")
                slines.append("".join(list(sline)))
                sline = []
        slines.append("".join(list(sline + ["\n"])))
        # --> write file
        sfid = open(self.startup_fn, "w")
        sfid.writelines(slines)
        sfid.close()

        print("Wrote Occam2D startup file to {0}".format(self.startup_fn))
