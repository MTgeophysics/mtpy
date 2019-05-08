#!/usr/bin/env python

"""
Helper functions for data interpolation.

The functions deal with the interpolation of
calibration files (only LEMI coils so far), non evenly sampled data, ...



@UofA, 2014
(LK)

"""

#=================================================================


import numpy as np
import math
import cmath
import sys


def interpolate_instrumentresponse(
        freq, instrument_response, instr_type='lemi'):
    """
    interpolate instrument-response style files.
    Wrapper for different calls, varying by type.

    Current state: only LEMI-coils implemented

    """

    if instr_type.lower() == 'lemi':
        return interpolate_lemi_coils_response(freq, instrument_response)

    else:
        print('\n\tERROR - instrument type', instr_type, ' not implemented yet\n')
        sys.exit()


def interpolate_lemi_coils_response(freq, instrument_response):
    """
        Find the best interpolated value for a given frequency from an instrument
        response array with LEMI characteristics.

        Input:
        - frequency of interest
        - instrument response array (dimension Nx3): frequency, real, imaginary

        Output:
        - complex valued instrument response for the given frequency

        LEMI characteristics:
        * frequencies <0.4 Hz:
          linear interpolation in LogLog scale (applied to real/imaginary values)
        * frequencies >=0.4 Hz < 1000 Hz:
          linear interpolation in LogLin scale of absolute/phase angle of the complex values
        * frequencies >= 1000 Hz:
          linear interpolation in LogLin scale of real/imaginary parts of the complex values

    """

    # find the value closest to the current freq, assume it's lower
    closest_lower = np.abs(freq - instrument_response[:, 0]).argmin()

    # in case the closest frequency value is not lower but higher,
    # take the freq value below as lower bound for the interval:
    if instrument_response[closest_lower, 0] > freq:
        closest_lower -= 1

    # define the interval:
    instrfreq1 = instrument_response[closest_lower, 0]
    instrfreq2 = instrument_response[closest_lower + 1, 0]

    # take the interval values:
    realval1 = instrument_response[closest_lower, 1]
    realval2 = instrument_response[closest_lower + 1, 1]
    imagval1 = instrument_response[closest_lower, 2]
    imagval2 = instrument_response[closest_lower + 1, 2]

    # for linear interpolation in abs/angle instead of real/imag:
    absval1 = np.abs(np.complex(realval1, imagval1))
    phival1 = np.angle(np.complex(realval1, imagval1)) / np.pi * 180
    absval2 = np.abs(np.complex(realval2, imagval2))
    phival2 = np.angle(np.complex(realval2, imagval2)) / np.pi * 180

    # interpolate real and imaginary part independently in log-space:
    logfreq1 = np.log(instrfreq1)
    logfreq2 = np.log(instrfreq2)

    loginterval = logfreq2 - logfreq1
    logfreq = np.log(freq)
    weight = (logfreq2 - logfreq) / loginterval

    # for low frequencies take the log of the values to get into loglog space:
    if freq <= 0.4:
        logrealval1 = np.log(realval1)
        logrealval2 = np.log(realval2)
        logimagval1 = np.log(imagval1)
        logimagval2 = np.log(imagval2)

        interpval_real = np.exp(weight * logrealval1 +
                                (1 - weight) * logrealval2)
        interpval_imag = np.exp(weight * logimagval1 +
                                (1 - weight) * logimagval2)

    # nominal range goes up to 1000 Hz
    elif 0.4 < freq <= 1000:
        # linear interpolation on res phase instead of real/imag
        interpval_abs = weight * absval1 + (1 - weight) * absval2
        interpval_phi = weight * phival1 + (1 - weight) * phival2
        interpval_real = np.real(
            cmath.rect(
                interpval_abs,
                interpval_phi /
                180. *
                np.pi))
        interpval_imag = np.imag(
            cmath.rect(
                interpval_abs,
                interpval_phi /
                180. *
                np.pi))

    else:
        interpval_real = weight * realval1 + (1 - weight) * realval2
        interpval_imag = weight * imagval1 + (1 - weight) * imagval2

    interpolated_value = np.complex(interpval_real, interpval_imag)

    return interpolated_value
