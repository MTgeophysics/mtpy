#!/usr/bin/env python

"""
mtpy/mtpy/analysis/niblettbostick.py

Contains functions for the calculation of the Niblett-Bostick transformation of
impedance tensors.

The methods follow
- Niblett
- Bostick
- Jones
- J. RODRIGUEZ, F.J. ESPARZA, E. GOMEZ-TREVINO

Niblett-Bostick transformations are possible in 1D and 2D.

@UofA, 2013 (LK)

Updated 2022-09 JP

"""
# =============================================================================
# Imports
# =============================================================================
import numpy as np
import scipy.interpolate as spi

from mtpy.utils import MU0

# =============================================================================

NB_SCALE_PARAMETER = 2.0 * np.pi * MU0


def calculate_niblett_bostick_depth(resistivity, period):
    """
    Use the Niblett-Bostick approximation for depth of penetration in meters

    :param resistivity: DESCRIPTION
    :type resistivity: TYPE
    :return: DESCRIPTION
    :rtype: TYPE

    """

    return np.sqrt(resistivity * period / NB_SCALE_PARAMETER)


def calculate_niblett_bostick_resistivity_weidelt(resistivity, phase):
    """
    Convert a period-dependent pair of resistivity/phase (Ohm meters/rad)
    into resistivity/depth (Ohm meters/meters)

    The conversion uses the simplified transformation without derivatives.

    :param resistivity: DESCRIPTION
    :type resistivity: TYPE
    :param phase: DESCRIPTION
    :type phase: TYPE
    :return: DESCRIPTION
    :rtype: TYPE

    """

    return resistivity * ((np.pi / 2) * np.deg2rad(phase % 90) - 1)


def calculate_niblett_bostick_resistivity_derivatives(resistivity, period):
    """

    Convert a period-dependent pair of resistivity/phase (Ohm meters/rad)
    into resistivity/depth (Ohm meters/meters)

    The conversion uses derivatives.

    :param resistivity: DESCRIPTION
    :type resistivity: TYPE
    :param period: DESCRIPTION
    :type period: TYPE
    :return: DESCRIPTION
    :rtype: TYPE

    """

    log_period = np.log10(period)
    log_resistivity = np.log10(resistivity)
    m = np.gradient(log_resistivity, log_period, edge_order=2)

    # bostick resistivity only valid for -1 < m < 1
    m[m > 1] = np.nan
    m[m < -1] = np.nan

    return resistivity * (1.0 + m) / (1.0 - m)


def calculate_depth_sensitivity(depth, period, rho=100):
    """
    compute sensitivty S(z,sigma, omega)= -kz*exp(-2*kz).
    The result is independent of sigma and freq.
    :param z:
    :param sigma_conduct:
    :param freq:
    :return: the sensitivity vslue
    """

    omega = 2 * np.pi / period

    k = np.sqrt((0.0 + 1j) * omega * MU0 * rho)

    # same as delta=sqrt(2/mu0*sigma*omega)
    p = 1 / np.real(k)  # It is the

    zp = depth / 1000 * p  # zp is normalized Z
    sensitivity = np.abs(-k * zp * np.exp(-2 * k * zp))

    return sensitivity


def calculate_depth_of_investigation(z_object):
    """
    Determine an array of Z_nb (depth dependent Niblett-Bostick transformed Z)
    from the 1D and 2D parts of an impedance tensor array Z.

    The calculation of the Z_nb needs 6 steps:

        1) Determine the dimensionality of the Z(T), discard all 3D parts
        2) Rotate all Z(T) to TE/TM setup (T_parallel/T_ortho)
        3) Transform every component individually by Niblett-Bostick
        4) collect the respective 2 components each for equal/similar depths
        5) interprete them as TE_nb/TM_nb
        6) set up Z_nb(depth)

    If 1D layers occur inbetween 2D layers, the strike angle is undefined therein.
    We take an - arbitrarily chosen - linear interpolation of strike angle for
    these layers, with the values varying between the angles of the bounding
    upper and lower 2D layers (linearly w.r.t. the periods).

    Use the output for instance for the determination of
    NB-transformed phase tensors.

    Note:
        No propagation of errors implemented yet!

    Arguments
    -------------
        *z_object* : mtpy.core.z object


    Returns
    ------------------
        *depth_array* : np.ndarray(num_periods,
                                   dtype=['period', 'depth_min', 'depth_max',
                                          'resistivity_min', 'resistivity_max'])
                        numpy structured array with keywords.
                            - period    --> period in s
                            - depth_min --> minimum depth estimated (m)
                            - depth_max --> maximum depth estimated (m)
                            - resistivity_min --> minimum resistivity estimated (Ohm-m)
                            - resistivity_max --> maximum resistivity estimated (Ohm-m)

    Example
    ------------
        >>> import mtpy.analysis.niblettbostick as nb
        >>> depth_array = nb.calculate_znb(z_object=z1)
        >>> # plot the results
        >>> import matplotlib.pyplot as plt
        >>> fig = plt.figure()
        >>> ax = fig.add_subplot(1,1,1)
        >>> ax.semilogy(depth_array['depth_min'], depth_array['period'])
        >>> ax.semilogy(depth_array['depth_max'], depth_array['period'])
        >>> plt.show()

    """

    if z_object.z.shape[0] > 1:

        dimensions = z_object.estimate_dimensionality()
        angles = z_object.phase_tensor.azimuth

        # reduce actual Z by the 3D layers:
        angles_2d = np.nan_to_num(angles[np.where(dimensions != 3)])
        periods_2d = z_object.period[np.where(dimensions != 3)]

        # interperpolate strike angle onto all periods
        # make a function for strike using only 2d angles
        strike_interp = spi.interp1d(
            periods_2d, angles_2d, bounds_error=False, fill_value=0
        )
        strike_angles = strike_interp(z_object.period)

        # rotate z to be along the interpolated strike angles
        z_object.rotate(strike_angles)

    depth_array = np.zeros(
        z_object.period.shape[0],
        dtype=[
            ("period", float),
            ("depth_xy", float),
            ("depth_yx", float),
            ("depth_det", float),
            ("depth_min", float),
            ("depth_max", float),
            ("resistivity_xy", float),
            ("resistivity_yx", float),
            ("resistivity_det", float),
            ("resistivity_min", float),
            ("resistivity_max", float),
        ],
    )

    depth_array["period"][:] = z_object.period

    for comp in ["xy", "yx", "det"]:
        res = getattr(z_object, f"res_{comp}")

        depth_array[f"depth_{comp}"][:] = calculate_niblett_bostick_depth(
            res, z_object.period
        )

        if z_object.z.shape[0] > 1:
            depth_array[f"resistivity_{comp}"][
                :
            ] = calculate_niblett_bostick_resistivity_derivatives(
                res, z_object.period
            )
        else:
            phase = getattr(z_object, f"phase_{comp}")
            depth_array[f"resistivity_{comp}"][
                :
            ] = calculate_niblett_bostick_resistivity_weidelt(res, phase)

    for x in ["depth", "resistivity"]:
        d = np.array(
            [
                depth_array[f"{x}_det"],
                depth_array[f"{x}_xy"],
                depth_array[f"{x}_yx"],
            ]
        )
        if np.all(np.isnan(d)):
            depth_array[f"{x}_min"] = np.nan
            depth_array[f"{x}_max"] = np.nan
            continue

        with np.warnings.catch_warnings():
            np.warnings.filterwarnings(
                "ignore", r"All-NaN (slice|axis) encountered"
            )
            depth_array[f"{x}_min"] = np.nanmin(d, axis=0)
            depth_array[f"{x}_max"] = np.nanmax(d, axis=0)

    return depth_array
