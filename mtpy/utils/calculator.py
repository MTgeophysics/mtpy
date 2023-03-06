#!/usr/bin/env python

"""
mtpy/utils/calculator.py

Helper functions for standard calculations, e.g. error propagation


@UofA, 2013
(LK)

"""

# =================================================================


import numpy as np
import math, cmath

import mtpy.utils.exceptions as MTex


# =================================================================

# define uncertainty for differences between time steps
epsilon = 1e-9
# magnetic permeability in free space in H/m (=Vs/Am)
mu0 = 4e-7 * np.pi

# =================================================================


def centre_point(xarray, yarray):
    """
    get the centre point of arrays of x and y values
    """
    return (xarray.max() + xarray.min()) / 2.0, (
        yarray.max() + yarray.min()
    ) / 2.0


def roundsf(number, sf):
    """
    round a number to a specified number of significant figures (sf)
    """
    # can't have < 1 s.f.
    sf = max(sf, 1.0)
    rounding = int(np.ceil(-np.log10(number) + sf - 1.0))

    return np.round(number, rounding)


def get_period_list(
    period_min, period_max, periods_per_decade, include_outside_range=True
):
    """
    get a list of values (e.g. periods), evenly spaced in log space and
    including values on multiples of 10

    :returns:
        numpy array containing list of values

    :inputs:
        period_min = minimum period
        period_max = maximum period
        periods_per_decade = number of periods per decade
        include_outside_range = option whether to start and finish the period
                                list just inside or just outside the bounds
                                specified by period_min and period_max
                                default True

    """

    log_period_min = np.log10(period_min)
    log_period_max = np.log10(period_max)

    # check if log_period_min is a whole number
    if log_period_min % 1 > 0:
        # list of periods, around the minimum period, that will be present in specified
        # periods per decade
        aligned_logperiods_min = np.linspace(
            np.floor(log_period_min),
            np.ceil(log_period_min),
            periods_per_decade + 1,
        )
        lpmin_diff = log_period_min - aligned_logperiods_min
        # index of starting period, smallest value > 0
        if include_outside_range:
            spimin = np.where(lpmin_diff > 0)[0][-1]
        else:
            spimin = np.where(lpmin_diff < 0)[0][0]
        start_period = aligned_logperiods_min[spimin]
    else:
        start_period = log_period_min

    if log_period_max % 1 > 0:
        # list of periods, around the maximum period, that will be present in specified
        # periods per decade
        aligned_logperiods_max = np.linspace(
            np.floor(log_period_max),
            np.ceil(log_period_max),
            periods_per_decade + 1,
        )
        lpmax_diff = log_period_max - aligned_logperiods_max
        # index of starting period, smallest value > 0
        if include_outside_range:
            spimax = np.where(lpmax_diff < 0)[0][0]
        else:
            spimax = np.where(lpmax_diff > 0)[0][-1]
        stop_period = aligned_logperiods_max[spimax]
    else:
        stop_period = log_period_max

    return np.logspace(
        start_period,
        stop_period,
        int((stop_period - start_period) * periods_per_decade + 1),
    )


def nearest_index(val, array):
    """
    find the index of the nearest value in the array
    :param val: the value to search for
    :param array: the array to search in

    :return: index: integer describing position of nearest value in array

    """
    # absolute difference between value and array
    diff = np.abs(array - val)

    return np.where(diff == min(diff))[0][0]


def make_log_increasing_array(
    z1_layer, target_depth, n_layers, increment_factor=0.999
):
    """
    create depth array with log increasing cells, down to target depth,
    inputs are z1_layer thickness, target depth, number of layers (n_layers)
    """

    # make initial guess for maximum cell thickness
    max_cell_thickness = target_depth
    # make initial guess for log_z
    log_z = np.logspace(
        np.log10(z1_layer), np.log10(max_cell_thickness), num=n_layers
    )
    counter = 0

    while np.sum(log_z) > target_depth:
        max_cell_thickness *= increment_factor
        log_z = np.logspace(
            np.log10(z1_layer), np.log10(max_cell_thickness), num=n_layers
        )
        counter += 1
        if counter > 1e6:
            break

    return log_z


def invertmatrix_incl_errors(inmatrix, inmatrix_error=None):

    if inmatrix is None:
        raise MTex.MTpyError_input_arguments("Matrix must be defined")

    if (inmatrix_error is not None) and (
        inmatrix.shape != inmatrix_error.shape
    ):
        raise MTex.MTpyError_input_arguments(
            "Matrix and err-matrix shapes do not match: %s - %s"
            % (str(inmatrix.shape), str(inmatrix_error.shape))
        )

    if inmatrix.shape[-2] != inmatrix.shape[-1]:
        raise MTex.MTpyError_input_arguments("Matrices must be square!")

    if (inmatrix_error is not None) and (
        inmatrix_error.shape[-2] != inmatrix_error.shape[-1]
    ):
        raise MTex.MTpyError_input_arguments("Matrices must be square!")

    dim = inmatrix.shape[-1]

    det = np.linalg.det(inmatrix)

    if det == 0:
        raise MTex.MTpyError_input_arguments(
            "Matrix is singular - I cannot invert that!"
        )

    inv_matrix = np.zeros_like(inmatrix)

    if dim != 2:
        raise MTex.MTpyError_input_arguments("Only 2D matrices supported yet")

    inv_matrix = np.linalg.inv(inmatrix)

    inv_matrix_err = None

    if inmatrix_error is not None:
        inmatrix_error = np.real(inmatrix_error)
        inv_matrix_err = np.zeros_like(inmatrix_error)

        for i in range(2):
            for j in range(2):
                # looping over the entries of the error matrix
                err = 0.0
                for k in range(2):
                    for l in range(2):
                        # each entry has 4 summands

                        err += np.abs(
                            -inv_matrix[i, k]
                            * inv_matrix[l, j]
                            * inmatrix_error[k, l]
                        )

                inv_matrix_err[i, j] = err

    return inv_matrix, inv_matrix_err


def rhophi2z(rho, phi, freq):
    """
    Convert impedance-style information given in Rho/Phi format into
    complex valued Z.

    Input:
    rho - 2x2 array (real) - in Ohm m
    phi - 2x2 array (real) - in degrees
    freq - scalar - frequency in Hz

    Output:
    Z - 2x2 array (complex)
    """

    try:
        if rho.shape != (2, 2) or phi.shape != (2, 2):
            raise
        if not (
            rho.dtype in ["float", "int"] and phi.dtype in ["float", "int"]
        ):
            raise

    except:
        raise MTex.MTpyError_input_arguments(
            "ERROR - arguments must be two 2x2 arrays (real)"
        )

    z = np.zeros((2, 2), "complex")
    for i in range(2):
        for j in range(2):
            abs_z = np.sqrt(5 * freq * rho[i, j])
            z[i, j] = cmath.rect(abs_z, math.radians(phi[i, j]))

    return z


def compute_determinant_error(
    z_array, z_err_array, method="theoretical", repeats=1000
):
    """
    compute the error of the determinant of z using a stochastic method
    seed random z arrays with a normal distribution around the input array

    :param z_array: z (impedance) array containing real and imaginary values
    :param z_err_array: impedance error array containing real values,
                        in MT we assume the real and imag errors are the same
    :param method: method to use, theoretical calculation or stochastic

    :return: error: array of real values with same shape as z_err_array
                    representing the error in the determinant of Z
    :return: error_sqrt: array of real values with same shape as z_err_array
                    representing the error in the (determinant of Z)**0.5

    """
    if method == "stochatic":
        arraylist = []

        for r in range(repeats):
            errmag = np.random.normal(
                loc=0, scale=z_err_array, size=z_array.shape
            )
            arraylist = np.append(arraylist, z_array + errmag * (1.0 + 1j))

        arraylist = arraylist.reshape(repeats, z_array.shape[0], 2, 2)
        detlist = np.linalg.det(arraylist)
        error = np.std(detlist, axis=0)

    else:
        error = np.abs(
            z_err_array[:, 0, 0] * np.abs(z_array[:, 1, 1])
            + z_err_array[:, 1, 1] * np.abs(z_array[:, 0, 0])
            - z_err_array[:, 0, 1] * np.abs(z_array[:, 1, 0])
            - z_err_array[:, 1, 0] * np.abs(z_array[:, 0, 1])
        )

    return error


def propagate_error_polar2rect(r, r_error, phi, phi_error):
    """
    Find error estimations for the transformation from polar to cartesian
    coordinates.

    Uncertainties in polar representation define a section of an annulus.
    Find the 4 corners of this section and additionally the outer boundary
    point, which is defined by phi = phi0, rho = rho0 + sigma rho.
    The cartesian "box" defining the uncertainties in x,y is the outer
    bound around the annulus section, defined by the four outermost
    points. So check the four corners as well as the outer boundary edge
    of the section to find the extrema in x znd y. These give you the
    sigma_x/y.

    """

    corners = [
        (
            np.real(cmath.rect(r - r_error, phi - phi_error)),
            np.imag(cmath.rect(r - r_error, phi - phi_error)),
        ),
        (
            np.real(cmath.rect(r + r_error, phi - phi_error)),
            np.imag(cmath.rect(r + r_error, phi - phi_error)),
        ),
        (
            np.real(cmath.rect(r + r_error, phi + phi_error)),
            np.imag(cmath.rect(r + r_error, phi + phi_error)),
        ),
        (
            np.real(cmath.rect(r - r_error, phi + phi_error)),
            np.imag(cmath.rect(r - r_error, phi + phi_error)),
        ),
        (
            np.real(cmath.rect(r + r_error, phi)),
            np.imag(cmath.rect(r + r_error, phi)),
        ),
    ]

    lo_x = [i[0] for i in corners]
    lo_y = [i[1] for i in corners]

    point = (np.real(cmath.rect(r, phi)), np.imag(cmath.rect(r, phi)))
    lo_xdiffs = [abs(point[0] - i) for i in lo_x]
    lo_ydiffs = [abs(point[1] - i) for i in lo_y]

    xerr = max(lo_xdiffs)
    yerr = max(lo_ydiffs)

    return xerr, yerr


def propagate_error_rect2polar(x, x_error, y, y_error):

    # x_error, y_error define a  rectangular uncertainty box

    # rho error is the difference between the closest and furthest point of the box (w.r.t. the origin)
    # approximation: just take corners and midpoint of edges
    lo_points = [
        (x + x_error, y),
        (x - x_error, y),
        (x, y - y_error),
        (x, y + y_error),
        (x - x_error, y - y_error),
        (x + x_error, y - y_error),
        (x + x_error, y + y_error),
        (x - x_error, y + y_error),
    ]

    # check, if origin is within the box:
    origin_in_box = False
    if x_error >= np.abs(x) and y_error >= np.abs(y):
        origin_in_box = True

    try:
        lo_polar_points = [cmath.polar(complex(*i)) for i in lo_points]
    except OverflowError:
        lo_polar_points = []
        for i in lo_points:
            real, imag = i
            if abs(real) > 1e32:
                real = 1e32
            if abs(imag) > 1e32:
                imag = 1e32
            lo_polar_points.append(cmath.polar(complex(real, imag)))

    lo_rho = [i[0] for i in lo_polar_points]
    lo_phi = [math.degrees(i[1]) % 360 for i in lo_polar_points]

    rho_err = 0.5 * (max(lo_rho) - min(lo_rho))
    phi_err = 0.5 * (max(lo_phi) - min(lo_phi))

    if (270 < max(lo_phi) < 360) and (0 < min(lo_phi) < 90):
        tmp1 = [i for i in lo_phi if (0 < i < 90)]
        tmp4 = [i for i in lo_phi if (270 < i < 360)]
        phi_err = 0.5 * ((max(tmp1) - min(tmp4)) % 360)

    if phi_err > 180:
        phi_err = (-phi_err) % 360

    if origin_in_box is True:
        # largest rho:
        rho_err = 2 * rho_err + min(lo_rho)
        # maximum angle uncertainty:
        phi_err = 180.0

    return rho_err, phi_err


def z_error2r_phi_error(z_real, z_imag, error):
    """
    Error estimation from rectangular to polar coordinates.

    By standard error propagation, relative error in resistivity is
    2*relative error in z amplitude.

    Uncertainty in phase (in degrees) is computed by defining a circle around
    the z vector in the complex plane. The uncertainty is the absolute angle
    between the vector to (x,y) and the vector between the origin and the
    tangent to the circle.

    :returns:
        tuple containing relative error in resistivity, absolute error in phase

    :inputs:
        z_real = real component of z (real number or array)
        z_imag = imaginary component of z (real number or array)
        error = absolute error in z (real number or array)

    """

    z_amp = np.abs(z_real + 1j * z_imag)

    if z_amp == 0:
        z_rel_err = 0.0
    else:
        z_rel_err = error / z_amp

    res_rel_err = 2.0 * z_rel_err

    # if the relative error of the amplitude is >=100% that means that the relative
    # error of the resistivity is 200% - that is then equivalent to an uncertainty
    # in the phase angle of 90 degrees:
    if np.iterable(z_real):
        phi_err = np.degrees(np.arctan(z_rel_err))
        phi_err[res_rel_err > 1.0] = 90.0

    else:
        if res_rel_err > 1.0:
            phi_err = 90
        else:
            phi_err = np.degrees(np.arctan(z_rel_err))

    return res_rel_err, phi_err


def old_z_error2r_phi_error(x, x_error, y, y_error):
    """
    Error estimation from rect to polar, but with small variation needed for
    MT: the so called 'relative phase error' is NOT the relative phase error,
    but the ABSOLUTE uncertainty in the angle that corresponds to the relative
    error in the amplitude.

    So, here we calculate the transformation from rect to polar coordinates,
    esp. the absolute/length of the value. Then we find the uncertainty in
    this length and calculate the relative error of this. The relative error of
    the resistivity will be double this value, because it's calculated by taking
    the square of this length.

    The relative uncertainty in length defines a circle around (x,y)
    (APPROXIMATION!). The uncertainty in phi is now the absolute of the
    angle beween the vector to (x,y) and the origin-vector tangential to the
    circle.
    BUT....since the phase angle uncertainty is interpreted with regard to
    the resistivity and not the Z-amplitude, we have to look at the square of
    the length, i.e. the relative error in question has to be halfed to get
    the correct relationship between resistivity and phase errors!!.

    """

    # x_error, y_error define a  rectangular uncertainty box

    # rho error is the difference between the closest and furthest point of the box (w.r.t. the origin)
    # approximation: just take corners and midpoint of edges
    lo_points = [
        (x + x_error, y),
        (x - x_error, y),
        (x, y - y_error),
        (x, y + y_error),
        (x - x_error, y - y_error),
        (x + x_error, y - y_error),
        (x + x_error, y + y_error),
        (x - x_error, y + y_error),
    ]

    lo_polar_points = [cmath.polar(np.complex(*i)) for i in lo_points]

    lo_rho = [i[0] for i in lo_polar_points]

    # uncertainty in amplitude is defined by half the diameter of the box around x,y
    rho_err = 0.5 * (max(lo_rho) - min(lo_rho))

    rho = cmath.polar(np.complex(x, y))[0]
    try:
        rel_error_rho = rho_err / rho
    except:
        rel_error_rho = 0.0

    # if the relative error of the amplitude is >=100% that means that the relative
    # error of the resistivity is 200% - that is then equivalent to an uncertainty
    # in the phase angle of 90 degrees:
    if rel_error_rho > 1.0:
        phi_err = 90
    else:
        phi_err = np.degrees(np.arcsin(rel_error_rho))

    return rho_err, phi_err


# rotation:
# 1. rotation positive in clockwise direction
# 2. orientation of new X-axis X' given by rotation angle
# 3. express contents of Z/tipper (points P) in this new system (points P')
# 4. rotation for points calculated as P' = ([cos , sin ],[-sin, cos]) * P <=> P' = R * P
# 5. => B' = R * B and E' = R * E
# (Rt is the inverse rotation matrix)
# 6. E = Z * B => Rt * E' = Z * Rt * B' => E' = (R*Z*Rt) * B' => Z' = (R*Z*Rt)

# 7. Bz = T * B => Bz = T * Rt * B' => T' = (T * Rt)
# (Since Tipper is a row vector)
# 7.a for general column vectors v:  v' = R * v

# Rotation of the uncertainties:
# a) rotate Z into Z'
# b) use propagation of errors on Z' to obtain the rotated Z'err
# That is NOT the same as the rotated error matrix Zerr (although the result is similar)


def rotate_matrix_with_errors(in_matrix, angle, error=None):
    """

    Rotate a matrix including errors clockwise given an angle in degrees.

    :param in_matrix: A n x 2 x 2  matrix to rotate
    :type inmatrix: np.ndarray

    :param angle: Angle to rotate by assuming clockwise positive from
                 0 = north
    :type angle: float

    :param error: A n x 2 x 2 matrix of associated errors,
                        defaults to None
    :type error: np.ndarray, optional

    :raises MTex: If input array is incorrect

    :return: rotated matrix
    :rtype: np.ndarray

    :return: rotated matrix errors
    :rtype: np.ndarray

    """

    if in_matrix is None:
        raise MTex.MTpyError_inputarguments("Matrix must be defined")

    if (error is not None) and (in_matrix.shape != error.shape):
        msg = "matricies are not the same shape in_matrix={0}, err={1}".format(
            in_matrix.shape, error.shape
        )
        raise MTex.MTpyError_inputarguments(msg)

    try:
        phi = np.deg2rad(float(angle) % 360)
    except TypeError:
        raise MTex.MTpyError_inputarguments('"Angle" must be a float')

    cphi = np.cos(phi)
    sphi = np.sin(phi)

    # clockwise rotation matrix is given by [[cos, -sin], [sin, cos]]
    rot_mat = np.array([[cphi, -sphi], [sphi, cphi]])
    rotated_matrix = np.dot(np.dot(rot_mat, in_matrix), np.linalg.inv(rot_mat))

    err_mat = None
    if error is not None:
        err_orig = np.real(error)
        err_mat = np.zeros_like(error)

        # standard propagation of errors:
        err_mat[0, 0] = np.sqrt(
            (cphi**2 * err_orig[0, 0]) ** 2
            + (cphi * sphi * err_orig[0, 1]) ** 2
            + (cphi * sphi * err_orig[1, 0]) ** 2
            + (sphi**2 * err_orig[1, 1]) ** 2
        )
        err_mat[0, 1] = np.sqrt(
            (cphi**2 * err_orig[0, 1]) ** 2
            + (cphi * sphi * err_orig[1, 1]) ** 2
            + (cphi * sphi * err_orig[0, 0]) ** 2
            + (sphi**2 * err_orig[1, 0]) ** 2
        )
        err_mat[1, 0] = np.sqrt(
            (cphi**2 * err_orig[1, 0]) ** 2
            + (cphi * sphi * err_orig[1, 1]) ** 2
            + (cphi * sphi * err_orig[0, 0]) ** 2
            + (sphi**2 * err_orig[0, 1]) ** 2
        )
        err_mat[1, 1] = np.sqrt(
            (cphi**2 * err_orig[1, 1]) ** 2
            + (cphi * sphi * err_orig[0, 1]) ** 2
            + (cphi * sphi * err_orig[1, 0]) ** 2
            + (sphi**2 * err_orig[0, 0]) ** 2
        )

    return rotated_matrix, err_mat


def rotate_vector_with_errors(in_vector, angle, error=None):
    """

    Rotate a vector including errors clockwise given an angle in degrees.

    :param in_matrix: A n x 1 x 2  vector to rotate
    :type invector: np.ndarray

    :param angle: Angle to rotate by assuming clockwise positive from
                 0 = north
    :type angle: float

    :param error: A n x 1 x 2 vector of associated errors,
                        defaults to None
    :type error: np.ndarray, optional

    :raises MTex: If input array is incorrect

    :return: rotated vector
    :rtype: np.ndarray

    :return: rotated vector errors
    :rtype: np.ndarray

    """

    if in_vector is None:
        raise MTex.MTpyError_inputarguments(
            "Vector AND error-vector must" + " be defined"
        )

    if (error is not None) and (in_vector.shape != error.shape):
        msg = "matricies are not the same shape in_vector={0}, err={1}".format(
            in_vector.shape, error.shape
        )
        raise MTex.MTpyError_inputarguments(msg)

    try:
        phi = np.deg2rad(float(angle) % 360)
    except TypeError:
        raise MTex.MTpyError_inputarguments('"Angle" must be a float')

    cphi = np.cos(phi)
    sphi = np.sin(phi)

    rot_mat = np.array([[cphi, -sphi], [sphi, cphi]])

    if in_vector.shape == (1, 2):
        rotated_vector = np.dot(in_vector, np.linalg.inv(rot_mat))
    else:
        rotated_vector = np.dot(rot_mat, in_vector)

    err_vector = None
    if error is not None:
        err_vector = np.zeros_like(error)

        if error.shape == (1, 2):
            err_vector = np.dot(error, np.abs(np.linalg.inv(rot_mat)))
        else:
            err_vector = np.dot(np.abs(rot_mat), error)

    return rotated_vector, err_vector


def multiplymatrices_incl_errors(
    inmatrix1, inmatrix2, inmatrix1_error=None, inmatrix2_error=None
):

    if inmatrix1 is None or inmatrix2 is None:
        raise MTex.MTpyError_input_arguments(
            "ERROR - two 2x2 arrays needed as input"
        )

    if inmatrix1.shape != inmatrix2.shape:
        raise MTex.MTpyError_input_arguments(
            "ERROR - two 2x2 arrays with same dimensions needed as input"
        )

    prod = np.array(np.dot(np.matrix(inmatrix1), np.matrix(inmatrix2)))

    if (inmatrix1_error is None) or (inmatrix1_error is None):
        return prod, None

    var = np.zeros((2, 2))
    var[0, 0] = (
        (inmatrix1_error[0, 0] * inmatrix2[0, 0]) ** 2
        + (inmatrix1_error[0, 1] * inmatrix2[1, 0]) ** 2
        + (inmatrix2_error[0, 0] * inmatrix1[0, 0]) ** 2
        + (inmatrix2_error[1, 0] * inmatrix1[0, 1]) ** 2
    )
    var[0, 1] = (
        (inmatrix1_error[0, 0] * inmatrix2[0, 1]) ** 2
        + (inmatrix1_error[0, 1] * inmatrix2[1, 1]) ** 2
        + (inmatrix2_error[0, 1] * inmatrix1[0, 0]) ** 2
        + (inmatrix2_error[1, 1] * inmatrix1[0, 1]) ** 2
    )
    var[1, 0] = (
        (inmatrix1_error[1, 0] * inmatrix2[0, 0]) ** 2
        + (inmatrix1_error[1, 1] * inmatrix2[1, 0]) ** 2
        + (inmatrix2_error[0, 0] * inmatrix1[1, 0]) ** 2
        + (inmatrix2_error[1, 0] * inmatrix1[1, 1]) ** 2
    )
    var[1, 1] = (
        (inmatrix1_error[1, 0] * inmatrix2[0, 1]) ** 2
        + (inmatrix1_error[1, 1] * inmatrix2[1, 1]) ** 2
        + (inmatrix2_error[0, 1] * inmatrix1[1, 0]) ** 2
        + (inmatrix2_error[1, 1] * inmatrix1[1, 1]) ** 2
    )

    return prod, np.sqrt(var)


def reorient_data2D(x_values, y_values, x_sensor_angle=0, y_sensor_angle=90):
    """
    Re-orient time series data of a sensor pair, which has not been in default (x=0, y=90) orientation.

    Input:
    - x-values - Numpy array
    - y-values - Numpy array
    Note: same length for both! - If not, the shorter length is taken

    Optional:
    - Angle of the x-sensor - measured in degrees, clockwise from North (0)
    - Angle of the y-sensor - measured in degrees, clockwise from North (0)

    Output:
    - corrected x-values (North)
    - corrected y-values (East)
    """

    x_values = np.array(x_values)
    y_values = np.array(y_values)

    try:
        if x_values.dtype not in ["complex", "float", "int"]:
            raise
        if len(x_values) != len(y_values):
            raise
    except:
        raise MTex.MTpyError_input_arguments(
            "ERROR - both input arrays must be of same length"
        )

    if len(x_values) != len(y_values):
        l = min(len(x_values), len(y_values))
        x_values = x_values[:l]
        y_values = y_values[:l]

    in_array = np.zeros((len(x_values), 2), x_values.dtype)

    in_array[:, 0] = x_values
    in_array[:, 1] = y_values

    try:
        x_angle = math.radians(x_sensor_angle)
        y_angle = math.radians(y_sensor_angle)
    except:
        raise MTex.MTpyError_input_arguments(
            "ERROR - both angles must be of type int or float"
        )

    T = np.matrix(
        [
            [np.real(cmath.rect(1, x_angle)), np.imag(cmath.rect(1, x_angle))],
            [np.real(cmath.rect(1, y_angle)), np.imag(cmath.rect(1, y_angle))],
        ]
    )

    try:
        new_array = np.dot(in_array, T.I)
    except:
        raise MTex.MTpyError_input_arguments(
            "ERROR - angles must define independent axes to span 2D"
        )

    # print new_array.shape

    return new_array[:, 0], new_array[:, 1]
