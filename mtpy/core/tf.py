# -*- coding: utf-8 -*-
"""
TransferFunction
====================

    * Deals with transfer functions and their covariance if given
    * Z and Tipper will be based on TransferFunction
    
Created on Tue Feb  2 10:49:18 2021

:copyright: 
    Jared Peacock (jpeacock@usgs.gov)

:license: MIT

"""
# =============================================================================
# Imports
# =============================================================================
import numpy as np

from mtpy.utils.mtpy_logger import get_mtpy_logger

# =============================================================================


class Index:
    """
    Container for transfer function index values for each channel

    default is:
        * hx = 0
        * hy = 1
        * hz = 2
        * ex = 3
        * ey = 4

    """

    def __init__(self, **kwargs):
        self.hx = 0
        self.hy = 1
        self.hz = 2
        self.ex = 3
        self.ey = 4

        for k, v in kwargs.items():
            setattr(self, k, v)


class Angles:
    """
    Container for channel angles

    default is:
        * hx = 0
        * hy = 90
        * hz = 90
        * ex = 0
        * ey = 90

    """

    def __init__(self, **kwargs):
        self.hx = 0
        self.hy = 90
        self.hz = 90
        self.ex = 0
        self.ey = 90

        for k, v in kwargs.items():
            setattr(self, k, v)


class TransferFunction:
    """
    Holds transfer function information and their covariances if given.

    Predicted channels are output channels or the electromagnetic response
    channels, namely Ex, Ey, Hz.

    Predictor channels are the input channes or the electromagnetic source
    channels, namely Hx, Hy.  

    :param tf: Transfer function of shape (n_periods, n_predicted, n_predictor)
    :type tf: np.ndarray

    :param sigma_e: Residual covariance matrix with shape 
    (n_periods, n_predicted_channels, n_predicted_channels) namely ex, ey, hz
    :type sigma_e: np.ndarray

    :param sigma_s: inverse coherent signal power with shape
    (n_periods, n_predictor_channels, n_predictor_channels) namely hx, hy. 
    For MT should always be a (nf, 2, 2)
    :type sigma_s: np.ndarray

    :param periods: periods array should have shape (n_periods), if shape is 
    (1, n_periods) or (n_periods, 1) it will be flattened, raises ValueError
    if given incorrect shape.
    :type periods: np.ndarray

    :param z_err: the error of the impedance tensor elements sqrt(variance), 
    should be shape (n_periods, 2, 2)
    :type z_err: np.ndarray

    :param t_err: the error of the induction vector elements sqrt(variance),
    should be shape (n_periods, 1, 2)
    :type t_err: np.ndarray


    """

    def __init__(
        self, tf=None, sigma_s=None, sigma_e=None, periods=None, z_err=None, t_err=None
    ):
        self.logger = get_mtpy_logger(f"{__name__}.{self.__class__.__name__}")
        self._tf = None
        self._sigma_e = None
        self._sigma_s = None
        self._z_err = None
        self._t_err = None
        self._periods = None
        self._z = None
        self._tipper = None

        self.index = Index()
        self.angles = Angles()

    @property
    def periods(self):
        return self._periods

    @periods.setter
    def periods(self, value):
        value = np.array(value, dtype=np.float)
        if len(value.shape) > 1:
            if len(value.shape) == 2:
                if 1 in value.shape:
                    self.logger.debug(f"Flattened period array of shape {value.shape}")
                    self._periods = value.flatten()
                    return

            msg = f"Input must be 1-D, not shape {value.shape}"
            self.logger.error(msg)
            raise ValueError(msg)
        else:
            self._periods = value

    @property
    def tf(self):
        """

        :return: transfer function
        :rtype: np.ndarray

        """
        return self._tf

    @tf.setter
    def tf(self, value):
        """
        set transfer function

        :param value: np.ndarray with shape (n_periods, n_predicted, n_predictor)
        with a complex data type.
        :type value: np.ndarray((n_periods, n_predicted, n_predictor), dtype=np.complex)

        """

        self._tf = self._validate_input_array(
            value, err_msg="Input must be (n_periods, n_predicted, n_predictor)"
        )

    @property
    def sigma_e(self):
        """

        :return: residual covariance matrix
        :rtype: np.ndarray((n_periods, n_predicted, n_predicted), dtype=np.complex)

        """
        return self._sigma_e

    @sigma_e.setter
    def sigma_e(self, value):
        """
        Set the residual covariance matrix

        :param sigma_e: Residual covariance matrix with shape 
        (n_periods, n_predicted_channels, n_predicted_channels) namely ex, ey, hz
        :type sigma_e: np.ndarray((n_periods, n_predicted, n_predicted), dtype=np.complex)

        """

        self._sigma_e = self._validate_input_array(
            value,
            err_msg="Input must be (n_periods, n_predicted, n_predicted)",
            shapes=True,
        )

    @property
    def sigma_s(self):
        """

        :return: inverse coherent signal power 
        :rtype: np.ndarray((n_periods, n_predictor, n_predictor), dtype=np.complex)

        """
        return self._sigma_s

    @sigma_s.setter
    def sigma_s(self, value):
        """
        Set the residual covariance matrix

        :param sigma_e: inverse coherent signal power 
        (n_periods, n_predicted_channels, n_predicted_channels) namely ex, ey, hz
        :type sigma_e: np.ndarray((n_periods, n_predicted, n_predicted), dtype=np.complex)

        """

        self._sigma_e = self._validate_input_array(
            value,
            err_msg="Input must be (n_periods, n_predictor, n_predictor)",
            shapes=True,
        )

    @property
    def z(self):
        """
        Return impedance
        
        :return: DESCRIPTION
        :rtype: TYPE

        """

        if self.tf is not None:
            return self._calculate_impedance()

    def _validate_input_array(
        self, value, dtype="complex", err_msg=None, test_square=False, test_shape=None
    ):
        """
        Make sure the input array has the correct shape and dtype.
        :param array: DESCRIPTION
        :type array: TYPE
        :return: DESCRIPTION
        :rtype: TYPE

        """
        dtype_dict = {"complex": np.complex_, "real": np.float_, "imaginary": np.float_}
        value = np.array(value, dtype=dtype_dict[dtype])

        if len(value.shape) == 2:
            self.logger.debug(
                "Assuming input is an estimate at a single period"
                + f" reshaping to (1, {value.shape[0]}, {value.shape[1]}"
            )
            value = value.reshape((1, value.shape[0], value.shape[1]))

        elif len(value.shape) == 3:
            if test_square:
                if test_shape is None:
                    if value.shape[1] != value.shape[2]:
                        msg = f"Must be a square matrix not {value.shape}"
                        self.logger.error(msg)
                        raise ValueError(msg)
            if test_shape:
                if len(test_shape) == 2:
                    if not value.shape[1:] != test_shape:
                        msg = (
                            "Input array must have shape "
                            + f"(n_periods, {test_shape[0]}, {test_shape[1]}) "
                            + f"not {value.shape}"
                        )
                        self.logger.error(msg)
                        raise ValueError(msg)
        else:
            msg = f"{err_msg}, not shape {value.shape}"
            self.logger.error(msg)
            raise ValueError(msg)

        if self.periods is not None:
            if value.shape[0] != self.periods.size:
                self.logger.warning(
                    f"Number of periods {self.periods.size} is "
                    + "Not the same as number of tf entries "
                    + f"{value.shape[0]}"
                )
        return value

    def _make_rotation_matrices(
        self, ch1_angle, ch2_angle, ch1_index, ch2_index, rot_angle
    ):
        """
        make rotation matrix
        
        :param angle_01: angle of channel 1
        :type angle_01: float in degrees
        :param angle_02: angle of channel 2
        :type angle_02: float in degrees
        :param rot_angle: angle of rotation
        :type rot_angle: float in degrees
        :param ch1_index: index of channel 1
        
        :return: rotation matrix
        """

        u = np.eye(2, 2)
        u[ch1_index, ch1_index] = np.cos(np.deg2rad(ch1_angle - rot_angle))
        u[ch1_index, ch2_index] = np.sin(np.deg2rad(ch1_angle - rot_angle))
        u[ch2_index, ch1_index] = np.cos(np.deg2rad(ch2_angle - rot_angle))
        u[ch2_index, ch2_index] = np.sin(np.deg2rad(ch2_angle - rot_angle))

        return u

    def calculate_impedance(self, angle=0.0):
        """
        calculate the impedances from the transfer functions
        """

        # check to see if there are actually electric fields in the TFs
        if not hasattr(self, "ex") or not hasattr(self, "ey"):
            msg = (
                "Cannot return apparent resistivity and phase "
                "data because these TFs do not contain electric "
                "fields as a predicted channel."
            )
            self.logger.error(msg)
            raise ValueError(msg)

        # transform the TFs first...
        # build transformation matrix for predictor channels
        #    (horizontal magnetic fields)
        hx_index = self.hx.index
        hy_index = self.hy.index
        u = np.eye(2, 2)
        u[hx_index, hx_index] = np.cos(np.deg2rad(self.hx.azimuth - angle))
        u[hx_index, hy_index] = np.sin(np.deg2rad(self.hx.azimuth - angle))
        u[hy_index, hx_index] = np.cos(np.deg2rad(self.hy.azimuth - angle))
        u[hy_index, hy_index] = np.sin(np.deg2rad(self.hy.azimuth - angle))
        u = np.linalg.inv(u)

        # build transformation matrix for predicted channels (electric fields)
        ex_index = self.ex.index
        ey_index = self.ey.index
        v = np.eye(self.transfer_functions.shape[1], self.transfer_functions.shape[1])
        v[ex_index - 2, ex_index - 2] = np.cos(np.deg2rad(self.ex.azimuth - angle))
        v[ey_index - 2, ex_index - 2] = np.sin(np.deg2rad(self.ex.azimuth - angle))
        v[ex_index - 2, ey_index - 2] = np.cos(np.deg2rad(self.ey.azimuth - angle))
        v[ey_index - 2, ey_index - 2] = np.sin(np.deg2rad(self.ey.azimuth - angle))

        # matrix multiplication...
        rotated_transfer_functions = np.matmul(
            v, np.matmul(self.transfer_functions, u.T)
        )
        rotated_sigma_s = np.matmul(u, np.matmul(self.sigma_s, u.T))
        rotated_sigma_e = np.matmul(v, np.matmul(self.sigma_e, v.T))

        # now pull out the impedance tensor
        z = np.zeros((self.num_freq, 2, 2), dtype=np.complex64)
        z[:, 0, 0] = rotated_transfer_functions[:, ex_index - 2, hx_index]  # Zxx
        z[:, 0, 1] = rotated_transfer_functions[:, ex_index - 2, hy_index]  # Zxy
        z[:, 1, 0] = rotated_transfer_functions[:, ey_index - 2, hx_index]  # Zyx
        z[:, 1, 1] = rotated_transfer_functions[:, ey_index - 2, hy_index]  # Zyy

        # and the variance information
        var = np.zeros((self.num_freq, 2, 2))
        var[:, 0, 0] = np.real(
            rotated_sigma_e[:, ex_index - 2, ex_index - 2]
            * rotated_sigma_s[:, hx_index, hx_index]
        )
        var[:, 0, 1] = np.real(
            rotated_sigma_e[:, ex_index - 2, ex_index - 2]
            * rotated_sigma_s[:, hy_index, hy_index]
        )
        var[:, 1, 0] = np.real(
            rotated_sigma_e[:, ey_index - 2, ey_index - 2]
            * rotated_sigma_s[:, hx_index, hx_index]
        )
        var[:, 1, 1] = np.real(
            rotated_sigma_e[:, ey_index - 2, ey_index - 2]
            * rotated_sigma_s[:, hy_index, hy_index]
        )

        error = np.sqrt(var)

        return z_object

    def calculate_impedance_err(self, angle=0.0):
        """
        Calculate impedance error
        
        """

    def calculate_tippers(self, angle=0.0):
        """
        calculate induction vectors
        """

        # check to see if there is a vertical magnetic field in the TFs
        if self.hz is None:
            raise ZMMError(
                "Cannot return tipper data because the TFs do not "
                "contain the vertical magnetic field as a "
                "predicted channel."
            )

        # transform the TFs first...
        # build transformation matrix for predictor channels
        #    (horizontal magnetic fields)
        hx_index = self.hx.index
        hy_index = self.hy.index
        u = np.eye(2, 2)
        u[hx_index, hx_index] = np.cos(np.deg2rad(self.hx.azimuth - angle))
        u[hx_index, hy_index] = np.sin(np.deg2rad(self.hx.azimuth - angle))
        u[hy_index, hx_index] = np.cos(np.deg2rad(self.hy.azimuth - angle))
        u[hy_index, hy_index] = np.sin(np.deg2rad(self.hy.azimuth - angle))
        u = np.linalg.inv(u)

        # don't need to transform predicated channels (assuming no tilt in Hz)
        hz_index = self.hz.index
        v = np.eye(self.transfer_functions.shape[1], self.transfer_functions.shape[1])

        # matrix multiplication...
        rotated_transfer_functions = np.matmul(
            v, np.matmul(self.transfer_functions, u.T)
        )
        rotated_sigma_s = np.matmul(u, np.matmul(self.sigma_s, u.T))
        rotated_sigma_e = np.matmul(v, np.matmul(self.sigma_e, v.T))

        # now pull out tipper information
        tipper = np.zeros((self.num_freq, 2), dtype=np.complex64)
        tipper[:, 0] = rotated_transfer_functions[:, hz_index - 2, hx_index]  # Tx
        tipper[:, 1] = rotated_transfer_functions[:, hz_index - 2, hy_index]  # Ty

        # and the variance/error information
        var = np.zeros((self.num_freq, 2))
        var[:, 0] = np.real(
            rotated_sigma_e[:, hz_index - 2, hz_index - 2]
            * rotated_sigma_s[:, hx_index, hx_index]
        )  # Tx
        var[:, 1] = np.real(
            rotated_sigma_e[:, hz_index - 2, hz_index - 2]
            * rotated_sigma_s[:, hy_index, hy_index]
        )  # Ty
        error = np.sqrt(var)

        tipper = tipper.reshape((self.num_freq, 1, 2))
        error = error.reshape((self.num_freq, 1, 2))

        tipper_obj = mtz.Tipper(tipper, error, self.frequencies)

        return tipper_obj
