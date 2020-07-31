#!/usr/bin/env python

"""
.. module:: Z
   :synopsis: Deal with MT responses Z and Tipper

.. moduleauthor:: Jared Peacock <jpeacock@usgs.gov>
.. moduleauthor:: Lars Krieger
"""

import cmath
import copy
import math

# =================================================================
import numpy as np

import mtpy.utils.calculator as MTcc
import mtpy.utils.exceptions as MTex
from mtpy.utils.mtpylog import MtPyLog


# get a logger object for this module, using the utility class MtPyLog to
# config the logger
# _logger = MtPyLog.get_mtpy_logger(__name__)


# ==============================================================================
# Resistivity and phase object
# ==============================================================================
class ResPhase(object):
    """
    resistivity and phase container
    """

    def __init__(self, z_array=None, z_err_array=None, freq=None, **kwargs):
        self._logger = MtPyLog.get_mtpy_logger(self.__class__.__name__)

        self._z = z_array
        self._z_err = z_err_array

        self._resistivity = None
        self._phase = None

        self._resistivity_err = None
        self._phase_err = None

        self.freq = freq

        for key in kwargs:
            setattr(self, key, kwargs[key])

    @property
    def resistivity(self):
        return self._resistivity

    @resistivity.setter
    def resistivity(self, res_array):
        self._resistivity = res_array

    @property
    def resistivity_err(self):
        return self._resistivity_err

    @resistivity_err.setter
    def resistivity_err(self, res_err_array):
        self._resistivity_err = res_err_array

    @property
    def phase(self):
        return self._phase

    @phase.setter
    def phase(self, phase_array):
        self._phase = phase_array

    @property
    def phase_err(self):
        return self._phase_err

    @phase_err.setter
    def phase_err(self, phase_err_array):
        self._phase_err = phase_err_array

    def compute_resistivity_phase(self, z_array=None, z_err_array=None,
                                  freq=None):
        """
        compute resistivity and phase from z and z_err
        """

        if z_array is not None:
            self._z = z_array
        if z_err_array is not None:
            self._z_err = z_err_array
        if freq is not None:
            self.freq = freq

        #if self._z is None or self._z_err is None or self.freq is None: #The _z_err can be None!!!
        if self._z is None or self.freq is None:
            raise MT_Z_Error('Values are None, check _z, _z_err, freq')

        self._resistivity = np.apply_along_axis(lambda x: np.abs(x) ** 2 / self.freq * 0.2,
                                                0, self._z)
        self._phase = np.rad2deg(np.angle(self._z))

        self._resistivity_err = np.zeros_like(self._resistivity, dtype=np.float)
        self._phase_err = np.zeros_like(self._phase, dtype=np.float)

        # calculate resistivity and phase
        if self._z_err is not None:
            for idx_f in range(self.freq.size):
                for ii in range(2):
                    for jj in range(2):
#                        r_err, phi_err = MTcc.z_error2r_phi_error(
#                            np.real(self._z[idx_f, ii, jj]),
#                            self._z_err[idx_f, ii, jj],
#                            np.imag(self._z[idx_f, ii, jj]),
#                            self._z_err[idx_f, ii, jj])

                        r_err, phi_err = MTcc.z_error2r_phi_error(
                                self._z[idx_f, ii, jj].real,
                                self._z[idx_f, ii, jj].imag,
                                self._z_err[idx_f, ii, jj])
                        self._resistivity_err[idx_f, ii, jj] = \
                            self._resistivity[idx_f, ii, jj] * r_err
#                        self._resistivity_err[idx_f, ii, jj] = \
#                            0.4 * np.abs(self._z[idx_f, ii, jj]) / \
#                            self.freq[idx_f] * r_err
                        self._phase_err[idx_f, ii, jj] = phi_err

    def set_res_phase(self, res_array, phase_array, freq, res_err_array=None,
                      phase_err_array=None):
        """
        Set values for resistivity (res - in Ohm m) and phase
        (phase - in degrees), including error propagation.


        :param res_array: resistivity array in Ohm-m
        :type res_array: np.ndarray(num_freq, 2, 2)

        :param phase_array: phase array in degrees
        :type phase_array: np.ndarray(num_freq, 2, 2)

        :param freq: frequency array in Hz
        :type freq: np.ndarray(num_freq)

        :param res_err_array: resistivity error array in Ohm-m
        :type res_err_array: np.ndarray(num_freq, 2, 2)

        :param phase_err_array: phase error array in degrees
        :type phase_err_array: np.ndarray(num_freq, 2, 2)


        """

        print('Resetting z and z_err')

        self._resistivity = res_array
        self._phase = phase_array
        self.freq = freq
        self._resistivity_err = res_err_array
        self._phase_err = phase_err_array

        # assert real array:
        if np.linalg.norm(np.imag(res_array)) != 0:
            raise MTex.MTpyError_inputarguments('Error - array "res" is not' + \
                                                'real valued !')

        if np.linalg.norm(np.imag(phase_array)) != 0:
            raise MTex.MTpyError_inputarguments('Error - array "phase" is' + \
                                                'not real valued !')

        abs_z = np.sqrt(5.0 * self.freq * (self.resistivity.T)).T
        self._z = abs_z * np.exp(1j * np.radians(self.phase))

        self._z_err = np.zeros_like(self._z, dtype=np.float)
        # ---------------------------
        # error propagation:
        if self._resistivity_err is None or self._phase_err is None:
            return

        for idx_f in range(self.freq.shape[0]):
            for ii in range(2):
                for jj in range(2):
                    abs_z = np.sqrt(5 * self.freq[idx_f] * \
                                    self.resistivity[idx_f, ii, jj])
                    rel_error_res = self.resistivity_err[idx_f, ii, jj] / \
                                    self.resistivity[idx_f, ii, jj]
                    # relative error varies by a factor of 0.5, which is the
                    # exponent in the relation between them:
                    abs_z_error = 0.5 * abs_z * rel_error_res

                    self._z_err[idx_f, ii, jj] = max(MTcc.propagate_error_polar2rect(
                        abs_z,
                        abs_z_error,
                        self.phase[idx_f, ii, jj],
                        self.phase_err[idx_f, ii, jj]))

    @property
    def res_xx(self):
        return self._resistivity[:, 0, 0]

    @property
    def res_xy(self):
        return self._resistivity[:, 0, 1]

    @property
    def res_yx(self):
        return self._resistivity[:, 1, 0]

    @property
    def res_yy(self):
        return self._resistivity[:, 1, 1]

    @property
    def phase_xx(self):
        return self._phase[:, 0, 0]

    @property
    def phase_xy(self):
        return self._phase[:, 0, 1]

    @property
    def phase_yx(self):
        return self._phase[:, 1, 0]

    @property
    def phase_yy(self):
        return self._phase[:, 1, 1]

    @property
    def res_err_xx(self):
        return self._resistivity_err[:, 0, 0]

    @property
    def res_err_xy(self):
        return self._resistivity_err[:, 0, 1]

    @property
    def res_err_yx(self):
        return self._resistivity_err[:, 1, 0]

    @property
    def res_err_yy(self):
        return self._resistivity_err[:, 1, 1]

    @property
    def phase_err_xx(self):
        return self._phase_err[:, 0, 0]

    @property
    def phase_err_xy(self):
        return self._phase_err[:, 0, 1]

    @property
    def phase_err_yx(self):
        return self._phase_err[:, 1, 0]

    @property
    def phase_err_yy(self):
        return self._phase_err[:, 1, 1]

    # calculate determinant values
    @property
    def _zdet(self):
        return np.array([np.linalg.det(zz) ** .5 for zz in self._z])

    @property
    def _zdet_var(self):
        if self._z_err is not None:
            return np.array([abs(np.linalg.det(zzv)) ** .5 for zzv in self._z_err])
        else:
            return np.ones_like(self._zdet, dtype=np.float)

    @property
    def phase_det(self):
        return np.arctan2(self._zdet.imag, self._zdet.real) * (180 / np.pi)

    @property
    def phase_det_err(self):
        return np.arcsin(self._zdet_var / abs(self._zdet)) * (180 / np.pi)

    @property
    def res_det(self):
        return 0.2 * (1. / self.freq) * abs(self._zdet) ** 2

    @property
    def res_det_err(self):
        return 0.2 * (1. / self.freq) * np.abs(self._zdet + self._zdet_var) ** 2 - self.res_det


# ==============================================================================
# Impedance Tensor Class
# ==============================================================================
class Z(ResPhase):
    """
    Z class - generates an impedance tensor (Z) object.

    Z is a complex array of the form (n_freq, 2, 2),
    with indices in the following order:

        - Zxx: (0,0)
        - Zxy: (0,1)
        - Zyx: (1,0)
        - Zyy: (1,1)

    All errors are given as standard deviations (sqrt(VAR))

    :param z_array: array containing complex impedance values
    :type z_array: numpy.ndarray(n_freq, 2, 2)


    :param z_err_array: array containing error values (standard deviation)
                        of impedance tensor elements
    :type z_err_array: numpy.ndarray(n_freq, 2, 2)

    :param freq: array of frequency values corresponding to impedance tensor
                 elements.
    :type freq: np.ndarray(n_freq)


    =============== ===========================================================
    Attributes      Description
    =============== ===========================================================
    freq             array of frequencies corresponding to elements of z
    rotation_angle   angle of which data is rotated by
    z                impedance tensor
    z_err             estimated errors of impedance tensor
    resistivity      apparent resisitivity estimated from z in Ohm-m
    resistivity_err  apparent resisitivity error
    phase            impedance phase (deg)
    phase_err        error in impedance phase
    =============== ===========================================================

    =================== =======================================================
    Methods             Description
    =================== =======================================================
    det                  calculates determinant of z with errors
    invariants           calculates the invariants of z
    inverse              calculates the inverse of z
    remove_distortion    removes distortion given a distortion matrix
    remove_ss            removes static shift by assumin Z = S * Z_0
    norm                 calculates the norm of Z
    only1d               zeros diagonal components and computes
                         the absolute valued mean of the off-diagonal
                         components.
    only2d               zeros diagonal components
    res_phase            computes resistivity and phase
    rotate               rotates z positive clockwise, angle assumes
                         North is 0.
    set_res_phase        recalculates z and z_err, needs attribute freq
    skew                 calculates the invariant skew (off diagonal trace)
    trace                calculates the trace of z
    =================== =======================================================

    :Example: ::

        >>> import mtpy.core.z as mtz
        >>> import numpy as np
        >>> z_test = np.array([[0+0j, 1+1j], [-1-1j, 0+0j]])
        >>> z_object = mtz.Z(z_array=z_test, freq=[1])
        >>> z_object.rotate(45)
        >>> z_object.resistivity


    """

    def __init__(self, z_array=None, z_err_array=None, freq=None):
        """
        Initialise an instance of the Z class.

        :param z_array: array containing complex impedance values
        :type z_array: numpy.ndarray(n_freq, 2, 2)

        :param z_err_array: array containing error values (standard deviation)
                            of impedance tensor elements
        :type z_err_array: numpy.ndarray(n_freq, 2, 2)

        :param freq: array of frequency values corresponding to impedance
                     tensor elements.
        :type freq: np.ndarray(n_freq)

        Initialises the attributes with None
        """

        ResPhase.__init__(self)

        self._z = z_array
        self._z_err = z_err_array
        self._freq = freq

        if z_array is not None:
            if len(z_array.shape) == 2 and z_array.shape == (2, 2):
                if z_array.dtype in ['complex', 'float', 'int']:
                    self._z = np.zeros((1, 2, 2), 'complex')
                    self._z[0] = z_array

        if z_err_array is not None:
            if len(z_err_array.shape) == 2 and z_err_array.shape == (2, 2):
                if z_err_array.dtype in ['complex', 'float', 'int']:
                    self._z_err = np.zeros((1, 2, 2), 'complex')
                    self._z_err[0] = z_err_array

        self.rotation_angle = 0.
        if self._z is not None:
            self.rotation_angle = np.zeros((len(self._z)))

        if self._z is not None:
            self.compute_resistivity_phase()

    # ---frequency-------------------------------------------------------------
    @property
    def freq(self):
        """
        Frequencies for each impedance tensor element

        Units are Hz.
        """
        return self._freq

    @freq.setter
    def freq(self, freq_arr):
        """
        Set the array of freq.

        :param freq_arr: array of frequnecies (Hz)
        :type freq_arr: np.ndarray
        """

        if freq_arr is not None:
            self._freq = np.array(freq_arr)
        else:
            return None

        if self.z is not None:
            if len(self.z.shape) == 3:
                if self._freq.size != len(self.z):
                    self._logger.warn('length of freq list/array not correct'
                                      '({0} instead of {1})'.format(self._freq.size,
                                                                    len(self.z)))
                    return
                else:
                    try:
                        self.compute_resistivity_phase()
                    except IndexError:
                        print('Need to input frequency array')

    # ----impedance tensor -----------------------------------------------------
    @property
    def z(self):
        """
        Impedance tensor

        np.ndarray(nfreq, 2, 2)
        """
        return self._z

    @z.setter
    def z(self, z_array):
        """
        Set the attribute 'z'.


        :param z_array: complex impedance tensor array
        :type z_array: np.ndarray(nfreq, 2, 2)

        Test for shape, but no test for consistency!

        Nulling the rotation_angle
        """

        try:
            if len(z_array.shape) == 3 and z_array.shape[1:3] == (2, 2):
                if z_array.dtype in ['complex', 'float', 'int']:
                    self._z = z_array
        except IndexError:
            try:
                if len(z_array.shape) == 2 and z_array.shape == (2, 2):
                    if z_array.dtype in ['complex', 'float', 'int']:
                        self._z = np.zeros((1, 2, 2), 'complex')
                        self._z[0] = z_array
            except IndexError:
                self._logger.error('provided Z array does not have correct dimensions- Z unchanged')

        if isinstance(self.rotation_angle, float):
            self.rotation_angle = np.repeat(self.rotation_angle,
                                            len(self._z))

        # for consistency recalculate resistivity and phase
        if self._z is not None and self._z_err is not None:
            try:
                self.compute_resistivity_phase()
            except IndexError:
                self._logger.error('Need to input frequency array')

    # ----impedance error-----------------------------------------------------
    @property
    def z_err(self):
        return self._z_err

    @z_err.setter
    def z_err(self, z_err_array):
        """
        Set the attribute z_err

        :param z_err_array: error of impedance tensor array as standard
                            deviation
        :type z_err_array: np.ndarray(nfreq, 2, 2)
        """
        if z_err_array.shape != self.z.shape:
            self._logger.warn('z_err_array shape {0} is not same shape as z {1}'.format(
                z_err_array.shape, self.z.shape))
        self._z_err = z_err_array

        # for consistency recalculate resistivity and phase
        if self._z_err is not None and self._z is not None:
            try:
                self.compute_resistivity_phase()
            except IndexError:
                self._logger.error('Need to input frequency array')

    @property
    def inverse(self):
        """
            Return the inverse of Z.

            (no error propagtaion included yet)

        """

        if self.z is None:
            self._logger.warn('z array is "None" - I cannot invert that')
            return

        inverse = copy.copy(self.z)
        for idx_f in range(len(inverse)):
            try:
                inverse[idx_f, :, :] = np.array((np.matrix(self.z[idx_f, :, :])).I)
            except:
                raise MTex.MTpyError_Z('The {0}ith impedance'.format(idx_f + 1) + \
                                       'tensor cannot be inverted')

        return inverse

    def rotate(self, alpha):
        """
        Rotate Z array by angle alpha.

        Rotation angle must be given in degrees. All angles are referenced
        to geographic North, positive in clockwise direction.
        (Mathematically negative!)

        In non-rotated state, X refs to North and Y to East direction.

        Updates the attributes
            - *z*
            - *z_err*
            - *zrot*
            - *resistivity*
            - *phase*
            - *resistivity_err*
            - *phase_err*

        """

        if self.z is None:
            self._logger.warn('Z array is "None" - I cannot rotate that')
            return

        # check for iterable list/set of angles - if so, it must have length
        # 1 or same as len(tipper):
        if np.iterable(alpha) == 0:
            try:
                degreeangle = float(alpha % 360)
            except ValueError:
                self._logger.error('"Angle" must be a valid number (in degrees)')
                return

            # make an n long list of identical angles
            lo_angles = [degreeangle for ii in self.z]
        else:
            if len(alpha) == 1:
                try:
                    degreeangle = float(alpha % 360)
                except ValueError:
                    self._logger.error('"Angle" must be a valid number (in degrees)')
                    return
                # make an n long list of identical angles
                lo_angles = [degreeangle for ii in self.z]
            else:
                try:
                    lo_angles = [float(ii % 360) for ii in alpha]
                except ValueError:
                    self._logger.error('"Angles" must be valid numbers (in degrees)')
                    return

        self.rotation_angle = np.array([(oldangle + lo_angles[ii]) % 360
                                        for ii, oldangle in enumerate(self.rotation_angle)])

        if len(lo_angles) != len(self.z):
            self._logger.warn('Wrong number of "angles" - I need {0}'.format(len(self.z)))
            # self.rotation_angle = 0.
            return

        z_rot = copy.copy(self.z)
        z_err_rot = copy.copy(self.z_err)

        for idx_freq in range(len(self.z)):

            angle = lo_angles[idx_freq]
            if np.isnan(angle):
                angle = 0.

            if self.z_err is not None:
                z_rot[idx_freq], z_err_rot[idx_freq] = \
                    MTcc.rotatematrix_incl_errors(self.z[idx_freq, :, :],
                                                  angle,
                                                  self.z_err[idx_freq, :, :])
            else:
                z_rot[idx_freq], z_err_rot = \
                    MTcc.rotatematrix_incl_errors(self.z[idx_freq, :, :],
                                                  angle)

        self.z = z_rot
        if self.z_err is not None:
            self.z_err = z_err_rot

        # for consistency recalculate resistivity and phase
        self.compute_resistivity_phase()

    def remove_ss(self, reduce_res_factor_x=1., reduce_res_factor_y=1.):
        """
        Remove the static shift by providing the respective correction factors
        for the resistivity in the x and y components.
        (Factors can be determined by using the "Analysis" module for the
        impedance tensor)

        Assume the original observed tensor Z is built by a static shift S
        and an unperturbated "correct" Z0 :

             * Z = S * Z0

        therefore the correct Z will be :
            * Z0 = S^(-1) * Z


        :param reduce_res_factor_x: static shift factor to be applied to x
                                    components (ie z[:, 0, :]).  This is
                                    assumed to be in resistivity scale
        :type reduce_res_factor_x: float or iterable list or array

        :param reduce_res_factor_y: static shift factor to be applied to y
                                    components (ie z[:, 1, :]).  This is
                                    assumed to be in resistivity scale
        :type reduce_res_factor_y: float or iterable list or array

        :returns: static shift matrix,
        :rtype: np.ndarray ((2, 2))

        :returns: corrected Z
        :rtype: mtpy.core.z.Z

        .. note:: The factors are in resistivity scale, so the
                  entries of  the matrix "S" need to be given by their
                  square-roots!

        """

        # check for iterable list/set of reduce_res_factor_x - if so, it must
        # have length 1 or same as len(z):
        if np.iterable(reduce_res_factor_x) == 0:
            try:
                x_factor = float(reduce_res_factor_x)
            except ValueError:
                self._logger.error('reduce_res_factor_x must be a valid numbers')
                return

            lo_x_factors = np.repeat(x_factor, len(self.z))
        elif len(reduce_res_factor_x) == 1:
            try:
                x_factor = float(reduce_res_factor_x)
            except ValueError:
                self._logger.error('reduce_res_factor_x must be a valid numbers')
                return
            lo_x_factors = np.repeat(x_factor, len(self.z))
        else:
            try:
                lo_x_factors = np.repeat(x_factor,
                                         len(reduce_res_factor_x))
            except ValueError:
                self._logger.error('"reduce_res_factor_x" must be valid numbers')
                return

        if len(lo_x_factors) != len(self.z):
            self._logger.error('Wrong number Number of reduce_res_factor_x - need {0}'.format(len(self.z)))
            return

        # check for iterable list/set of reduce_res_factor_y - if so,
        # it must have length 1 or same as len(z):
        if np.iterable(reduce_res_factor_y) == 0:
            try:
                y_factor = float(reduce_res_factor_y)
            except ValueError:
                self._logger.error('"reduce_res_factor_y" must be a valid numbers')
                return

            lo_y_factors = np.repeat(y_factor, len(self.z))
        elif len(reduce_res_factor_y) == 1:
            try:
                y_factor = float(reduce_res_factor_y)
            except ValueError:
                self._logger.error('"reduce_res_factor_y" must be a valid numbers')
                return
            lo_y_factors = np.repeat(y_factor, len(self.z))
        else:
            try:
                lo_y_factors = np.repeat(y_factor,
                                         len(reduce_res_factor_y))
            except ValueError:
                self._logger.error('"reduce_res_factor_y" must be valid numbers')
                return

        if len(lo_y_factors) != len(self.z):
            self._logger.error('Wrong number Number of "reduce_res_factor_y"' + \
                               '- need {0} '.format(len(self.z)))
            return

        z_corrected = copy.copy(self.z)
        static_shift = np.zeros((len(self.z), 2, 2))

        for idx_f in range(len(self.z)):
            # correct for x-direction
            z_corrected[idx_f, 0, :] = self.z[idx_f, 0, :] / \
                                       np.sqrt(lo_x_factors[idx_f])
            # correct for y-direction
            z_corrected[idx_f, 1, :] = self.z[idx_f, 1, :] / \
                                       np.sqrt(lo_y_factors[idx_f])
            # make static shift array
            static_shift[idx_f, 0, 0] = np.sqrt(lo_x_factors[idx_f])
            static_shift[idx_f, 1, 1] = np.sqrt(lo_y_factors[idx_f])

        return static_shift, z_corrected

    def remove_distortion(self, distortion_tensor, distortion_err_tensor=None):
        """
        Remove distortion D form an observed impedance tensor Z to obtain
        the uperturbed "correct" Z0:
        Z = D * Z0

        Propagation of errors/uncertainties included


        :param distortion_tensor: real distortion tensor as a 2x2
        :type distortion_tensor: np.ndarray(2, 2, dtype=real)


        :param distortion_err_tensor: default is None
        :type distortion_err_tensor: np.ndarray(2, 2, dtype=real),

			:returns: input distortion tensor
        :rtype: np.ndarray(2, 2, dtype='real')

			:returns: impedance tensor with distorion removed
        :rtype: np.ndarray(num_freq, 2, 2, dtype='complex')


			:returns: impedance tensor error after distortion is removed
        :rtype: np.ndarray(num_freq, 2, 2, dtype='complex')


    		:Example: ::

    			>>> import mtpy.core.z as mtz
    			>>> distortion = np.array([[1.2, .5],[.35, 2.1]])
    			>>> d, new_z, new_z_err = z_obj.remove_distortion(distortion)
        """

        if distortion_err_tensor is None:
            distortion_err_tensor = np.zeros_like(distortion_tensor)
        # for all freq, calculate D.Inverse, then obtain Z0 = D.I * Z
        try:
            if not (len(distortion_tensor.shape) in [2, 3]) and \
                    (len(distortion_err_tensor.shape) in [2, 3]):
                raise ValueError('Shape not the same')
            if len(distortion_tensor.shape) == 3 or \
                            len(distortion_err_tensor.shape) == 3:
                self._logger.info('Distortion is not time-dependent - take only first' + \
                                  'of given distortion tensors')
                try:
                    distortion_tensor = distortion_tensor[0]
                    distortion_err_tensor = distortion_err_tensor[0]
                except IndexError:
                    raise ValueError('distortion tensor the wrong shape')

            if not (distortion_tensor.shape == (2, 2)) and \
                    (distortion_err_tensor.shape == (2, 2)):
                raise ValueError('Shape not the same')

            distortion_tensor = np.matrix(np.real(distortion_tensor))

        except ValueError:
            raise MTex.MTpyError_Z('The array provided is not a proper' + \
                                   'distortion tensor')

        try:
            DI = distortion_tensor.I
        except np.linalg.LinAlgError:
            raise MTex.MTpyError_Z('The provided distortion tensor is' + \
                                   'singular - I cannot invert that!')

        # propagation of errors (using 1-norm) - step 1 - inversion of D:
        DI_err = np.zeros_like(distortion_err_tensor)

        # todo :include error on  determinant!!
        # D_det = np.linalg.det(distortion_tensor)

        dummy, DI_err = MTcc.invertmatrix_incl_errors(distortion_tensor,
                                                      distortion_err_tensor)

        # propagation of errors - step 2 - product of D.inverse and Z;
        # D.I * Z, making it 4 summands for each component:
        z_corrected = np.zeros_like(self.z)
        z_corrected_err = np.zeros_like(self.z_err)

        for idx_f in range(len(self.z)):
            z_corrected[idx_f] = np.array(np.dot(DI, np.matrix(self.z[idx_f])))
            for ii in range(2):
                for jj in range(2):
                    z_corrected_err[idx_f, ii, jj] = np.sum(np.abs(
                        np.array([DI_err[ii, 0] * \
                                  self.z[idx_f, 0, jj], \
                                  DI[ii, 0] * \
                                  self.z_err[idx_f, 0, jj], \
                                  DI_err[ii, 1] * \
                                  self.z[idx_f, 1, jj], \
                                  DI[ii, 1] * \
                                  self.z_err[idx_f, 1, jj]])))

        return distortion_tensor, z_corrected, z_corrected_err

    def _compute_det_variance(self):
        """
        compute the variance of the determinant of Z, 
        """

    @property
    def only_1d(self):
        """
        Return Z in 1D form.

        If Z is not 1D per se, the diagonal elements are set to zero,
        the off-diagonal elements keep their signs, but their absolute
        is set to the mean of the original Z off-diagonal absolutes.
        """

        z1d = copy.copy(self.z)

        for ii in range(len(z1d)):
            z1d[ii, 0, 0] = 0
            z1d[ii, 1, 1] = 0
            sign01 = np.sign(z1d[ii, 0, 1])
            sign10 = np.sign(z1d[ii, 1, 0])
            mean1d = 0.5 * (z1d[ii, 1, 0] + z1d[ii, 0, 1])
            z1d[ii, 0, 1] = sign01 * mean1d
            z1d[ii, 1, 0] = sign10 * mean1d

        return z1d

    @property
    def only_2d(self):
        """
        Return Z in 2D form.

        If Z is not 2D per se, the diagonal elements are set to zero.
        """

        z2d = copy.copy(self.z)

        for ii in range(len(z2d)):
            z2d[ii, 0, 0] = 0
            z2d[ii, 1, 1] = 0

        return z2d

    @property
    def trace(self):
        """
        Return the trace of Z

        :returns: Trace(z)
        :rtype: np.ndarray(nfreq, 2, 2)

        """

        tr = np.array([np.trace(ii) for ii in self.z])

        return tr

    @property
    def trace_err(self):
        """
        Return the trace of Z

        :returns: Trace(z)
        :rtype: np.ndarray(nfreq, 2, 2)

        """

        tr_err = None
        if self.z_err is not None:
            tr_err = np.zeros_like(self.trace, dtype=np.float)
            tr_err[:] = self.z_err[:, 0, 0] + self.z_err[:, 1, 1]

        return tr_err

    @property
    def skew(self):
        """
        Returns the skew of Z as defined by Z[0, 1] + Z[1, 0]

        .. note:: This is not the MT skew, but simply the linear algebra skew


        :returns: skew
        :rtype: np.ndarray(nfreq, 2, 2)
        """

        skew = np.array([ii[0, 1] - ii[1, 0] for ii in self.z])

        return skew

    @property
    def skew_err(self):
        """
        Returns the skew error of Z as defined by Z_err[0, 1] + Z_err[1, 0]

        .. note:: This is not the MT skew, but simply the linear algebra skew

        :returns: skew_err
        :rtype: np.ndarray(nfreq, 2, 2)
        """

        skew_err = None
        if self.z_err is not None:
            skew_err = np.zeros_like(self.skew, dtype=np.float)
            skew_err[:] = self.z_err[:, 0, 1] + self.z_err[:, 1, 0]

        return skew_err

    @property
    def det(self):
        """
        Return the determinant of Z

        :returns: det_Z
        :rtype: np.ndarray(nfreq)
        """

        det_Z = np.array([np.linalg.det(ii) for ii in self.z])

        return det_Z

    @property
    def det_err(self):
        """
        Return the determinant of Z error

        :returns: det_Z_err
        :rtype: np.ndarray(nfreq)
        """
        det_Z_err = None
        if self.z_err is not None:
            det_Z_err = np.zeros_like(self.det, dtype=np.float)
            # components of the impedance tensor are not independent variables
            # so can't use standard error propagation
            # calculate manually:
            # difference of determinant of z + z_err and z - z_err then divide by 2
            det_Z_err[:] = np.abs(np.linalg.det(self.z + self.z_err) - \
                                  np.linalg.det(self.z - self.z_err)) / 2.
        # det_Z_err[:] = np.abs(self.z[:, 1, 1] * self.z_err[:, 0, 0]) +\
        #                           np.abs(self.z[:, 0, 0] * self.z_err[:, 1, 1]) +\
        #                           np.abs(self.z[:, 0, 1] * self.z_err[:, 1, 0]) +\
        #                           np.abs(self.z[:, 1, 0] * self.z_err[:, 0, 1])

        return det_Z_err

    @property
    def norm(self):
        """
        Return the 2-/Frobenius-norm of Z

        :returns: norm
        :rtype: np.ndarray(nfreq)
        """

        norm = np.array([np.linalg.norm(ii) for ii in self.z])

        return norm

    @property
    def norm_err(self):
        """
        Return the 2-/Frobenius-norm of Z  error

        :returns: norm_err
        :rtype: np.ndarray(nfreq)
        """
        norm_err = None

        if self.z_err is not None:
            norm_err = np.zeros_like(self.norm, dtype=np.float)
            for idx, z_tmp in enumerate(self.z):
                value = self.norm[idx]
                error_matrix = self.z_err[idx]
                radicand = 0.
                for ii in range(2):
                    for jj in range(2):
                        radicand += (error_matrix[ii, jj] * \
                                     np.real(z_tmp[ii, jj])) ** 2
                        radicand += (error_matrix[ii, jj] * \
                                     np.imag(z_tmp[ii, jj])) ** 2

                norm_err[idx] = 1. / value * np.sqrt(radicand)

        return norm_err

    @property
    def invariants(self):
        """
        Return a dictionary of Z-invariants.

        Contains
    		-----------
    			* z1
    			* det
    			* det_real
    			* det_imag
    			* trace
    			* skew
    			* norm
    			* lambda_plus/minus,
    			* sigma_plus/minus
        """

        invariants_dict = {}

        z1 = (self.z[:, 0, 1] - self.z[:, 1, 0]) / 2.
        invariants_dict['z1'] = z1

        invariants_dict['det'] = self.det[0]

        det_real = np.array([np.linalg.det(ii) for ii in np.real(self.z)])
        invariants_dict['det_real'] = det_real

        det_imag = np.array([np.linalg.det(ii) for ii in np.imag(self.z)])
        invariants_dict['det_imag'] = det_imag

        invariants_dict['trace'] = self.trace

        invariants_dict['skew'] = self.skew

        invariants_dict['norm'] = self.norm

        invariants_dict['lambda_plus'] = z1 + np.sqrt(z1 * z1 / self.det)

        invariants_dict['lambda_minus'] = z1 - np.sqrt(z1 * z1 / self.det)

        invariants_dict['sigma_plus'] = 0.5 * self.norm ** 2 + \
                                        np.sqrt(0.25 * self.norm ** 4) + \
                                        np.abs(self.det ** 2)

        invariants_dict['sigma_minus'] = 0.5 * self.norm ** 2 - \
                                         np.sqrt(0.25 * self.norm ** 4) + \
                                         np.abs(self.det ** 2)

        return invariants_dict


# ==============================================================================
# errors
# ==============================================================================
class MT_Z_Error(Exception):
    pass


# ======================================================================
#                               TIPPER
# ======================================================================

class Tipper(object):
    """
    Tipper class --> generates a Tipper-object.

    Errors are given as standard deviations (sqrt(VAR))

    :param tipper_array: tipper array in the shape of [Tx, Ty]
                         *default* is None
    :type tipper_array: np.ndarray((nf, 1, 2), dtype='complex')


    :param tipper_err_array: array of estimated tipper errors
                             in the shape of [Tx, Ty].
                             Must be the same shape as tipper_array.
                             *default* is None
    :type tipper_err_array: np.ndarray((nf, 1, 2))


    :param freq: array of frequencies corresponding to the tipper elements.
                 Must be same length as tipper_array.
                 *default* is None
    :type freq: np.ndarray(nf)


    =============== ===========================================================
    Attributes      Description
    =============== ===========================================================
    freq            array of frequencies corresponding to elements of z
    rotation_angle  angle of which data is rotated by

    tipper          tipper array
    tipper_err       tipper error array
    =============== ===========================================================

    =============== ===========================================================
    Methods         Description
    =============== ===========================================================
    mag_direction   computes magnitude and direction of real and imaginary
                    induction arrows.
    amp_phase       computes amplitude and phase of Tx and Ty.
    rotate          rotates the data by the given angle
    =============== ===========================================================
    """

    def __init__(self, tipper_array=None, tipper_err_array=None,
                 freq=None):
        """
        initialize
        """
        self._logger = MtPyLog.get_mtpy_logger(self.__class__.__name__)
        self._tipper = tipper_array
        self._tipper_err = tipper_err_array
        self._freq = freq

        self.rotation_angle = 0.
        if self.tipper is not None:
            self.rotation_angle = np.zeros((len(self.tipper)))

        self._amplitude = None
        self._amplitude_err = None
        self._phase = None
        self._phase_err = None

        self._mag_real = None
        self._mag_imag = None
        self._angle_real = None
        self._angle_imag = None
        self._mag_err = None
        self._angle_err = None

        if self._tipper is not None and self._freq is not None:
            self.compute_amp_phase()
            self.compute_mag_direction()

    # ==========================================================================
    # Define get/set and properties
    # ==========================================================================
    # ----freq----------------------------------------------------------
    @property
    def freq(self):
        return self._freq

    @freq.setter
    def freq(self, freq_arr):
        """
        Set the array of freq.

        :param freq_arr: array of frequnecies (Hz)
        :type freq_arr: np.ndarray(num_frequencies)
        """
        if freq_arr is not None:
            self._freq = np.array(freq_arr)

        if self._freq.size is not len(self.tipper):
            self._logger.info('length of freq list/array not correct' + \
                              ' (%ii instead of %ii)' % (self._freq.size, len(self.tipper)))
            return

        # for consistency recalculate amplitude and phase
        self.compute_amp_phase()

    # ---tipper--------------------------------------------------------------
    @property
    def tipper(self):
        return self._tipper

    @tipper.setter
    def tipper(self, tipper_array):
        """
        Set the attribute *tipper*

        :param tipper_array: tipper array in the shape of [Tx, Ty]
                             *default* is None
        :type tipper_array: np.ndarray((nf, 1, 2), dtype='complex')
        """

        # check to see if the new tipper array is the same shape as the old
        if self._tipper is not None and self._tipper.shape != tipper_array.shape:
            raise MT_Z_Error('Shape of new "tipper" array does not match old' + \
                             'new shape {0} != old shape {1}'.format(tipper_array.shape,
                                                                     self._tipper.shape) + \
                             '\n***Make new Tipper object***')
        if tipper_array is not None:
            if len(tipper_array.shape) == 3 and tipper_array.shape[1:3] == (1, 2):
                if tipper_array.dtype in ['complex', 'float', 'int']:
                    self._tipper = tipper_array

        # neeed to set the rotation angle such that it is an array
        if self.rotation_angle is float:
            self.rotation_angle = np.repeat(self.rotation_angle,
                                            len(self._tipper))

            # for consistency recalculate mag and angle
        self.compute_mag_direction()

        # for consistency recalculate amplitude and phase
        self.compute_amp_phase()

    # ----tipper error---------------------------------------------------------
    @property
    def tipper_err(self):
        return self._tipper_err

    @tipper_err.setter
    def tipper_err(self, tipper_err_array):
        """
        Set the attribute *tipper_err*.

        :param tipper_err_array: array of estimated tipper errors
                                 in the shape of [Tx, Ty].
                                 Must be the same shape as tipper_array.
                                 *default* is None
        :type tipper_err_array: np.ndarray((nf, 1, 2))
        """
        if self.tipper_err is not None and \
                (self._tipper_err.shape != tipper_err_array.shape):
            raise MT_Z_Error('Shape of new "tipper_err" array does not match old' + \
                             'new shape {0} != old shape {1}'.format(tipper_err_array.shape),
                             self._tipper_err.shape)

        # make sure the input array is of required shape
        if tipper_err_array is not None:
            if len(tipper_err_array.shape) == 3 and \
                            tipper_err_array.shape[1:3] == (1, 2):
                if tipper_err_array.dtype in ['float', 'int']:
                    self._tipper_err = tipper_err_array

                    assert self._tipper_err.shape == self._tipper.shape

        # for consistency recalculate mag and angle
        self.compute_mag_direction()

        # for consistency recalculate amplitude and phase
        self.compute_amp_phase()

    # ----amplitude and phase
    def compute_amp_phase(self):
        """
        Sets attributes:
    			* *amplitude*
    			* *phase*
    			* *amplitude_err*
    			* *phase_err*

        values for resistivity are in in Ohm m and phase in degrees.
        """

        if self.tipper is None:
            # logging.error( 'tipper array is None - cannot calculate rho/phi')
            # print 'tipper array is None - cannot calculate rho/phi'
            return None

        self._amplitude_err = None
        self._phase_err = None
        if self.tipper_err is not None:
            self._amplitude_err = np.zeros(self.tipper_err.shape)
            self._phase_err = np.zeros(self.tipper_err.shape)

        self._amplitude = np.abs(self.tipper)
        self._phase = np.rad2deg(np.angle(self.tipper))

        if self.tipper_err is not None:
            for idx_f in range(len(self.tipper)):
                for jj in range(2):
                    if self.tipper_err is not None:
                        if type(self.tipper) == np.ma.core.MaskedArray:
                            if self.tipper.mask[idx_f, 0, jj]:
                                continue
                        r_err, phi_err = MTcc.propagate_error_rect2polar(
                            np.real(self.tipper[idx_f, 0, jj]),
                            self.tipper_err[idx_f, 0, jj],
                            np.imag(self.tipper[idx_f, 0, jj]),
                            self.tipper_err[idx_f, 0, jj])

                        self.amplitude_err[idx_f, 0, jj] = r_err
                        self._phase_err[idx_f, 0, jj] = phi_err

    def set_amp_phase(self, r_array, phi_array):
        """
        Set values for amplitude(r) and argument (phi - in degrees).

        Updates the attributes:
    			* tipper
    			* tipper_err

        """

        if self.tipper is not None:

            tipper_new = copy.copy(self.tipper)

            if self.tipper.shape != r_array.shape:
                self._logger.error('Error - shape of "r" array does not match shape of ' + \
                                   'tipper array: %s ; %s' % (str(r_array.shape),
                                                              str(self.tipper.shape)))
                return

            if self.tipper.shape != phi_array.shape:
                self._logger.error('Error - shape of "phi" array does not match shape of ' + \
                                   'tipper array: %s ; %s' % (str(phi_array.shape),
                                                              str(self.tipper.shape)))
                return
        else:

            tipper_new = np.zeros(r_array.shape, 'complex')

            if r_array.shape != phi_array.shape:
                self._logger.error('Error - shape of "phi" array does not match shape ' + \
                                   'of "r" array: %s ; %s' % (str(phi_array.shape),
                                                              str(r_array.shape)))
                return

        # assert real array:
        if np.linalg.norm(np.imag(r_array)) != 0:
            self._logger.error('Error - array "r" is not real valued !')
            return
        if np.linalg.norm(np.imag(phi_array)) != 0:
            self._logger.error('Error - array "phi" is not real valued !')
            return

        for idx_f in range(len(r_array)):
            for jj in range(2):
                tipper_new[idx_f, 0, jj] = cmath.rect(r_array[idx_f, 0, jj],
                                                      math.radians(phi_array[idx_f, 0, jj]))

        self.tipper = tipper_new

        # for consistency recalculate amplitude and phase
        self.compute_amp_phase()
        self.compute_mag_direction()

    # ---------------------------------
    # properties
    @property
    def amplitude(self):
        return self._amplitude

    @property
    def phase(self):
        return self._phase

    @property
    def amplitude_err(self):
        return self._amplitude_err

    @property
    def phase_err(self):
        return self._phase_err

    # ----magnitude and direction----------------------------------------------
    def compute_mag_direction(self):
        """
        computes the magnitude and direction of the real and imaginary
        induction vectors.
        """

        if self.tipper is None:
            return None
        self._mag_real = np.sqrt(self.tipper[:, 0, 0].real ** 2 + \
                                 self.tipper[:, 0, 1].real ** 2)
        self._mag_imag = np.sqrt(self.tipper[:, 0, 0].imag ** 2 +
                                 self.tipper[:, 0, 1].imag ** 2)

        self._mag_err = None
        self._angle_err = None
        # get the angle, need to make both parts negative to get it into the
        # parkinson convention where the arrows point towards the conductor

        self._angle_real = np.rad2deg(np.arctan2(-self.tipper[:, 0, 1].real,
                                                 -self.tipper[:, 0, 0].real))

        self._angle_imag = np.rad2deg(np.arctan2(-self.tipper[:, 0, 1].imag,
                                                 -self.tipper[:, 0, 0].imag))

        ## estimate error: THIS MAYBE A HACK
        if self.tipper_err is not None:
            self._mag_err = np.sqrt(self.tipper_err[:, 0, 0] ** 2 + \
                                    self.tipper_err[:, 0, 1] ** 2)
            self._angle_err = np.rad2deg(np.arctan2(self.tipper_err[:, 0, 0],
                                                    self.tipper_err[:, 0, 1])) % 45

    def set_mag_direction(self, mag_real, ang_real, mag_imag, ang_imag):
        """
        computes the tipper from the magnitude and direction of the real
        and imaginary components.

        Updates tipper

        No error propagation yet
        """

        self.tipper[:, 0, 0].real = np.sqrt((mag_real ** 2 * np.arctan(ang_real) ** 2) / \
                                            (1 - np.arctan(ang_real) ** 2))

        self.tipper[:, 0, 1].real = np.sqrt(mag_real ** 2 / \
                                            (1 - np.arctan(ang_real) ** 2))

        self.tipper[:, 0, 0].imag = np.sqrt((mag_imag ** 2 * np.arctan(ang_imag) ** 2) / \
                                            (1 - np.arctan(ang_imag) ** 2))

        self.tipper[:, 0, 1].imag = np.sqrt(mag_imag ** 2 / \
                                            (1 - np.arctan(ang_imag) ** 2))
        # for consistency recalculate mag and angle
        self.compute_mag_direction()
        self.compute_amp_phase()

    @property
    def mag_real(self):
        return self._mag_real

    @property
    def mag_imag(self):
        return self._mag_imag

    @property
    def angle_real(self):
        return self._angle_real

    @property
    def angle_imag(self):
        return self._angle_imag

    @property
    def mag_err(self):
        return self._mag_err

    @property
    def angle_err(self):
        return self._angle_err

    # ----rotate---------------------------------------------------------------
    def rotate(self, alpha):
        """
        Rotate  Tipper array.

        Rotation angle must be given in degrees. All angles are referenced
        to geographic North=0, positive in clockwise direction.
        (Mathematically negative!)

        In non-rotated state, 'X' refs to North and 'Y' to East direction.

        Updates the attributes:
    			* *tipper*
    			* *tipper_err*
    			* *rotation_angle*

        """

        if self.tipper is None:
            self._logger.error('tipper array is "None" - I cannot rotate that')
            return

        # check for iterable list/set of angles - if so, it must have length 1
        # or same as len(tipper):
        if np.iterable(alpha) == 0:
            try:
                degreeangle = float(alpha % 360)
            except ValueError:
                self._logger.error('"Angle" must be a valid number (in degrees)')
                return

            # make an n long list of identical angles
            lo_angles = [degreeangle for ii in self.tipper]
        elif len(alpha) == 1:
            try:
                degreeangle = float(alpha % 360)
            except ValueError:
                self._logger.error('"Angle" must be a valid number (in degrees)')
                return
            # make an n long list of identical angles
            lo_angles = [degreeangle for ii in self.tipper]
        else:
            try:
                lo_angles = [float(ii % 360) for ii in alpha]
            except ValueError:
                self._logger.error('"Angles" must be valid numbers (in degrees)')
                return

        self.rotation_angle = np.array([(oldangle + lo_angles[ii]) % 360
                                        for ii, oldangle in enumerate(self.rotation_angle)])

        if len(lo_angles) != len(self.tipper):
            self._logger.error('Wrong number Number of "angles" - need %ii ' % (len(self.tipper)))
            self.rotation_angle = 0.
            return

        tipper_rot = copy.copy(self.tipper)
        tipper_err_rot = copy.copy(self.tipper_err)

        for idx_freq in range(len(tipper_rot)):
            angle = lo_angles[idx_freq]

            if self.tipper_err is not None:
                tipper_rot[idx_freq], tipper_err_rot[idx_freq] = \
                    MTcc.rotatevector_incl_errors(self.tipper[idx_freq, :, :],
                                                  angle,
                                                  self.tipper_err[idx_freq, :, :])
            else:
                tipper_rot[idx_freq], tipper_err_rot = \
                    MTcc.rotatevector_incl_errors(self.tipper[idx_freq, :, :],
                                                  angle)

        self.tipper = tipper_rot
        self.tipper_err = tipper_err_rot

        # for consistency recalculate mag and angle
        self.compute_mag_direction()

        # for consistency recalculate amplitude and phase
        self.compute_amp_phase()


# ------------------------
def correct4sensor_orientation(Z_prime, Bx=0, By=90, Ex=0, Ey=90,
                               Z_prime_error=None):
    """
    Correct a Z-array for wrong orientation of the sensors.

    Assume, E' is measured by sensors orientated with the angles
        E'x: a
        E'y: b

    Assume, B' is measured by sensors orientated with the angles
        B'x: c
        B'y: d

    With those data, one obtained the impedance tensor Z':
        E' = Z' * B'

    Now we define change-of-basis matrices T,U so that
        E = T * E'
        B = U * B'

    =>   T contains the expression of the E'-basis in terms of E
    (the standard basis)
    and  U contains the expression of the B'-basis in terms of B
    (the standard basis)
    The respective expressions for E'x-basis vector and E'y-basis
    vector are the columns of T.
    The respective expressions for B'x-basis vector and B'y-basis
    vector are the columns of U.

    We obtain the impedance tensor in default coordinates as:

    E' = Z' * B' => T^(-1) * E = Z' * U^(-1) * B
                 => E = T * Z' * U^(-1) * B
                 => Z = T * Z' * U^(-1)

    :param Z_prime: impedance tensor to be adjusted
    :dtype Z_prime: np.ndarray(num_freq, 2, 2, dtype='complex')


    :param Bx: orientation of Bx relative to geographic north (0)
				   *default* is 0
    :type Bx: float (angle in degrees)

    :param By:
    :type By: float (angle in degrees)
		         orientation of By relative to geographic north (0)
				 *default* is 90

    :param Ex: orientation of Ex relative to geographic north (0)
				   *default* is 0
    :type Ex: float (angle in degrees)

    :param Ey: orientation of Ey relative to geographic north (0)
				  *default* is 90
    :type Ey: float (angle in degrees)

    :param Z_prime_error: impedance tensor error (std)
            					 *default* is None
    :type Z_prime_error: np.ndarray(Z_prime.shape)

    :returns: adjusted impedance tensor
    :rtype: np.ndarray(Z_prime.shape, dtype='complex')

    :returns: impedance tensor standard deviation in
					default orientation
    :rtype: np.ndarray(Z_prime.shape, dtype='real')
    """
    try:
        if len(Z_prime.shape) != 2:
            raise
        if Z_prime.shape != (2, 2):
            raise

        if Z_prime.dtype not in ['complex', 'float', 'int']:
            raise

        Z_prime = np.matrix(Z_prime)

    except:
        raise MTex.MTpyError_inputarguments('ERROR - Z array not valid!' + \
                                            'Must be 2x2 complex array')

    if Z_prime_error is not None:
        try:
            if len(Z_prime_error.shape) != 2:
                raise
            if Z_prime_error.shape != (2, 2):
                raise

            if Z_prime_error.dtype not in ['float', 'int']:
                raise

        except:
            raise MTex.MTpyError_inputarguments('ERROR - Z-error array not' + \
                                                'valid! Must be 2x2 real array')

    T = np.matrix(np.zeros((2, 2)))
    U = np.matrix(np.zeros((2, 2)))

    dummy1 = cmath.rect(1, math.radians(Ex))

    T[0, 0] = np.real(dummy1)
    T[1, 0] = np.imag(dummy1)
    dummy2 = cmath.rect(1, math.radians(Ey))
    T[0, 1] = np.real(dummy2)
    T[1, 1] = np.imag(dummy2)

    dummy3 = cmath.rect(1, math.radians(Bx))
    U[0, 0] = np.real(dummy3)
    U[1, 0] = np.imag(dummy3)
    dummy4 = cmath.rect(1, math.radians(By))
    U[0, 1] = np.real(dummy4)
    U[1, 1] = np.imag(dummy4)

    try:
        z_arr = np.array(np.dot(T, np.dot(Z_prime, U.I)))
    except:
        raise MTex.MTpyError_inputarguments("ERROR - Given angles do not" + \
                                            "define basis for 2 dimensions - cannot convert Z'")

    z_err_arr = copy.copy(Z_prime_error)

    # TODO: calculate error propagation

    return z_arr, z_err_arr
