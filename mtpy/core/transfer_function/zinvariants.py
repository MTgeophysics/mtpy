# -*- coding: utf-8 -*-
"""
Created on Wed May 08 09:40:42 2013

Originally written by Stephan Thiel in Matlab 2005
translated to Python by Lars Krieger

Revised by J. Peacock 2023 to fit with version 2.
"""

# =============================================================================
# Imports
# =============================================================================
import copy
import numpy as np

# =============================================================================


class ZInvariants:
    """
    Calculates invariants from Weaver et al. [2000, 2003].  At the moment it
    does not calculate the error for each invariant, only the strike.

    :type z: complex np.array(nf,2,2)
    :param z: impedance tensor array



    Further reading
    ----------------

        Weaver, J. T., Agarwal, A. K., Lilley, F. E. M., 2000,
           Characterization of the magnetotelluric tensor in terms of its
           invariants, Geophysical Journal International, 141, 321--336.

        Weaver, J. T., Agarwal, A. K., Lilley, F. E. M., 2003,
            The relationship between the magnetotelluric tensor invariants and
            the phase tensor of Caldwell, Bibby and Brown,
            presented at 3D Electromagnetics III, ASEG, paper 43.

        Lilley, F. E. M, 1998, Magnetotelluric tensor dcomposition: 1: Theory
            for a basic procedure, Geophysics, 63, 1885--1897.

        Lilley, F. E. M, 1998, Magnetotelluric tensor dcomposition: 2: Examples
            of a basic procedure, Geophysics, 63, 1898--1907.

        Szarka, L. and Menvielle, M., 1997, Analysis of rotational invariants
            of the magnetotelluric impedance tensor, Geophysical Journal
            International, 129, 133--142.

    """

    def __init__(
        self,
        z=None,
    ):

        self.z = z

    def __str__(self):
        return f"Weaver Invariants \n\tHas Impedance:  {self.has_impedance()}"

    def __repr__(self):
        return self.__str__()

    def has_impedance(self):
        if np.all(self.z == 0):
            return False
        return True

    def _zero_to_nan(self, array):
        """
         convert zeros to nans
        :param array: DESCRIPTION
        :type array: TYPE
        :return: DESCRIPTION
        :rtype: TYPE

        """
        array[np.where(array == 0)] = np.nan
        return array

    @property
    def _x1(self):
        if self.has_impedance():
            return self._zero_to_nan(
                0.5 * (self.z[:, 0, 0].real + self.z[:, 1, 1].real)
            )

    @property
    def _x2(self):
        if self.has_impedance():
            return self._zero_to_nan(
                0.5 * (self.z[:, 0, 1].real + self.z[:, 1, 0].real)
            )

    @property
    def _x3(self):
        if self.has_impedance():
            return self._zero_to_nan(
                0.5 * (self.z[:, 0, 0].real - self.z[:, 1, 1].real)
            )

    @property
    def _x4(self):
        if self.has_impedance():
            return self._zero_to_nan(
                0.5 * (self.z[:, 0, 1].real - self.z[:, 1, 0].real)
            )

    @property
    def _e1(self):
        if self.has_impedance():
            return self._zero_to_nan(
                0.5 * (self.z[:, 0, 0].imag + self.z[:, 1, 1].imag)
            )

    @property
    def _e2(self):
        if self.has_impedance():
            return self._zero_to_nan(
                0.5 * (self.z[:, 0, 1].imag + self.z[:, 1, 0].imag)
            )

    @property
    def _e3(self):
        if self.has_impedance():
            return self._zero_to_nan(
                0.5 * (self.z[:, 0, 0].imag - self.z[:, 1, 1].imag)
            )

    @property
    def _e4(self):
        if self.has_impedance():
            return self._zero_to_nan(
                0.5 * (self.z[:, 0, 1].imag - self.z[:, 1, 0].imag)
            )

    @property
    def _ex(self):
        if self.has_impedance():
            ex = (
                self._x1 * self._e1
                - self._x2 * self._e2
                - self._x3 * self._e3
                + self._x4 * self._e4
            )
            ex[np.where(ex == 0)] = np.nan

            return self._zero_to_nan(ex)

    @property
    def normalizing_real(self):
        """inv 1"""
        if self.has_impedance():
            return np.sqrt(self._x4**2 + self._x1**2)

    @property
    def normalizing_imag(self):
        """inv 2"""
        if self.has_impedance():
            return np.sqrt(self._e4**2 + self._e1**2)

    @property
    def anisotropic_real(self):
        """inv 3"""
        if self.has_impedance():
            return (
                np.sqrt(self._x2**2 + self._x3**2) / self.normalizing_real
            )

    @property
    def anisotropic_imag(self):
        """inv 4"""
        if self.has_impedance():
            return (
                np.sqrt(self._e2**2 + self._e3**2) / self.normalizing_imag
            )

    @property
    def electric_twist(self):
        """inv 5"""
        if self.has_impedance():
            return (self._x4 * self._e1 + self._x1 * self._e4) / (
                self.normalizing_real * self.normalizing_imag
            )

    @property
    def phase_distortion(self):
        """inv 6"""
        if self.has_impedance():
            return (self._x4 * self._e1 - self._x1 * self._e4) / (
                self.normalizing_real * self.normalizing_imag
            )

    @property
    def dimensionality(self):
        """q"""
        if self.has_impedance():
            return np.sqrt(
                (
                    (self._x1 * self._e2 - self._x2 * self._e1) / self._ex
                    - (self._x3 * self._e4 - self._x4 * self._e3) / self._ex
                )
                ** 2
                + (
                    (self._x1 * self._e3 - self._x3 * self._e1) / self._ex
                    + (self._x2 * self._e4 - self._x4 * self._e2) / self._ex
                )
                ** 2
            )

    @property
    def structure_3d(self):
        """inv 7"""
        if self.has_impedance():
            return (
                (self._x4 * self._e1 - self._x1 * self._e4) / self._ex
                - (self._x2 * self._e3 - self._x3 * self._e2) / self._ex
            ) / self.dimensionality

    @property
    def strike(self):
        if self.has_impedance():
            return (
                0.5
                * np.rad2deg(
                    np.arctan2(
                        (self._x1 * self._e2 - self._x2 * self._e1) / self._ex
                        - (self._x3 * self._e4 - self._x4 * self._e3)
                        / self._ex,
                        (self._x1 * self._e3 - self._x3 * self._e1) / self._ex
                        + (self._x2 * self._e4 - self._x4 * self._e2)
                        / self._ex,
                    )
                )
                % 360
            )

    @property
    def strike_error(self):
        if self.has_impedance():
            return np.rad2deg(abs(0.5 * np.arcsin(self.structure_3d)))
