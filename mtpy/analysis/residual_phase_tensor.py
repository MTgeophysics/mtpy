# -*- coding: utf-8 -*-
"""
Created on Wed Feb 22 14:13:57 2023

@author: jpeacock
"""

# =============================================================================
# Imports
# =============================================================================
import numpy as np
from loguru import logger

from mtpy.core.transfer_function.pt import PhaseTensor
import mtpy.utils.calculator as MTcc

# =============================================================================
class ResidualPhaseTensor:
    """
    PhaseTensor class - generates a Phase Tensor (PT) object DeltaPhi
    DeltaPhi = 1 - Phi1^-1*Phi2
    """

    def __init__(
        self, pt_object1=None, pt_object2=None, residual_type="heise"
    ):
        """
        Initialise an instance of the ResidualPhaseTensor class.
        Optional input:
        pt_object1 : instance of the PhaseTensor class
        pt_object2 : instance of the PhaseTensor class
        Initialise the attributes with None
        """
        self.logger = logger

        self.residual_pt = None
        self.rpt = None
        self.rpt_error = None
        self.pt1 = None
        self.pt2 = None
        self.pt1_error = None
        self.pt2_error = None
        self.frequency = None
        self.residual_type = residual_type

        if pt_object1 is not None or pt_object2 is not None:
            if not (
                (
                    isinstance(pt_object1, PhaseTensor)
                    and isinstance(pt_object2, PhaseTensor)
                )
            ):
                self.logger.warning(type(pt_object1), type(pt_object2))
                raise TypeError(
                    "ERROR - arguments must be instances "
                    "of the PhaseTensor class"
                )

            self.pt1 = pt_object1
            self.pt2 = pt_object2

            self.frequency = self.pt1.frequency
            self.compute_residual_pt()

    def compute_residual_pt(self):
        """
        Read in two instance of the MTpy PhaseTensor class.
        Update attributes:
        rpt, rpt_error, _self.pt1, _self.pt2, _self.pt1_error, _self.pt2_error
        """

        # --> compute residual phase tensor
        if self.pt1 is not None and self.pt2 is not None:
            if self.pt1.pt.dtype not in [float, int]:
                raise ValueError
            if self.pt2.pt.dtype not in [float, int]:
                raise ValueError
            if not self.pt1.pt.shape == self.pt2.pt.shape:
                raise TypeError("PT arrays not the same shape")
            if not len(self.pt1.pt.shape) in [2, 3]:
                raise TypeError("PT array is not a valid shape")
            if self.residual_type == "heise":
                if len(self.pt1.pt.shape) == 3:
                    self.rpt = np.zeros_like(self.pt1.pt)

                    for idx in range(len(self.pt1.pt)):
                        with np.errstate(divide="ignore", invalid="ignore"):
                            try:
                                self.rpt[idx] = np.eye(2) - 0.5 * (
                                    np.dot(
                                        np.matrix(self.pt1.pt[idx]).I,
                                        np.matrix(self.pt2.pt[idx]),
                                    )
                                    + np.dot(
                                        np.matrix(self.pt2.pt[idx]),
                                        np.matrix(self.pt1.pt[idx]).I,
                                    )
                                )
                            except np.linalg.LinAlgError:
                                self.rpt[idx] = np.zeros((2, 2))

                else:
                    self.rpt = np.zeros((1, 2, 2))
                    try:
                        with np.errstate(divide="ignore", invalid="ignore"):
                            self.rpt[idx] = np.eye(2) - 0.5 * (
                                np.dot(
                                    np.matrix(self.pt2.pt[idx]).I,
                                    np.matrix(self.pt1.pt[idx]),
                                )
                                + np.dot(
                                    np.matrix(self.pt1.pt[idx]),
                                    np.matrix(self.pt2.pt[idx]).I,
                                )
                            )

                    except np.linalg.LinAlgError:
                        pass

            elif self.residual_type == "booker":
                self.rpt = self.pt1.pt - self.pt2.pt

        else:
            self.logger.warning(
                "Could not determine ResPT - both PhaseTensor objects must"
                "contain PT arrays of the same shape"
            )

        # --> compute residual error

        if self.pt1.pt_error is not None and self.pt2.pt_error is not None:
            self.rpt_error = np.zeros(self.rpt.shape)
            try:
                if (self.pt1.pt_error.dtype not in [float, int]) or (
                    self.pt2.pt_error.dtype not in [float, int]
                ):
                    raise ValueError
                if not self.pt1.pt_error.shape == self.pt2.pt_error.shape:
                    raise ValueError
                if not len(self.pt1.pt_error.shape) in [2, 3]:
                    raise ValueError
                if self.rpt_error is not None:
                    if self.rpt_error.shape != self.pt1.pt_error.shape:
                        raise ValueError
                if self.residual_type == "heise":
                    if len(self.pt1.pt_error.shape) == 3:
                        self.rpt_error = np.zeros((len(self.pt1.pt), 2, 2))

                        with np.errstate(divide="ignore", invalid="ignore"):
                            for idx in range(len(self.pt1.pt_error)):
                                matrix1 = self.pt1.pt[idx]
                                matrix1error = self.pt1.pt_error[idx]
                                try:
                                    (
                                        matrix2,
                                        matrix2error,
                                    ) = MTcc.invertmatrix_incl_errors(
                                        self.pt2.pt[idx],
                                        inmatrix_error=self.pt2.pt_error[idx],
                                    )

                                    (
                                        summand1,
                                        error1,
                                    ) = MTcc.multiplymatrices_incl_errors(
                                        matrix2,
                                        matrix1,
                                        inmatrix1_error=matrix2error,
                                        inmatrix2_error=matrix1error,
                                    )
                                    (
                                        summand2,
                                        error2,
                                    ) = MTcc.multiplymatrices_incl_errors(
                                        matrix1,
                                        matrix2,
                                        inmatrix1_error=matrix1error,
                                        inmatrix2_error=matrix2error,
                                    )
                                    self.rpt_error[idx] = np.sqrt(
                                        0.25 * error1**2 + 0.25 * error2**2
                                    )
                                except ValueError:
                                    self.rpt_error[idx] = 1e10

                    else:
                        self.rpt_error = np.zeros((1, 2, 2))
                        try:
                            with np.errstate(
                                divide="ignore", invalid="ignore"
                            ):
                                self.rpt_error[0] = np.eye(2) - 0.5 * np.array(
                                    np.dot(
                                        np.matrix(self.pt2.pt).I,
                                        np.matrix(self.pt1.pt),
                                    )
                                    + np.dot(
                                        np.matrix(self.pt1.pt),
                                        np.matrix(self.pt2.pt).I,
                                    )
                                )
                                matrix1 = self.pt1.pt
                                matrix1error = self.pt1.pt_error
                                (
                                    matrix2,
                                    matrix2error,
                                ) = MTcc.invertmatrix_incl_errors(
                                    self.pt2.pt,
                                    inmatrix_error=self.pt2.pt_error,
                                )

                                (
                                    summand1,
                                    error1,
                                ) = MTcc.multiplymatrices_incl_errors(
                                    matrix2,
                                    matrix1,
                                    inmatrix1_error=matrix2error,
                                    inmatrix2_error=matrix1error,
                                )
                                (
                                    summand2,
                                    error2,
                                ) = MTcc.multiplymatrices_incl_errors(
                                    matrix1,
                                    matrix2,
                                    inmatrix1_error=matrix1error,
                                    inmatrix2_error=matrix2error,
                                )

                                self.rpt_error = np.sqrt(
                                    0.25 * error1**2 + 0.25 * error2**2
                                )
                        except ValueError:
                            self.rpt_error[idx] = 1e10

                elif self.residual_type == "booker":
                    self.rpt_error = self.pt1.pt_error + self.pt2.pt_error

            except ValueError:
                raise TypeError(
                    "ERROR - both PhaseTensor objects must"
                    "contain PT-erroror arrays of the same shape"
                )

        else:
            self.logger.warning(
                "Could not determine Residual PT uncertainties - both"
                " PhaseTensor objects must contain PT-error arrays of the"
                "same shape"
            )

        # --> make a pt object that is the residual phase tensor
        with np.errstate(divide="ignore", invalid="ignore"):
            self.residual_pt = PhaseTensor(
                pt=self.rpt,
                pt_error=self.rpt_error,
                frequency=self.frequency,
            )

    def read_pts(self, pt1, pt2, pt1_error=None, pt2error=None):
        """
        Read two PT arrays and calculate the ResPT array (incl. uncertainties).
        Input:
        - 2x PT array
        Optional:
        - 2x pt_erroror array

        """

        try:
            if self.pt1.shape != self.pt2.shape:
                raise

        except:
            raise TypeError(
                "ERROR - could not build ResPT array from given PT arrays - check shapes! "
            )
        # TODO - check arrays here:

        pt_o1 = PhaseTensor(pt_array=self.pt1, pt_error_array=self.pt1_error)
        pt_o2 = PhaseTensor(pt_array=self.pt2, pt_error_array=self.pt2_error)

        self.compute_residual_pt(pt_o1, pt_o2)

    def set_rpt(self, rpt_array):
        """
        Set the attribute 'rpt' (ResidualPhaseTensor array).
        Input:
        ResPT array
        Test for shape, but no test for consistency!
        """
        if (self.rpt is not None) and (self.rpt.shape != rpt_array.shape):
            self.logger.erroror(
                'Shape of "ResPT" array does not match shape of existing rpt array: %s ; %s'
                % (str(rpt_array.shape), str(self.rpt.shape))
            )
            return

        self.rpt = rpt_array

        # --> make a pt object that is the residual phase tensor
        self.residual_pt = PhaseTensor(
            pt_array=self.rpt,
            pt_error_array=self.rpt_error,
            freq=self.frequency,
        )

    def set_rpt_error(self, rpt_error_array):
        """
        Set the attribute 'rpt_error' (ResidualPhaseTensor-erroror array).
        Input:
        ResPT-erroror array
        Test for shape, but no test for consistency!
        """
        if (self.rpt_error is not None) and (
            self.rpt_error.shape != rpt_error_array.shape
        ):
            self.logger.erroror(
                'Error - shape of "ResPT-erroror" array does not match shape of existing rpt_error array: %s ; %s'
                % (str(rpt_error_array.shape), str(self.rpt_error.shape))
            )
            return

        self.rpt_error = rpt_error_array

        # --> make a pt object that is the residual phase tensor
        self.residual_pt = PhaseTensor(
            pt_array=self.rpt,
            pt_error_array=self.rpt_error,
            freq=self.frequency,
        )
