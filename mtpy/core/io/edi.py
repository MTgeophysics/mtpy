# -*- coding: utf-8 -*-
"""
.. module:: EDI
   :synopsis: Deal with EDI files.  The Edi class can read and write an .edi
             file, the 'standard format' of magnetotellurics.  Each section
             of the .edi file is given its own class, so the elements of each
             section are attributes for easy access.

.. moduleauthor:: Jared Peacock <jpeacock@usgs.gov>
"""

# ==============================================================================
#  Imports
# ==============================================================================
import os
import numpy as np
import logging
from pathlib import Path

import mtpy.utils.gis_tools as gis_tools
import mtpy.utils.exceptions as MTex
import mtpy.utils.filehandling as MTfh
from mtpy.utils.mttime import MTime, get_now_utc
import mtpy.core.z as MTz
from mtpy.core import metadata
from mtpy.core import mt
from mtpy import __version__

import scipy.stats.distributions as ssd

tab = " " * 4
# ==============================================================================
# EDI Class
# ==============================================================================
class Edi(object):
    """
    This class is for .edi files, mainly reading and writing.  Has been tested
    on Winglink and Phoenix output .edi's, which are meant to follow the
    archaic EDI format put forward by SEG. Can read impedance, Tipper and/or
    spectra data.

    The Edi class contains a class for each major section of the .edi file.

    Frequency and components are ordered from highest to lowest frequency.

    :param fn: full path to .edi file to be read in.
                  *default* is None. If an .edi file is input, it is
                  automatically read in and attributes of Edi are filled
    :type fn: string

    ===================== =====================================================
    Methods               Description
    ===================== =====================================================
    read_edi_file         Reads in an edi file and populates the associated
                          classes and attributes.
    write_edi_file        Writes an .edi file following the EDI format given
                          the apporpriate attributes are filled.  Writes out
                          in impedance and Tipper format.
    _read_data            Reads in the impedance and Tipper blocks, if the
                          .edi file is in 'spectra' format, read_data converts
                          the data to impedance and Tipper.
    _read_mt              Reads impedance and tipper data from the appropriate
                          blocks of the .edi file.
    _read_spectra         Reads in spectra data and converts it to impedance
                          and Tipper data.
    ===================== =====================================================

    ===================== ========================================== ==========
    Attributes            Description                                default
    ===================== ========================================== ==========
    Data             DataSection class, contains basic
                          information on the data collected and in
                          whether the data is in impedance or
                          spectra.
    Measurement    DefineMeasurement class, contains
                          information on how the data was
                          collected.
    fn                full path to edi file read in              None
    Header                Header class, contains metadata on
                          where, when, and who collected the data
    Info                  Information class, contains information
                          on how the data was processed and how the
                          transfer functions where estimated.
    Tipper                mtpy.core.z.Tipper class, contains the
                          tipper data
    Z                     mtpy.core.z.Z class, contains the
                          impedance data
    _block_len            number of data in one line.                6
    _data_header_str      header string for each of the data         '>!****{0}****!'
                          section
    _num_format           string format of data.                     ' 15.6e'
    _t_labels             labels for tipper blocks
    _z_labels             labels for impedance blocks
    ===================== ========================================== ==========

    :Change Latitude: ::

        >>> import mtpy.core.edi as mtedi
        >>> edi_obj = mtedi.Edi(fn=r"/home/mt/mt01.edi")
        >>> # change the latitude
        >>> edi_obj.header.lat = 45.7869
        >>> new_edi_fn = edi_obj.write_edi_file()
    """

    def __init__(self, fn=None):
        self.logger = logging.getLogger(f"{__name__}.{self.__class__.__name__}")
        self._fn = None
        self._edi_lines = None

        self.fn = fn

        self.Header = Header()
        self.Info = Information()
        self.Measurement = DefineMeasurement()
        self.Data = DataSection()
        self.Z = MTz.Z()
        self.Tipper = MTz.Tipper()

        self._z_labels = [
            ["zxxr", "zxxi", "zxx.var"],
            ["zxyr", "zxyi", "zxy.var"],
            ["zyxr", "zyxi", "zyx.var"],
            ["zyyr", "zyyi", "zyy.var"],
        ]

        self._t_labels = [
            ["txr.exp", "txi.exp", "txvar.exp"],
            ["tyr.exp", "tyi.exp", "tyvar.exp"],
        ]

        self._data_header_str = ">!****{0}****!\n"

        self._num_format = " 15.6e"
        self._block_len = 6

        if self.fn is not None:
            self.read_edi_file()

    @property
    def fn(self):
        return self._fn

    @fn.setter
    def fn(self, fn):
        if fn is not None:
            self._fn = Path(fn)

    def read_edi_file(self, fn=None):
        """
        Read in an edi file and fill attributes of each section's classes.
        Including:
            * Header
            * Info
            * Measurement
            * Data
            * Z
            * Tipper

            .. note:: Automatically detects if data is in spectra format.  All
                  data read in is converted to impedance and Tipper.


        :param fn: full path to .edi file to be read in
                       *default* is None
        :type fn: string

        :Example: ::

            >>> import mtpy.core.Edi as mtedi
            >>> edi_obj = mtedi.Edi()
            >>> edi_obj.read_edi_file(fn=r"/home/mt/mt01.edi")

        """

        if fn is not None:
            self.fn = fn

        if self.fn is None:
            msg = "Must input EDI file to read"
            self.logger.error(msg)
            raise MTex.MTpyError_EDI(msg)

        if not self.fn.exists():
            msg = f"Cannot find EDI file: {self.fn}"
            self.logger.error(msg)
            raise MTex.MTpyError_EDI(msg)

        with open(self.fn, "r") as fid:
            self._edi_lines = _validate_edi_lines(fid.readlines())

        self.Header = Header(edi_lines=self._edi_lines)
        self.Info = Information(edi_lines=self._edi_lines)
        self.Measurement = DefineMeasurement(edi_lines=self._edi_lines)
        self.Data = DataSection(edi_lines=self._edi_lines)

        self._read_data()

        if self.Header.lat is None:
            self.Header.lat = self.Measurement.reflat
            self.logger.info(
                "Got latitude from reflat for {0}".format(self.Header.dataid)
            )
        if self.Header.lon is None:
            self.Header.lon = self.Measurement.reflon
            self.logger.info(
                "Got longitude from reflon for {0}".format(self.Header.dataid)
            )
        if self.Header.elev is None:
            self.Header.elev = self.Measurement.refelev
            self.logger.info(
                "Got elevation from refelev for {0}".format(self.Header.dataid)
            )

        self.logger.info("Read in edi file for station {0}".format(self.Header.dataid))

    def _read_data(self):
        """
        Read either impedance or spectra data depending on what the type is
        in the data section.
        """

        if self.fn is None:
            raise MTex.MTpyError_EDI("No edi file input, check fn")
        if self.fn.exists() is False:
            raise MTex.MTpyError_EDI("No edi file input, check fn")

        lines = self._edi_lines[self.Data.line_num :]

        if self.Data.data_type == "spectra":
            self.logger.info("Converting Spectra to Impedance and Tipper")
            self.logger.info(
                "Check to make sure input channel list is correct if the data looks incorrect"
            )
            if self.Data.nchan == 5:
                c_list = ["hx", "hy", "hz", "ex", "ey"]
            elif self.Data.nchan == 4:
                c_list = ["hx", "hy", "ex", "ey"]
            elif self.Data.nchan == 6:
                c_list = ["hx", "hy", "ex", "ey", "rhx", "rhy"]
            elif self.Data.nchan == 7:
                c_list = ["hx", "hy", "hz", "ex", "ey", "rhx", "rhy"]
            self._read_spectra(lines, comp_list=c_list)

        elif self.Data.data_type == "z":
            self._read_mt(lines)

    def _read_mt(self, data_lines):
        """
        Read in impedance and tipper data

        :param data_lines: list of data lines from the edi file
        :type data_lines: list
        """
        flip = False
        data_dict = {}
        data_find = False
        for line in data_lines:
            line = line.strip()
            if ">" in line and "!" not in line:
                line_list = line[1:].strip().split()
                if len(line_list) == 0:
                    continue
                key = line_list[0].lower()
                if key[0] == "z" or key[0] == "t" or key == "freq":
                    data_find = True
                    data_dict[key] = []
                else:
                    data_find = False

            elif data_find and ">" not in line and "!" not in line:
                d_lines = line.strip().split()
                for ii, dd in enumerate(d_lines):
                    # check for empty values and set them to 0, check for any
                    # other characters sometimes there are ****** for a null
                    # component
                    try:
                        d_lines[ii] = float(dd)
                        if d_lines[ii] == 1.0e32:
                            d_lines[ii] = 0.0
                    except ValueError:
                        d_lines[ii] = 0.0
                data_dict[key] += d_lines

        # fill useful arrays
        freq_arr = np.array(data_dict["freq"], dtype=np.float)
        z_arr = np.zeros((freq_arr.size, 2, 2), dtype=np.complex)
        z_err_arr = np.zeros((freq_arr.size, 2, 2), dtype=np.float)

        # fill impedance tensor
        if "zxxr" in data_dict.keys():
            z_arr[:, 0, 0] = (
                np.array(data_dict["zxxr"]) + np.array(data_dict["zxxi"]) * 1j
            )
            z_err_arr[:, 0, 0] = np.array(data_dict["zxx.var"]) ** 0.5
        if "zxyr" in data_dict.keys():
            z_arr[:, 0, 1] = (
                np.array(data_dict["zxyr"]) + np.array(data_dict["zxyi"]) * 1j
            )
            z_err_arr[:, 0, 1] = np.array(data_dict["zxy.var"]) ** 0.5
        if "zyxr" in data_dict.keys():
            z_arr[:, 1, 0] = (
                np.array(data_dict["zyxr"]) + np.array(data_dict["zyxi"]) * 1j
            )
            z_err_arr[:, 1, 0] = np.array(data_dict["zyx.var"]) ** 0.5
        if "zyyr" in data_dict.keys():
            z_arr[:, 1, 1] = (
                np.array(data_dict["zyyr"]) + np.array(data_dict["zyyi"]) * 1j
            )
            z_err_arr[:, 1, 1] = np.array(data_dict["zyy.var"]) ** 0.5

        # check for order of frequency, we want high togit  low
        if freq_arr[0] < freq_arr[1]:
            self.logger.info("Ordered arrays to be arranged from high to low frequency")
            freq_arr = freq_arr[::-1]
            z_arr = z_arr[::-1]
            z_err_arr = z_err_arr[::-1]
            flip = True

        # set the attributes as private variables to avoid redundant estimation
        # of res and phase
        self.Z._freq = freq_arr
        self.Z._z = z_arr
        self.Z._z_err = z_err_arr

        try:
            self.Z.rotation_angle = np.array(data_dict["zrot"])
        except KeyError:
            self.Z.rotation_angle = np.zeros_like(freq_arr)

        # compute resistivity and phase
        self.Z.compute_resistivity_phase()

        # fill tipper data if there it exists
        tipper_arr = np.zeros((freq_arr.size, 1, 2), dtype=np.complex)
        tipper_err_arr = np.zeros((freq_arr.size, 1, 2), dtype=np.float)

        try:
            self.Tipper.rotation_angle = np.array(data_dict["trot"])
        except KeyError:
            try:
                self.Tipper.rotation_angle = np.array(data_dict["zrot"])
            except KeyError:
                self.Tipper.rotation_angle = np.zeros_like(freq_arr)

        if "txr.exp" in list(data_dict.keys()):
            tipper_arr[:, 0, 0] = (
                np.array(data_dict["txr.exp"]) + np.array(data_dict["txi.exp"]) * 1j
            )
            tipper_arr[:, 0, 1] = (
                np.array(data_dict["tyr.exp"]) + np.array(data_dict["tyi.exp"]) * 1j
            )

            tipper_err_arr[:, 0, 0] = np.array(data_dict["txvar.exp"]) ** 0.5
            tipper_err_arr[:, 0, 1] = np.array(data_dict["tyvar.exp"]) ** 0.5

            if flip:
                tipper_arr = tipper_arr[::-1]
                tipper_err_arr = tipper_err_arr[::-1]

        else:
            self.logger.info("Could not find any Tipper data.")

        self.Tipper._freq = freq_arr
        self.Tipper._tipper = tipper_arr
        self.Tipper._tipper_err = tipper_err_arr
        self.Tipper.compute_amp_phase()
        self.Tipper.compute_mag_direction()

    def _read_spectra(
        self, data_lines, comp_list=["hx", "hy", "hz", "ex", "ey", "rhx", "rhy"]
    ):
        """
        Read in spectra data and convert to impedance and Tipper.

        :param data_lines: list of lines from edi file
        :type data_lines: list

        :param comp_list: list of components that correspond to the columns
                          of the spectra data.
        :type comp_list: list
        """

        data_dict = {}
        avgt_dict = {}
        data_find = False
        for line in data_lines:
            if line.lower().find(">spectra") == 0 and line.find("!") == -1:
                line_list = _validate_str_with_equals(line)
                data_find = True

                # frequency will be the key
                try:
                    key = float(
                        [
                            ss.split("=")[1]
                            for ss in line_list
                            if ss.lower().find("freq") == 0
                        ][0]
                    )
                    data_dict[key] = []
                    avgt = float(
                        [
                            ss.split("=")[1]
                            for ss in line_list
                            if ss.lower().find("avgt") == 0
                        ][0]
                    )
                    avgt_dict[key] = avgt
                except ValueError:
                    self.logger.info("did not find frequency key")

            elif data_find and line.find(">") == -1 and line.find("!") == -1:
                data_dict[key] += [float(ll) for ll in line.strip().split()]

            elif line.find(">spectra") == -1:
                data_find = False

        # get an object that contains the indices for each component
        cc = index_locator(comp_list)

        freq_arr = np.array(sorted(list(data_dict.keys()), reverse=True))

        z_arr = np.zeros((len(list(data_dict.keys())), 2, 2), dtype=np.complex)
        t_arr = np.zeros((len(list(data_dict.keys())), 1, 2), dtype=np.complex)

        z_err_arr = np.zeros_like(z_arr, dtype=np.float)
        t_err_arr = np.zeros_like(t_arr, dtype=np.float)

        for kk, key in enumerate(freq_arr):
            spectra_arr = np.reshape(
                np.array(data_dict[key]), (len(comp_list), len(comp_list))
            )

            # compute cross powers
            s_arr = np.zeros_like(spectra_arr, dtype=np.complex)
            for ii in range(s_arr.shape[0]):
                for jj in range(ii, s_arr.shape[0]):
                    if ii == jj:
                        s_arr[ii, jj] = spectra_arr[ii, jj]
                    else:
                        # minus sign for complex conjugation
                        # original spectra data are of form <A,B*>, but we need
                        # the order <B,A*>...
                        # this is achieved by complex conjugation of the
                        # original entries
                        s_arr[ii, jj] = np.complex(
                            spectra_arr[jj, ii], -spectra_arr[ii, jj]
                        )
                        # keep complex conjugated entries in the lower
                        # triangular matrix:
                        s_arr[jj, ii] = np.complex(
                            spectra_arr[jj, ii], spectra_arr[ii, jj]
                        )

            # use formulas from Bahr/Simpson to convert the Spectra into Z
            # the entries of S are sorted like
            # <X,X*>  <X,Y*>  <X,Z*>  <X,En*>  <X,Ee*>  <X,Rx*>  <X,Ry*>
            #         <Y,Y*>  <Y,Z*>  <Y,En*>  <Y,Ee*>  <Y,Rx*>  <Y,Ry*>
            # .....

            z_arr[kk, 0, 0] = (
                s_arr[cc.ex, cc.rhx] * s_arr[cc.hy, cc.rhy]
                - s_arr[cc.ex, cc.rhy] * s_arr[cc.hy, cc.rhx]
            )
            z_arr[kk, 0, 1] = (
                s_arr[cc.ex, cc.rhy] * s_arr[cc.hx, cc.rhx]
                - s_arr[cc.ex, cc.rhx] * s_arr[cc.hx, cc.rhy]
            )
            z_arr[kk, 1, 0] = (
                s_arr[cc.ey, cc.rhx] * s_arr[cc.hy, cc.rhy]
                - s_arr[cc.ey, cc.rhy] * s_arr[cc.hy, cc.rhx]
            )
            z_arr[kk, 1, 1] = (
                s_arr[cc.ey, cc.rhy] * s_arr[cc.hx, cc.rhx]
                - s_arr[cc.ey, cc.rhx] * s_arr[cc.hx, cc.rhy]
            )

            z_arr[kk] /= (
                s_arr[cc.hx, cc.rhx] * s_arr[cc.hy, cc.rhy]
                - s_arr[cc.hx, cc.rhy] * s_arr[cc.hy, cc.rhx]
            )

            # compute error only if scipy package exists
            # 68% Quantil of the Fisher distribution:
            z_det = np.real(
                s_arr[cc.hx, cc.hx] * s_arr[cc.hy, cc.hy]
                - np.abs(s_arr[cc.hx, cc.hy] ** 2)
            )

            sigma_quantil = ssd.f.ppf(0.68, 4, avgt_dict[key] - 4)

            ## 1) Ex
            a = (
                s_arr[cc.ex, cc.hx] * s_arr[cc.hy, cc.hy]
                - s_arr[cc.ex, cc.hy] * s_arr[cc.hy, cc.hx]
            )
            b = (
                s_arr[cc.ex, cc.hy] * s_arr[cc.hx, cc.hx]
                - s_arr[cc.ex, cc.hx] * s_arr[cc.hx, cc.hy]
            )
            a /= z_det
            b /= z_det

            psi_squared = np.real(
                1.0
                / s_arr[cc.ex, cc.ex].real
                * (a * s_arr[cc.hx, cc.ex] + b * s_arr[cc.hy, cc.ex])
            )
            epsilon_squared = 1.0 - psi_squared

            scaling = (
                sigma_quantil
                * 4
                / (avgt_dict[key] - 4.0)
                * epsilon_squared
                / z_det
                * s_arr[cc.ex, cc.ex].real
            )
            z_err_arr[kk, 0, 0] = np.sqrt(abs(scaling * s_arr[cc.hy, cc.hy].real))
            z_err_arr[kk, 0, 1] = np.sqrt(abs(scaling * s_arr[cc.hx, cc.hx].real))

            # 2) EY
            a = (
                s_arr[cc.ey, cc.hx] * s_arr[cc.hy, cc.hy]
                - s_arr[cc.ey, cc.hy] * s_arr[cc.hy, cc.hx]
            )
            b = (
                s_arr[cc.ey, cc.hy] * s_arr[cc.hx, cc.hx]
                - s_arr[cc.ey, cc.hx] * s_arr[cc.hx, cc.hy]
            )
            a /= z_det
            b /= z_det

            psi_squared = np.real(
                1.0
                / np.real(s_arr[cc.ey, cc.ey])
                * (a * s_arr[cc.hx, cc.ey] + b * s_arr[cc.hy, cc.ey])
            )
            epsilon_squared = 1.0 - psi_squared

            scaling = (
                sigma_quantil
                * 4
                / (avgt_dict[key] - 4.0)
                * epsilon_squared
                / z_det
                * s_arr[cc.ey, cc.ey].real
            )
            z_err_arr[kk, 1, 0] = np.sqrt(abs(scaling * s_arr[cc.hy, cc.hy].real))
            z_err_arr[kk, 1, 1] = np.sqrt(abs(scaling * s_arr[cc.hx, cc.hx].real))

            # if HZ information is present:
            if len(comp_list) > 5:
                t_arr[kk, 0, 0] = (
                    s_arr[cc.hz, cc.rhx] * s_arr[cc.hy, cc.rhy]
                    - s_arr[cc.hz, cc.rhy] * s_arr[cc.hy, cc.rhx]
                )
                t_arr[kk, 0, 1] = (
                    s_arr[cc.hz, cc.rhy] * s_arr[cc.hx, cc.rhx]
                    - s_arr[cc.hz, cc.rhx] * s_arr[cc.hx, cc.rhy]
                )

                t_arr[kk] /= (
                    s_arr[cc.hx, cc.rhx] * s_arr[cc.hy, cc.rhy]
                    - s_arr[cc.hx, cc.rhy] * s_arr[cc.hy, cc.rhx]
                )

                a = (
                    s_arr[cc.hz, cc.hx] * s_arr[cc.hy, cc.hy]
                    - s_arr[cc.hz, cc.hy] * s_arr[cc.hy, cc.hx]
                )
                b = (
                    s_arr[cc.hz, cc.hy] * s_arr[cc.hx, cc.hx]
                    - s_arr[cc.hz, cc.hx] * s_arr[cc.hx, cc.hy]
                )
                a /= z_det
                b /= z_det

                psi_squared = np.real(
                    1.0
                    / s_arr[cc.hz, cc.hz].real
                    * (a * s_arr[cc.hx, cc.hz] + b * s_arr[cc.hy, cc.hz])
                )
                epsilon_squared = 1.0 - psi_squared

                scaling = (
                    sigma_quantil
                    * 4
                    / (avgt_dict[key] - 4.0)
                    * epsilon_squared
                    / z_det
                    * s_arr[cc.hz, cc.hz].real
                )
                t_err_arr[kk, 0, 0] = np.sqrt(abs(scaling * s_arr[cc.hy, cc.hy].real))
                t_err_arr[kk, 0, 1] = np.sqrt(abs(scaling * s_arr[cc.hx, cc.hx].real))

        # check for nans
        z_err_arr = np.nan_to_num(z_err_arr)
        t_err_arr = np.nan_to_num(t_err_arr)

        z_err_arr[np.where(z_err_arr == 0.0)] = 1.0
        t_err_arr[np.where(t_err_arr == 0.0)] = 1.0

        # be sure to fill attributes
        self.Z.freq = freq_arr
        self.Z.z = z_arr
        self.Z.z_err = z_err_arr
        self.Z.rotation_angle = np.zeros_like(freq_arr)
        self.Z.compute_resistivity_phase()

        self.Tipper.tipper = t_arr
        self.Tipper.tipper_err = t_err_arr
        self.Tipper.freq = freq_arr
        self.Tipper.rotation_angle = np.zeros_like(freq_arr)
        self.Tipper.compute_amp_phase()
        self.Tipper.compute_mag_direction()

    def write_edi_file(
        self, new_edi_fn=None, longitude_format="LON", latlon_format="dms"
    ):
        """
        Write a new edi file from either an existing .edi file or from data
        input by the user into the attributes of Edi.


        :param new_edi_fn: full path to new edi file.
                           *default* is None, which will write to the same
                           file as the input .edi with as:
                           r"/home/mt/mt01_1.edi"
        :type new_edi_fn: string
        :param longitude_format:  whether to write longitude as LON or LONG. 
                                  options are 'LON' or 'LONG', default 'LON'
        :type longitude_format:  string
        :param latlon_format:  format of latitude and longitude in output edi,
                               degrees minutes seconds ('dms') or decimal 
                               degrees ('dd')
        :type latlon_format:  string

        :returns: full path to new edi file
        :rtype: string

        :Example: ::

            >>> import mtpy.core.edi as mtedi
            >>> edi_obj = mtedi.Edi(fn=r"/home/mt/mt01/edi")
            >>> edi_obj.Header.dataid = 'mt01_rr'
            >>> n_edi_fn = edi_obj.write_edi_file()
        """

        if new_edi_fn is None:
            if self.fn is not None:
                new_edi_fn = self.fn
            else:
                new_edi_fn = os.path.join(
                    os.getcwd(), "{0}.edi".format(self.Header.dataid)
                )
        new_edi_fn = MTfh.make_unique_filename(new_edi_fn)

        if self.Header.dataid is None:
            self.read_edi_file()

        # write lines
        header_lines = self.Header.write_header(
            longitude_format=longitude_format, latlon_format=latlon_format
        )
        info_lines = self.Info.write_info()
        define_lines = self.Measurement.write_measurement(
            longitude_format=longitude_format, latlon_format=latlon_format
        )
        dsect_lines = self.Data.write_Data(
            over_dict={"nfreq": len(self.Z.freq)}
        )

        # write out frequencies
        freq_lines = [self._data_header_str.format("frequencies".upper())]
        freq_lines += self._write_data_block(self.Z.freq, "freq")

        # write out rotation angles
        zrot_lines = [self._data_header_str.format("impedance rotation angles".upper())]
        zrot_lines += self._write_data_block(self.Z.rotation_angle, "zrot")

        # write out data only impedance and tipper
        z_data_lines = [self._data_header_str.format("impedances".upper())]
        self.Z.z = np.nan_to_num(self.Z.z)
        self.Z.z_err = np.nan_to_num(self.Z.z_err)
        self.Tipper.tipper = np.nan_to_num(self.Tipper.tipper)
        self.Tipper.tipper_err = np.nan_to_num(self.Tipper.tipper_err)
        for ii in range(2):
            for jj in range(2):
                z_lines_real = self._write_data_block(
                    self.Z.z[:, ii, jj].real, self._z_labels[2 * ii + jj][0]
                )
                z_lines_imag = self._write_data_block(
                    self.Z.z[:, ii, jj].imag, self._z_labels[2 * ii + jj][1]
                )
                z_lines_var = self._write_data_block(
                    self.Z.z_err[:, ii, jj] ** 2.0, self._z_labels[2 * ii + jj][2]
                )

                z_data_lines += z_lines_real
                z_data_lines += z_lines_imag
                z_data_lines += z_lines_var

        if self.Tipper.tipper is not None and np.all(self.Tipper.tipper == 0):
            trot_lines = [""]
            t_data_lines = [""]
        else:
            try:
                # write out rotation angles
                trot_lines = [
                    self._data_header_str.format("tipper rotation angles".upper())
                ]
                if isinstance(self.Tipper.rotation_angle, float):
                    trot = np.repeat(self.Tipper.rotation_angle, self.Tipper.freq.size)
                else:
                    trot = self.Tipper.rotation_angle
                trot_lines += self._write_data_block(np.array(trot), "trot")

                # write out tipper lines
                t_data_lines = [self._data_header_str.format("tipper".upper())]
                for jj in range(2):
                    t_lines_real = self._write_data_block(
                        self.Tipper.tipper[:, 0, jj].real, self._t_labels[jj][0]
                    )
                    t_lines_imag = self._write_data_block(
                        self.Tipper.tipper[:, 0, jj].imag, self._t_labels[jj][1]
                    )
                    t_lines_var = self._write_data_block(
                        self.Tipper.tipper_err[:, 0, jj] ** 2.0, self._t_labels[jj][2]
                    )

                    t_data_lines += t_lines_real
                    t_data_lines += t_lines_imag
                    t_data_lines += t_lines_var
            except AttributeError:
                trot_lines = [""]
                t_data_lines = [""]

        edi_lines = (
            header_lines
            + info_lines
            + define_lines
            + dsect_lines
            + freq_lines
            + zrot_lines
            + z_data_lines
            + trot_lines
            + t_data_lines
            + [">END"]
        )

        with open(new_edi_fn, "w") as fid:
            fid.write("".join(edi_lines))

        self.logger.info("Wrote {0}".format(new_edi_fn))
        return new_edi_fn

    def _write_data_block(self, data_comp_arr, data_key):
        """
        Write a data block

        :param data_comp_arr: array of data components
        :type data_comp_arr: np.ndarray

        :param data_key: the component to write out
        :type data_key: string

        :returns: list of lines to write to edi file
        :rtype: list
        """
        if data_key.lower().find("z") >= 0 and data_key.lower() not in ["zrot", "trot"]:
            block_lines = [
                ">{0} ROT=ZROT // {1:.0f}\n".format(
                    data_key.upper(), data_comp_arr.size
                )
            ]
        elif data_key.lower().find("t") >= 0 and data_key.lower() not in [
            "zrot",
            "trot",
        ]:
            block_lines = [
                ">{0} ROT=TROT // {1:.0f}\n".format(
                    data_key.upper(), data_comp_arr.size
                )
            ]
        elif data_key.lower() == "freq":
            block_lines = [
                ">{0} // {1:.0f}\n".format(data_key.upper(), data_comp_arr.size)
            ]

        elif data_key.lower() in ["zrot", "trot"]:
            block_lines = [
                ">{0} // {1:.0f}\n".format(data_key.upper(), data_comp_arr.size)
            ]

        else:
            raise MTex.MTpyError_EDI("Cannot write block for {0}".format(data_key))

        for d_index, d_comp in enumerate(data_comp_arr, 1):
            if d_comp == 0.0 and data_key.lower() not in ["zrot", "trot"]:
                d_comp = float(self.Header.empty)
            # write the string in the specified format
            num_str = "{0:{1}}".format(d_comp, self._num_format)

            # check to see if a new line is needed
            if d_index % self._block_len == 0:
                num_str += "\n"
            # at the end of the block add a return
            if d_index == data_comp_arr.size:
                num_str += "\n"

            block_lines.append(num_str)

        return block_lines

    # -----------------------------------------------------------------------
    # set a few important properties
    # --> Latitude
    @property
    def lat(self):
        """latitude in decimal degrees"""
        return self.Header.lat

    @lat.setter
    def lat(self, input_lat):
        """set latitude and make sure it is converted to a float"""

        self.Header.lat = gis_tools.assert_lat_value(input_lat)
        self.logger.info(
            "Converted input latitude to decimal degrees: {0: .6f}".format(
                self.Header.lat
            )
        )

    # --> Longitude
    @property
    def lon(self):
        """longitude in decimal degrees"""
        return self.Header.lon

    @lon.setter
    def lon(self, input_lon):
        """set latitude and make sure it is converted to a float"""
        self.Header.lon = gis_tools.assert_lon_value(input_lon)
        self.logger.info(
            "Converted input longitude to decimal degrees: {0: .6f}".format(
                self.Header.lon
            )
        )

    # --> Elevation
    @property
    def elev(self):
        """Elevation in elevation units"""
        return self.Header.elev

    @elev.setter
    def elev(self, input_elev):
        """set elevation and make sure it is converted to a float"""
        self.Header.elev = gis_tools.assert_elevation_value(input_elev)

    # --> station
    @property
    def station(self):
        """station name"""
        return self.Header.dataid

    @station.setter
    def station(self, new_station):
        """station name"""
        if not isinstance(new_station, str):
            new_station = f"{new_station}".replace(r"/", "")
        self.Header.dataid = new_station
        self.Data.sectid = new_station

    @property
    def survey_metadata(self):
        sm = metadata.Survey()
        sm.project = self.Header.project
        sm.survey_id = self.Header.survey
        sm.acquired_by.author = self.Header.acqby
        sm.geographic_name = self.Header.loc

        for key, value in self.Info.info_dict.items():
            key = key.lower()
            if key in ["project"]:
                setattr(sm, key, value)
            if key in ["survey"]:
                sm.survey_id = value

        return sm

    @property
    def station_metadata(self):
        sm = metadata.Station()
        sm.run_list.append(metadata.Run(id=f"{self.station}a"))
        sm.id = self.station
        sm.data_type = "MT"
        sm.channels_recorded = self.Measurement.channels_recorded
        # location
        sm.location.latitude = self.lat
        sm.location.longitude = self.lon
        sm.location.elevation = self.elev
        sm.location.datum = self.Header.datum
        sm.location.declination.value = self.Header.declination
        sm.orientation.reference_frame = self.Header.coordinate_system.split()[0]
        # provenance
        sm.acquired_by.author = self.Header.acqby
        sm.provenance.creation_time = self.Header.filedate
        sm.provenance.submitter.author = self.Header.fileby
        sm.provenance.software.name = self.Header.fileby
        sm.provenance.software.version = self.Header.progvers
        sm.transfer_function.processed_date = self.Header.filedate
        # dates
        if self.Header.acqdate is not None:
            sm.time_period.start = self.Header.acqdate

        for key, value in self.Info.info_dict.items():
            key = key.lower()
            if "processing" in key:
                key = key.split("processing")[1]
                if key in ["software"]:
                    sm.transfer_function.software.name = value
                elif key in ["tag"]:
                    if value.count(",") > 0:
                        sm.transfer_function.remote_references = value.split(",")
                    else:
                        sm.transfer_function.remote_references = value.split()
            
            elif key == "processedby":
                sm.transfer_function.processed_by.author = value
            
            elif key == "runlist":
                if value.count(",") > 0:
                    runs = value.split(",")
                else:
                    runs = value.split()
                sm.run_list = []
                for rr in runs:
                    sm.run_list.append(metadata.Run(id=rr))
                sm.transfer_function.runs_processed = runs
               
            elif key == "sitename":
                sm.geographic_name = value
            elif key == "signconvention":
                sm.transfer_function.sign_convention = value
            if 'mtft' in key or 'emtf' in key or 'mtedit' in key:
                sm.transfer_function.processing_parameters.append(f"{key}={value}")

        if self.Header.filedate is not None:
            sm.transfer_function.processed_date = self.Header.filedate
            
        # make any extra information in info list into a comment
        sm.comments = '\n'.join(self.Info.info_list)
        
        # add information to runs
        for rr in sm.run_list:
            rr.ex = self.ex_metadata
            rr.ey = self.ey_metadata
            rr.hx = self.hx_metadata
            rr.hy = self.hy_metadata
            if self.hz_metadata.component in ['hz']:
                rr.hz = self.hz_metadata
            if self.rrhx_metadata.component in ['rrhx']:
                rr.rrhx = self.rrhx_metadata
            if self.rrhy_metadata.component in ['rrhy']:
                rr.rrhy = self.rrhy_metadata
                

        return sm
    
    def _get_electric_metadata(self, comp):
        """
        get electric information from the various metadata
        """
        comp = comp.lower()
        electric = metadata.Electric()
        if hasattr(self.Measurement, f"meas_{comp}"):
            meas = getattr(self.Measurement, f"meas_{comp}")
            electric.dipole_length = meas.dipole_length
            electric.measurement_azimuth = meas.azimuth
            electric.component = meas.chtype
            electric.channel_number = meas.channel_number
            electric.negative.x = meas.x
            electric.positive.x2 = meas.x2
            electric.negative.y = meas.y
            electric.positive.y2 = meas.y2
            for k, v in self.Info.info_dict.items():
                if f'{comp}.' in k:
                    key = k.split(f'{comp}.')[1].strip()
                    if key  == 'manufacturer':
                        electric.negative.manufacturer = v
                        electric.positive.manufacturer = v
                    if key == 'type':
                        electric.negative.type = v
                        electric.positive.type = v
                        
            if electric.positive.x2 == 0 and electric.positive.y2 == 0.0 and\
                electric.negative.x == 0 and electric.negative.y == 0.0:
                electric.positive.x2 = electric.dipole_length * np.cos(np.deg2rad(meas.azimuth))
                electric.positive.y2 = electric.dipole_length * np.sin(np.deg2rad(meas.azimuth))

        return electric
        

    @property
    def ex_metadata(self):
        return self._get_electric_metadata('ex')

    @property
    def ey_metadata(self):
        return self._get_electric_metadata('ey')
    
    def _get_magnetic_metadata(self, comp):
        """
        
        get magnetic metadata from the various sources
        
        :param comp: DESCRIPTION
        :type comp: TYPE
        :return: DESCRIPTION
        :rtype: TYPE

        """
        
        magnetic = metadata.Magnetic()
        if hasattr(self.Measurement, f"meas_{comp}"):
            meas = getattr(self.Measurement, f"meas_{comp}")
            magnetic.measurement_azimuth = meas.azm
            magnetic.component = meas.chtype
            magnetic.channel_number = meas.channel_number
            magnetic.location.x = meas.x
            magnetic.location.y = meas.y
            try:
                magnetic.sensor.id = meas.meas_magnetic.sensor
            except AttributeError:
                pass
            for k, v in self.Info.info_dict.items():
                if f'{comp}.' in k:
                    key = k.split(f'{comp}.')[1].strip()
                    if key  == 'manufacturer':
                        magnetic.sensor.manufacturer = v
                    if key == 'type':
                        magnetic.sensor.type = v

        return magnetic

    @property
    def hx_metadata(self):
        return self._get_magnetic_metadata('hx')

    @property
    def hy_metadata(self):
        return self._get_magnetic_metadata('hy')

    @property
    def hz_metadata(self):
        return self._get_magnetic_metadata('hz')
    
    @property
    def rrhx_metadata(self):
        return self._get_magnetic_metadata('rrhx')
    @property
    def rrhy_metadata(self):
        return self._get_magnetic_metadata('rrhy')
    


# ==============================================================================
# Index finder
# ==============================================================================
class index_locator(object):
    def __init__(self, component_list):
        self.ex = None
        self.ey = None
        self.hx = None
        self.hy = None
        self.hz = None
        self.rhx = None
        self.rhy = None
        self.rhz = None
        for ii, comp in enumerate(component_list):
            setattr(self, comp, ii)
        if self.rhx is None:
            self.rhx = self.hx
        if self.rhy is None:
            self.rhy = self.hy


# ==============================================================================
#  Header object
# ==============================================================================
class Header(object):
    """
    Header class contains all the information in the header section of the .edi
    file. A typical header block looks like::

        >HEAD

            ACQBY=None
            ACQDATE=None
            DATAID=par28ew
            ELEV=0.000
            EMPTY=1e+32
            FILEBY=WG3DForward
            FILEDATE=2016/04/11 19:37:37 UTC
            LAT=-30:12:49
            LOC=None
            LON=139:47:50
            PROGDATE=2002-04-22
            PROGVERS=WINGLINK EDI 1.0.22
            COORDINATE SYSTEM = GEOGRAPHIC NORTH
            DECLINATION = 10.0


    :param fn: full path to .edi file to be read in.
                   *default* is None. If an .edi file is input attributes
                   of Header are filled.
    :type fn: string

    Many of the attributes are needed in the .edi file.  They are marked with
    a yes for 'In .edi'

    ============== ======================================= ======== ===========
    Attributes     Description                             Default  In .edi
    ============== ======================================= ======== ===========
    acqby          Acquired by                             None     yes
    acqdate        Acquired date (YYYY-MM-DD)              None     yes
    coordinate     [ geographic | geomagnetic ]            None     yes 
    dataid         Station name, should be a string        None     yes
    declination    geomagnetic declination                 None     yes 
    fn         Full path to .edi file                  None     no
    elev           Elevation of station (m)                None     yes
    empty          Value for missing data                  1e32     yes
    fileby         File written by                         None     yes
    filedate       Date the file is written (YYYY-MM-DD)   None     yes
    header_list    List of header lines                    None     no
    lat            Latitude of station [1]_                None     yes
    loc            Location name where station was         None     yes
                   collected
    lon            Longitude of station [1]_               None     yes
    phoenix_edi    [ True | False ] if phoenix .edi format False    no
    progdate       Date of program version to write .edi   None     yes
    progvers       Version of program writing .edi         None     yes
    stdvers        Standard version                        None     yes
    units          Units of distance                       m        yes
    _header_keys   list of metadata input into .edi        [2]_
                   header block.                                    no
    ============== ======================================= ======== ===========

    .. [1] Internally everything is converted to decimal degrees.  Output is
          written as HH:MM:SS.ss so Winglink can read them in.
    .. [2] If you want to change what metadata is written into the .edi file
           change the items in _header_keys.  Default attributes are:
               * acqby
               * acqdate
               * coordinate_system
               * dataid
               * declination
               * elev
               * fileby
               * lat
               * loc
               * lon
               * filedate
               * empty
               * progdate
               * progvers


    ====================== ====================================================
    Methods                Description
    ====================== ====================================================
    get_header_list        get header lines from edi file
    read_header            read in header information from header_lines
    write_header           write header lines, returns a list of lines to write
    ====================== ====================================================

    :Read Header: ::

        >>> import mtpy.core.edi as mtedi
        >>> header_obj = mtedi.Header(fn=r"/home/mt/mt01.edi")

    """

    def __init__(self, fn=None, **kwargs):
        self.logger = logging.getLogger(f"{__name__}.{self.__class__.__name__}")
        self._fn = None
        self.fn = fn
        self.edi_lines = None
        self.dataid = None
        self.acqby = None
        self.fileby = "MTpy"
        self._acqdate = MTime()
        self._enddate = None
        self._filedate = MTime()
        self.loc = None
        self._lat = None
        self._lon = None
        self._elev = None
        self.units = "[mV/km]/[nT]"
        self.empty = 1e32
        self.progvers = __version__
        self._progdate = MTime("2020-11-10")
        self.progname = "MTpy"
        self.project = None
        self.survey = None
        self.coordinate_system = "Geographic North"
        self.declination = None
        self.datum = "WGS84"
        self.phoenix_edi = False
        self.stdvers = "SEG 1.0"
        self.state = None
        self.country = None

        self.header_list = None

        self._header_keys = [
            "acqby",
            "acqdate",
            "dataid",
            "elev",
            "fileby",
            "lat",
            "loc",
            "lon",
            "filedate",
            "empty",
            "progname",
            "progdate",
            "progvers",
            "coordinate_system",
            "declination",
            "datum",
            "project",
            "survey",
            "units",
            "stdvers",
        ]
        
        self._optional_keys = ["enddate", "state", "country"]

        for key in list(kwargs.keys()):
            setattr(self, key, kwargs[key])

        if self.fn is not None or self.edi_lines is not None:
            self.read_header()

    def __str__(self):
        return "".join(self.write_header())

    def __repr__(self):
        return self.__str__()
    
    @property
    def fn(self):
        return self._fn
    
    @fn.setter
    def fn(self, value):
        if value is None:
            self._fn = None
            return 
        self._fn = Path(value) 
        if self._fn.exists():
            self.read_header()
        
    @property
    def lat(self):
        return self._lat
    
    @lat.setter
    def lat(self, value):
        self._lat = gis_tools.assert_lat_value(value)
        
    @property
    def lon(self):
        return self._lon
    
    @lon.setter
    def lon(self, value):
        self._lon = gis_tools.assert_lon_value(value)
        
    @property
    def elev(self):
        return self._elev
    
    @elev.setter
    def elev(self, value):
        self._elev = gis_tools.assert_elevation_value(value)
        
    @property
    def acqdate(self):
        return self._acqdate.date
    
    @acqdate.setter
    def acqdate(self, value):
        self._acqdate = MTime(value)
        
    @property
    def enddate(self):
        if self._enddate is not None:
            return self._enddate.date
    
    @enddate.setter
    def enddate(self, value):
        self._enddate = MTime(value)
        
    @property
    def filedate(self):
        return self._filedate.date
    
    @filedate.setter
    def filedate(self, value):
        self._filedate = MTime(value)
        
    @property
    def progdate(self):
        return self._progdate.date
    
    @progdate.setter
    def progdate(self, value):
        self._progdate = MTime(value)
            

    def get_header_list(self):
        """
        Get the header information from the .edi file in the form of a list,
        where each item is a line in the header section.
        """

        if self.fn is None and self.edi_lines is None:
            self.logger.info("No edi file to read.")
            return

        header_list = []
        head_find = False

        # read in file line by line
        if self.fn is not None:
            if not self.fn.exists():
                msg = f"Could not find {self.fn}"
                self.logger.error(msg)
                raise IOError(msg)
            with open(self.fn, "r") as fid:
                self.edi_lines = _validate_edi_lines(fid.readlines())

        # read in list line by line and then truncate
        for line in self.edi_lines:
            # check for header label
            if ">" in line and "head" in line.lower():
                head_find = True
            # if the header line has been found then the next >
            # should be the next section so stop
            elif ">" in line:
                if head_find is True:
                    break
                else:
                    pass
            # get the header information into a list
            elif head_find:
                # skip any blank lines
                if len(line.strip()) > 2:
                    line = line.strip().replace('"', "")
                    h_list = line.split("=")
                    if len(h_list) == 2:
                        key = h_list[0].strip()
                        value = h_list[1].strip()
                        header_list.append("{0}={1}".format(key, value))
                        
        return header_list

    def read_header(self, header_list=None):
        """
        read a header information from either edi file or a list of lines
        containing header information.

        :param header_list: should be read from an .edi file or input as
                            ['key_01=value_01', 'key_02=value_02']
        :type header_list: list

        :Input header_list: ::

            >>> h_list = ['lat=36.7898', 'lon=120.73532', 'elev=120.0', ...
            >>>           'dataid=mt01']
            >>> import mtpy.core.edi as mtedi
            >>> header = mtedi.Header()
            >>> header.read_header(h_list)

        """

        if header_list is not None:
            self.header_list = self._validate_header_list(header_list)

        if self.header_list is None and self.fn is None and self.edi_lines is None:
            self.logger.info("Nothing to read. header_list and fn are None")

        elif self.fn is not None or self.edi_lines is not None:
            self.header_list = self.get_header_list()

        for h_line in self.header_list:
            h_list = h_line.split("=")
            key = h_list[0].lower()
            value = h_list[1]
            # test if its a phoenix formated .edi file
            if key in ["progvers"]:
                if value.lower().find("mt-editor") != -1:
                    self.phoenix_edi = True
            if key in "latitude":
                key = "lat"
            elif key in "longitude":
                key = "lon"
            elif key in "elevation":
                key = "elev"
            elif key in "location":
                key = "loc"

            setattr(self, key, value)
            
            # be sure to pass any uncommon keys through to new file
            if key not in self._header_keys:
                self._optional_keys.append(key)

    def write_header(
        self, header_list=None, longitude_format="LON", latlon_format="dms"
    ):
        """
        Write header information to a list of lines.


        :param header_list: should be read from an .edi file or input as
                            ['key_01=value_01', 'key_02=value_02']
        :type header_list: list
        :param longitude_format:  whether to write longitude as LON or LONG. 
                                  options are 'LON' or 'LONG', default 'LON'
        :type longitude_format:  string
        :param latlon_format:  format of latitude and longitude in output edi,
                               degrees minutes seconds ('dms') or decimal 
                               degrees ('dd')
        :type latlon_format:  string

        :returns header_lines: list of lines containing header information
                               will be of the form::

                               ['>HEAD\n',
                                '    key_01=value_01\n']
                                if None is input then reads from input .edi
                                file or uses attribute information to write
                                metadata.

        """

        if header_list is not None:
            self.read_header(header_list)

        if self.header_list is None and self.fn is not None:
            self.header_list = self.get_header_list()

        header_lines = [">HEAD\n"]
        for key in sorted(self._header_keys + self._optional_keys):
            value = getattr(self, key)
            if key in self._optional_keys and value is None:
                continue
            if key in ["lat", "lon"] and value is not None:
                if latlon_format.upper() == "DD":
                    value = "%.6f" % value
                else:
                    value = gis_tools.convert_position_float2str(value)
            if key in ["elev", "declination"] and value is not None:
                try:
                    value = "{0:.3f}".format(value)
                except ValueError:
                    raise Exception("value error for key elev or declination")
                    # value = '0.000'

            if key in ["filedate"]:
                value = get_now_utc()
            #
            if key == "lon":
                if longitude_format == "LONG":
                    key = "long"
            if isinstance(value, list):
                value = ",".join(value)

            header_lines.append(f"{tab}{key.upper()}={value}\n")
            
            
        header_lines.append("\n")
        return header_lines

    def _validate_header_list(self, header_list):
        """
        make sure the input header list is valid

        returns a validated header list
        """

        if header_list is None:
            self.logger.info("No header information to read")
            return None

        new_header_list = []
        for h_line in header_list:
            h_line = h_line.strip().replace('"', "")
            if len(h_line) > 1:
                h_list = h_line.split("=")
                if len(h_list) == 2:
                    key = h_list[0].strip().lower()
                    value = h_list[1].strip()
                    new_header_list.append("{0}={1}".format(key, value))

        return new_header_list


# ==============================================================================
# Info object
# ==============================================================================
class Information(object):
    """
    Contain, read, and write info section of .edi file

    not much to really do here, but just keep it in the same format that it is
    read in as, except if it is in phoenix format then split the two paragraphs
    up so they are sequential.

    """

    def __init__(self, fn=None, edi_lines=None):
        self.logger = logging.getLogger(f"{__name__}.{self.__class__.__name__}")
        self._fn = None
        self.fn = fn
        self.edi_lines = edi_lines
        self.info_list = []
        self.info_dict = {}

        if self.fn is not None or self.edi_lines is not None:
            self.read_info()
            
    @property
    def fn(self):
        return self._fn
    
    @fn.setter
    def fn(self, value):
        if value is None:
            self._fn = None
            return 
        self._fn = Path(value) 
        if self._fn.exists():
            self.read_info()
                
    def __str__(self):
        return "".join(self.write_info())

    def __repr__(self):
        return self.__str__()

    def get_info_list(self):
        """
        get a list of lines from the info section
        """

        if self.fn is None and self.edi_lines is None:
            self.logger.info("No EDI file input")
            return

        self.info_list = []
        info_find = False
        phoenix_file = False
        phoenix_list_02 = []

        if self.fn is not None:
            if os.path.isfile(self.fn) is False:
                msg = f"Could not find EDI file: {self.fn}"
                self.logger.info(msg)
                return

            with open(self.fn, "r") as fid:
                self.edi_lines = _validate_edi_lines(fid.readlines())

        for line in self.edi_lines:
            if ">" in line and "info" in line.lower():
                info_find = True
            elif ">" in line:
                # need to check for xml type formating
                if "<" in line:
                    pass
                else:
                    if info_find is True:
                        break
                    else:
                        pass
            elif info_find:
                if line.lower().find("run information") >= 0:
                    phoenix_file = True
                if phoenix_file and len(line) > 40:
                    self.info_list.append(line[0:37].strip())
                    phoenix_list_02.append(line[38:].strip())
                else:
                    if len(line.strip()) > 1:
                        self.info_list.append(line.strip())

        self.info_list += phoenix_list_02
        # validate the information list
        self.info_list = self._validate_info_list(self.info_list)

    def read_info(self, info_list=None):
        """
        read information section of the .edi file
        """

        if info_list is not None:
            self.info_list = self._validate_info_list(info_list)

        elif self.fn is not None or self.edi_lines is not None:
            self.get_info_list()

        self.info_dict = {}
        colon_find = None
        equals_find = None
        # make info items attributes of Information
        for ll in self.info_list:
            l_list = [None, ""]
            # phoenix has lat an lon information in the notes but separated by 
            # a space instead of an = or :
            if "lat" in ll.lower() or "lon" in ll.lower() or "lng" in ll.lower():
                l_list = ll.split()
                if len(l_list) == 2:
                    self.info_dict[l_list[0]] = l_list[1]
                    continue
                elif len(l_list) == 4:
                    self.info_dict[l_list[0]] = l_list[1]
                    self.info_dict[l_list[2]] = l_list[3]
                    continue
                elif len(l_list) == 6:
                    self.info_dict[l_list[0]] = l_list[1] + l_list[2]
                    self.info_dict[l_list[3]] = l_list[4] + l_list[5]
                    continue
                
            # need to check if there is an = or : seperator, which ever
            # comes first is assumed to be the delimiter
            if ll.find(":") > 0:
                colon_find = ll.find(':')
            if ll.find("=") > 0:
                equals_find = ll.find('=')
            if colon_find is not None and equals_find is not None:
                if ll.find(":") < ll.find("="):
                    l_list = ll.split(":")
                else:
                    l_list = ll.split("=")
            elif colon_find is not None:
                l_list = ll.split(":")
            elif equals_find is not None:
                l_list = ll.split("=")
            else:
                l_list[0] = ll
            if l_list[0] is not None and len(l_list) > 1:
                l_key = l_list[0]
                l_value = l_list[1].strip()
                    
                self.info_dict[l_key] = l_value.replace('"', "")
            else: 
                self.info_dict[l_list[0]] = None

        if self.info_list is None:
            self.logger.info("Could not read information")
            return

    def write_info(self, info_list=None):
        """
        write out information
        """

        if info_list is not None:
            self.info_list = self._validate_info_list(info_list)

        info_lines = [">INFO\n"]
        for line in self.info_list:
            info_lines.append("{0}{1}\n".format(tab, line))

        return info_lines

    def _validate_info_list(self, info_list):
        """
        check to make sure the info list input is valid, really just checking
        for Phoenix format where they put two columns in the file and remove
        any blank lines and the >info line
        """

        new_info_list = []
        for line in info_list:
            # get rid of empty lines
            lt = str(line).strip()
            if len(lt) > 1:
                if ">" in line:
                    pass
                else:
                    new_info_list.append(line.strip())

        return new_info_list


# ==============================================================================
#  Define measurement class
# ==============================================================================
class DefineMeasurement(object):
    """
    DefineMeasurement class holds information about the measurement.  This
    includes how each channel was setup.  The main block contains information
    on the reference location for the station.  This is a bit of an archaic
    part and was meant for a multiple station .edi file.  This section is also
    important if you did any forward modeling with Winglink cause it only gives
    the station location in this section.  The other parts are how each channel
    was collected.  An example define measurement section looks like::

        >=DEFINEMEAS

            MAXCHAN=7
            MAXRUN=999
            MAXMEAS=9999
            UNITS=M
            REFTYPE=CART
            REFLAT=-30:12:49.4693
            REFLONG=139:47:50.87
            REFELEV=0

        >HMEAS ID=1001.001 CHTYPE=HX X=0.0 Y=0.0 Z=0.0 AZM=0.0
        >HMEAS ID=1002.001 CHTYPE=HY X=0.0 Y=0.0 Z=0.0 AZM=90.0
        >HMEAS ID=1003.001 CHTYPE=HZ X=0.0 Y=0.0 Z=0.0 AZM=0.0
        >EMEAS ID=1004.001 CHTYPE=EX X=0.0 Y=0.0 Z=0.0 X2=0.0 Y2=0.0
        >EMEAS ID=1005.001 CHTYPE=EY X=0.0 Y=0.0 Z=0.0 X2=0.0 Y2=0.0
        >HMEAS ID=1006.001 CHTYPE=HX X=0.0 Y=0.0 Z=0.0 AZM=0.0
        >HMEAS ID=1007.001 CHTYPE=HY X=0.0 Y=0.0 Z=0.0 AZM=90.0

    :param fn: full path to .edi file to read in.
    :type fn: string

    ================= ==================================== ======== ===========
    Attributes        Description                          Default  In .edi
    ================= ==================================== ======== ===========
    fn            Full path to edi file read in        None     no
    maxchan           Maximum number of channels measured  None     yes
    maxmeas           Maximum number of measurements       9999     yes
    maxrun            Maximum number of measurement runs   999      yes
    meas_####         HMeasurement or EMEasurment object   None     yes
                      defining the measurement made [1]_
    refelev           Reference elevation (m)              None     yes
    reflat            Reference latitude [2]_              None     yes
    refloc            Reference location                   None     yes
    reflon            Reference longituted [2]_            None     yes
    reftype           Reference coordinate system          'cart'   yes
    units             Units of length                      m        yes
    _define_meas_keys Keys to include in define_measurment [3]_     no
                      section.
    ================= ==================================== ======== ===========

    .. [1] Each channel with have its own define measurement and depending on
           whether it is an E or H channel the metadata will be different.
           the #### correspond to the channel number.
    .. [2] Internally everything is converted to decimal degrees.  Output is
          written as HH:MM:SS.ss so Winglink can read them in.
    .. [3] If you want to change what metadata is written into the .edi file
           change the items in _header_keys.  Default attributes are:
               * maxchan
               * maxrun
               * maxmeas
               * reflat
               * reflon
               * refelev
               * reftype
               * units

    """

    def __init__(self, fn=None, edi_lines=None):
        self.logger = logging.getLogger(f"{__name__}.{self.__class__.__name__}")
        self._fn = None
        self.fn = fn
        self.edi_lines = edi_lines
        self.measurement_list = None

        self.maxchan = None
        self.maxmeas = 7
        self.maxrun = 999
        self.refelev = None
        self.reflat = None
        self.reflon = None
        self.reftype = "cartesian"
        self.units = "m"

        self._define_meas_keys = [
            "maxchan",
            "maxrun",
            "maxmeas",
            "reflat",
            "reflon",
            "refelev",
            "reftype",
            "units",
        ]

        if self.edi_lines is not None:
            self.read_measurement()
            
    def __str__(self):
        return "".join(self.write_measurement())
    
    def __repr__(self):
        return self.__str__()
    
    @property
    def fn(self):
        return self._fn
    
    @fn.setter
    def fn(self, value):
        if value is None:
            self._fn = None
            return 
        self._fn = Path(value) 
        if self._fn.exists():
            self.read_measurement()

    def get_measurement_lists(self):
        """
        get measurement list including measurement setup
        """
        if self.fn is None and self.edi_lines is None:
            self.logger.info("No edi file input, check fn attribute")
            return

        self.measurement_list = []
        meas_find = False
        count = 0

        if self.fn is not None:
            if os.path.isfile(self.fn) is False:
                self.logger.info("Could not find {0}, check path".format(self.fn))
                return

            with open(self.fn, "r") as fid:
                self.edi_lines = _validate_edi_lines(fid.readlines())

        for line in self.edi_lines:
            if ">=" in line and "definemeas" in line.lower():
                meas_find = True
            elif ">=" in line:
                if meas_find is True:
                    break
            elif meas_find is True and ">" not in line:
                line = line.strip()
                if len(line) > 2:
                    if count > 0:
                        line_list = _validate_str_with_equals(line)
                        for ll in line_list:
                            ll_list = ll.split("=")
                            key = ll_list[0].lower()
                            value = ll_list[1]
                            self.measurement_list[-1][key] = value
                    else:
                        self.measurement_list.append(line.strip())

            # look for the >XMEAS parts
            elif ">" in line and meas_find:
                if line.find("!") > 0:
                    pass
                else:
                    count += 1
                    line_list = _validate_str_with_equals(line)
                    m_dict = {}
                    for ll in line_list:
                        ll_list = ll.split("=")
                        key = ll_list[0].lower()
                        value = ll_list[1]
                        m_dict[key] = value
                    self.measurement_list.append(m_dict)

    def read_measurement(self, measurement_list=None):
        """
        read the define measurment section of the edi file

        should be a list with lines for:
            - maxchan
            - maxmeas
            - maxrun
            - refelev
            - reflat
            - reflon
            - reftype
            - units
            - dictionaries for >XMEAS with keys:
                - id
                - chtype
                - x
                - y
                - axm
                -acqchn

        """

        if measurement_list is not None:
            self.measurement_list = measurement_list

        elif self.fn is not None or self.edi_lines is not None:
            self.get_measurement_lists()

        if self.measurement_list is None:
            self.logger.info("Nothing to read, check fn or measurement_list attributes")
            return

        for line in self.measurement_list:
            if isinstance(line, str):
                line_list = line.split("=")
                key = line_list[0].lower()
                value = line_list[1].strip()
                if key in "reflatitude":
                    key = "reflat"
                    value = gis_tools.assert_lat_value(value)
                elif key in "reflongitude":
                    key = "reflon"
                    value = gis_tools.assert_lon_value(value)
                elif key in "refelevation":
                    key = "refelev"
                    value = gis_tools.assert_elevation_value(value)
                elif key in "maxchannels":
                    key = "maxchan"
                    try:
                        value = int(value)
                    except ValueError:
                        value = 0
                elif key in "maxmeasurements":
                    key = "maxmeas"
                    try:
                        value = int(value)
                    except ValueError:
                        value = 0
                elif key in "maxruns":
                    key = "maxrun"
                    try:
                        value = int(value)
                    except ValueError:
                        value = 0
                setattr(self, key, value)

            elif isinstance(line, dict):
                key = "meas_{0}".format(line["chtype"].lower())
                if key[4:].find("h") >= 0:
                    value = HMeasurement(**line)
                elif key[4:].find("e") >= 0:
                    value = EMeasurement(**line)
                if hasattr(self, key):
                    key = key.replace('_', '_rr')
                    try:
                        value.chtype = f'RR{value.chtype}'
                    except AttributeError:
                        pass
                setattr(self, key, value)

    def write_measurement(
        self, measurement_list=None, longitude_format="LON", latlon_format="dd"
    ):
        """
        write the define measurement block as a list of strings
        """

        if measurement_list is not None:
            self.read_measurement(measurement_list=measurement_list)

        measurement_lines = ["\n>=DEFINEMEAS\n"]
        for key in self._define_meas_keys:
            value = getattr(self, key)
            if key == "reflat" or key == "reflon":
                if latlon_format.upper() == "DD":
                    value = "%.6f" % value
                else:
                    value = gis_tools.convert_position_float2str(value)
            elif key == "refelev":
                value = "{0:.3f}".format(gis_tools.assert_elevation_value(value))
            if key.upper() == "REFLON":
                if longitude_format == "LONG":
                    key += "G"
            measurement_lines.append("{0}{1}={2}\n".format(tab, key.upper(), value))
        measurement_lines.append("\n")

        # need to write the >XMEAS type, but sort by channel number
        m_key_list = [
            (kk.strip(), float(self.__dict__[kk].id))
            for kk in list(self.__dict__.keys())
            if kk.find("meas_") == 0
        ]
        if len(m_key_list) == 0:
            self.logger.info("No XMEAS information.")
        else:
            # need to sort the dictionary by chanel id
            chn_count = 1
            for x_key in sorted(m_key_list, key=lambda x: x[1]):
                x_key = x_key[0]
                m_obj = getattr(self, x_key)
                if m_obj.chtype is not None:
                    m_obj.chtype = str(m_obj.chtype)
                    if m_obj.chtype.lower().find("h") >= 0:
                        head = "hmeas"
                    elif m_obj.chtype.lower().find("e") >= 0:
                        head = "emeas"
                    else:
                        head = "None"
                else:
                    head = "None"

                m_list = [">{0}".format(head.upper())]

                for mkey, mfmt in zip(m_obj._kw_list, m_obj._fmt_list):
                    if mkey == "acqchan":
                        if (
                            getattr(m_obj, mkey) is None
                            or getattr(m_obj, mkey) == "None"
                        ):
                            setattr(m_obj, mkey, chn_count)
                            chn_count += 1

                    try:
                        m_list.append(
                            " {0}={1:{2}}".format(
                                mkey.upper(), getattr(m_obj, mkey), mfmt
                            )
                        )
                    except ValueError:
                        m_list.append(" {0}={1:{2}}".format(mkey.upper(), 0.0, mfmt))

                m_list.append("\n")
                measurement_lines.append("".join(m_list))

        return measurement_lines

    def get_measurement_dict(self):
        """
        get a dictionary for the xmeas parts
        """
        meas_dict = {}
        for key in list(self.__dict__.keys()):
            if key.find("meas_") == 0:
                meas_attr = getattr(self, key)
                meas_key = meas_attr.chtype
                meas_dict[meas_key] = meas_attr

        return meas_dict
    
    def from_metadata(self, channel):
        """
        create a measurement class from metadata
        
        :param channel: DESCRIPTION
        :type channel: TYPE
        :return: DESCRIPTION
        :rtype: TYPE

        """
        
        if 'e' in channel.component:
            meas = EMeasurement(
                **{
                    "x": channel.negative.x,
                    "x2": channel.positive.x2,
                    "y": channel.negative.y,
                    "y2": channel.positive.y2,
                    "chtype": channel.component,
                    "id": channel.channel_number,
                    "acqchan": channel.channel_number,
                }
            )
            setattr(self, f"meas_{channel.component.lower()}", meas)
            
        if 'h' in channel.component:
            meas = HMeasurement(
                **{
                    "x": channel.location.x,
                    "y": channel.location.y,
                    "azm": channel.measurement_azimuth,
                    "chtype": channel.component,
                    "id": channel.channel_number,
                    "acqchan": channel.channel_number,
                }
            )
            setattr(self, f"meas_{channel.component.lower()}", meas)

    @property
    def channels_recorded(self):
        """ Get the channels recorded """

        return [cc.lower() for cc in self.get_measurement_dict().keys()]


# ==============================================================================
# magnetic measurements
# ==============================================================================
class HMeasurement(object):
    """
    HMeasurement contains metadata for a magnetic field measurement

    ====================== ====================================================
    Attributes             Description
    ====================== ====================================================
    id                     Channel number
    chtype                 [ HX | HY | HZ | RHX | RHY ]
    x                      x (m) north from reference point (station)
    y                      y (m) east from reference point (station)
    azm                    angle of sensor relative to north = 0
    acqchan                name of the channel acquired usually same as chtype
    ====================== ====================================================

    :Fill Metadata: ::

        >>> import mtpy.core.edi as mtedi
        >>> h_dict = {'id': '1', 'chtype':'hx', 'x':0, 'y':0, 'azm':0}
        >>> h_dict['acqchn'] = 'hx'
        >>> hmeas = mtedi.HMeasurement(**h_dict)
    """

    def __init__(self, **kwargs):

        self._kw_list = ["id", "chtype", "x", "y", "azm", "acqchan"]
        self._fmt_list = ["<4.0f", "<3", "<4.1f", "<4.1f", "<4.1f", "<4"]
        for key, fmt in zip(self._kw_list, self._fmt_list):
            if "f" in fmt:
                setattr(self, key.lower(), 0.0)
            else:
                setattr(self, key.lower(), "0.0")

        for key, value in kwargs.items():
            if value is None:
                value = 0.0
            try:
                setattr(self, key.lower(), float(value))
            except ValueError:
                setattr(self, key.lower(), value)
                
    def __str__(self):
        return "\n".join(
            [f"{k} = {v}" for k, v in self.__dict__.items() if k[0] != "_"]
        )

    def __repr__(self):
        return self.__str__()

    @property
    def channel_number(self):
        if not isinstance(self.acqchan, (int, float)):
            try:
                return [int("".join(i for i in self.acqchan if i.isdigit()))][0]
            except IndexError:
                return 0
        return self.acqchan


# ==============================================================================
# electric measurements
# ==============================================================================
class EMeasurement(object):
    """
    EMeasurement contains metadata for an electric field measurement

    ====================== ====================================================
    Attributes             Description
    ====================== ====================================================
    id                     Channel number
    chtype                 [ EX | EY ]
    x                      x (m) north from reference point (station) of one
                           electrode of the dipole
    y                      y (m) east from reference point (station) of one
                           electrode of the dipole
    x2                     x (m) north from reference point (station) of the
                           other electrode of the dipole
    y2                     y (m) north from reference point (station) of the
                           other electrode of the dipole
    acqchan                name of the channel acquired usually same as chtype
    ====================== ====================================================

    :Fill Metadata: ::

        >>> import mtpy.core.edi as mtedi
        >>> e_dict = {'id': '1', 'chtype':'ex', 'x':0, 'y':0, 'x2':50, 'y2':50}
        >>> e_dict['acqchn'] = 'ex'
        >>> emeas = mtedi.EMeasurement(**e_dict)
    """

    def __init__(self, **kwargs):
        self._logger = logging.getLogger(f"{__name__}.{self.__class__.__name__}")
        self._kw_list = ["id", "chtype", "x", "y", "x2", "y2", "acqchan"]
        self._fmt_list = ["<4.4g", "<3", "<4.1f", "<4.1f", "<4.1f", "<4.1f", "<4"]
        for key, fmt in zip(self._kw_list, self._fmt_list):
            if "f" in fmt:
                setattr(self, key.lower(), 0.0)
            else:
                setattr(self, key.lower(), "0.0")

        for key, value in kwargs.items():
            if value is None:
                value = 0.0
            try:
                setattr(self, key.lower(), float(value))
            except ValueError:
                setattr(self, key.lower(), value)

    def __str__(self):
        return "\n".join(
            [f"{k} = {v}" for k, v in self.__dict__.items() if k[0] != "_"]
        )

    def __repr__(self):
        return self.__str__()

    @property
    def dipole_length(self):
        if hasattr(self, "z"):
            return np.sqrt(
                (self.x2 - self.x) ** 2 + (self.y2 - self.y) ** 2 + self.z ** 2
            )
        else:
            return np.sqrt((self.x2 - self.x) ** 2 + (self.y2 - self.y) ** 2)

    @property
    def azimuth(self):
        if hasattr(self, "azm"):
            return self.azm
        try:
            return np.rad2deg(np.arctan2((self.y2 - self.y), (self.x2 - self.x)))
        except ZeroDivisionError:
            return 0.0

    @property
    def channel_number(self):
        if not isinstance(self.acqchan, (int, float)):
            try:
                return [int("".join(i for i in self.acqchan if i.isdigit()))][0]
            except IndexError:
                return 0
        return self.acqchan


# ==============================================================================
# data section
# ==============================================================================
class DataSection(object):
    """
    DataSection contains the small metadata block that describes which channel
    is which.  A typical block looks like::

        >=MTSECT

            ex=1004.001
            ey=1005.001
            hx=1001.001
            hy=1002.001
            hz=1003.001
            nfreq=14
            sectid=par28ew
            nchan=None
            maxblks=None


    :param fn: full path to .edi file to read in.
    :type fn: string


    ================= ==================================== ======== ===========
    Attributes        Description                          Default  In .edi
    ================= ==================================== ======== ===========
    ex                ex channel id number                 None     yes
    ey                ey channel id number                 None     yes
    hx                hx channel id number                 None     yes
    hy                hy channel id number                 None     yes
    hz                hz channel id number                 None     yes
    nfreq             number of frequencies                None     yes
    sectid            section id, should be the same
                      as the station name -> Header.dataid None     yes
    maxblks           maximum number of data blocks        None     yes
    nchan             number of channels                   None     yes
    _kw_list          list of key words to put in metadata [1]_     no
    ================= ==================================== ======== ===========

    .. [1] Changes these values to change what is written to edi file
    """

    def __init__(self, fn=None, edi_lines=None):
        """
        writing the EDI files MTSECT
        :param fn:
        :param edi_lines:
        """
        self.logger = logging.getLogger(f"{__name__}.{self.__class__.__name__}")
        self.fn = fn
        self.edi_lines = edi_lines

        self.data_type = "z"
        self.line_num = 0
        self.Data_list = None

        self.nfreq = None
        self.sectid = None
        self.nchan = None
        self.maxblks = None
        self.ex = None
        self.ey = None
        self.hx = None
        self.hy = None
        self.hz = None

        self._kw_list = [
            "nfreq",
            "sectid",
            "nchan",
            "maxblks",
            "ex",
            "ey",
            "hx",
            "hy",
            "hz",
        ]

        if self.fn is not None or self.edi_lines is not None:
            self.read_Data()
            
    def __str__(self):
        return ''.join(self.write_Data())
    
    def __repr__(self):
        return self.__str__()

    def get_Data(self):
        """
        read in the data of the file, will detect if reading spectra or
        impedance.
        """

        if self.fn is None and self.edi_lines is None:
            raise MTex.MTpyError_EDI("No edi file to read. Check fn")

        self.Data_list = []
        Data_find = False

        if self.fn is not None:
            if os.path.isfile(self.fn) is False:
                raise MTex.MTpyError_EDI(
                    "Could not find {0}. Check path".format(self.fn)
                )
            with open(self.fn) as fid:
                self.edi_lines = _validate_edi_lines(fid.readlines())

        for ii, line in enumerate(self.edi_lines):
            if ">=" in line and "sect" in line.lower():
                Data_find = True
                self.line_num = ii
                if line.lower().find("spect") > 0:
                    self.data_type = "spectra"
                elif line.lower().find("mt") > 0:
                    self.data_type = "z"
            elif ">" in line and Data_find is True:
                self.line_num = ii
                break

            elif Data_find:
                if len(line.strip()) > 2:
                    self.Data_list.append(line.strip())

    def read_Data(self, Data_list=None):
        """
        read data section
        """

        if Data_list is not None:
            self.Data_list = Data_list

        elif self.fn is not None or self.edi_lines is not None:
            self.get_Data()

        for d_line in self.Data_list:
            d_list = d_line.split("=")
            if len(d_list) > 1:
                key = d_list[0].lower()
                try:
                    value = int(d_list[1].strip())
                except ValueError:
                    value = d_list[1].strip().replace('"', "")

                setattr(self, key, value)

    def write_Data(self, Data_list=None, over_dict=None):
        """
        write a data section
        """

        # FZ: need to modify the nfreq (number of freqs), when re-writing effective EDI files)
        if over_dict is not None:
            for akey in list(over_dict.keys()):
                self.__setattr__(akey, over_dict[akey])

        if Data_list is not None:
            self.read_Data(Data_list)

        self.data_type = "z"
        self.logger.info("Writing out data a impedances")

        Data_lines = ["\n>=mtsect\n".upper()]

        for key in self._kw_list[0:4]:
            Data_lines.append(
                "{0}{1}={2}\n".format(tab, key.upper(), getattr(self, key))
            )

        # need to sort the list so it is descending order by channel number
        ch_list = [(key.upper(), getattr(self, key)) for key in self._kw_list[4:]]
        # ch_list = sorted(ch_list, key=lambda x: x[1])  #FZ: x[1] can be None, not working for Py3
        ch_list2 = sorted(ch_list, key=lambda x: x[0])

        for ch in ch_list2:
            Data_lines.append("{0}{1}={2}\n".format(tab, ch[0], ch[1]))

        Data_lines.append("\n")

        return Data_lines


def _validate_str_with_equals(input_string):
    """
    make sure an input string is of the format {0}={1} {2}={3} {4}={5} ...
    Some software programs put spaces after the equals sign and that's not
    cool.  So we make the string into a readable format

    :param input_string: input string from an edi file
    :type input_string: string

    :returns line_list: list of lines as ['key_00=value_00',
                                          'key_01=value_01']
    :rtype line_list: list
    """
    input_string = input_string.strip()
    # remove the first >XXXXX
    if ">" in input_string:
        input_string = input_string[input_string.find(" ") :]

    # check if there is a // at the end of the line
    if input_string.find("//") > 0:
        input_string = input_string[0 : input_string.find("//")]

    # split the line by =
    l_list = input_string.strip().split("=")

    # split the remaining strings
    str_list = []
    for line in l_list:
        s_list = line.strip().split()
        for l_str in s_list:
            str_list.append(l_str.strip())

    # probably not a good return
    if len(str_list) % 2 != 0:
        # _logger.info(
        #     'The number of entries in {0} is not even'.format(str_list))
        return str_list

    line_list = [
        "{0}={1}".format(str_list[ii], str_list[ii + 1])
        for ii in range(0, len(str_list), 2)
    ]

    return line_list


def _validate_edi_lines(edi_lines):
    """
    check for carriage returns or hard returns

    :param edi_lines: list of edi lines
    :type edi_lines: list

    :returns: list of edi lines
    :rtype: list
    """

    if len(edi_lines) == 1:
        edi_lines = edi_lines[0].replace("\r", "\n").split("\n")
        if len(edi_lines) > 1:
            return edi_lines
        else:
            raise ValueError("*** EDI format not correct check file ***")
    else:
        return edi_lines


def read_edi(fn):
    """
    
    Read an edi file and return a :class:`mtpy.core.mt.MT` object
    
    :param fn: DESCRIPTION
    :type fn: TYPE
    :return: DESCRIPTION
    :rtype: TYPE

    """

    edi_obj = Edi()
    edi_obj.read_edi_file(fn)

    mt_obj = mt.MT()

    for attr in [
        "Z",
        "Tipper",
        "survey_metadata",
        "station_metadata",
    ]:
        setattr(mt_obj, attr, getattr(edi_obj, attr))

    return mt_obj


def write_edi(mt_object):
    """
    Write an edi file from an :class:`mtpy.core.mt.MT` object
    
    :param mt_object: DESCRIPTION
    :type mt_object: TYPE
    :return: DESCRIPTION
    :rtype: TYPE

    """

    if not isinstance(mt_object, mt.MT):
        raise ValueError("Input must be an mtpy.core.mt.MT object")

    edi_obj = Edi()
    edi_obj.Z = mt_object.Z
    edi_obj.Tipper = mt_object.Tipper
    
    ### fill header information from survey
    edi_obj.Header.survey = mt_object.survey_metadata.survey_id
    edi_obj.Header.project = mt_object.survey_metadata.project
    edi_obj.Header.loc = mt_object.survey_metadata.geographic_name

    ### fill header information from station
    edi_obj.Header.acqby = mt_object.station_metadata.acquired_by.author
    edi_obj.Header.acqdate = mt_object.station_metadata.time_period.start_date
    edi_obj.Header.coordinate_system = (
        mt_object.station_metadata.orientation.reference_frame
    )
    edi_obj.Header.dataid = mt_object.station
    edi_obj.Header.declination = mt_object.station_metadata.location.declination.value
    edi_obj.Header.elev = mt_object.elevation
    edi_obj.Header.fileby = mt_object.station_metadata.provenance.submitter.author
    edi_obj.Header.filedate = mt_object.station_metadata.provenance.creation_time
    edi_obj.Header.lat = mt_object.latitude
    edi_obj.Header.lon = mt_object.longitude
    edi_obj.Header.datum = mt_object.station_metadata.location.datum
    edi_obj.Header.stdvers = "SEG 1.0"
    edi_obj.Header.units = mt_object.station_metadata.transfer_function.units
    
    ### write notes
    # write transfer function info first
    for k, v in mt_object.station_metadata.transfer_function.to_dict(single=True).items():
        if not v in [None]:
            if k in ["processing_parameters"]:
                for item in v:
                    edi_obj.Info.info_list.append(item.replace('=', ' = '))
            else:
                edi_obj.Info.info_list.append(f"{k} = {v}")
                
    # write comments, which would be anything in the info section from an edi
    edi_obj.Info.info_list += mt_object.station_metadata.comments.split('\n')
            
    # write field notes
    for run in mt_object.station_metadata.run_list:
        write_dict = dict([(comp, False) for comp in ['ex', 'ey', 'hx', 'hy', 'hz',
                                                      'temperature', 'rrhx', 'rrhy']])
        for cc in write_dict.keys():
            if getattr(run, cc).component is not None:
                write_dict[cc] = True
                                                       
        r_dict = run.to_dict(single=True)
        
        for rk, rv in r_dict.items():
            if rv not in [None]:
                if rk[0:2] in ['ex', 'ey', 'hx', 'hy', 'hz', 'te', 'rr']:
                    if rk[0:2] == 'te':
                        comp = 'temperature'
                    elif rk[0:2] == 'rr':
                        comp = rk[0:4]
                    else:
                        comp = rk[0:2]
                    if write_dict[comp] is False:
                        continue
                    skip_list = [f"{comp}.{ff}" for ff in ['filter.name',
                                                           'filter.applied',
                                                           'time_period.start', 
                                                           'time_period.end',
                                                           'location.elevation',
                                                           'location.latitude',
                                                           'location.longitude',
                                                           'positive.latitude',
                                                           'positive.longitude',
                                                           'positive.elevation',
                                                           'negative.latitude',
                                                           'negative.longitude',
                                                           'negative.elevation']]
                    if rk not in skip_list:
                        edi_obj.Info.info_list.append(f"{run.id}.{rk} = {rv}")
                else:
                    edi_obj.Info.info_list.append(f"{run.id}.{rk} = {rv}")

    ### fill measurement
    edi_obj.Measurement.refelev = mt_object.elevation
    edi_obj.Measurement.reflat = mt_object.latitude
    edi_obj.Measurement.reflon = mt_object.longitude
    edi_obj.Measurement.maxchan = len(mt_object.station_metadata.channels_recorded)
    for comp in ['ex', 'ey', 'hx', 'hy', 'hz', 'rrhx', 'rrhy']:
        try:
            edi_obj.Measurement.from_metadata(getattr(mt_object, f"{comp}_metadata"))
        except AttributeError:
            edi_obj.logger.debug(f"Did not find information on {comp}")
                                             
    # input data section
    edi_obj.Data.data_type = mt_object.station_metadata.data_type
    edi_obj.Data.nfreq = mt_object.Z.z.shape[0]
    edi_obj.Data.sectid = mt_object.station
    edi_obj.Data.nchan = 5
    if np.all(mt_object.Tipper.tipper == 0) == True:
        edi_obj.Data.nchan = 4
        
    edi_obj.Data.maxblks = 999
    for comp in ['ex', 'ey', 'hx', 'hy', 'hz']:
        if hasattr(edi_obj.Measurement, f"meas_{comp}"):
            setattr(edi_obj.Data, comp, getattr(edi_obj.Measurement, f'meas_{comp}').acqchan)


    edi_obj.write_edi_file()
    
    return edi_obj
