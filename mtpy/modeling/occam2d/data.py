# -*- coding: utf-8 -*-
"""
Created on Tue Mar  7 19:01:14 2023

@author: jpeacock
"""
# =============================================================================
# Imports
# =============================================================================
from pathlib import Path

import numpy as np
import pandas as pd
from loguru import logger

from mtpy.core.mt_dataframe import MTDataFrame

# =============================================================================


class Data:
    """
    Reads and writes data files and more.

    Inherets Profile, so the intended use is to use Data to project stations
    onto a profile, then write the data file.

    ===================== =====================================================
    Model Modes           Description
    ===================== =====================================================
    1 or log_all          Log resistivity of TE and TM plus Tipper
    2 or log_te_tip       Log resistivity of TE plus Tipper
    3 or log_tm_tip       Log resistivity of TM plus Tipper
    4 or log_te_tm        Log resistivity of TE and TM
    5 or log_te           Log resistivity of TE
    6 or log_tm           Log resistivity of TM
    7 or all              TE, TM and Tipper
    8 or te_tip           TE plus Tipper
    9 or tm_tip           TM plus Tipper
    10 or te_tm           TE and TM mode
    11 or te              TE mode
    12 or tm              TM mode
    13 or tip             Only Tipper
    ===================== =====================================================


    **data** : is a list of dictioinaries containing the data for each station.
               keys include:
                   * 'station' -- name of station
                   * 'offset' -- profile line offset
                   * 'te_res' -- TE resisitivity in linear scale
                   * 'tm_res' -- TM resistivity in linear scale
                   * 'te_phase' -- TE phase in degrees
                   * 'tm_phase' --  TM phase in degrees in first quadrant
                   * 're_tip' -- real part of tipper along profile
                   * 'im_tip' -- imaginary part of tipper along profile

               each key is a np.ndarray(2, num_freq)
               index 0 is for data
               index 1 is for error

    ===================== =====================================================
    Key Words/Attributes  Desctription
    ===================== =====================================================
    _data_header          header line in data file
    _data_string          full data string
    _profile_generated    [ True | False ] True if profile has already been
                          generated.
    _rotate_to_strike     [ True | False ] True to rotate data to strike
                          angle.  *default* is True
    data                  list of dictionaries of data for each station.
                          see above
    data_fn               full path to data file
    data_list             list of lines to write to data file
    edi_list              list of mtpy.core.mt instances for each .edi file
                          read
    edi_path              directory path where .edi files are
    edi_type              [ 'z' | 'spectra' ] for .edi format
    elevation_model       model elevation np.ndarray(east, north, elevation)
                          in meters
    elevation_profile     elevation along profile np.ndarray (x, elev) (m)
    fn_basename           data file basename *default* is OccamDataFile.dat
    freq                  list of frequencies to use for the inversion
    freq_max              max frequency to use in inversion. *default* is None
    freq_min              min frequency to use in inversion. *default* is None
    freq_num              number of frequencies to use in inversion
    geoelectric_strike    geoelectric strike angle assuming N = 0, E = 90
    masked_data           similar to data, but any masked points are now 0
    mode_dict             dictionary of model modes to chose from
    model_mode            model mode to use for inversion, see above
    num_edi               number of stations to invert for
    occam_dict            dictionary of occam parameters to use internally
    occam_format          occam format of data file.
                          *default* is OCCAM2MTDATA_1.0
    phase_te_err          percent error in phase for TE mode. *default* is 5
    phase_tm_err          percent error in phase for TM mode. *default* is 5
    profile_angle         angle of profile line realtive to N = 0, E = 90
    profile_line          m, b coefficients for mx+b definition of profile line
    res_te_err            percent error in resistivity for TE mode.
                          *default* is 10
    res_tm_err            percent error in resistivity for TM mode.
                          *default* is 10
    error_type            'floor' or 'absolute' - option to set error as floor
                          (i.e. maximum of the data error and a specified value)
                          or a straight value
    save_path             directory to save files to
    station_list          list of station for inversion
    station_locations     station locations along profile line
    tipper_err            percent error in tipper. *default* is 5
    title                 title in data file.
    ===================== =====================================================

    =========================== ===============================================
    Methods                     Description
    =========================== ===============================================
    _fill_data                  fills the data array that is described above
    _get_data_list              gets the lines to write to data file
    _get_frequencies            gets frequency list to invert for
    get_profile_origin          get profile origin in UTM coordinates
    mask_points                 masks points in data picked from
                                plot_mask_points
    plot_mask_points            plots data responses to interactively mask
                                data points.
    plot_resonse                plots data/model responses, returns
                                PlotResponse data type.
    read_data_file              read in existing data file and fill appropriate
                                attributes.
    write_data_file             write a data file according to Data attributes
    =========================== ===============================================

    :Example Write Data File: ::
        >>> import mtpy.modeling.occam2d as occam2d
        >>> edipath = r"/home/mt/edi_files"
        >>> slst = ['mt{0:03}'.format(ss) for ss in range(1, 20)]
        >>> ocd = occam2d.Data(edi_path=edipath, station_list=slst)
        >>> # model just the tm mode and tipper
        >>> ocd.model_mode = 3
        >>> ocd.save_path = r"/home/occam/Line1/Inv1"
        >>> ocd.write_data_file()
        >>> # mask points
        >>> ocd.plot_mask_points()
        >>> ocd.mask_points()

    """

    def __init__(self, dataframe=None, center_point=None, **kwargs):

        self.logger = logger
        self.dataframe = dataframe
        self.data_fn = None
        self.fn_basename = "OccamDataFile.dat"
        self.save_path = Path()
        self.freq = None
        self.interpolate_freq = None
        self.model_mode = "1"
        self.data = None
        self.data_list = None

        self.res_te_err = 10
        self.res_tm_err = 10
        self.phase_te_err = 5
        self.phase_tm_err = 5
        self.tipper_err = 10
        self.error_type = "floor"

        self.freq_min = None
        self.freq_max = None
        self.freq_num = None
        self.freq_tol = None

        self.occam_format = "OCCAM2MTDATA_1.0"
        self.title = "MTpy-OccamDatafile"
        self.edi_type = "z"
        self.masked_data = None

        self.occam_dict = {
            "1": "log_te_res",
            "2": "te_phase",
            "3": "re_tip",
            "4": "im_tip",
            "5": "log_tm_res",
            "6": "tm_phase",
            "9": "te_res",
            "10": "tm_res",
        }

        self.mode_dict = {
            "log_all": [1, 2, 3, 4, 5, 6],
            "log_te_tip": [1, 2, 3, 4],
            "log_tm_tip": [5, 6, 3, 4],
            "log_te_tm": [1, 2, 5, 6],
            "log_te": [1, 2],
            "log_tm": [5, 6],
            "all": [9, 2, 3, 4, 10, 6],
            "te_tip": [9, 2, 3, 4],
            "tm_tip": [10, 6, 3, 4],
            "te_tm": [9, 2, 10, 6],
            "te": [9, 2],
            "tm": [10, 6],
            "tip": [3, 4],
            "1": [1, 2, 3, 4, 5, 6],
            "2": [1, 2, 3, 4],
            "3": [5, 6, 3, 4],
            "4": [1, 2, 5, 6],
            "5": [1, 2],
            "6": [5, 6],
            "7": [9, 2, 3, 4, 10, 6],
            "8": [9, 2, 3, 4],
            "9": [10, 6, 3, 4],
            "10": [9, 2, 10, 6],
            "11": [9, 2],
            "12": [10, 6],
            "13": [3, 4],
        }

        self._data_string = "{0:^6}{1:^6}{2:^6} {3: >8} {4: >8}\n"
        self._data_header = "{0:<6}{1:<6}{2:<6} {3:<8} {4:<8}\n".format(
            "SITE", "FREQ", "TYPE", "DATUM", "ERROR"
        )

        for key, value in kwargs.items():
            setattr(self, key, value)

    @property
    def dataframe(self):
        return self._mt_dataframe.dataframe

    @dataframe.setter
    def dataframe(self, df):
        """
        Set dataframe to an MTDataframe
        :param df: DESCRIPTION
        :type df: TYPE
        :return: DESCRIPTION
        :rtype: TYPE

        """

        if df is None:
            self._mt_dataframe = MTDataFrame()

        elif isinstance(df, (pd.DataFrame, MTDataFrame, np.ndarray)):
            self._mt_dataframe = MTDataFrame(df)

        else:
            raise TypeError(
                f"Input must be a dataframe or MTDataFrame object not {type(df)}"
            )

    def read_data_file(self, data_fn=None):
        """
        Read in an existing data file and populate appropriate attributes
            * data
            * data_list
            * freq
            * station_list
            * station_locations

        Arguments:
        -----------
            **data_fn** : string
                          full path to data file
                          *default* is None and set to save_path/fn_basename

        :Example: ::

            >>> import mtpy.modeling.occam2d as occam2d
            >>> ocd = occam2d.Data()
            >>> ocd.read_data_file(r"/home/Occam2D/Line1/Inv1/Data.dat")

        """

        if data_fn is not None:
            self.data_fn = data_fn

        if not self.data_fn.is_file():
            raise ValueError("Could not find {0}".format(self.data_fn))
        if self.data_fn is None:
            raise ValueError("data_fn is None, input filename")

        self.save_path = self.data_fn.parent

        print(f"Reading from {self.data_fn}")

        with open(self.data_fn, "r") as dfid:
            dlines = dfid.readlines()

        # get format of input data
        self.occam_format = dlines[0].strip().split(":")[1].strip()

        # get title
        title_str = dlines[1].strip().split(":")[1].strip()

        title_list = title_str.split(",")
        self.title = title_list[0]

        # get strike angle and profile angle
        if len(title_list) > 1:
            for t_str in title_list[1:]:
                t_list = t_str.split("=")
                if len(t_list) > 1:
                    key = t_list[0].strip().lower().replace(" ", "_")
                    if key == "profile":
                        key = "profile_angle"
                    elif key == "strike":
                        key = "geoelectric_strike"
                    value = t_list[1].split("deg")[0].strip()
                    print("    {0} = {1}".format(key, value))
                    try:
                        setattr(self, key, float(value))
                    except ValueError:
                        setattr(self, key, value)

        # get number of sites
        nsites = int(dlines[2].strip().split(":")[1].strip())
        print("    {0} = {1}".format("number of sites", nsites))

        # get station names
        self.station_list = np.array(
            [dlines[ii].strip() for ii in range(3, nsites + 3)]
        )

        # get offsets in meters
        self.station_locations = np.array(
            [
                float(dlines[ii].strip())
                for ii in range(4 + nsites, 4 + 2 * nsites)
            ]
        )

        # get number of frequencies
        nfreq = int(dlines[4 + 2 * nsites].strip().split(":")[1].strip())
        print("    {0} = {1}".format("number of frequencies", nfreq))

        # get frequencies
        self.freq = np.array(
            [
                float(dlines[ii].strip())
                for ii in range(5 + 2 * nsites, 5 + 2 * nsites + nfreq)
            ]
        )

        # get periods
        self.period = 1.0 / self.freq

        # -----------get data-------------------
        # set zero array size the first row will be the data and second the
        # error
        asize = (2, self.freq.shape[0])

        # make a list of dictionaries for each station.
        self.data = [
            {
                "station": station,
                "offset": offset,
                "te_phase": np.zeros(asize),
                "tm_phase": np.zeros(asize),
                "re_tip": np.zeros(asize),
                "im_tip": np.zeros(asize),
                "te_res": np.zeros(asize),
                "tm_res": np.zeros(asize),
            }
            for station, offset in zip(
                self.station_list, self.station_locations
            )
        ]

        self.data_list = dlines[7 + 2 * nsites + nfreq :]
        for line in self.data_list:
            try:
                station, freq, comp, odata, oerr = line.split()
                # station index -1 cause python starts at 0
                ss = int(station) - 1

                # frequency index -1 cause python starts at 0
                ff = int(freq) - 1
                # data key
                key = self.occam_dict[comp]

                # put into array
                if int(comp) == 1 or int(comp) == 5:
                    self.data[ss][key[4:]][0, ff] = 10 ** float(odata)
                    # error
                    self.data[ss][key[4:]][1, ff] = float(oerr) * np.log(10)
                else:
                    self.data[ss][key][0, ff] = float(odata)
                    # error
                    self.data[ss][key][1, ff] = float(oerr)
            except ValueError:
                print("Could not read line {0}".format(line))

    def _get_frequencies(self):
        """
        from the list of edi's get a frequency list to invert for.

        Uses Attributes:
        ------------
            **freq_min** : float (Hz)
                           minimum frequency to invert for.
                           *default* is None and will use the data to find min

            **freq_max** : float (Hz)
                           maximum frequency to invert for
                           *default* is None and will use the data to find max

            **freq_num** : int
                           number of frequencies to invert for
                           *default* is None and will use the data to find num
        """
        if self.interpolate_freq:
            if self.freq is not None:
                return

        if self.freq is not None:
            return

        # get all frequencies from all edi files
        lo_all_freqs = []
        for edi in self.edi_list:
            lo_all_freqs.extend(list(edi.Z.freq))

        # sort all frequencies so that they are in descending order,
        # use set to remove repeats and make an array
        all_freqs = np.array(sorted(list(set(lo_all_freqs)), reverse=True))

        # --> get min and max values if none are given
        if (
            (self.freq_min is None)
            or (self.freq_min < all_freqs.min())
            or (self.freq_min > all_freqs.max())
        ):
            self.freq_min = all_freqs.min()

        if (
            (self.freq_max is None)
            or (self.freq_max > all_freqs.max())
            or (self.freq_max < all_freqs.min())
        ):
            self.freq_max = all_freqs.max()

        # --> get all frequencies within the given range
        self.freq = all_freqs[
            np.where(
                (all_freqs >= self.freq_min) & (all_freqs <= self.freq_max)
            )
        ]

        if len(self.freq) == 0:
            raise ValueError(
                "No frequencies in user-defined interval "
                "[{0}, {1}]".format(self.freq_min, self.freq_max)
            )

        # check, if frequency list is longer than given max value
        if self.freq_num is not None:
            if int(self.freq_num) < self.freq.shape[0]:
                print(
                    (
                        "Number of frequencies exceeds freq_num "
                        "{0} > {1} ".format(self.freq.shape[0], self.freq_num)
                        + "Trimming frequencies to {0}".format(self.freq_num)
                    )
                )

                excess = self.freq.shape[0] / float(self.freq_num)
                if excess < 2:
                    offset = 0
                else:
                    stepsize = (self.freq.shape[0] - 1) / self.freq_num
                    offset = stepsize / 2.0
                indices = np.array(
                    np.around(
                        np.linspace(
                            offset,
                            self.freq.shape[0] - 1 - offset,
                            self.freq_num,
                        ),
                        0,
                    ),
                    dtype="int",
                )
                if indices[0] > (self.freq.shape[0] - 1 - indices[-1]):
                    indices -= 1
                self.freq = self.freq[indices]

    def _fill_data(self):
        """
        Read all Edi files.
        Create a profile
        rotate impedance and tipper
        Extract frequencies.

        Collect all information sorted according to occam specifications.

        Data of Z given in muV/m/nT = km/s
        Error is assumed to be 1 stddev.
        """

        # create a profile line, this sorts the stations by offset and rotates
        # data.
        self.generate_profile()
        self.plot_profile()

        # --> get frequencies to invert for
        self._get_frequencies()

        # set zero array size the first row will be the data and second the
        # error
        asize = (2, self.freq.shape[0])

        # make a list of dictionaries for each station.
        self.data = [
            {
                "station": station,
                "offset": offset,
                "te_phase": np.zeros(asize),
                "tm_phase": np.zeros(asize),
                "re_tip": np.zeros(asize),
                "im_tip": np.zeros(asize),
                "te_res": np.zeros(asize),
                "tm_res": np.zeros(asize),
            }
            for station, offset in zip(
                self.station_list, self.station_locations
            )
        ]

        # loop over mt object in edi_list and use a counter starting at 1
        # because that is what occam starts at.
        for s_index, edi in enumerate(self.edi_list):

            if self.interpolate_freq:
                station_freq = edi.Z.freq
                interp_freq = self.freq[
                    np.where(
                        (self.freq >= station_freq.min())
                        & (self.freq <= station_freq.max())
                    )
                ]
                # interpolate data onto given frequency list
                z_interp, t_interp = edi.interpolate(interp_freq)
                z_interp.compute_resistivity_phase()

                rho = z_interp.resistivity
                phi = z_interp.phase
                rho_err = z_interp.resistivity_err
                if t_interp is not None:
                    tipper = t_interp.tipper
                    tipper_err = t_interp.tipper_err
                else:
                    tipper = None
                    tipper_err = None
                # update station freq, as we've now interpolated new z values
                # for the station
                station_freq = self.freq[
                    np.where(
                        (self.freq >= station_freq.min())
                        & (self.freq <= station_freq.max())
                    )
                ]
            else:
                station_freq = edi.Z.freq
                rho = edi.Z.resistivity
                rho_err = edi.Z.resistivity_err
                phi = edi.Z.phase
                tipper = edi.Tipper.tipper
                tipper_err = edi.Tipper.tipper_err

            self.data[s_index]["station"] = edi.station
            self.data[s_index]["offset"] = edi.offset

            for freq_num, frequency in enumerate(self.freq):
                if self.freq_tol is not None:
                    try:
                        f_index = np.where(
                            (station_freq >= frequency * (1 - self.freq_tol))
                            & (station_freq <= frequency * (1 + self.freq_tol))
                        )[0][0]

                    except IndexError:
                        f_index = None
                else:
                    # skip, if the listed frequency is not available for the
                    # station
                    if frequency in station_freq:
                        # find the respective frequency index for the station
                        f_index = np.abs(station_freq - frequency).argmin()
                    else:
                        f_index = None

                if f_index == None:
                    continue

                # --> get te resistivity
                self.data[s_index]["te_res"][0, freq_num] = rho[f_index, 0, 1]
                # compute error
                if rho[f_index, 0, 1] != 0.0:
                    # --> get error from data
                    if (self.res_te_err is None) or (
                        self.error_type == "floor"
                    ):
                        error_val = np.abs(rho_err[f_index, 0, 1])
                        if error_val > rho[f_index, 0, 1]:
                            error_val = rho[f_index, 0, 1]
                            # set error floor if desired
                        if self.error_type == "floor":
                            error_val = max(
                                error_val,
                                rho[f_index, 0, 1] * self.res_te_err / 100.0,
                            )

                        self.data[s_index]["te_res"][1, freq_num] = error_val
                        # --> set generic error
                    else:
                        self.data[s_index]["te_res"][1, freq_num] = (
                            self.res_te * self.res_te_err / 100.0
                        )

                # --> get tm resistivity
                self.data[s_index]["tm_res"][0, freq_num] = rho[f_index, 1, 0]
                # compute error
                if rho[f_index, 1, 0] != 0.0:
                    # --> get error from data
                    if (self.res_tm_err is None) or (
                        self.error_type == "floor"
                    ):
                        error_val = np.abs(rho_err[f_index, 1, 0])
                        if error_val > rho[f_index, 1, 0]:
                            error_val = rho[f_index, 1, 0]
                        if self.error_type == "floor":
                            error_val = max(
                                error_val,
                                rho[f_index, 1, 0] * self.res_tm_err / 100.0,
                            )
                        self.data[s_index]["tm_res"][1, freq_num] = error_val
                    # --> set generic error
                    else:
                        self.data[s_index]["tm_res"][1, freq_num] = (
                            self.res_tm * self.res_tm_err / 100.0
                        )

                # --> get te phase
                # be sure the phase is positive and in the first quadrant
                phase_te = phi[f_index, 0, 1] % 180

                if (phase_te < 0) or (phase_te > 90):
                    phase_te = 0
                    self.data[s_index]["te_res"][0, freq_num] = 0

                # assign phase to array
                self.data[s_index]["te_phase"][0, freq_num] = phase_te

                # compute error
                # if phi[f_index, 0, 1] != 0.0:
                # --> get error from data
                if (self.phase_te_err is None) or (self.error_type == "floor"):
                    error_val = np.degrees(
                        np.arcsin(
                            min(
                                0.5
                                * rho_err[f_index, 0, 1]
                                / rho[f_index, 0, 1],
                                1.0,
                            )
                        )
                    )
                    if self.error_type == "floor":
                        error_val = max(
                            error_val, (self.phase_te_err / 100.0) * 57.0 / 2.0
                        )
                    self.data[s_index]["te_phase"][1, freq_num] = error_val
                # --> set generic error floor
                else:
                    self.data[s_index]["te_phase"][1, freq_num] = (
                        (self.phase_te_err / 100.0) * 57.0 / 2.0
                    )

                # --> get tm phase and be sure it's positive and in the first quadrant
                phase_tm = phi[f_index, 1, 0] % 180

                if (phase_tm < 0) or (phase_tm > 90):
                    phase_tm = 0
                    self.data[s_index]["tm_res"][0, freq_num] = 0

                # assign phase to array
                self.data[s_index]["tm_phase"][0, freq_num] = phase_tm

                # compute error
                # if phi[f_index, 1, 0] != 0.0:
                # --> get error from data
                if (self.phase_tm_err is None) or (self.error_type == "floor"):
                    error_val = np.degrees(
                        np.arcsin(
                            min(
                                0.5
                                * rho_err[f_index, 1, 0]
                                / rho[f_index, 1, 0],
                                1.0,
                            )
                        )
                    )
                    if self.error_type == "floor":
                        error_val = max(
                            error_val, (self.phase_tm_err / 100.0) * 57.0 / 2.0
                        )
                    self.data[s_index]["tm_phase"][1, freq_num] = error_val
                # --> set generic error floor
                else:
                    self.data[s_index]["tm_phase"][1, freq_num] = (
                        (self.phase_tm_err / 100.0) * 57.0 / 2.0
                    )

                # --> get Tipper
                if tipper is not None:
                    self.data[s_index]["re_tip"][0, freq_num] = tipper[
                        f_index, 0, 1
                    ].real
                    self.data[s_index]["im_tip"][0, freq_num] = tipper[
                        f_index, 0, 1
                    ].imag

                    # get error
                    if (self.tipper_err is not None) or (
                        self.error_type == "floor"
                    ):
                        error_val = self.tipper_err / 100.0
                        if self.error_type == "floor":
                            error_val = max(
                                error_val, tipper_err[f_index, 0, 1]
                            )
                        self.data[s_index]["re_tip"][1, freq_num] = error_val
                        self.data[s_index]["im_tip"][1, freq_num] = error_val
                    else:
                        self.data[s_index]["re_tip"][1, freq_num] = (
                            tipper[f_index, 0, 1].real
                            / tipper_err[f_index, 0, 1]
                        )
                        self.data[s_index]["im_tip"][1, freq_num] = (
                            tipper[f_index, 0, 1].imag
                            / tipper_err[f_index, 0, 1]
                        )

    def _get_data_list(self):
        """
        Get all the data needed to write a data file.

        """

        self.data_list = []
        for ss, sdict in enumerate(self.data, 1):
            for ff in range(self.freq.shape[0]):
                for mmode in self.mode_dict[self.model_mode]:
                    # log(te_res)
                    if mmode == 1:
                        if sdict["te_res"][0, ff] != 0.0:
                            dvalue = np.log10(sdict["te_res"][0, ff])
                            derror = (
                                sdict["te_res"][1, ff] / sdict["te_res"][0, ff]
                            ) / np.log(10.0)
                            dstr = "{0:.4f}".format(dvalue)
                            derrstr = "{0:.4f}".format(derror)
                            line = self._data_string.format(
                                ss, ff + 1, mmode, dstr, derrstr
                            )
                            self.data_list.append(line)
                    # te_res
                    if mmode == 9:
                        if sdict["te_res"][0, ff] != 0.0:
                            dvalue = sdict["te_res"][0, ff]
                            derror = sdict["te_res"][1, ff]
                            dstr = "{0:.4f}".format(dvalue)
                            derrstr = "{0:.4f}".format(derror)
                            line = self._data_string.format(
                                ss, ff + 1, mmode, dstr, derrstr
                            )
                            self.data_list.append(line)

                    # te_phase
                    if mmode == 2:
                        if sdict["te_phase"][0, ff] != 0.0:
                            dvalue = sdict["te_phase"][0, ff]
                            derror = sdict["te_phase"][1, ff]
                            dstr = "{0:.4f}".format(dvalue)
                            derrstr = "{0:.4f}".format(derror)
                            line = self._data_string.format(
                                ss, ff + 1, mmode, dstr, derrstr
                            )
                            self.data_list.append(line)
                            # log(tm_res)
                    if mmode == 5:
                        if sdict["tm_res"][0, ff] != 0.0:
                            dvalue = np.log10(sdict["tm_res"][0, ff])
                            (
                                sdict["tm_res"][1, ff] / sdict["tm_res"][0, ff]
                            ) / np.log(10)
                            dstr = "{0:.4f}".format(dvalue)
                            derrstr = "{0:.4f}".format(derror)
                            line = self._data_string.format(
                                ss, ff + 1, mmode, dstr, derrstr
                            )
                            self.data_list.append(line)

                    # tm_res
                    if mmode == 10:
                        if sdict["tm_res"][0, ff] != 0.0:
                            dvalue = sdict["tm_res"][0, ff]
                            derror = sdict["tm_res"][1, ff]
                            dstr = "{0:.4f}".format(dvalue)
                            derrstr = "{0:.4f}".format(derror)
                            line = self._data_string.format(
                                ss, ff + 1, mmode, dstr, derrstr
                            )
                            self.data_list.append(line)

                    # tm_phase
                    if mmode == 6:
                        if sdict["tm_phase"][0, ff] != 0.0:
                            dvalue = sdict["tm_phase"][0, ff]
                            derror = sdict["tm_phase"][1, ff]
                            dstr = "{0:.4f}".format(dvalue)
                            derrstr = "{0:.4f}".format(derror)
                            line = self._data_string.format(
                                ss, ff + 1, mmode, dstr, derrstr
                            )
                            self.data_list.append(line)

                    # Re_tip
                    if mmode == 3:
                        if sdict["re_tip"][0, ff] != 0.0:
                            dvalue = sdict["re_tip"][0, ff]
                            derror = sdict["re_tip"][1, ff]
                            dstr = "{0:.4f}".format(dvalue)
                            derrstr = "{0:.4f}".format(derror)
                            line = self._data_string.format(
                                ss, ff + 1, mmode, dstr, derrstr
                            )
                            self.data_list.append(line)

                    # Im_tip
                    if mmode == 4:
                        if sdict["im_tip"][0, ff] != 0.0:
                            dvalue = sdict["im_tip"][0, ff]
                            derror = sdict["im_tip"][1, ff]
                            dstr = "{0:.4f}".format(dvalue)
                            derrstr = "{0:.4f}".format(derror)
                            line = self._data_string.format(
                                ss, ff + 1, mmode, dstr, derrstr
                            )
                            self.data_list.append(line)

    def mask_from_datafile(self, mask_datafn):
        """
        reads a separate data file and applies mask from this data file.
        mask_datafn needs to have exactly the same frequencies, and station names
        must match exactly.

        """
        ocdm = Data()
        ocdm.read_data_file(mask_datafn)
        # list of stations, in order, for the mask_datafn and the input data
        # file
        ocdm_stlist = [ocdm.data[i]["station"] for i in range(len(ocdm.data))]
        ocd_stlist = [self.data[i]["station"] for i in range(len(self.data))]

        for i_ocd, stn in enumerate(ocd_stlist):
            i_ocdm = ocdm_stlist.index(stn)
            for dmode in [
                "te_res",
                "tm_res",
                "te_phase",
                "tm_phase",
                "im_tip",
                "re_tip",
            ]:

                for i in range(len(self.freq)):
                    if self.data[i_ocdm][dmode][0][i] == 0:
                        self.data[i_ocd][dmode][0][i] = 0.0
        self.fn_basename = (
            self.fn_basename[:-4] + "Masked" + self.fn_basename[-4:]
        )
        self.write_data_file()

    def write_data_file(self, data_fn=None):
        """
        Write a data file.

        Arguments:
        -----------
            **data_fn** : string
                          full path to data file.
                          *default* is save_path/fn_basename

        If there data is None, then _fill_data is called to create a profile,
        rotate data and get all the necessary data.  This way you can use
        write_data_file directly without going through the steps of projecting
        the stations, etc.

        :Example: ::
            >>> edipath = r"/home/mt/edi_files"
            >>> slst = ['mt{0:03}'.format(ss) for ss in range(1, 20)]
            >>> ocd = occam2d.Data(edi_path=edipath, station_list=slst)
            >>> ocd.save_path = r"/home/occam/line1/inv1"
            >>> ocd.write_data_file()

        """

        if self.data is None:
            self._fill_data()

        # get the appropriate data to write to file
        self._get_data_list()

        if data_fn is not None:
            self.data_fn = data_fn
        else:
            if self.save_path is None:
                self.save_path = Path()
            if not self.save_path.exists():
                self.save_path.mkdir()

            self.data_fn = self.save_path.joinpath(self.fn_basename)

        data_lines = []

        # --> header line
        data_lines.append("{0:<18}{1}\n".format("FORMAT:", self.occam_format))

        # --> title line
        if self.profile_angle is None:
            self.profile_angle = 0
        if self.geoelectric_strike is None:
            self.geoelectric_strike = 0.0
        t_str = "{0}, Profile={1:.1f} deg, Strike={2:.1f} deg".format(
            self.title, self.profile_angle, self.geoelectric_strike
        )
        data_lines.append("{0:<18}{1}\n".format("TITLE:", t_str))

        # --> sites
        data_lines.append("{0:<18}{1}\n".format("SITES:", len(self.data)))
        for sdict in self.data:
            data_lines.append("   {0}\n".format(sdict["station"]))

        # --> offsets
        data_lines.append("{0:<18}\n".format("OFFSETS (M):"))
        for sdict in self.data:
            data_lines.append("   {0:.1f}\n".format(sdict["offset"]))
        # --> frequencies
        data_lines.append(
            "{0:<18}{1}\n".format("FREQUENCIES:", self.freq.shape[0])
        )
        for ff in self.freq:
            data_lines.append("   {0:<10.6e}\n".format(ff))

        # --> data
        data_lines.append(
            "{0:<18}{1}\n".format("DATA BLOCKS:", len(self.data_list))
        )
        data_lines.append(self._data_header)
        data_lines += self.data_list

        dfid = open(self.data_fn, "w")
        dfid.writelines(data_lines)
        dfid.close()

        print("Wrote Occam2D data file to {0}".format(self.data_fn))

    def get_profile_origin(self):
        """
        get the origin of the profile in real world coordinates

        Author: Alison Kirkby (2013)

        NEED TO ADAPT THIS TO THE CURRENT SETUP.
        """

        x, y = self.easts, self.norths
        x1, y1 = x[0], y[0]
        [m, c1] = self.profile
        x0 = (y1 + (1.0 / m) * x1 - c1) / (m + (1.0 / m))
        y0 = m * x0 + c1
        self.profile_origin = [x0, y0]
