# -*- coding: utf-8 -*-
"""
.. module:: JFile
   :synopsis: Deal with J-Files of the format propsed by Alan Jones 

.. moduleauthor:: Jared Peacock <jpeacock@usgs.gov>
"""

# ==============================================================================
from pathlib import Path
import numpy as np

import mtpy.core.z as mtz
from mtpy.utils.mtpy_logger import get_mtpy_logger

from mt_metadata.transfer_functions import tf as metadata
from mt_metadata.utils.mttime import MTime

# ==============================================================================
# Class to read j_file
# ==============================================================================
class JFile(object):
    """
    be able to read and write a j-file
    """

    def __init__(self, j_fn=None):
        self.logger = get_mtpy_logger(f"{__name__}.{self.__class__.__name__}")
        self._jfn = None
        self.j_fn = j_fn
        self.header_dict = None
        self.metadata_dict = None
        self.Z = None
        self.Tipper = None
        self.station = None

        if self.j_fn is not None:
            self.read_j_file()
            
    def __str__(self):
        lines = [f"Station: {self.station}", "-" * 50]
        lines.append(f"\tSurvey:        {self.survey_metadata.survey_id}")
        lines.append(f"\tProject:       {self.survey_metadata.project}")
        lines.append(f"\tAcquired by:   {self.station_metadata.acquired_by.author}")
        lines.append(f"\tAcquired date: {self.station_metadata.time_period.start_date}")
        lines.append(f"\tLatitude:      {self.latitude:.3f}")
        lines.append(f"\tLongitude:     {self.longitude:.3f}")
        lines.append(f"\tElevation:     {self.elevation:.3f}")
        if self.Z.z is not None:
            lines.append("\tImpedance:     True")
        else:
            lines.append("\tImpedance:     False")
        if self.Tipper.tipper is not None:
            lines.append("\tTipper:        True")
        else:
            lines.append("\tTipper:        False")

        if self.Z.z is not None:
            lines.append(f"\tN Periods:     {len(self.Z.freq)}")

            lines.append("\tPeriod Range:")
            lines.append(f"\t\tMin:   {self.periods.min():.5E} s")
            lines.append(f"\t\tMax:   {self.periods.max():.5E} s")

            lines.append("\tFrequency Range:")
            lines.append(f"\t\tMin:   {self.frequencies.max():.5E} Hz")
            lines.append(f"\t\tMax:   {self.frequencies.min():.5E} Hz")

        return "\n".join(lines)

    def __repr__(self):
        lines = []
        lines.append(f"station='{self.station}'")
        lines.append(f"latitude={self.latitude:.2f}")
        lines.append(f"longitude={self.longitude:.2f}")
        lines.append(f"elevation={self.elevation:.2f}")

        return f"MT( {(', ').join(lines)} )"

    @property
    def j_fn(self):
        return self._jfn

    @j_fn.setter
    def j_fn(self, value):
        """
        set file name
        :param value: DESCRIPTION
        :type value: TYPE
        :return: DESCRIPTION
        :rtype: TYPE

        """
        self._jfn = Path(value)
        
    @property
    def latitude(self):
        return self.metadata_dict["latitude"]
    
    @latitude.setter
    def latitude(self, value):
        self.metadata_dict["latitude"] = value
        
    @property
    def longitude(self):
        return self.metadata_dict["longitude"]
    
    @longitude.setter
    def longitude(self, value):
        self.metadata_dict["longitude"] = value
        
    @property
    def elevation(self):
        return self.metadata_dict["elevation"]
    
    @elevation.setter
    def elevation(self, value):
        self.metadata_dict["elevation"] = value
        
    @property
    def periods(self):
        if not np.all(self.Z.z == 0) and self.Z is not None:
            return 1./self.Z.freq
        if not np.all(self.Tipper.tipper == 0) and self.Tipper is not None:
            return 1./self.Tipper.freq
        
    @property
    def frequencies(self):
        if not np.all(self.Z.z == 0) and self.Z is not None:
            return self.Z.freq
        if not np.all(self.Tipper.tipper == 0) and self.Tipper is not None:
            return self.Tipper.freq

    def _validate_j_file(self):
        """
        change the lat, lon, elev lines to something machine readable,
        if they are not.
        """

        if not self.j_fn.exists():
            msg = f"Could not find {self.j_fn}, check path"
            self.logger.error(msg)
            raise NameError(msg)

        with open(self.j_fn, "r", errors="replace") as fid:
            j_lines = fid.readlines()

        for variable in ["lat", "lon", "elev"]:
            for ii, line in enumerate(j_lines):
                if variable in line.lower():
                    name = line.split("=")[0]
                    try:
                        value = float(line.split("=")[1].strip())
                    except ValueError:
                        value = 0.0
                        self.logger.debug(f"Changed {name[1:]} to 0.0")
                    j_lines[ii] = "{0} = {1}\n".format(name, value)
                    break

        return j_lines

    def _read_header_line(self, line):
        """
        read a header line
        """
        line = " ".join(line[1:].strip().split())

        new_line = ""

        # need to restructure the string so its readable, at least the way
        # that birrp outputs the file
        e_find = 0
        for ii in range(len(line)):
            if line[ii] == "=":
                e_find = ii
                new_line += line[ii]
            elif line[ii] == " ":
                if abs(e_find - ii) == 1:
                    pass
                else:
                    new_line += ","
            else:
                new_line += line[ii]

        # now that we have a useful line, split it into its parts
        line_list = new_line.split(",")

        # try to split up the parts into a key=value setup
        # and try to make the values floats if they can be
        l_dict = {}
        key = "null"
        for ll in line_list:
            ll_list = ll.split("=")
            if len(ll_list) == 1:
                continue

            # some times there is just a list of numbers, need a way to read
            # that.
            if len(ll_list) != 2:
                if type(l_dict[key]) is not list:
                    l_dict[key] = list([l_dict[key]])
                try:
                    l_dict[key].append(float(ll))
                except ValueError:
                    l_dict[key].append(ll)
            else:
                key = ll_list[0]
                try:
                    value = float(ll_list[1])
                except ValueError:
                    value = ll_list[1]

                l_dict[key] = value

        return l_dict

    def read_header(self, j_lines=None):
        """
        Parsing the header lines of a j-file to extract processing information.
    
        Input:
        - j-file as list of lines (output of readlines())
    
        Output:
        - Dictionary with all parameters found

        """
        if j_lines is None:
            j_lines = self._validate_j_file()
        header_lines = [j_line for j_line in j_lines if "#" in j_line]
        header_dict = {"title": header_lines[0][1:].strip()}

        fn_count = 0
        theta_count = 0
        # put the information into a dictionary
        for h_line in header_lines[1:]:
            h_dict = self._read_header_line(h_line)
            for key in list(h_dict.keys()):
                if key == "filnam":
                    h_key = "{0}_{1:02}".format(key, fn_count)
                    fn_count += 1
                elif key == "nskip" or key == "nread":
                    h_key = "{0}_{1:02}".format(key, fn_count - 1)

                # if its the line of angles, put them all in a list with a unique key
                elif key == "theta1":
                    h_key = "{0}_{1:02}".format(key, theta_count)
                    theta_count += 1

                elif key == "theta2" or key == "phi":
                    h_key = "{0}_{1:02}".format(key, theta_count - 1)
                else:
                    h_key = key

                header_dict[h_key] = h_dict[key]

        self.header_dict = header_dict

    def read_metadata(self, j_lines=None, j_fn=None):
        """
        read in the metadata of the station, or information of station 
        logistics like: lat, lon, elevation
        
        Not really needed for a birrp output since all values are nan's
        """
        if j_lines is None:
            j_lines = self._validate_j_file()

        metadata_lines = [j_line for j_line in j_lines if ">" in j_line]

        metadata_dict = {}
        for m_line in metadata_lines:
            m_list = m_line.strip().split("=")
            m_key = m_list[0][1:].strip().lower()
            try:
                m_value = float(m_list[0].strip())
            except ValueError:
                m_value = 0.0

            metadata_dict[m_key] = m_value

        self.metadata_dict = metadata_dict

    def read_j_file(self, j_fn=None):
        """
        read_j_file will read in a *.j file output by BIRRP (better than reading lots of *.<k>r<l>.rf files)
    
        Input:
        j-filename
    
        Output: 4-tuple
        - periods : N-array
        - Z_array : 2-tuple - values and errors
        - tipper_array : 2-tuple - values and errors
        - processing_dict : parsed processing parameters from j-file header
    
        """
        # read data
        z_index_dict = {"zxx": (0, 0), "zxy": (0, 1), "zyx": (1, 0), "zyy": (1, 1)}
        t_index_dict = {"tzx": (0, 0), "tzy": (0, 1)}

        if j_fn is not None:
            self.j_fn = j_fn

        self.logger.debug(f"Reading {self.j_fn}")

        j_line_list = self._validate_j_file()

        self.read_header(j_lines=j_line_list)
        self.read_metadata(j_lines=j_line_list)

        data_lines = [
            j_line for j_line in j_line_list if not ">" in j_line and not "#" in j_line
        ][:]

        self.station = data_lines[0].strip()

        # sometimes birrp outputs some missing periods, so the best way to deal with
        # this that I could come up with was to get things into dictionaries with
        # key words that are the period values, then fill in Z and T from there
        # leaving any missing values as 0

        # make empty dictionary that have keys as the component
        z_dict = dict([(z_key, {}) for z_key in list(z_index_dict.keys())])
        t_dict = dict([(t_key, {}) for t_key in list(t_index_dict.keys())])
        for d_line in data_lines[1:]:
            # check to see if we are at the beginning of a component block, if so
            # set the dictionary key to that value
            if "z" in d_line.lower():
                d_key = d_line.strip().split()[0].lower()
            # if we are at the number of periods line, skip it
            elif len(d_line.strip().split()) == 1 and "r" not in d_line.lower():
                continue
            elif "r" in d_line.lower():
                break
            # get the numbers into the correct dictionary with a key as period and
            # for now we will leave the numbers as a list, which we will parse later
            else:
                # split the line up into each number
                d_list = d_line.strip().split()

                # make a copy of the list to be sure we don't rewrite any values,
                # not sure if this is necessary at the moment
                d_value_list = list(d_list)
                for d_index, d_value in enumerate(d_list):
                    # check to see if the column number can be converted into a float
                    # if it can't, then it will be set to 0, which is assumed to be
                    # a masked number when writing to an .edi file

                    try:
                        d_value = float(d_value)
                        # need to check for masked points represented by
                        # birrp as -999, apparently
                        if d_value == -999 or np.isnan(d_value):
                            d_value_list[d_index] = 0.0
                        else:
                            d_value_list[d_index] = d_value
                    except ValueError:
                        d_value_list[d_index] = 0.0

                # put the numbers in the correct dictionary as:
                # key = period, value = [real, imaginary, error]
                if d_key in list(z_index_dict.keys()):
                    z_dict[d_key][d_value_list[0]] = d_value_list[1:4]
                elif d_key in list(t_index_dict.keys()):
                    t_dict[d_key][d_value_list[0]] = d_value_list[1:4]

        # --> now we need to get the set of periods for all components
        # check to see if there is any tipper data output

        all_periods = []
        for z_key in list(z_index_dict.keys()):
            for f_key in list(z_dict[z_key].keys()):
                all_periods.append(f_key)

        if len(list(t_dict["tzx"].keys())) == 0:
            self.logger.info(f"Could not find any Tipper data in {self.j_fn}")
            find_tipper = False

        else:
            for t_key in list(t_index_dict.keys()):
                for f_key in list(t_dict[t_key].keys()):
                    all_periods.append(f_key)
            find_tipper = True

        all_periods = np.array(sorted(list(set(all_periods))))
        all_periods = all_periods[np.nonzero(all_periods)]
        num_per = len(all_periods)

        # fill arrays using the period key from all_periods
        z_arr = np.zeros((num_per, 2, 2), dtype=np.complex)
        z_err_arr = np.zeros((num_per, 2, 2), dtype=np.float)

        t_arr = np.zeros((num_per, 1, 2), dtype=np.complex)
        t_err_arr = np.zeros((num_per, 1, 2), dtype=np.float)

        for p_index, per in enumerate(all_periods):
            for z_key in sorted(z_index_dict.keys()):
                kk = z_index_dict[z_key][0]
                ll = z_index_dict[z_key][1]
                try:
                    z_value = z_dict[z_key][per][0] + 1j * z_dict[z_key][per][1]
                    z_arr[p_index, kk, ll] = z_value
                    z_err_arr[p_index, kk, ll] = z_dict[z_key][per][2]
                except KeyError:
                    self.logger.debug(f"No value found for period {per:.4g}")
                    self.logger.debug(f"For component {z_key}")
            if find_tipper is True:
                for t_key in sorted(t_index_dict.keys()):
                    kk = t_index_dict[t_key][0]
                    ll = t_index_dict[t_key][1]
                    try:
                        t_value = t_dict[t_key][per][0] + 1j * t_dict[t_key][per][1]
                        t_arr[p_index, kk, ll] = t_value
                        t_err_arr[p_index, kk, ll] = t_dict[t_key][per][2]
                    except KeyError:
                        self.logger.debug(f"No value found for period {per:.4g}")
                        self.logger.debug(f"For component {t_key}")

        # put the results into mtpy objects
        freq = 1.0 / all_periods
        z_arr[np.where(z_arr == np.inf)] = 0 + 0j
        t_arr[np.where(t_arr == np.inf)] = 0 + 0j
        z_err_arr[np.where(z_err_arr == np.inf)] = 10 ** 6
        t_err_arr[np.where(t_err_arr == np.inf)] = 10 ** 6

        self.Z = mtz.Z(z_arr, z_err_arr, freq)
        self.Tipper = mtz.Tipper(t_arr, t_err_arr, freq)

    @property
    def station_metadata(self):
        sm = metadata.Station()
        r1 = metadata.Run(id=f"{self.station}a")

        if not np.all(self.Z.z == 0):
            r1.ex = metadata.Electric(component="ex", channel_id=1)
            r1.ey = metadata.Electric(component="ey", channel_id=2)
            r1.hx = metadata.Magnetic(component="hx", channel_id=3)
            r1.hy = metadata.Magnetic(component="hy", channel_id=4)

        if not np.all(self.Tipper.tipper == 0):
            r1.hz = metadata.Magnetic(component="hz", channel_id=5)

        sm.run_list.append(r1)
        sm.id = self.station
        sm.data_type = "MT"

        sm.location.latitude = self.metadata_dict["latitude"]
        sm.location.longitude = self.metadata_dict["longitude"]
        sm.location.elevation = self.metadata_dict["elevation"]

        # provenance
        sm.provenance.software.name = "BIRRP"
        sm.provenance.software.version = "5"
        sm.transfer_function.processed_date = MTime(self.j_fn.stat().st_ctime).iso_str
        sm.transfer_function.runs_processed = sm.run_names
        # add birrp parameters
        for key, value in self.header_dict.items():
            sm.transfer_function.processing_parameters.append(f"{key} = {value}")

        return sm

    @property
    def survey_metadata(self):
        sm = metadata.Survey()

        return sm


def read_jfile(j_fn):
    """
    Read a .j file output by BIRRP
    
    :param j_fn: full path to j file
    :type j_fn: string or :class:`pathlib.Path`
    
    """

    from mtpy.core import mt

    j_obj = JFile(j_fn)

    mt_obj = mt.MT()
    mt_obj._fn = j_fn

    for attr in [
        "Z",
        "Tipper",
        "survey_metadata",
        "station_metadata",
    ]:
        setattr(mt_obj, attr, getattr(j_obj, attr))

    # need to set latitude to compute UTM coordinates to make sure station
    # location is estimated for ModEM
    mt_obj.latitude = j_obj.station_metadata.location.latitude

    return mt_obj


def write_jfile(mt_obj, fn=None):
    """
    
    :param mt_obj: DESCRIPTION
    :type mt_obj: TYPE
    :param fn: DESCRIPTION, defaults to None
    :type fn: TYPE, optional
    :return: DESCRIPTION
    :rtype: TYPE

    """

    raise IOError("write_jfile not implemented yet.")
