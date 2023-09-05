"""
==================
ModEM
==================

# Generate files for ModEM

# revised by JP 2017
# revised by AK 2017 to bring across functionality from ak branch
# revised by JP 2021 adding functionality and updating.
# revised by JP 2022 to work with new structure of a central object

"""
# =============================================================================
# Imports
# =============================================================================
import numpy as np
from pathlib import Path
import pandas as pd
from loguru import logger

from mtpy.core.mt_dataframe import MTDataFrame
from mtpy.core.mt_location import MTLocation
from mtpy.modeling.errors import ModelErrors

# =============================================================================
class Data:
    """
    Data will read and write .dat files for ModEM and convert a WS data file
    to ModEM format.

    ..note: :: the data is interpolated onto the given periods such that all
               stations invert for the same periods.  The interpolation is
               a linear interpolation of each of the real and imaginary parts
               of the impedance tensor and induction tensor.
               See mtpy.core.mt.MT.interpolate for more details
    
    :param edi_list: list of edi files to read

    ====================== ====================================================
    Attributes              Description
    ====================== ====================================================
    _dtype                 internal variable defining the data type of
                           data_array
    _logger                python logging object that put messages in logging
                           format defined in logging configure file, see MtPyLog
                           for more information
    _t_shape               internal variable defining shape of tipper array in
                           _dtype
    _z_shape               internal variable defining shape of Z array in
                           _dtype
    center_position        (east, north, evel) for center point of station
                           array.  All stations are relative to this location
                           for plotting purposes.
    comp_index_dict        dictionary for index values of component of Z and T
    station_locations      Stations object
    data_array             numpy.ndarray (num_stations) structured to store
                           data.  keys are:
                               * station --> station name
                               * lat --> latitude in decimal degrees
                               * lon --> longitude in decimal degrees
                               * elev --> elevation (m)
                               * rel_east -- > relative east location to
                                               center_position (m)
                               * rel_north --> relative north location to
                                               center_position (m)
                               * east --> UTM east (m)
                               * north --> UTM north (m)
                               * zone --> UTM zone
                               * z --> impedance tensor array with shape
                                       (num_freq, 2, 2)
                               * z_err --> impedance tensor error array with
                                       shape (num_freq, 2, 2)
                               * tip --> Tipper array with shape
                                       (num_freq, 1, 2)
                               * tipperr --> Tipper array with shape
                                       (num_freq, 1, 2)
    data_fn                full path to data file
    data_period_list       period list from all the data
    edi_list               list of full paths to edi files
    error_type_tipper      [ 'abs' | 'floor' ]
                           *default* is 'abs'
    error_type_z           [ 'egbert' | 'mean_od' | 'eigen' | 'median']
                           *default* is 'egbert_floor'
                                * add '_floor' to any of the above to set the
                                  error as an error floor, otherwise all
                                  components are give weighted the same

                                * 'egbert'  sets error to
                                            error_value_z * sqrt(abs(zxy*zyx))
                                * 'mean_od' sets error to
                                            error_value_z * mean([Zxy, Zyx])
                                            (non zeros)
                                * 'eigen'   sets error to
                                            error_value_z * eigenvalues(Z[ii])
                                * 'median'  sets error to
                                            error_value_z * median([Zxx, Zxy, Zyx, Zyy])
                                            (non zeros)
                           A 2x2 numpy array of error_type_z can be specified to
                           explicitly set the error_type_z for each component.

    error_value_z          percentage to multiply Z by to set error
                           *default* is 5 for 5% of Z as error
                           A 2x2 numpy array of values can be specified to
                           explicitly set the error_value_z for each component.

    error_value_tipper     absolute error between 0 and 1.
    fn_basename            basename of data file. *default* is 'ModEM_Data.dat'
    formatting             ['1' | '2'], format of the output data file, *default* is '1'
    header_strings         strings for header of data file following the format
                           outlined in the ModEM documentation
    inv_comp_dict          dictionary of inversion components
    inv_mode               inversion mode, options are: *default* is '1'
                               * '1' --> for 'Full_Impedance' and
                                             'Full_Vertical_Components'
                               * '2' --> 'Full_Impedance'
                               * '3' --> 'Off_Diagonal_Impedance' and
                                         'Full_Vertical_Components'
                               * '4' --> 'Off_Diagonal_Impedance'
                               * '5' --> 'Full_Vertical_Components'
                               * '6' --> 'Full_Interstation_TF'
                               * '7' --> 'Off_Diagonal_Rho_Phase'

    inv_mode_dict          dictionary for inversion modes
    max_num_periods        maximum number of periods
    model_epsg             epsg code for model projection, provide this to
                           project model to non-utm coordinates. Find the epsg
                           code for your projection on
                           http://spatialreference.org/ref/ or google search
                           epsg "your projection"
    model_utm_zone         alternative to model_epsg, choose a utm zone to
                           project all sites to (e.g. '55S')
    mt_dict                dictionary of mtpy.core.mt.MT objects with keys
                           being station names
    period_buffer          float or int
                           if specified, apply a buffer so that interpolation doesn't
                           stretch too far over periods
    period_dict            dictionary of period index for period_list
    period_list            list of periods to invert for
    period_max             maximum value of period to invert for
    period_min             minimum value of period to invert for
    period_buffer          buffer so that interpolation doesn't stretch too far
                              over periods. Provide a float or integer factor, 
                              greater than which interpolation will not stretch.
                              e.g. 1.5 means only interpolate to a maximum of
                              1.5 times each side of each frequency value
    rotate_angle           Angle to rotate data to assuming 0 is N and E is 90
    save_path              path to save data file to
    units                  [ [V/m]/[T] | [mV/km]/[nT] | Ohm ] units of Z
                           *default* is [mV/km]/[nT]
    wave_sign_impedance    [ + | - ] sign of time dependent wave.
                           *default* is '+' as positive downwards.
    wave_sign_tipper       [ + | - ] sign of time dependent wave.
                           *default* is '+' as positive downwards.
    ====================== ====================================================


    :Example 1 --> create inversion period list: ::

        >>> from pathlib import Path
        >>> import mtpy.modeling.modem as modem
        >>> edi_path = Path(r"/home/mt/edi_files")
        >>> edi_list = list(edi_path.glob("*.edi"))
        >>> md = modem.Data(edi_list, period_min=.1, period_max=300,\
        >>> ...             max_num_periods=12)
        >>> md.write_data_file(save_path=r"/home/modem/inv1")
        >>> md
        

    :Example 2 --> set inverions period list from data: ::

        >>> md = modem.Data(edi_list)
        >>> #get period list from an .edi file
        >>> inv_period_list = 1./md.mt_dict["mt01"].Z.freq
        >>> #invert for every third period in inv_period_list
        >>> inv_period_list = inv_period_list[np.arange(0, len(inv_period_list, 3))]
        >>> md.period_list = inv_period_list
        >>> md.write_data_file(save_path=r"/home/modem/inv1")

    :Example 3 --> change error values: ::

        >>> mdr.error_type = 'floor'
        >>> mdr.error_floor = 10
        >>> mdr.error_tipper = .03
        >>> mdr.write_data_file(save_path=r"/home/modem/inv2")

    :Example 4 --> change inversion type: ::

        >>> mdr.inv_mode = '3'
        >>> mdr.write_data_file(save_path=r"/home/modem/inv2")

    :Example 5 --> rotate data: ::

        >>> md.rotation_angle = 60
        >>> md.write_data_file(save_path=r"/home/modem/Inv1")
        >>> # or
        >>> md.write_data_file(save_path=r"/home/modem/Inv1", \
                               rotation_angle=60)


    """

    def __init__(self, dataframe=None, center_point=None, **kwargs):

        self.logger = logger

        self.dataframe = dataframe

        if center_point is None:
            self.center_point = MTLocation()
        else:
            self.center_point = center_point

        self.wave_sign_impedance = "+"
        self.wave_sign_tipper = "+"
        self.z_units = "[mV/km]/[nT]"
        self.t_units = ""
        self.inv_mode = "1"
        self.formatting = "1"
        self.rotation_angle = 0

        self.z_model_error = ModelErrors(
            error_value=5,
            error_type="geometric_mean",
            floor=True,
            mode="impedance",
        )
        self.t_model_error = ModelErrors(
            error_value=0.02,
            error_type="absolute",
            floor=True,
            mode="tipper",
        )

        self.fn_basename = "ModEM_Data.dat"
        self.save_path = Path.cwd()

        self.topography = True

        self.inv_mode_dict = {
            "1": ["Full_Impedance", "Full_Vertical_Components"],
            "2": ["Full_Impedance"],
            "3": ["Off_Diagonal_Impedance", "Full_Vertical_Components"],
            "4": ["Off_Diagonal_Impedance"],
            "5": ["Full_Vertical_Components"],
            "6": ["Full_Interstation_TF"],
            "7": ["Off_Diagonal_Rho_Phase"],
        }
        self.inv_comp_dict = {
            "Full_Impedance": ["zxx", "zxy", "zyx", "zyy"],
            "Off_Diagonal_Impedance": ["zxy", "zyx"],
            "Full_Vertical_Components": ["tzx", "tzy"],
        }

        self.header_string = " ".join(
            [
                "# Period(s)",
                "Code",
                "GG_Lat",
                "GG_Lon",
                "X(m)",
                "Y(m)",
                "Z(m)",
                "Component",
                "Real",
                "Imag",
                "Error",
            ]
        )

        self._df_keys = [
            "period",
            "station",
            "latitude",
            "longitude",
            "model_north",
            "model_east",
            "model_elevation",
            "comp",
            "real",
            "imag",
            "error",
        ]

        for key, value in kwargs.items():
            setattr(self, key, value)

    def __str__(self):
        lines = ["ModEM Data Object:"]
        if self.dataframe is not None:
            lines += [
                f"\tNumber of impedance stations: {self.get_n_stations('impedance')}"
            ]
            lines += [
                f"\tNumber of tipper stations: {self.get_n_stations('vertical')}"
            ]
            lines += [
                f"\tNumber of phase tensor stations: {self.get_n_stations('phase_tensor')}"
            ]
            lines += [f"\tNumber of periods:  {self.n_periods}"]
            lines += ["\tPeriod range (s):  "]
            lines += [f"\t\tMin: {self.period.min():.5g}"]
            lines += [f"\t\tMax: {self.period.max():.5g}"]
            lines += [f"\tRotation angle:     {self.rotation_angle}"]
            lines += ["\tData center:        "]
            lines += [
                f"\t\tLatitude:  {self.center_point.latitude:>8.4f} deg "
                f"\tNorthing: {self.center_point.north:.4f} m"
            ]
            lines += [
                f"\t\tLongitude: {self.center_point.longitude:>8.4f} deg "
                f"\tEasting: {self.center_point.east:.4f} m"
            ]
            lines += [
                f"\t\tDatum epsg: {self.center_point.datum_epsg}"
                f"\t\t\tUTM epsg:   {self.center_point.utm_epsg}"
            ]
            lines += [f"\t\tElevation:  {self.center_point.elevation:.1f} m"]

            lines += [
                f"\tImpedance data:     {self.dataframe.zxy.mean() != 0.0}"
            ]
            lines += [
                f"\tTipper data:        {self.dataframe.tzx.mean() != 0.0}"
            ]
            lines += [
                f"\tInversion Mode:   {', '.join(self.inv_mode_dict[self.inv_mode])}"
            ]

        return "\n".join(lines)

    def __repr__(self):
        return self.__str__()

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

        self._mt_dataframe.dataframe.reset_index(drop=True, inplace=True)

    @property
    def model_parameters(self):
        params = {
            "wave_sign_impedance": self.wave_sign_impedance,
            "wave_sign_tipper": self.wave_sign_tipper,
            "z_units": self.z_units,
            "t_units": self.t_units,
            "inv_mode": self.inv_mode,
            "formatting": self.formatting,
            "data_filename": self.data_filename,
            "topography": self.topography,
            "rotation_angle": self.rotation_angle,
            "center_point.latitude": self.center_point.latitude,
            "center_point.longitue": self.center_point.longitude,
            "center_point.elevation": self.center_point.elevation,
            "center_point.utm_epsg": self.center_point.utm_epsg,
            "center_point.datum_epsg": self.center_point.datum_epsg,
        }

        for key, value in self.z_model_error.error_parameters.items():
            params[f"z_model_error.{key}"] = value
        for key, value in self.t_model_error.error_parameters.items():
            params[f"t_model_error.{key}"] = value

        return params

    @property
    def data_filename(self):
        return self.save_path.joinpath(self.fn_basename)

    @data_filename.setter
    def data_filename(self, value):
        if value is not None:
            value = Path(value)
            if value.parent == Path("."):
                self.fn_basename = value.name
            else:
                self.save_path = value.parent
                self.fn_basename = value.name

    @property
    def period(self):
        if self.dataframe is not None:
            return np.sort(self.dataframe.period.unique())

    def get_n_stations(self, mode):
        if self.dataframe is not None:
            if "impedance" in mode.lower():
                return (
                    self.dataframe.loc[
                        (self.dataframe.zxx != 0)
                        | (self.dataframe.zxy != 0)
                        | (self.dataframe.zyx != 0)
                        | (self.dataframe.zyy != 0),
                        "station",
                    ]
                    .unique()
                    .size
                )
            elif "vertical" in mode.lower():
                return (
                    self.dataframe.loc[
                        (self.dataframe.tzx != 0) | (self.dataframe.tzy != 0),
                        "station",
                    ]
                    .unique()
                    .size
                )
            elif "phase_tensor" in mode.lower():
                return (
                    self.dataframe.loc[
                        (self.dataframe.ptxx != 0)
                        | (self.dataframe.ptxy != 0)
                        | (self.dataframe.ptyx != 0)
                        | (self.dataframe.ptyy != 0),
                        "station",
                    ]
                    .unique()
                    .size
                )

    @property
    def n_periods(self):
        return self.period.size

    def _get_components(self):
        """
        get components to write out
        """

        comps = []
        for inv_modes in self.inv_mode_dict[self.inv_mode]:
            comps += self.inv_comp_dict[inv_modes]

        return comps

    def _get_header_string(self, error_type, error_value):
        """
        Create the header strings

        # Created using MTpy calculated egbert_floor error of 5% data rotated 0.0_deg
        clockwise from N

        :param error_type: The method to calculate the errors
        :type error_type: string
        :param error_value: value of error or error floor
        :type error_value: float
        :param rotation_angle: angle data have been rotated by
        :type rotation_angle: float

        """

        h_str = []
        if np.atleast_1d(error_type).ndim == 2:
            h_str = (
                f"# Creating_software: MTpy v2, "
                f"error: [{error_type[0, 0]}, {error_type[0, 1]}, "
                f"{error_type[1, 0]}, {error_type[1, 1]}], "
            )
        else:
            h_str = f"# Creating_software: MTpy v2, error: {error_type}, "

        if np.atleast_1d(error_value).ndim == 2:
            h_str += (
                f"error floors of {error_value[0, 0]:.0f}%, "
                f"{error_value[0, 1]:.0f}%, "
                f"{error_value[1, 0]:.0f}%, "
                f"{error_value[1, 1]:.0f}%, "
                f"data rotated {self.rotation_angle:.1f}_deg clockwise from N, "
                f"{self.center_point.utm_crs}"
            )

        else:
            if error_value > 1:
                fmt = ".0f"
                units = "%"
            elif error_value < 1:
                fmt = ".2f"
                units = ""
            h_str += (
                f"error_value: {error_value:{fmt}}{units}, data_rotation: "
                f"{self.rotation_angle:.1f} deg clockwise, "
                f"model_{self.center_point.utm_crs}"
            )

        return h_str

    def _write_header(self, mode):
        """ """
        d_lines = []
        if "impedance" in mode.lower():
            d_lines.append(
                self._get_header_string(
                    self.z_model_error.error_type,
                    self.z_model_error.error_value,
                )
            )
            d_lines.append(self.header_string)
            d_lines.append(f"> {mode}")
            d_lines.append(f"> exp({self.wave_sign_impedance}i\omega t)")
            d_lines.append(f"> {self.z_units}")

        elif "vertical" in mode.lower():
            d_lines.append(
                self._get_header_string(
                    self.t_model_error.error_type,
                    self.t_model_error.error_value,
                )
            )
            d_lines.append(self.header_string)
            d_lines.append(f"> {mode}")
            d_lines.append(f"> exp({self.wave_sign_tipper}i\omega t)")
            d_lines.append(f"> [{self.t_units}]")

        d_lines.append(
            f"> {self.rotation_angle:.3g}"
        )  # orientation, need to add at some point
        if self.topography:
            d_lines.append(
                f"> {self.center_point.latitude:>10.6f} "
                f"{self.center_point.longitude:>10.6f} "
                f"{self.center_point.model_elevation:>10.2f}"
            )
        else:
            d_lines.append(
                f"> {self.center_point.latitude:>10.6f} "
                f"{self.center_point.longitude:>10.6f}"
            )

        n_stations = self.get_n_stations(mode)
        d_lines.append(f"> {self.n_periods} {n_stations}")

        return d_lines

    def _write_comp(self, row, comp):
        """
        write a single row

        :param row: DESCRIPTION
        :type row: TYPE
        :param comp: DESCRIPTION
        :type comp: TYPE
        :return: DESCRIPTION
        :rtype: TYPE

        """

        value = np.nan_to_num(getattr(row, comp))
        err = getattr(row, f"{comp}_model_error")

        if (
            value.real != 0.0
            and value.imag != 0.0
            and value.real != 1e32
            and value.imag != 1e32
        ):
            if self.formatting == "1":
                per = f"{row.period:<12.5e}"
                sta = f"{row.station:>7}"
                lat = f"{row.latitude:> 9.3f}"
                lon = f"{row.longitude:> 9.3f}"
                eas = f"{row.model_east:> 12.3f}"
                nor = f"{row.model_north:> 12.3f}"
                if self.topography:
                    ele = f"{row.model_elevation:> 12.3f}"
                else:
                    ele = f"{0:> 12.3f}"
                if comp[1].lower() == "z":
                    comp = comp.replace("z", "")
                com = f"{comp:>4}".upper()
                if self.z_units.lower() == "ohm":
                    rea = f"{value.real / 796.:> 14.6e}"
                    ima = f"{value.imag / 796.:> 14.6e}"
                elif self.z_units.lower() not in (
                    "[v/m]/[t]",
                    "[mv/km]/[nt]",
                ):
                    raise ValueError(f"Unsupported unit '{self.z_units}'")
                else:
                    rea = f"{value.real:> 14.6e}"
                    ima = f"{value.imag:> 14.6e}"

            elif self.formatting == "2":
                per = f"{row.period:<14.6e}"
                sta = f"{row.station:>10}"
                lat = f"{row.latitude:> 14.6f}"
                lon = f"{row.longitude:> 14.6f}"
                eas = f"{row.model_east:> 15.3f}"
                nor = f"{row.model_north:> 15.3f}"
                if self.topography:
                    ele = f"{row.model_elevation:> 10.3f}"
                else:
                    ele = f"{0:> 10.3f}"
                if comp[1].lower() == "z":
                    comp = comp.replace("z", "")
                com = f"{comp:>4}".upper()
                if self.z_units.lower() == "ohm":
                    rea = f"{value.real / 796.:> 17.6e}"
                    ima = f"{value.imag / 796.:> 17.6e}"
                elif self.z_units.lower() not in (
                    "[v/m]/[t]",
                    "[mv/km]/[nt]",
                ):
                    raise ValueError(f"Unsupported unit '{self.z_units}'")
                else:
                    rea = f"{value.real:> 17.6e}"
                    ima = f"{value.imag:> 17.6e}"

            else:
                raise NotImplementedError(
                    f"format {self.formatting} ({type(self.formatting)}) is "
                    "not supported."
                )

            if np.isinf(err) or np.isnan(err):
                err = 10 ** (
                    np.floor(np.log10(abs(max([float(rea), float(ima)]))))
                )
            abs_err = f"{err:> 14.6e}"

            return "".join(
                [
                    per,
                    sta,
                    lat,
                    lon,
                    nor,
                    eas,
                    ele,
                    com,
                    rea,
                    ima,
                    abs_err,
                ]
            )

    def _check_for_errors_of_zero(self):
        """
        Need to check for any zeros in the error values which can prevent
        ModEM from running.

        :return: DESCRIPTION
        :rtype: TYPE

        """

        ## check for zeros in model error
        for comp in ["zxx", "zxy", "zyx", "zyy", "tzx", "tzy"]:
            find_zeros = np.where(self.dataframe[f"{comp}_model_error"] == 0)[
                0
            ]
            if find_zeros.shape[0] > 0:
                if comp in ["zxx", "zxy", "zyx", "zyy"]:
                    error_percent = self.z_model_error.error_value
                elif "t" in comp:
                    error_percent = self.t_model_error.error_value

                self.logger.warning(
                    f"Found errors with values of 0 in {comp} "
                    f"{len(find_zeros)} times. Setting error as {comp} x "
                    f"{error_percent}."
                )

                self.dataframe.loc[
                    find_zeros.tolist(), f"{comp}_model_error"
                ] = (
                    abs(self.dataframe[f"{comp}"].iloc[list(find_zeros)])
                    * error_percent
                )

    def _check_for_too_small_errors(self, tol=0.02):
        """
        Check for too small of errors relative to the error floor
        """

        for comp in ["zxx", "zxy", "zyx", "zyy", "tzx", "tzy"]:
            find_small = np.where(
                self.dataframe[f"{comp}_model_error"]
                / abs(self.dataframe[comp])
                < tol
            )[0]
            if find_small.shape[0] > 0:

                if "z" in comp:
                    error_percent = self.z_model_error.error_value
                elif "t" in comp:
                    error_percent = self.t_model_error.error_value

                self.logger.warning(
                    f"Found errors with values less than {tol} in {comp} "
                    f"{len(find_small)} times. Setting error as {comp} x "
                    f"{error_percent}."
                )
                self.dataframe.loc[
                    find_small.tolist(), f"{comp}_model_error"
                ] = (
                    abs(self.dataframe[f"{comp}"].iloc[list(find_small)])
                    * error_percent
                )

    def write_data_file(
        self,
        file_name=None,
        save_path=None,
        fn_basename=None,
        elevation=False,
    ):
        """
        
        :param save_path: full directory to save file to, defaults to None
        :type save_path: string or Path, optional
        :param fn_basename: Basename of the saved file, defaults to None
        :type fn_basename: string, optional
        :param elevation: If True adds in elevation from 'rel_elev' column in data
         array, defaults to False
        :type elevation: boolean, optional

        :raises NotImplementedError: If the inversion mode is not supported
        :raises ValueError: :class:`mtpy.utils.exceptions.ValueError` if a parameter
        is missing
        :return: full path to data file
        :rtype: Path

        .. code-block::
            :linenos:

            >>> from pathlib import Path
            >>> import mtpy.modeling.modem as modem
            >>> edi_path = Path(r"/home/mt/edi_files")
            >>> edi_list = list(edi_path.glob("*.ed"))
            >>> md = modem.Data(edi_list, period_min=.1, period_max=300,\
            >>> ...             max_num_periods=12)
            >>> md.write_data_file(save_path=r"/home/modem/inv1")
            /home/modem/inv1/ModemDataFile.dat
            
        """

        if self.dataframe is None:
            raise ValueError(
                "A DataFrame needs to be present to write a ModEM data file"
            )

        if file_name is not None:
            self.data_filename = file_name

        if save_path is not None:
            self.save_path = Path(save_path)
        if fn_basename is not None:
            self.data_filename = fn_basename

        self._check_for_errors_of_zero()
        self._check_for_too_small_errors()

        for inv_mode in self.inv_mode_dict[self.inv_mode]:
            if "impedance" in inv_mode.lower():
                z_lines = self._write_header(inv_mode)

            elif "vertical" in inv_mode.lower():
                t_lines = self._write_header(inv_mode)

            else:
                # maybe error here
                raise NotImplementedError(
                    f"inv_mode {inv_mode} is not supported yet"
                )

        comps = self._get_components()
        # Iterate over stations and sort by period
        for station in self.dataframe.station.unique():
            sdf = self.dataframe.loc[self.dataframe.station == station]
            sdf.sort_values("period")

            for row in sdf.itertuples():
                for comp in comps:
                    d_line = self._write_comp(row, comp)
                    if d_line is None:
                        continue

                    if comp.startswith("z"):
                        z_lines.append(d_line)
                    elif comp.startswith("t"):
                        t_lines.append(d_line)

        with open(self.data_filename, "w") as dfid:
            dfid.write("\n".join(z_lines + t_lines))

        self.logger.info(
            "Wrote ModEM data file to {0}".format(self.data_filename)
        )
        return self.data_filename

    def _read_header(self, header_lines):
        """
        Read header lines

        :param header_lines: DESCRIPTION
        :type header_lines: TYPE
        :return: DESCRIPTION
        :rtype: TYPE

        """

        mode = None
        inv_list = []
        header_list = []
        metadata_list = []
        n_periods = 0
        n_stations = 0
        self.center_point = MTLocation()
        for hline in header_lines:
            if hline.find("#") == 0:
                if "period" not in hline.lower():
                    header_list.append(hline.strip())
            elif hline.find(">") == 0:
                # modem outputs only 7 characters for the lat and lon
                # if there is a negative they merge together, need to split
                # them up
                hline = hline.replace("-", " -")
                metadata_list.append(hline[1:].strip())
                if hline.lower().find("ohm") > 0:
                    self.z_units = "ohm"
                    continue
                elif hline.lower().find("mv") > 0:
                    self.z_units = "[mV/km]/[nT]"
                    continue
                elif hline.lower().find("vertical") > 0:
                    mode = "vertical"
                    inv_list.append("Full_Vertical_Components")
                    continue
                elif hline.lower().find("impedance") > 0:
                    mode = "impedance"
                    inv_list.append("Full_Impedance")
                    continue

                if hline.find("exp") > 0:
                    if mode in ["impedance"]:
                        self.wave_sign_impedance = hline[hline.find("(") + 1]
                    elif mode in ["vertical"]:
                        self.wave_sign_tipper = hline[hline.find("(") + 1]

                elif (
                    len(hline[1:].strip().split()) >= 2
                    and hline.count(".") > 0
                ):
                    value_list = [
                        float(value) for value in hline[1:].strip().split()
                    ]
                    if value_list[0] != 0.0:
                        self.center_point.latitude = value_list[0]
                    if value_list[1] != 0.0:
                        self.center_point.longitude = value_list[1]
                    try:
                        self.center_point.elevation = value_list[2]
                    except IndexError:
                        self.center_point.elevation = 0.0
                        self.logger.debug(
                            "Did not find center elevation in data file"
                        )
                elif len(hline[1:].strip().split()) < 2:
                    try:
                        self.rotation_angle = float(hline[1:].strip())
                    except ValueError:
                        continue
                elif len(hline[1:].strip().split()) == 2:
                    n_periods = int(hline[1:].strip().split()[0])
                    n_stations = int(hline[1:].strip().split()[1])

        for head_line, inv_mode in zip(header_list, inv_list):
            self._parse_header_line(head_line, inv_mode)

        self._get_inv_mode(inv_list)

        return n_periods, n_stations

    def _parse_header_line(self, header_line, mode):
        """
        Parse header line

        """

        if header_line == self.header_string:
            return

        item_dict = {
            "error": "error_type",
            "error_value": "error_value",
            "data_rotation": "rotation_angle",
            "model_epsg": "center_point.utm_epsg",
        }

        if header_line.count(",") > 0:
            header_list = header_line.split(",")
        else:
            header_list = header_line.split()

        if "impedance" in mode.lower():
            obj = self.z_model_error

        elif "vertical" in mode.lower():
            obj = self.t_model_error

        for ii, item in enumerate(header_list):
            item = item.lower()
            if item.count(":") > 0:
                item_list = [k.strip() for k in item.split(":")]
                if len(item_list) == 2:
                    key = item_list[0]
                    value = item_list[1].replace("%", "").split()[0]
                    if key in ["error_value", "data_rotation"]:
                        try:
                            value = float(value)
                        except ValueError:
                            pass
                    try:
                        if key in ["model_epsg"]:
                            setattr(self.center_point, "utm_epsg", value)
                        elif "error" in key:

                            setattr(
                                obj,
                                item_dict[key],
                                value,
                            )
                        else:
                            setattr(self, item_dict["key"], value)
                    except KeyError:
                        continue

            ## Older files
            else:
                if item in ["calculated"]:
                    value = header_list[ii + 1]

                    if "floor" in value:
                        setattr(obj, "floor", True)
                        value = value.replace("_floor", "")
                    setattr(obj, "error_type", value)

                if item in ["of"]:
                    value = float(header_list[ii + 1].replace("%", ""))
                    setattr(obj, item_dict["error_value"], value)

                if "deg" in item:
                    setattr(
                        self,
                        item_dict["data_rotation"],
                        float(item.split("_")[0]),
                    )

    def _get_rotation_angle(self, header_line):
        # try to find rotation angle
        h_list = header_line.split()
        for hh, h_str in enumerate(h_list):
            if h_str.find("_deg") > 0:
                try:
                    self.rotation_angle = float(h_str[0 : h_str.find("_deg")])
                except ValueError:
                    pass

    def _get_inv_mode(self, inv_list):
        # find inversion mode
        for inv_key in list(self.inv_mode_dict.keys()):
            inv_mode_list = self.inv_mode_dict[inv_key]
            if len(inv_mode_list) != inv_list:
                continue
            else:
                tf_arr = np.zeros(len(inv_list), dtype=np.bool)

                for tf, data_inv in enumerate(inv_list):
                    if data_inv in self.inv_mode_dict[inv_key]:
                        tf_arr[tf] = True

                if np.alltrue(tf_arr):
                    self.inv_mode = inv_key
                    break

    def _read_line(self, line):
        """
        read a single line
        :param line: DESCRIPTION
        :type line: TYPE
        :return: DESCRIPTION
        :rtype: TYPE

        .. note:: Pandas Groupby does not play nice with complex numbers so
         we will be keeping the real and imaginary part separate for now.

        """

        line_dict = dict(
            [(key, value) for key, value in zip(self._df_keys, line.split())]
        )
        for key in [
            "period",
            "latitude",
            "longitude",
            "model_east",
            "model_north",
            "model_elevation",
            "real",
            "imag",
            "error",
        ]:
            line_dict[key] = float(line_dict[key])

        comp = line_dict.pop("comp").lower()
        if comp.startswith("t"):
            comp = comp.replace("t", "tz")
        line_dict[f"{comp}_real"] = line_dict.pop("real")
        line_dict[f"{comp}_imag"] = line_dict.pop("imag")
        line_dict[f"{comp}_model_error"] = line_dict.pop("error")
        if line_dict[f"{comp}_model_error"] > 1e10:
            line_dict[f"{comp}_model_error"] = np.nan

        return line_dict

    def read_data_file(self, data_fn):
        """

        :param data_fn: full path to data file name
        :type data_fn: string or Path
        :raises ValueError: If cannot compute component

        Fills attributes:
            * data_array
            * period_list
            * mt_dict

        .. code-block::

            >>> md = Data()
            >>> md.read_data_file(r"/home/modem_data.dat")
            >>> md
            ModEM Data Object:
                Number of stations: 169
                Number of periods:  22
                Period range:
                        Min: 0.01 s
                        Max: 15230.2 s
                Rotation angle:     0.0
                Data center:
                         latitude:  39.6351 deg
                         longitude: -119.8039 deg
                         Elevation: 0.0 m
                         Easting:   259368.9746 m
                         Northing:  4391021.1981 m
                         UTM zone:  11S
                Model EPSG:         None
                Model UTM zone:     None
                Impedance data:     True
                Tipper data:        True


        """

        self.data_filename = Path(data_fn)

        if self.data_filename is None:
            raise ValueError("data_fn is None, enter a data file to read.")
        elif not self.data_filename.is_file():
            raise ValueError(
                "Could not find {0}, check path".format(self.data_filename)
            )

        self.center_point = MTLocation()

        # open file get lines
        with open(self.data_filename, "r") as dfid:
            dlines = dfid.readlines()

        # read header information
        n_periods, n_stations = self._read_header(
            [line for line in dlines if ">" in line or "#" in line]
        )

        # create a list of dictionaries to make into a pandas dataframe
        entries = []
        for dline in dlines:
            if "#" in dline or ">" in dline:
                continue

            elif len(dline.split()) == len(self._df_keys):
                line_dict = self._read_line(dline)
                entries.append(line_dict)

        full_df = pd.DataFrame(entries)

        # group by period and station so that there is 1 row per period per station
        combined_df = full_df.groupby(
            ["station", "period"], as_index=False
        ).first()

        # combine real and imaginary
        cols = [c.split("_")[0] for c in combined_df.columns if "real" in c]
        for col in cols:
            combined_df[col] = (
                combined_df[f"{col}_real"] + 1j * combined_df[f"{col}_imag"]
            )
            combined_df.drop(
                [f"{col}_real", f"{col}_imag"], axis=1, inplace=True
            )

        return MTDataFrame(combined_df)

    def fix_data_file(self, fn=None, n=3):
        """
        A newer compiled version of Modem outputs slightly different headers
        This aims to convert that into the older format

        :param fn: DESCRIPTION, defaults to None
        :type fn: TYPE, optional
        :param n: DESCRIPTION, defaults to 3
        :type n: TYPE, optional
        :return: DESCRIPTION
        :rtype: TYPE

        """
        if fn:
            self.data_filename = Path(fn)
        with self.data_filename.open() as fid:
            lines = fid.readlines()

        def fix_line(line_list):
            return (
                " ".join("".join(line_list).replace("\n", "").split()) + "\n"
            )

        h1 = fix_line(lines[0:n])
        h2 = fix_line(lines[n : 2 * n])

        find = None
        for index, line in enumerate(lines[2 * n + 1 :], start=2 * n + 1):
            if line.find("#") >= 0:
                find = index
                break

        if find is not None:
            h3 = fix_line(lines[find : find + n])
            h4 = fix_line(lines[find + n : find + 2 * n])

            new_lines = (
                [h1, h2]
                + lines[2 * n : find]
                + [h3, h4]
                + lines[find + 2 * n :]
            )
        else:
            new_lines = [h1, h2] + lines[2 * n :]

        with self.data_filename.open("w") as fid:
            fid.writelines(new_lines)

        return self.data_filename
