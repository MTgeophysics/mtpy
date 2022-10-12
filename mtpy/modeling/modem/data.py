"""
==================
ModEM
==================

# Generate files for ModEM

# revised by JP 2017
# revised by AK 2017 to bring across functionality from ak branch
# revised by JP 2021 adding functionality and updating.

"""
# =============================================================================
# Imports
# =============================================================================
import numpy as np
from pathlib import Path

import matplotlib.pyplot as plt

from mtpy.core import mt as mt
from mtpy.core import z as mtz

# from mtpy.modeling import ws3dinv as ws
from mtpy.utils.mtpy_logger import get_mtpy_logger

# =============================================================================
class Data(object):
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

    def __init__(self, df, **kwargs):

        self.logger = get_mtpy_logger(f"{__name__}.{self.__class__.__name__}")

        self.dataframe = df

        self.wave_sign_impedance = "+"
        self.wave_sign_tipper = "+"
        self.z_units = "[mV/km]/[nT]"
        self.t_units = ""
        self.inv_mode = "1"
        self.formatting = "1"

        self.error_type_z = "geometric_mean"
        self.error_value_z = 5
        self.rotation_angle = 0
        self.error_type_tipper = "absolute"
        self.error_value_tipper = 0.02

        self.data_fn = "ModEM_Data.dat"
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

    def __str__(self):
        lines = ["ModEM Data Object:"]
        if self.dataframe is not None:
            lines += [f"\tNumber of stations: {self.n_stations}"]
            lines += [f"\tNumber of periods:  {self.n_periods}"]
            lines += ["\tPeriod range:  "]
            lines += [f"\t\tMin: {self.period.min()} s"]
            lines += [f"\t\tMax: {self.period.max()} s"]
            # lines += [f"\tRotation angle:     {self.rotation_angle}"]
            # lines += ["\tData center:        "]
            # lines += [f"\t\t latitude:  {self.center_point.lat[0]:.4f} deg"]
            # lines += [f"\t\t longitude: {self.center_point.lon[0]:.4f} deg"]
            # lines += [f"\t\t Elevation: {self.center_point.elev[0]:.1f} m"]
            # lines += [f"\t\t Easting:   {self.center_point.east[0]:.4f} m"]
            # lines += [f"\t\t Northing:  {self.center_point.north[0]:.4f} m"]
            # lines += [f"\t\t UTM zone:  {self.center_point.zone[0]}"]
            # lines += [f"\tModel EPSG:         {self.model_epsg}"]
            # lines += [f"\tModel UTM zone:     {self.model_utm_zone}"]
            # lines += [
            #     f"\tImpedance data:     {self.data_array['z'].mean() != 0.0}"
            # ]
            # lines += [
            #     f"\tTipper data:        {self.data_array['tip'].mean() != 0.0}"
            # ]
        return "\n".join(lines)

    def __repr__(self):
        return self.__str__()

    @property
    def period(self):
        if self.dataframe is not None:
            return np.sort(self.dataframe.period.unique())

    @property
    def n_stations(self):
        if self.dataframe is not None:
            return self.dataframe.station.unique().size

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

    @staticmethod
    def get_header_string(error_type, error_value, rotation_angle):
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
                f"# Created using MTpy v2 calculated "
                f"{error_type[0, 0]}, {error_type[0, 1]}, "
                f"{error_type[1, 0]}, {error_type[1, 1]} "
            )
        else:
            h_str = f"# Created using MTpy v2 calculated {error_type} "

        if np.atleast_1d(error_value).ndim == 2:
            h_str += (
                f"error floors of {error_value[0, 0]:.0f}%, "
                f"{error_value[0, 1]:.0f}%, "
                f"{error_value[1, 0]:.0f}%, "
                f"{error_value[1, 1]:.0f}%, "
                f"data rotated {rotation_angle:.1f}_deg clockwise from N"
            )

        else:
            h_str += (
                f"error of {error_value:.0f}% data rotated "
                f"{rotation_angle:.1f}_deg clockwise from N"
            )

        return h_str

    def _write_header(self, mode, center_point):
        """ """
        d_lines = []
        if "impedance" in mode.lower():
            d_lines.append(
                self.get_header_string(
                    self.error_type_z,
                    self.error_value_z,
                    self.rotation_angle,
                )
            )
            d_lines.append(self.header_string)
            d_lines.append(f"> {mode}")
            d_lines.append(f"> exp({self.wave_sign_impedance}i\omega t)")
            d_lines.append(f"> {self.z_units}")
        elif "vertical" in mode.lower():
            d_lines.append(
                self.get_header_string(
                    self.error_type_tipper,
                    self.error_value_tipper,
                    self.rotation_angle,
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
                f"> {center_point.latitude:>10.6f} "
                f"{center_point.longitude:>10.6f} "
                f"{center_point.model_elevation:>10.2f}"
            )
        else:
            d_lines.append(
                f"> {center_point.latitude:>10.6f} "
                f"{center_point.longitude:>10.6f}"
            )
        d_lines.append(f"> {self.n_periods} {self.n_stations}")

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

        value = getattr(row, comp)
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
                com = f"{comp:>4}"
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
                com = f"{comp:>4}"
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

    def write_data_file(
        self,
        center_point,
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
        if save_path is not None:
            self.save_path = Path(save_path)
        if fn_basename is not None:
            self.data_fn = fn_basename

        self.data_fn = Path(self.save_path, self.data_fn)

        for inv_mode in self.inv_mode_dict[self.inv_mode]:
            if "impedance" in inv_mode.lower():
                z_lines = self._write_header(inv_mode, center_point)

            elif "vertical" in inv_mode.lower():
                t_lines = self._write_header(inv_mode, center_point)

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

                    if comp.startswith("z"):
                        z_lines.append(d_line)
                    elif comp.startswith("t"):
                        t_lines.append(d_line)
        with open(self.data_fn, "w") as dfid:
            dfid.write("\n".join(z_lines + t_lines))

        self.logger.info("Wrote ModEM data file to {0}".format(self.data_fn))
        return self.data_fn

    # def convert_ws3dinv_data_file(
    #     self,
    #     ws_data_fn,
    #     station_fn=None,
    #     save_path=None,
    #     fn_basename="ws_data_file.dat",
    # ):
    #     """
    #     convert a ws3dinv data file into ModEM format

    #     Arguments:
    #     ------------
    #         **ws_data_fn** : string
    #                          full path to WS data file

    #         **station_fn** : string
    #                          full path to station info file output by
    #                          mtpy.modeling.ws3dinv. Or you can create one using
    #                          mtpy.modeling.ws3dinv.WSStation

    #         **save_path** : string
    #                         directory path to save data file to.
    #                         *default* is cwd

    #         **fn_basename** : string
    #                           basename to save data file as
    #                           *default* is 'ModEM_Data.dat'

    #     Outputs:
    #     -----------
    #         **data_fn** : string
    #                       full path to created data file

    #     :Example: ::

    #         >>> import mtpy.modeling.modem as modem
    #         >>> mdr = modem.Data()
    #         >>> mdr.convert_ws3dinv_data_file(r"/home/ws3dinv/inv1/WSData.dat",
    #                 station_fn=r"/home/ws3dinv/inv1/WS_Station_Locations.txt")
    #     """

    #     if not Path(ws_data_fn).is_file():
    #         raise ValueError(
    #             "Did not find {0}, check path".format(ws_data_fn)
    #         )

    #     if save_path is not None:
    #         save_path = Path(save_path)
    #     else:
    #         save_path = self.save_path

    #     # --> get data from data file
    #     wsd = ws.WSData()
    #     wsd.read_data_file(ws_data_fn, station_fn=station_fn)

    #     ns = wsd.data["station"].shape[0]
    #     nf = wsd.period_list.shape[0]

    #     self.period_list = wsd.period_list.copy()
    #     self.data_array = np.zeros(
    #         ns, dtype=self.make_dtype((nf, 2, 2), (nf, 1, 2))
    #     )

    #     # --> fill data array
    #     for ii, d_arr in enumerate(wsd.data):
    #         self.data_array[ii]["station"] = d_arr["station"]
    #         self.data_array[ii]["rel_east"] = d_arr["east"]
    #         self.data_array[ii]["rel_north"] = d_arr["north"]
    #         self.data_array[ii]["z"][:] = d_arr["z_data"]
    #         self.data_array[ii]["z_err"][:] = (
    #             d_arr["z_data_err"].real * d_arr["z_err_map"].real
    #         )
    #         self.data_array[ii]["station"] = d_arr["station"]
    #         self.data_array[ii]["lat"] = 0.0
    #         self.data_array[ii]["lon"] = 0.0
    #         self.data_array[ii]["rel_east"] = d_arr["east"]
    #         self.data_array[ii]["rel_north"] = d_arr["north"]
    #         self.data_array[ii]["elev"] = 0.0

    #     # need to change the inversion mode to be the same as the ws_data file
    #     if self.data_array["z"].all() == 0.0:
    #         if self.data_array["tip"].all() == 0.0:
    #             self.inv_mode = "4"
    #         else:
    #             self.inv_mode = "3"
    #     else:
    #         if self.data_array["tip"].all() == 0.0:
    #             self.inv_mode = "2"
    #         else:
    #             self.inv_mode = "1"

    #     # -->write file
    #     self.write_data_file()

    # def convert_modem_to_ws(
    #     self, data_fn=None, ws_data_fn=None, error_map=[1, 1, 1, 1]
    # ):
    #     """
    #     convert a ModEM data file to WS format.

    #     Arguments
    #     -------------
    #         **data_fn** : string
    #                      full path to modem data file to convert

    #         **ws_data_fn** : string
    #                          full path to write ws format data file

    #         **error_map** : [zxx, zxy, zyx, zyy] floats
    #                         error map that ws uses, weights for each component
    #                         *default* is [1, 1, 1, 1] for equal weighting
    #     Returns
    #     ------------
    #         **ws_data_fn** : string
    #                          full path of ws data file

    #         **ws_station_fn** : string
    #                             full path to ws station file

    #     Example
    #     -----------
    #         :Convert ModEM data file to WS: ::
    #             >>> import mtpy.modeling.modem as modem
    #             >>> md = modem.Data()
    #             >>> md.convert_modem_to_ws(data_fn=r"/home/mt/modem/data.dat")
    #     """

    #     if self.data_fn is not None:
    #         self.read_data_file(data_fn)

    #     if ws_data_fn is None:
    #         save_path = self.save_path
    #         ws_data_fn = Path(save_path, "WS_Data.dat")

    #     else:
    #         save_path = Path(ws_data_fn).parent

    #     station_info = ws.WSStation()
    #     station_info.east = self.data_array["rel_east"]
    #     station_info.north = self.data_array["rel_north"]
    #     station_info.names = self.data_array["station"]
    #     station_info.elev = self.data_array["elev"]
    #     station_info.save_path = save_path
    #     station_info.write_station_file()

    #     ws_data = ws.WSData()
    #     ws_data.period_list = self.period_list.copy()
    #     ws_data.z_err_map = error_map
    #     ws_data.z_err = "data"
    #     z_shape = (self.period_list.size, 2, 2)
    #     data_dtype = [
    #         ("station", "|U50"),
    #         ("east", np.float),
    #         ("north", np.float),
    #         ("z_data", (np.complex, z_shape)),
    #         ("z_data_err", (np.complex, z_shape)),
    #         ("z_err_map", (np.complex, z_shape)),
    #     ]
    #     ws_data.data = np.zeros(
    #         self.data_array["station"].size, dtype=data_dtype
    #     )
    #     ws_data.data["station"][:] = self.data_array["station"]
    #     ws_data.data["east"] = self.data_array["rel_east"]
    #     ws_data.data["north"] = self.data_array["rel_north"]
    #     ws_data.data["z_data"][:, :, :] = self.data_array["z"]
    #     ws_data.data["z_data_err"][:, :, :] = self.data_array["z_err"] * (
    #         1 + 1j
    #     )
    #     ws_data.data["z_err_map"][:, :, :] = np.array([[1, 1], [1, 1]])

    #     ws_data.write_data_file(save_path=save_path, data_fn=ws_data_fn)

    #     return ws_data.data_fn, station_info.station_fn

    def read_data_file(self, data_fn, center_utm=None):
        """

        :param data_fn: full path to data file name
        :type data_fn: string or Path
        :param center_utm: option to provide real world coordinates of the center of
        the grid for putting the data and model back into utm/grid coordinates,
        format [east_0, north_0, z_0], defaults to None
        :type center_utm: list or tuple, optional
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

        self.data_fn = Path(data_fn)
        self.save_path = self.data_fn.parent
        self.fn_basename = self.data_fn.name

        if self.data_fn is None:
            raise ValueError("data_fn is None, enter a data file to read.")
        elif not self.data_fn.is_file():
            raise ValueError(
                "Could not find {0}, check path".format(self.data_fn)
            )

        with open(self.data_fn, "r") as dfid:
            dlines = dfid.readlines()

        header_list = []
        metadata_list = []
        data_list = []
        period_list = []
        station_list = []
        read_impedance = False
        read_tipper = False
        inv_list = []
        for dline in dlines:
            if dline.find("#") == 0:
                header_list.append(dline.strip())
            elif dline.find(">") == 0:
                # modem outputs only 7 characters for the lat and lon
                # if there is a negative they merge together, need to split
                # them up
                dline = dline.replace("-", " -")
                metadata_list.append(dline[1:].strip())
                if dline.lower().find("ohm") > 0:
                    self.units = "ohm"
                    continue
                elif dline.lower().find("mv") > 0:
                    self.units = "[mV/km]/[nT]"
                    continue
                elif dline.lower().find("vertical") > 0:
                    read_tipper = True
                    read_impedance = False
                    inv_list.append("Full_Vertical_Components")
                    continue
                elif dline.lower().find("impedance") > 0:
                    read_impedance = True
                    read_tipper = False
                    inv_list.append("Full_Impedance")
                    continue
                if dline.find("exp") > 0:
                    if read_impedance is True:
                        self.wave_sign_impedance = dline[dline.find("(") + 1]
                    elif read_tipper is True:
                        self.wave_sign_tipper = dline[dline.find("(") + 1]
                elif len(dline[1:].strip().split()) >= 2:
                    if dline.find(".") > 0:
                        value_list = [
                            float(value) for value in dline[1:].strip().split()
                        ]
                        if value_list[0] != 0.0:
                            self._center_lat = value_list[0]
                        if value_list[1] != 0.0:
                            self._center_lon = value_list[1]
                        try:
                            self._center_elev = value_list[2]
                        except IndexError:
                            self._center_elev = 0.0
                            self.logger.debug(
                                "Did not find center elevation in data file"
                            )
                elif len(dline[1:].strip().split()) < 2:
                    try:
                        self.rotation_angle = float(dline[1:].strip())
                    except ValueError:
                        continue

            else:
                dline_list = dline.strip().split()
                if len(dline_list) == 11:
                    for ii, d_str in enumerate(dline_list):
                        if ii != 1:
                            try:
                                dline_list[ii] = float(d_str.strip())
                            except ValueError:
                                pass
                        # be sure the station name is a string
                        else:
                            dline_list[ii] = d_str.strip()
                    period_list.append(dline_list[0])
                    station_list.append(dline_list[1])

                    data_list.append(dline_list)

        # try to find rotation angle
        h_list = header_list[0].split()
        for hh, h_str in enumerate(h_list):
            if h_str.find("_deg") > 0:
                try:
                    self._rotation_angle = float(h_str[0 : h_str.find("_deg")])
                except ValueError:
                    pass

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

        self.period_list = np.array(sorted(set(period_list)))
        station_list = sorted(set(station_list))

        # make a period dictionary to with key as period and value as index
        period_dict = dict(
            [(per, ii) for ii, per in enumerate(self.period_list)]
        )

        data_dict = {}
        z_dummy = np.ones((len(self.period_list), 2, 2), dtype="complex")
        t_dummy = np.ones((len(self.period_list), 1, 2), dtype="complex")

        index_dict = {
            "zxx": (0, 0),
            "zxy": (0, 1),
            "zyx": (1, 0),
            "zyy": (1, 1),
            "tx": (0, 0),
            "ty": (0, 1),
        }

        # dictionary for true false if station data (lat, lon, elev, etc)
        # has been filled already so we don't rewrite it each time
        tf_dict = {}
        for station in station_list:
            data_dict[station] = mt.MT()
            data_dict[station].period = self.period_list
            data_dict[station].impedance = z_dummy.copy()
            data_dict[station].impedance_error = z_dummy.copy().real
            data_dict[station].tipper = t_dummy.copy()
            data_dict[station].tipper_error = t_dummy.copy().real

            # make sure that the station data starts out with false to fill
            # the data later
            tf_dict[station] = False

        # fill in the data for each station
        for dd in data_list:
            # get the period index from the data line
            p_index = period_dict[dd[0]]
            # get the component index from the data line
            ii, jj = index_dict[dd[7].lower()]

            # if the station data has not been filled yet, fill it
            if not tf_dict[dd[1]]:
                data_dict[dd[1]].latitude = dd[2]
                data_dict[dd[1]].longitude = dd[3]
                data_dict[dd[1]].grid_north = dd[4]
                data_dict[dd[1]].grid_east = dd[5]
                data_dict[dd[1]].grid_elev = dd[6]
                data_dict[dd[1]].elevation = dd[6]
                data_dict[dd[1]].station = dd[1]
                tf_dict[dd[1]] = True
            # fill in the impedance tensor with appropriate values
            if dd[7].find("Z") == 0:
                z_err = dd[10]
                if self.wave_sign_impedance == "+":
                    z_value = dd[8] + 1j * dd[9]
                elif self.wave_sign_impedance == "-":
                    z_value = dd[8] - 1j * dd[9]
                else:
                    raise ValueError(
                        'Incorrect wave sign "{}" (impedance)'.format(
                            self.wave_sign_impedance
                        )
                    )

                if self.units.lower() == "ohm":
                    z_value *= 796.0
                    z_err *= 796.0
                elif self.units.lower() not in ("[v/m]/[t]", "[mv/km]/[nt]"):
                    raise ValueError('Unsupported unit "{}"'.format(self.units))

                data_dict[dd[1]].impedance[p_index, ii, jj] = z_value
                data_dict[dd[1]].impedance_error[p_index, ii, jj] = z_err
            # fill in tipper with appropriate values
            elif dd[7].find("T") == 0:
                if self.wave_sign_tipper == "+":
                    data_dict[dd[1]].tipper[p_index, ii, jj] = (
                        dd[8] + 1j * dd[9]
                    )
                elif self.wave_sign_tipper == "-":
                    data_dict[dd[1]].tipper[p_index, ii, jj] = (
                        dd[8] - 1j * dd[9]
                    )
                else:
                    raise ValueError(
                        'Incorrect wave sign "{}" (tipper)'.format(
                            self.wave_sign_tipper
                        )
                    )
                data_dict[dd[1]].tipper_error[p_index, ii, jj] = dd[10]

        # make mt_dict an attribute for easier manipulation later
        self.mt_dict = data_dict

        ns = len(list(self.mt_dict.keys()))
        nf = len(self.period_list)
        self.data_array = np.zeros(
            ns, dtype=self.make_dtype((nf, 2, 2), (nf, 1, 2))
        )

        # Be sure to caclulate invariants and phase tensor for each station
        for ii, s_key in enumerate(sorted(self.mt_dict.keys())):
            mt_obj = self.mt_dict[s_key]

            # self.mt_dict[s_key].zinv.compute_invariants()
            self.mt_dict[s_key].pt.set_z_object(mt_obj.Z)
            self.mt_dict[s_key].Tipper.compute_amp_phase()
            self.mt_dict[s_key].Tipper.compute_mag_direction()

            self.data_array[ii]["station"] = mt_obj.station
            self.data_array[ii]["lat"] = mt_obj.latitude
            self.data_array[ii]["lon"] = mt_obj.longitude
            self.data_array[ii]["east"] = mt_obj.east
            self.data_array[ii]["north"] = mt_obj.north
            self.data_array[ii]["zone"] = mt_obj.utm_zone
            self.data_array[ii]["elev"] = mt_obj.elevation
            self.data_array[ii]["rel_elev"] = mt_obj.grid_elev
            self.data_array[ii]["rel_east"] = mt_obj.grid_east
            self.data_array[ii]["rel_north"] = mt_obj.grid_north

            self.data_array[ii]["z"][:] = mt_obj.Z.z
            self.data_array[ii]["z_err"][:] = mt_obj.Z.z_err
            self.data_array[ii]["z_inv_err"][:] = mt_obj.Z.z_err

            self.data_array[ii]["tip"][:] = mt_obj.Tipper.tipper
            self.data_array[ii]["tip_err"][:] = mt_obj.Tipper.tipper_err
            self.data_array[ii]["tip_inv_err"][:] = mt_obj.Tipper.tipper_err

        # option to provide real world coordinates in eastings/northings
        # (ModEM data file contains real world center in lat/lon but projection
        # is not provided so utm is assumed, causing errors when points cross
        # utm zones. And lat/lon cut off to 3 d.p. causing errors in smaller
        # areas)
        if center_utm is not None:
            self.data_array["east"] = (
                self.data_array["rel_east"] + center_utm[0]
            )
            self.data_array["north"] = (
                self.data_array["rel_north"] + center_utm[1]
            )

        if np.all(self.data_array["rel_elev"] == 0):
            self.topography = False

    def get_parameters(self):
        """
        get important parameters for documentation
        """

        parameter_list = [
            "error_type_z",
            "error_value_z",
            "error_type_tipper",
            "error_value_tipper",
            "wave_sign_impedance",
            "wave_sign_tipper",
            "rotation_angle",
            "save_path",
        ]

        parameter_dict = {}
        for parameter in parameter_list:
            key = "data.{0}".format(parameter)
            parameter_dict[key] = getattr(self, parameter)

        parameter_dict["data.period_min"] = self.period_list.min()
        parameter_dict["data.period_max"] = self.period_list.max()
        parameter_dict["data.period_num"] = self.period_list.size

        parameter_dict["data.inv_mode"] = self.inv_mode_dict[self.inv_mode]
        parameter_dict[
            "data.num_stations"
        ] = self.station_locations.station.size
        parameter_dict["data.center_point_ll"] = (
            self.center_point.lat[0],
            self.center_point.lon[0],
        )

        parameter_dict["data.center_point_utm"] = (
            self.center_point.north[0],
            self.center_point.east[0],
            self.center_point.zone[0],
        )
        return parameter_dict

    def estimate_starting_rho(self):
        """
        Estimate starting resistivity from the data.
        Creates a plot of the mean and median apparent resistivity values.

        :return: array of the median rho per period
        :rtype: np.ndarray(n_periods)
        :return: array of the mean rho per period
        :rtype: np.ndarray(n_periods)

        >>> d = Data()
        >>> d.read_data_file(r"example/data.dat")
        >>> rho_median, rho_mean = d.estimate_starting_rho()

        """
        rho = np.zeros((self.data_array.shape[0], self.period_list.shape[0]))

        for ii, d_arr in enumerate(self.data_array):
            z_obj = mtz.Z(d_arr["z"], freq=1.0 / self.period_list)
            rho[ii, :] = z_obj.res_det

        mean_rho = np.apply_along_axis(
            lambda x: x[np.nonzero(x)].mean(), 0, rho
        )
        median_rho = np.apply_along_axis(
            lambda x: np.median(x[np.nonzero(x)]), 0, rho
        )

        fig = plt.figure()

        ax = fig.add_subplot(1, 1, 1)
        (l1,) = ax.loglog(
            self.period_list, mean_rho, lw=2, color=(0.75, 0.25, 0)
        )
        (l2,) = ax.loglog(
            self.period_list, median_rho, lw=2, color=(0, 0.25, 0.75)
        )

        ax.loglog(
            self.period_list,
            np.repeat(mean_rho.mean(), self.period_list.size),
            ls="--",
            lw=2,
            color=(0.75, 0.25, 0),
        )
        ax.loglog(
            self.period_list,
            np.repeat(np.median(median_rho), self.period_list.size),
            ls="--",
            lw=2,
            color=(0, 0.25, 0.75),
        )

        ax.set_xlabel("Period (s)", fontdict={"size": 12, "weight": "bold"})
        ax.set_ylabel(
            "Resistivity (Ohm-m)", fontdict={"size": 12, "weight": "bold"}
        )

        ax.legend(
            [l1, l2],
            [
                "Mean = {0:.1f}".format(mean_rho.mean()),
                "Median = {0:.1f}".format(np.median(median_rho)),
            ],
            loc="upper left",
        )
        ax.grid(which="both", ls="--", color=(0.75, 0.75, 0.75))

        plt.show()

        return median_rho, mean_rho

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
            self.data_fn = Path(fn)
        with self.data_fn.open() as fid:
            lines = fid.readlines()

        def fix_line(line_list):
            return " ".join("".join(line_list).replace("\n", "").split()) + "\n"

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

        with self.data_fn.open("w") as fid:
            fid.writelines(new_lines)

        return self.data_fn
