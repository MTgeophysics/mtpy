"""
==================
ModEM
==================

# Generate files for ModEM

# revised by JP 2017
# revised by AK 2017 to bring across functionality from ak branch

"""
# =============================================================================
# Imports
# =============================================================================
from pathlib import Path
from copy import deepcopy
import numpy as np
import pandas as pd
from pyproj import CRS
import geopandas as gpd
from scipy import stats
from loguru import logger

from mtpy.core.mt_location import MTLocation

from pyevtk.hl import pointsToVTK

# =============================================================================


class MTStations:
    """
    Object to deal with station location and geographic projection.

    Geographic projections are done using pyproj.CRS objects.

    Takes in a list of :class:`mtpy.core.mt.MT` objects which are inherit
    :class:`mtpy.core.mt_location.MTLocation` objects, which deal with
    transformation of point data using pyproj.

    """

    def __init__(self, utm_epsg, datum_epsg=None, **kwargs):

        self.logger = logger

        self.dtype = dict(
            [
                ("survey", "U50"),
                ("station", "U50"),
                ("latitude", float),
                ("longitude", float),
                ("elevation", float),
                ("datum_epsg", "U6"),
                ("east", float),
                ("north", float),
                ("utm_epsg", "U6"),
                ("model_east", float),
                ("model_north", float),
                ("model_elevation", float),
                ("profile_offset", float),
            ]
        )
        self._datum_crs = CRS.from_epsg(4326)
        self._utm_crs = None
        self._center_lat = None
        self._center_lon = None
        self._center_elev = 0.0
        self.shift_east = 0
        self.shift_north = 0
        self.rotation_angle = 0.0
        self.mt_list = None
        self.utm_epsg = utm_epsg

        for key in list(kwargs.keys()):
            if hasattr(self, key):
                setattr(self, key, kwargs[key])

        if self.mt_list is not None:
            if len(self.mt_list) > 0:
                self.compute_relative_locations()
                self.station_locations

    def __str__(self):
        if self.mt_list is None:
            return ""
        elif len(self.mt_list) == 0:
            return ""

        fmt_dict = dict(
            [
                ("survey", "<8"),
                ("station", "<8"),
                ("latitude", "<10.4f"),
                ("longitude", "<10.4f"),
                ("elevation", "<8.2f"),
                ("model_east", "<13.2f"),
                ("model_north", "<13.2f"),
                ("model_elevation", "<13.2f"),
                ("profile_offset", "<13.2f"),
                ("east", "<12.2f"),
                ("north", "<12.2f"),
                ("utm_epsg", "<6"),
                ("datum_epsg", "<6"),
            ]
        )
        lines = ["  ".join([n for n in self.station_locations.columns])]
        lines.append("-" * 72)
        for row in self.station_locations.itertuples():
            l = []
            for key in self.station_locations.columns:
                l.append(f"{getattr(row, key):{fmt_dict[key]}}")
            lines.append("".join(l))

        lines.append("\nModel Center:")
        l = []
        for n in [
            "latitude",
            "longitude",
            "elevation",
            "east",
            "north",
            "utm_epsg",
        ]:
            l.append(f"{getattr(self.center_point, n):{fmt_dict[n]}}")
        lines.append("".join(l))

        lines.append("\nMean Values:")
        l = []
        for n in ["latitude", "longitude", "elevation", "east", "north"]:
            l.append(f"{self.station_locations[n].mean():{fmt_dict[n]}}")

        return "\n".join(lines)

    def __repr__(self):
        return self.__str__()

    def __eq__(self, other):
        if not isinstance(other, MTStations):
            raise TypeError(f"Can not compare {type(other)} to MTStations")

        if not (self.station_locations == other.station_locations).all().all():
            return False

        if not self.center_point == other.center_point:
            return False

        return True

    def __len__(self):
        if self.mt_list is None:
            return 0
        else:
            return len(self.mt_list)

    def copy(self):
        """
        create a deep copy of the MTStations object.

        .. note:: At the moment this is very slow because it is making a lot
         of deep copies.  Use sparingly.

        :return: deep copy of MTStation object
        :rtype: :class:`mtpy.core.mt_stations.MTStations`

        """

        if self.mt_list is not None:
            mt_list_copy = [m.copy() for m in self.mt_list]
        else:
            mt_list_copy = None
        copied = MTStations(None, mt_list=mt_list_copy)
        for key in [
            "utm_crs",
            "datum_crs",
            "_center_lat",
            "_center_lon",
            "_center_elev",
            "shift_east",
            "shift_north",
            "rotation_angle",
        ]:
            setattr(copied, key, deepcopy(getattr(self, key)))

        return copied

    @property
    def model_epsg(self):
        """

        :return: model epsg number from the model_crs object
        :rtype: int

        """
        return self.utm_epsg

    @model_epsg.setter
    def model_epsg(self, value):
        """

        :param value: EPSG number for the model
        :type value: integer or string

        """
        self.utm_epsg = value

    @property
    def utm_crs(self):
        """

        :return: UTM CRS object
        :rtype: :class:`pyproj.CRS`

        """
        if self._utm_crs is not None:
            return self._utm_crs

    @property
    def utm_name(self):
        """

        :return: UTM CRS name
        :rtype: string

        """
        if self._utm_crs is not None:
            return self._utm_crs.name

    @property
    def utm_epsg(self):
        """

        :return: UTM EPSG number
        :rtype: int

        """
        if self._utm_crs is not None:
            return self._utm_crs.to_epsg()

    @utm_epsg.setter
    def utm_epsg(self, value):
        """

        :param value: EPSG number
        :type value: int or str

        """
        self.utm_crs = value

    @property
    def utm_zone(self):
        """

        :return: UTM Zone number
        :rtype: str

        """
        if self._utm_crs is not None:
            return self._utm_crs.utm_zone

    @utm_crs.setter
    def utm_crs(self, value):
        """

        :param value: UTM CRS object, EPSG number, proj4 string
        :type value: :class:`pyproj.CRS`, int, str

        """
        if value in [None, "None", "none", "null"]:
            return

        self._utm_crs = CRS.from_user_input(value)
        if self.mt_list is not None:
            for mt_obj in self.mt_list:
                mt_obj.utm_crs = value

    @property
    def datum_crs(self):
        """

        :return: Datum CRS object
        :rtype: :class:`pyproj.CRS`

        """
        if self._datum_crs is not None:
            return self._datum_crs

    @property
    def datum_name(self):
        """

        :return: Datum well known name
        :rtype: str

        """
        if self._datum_crs is not None:
            return self._datum_crs.name

    @property
    def datum_epsg(self):
        """

        :return: Datum EPSG number
        :rtype: int

        """
        if self._datum_crs is not None:
            return self._datum_crs.to_epsg()

    @datum_epsg.setter
    def datum_epsg(self, value):
        """

        :param value: Datum EPSG number
        :type value: int or str

        """
        self.datum_crs = value

    @datum_crs.setter
    def datum_crs(self, value):
        """
        set the model epsg number an project east, north

        :param value: Datum CRS object, EPSG number, proj4 string
        :type value: :class:`pyproj.CRS`, int, str
        """
        if value in [None, "None", "none", "null"]:
            return

        self._datum_crs = CRS.from_user_input(value)
        if self.mt_list is not None:
            for mt_obj in self.mt_list:
                mt_obj.datum_crs = value

    @property
    def station_locations(self):
        """

        :return: Dataframe of station location information
        :rtype: :class:`pandas.DataFrame`

        """

        # make a structured array to put station location information into
        if self.mt_list is None:
            return

        entries = dict(
            [
                (col, np.zeros(len(self.mt_list), dtype))
                for col, dtype in self.dtype.items()
            ]
        )
        # get station locations in meters
        for ii, mt_obj in enumerate(self.mt_list):
            entries["survey"][ii] = mt_obj.survey
            entries["station"][ii] = mt_obj.station
            entries["latitude"][ii] = mt_obj.latitude
            entries["longitude"][ii] = mt_obj.longitude
            entries["elevation"][ii] = mt_obj.elevation
            entries["datum_epsg"][ii] = mt_obj.datum_epsg
            entries["east"][ii] = mt_obj.east
            entries["north"][ii] = mt_obj.north
            entries["utm_epsg"][ii] = mt_obj.utm_epsg
            entries["model_east"][ii] = mt_obj.model_east
            entries["model_north"][ii] = mt_obj.model_north
            entries["model_elevation"][ii] = mt_obj.model_elevation
            entries["profile_offset"][ii] = mt_obj.profile_offset

        station_df = pd.DataFrame(entries)
        self.datum_epsg = self._validate_epsg(station_df, key="datum")
        self.utm_epsg = self._validate_epsg(station_df, key="utm")

        return station_df

    def _validate_epsg(self, df, key="datum"):
        """
        Make sure that there is only one EPSG number for each of the Datum
        and UTM.  If there are more than one use the median value or the
        first in a unique list of EPSG numbers

        :param df: station_location dataframe
        :type df: :class:`pandas.DataFrame`
        :return: EPSG number
        :rtype: int

        """

        key = f"{key}_epsg"
        if len(df[key].unique()) > 1:
            epsg = int(df[key].median())
            self.logger.warning(
                f"Found more than one {key} number, using median EPSG number {epsg}"
            )
            return epsg
        else:
            if getattr(self, key) is None:
                return df[key].unique()[0]

    def compute_relative_locations(self):
        """
        Calculate model station locations relative to the center point in meters.

        Uses `mtpy.core.MTLocation.compute_model_location` to calculate the
        relative distance.

        Computes inplace.

        """

        for mt_obj in self.mt_list:
            mt_obj.compute_model_location(self.center_point)

    # make center point a get property, can't set it.
    @property
    def center_point(self):
        """
        calculate the center point from the given station locations

        If _center attributes are set, that is returned as the center point.

        Otherwise, looks for non-zero locations in E-N first, then Lat/Lon and
        estimates the center point as (max - min) / 2.

        :return: Center point
        :rtype: :class:`mtpy.core.MTLocation`

        """

        center_location = MTLocation()
        if self._center_lat is not None and self._center_lon is not None:
            self.logger.debug("assigning center from user set values")
            center_location.latitude = self._center_lat
            center_location.longitude = self._center_lon
            center_location.elevation = self._center_elev
            center_location.utm_epsg = self.utm_epsg
            center_location.model_east = center_location.east
            center_location.model_north = center_location.north
            center_location.model_elevation = self._center_elev

            return center_location

        else:
            center_location.datum_epsg = self.datum_epsg
            center_location.utm_epsg = self.utm_epsg
            st_df = self.station_locations.copy()

            st_en = st_df.loc[(st_df.east != 0) & (st_df.north != 0)]
            if st_en.empty:
                st_ll = st_df.loc[
                    (st_df.latitude != 0) & (st_df.longitude != 0)
                ]
                if st_ll.empty:
                    raise ValueError(
                        "Station locations are all 0 cannot find center."
                    )

                else:
                    self.logger.debug(
                        "locating center from latitude and longitude"
                    )
                    center_location.latitude = (
                        st_ll.latitude.max() + st_ll.latitude.min()
                    ) / 2
                    center_location.longitude = (
                        st_ll.longitude.max() + st_ll.longitude.min()
                    ) / 2

            else:
                self.logger.debug("locating center from UTM grid")
                center_location.east = (st_en.east.max() + st_en.east.min()) / 2
                center_location.north = (
                    st_en.north.max() + st_en.north.min()
                ) / 2

            center_location.model_east = center_location.east
            center_location.model_north = center_location.north
            center_location.model_elevation = self._center_elev

        return center_location

    def rotate_stations(self, rotation_angle):
        """
        Rotate stations model postions only assuming N is 0 and east is 90.

        .. note:: Computes in place and rotates according to already set
         rotation angle.  Therefore if the station locations have already been
         rotated the function will rotate the already rotate stations.  For
         example if you rotate the stations 15 degrees, then again by 20 degrees
         the resulting station locations will be 35 degrees rotated from the
         original locations.

        :param rotation_angle: rotation angle in degrees assuming N=0, E=90.
         Positive clockwise.
        :type rotation_angle: float

        """

        cos_ang = np.cos(np.deg2rad(rotation_angle))
        sin_ang = np.sin(np.deg2rad(rotation_angle))
        rot_matrix = np.array([[cos_ang, sin_ang], [-sin_ang, cos_ang]])

        for mt_obj in self.mt_list:
            coords = np.array(
                [
                    mt_obj.model_east,
                    mt_obj.model_north,
                ]
            )

            # rotate the relative station locations
            new_coords = np.array(np.dot(rot_matrix, coords))

            mt_obj.model_east = new_coords[0]
            mt_obj.model_north = new_coords[1]

        self.rotation_angle += rotation_angle

        self.logger.info(
            f"Rotated stations by {rotation_angle:.1f} deg clockwise from N. "
            f"Total rotation = {self.rotation_angle:.1f} deg."
        )

    def center_stations(self, model_obj):
        """
        Center station locations to the middle of cells, is useful for
        topography cause it reduces edge effects of stations close to cell edges.
        Recalculates rel_east, rel_north to center of model cell.

        :param model_obj: :class:`mtpy.modeling.Structured` object of the model
        :type model_obj: :class:`mtpy.modeling.modem.Model`


        """

        for mt_obj in self.mt_list:
            e_index = (
                np.where(model_obj.grid_east >= mt_obj.model_east)[0][0] - 1
            )
            n_index = (
                np.where(model_obj.grid_north >= mt_obj.model_north)[0][0] - 1
            )

            mt_obj.model_east = model_obj.grid_east[
                e_index : e_index + 2
            ].mean()
            mt_obj.model_north = model_obj.grid_north[
                n_index : n_index + 2
            ].mean()

    def project_stations_on_topography(
        self,
        model_object,
        air_resistivity=1e12,
        sea_resistivity=0.3,
        ocean_bottom=False,
    ):
        """
        Project stations on topography of a given model

        :param model_obj: :class:`mtpy.modeling.modem.Model` object of the model
        :type model_obj: :class:`mtpy.modeling.modem.Model`
        :param air_resistivity: resistivity value of air cells in the model
        :type air_resistivity:  float
        :param sea_resistivity: resistivity of sea
        :type sea_resistivity: float
        :param ocean_bottom: If True places stations at bottom of sea cells
        :type ocean_bottom: boolean

        Recaluclates rel_elev
        """

        # find index of each station on grid
        for mt_obj in self.mt_list:
            # relative locations of stations
            sx = mt_obj.model_east
            sy = mt_obj.model_north

            # indices of stations on model grid
            sxi = np.where(
                (sx <= model_object.grid_east[1:])
                & (sx > model_object.grid_east[:-1])
            )[0][0]

            syi = np.where(
                (sy <= model_object.grid_north[1:])
                & (sy > model_object.grid_north[:-1])
            )[0][0]

            # first, check if there are any air cells
            if np.any(
                model_object.res_model[syi, sxi] > 0.95 * air_resistivity
            ):
                szi = np.amin(
                    np.where(
                        (
                            model_object.res_model[syi, sxi]
                            < 0.95 * air_resistivity
                        )
                    )[0]
                )
            # otherwise place station at the top of the model
            else:
                szi = 0

            # JP: estimate ocean bottom stations if requested
            if ocean_bottom:
                if np.any(model_object.res_model[syi, sxi] <= sea_resistivity):
                    szi = np.amax(
                        np.where(
                            (
                                model_object.res_model[syi, sxi]
                                <= sea_resistivity
                            )
                        )[0]
                    )

            # get relevant grid point elevation
            topoval = model_object.grid_z[szi]

            # update elevation in station locations and data array, +1 m as
            # data elevation needs to be below the topography (as advised by Naser)
            mt_obj.model_elevation = topoval + 0.001

        # BM: After applying topography, center point of grid becomes
        #  highest point of surface model.
        self._center_elev = model_object.grid_z[0]

    def to_geopd(self):
        """
        create a geopandas dataframe

        :return: Geopandas DataFrame with points from latitude and longitude
        :rtype: :class:`geopandas.DataFrame`

        """

        gdf = gpd.GeoDataFrame(
            self.station_locations,
            geometry=gpd.points_from_xy(
                self.station_locations.longitude,
                self.station_locations.latitude,
            ),
            crs=self.center_point.datum_crs,
        )

        return gdf

    def to_shp(self, shp_fn):
        """
        Write a shape file of the station locations using geopandas which only takes
        in epsg numbers

        :param shp_fn: full path to new shapefile
        :type shp_fn: string

        """
        sdf = self.to_geopd()

        sdf.to_file(shp_fn)
        return shp_fn

    def to_csv(self, csv_fn, geometry=False):
        """
        Write a shape file of the station locations using geopandas which only takes
        in epsg numbers

        :param csv_fn: full path to new shapefile
        :type csv_fn: string

        """
        sdf = self.to_geopd()
        use_columns = list(sdf.columns)
        if not geometry:
            use_columns.remove("geometry")
        sdf.to_csv(csv_fn, index=False, columns=use_columns)

    def to_vtk(
        self,
        vtk_fn=None,
        vtk_save_path=None,
        vtk_fn_basename="ModEM_stations",
        geographic=False,
        shift_east=0,
        shift_north=0,
        shift_elev=0,
        units="km",
        coordinate_system="nez+",
    ):
        """

        :param vtk_save_path: directory to save vtk file to, defaults to None
        :type vtk_save_path: string or Path, optional
        :param vtk_fn_basename: filename basename of vtk file, note that .vtr
        extension is automatically added, defaults to "ModEM_stations"
        :type vtk_fn_basename: string, optional
        :param geographic: If true puts the grid on geographic coordinates based
        on the model_utm_zone, defaults to False
        :type geographic: boolean, optional
        :param shift_east: shift in east directions in meters, defaults to 0
        :type shift_east: float, optional
        :param shift_north: shift in north direction in meters, defaults to 0
        :type shift_north: float, optional
        :param shift_elev: shift in elevation + down in meters, defaults to 0
        :type shift_elev: float, optional
        :param units: Units of the spatial grid [ km | m | ft ], defaults to "km"
        :type units: string, optional
        :type : string
        :param coordinate_system: coordinate system for the station, either the
        normal MT right-hand coordinate system with z+ down or the sinister
        z- down [ nez+ | enz- ], defaults to nez+
        :return: full path to VTK file
        :rtype: Path

        Write VTK file
        >>> md.write_vtk_station_file(vtk_fn_basename="modem_stations")

        Write VTK file in geographic coordinates
        >>> md.write_vtk_station_file(vtk_fn_basename="modem_stations",
        >>> ...                       geographic=True)

        Write VTK file in geographic coordinates with z+ up
        >>> md.write_vtk_station_file(vtk_fn_basename="modem_stations",
        >>> ...                       geographic=True,
        >>> ...                       coordinate_system='enz-')

        """

        if isinstance(units, str):
            if units.lower() == "km":
                scale = 1.0 / 1000.00
            elif units.lower() == "m":
                scale = 1.0
            elif units.lower() == "ft":
                scale = 3.2808
        elif isinstance(units, (int, float)):
            scale = units

        if vtk_fn is None:
            if vtk_save_path is None:
                raise ValueError("Need to input vtk_save_path")
            vtk_fn = Path(vtk_save_path, vtk_fn_basename)
        else:
            vtk_fn = Path(vtk_fn)

        if vtk_fn.suffix != "":
            vtk_fn = vtk_fn.parent.joinpath(vtk_fn.stem)

        sdf = self.station_locations.copy()

        if not geographic:
            if coordinate_system == "nez+":
                vtk_x = (sdf.model_north + shift_north) * scale
                vtk_y = (sdf.model_east + shift_east) * scale
                vtk_z = (sdf.model_elevation + shift_elev) * scale
                extra = (sdf.model_elevation + shift_elev) * scale
            elif coordinate_system == "enz-":
                vtk_x = (sdf.model_north + shift_north) * scale
                vtk_y = (sdf.model_east + shift_east) * scale
                vtk_z = -1 * (sdf.model_elevation + shift_elev) * scale
                extra = -1 * (sdf.model_elevation + shift_elev) * scale

        else:
            if coordinate_system == "nez+":
                vtk_y = (sdf.north + shift_north) * scale
                vtk_x = (sdf.east + shift_east) * scale
                vtk_z = -1 * (sdf.elevation + shift_elev) * scale
                extra = -1 * (sdf.elevation + shift_elev)
            elif coordinate_system == "enz-":
                vtk_y = (sdf.north + shift_north) * scale
                vtk_x = (sdf.east + shift_east) * scale
                vtk_z = (sdf.elevation + shift_elev) * scale
                extra = sdf.elevation + shift_elev

        # write file
        pointsToVTK(
            vtk_fn.as_posix(),
            vtk_x.to_numpy(),
            vtk_y.to_numpy(),
            vtk_z.to_numpy(),
            data={"elevation": extra.to_numpy()},
        )

        self.logger.info(f"Wrote station VTK file to {vtk_fn}.vtu")
        return vtk_fn

    def generate_profile(self, units="deg"):
        """
        Estimate a profile from the data
        :return: DESCRIPTION
        :rtype: TYPE

        """

        if units == "deg":
            x = self.station_locations.longitude
            y = self.station_locations.latitude

        elif units == "m":
            if self.utm_crs is not None:
                x = self.station_locations.east
                y = self.station_locations.north
            else:
                raise ValueError("Must input a UTM CRS or EPSG")

        # check regression for 2 profile orientations:
        # horizontal (N=N(E)) or vertical(E=E(N))
        # use the one with the lower standard deviation
        profile1 = stats.linregress(x, y)
        profile2 = stats.linregress(y, x)
        # if the profile is rather E=E(N), the parameters have to converted
        # into N=N(E) form:
        if profile2.stderr < profile1.stderr:
            profile_line = {
                "slope": 1.0 / profile2.slope,
                "intercept": -profile2.intercept / profile2.slope,
            }
        else:
            profile_line = {
                "slope": profile1.slope,
                "intercept": profile1.intercept,
            }

        x1 = x.min()
        x2 = x.max()
        y1 = profile_line["slope"] * x1 + profile_line["intercept"]
        y2 = profile_line["slope"] * x2 + profile_line["intercept"]

        return x1, y1, x2, y2, profile_line

    def generate_profile_from_strike(self, strike, units="deg"):
        """
        Estimate a profile line from a given geoelectric strike

        :param units: DESCRIPTION, defaults to "deg"
        :type units: TYPE, optional
        :return: DESCRIPTION
        :rtype: TYPE

        """

        if units == "deg":
            x = self.station_locations.longitude
            y = self.station_locations.latitude

        elif units == "m":
            if self.utm_crs is not None:
                x = self.station_locations.east
                y = self.station_locations.north
            else:
                raise ValueError("Must input a UTM CRS or EPSG")

        profile_line = {"slope": np.arctan(np.deg2rad(90 - strike))}
        profile_line["intercept"] = y.min() - profile_line["slope"] * x.min()

        x1 = x.min()
        x2 = x.max()
        y1 = profile_line["slope"] * x1 + profile_line["intercept"]
        y2 = profile_line["slope"] * x2 + profile_line["intercept"]

        return x1, y1, x2, y2, profile_line

    def _extract_profile(self, x1, y1, x2, y2, radius):
        """
        extract stations along a profile line that lie with in the given
        radius

        :param point1: DESCRIPTION
        :type point1: TYPE
        :param point2: DESCRIPTION
        :type point2: TYPE
        :param radius: DESCRIPTION
        :type radius: TYPE
        :return: DESCRIPTION
        :rtype: TYPE

        """

        if np.abs(x2 - x1) < 100:
            if self.utm_crs is None:
                raise ValueError("Must input UTM CRS or EPSG.")
            point_1 = MTLocation(
                longitude=x1, latitude=y1, utm_crs=self.utm_crs
            )
            point_2 = MTLocation(
                longitude=x2, latitude=y2, utm_crs=self.utm_crs
            )
            x1 = point_1.east
            y1 = point_1.north
            x2 = point_2.east
            y2 = point_2.north

        if radius is None:
            radius = 1e12

        def distance(x, y):
            return np.abs(
                (x2 - x1) * (y1 - y) - (x1 - x) * (y2 - y1)
            ) / np.sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2)

        slope = (y2 - y1) / (x2 - x1)
        intersection = y1 - slope * x1

        profile_list = []
        offsets = []
        for mt_obj in self.mt_list:
            d = distance(mt_obj.east, mt_obj.north)

            if d <= radius:
                mt_obj.project_onto_profile_line(slope, intersection)
                profile_list.append(mt_obj)
                offsets.append(mt_obj.profile_offset)

        offsets = np.array(offsets)
        indexes = np.argsort(offsets)

        sorted_profile_list = []
        for index in indexes:
            profile_list[index].profile_offset -= offsets.min()
            sorted_profile_list.append(profile_list[index])

        return sorted_profile_list
