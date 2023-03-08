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
import numpy as np
import pandas as pd
from pyproj import CRS
import geopandas as gpd
from scipy import stats

from mtpy.core.mt_location import MTLocation
from mtpy.utils.mtpy_logger import get_mtpy_logger

try:
    from pyevtk.hl import pointsToVTK
except ImportError:
    print(
        "If you want to write a vtk file for 3d viewing, you need to install pyevtk"
    )
# =============================================================================


class MTStations:
    """
    station locations class

    """

    def __init__(self, utm_epsg, datum_epsg=None, **kwargs):

        self.logger = get_mtpy_logger(f"{__name__}.{self.__class__.__name__}")

        self.dtype = dict(
            [
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

        for key in list(kwargs.keys()):
            if hasattr(self, key):
                setattr(self, key, kwargs[key])

        if self.mt_list is not None:
            self.compute_relative_locations()

    def __str__(self):
        fmt_dict = dict(
            [
                ("station", "<8"),
                ("latitude", "<10.4f"),
                ("longitude", "<10.4f"),
                ("elevation", "<8.2f"),
                ("model_east", "<13.2f"),
                ("model_north", "<13.2f"),
                ("model_elevation", "<8.2f"),
                ("east", "<12.2f"),
                ("north", "<12.2f"),
                ("utm_epsg", "<6"),
            ]
        )
        lines = [
            "".join(
                [
                    f"{n.capitalize():<10}"
                    for n in self.station_locations.dtype.names
                ]
            )
        ]
        lines.append("-" * 72)
        for ss in self.station_locations:
            l = []
            for k in self.station_locations.dtype.names:
                l.append(f"{ss[k]:{fmt_dict[k]}}")
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
            l.append(f"{self.center_point[n][0]:{fmt_dict[n]}}")
        lines.append("".join(l))

        lines.append("\nMean Values:")
        l = []
        for n in ["latitude", "longitude", "elevation", "east", "north"]:
            l.append(f"{self.station_locations[n].mean():{fmt_dict[n]}}")

        return "\n".join(lines)

    def __repr__(self):
        return self.__str__()

    @property
    def model_epsg(self):
        return self.utm_epsg

    @model_epsg.setter
    def model_epsg(self, value):
        self.utm_epsg = value

    @property
    def utm_crs(self):
        if self._utm_crs is not None:
            return self._utm_crs

    @property
    def utm_name(self):
        if self._utm_crs is not None:
            return self._utm_crs.name

    @property
    def utm_epsg(self):
        if self._utm_crs is not None:
            return self._utm_crs.to_epsg()

    @utm_epsg.setter
    def utm_epsg(self, value):
        self.utm_crs = value

    @property
    def utm_zone(self):
        if self._utm_crs is not None:
            return self._utm_crs.utm_zone

    @utm_crs.setter
    def utm_crs(self, value):
        if value in [None, "None", "none", "null"]:
            return

        self._utm_crs = CRS.from_user_input(value)
        for mt_obj in self.mt_list:
            mt_obj.utm_crs = value

    @property
    def datum_crs(self):
        if self._datum_crs is not None:
            return self._datum_crs

    @property
    def datum_name(self):
        if self._datum_crs is not None:
            return self._datum_crs.name

    @property
    def datum_epsg(self):
        if self._datum_crs is not None:
            return self._datum_crs.to_epsg()

    @datum_epsg.setter
    def datum_epsg(self, value):
        self.datum_crs = value

    @datum_crs.setter
    def datum_crs(self, value):
        """
        set the model epsg number an project east, north
        """
        if value in [None, "None", "none", "null"]:
            return

        self._datum_crs = CRS.from_user_input(value)
        for mt_obj in self.mt_list:
            mt_obj.datum_crs = value

    @property
    def station_locations(self):
        """
        get station locations from a list of edi files

        Arguments
        -------------
            **input_list** : list
                             list of edi file names, or mt_objects


        Returns
        ------------
            * fills station_locations array

        """

        # make a structured array to put station location information into
        entries = dict(
            [
                (col, np.zeros(len(self.mt_list), dtype))
                for col, dtype in self.dtype.items()
            ]
        )
        # get station locations in meters
        for ii, mt_obj in enumerate(self.mt_list):
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
        self._validate_station_locations(station_df)

        return station_df

    def _validate_station_locations(self, df):

        if len(df.datum_epsg.unique()) > 1:
            self.datum_epsg = df.datum_epsg.median()
        else:
            if self.datum_epsg is None:
                self.datum_epsg = df.datum_epsg.unique()[0]

        if len(df.utm_epsg.unique()) > 1:
            self.utm_epsg = df.utm_epsg.median()
        else:
            if self.utm_epsg is None:
                self.utm_epsg = df.utm_epsg.unique()[0]

    def compute_relative_locations(self, shift_east=0, shift_north=0):
        """
        put station in a coordinate system relative to
        (shift_east, shift_north)
        (+) shift right or up
        (-) shift left or down

        """

        for mt_obj in self.mt_list:
            mt_obj.compute_model_location(self.center_point)

    # make center point a get property, can't set it.
    @property
    def center_point(self):
        """
        calculate the center point from the given station locations

        Returns
        -------------
            **center_location** : np.ndarray
                                  structured array of length 1
                                  dtype includes (east, north, zone, lat, lon)
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
            if np.all(self.station_locations.east == 0) or np.all(
                self.station_locations.north == 0
            ):
                if np.all(self.station_locations.latitude != 0) and np.all(
                    self.station_locations.longitude != 0
                ):
                    self.logger.debug(
                        "locating center from latitude and longitude"
                    )
                    center_location.latitude = (
                        self.station_locations.latitude.max()
                        + self.station_locations.latitude.min()
                    ) / 2
                    center_location.longitude = (
                        self.station_locations.longitude.max()
                        + self.station_locations.longitude.min()
                    ) / 2
                else:
                    raise ValueError(
                        "Station locations are all 0 cannot find center."
                    )

            else:
                self.logger.debug("locating center from UTM grid")
                center_location.east = (
                    self.station_locations.east.max()
                    + self.station_locations.east.min()
                ) / 2
                center_location.north = (
                    self.station_locations.north.max()
                    + self.station_locations.north.min()
                ) / 2

            center_location.model_east = center_location.east
            center_location.model_north = center_location.north
            center_location.model_elevation = self._center_elev

        return center_location

    def rotate_stations(self, rotation_angle):
        """
        Rotate stations assuming N is 0

        Arguments
        -------------
            **rotation_angle** : float
                                 angle in degrees assuming N is 0

        Returns
        -------------
            * refils rel_east and rel_north in station_locations.  Does this
              because you will still need the original locations for plotting
              later.

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

        self.logger.info(
            f"Rotated stations by {rotation_angle:.1f} deg clockwise from N"
        )

    def center_stations(self, model_obj):
        """
        Center station locations to the middle of cells, is useful for
        topography cause it reduces edge effects of stations close to cell edges.
        Recalculates rel_east, rel_north to center of model cell.

        :param model_obj: :class:`mtpy.modeling.modem.Model` object of the model
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
            print(f"{mt_obj.station}: {mt_obj.model_elevation:.2f}")

        # BM: After applying topography, center point of grid becomes
        #  highest point of surface model.
        self._center_elev = model_object.grid_z[0]

    def to_geopd(self):
        """
        create a geopandas dataframe

        """

        gdf = gpd.GeoDataFrame(
            self.station_locations,
            geometry=gpd.points_from_xy(
                self.station_locations.lon, self.station_locations.lat
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

    def write_vtk_station_file(
        self,
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

        if vtk_save_path is None:
            vtk_fn = self.save_path.joinpath(vtk_fn_basename)
        else:
            vtk_fn = Path(vtk_save_path, vtk_fn_basename)

        sdf = self.station_locations.copy()

        if not geographic:
            if coordinate_system == "nez+":
                vtk_x = (sdf.model_north + shift_north) * scale
                vtk_y = (sdf.model_east + shift_east) * scale
                vtk_z = (sdf.model_elev + shift_elev) * scale
                extra = (sdf.model_elev + shift_elev) * scale
            elif coordinate_system == "enz-":
                vtk_x = (sdf.model_north + shift_north) * scale
                vtk_y = (sdf.model_east + shift_east) * scale
                vtk_z = (sdf.model_elev + shift_elev) * scale
                extra = (sdf.model_elev + shift_elev) * scale

        else:
            if coordinate_system == "nez+":
                vtk_y = (sdf.north + shift_north) * scale
                vtk_x = (sdf.east + shift_east) * scale
                vtk_z = -1 * (sdf.elev + shift_elev) * scale
                extra = -1 * (sdf.elev + shift_elev)
            elif coordinate_system == "enz-":
                vtk_y = (sdf.north + shift_north) * scale
                vtk_x = (sdf.east + shift_east) * scale
                vtk_z = -1 * (sdf.elev + shift_elev) * scale
                extra = -1 * (sdf.elev + shift_elev)

        # write file
        pointsToVTK(
            vtk_fn.as_posix(), vtk_x, vtk_y, vtk_z, data={"elevation": extra}
        )

        self.logger.info("Wrote station VTK file to {0}".format(vtk_fn))
        return vtk_fn

    def _generate_profile(self, units="deg"):
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
                "slope": profile_line[0],
                "intercept": profile_line[1],
            }

        x1 = x.min()
        x2 = x.max()
        y1 = profile_line["slope"] * x1 + profile_line["intercept"]
        y2 = profile_line["slope"] * x2 + profile_line["intercept"]

        return x1, y1, x2, y2, profile_line

    def _generate_profile_from_strike(self, strike, units="deg"):
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

        profile_line = {"slope": strike}
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
