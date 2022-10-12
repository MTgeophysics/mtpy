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
import numpy as np
import pandas as pd
import geopandas as gpd

from mtpy.core.mt_location import MTLocation
from mtpy.utils.mtpy_logger import get_mtpy_logger

# =============================================================================


class MTStations:
    """
    station locations class

    """

    def __init__(self, mt_list, model_epsg, datum_epsg=None, **kwargs):

        self.logger = get_mtpy_logger(f"{__name__}.{self.__class__.__name__}")

        self.mt_list = mt_list

        self.dtype = dict(
            [
                ("station", "U50"),
                ("lat", float),
                ("lon", float),
                ("elev", float),
                ("datum_epsg", "U6"),
                ("east", float),
                ("north", float),
                ("utm_epsg", "U6"),
                ("model_east", float),
                ("model_north", float),
                ("model_elev", float),
            ]
        )
        self.model_epsg = model_epsg
        self.datum_epsg = datum_epsg
        self._center_lat = None
        self._center_lon = None
        self._center_elev = 0.0
        self.shift_east = 0
        self.shift_north = 0

        for key in list(kwargs.keys()):
            if hasattr(self, key):
                setattr(self, key, kwargs[key])

        self.calculate_rel_locations()

    def __str__(self):
        fmt_dict = dict(
            [
                ("station", "<8"),
                ("lat", "<10.4f"),
                ("lon", "<10.4f"),
                ("elev", "<8.2f"),
                ("rel_east", "<13.2f"),
                ("rel_north", "<13.2f"),
                ("rel_elev", "<8.2f"),
                ("east", "<12.2f"),
                ("north", "<12.2f"),
                ("zone", "<6"),
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
        for n in ["lat", "lon", "elev", "east", "north", "zone"]:
            l.append(f"{self.center_point[n][0]:{fmt_dict[n]}}")
        lines.append("".join(l))

        lines.append("\nMean Values:")
        l = []
        for n in ["lat", "lon", "elev", "east", "north"]:
            l.append(f"{self.station_locations[n].mean():{fmt_dict[n]}}")
        lines.append("".join(l) + f"{self.center_point.zone[0]:<6}")

        return "\n".join(lines)

    def __repr__(self):
        return self.__str__()

    @property
    def model_epsg(self):
        return self._model_epsg

    @model_epsg.setter
    def model_epsg(self, value):
        """
        set the model epsg number an project east, north
        """
        self._model_epsg = value
        if self._model_epsg is not None:
            for mt_obj in self.mt_list:
                mt_obj.utm_epsg = value

    @property
    def datum_epsg(self):
        return self._datum_epsg

    @datum_epsg.setter
    def datum_epsg(self, value):
        """
        set the model epsg number an project east, north
        """
        self._datum_epsg = value
        if self._datum_epsg is not None:
            for mt_obj in self.mt_list:
                mt_obj.datum_epsg = value

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
            entries["lat"][ii] = mt_obj.latitude
            entries["lon"][ii] = mt_obj.longitude
            entries["elev"][ii] = mt_obj.elevation
            entries["datum_epsg"][ii] = mt_obj.datum_epsg
            entries["east"][ii] = mt_obj.east
            entries["north"][ii] = mt_obj.north
            entries["utm_epsg"][ii] = mt_obj.utm_epsg
            entries["model_east"][ii] = mt_obj.model_east
            entries["model_north"][ii] = mt_obj.model_north
            entries["model_elev"][ii] = mt_obj.model_elevation

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
            self.model_epsg = df.utm_epsg.median()
        else:
            if self.model_epsg is None:
                self.model_epsg = df.utm_epsg.unique()[0]

    def calculate_rel_locations(self, shift_east=0, shift_north=0):
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
            center_location.utm_epsg = self.model_epsg
            center_location.model_east = center_location.east
            center_location.model_north = center_location.north
            center_location.model_elevation = self._center_elev

            return center_location

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

            center_location.datum_epsg = self.datum_epsg
            center_location.utm_epsg = self.model_epsg
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
        station_index_x = []
        station_index_y = []
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

            station_index_x.append(sxi)
            station_index_y.append(syi)

            # update elevation in station locations and data array, +1 m as
            # data elevation needs to be below the topography (as advised by Naser)
            mt_obj.model_elev = topoval + 0.001

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
