"""
==================
ModEM
==================

# Generate files for ModEM

# revised by JP 2017
# revised by AK 2017 to bring across functionality from ak branch

"""
import numpy as np
from mtpy.core import mt as mt
from mtpy.utils import gis_tools as gis_tools
# in module imports
from .exception import ModEMError

__all__ = ['Stations']


class Stations(object):
    """
    station locations class

    ..note:: If the survey steps across multiple UTM zones, then a
             distance will be added to the stations to place them in
             the correct location.  This distance is
             _utm_grid_size_north and _utm_grid_size_east.  You should
             these parameters to place the locations in the proper spot
             as grid distances and overlaps change over the globe.
             **This is not implemented yet**
    """

    def __init__(self, **kwargs):

        self.dtype = [('station', '|S10'),
                      ('lat', np.float),
                      ('lon', np.float),
                      ('elev', np.float),
                      ('rel_east', np.float),
                      ('rel_north', np.float),
                      ('rel_elev', np.float),
                      ('east', np.float),
                      ('north', np.float),
                      ('zone', 'S4')]
        self.station_locations = np.zeros(0, dtype=self.dtype)
        self.model_epsg = None
        self.model_utm_zone = None

        for key in list(kwargs.keys()):
            if hasattr(self, key):
                setattr(self, key, kwargs[key])

    ## --> define properties that can only be returned and not set
    @property
    def lat(self):
        return self.station_locations['lat']

    @property
    def lon(self):
        return self.station_locations['lon']

    @property
    def east(self):
        return self.station_locations['east']

    @property
    def north(self):
        return self.station_locations['north']

    @property
    def elev(self):
        return self.station_locations['elev']

    @property
    def rel_east(self):
        return self.station_locations['rel_east']

    @property
    def rel_north(self):
        return self.station_locations['rel_north']

    @property
    def rel_elev(self):
        return self.station_locations['rel_elev']

    @property
    def utm_zone(self):
        return self.station_locations['zone']

    @property
    def station(self):
        return self.station_locations['station']

    def _get_mt_objs_from_list(self, input_list):
        """
        get mt_objects from a list of files or mt_objects
        """

        if type(input_list) not in [list, np.ndarray]:
            raise ValueError('Input list needs to be type list, not {0}'.format(type(input_list)))

        if type(input_list[0]) is mt.MT:
            return input_list

        if type(input_list[0]) is str:
            if input_list[0].endswith('.edi'):
                return [mt.MT(fn) for fn in input_list]

            else:
                raise ModEMError('file {0} not supported yet'.format(input_list[0][-4:]))

    def get_station_locations(self, input_list):
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
        # print input_list
        mt_obj_list = self._get_mt_objs_from_list(input_list)

        # if station locations are not input read from the edi files
        if mt_obj_list is None:
            raise AttributeError('mt_obj_list is None, need to input a list of '
                                 'mt objects to read in.')

        n_stations = len(mt_obj_list)

        if n_stations == 0:
            raise ModEMError('No .edi files in edi_list, please check '
                             'file locations.')

        # make a structured array to put station location information into
        self.station_locations = np.zeros(n_stations,
                                          dtype=self.dtype)
        # get station locations in meters
        for ii, mt_obj in enumerate(mt_obj_list):
            self.station_locations[ii]['lat'] = mt_obj.lat
            self.station_locations[ii]['lon'] = mt_obj.lon
            self.station_locations[ii]['station'] = mt_obj.station
            self.station_locations[ii]['elev'] = mt_obj.elev

            if (self.model_epsg is not None) or (self.model_utm_zone is not None):
                east, north, utm_zone = gis_tools.project_point_ll2utm(mt_obj.lat,
                                                                       mt_obj.lon,
                                                                       utm_zone=self.model_utm_zone,
                                                                       epsg=self.model_epsg)
                self.station_locations[ii]['east'] = east
                self.station_locations[ii]['north'] = north
                self.station_locations[ii]['zone'] = utm_zone
            else:
                self.station_locations[ii]['east'] = mt_obj.east
                self.station_locations[ii]['north'] = mt_obj.north
                self.station_locations[ii]['zone'] = mt_obj.utm_zone

        # get relative station locations
        self.calculate_rel_locations()

    def calculate_rel_locations(self, shift_east=0, shift_north=0):
        """
        put station in a coordinate system relative to
        (shift_east, shift_north)
        (+) shift right or up
        (-) shift left or down

        """
        #
        #        #remove the average distance to get coordinates in a relative space
        #        self.station_locations['rel_east'] = self.east-self.east.mean()
        #        self.station_locations['rel_north'] = self.north-self.north.mean()
        #
        #        #translate the stations so they are relative to 0,0
        #        east_center = (self.rel_east.max()-np.abs(self.rel_east.min()))/2.
        #        north_center = (self.rel_north.max()-np.abs(self.rel_north.min()))/2.
        #
        #
        #        #remove the average distance to get coordinates in a relative space
        #        self.station_locations['rel_east'] -= east_center+shift_east
        #        self.station_locations['rel_north'] -= north_center+shift_north

        # translate the stations so they are relative to 0,0
        east_center = (self.east.max() + self.east.min()) / 2.
        north_center = (self.north.max() + self.north.min()) / 2.

        self.station_locations['rel_east'] = self.east - east_center
        self.station_locations['rel_north'] = self.north - north_center
        
        # BM: Before topograhy is applied to the model, the station
        #  elevation isn't relative to anything (according to 
        #  Data.project_stations_on_topography, station elevation is
        #  relevant to topography). So rel_elev and elev are the same.
        #  Once topography has been applied, rel_elev can be calcuated
        #  by calling Data.project_stations_on_topography.
        self.station_locations['rel_elev'] = self.elev

    # make center point a get method, can't set it.
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
        dtype = [('lat', np.float),
                 ('lon', np.float),
                 ('east', np.float),
                 ('north', np.float),
                 ('elev', np.float),
                 ('zone', 'S4')]
        center_location = np.recarray(1, dtype=dtype)
        #        AK - using the mean here but in get_relative_locations used (max + min)/2, why???

        #        center_point = np.array([self.east.mean(), self.north.mean()])
        #
        #        #translate the stations so they are relative to 0,0
        #        east_center = (self.rel_east.max()-np.abs(self.rel_east.min()))/2
        #        north_center = (self.rel_north.max()-np.abs(self.rel_north.min()))/2
        #
        #        center_point[0] -= east_center
        #        center_point[1] -= north_center
        #
        #        # calculate center point in lat, lon, easting, northing
        #        center_location['east'] = center_point[0]
        #        center_location['north'] = center_point[1]
        center_point = np.array([self.east.max() + self.east.min(),
                                 self.north.max() + self.north.min()]) / 2.
        center_location['east'] = center_point[0]
        center_location['north'] = center_point[1]

        center_location['zone'] = self.utm_zone[0]

        center_ll = gis_tools.project_point_utm2ll(float(center_point[0]),
                                                   float(center_point[1]),
                                                   self.utm_zone[0],
                                                   epsg=self.model_epsg)

        center_location['lat'] = center_ll[0]
        center_location['lon'] = center_ll[1]
        # BM: Because we are now writing center_point.elev to ModEm
        #  data file, we need to provide it.
        #  The center point elevation is the highest point of the
        #  model. Before topography is applied, this is the highest
        #  station. After it's applied, it's the highest point
        #  point of the surface model (this will be set by calling
        #  Data.project_stations_on_topography).
        center_location['elev'] = self.elev.max()

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
        rot_matrix = np.array([[cos_ang, sin_ang],
                               [-sin_ang, cos_ang]])

        coords = np.array([self.station_locations['rel_east'],
                           self.station_locations['rel_north']])

        # rotate the relative station locations
        new_coords = np.array(np.dot(rot_matrix, coords))


        self.station_locations['rel_east'] = new_coords[0, :]
        self.station_locations['rel_north'] = new_coords[1, :]
        

        print('Rotated stations by {0:.1f} deg clockwise from N'.format(
            rotation_angle))
        


    def check_utm_crossing(self):
        """
        If the stations cross utm zones, then estimate distance by computing
        distance on a sphere.
        """
        #
        #        latMid = (Lat1+Lat2 )/2.0;  // or just use Lat1 for slightly less accurate estimate
        #
        #
        #        m_per_deg_lat = 111132.954 - 559.822 * cos( 2.0 * latMid ) + 1.175 * cos( 4.0 * latMid);
        #        m_per_deg_lon = (3.14159265359/180 ) * 6367449 * cos ( latMid );
        #
        #        deltaLat = fabs(Lat1 - Lat2);
        #        deltaLon = fabs(Lon1 - Lon2);
        #
        #        dist_m = sqrt (  pow( deltaLat * m_per_deg_lat,2) + pow( deltaLon * m_per_deg_lon , 2) );
        #
        pass
