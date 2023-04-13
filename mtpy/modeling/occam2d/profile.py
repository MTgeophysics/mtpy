# -*- coding: utf-8 -*-
"""
Created on Tue Mar  7 18:08:55 2023

@author: jpeacock
"""
# =============================================================================
# Imports
# =============================================================================
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import scipy as sp
from scipy.stats import mode

import mtpy.analysis.geometry as MTgy
import mtpy.core.mt as mt
from mtpy.utils.calculator import centre_point
from mtpy.utils.gis_tools import get_epsg, project_point_ll2utm

# =============================================================================
class Profile:
    """
    Takes data from .edi files to create a profile line for 2D modeling.
    Can project the stations onto a profile that is perpendicular to strike
    or a given profile direction.

    If _rotate_to_strike is True, the impedance tensor and tipper are rotated
    to align with the geoelectric strike angle.

    If _rotate_to_strike is True and geoelectric_strike is not given,
    then it is calculated using the phase tensor.  First, 2D sections are
    estimated from the impedance tensor then the strike is estimated from the
    phase tensor azimuth + skew.  This angle is then used to project the
    stations perpendicular to the strike angle.

    If you want to project onto an angle not perpendicular to strike, give
    profile_angle and set _rotate_to_strike to False.  This will project
    the impedance tensor and tipper to be perpendicular with the
    profile_angle.

    Arguments:
    -----------

        **edi_path** : string
                       full path to edi files

        **station_list** : list of stations to create profile for if None is
                           given all .edi files in edi_path will be used.
                           .. note:: that algorithm assumes .edi files are
                                     named by station and it only looks for
                                     the station within the .edi file name
                                     it does not match exactly, so if you have
                                     .edi files with similar names there
                                     might be some problems.

        **geoelectric_strike** : float
                                 geoelectric strike direction in degrees
                                 assuming 0 is North and East is 90

        **profile_angle** : float
                            angle to project the stations onto a profile line
                            .. note:: the geoelectric strike angle and
                                      profile angle should be orthogonal for
                                      best results from 2D modeling.


    ======================= ===================================================
    **Attributes**          Description
    ======================= ===================================================
    edi_list                list of mtpy.core.mt.MT instances for each .edi
                            file found in edi_path
    elevation_model         numpy.ndarray(3, num_elevation_points) elevation
                            values for the profile line (east, north, elev)
    geoelectric_strike      geoelectric strike direction assuming N == 0
    profile_angle           angle of profile line assuming N == 0
    profile_line            (slope, N-intercept) of profile line
    _profile_generated      [ True | False ] True if profile has already been
                            generated
    edi_path                path to find .edi files
    station_list            list of stations to extract from edi_path
    num_edi                 number of edi files to create a profile for
    _rotate_to_strike       [ True | False] True to project the stations onto
                            a line that is perpendicular to geoelectric strike
                            also Z and Tipper are rotated to strike direction.
    ======================= ===================================================

    .. note:: change _rotate_to_strike to False if you want to project the
              stations onto a given profile direction.  This will rotate
              Z and Tipper to be orthogonal to this direction

    ======================= ===================================================
    Methods                 Description
    ======================= ===================================================
    generate_profile        generates a profile for the given stations
    plot_profile            plots the profile line along with original station
                            locations to compare.
    ======================= ===================================================

    :Example: ::

        >>> import mtpy.modeling.occam2d as occam
        >>> edi_path = r"/home/mt/edi_files"
        >>> station_list = ['mt{0:03}'.format(ss) for ss in range(0, 15)]
        >>> prof_line = occam.Profile(edi_path, station_list=station_list)
        >>> prof_line.plot_profile()
        >>> #if you want to project to a given strike
        >>> prof_line.geoelectric_strike = 36.7
        >>> prof_line.generate_profile()
        >>> prof_line.plot_profile()


    """

    def __init__(self, dataframe=None, **kwargs):

        self.station_list = None
        self.geoelectric_strike = None
        self.profile_angle = None
        self.model_epsg = None
        self._rotate_to_strike = True
        self.num_edi = 0
        self._profile_generated = False
        self.profile_line = None
        self.station_locations = None
        self.elevation_model = None
        self.elevation_profile = None
        self.estimate_elevation = True

    def _get_edi_list(self):
        """
        get a list of edi files that coorespond to the station list

        each element of the list is a mtpy.core.mt.MT object
        """
        if self.edi_path is not None:
            self.edi_list = []
            if self.station_list is not None:
                for station in self.station_list:
                    tmp_edi_list = []
                    for edi in os.listdir(self.edi_path):
                        # make a temporary list to find all station/edi matches
                        if edi.find(station) == 0 and edi[-3:] == "edi":
                            tmp_edi_list.append(edi)

                    if len(tmp_edi_list) == 0:
                        print(
                            "Didn't find edi file {} in directory {}".format(
                                edi, self.edi_path
                            )
                        )
                    # if only one matching edi use that one
                    elif len(tmp_edi_list) == 1:
                        edi_to_use = tmp_edi_list[0]
                    # if more than one matching edi then find exact match
                    else:
                        found_edi = False
                        for edifile in tmp_edi_list:
                            if edifile[:-4] == station:
                                found_edi = True
                                edi_to_use = edifile
                                break

                        # if no exact matches, raise an error
                        if not found_edi:
                            raise OccamInputError(
                                "Invalid station name in station list"
                            )

                    # append the correct edi
                    self.edi_list.append(
                        mt.MT(os.path.join(self.edi_path, edi_to_use))
                    )

            else:
                self.edi_list = [
                    mt.MT(os.path.join(self.edi_path, edi))
                    for edi in os.listdir(self.edi_path)
                    if edi[-3:] == "edi"
                ]
                self.station_list = [mtObj.station for mtObj in self.edi_list]
        elif self.edi_list is not None and not self.edi_list:
            # use existing edi list
            if self.station_list is not None:
                # filter edi list by station
                filtered_edi_list = []
                for edi in self.edi_list:
                    if edi.station in self.station_list:
                        filtered_edi_list.append(edi)
                self.edi_list = filtered_edi_list

        self.num_edi = len(self.edi_list)

        for edi in self.edi_list:
            if type(edi.Tipper.rotation_angle) is list:
                edi.Tipper.rotation_angle = np.array(edi.Tipper.rotation_angle)

    def generate_profile(self):
        """
        Generate linear profile by regression of station locations.

        If profile_angle is not None, then station are projected onto that
        line.  Else, the a geoelectric strike is calculated from the data
        and the stations are projected onto an angle perpendicular to the
        estimated strike direction.  If _rotate_to_strike is True, the
        impedance tensor and Tipper data are rotated to align with strike.
        Else, data is not rotated to strike.

        To project stations onto a given line, set profile_angle and
        _rotate_to_strike to False.  This will project the stations onto
        profile_angle and rotate the impedance tensor and tipper to be
        perpendicular to the profile_angle.
        """

        self._get_edi_list()

        strike_angles = np.zeros(self.num_edi)

        easts = np.zeros(self.num_edi)
        norths = np.zeros(self.num_edi)
        utm_zones = np.zeros(self.num_edi)

        if self.model_epsg is None:
            latlist = np.array([mtObj.latitude for mtObj in self.edi_list])
            lonlist = np.array([mtObj.longitude for mtObj in self.edi_list])
            lonc, latc = centre_point(lonlist, latlist)
            self.model_epsg = get_epsg(latc, lonc)

        for ii, edi in enumerate(self.edi_list):
            # find strike angles for each station if a strike angle is not
            # given
            if self.geoelectric_strike is None:
                try:
                    # check dimensionality to be sure strike is estimate for 2D
                    dim = MTgy.dimensionality(z_object=edi.Z)
                    # get strike for only those periods
                    gstrike = MTgy.strike_angle(edi.Z.z[np.where(dim == 2)])[
                        :, 0
                    ]
                    strike_angles[ii] = np.median(gstrike)
                except:
                    pass

            if self.model_epsg is not None:
                edi.east, edi.north, edi.utm_zone = project_point_ll2utm(
                    edi.latitude, edi.longitude, epsg=self.model_epsg
                )

            easts[ii] = edi.east
            norths[ii] = edi.north
            utm_zones[ii] = int(edi.utm_zone[:-1])

        if len(self.edi_list) == 0:
            raise IOError(
                "Could not find and .edi file in {0}".format(self.edi_path)
            )

        if self.geoelectric_strike is None:
            try:
                # might try mode here instead of mean
                self.geoelectric_strike = np.median(
                    strike_angles[np.nonzero(strike_angles)]
                )
            except:
                # empty list or so....
                # can happen, if everyhing is just 1D
                self.geoelectric_strike = 0.0

        # need to check the zones of the stations
        main_utmzone = mode(utm_zones)[0][0]

        for ii, zone in enumerate(utm_zones):
            if zone == main_utmzone:
                continue
            else:
                print(
                    (
                        "station {0} is out of main utm zone".format(
                            self.edi_list[ii].station
                        )
                        + " will not be included in profile"
                    )
                )

        # check regression for 2 profile orientations:
        # horizontal (N=N(E)) or vertical(E=E(N))
        # use the one with the lower standard deviation
        profile1 = sp.stats.linregress(easts, norths)
        profile2 = sp.stats.linregress(norths, easts)
        profile_line = profile1[:2]
        # if the profile is rather E=E(N), the parameters have to converted
        # into N=N(E) form:
        if profile2[4] < profile1[4]:
            profile_line = (1.0 / profile2[0], -profile2[1] / profile2[0])
        self.profile_line = profile_line
        # profile_line = sp.polyfit(lo_easts, lo_norths, 1)
        if self.profile_angle is None:
            self.profile_angle = (
                90 - (np.arctan(profile_line[0]) * 180 / np.pi)
            ) % 180
        # rotate Z according to strike angle,

        # if strike was explicitely given, use that value!

        # otherwise:
        # have 90 degree ambiguity in strike determination
        # choose strike which offers larger angle with profile
        # if profile azimuth is in [0,90].

        if self._rotate_to_strike is False:
            if 0 <= self.profile_angle < 90:
                if np.abs(self.profile_angle - self.geoelectric_strike) < 45:
                    self.geoelectric_strike += 90
            elif 90 <= self.profile_angle < 135:
                if self.profile_angle - self.geoelectric_strike < 45:
                    self.geoelectric_strike -= 90
            else:
                if self.profile_angle - self.geoelectric_strike >= 135:
                    self.geoelectric_strike += 90

        self.geoelectric_strike = self.geoelectric_strike % 180

        # rotate components of Z and Tipper to align with geoelectric strike
        # which should now be perpendicular to the profile strike
        if self._rotate_to_strike == True:
            self.profile_angle = self.geoelectric_strike + 90
            p1 = np.tan(np.deg2rad(90 - self.profile_angle))
            # need to project the y-intercept to the new angle
            p2 = (self.profile_line[0] - p1) * easts[0] + self.profile_line[1]
            self.profile_line = (p1, p2)

            for edi in self.edi_list:
                edi.Z.rotate(self.geoelectric_strike - edi.Z.rotation_angle)
                # rotate tipper to profile azimuth, not strike.
                try:
                    edi.Tipper.rotate(
                        (self.profile_angle - 90) % 180
                        - edi.Tipper.rotation_angle.mean()
                    )
                except AttributeError:
                    edi.Tipper.rotate(
                        (self.profile_angle - 90) % 180
                        - edi.Tipper.rotation_angle
                    )

            print("=" * 72)
            print(
                (
                    "Rotated Z and Tipper to align with "
                    "{0:+.2f} degrees E of N".format(self.geoelectric_strike)
                )
            )
            print(
                (
                    "Profile angle is "
                    "{0:+.2f} degrees E of N".format(self.profile_angle)
                )
            )
            print("=" * 72)
        else:
            for edi in self.edi_list:
                edi.Z.rotate(
                    (self.profile_angle - 90) % 180 - edi.Z.rotation_angle
                )
                # rotate tipper to profile azimuth, not strike.
                try:
                    edi.Tipper.rotate(
                        (self.profile_angle - 90) % 180
                        - edi.Tipper.rotation_angle.mean()
                    )
                except AttributeError:
                    edi.Tipper.rotate(
                        (self.profile_angle - 90) % 180
                        - edi.Tipper.rotation_angle
                    )

            print("=" * 72)
            print(
                (
                    "Rotated Z and Tipper to be perpendicular  with "
                    "{0:+.2f} profile angle".format(
                        (self.profile_angle - 90) % 180
                    )
                )
            )
            print(
                (
                    "Profile angle is "
                    "{0:+.2f} degrees E of N".format(self.profile_angle)
                )
            )
            print("=" * 72)

        # --> project stations onto profile line
        projected_stations = np.zeros((self.num_edi, 2))
        self.station_locations = np.zeros(self.num_edi)

        # create profile vector
        profile_vector = np.array([1, self.profile_line[0]])
        # be sure the amplitude is 1 for a unit vector
        profile_vector /= np.linalg.norm(profile_vector)

        for ii, edi in enumerate(self.edi_list):
            station_vector = np.array(
                [easts[ii], norths[ii] - self.profile_line[1]]
            )
            position = np.dot(profile_vector, station_vector) * profile_vector
            self.station_locations[ii] = np.linalg.norm(position)
            edi.offset = np.linalg.norm(position)
            edi.projected_east = position[0]
            edi.projected_north = position[1] + self.profile_line[1]
            projected_stations[ii] = [
                position[0],
                position[1] + self.profile_line[1],
            ]

        # set the first station to 0
        for edi in self.edi_list:
            edi.offset -= self.station_locations.min()
        self.station_locations -= self.station_locations.min()

        # Sort from West to East:
        index_sort = np.argsort(self.station_locations)
        if self.profile_angle == 0:
            # Exception: sort from North to South
            index_sort = np.argsort(norths)

        # sorting along the profile
        self.edi_list = [self.edi_list[ii] for ii in index_sort]
        self.station_locations = np.array(
            [self.station_locations[ii] for ii in index_sort]
        )

        if self.estimate_elevation == True:
            self.project_elevation()

        self._profile_generated = True

    def project_elevation(self, elevation_model=None):
        """
        projects elevation data into the profile

        Arguments:
        -------------
            **elevation_model** : np.ndarray(3, num_elevation_points)
                                  (east, north, elevation)
                                  for now needs to be in utm coordinates
                                  if None then elevation is taken from edi_list

        Returns:
        ----------
            **elevation_profile** :
        """
        self.elevation_model = elevation_model

        # --> get an elevation model for the mesh
        if self.elevation_model == None:
            self.elevation_profile = np.zeros((2, len(self.edi_list)))
            self.elevation_profile[0, :] = np.array(
                [ss for ss in self.station_locations]
            )
            self.elevation_profile[1, :] = np.array(
                [edi.elevation for edi in self.edi_list]
            )

        # --> project known elevations onto the profile line
        else:
            self.elevation_profile = np.zeros(
                (2, self.elevation_model.shape[1])
            )
            # create profile vector
            profile_vector = np.array([1, self.profile_line[0]])
            # be sure the amplitude is 1 for a unit vector
            profile_vector /= np.linalg.norm(profile_vector)
            for ii in range(self.elevation_model.shape[1]):
                east = self.elevation_model[0, ii]
                north = self.elevation_model[1, ii]
                elev = self.elevation_model[2, ii]
                elev_vector = np.array([east, north - self.profile_line[1]])
                position = np.dot(profile_vector, elev_vector) * profile_vector
                self.elevation_profile[0, ii] = np.linalg.norm(position)
                self.elevation_profile[1, ii] = elev

    def plot_profile(self, **kwargs):
        """
        Plot the projected profile line along with original station locations
        to make sure the line projected is correct.

        ===================== =================================================
        Key Words             Description
        ===================== =================================================
        fig_dpi               dots-per-inch resolution of figure
                              *default* is 300
        fig_num               number if figure instance
                              *default* is 'Projected Profile'
        fig_size              size of figure in inches (width, height)
                              *default* is [5, 5]
        fs                    [ float ] font size in points of axes tick labels
                              axes labels are fs+2
                              *default* is 6
        lc                    [ string | (r, g, b) ]color of profile line
                              (see matplotlib.line for options)
                              *default* is 'b' -- blue
        lw                    float, width of profile line in points
                              *default* is 1
        marker                [ string ] marker for stations
                              (see matplotlib.pyplot.plot) for options
        mc                    [ string | (r, g, b) ] color of projected
                              stations.  *default* is 'k' -- black
        ms                    [ float ] size of station marker
                              *default* is 5
        station_id            [min, max] index values for station labels
                              *default* is None
        ===================== =================================================

        :Example: ::
            >>> edipath = r"/home/mt/edi_files"
            >>> pr = occam2d.Profile(edi_path=edipath)
            >>> pr.generate_profile()
            >>> # set station labels to only be from 1st to 4th index
            >>> # of station name
            >>> pr.plot_profile(station_id=[0,4])

        """

        fig_num = kwargs.pop("fig_num", "Projected Profile")
        fig_size = kwargs.pop("fig_size", [5, 5])
        fig_dpi = kwargs.pop("fig_dpi", 300)
        marker = kwargs.pop("marker", "v")
        ms = kwargs.pop("ms", 5)
        mc = kwargs.pop("mc", "k")
        lc = kwargs.pop("lc", "b")
        lw = kwargs.pop("ls", 1)
        fs = kwargs.pop("fs", 6)
        station_id = kwargs.pop("station_id", None)

        plt.rcParams["figure.subplot.left"] = 0.12
        plt.rcParams["figure.subplot.right"] = 0.98
        plt.rcParams["font.size"] = fs

        if self._profile_generated is False:
            self.generate_profile()

        fig = plt.figure(fig_num, figsize=fig_size, dpi=fig_dpi)
        ax = fig.add_subplot(1, 1, 1, aspect="equal")

        for edi in self.edi_list:
            (m1,) = ax.plot(
                edi.projected_east,
                edi.projected_north,
                marker=marker,
                ms=ms,
                mfc=mc,
                mec=mc,
                color=lc,
            )

            (m2,) = ax.plot(
                edi.east,
                edi.north,
                marker=marker,
                ms=0.5 * ms,
                mfc=(0.6, 0.6, 0.6),
                mec=(0.6, 0.6, 0.6),
                color=lc,
            )

            if station_id is None:
                ax.text(
                    edi.projected_east,
                    edi.projected_north * 1.00025,
                    edi.station,
                    ha="center",
                    va="baseline",
                    fontdict={"size": fs, "weight": "bold"},
                )
            else:
                ax.text(
                    edi.projected_east,
                    edi.projected_north * 1.00025,
                    edi.station[station_id[0] : station_id[1]],
                    ha="center",
                    va="baseline",
                    fontdict={"size": fs, "weight": "bold"},
                )

        peasts = np.array([edi.projected_east for edi in self.edi_list])
        pnorths = np.array([edi.projected_north for edi in self.edi_list])
        easts = np.array([edi.east for edi in self.edi_list])
        norths = np.array([edi.north for edi in self.edi_list])

        ploty = sp.polyval(self.profile_line, easts)
        ax.plot(easts, ploty, lw=lw, color=lc)
        ax.set_title("Original/Projected Stations")
        ax.set_ylim(
            (
                min([norths.min(), pnorths.min()]) * 0.999,
                max([norths.max(), pnorths.max()]) * 1.001,
            )
        )
        ax.set_xlim(
            (
                min([easts.min(), peasts.min()]) * 0.98,
                max([easts.max(), peasts.max()]) * 1.02,
            )
        )
        ax.set_xlabel(
            "Easting (m)", fontdict={"size": fs + 2, "weight": "bold"}
        )
        ax.set_ylabel(
            "Northing (m)", fontdict={"size": fs + 2, "weight": "bold"}
        )
        ax.grid(alpha=0.5)
        ax.legend(
            [m1, m2],
            ["Projected", "Original"],
            loc="upper left",
            prop={"size": fs},
        )

        plt.show()
