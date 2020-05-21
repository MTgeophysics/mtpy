# -*- coding: utf-8 -*-
"""
Create shape files for phase tensor ellipses.
https://pcjericks.github.io/py-gdalogr-cookbook/vector_layers.html#create-a-new-shapefile-and-add-data

Created on Sun Apr 13 12:32:16 2014

@author: jrpeacock
"""
import mtpy.modeling.modem
from mtpy.utils.gis_tools import project_point_ll2utm, get_utm_zone

try:
    from osgeo import ogr, gdal, osr
except ImportError:
    raise ImportError('Did not find GDAL, be sure it is installed correctly and '
          'all the paths are correct')
import numpy as np
import os
import mtpy.core.mt as mt
import mtpy.modeling.modem as modem
import mtpy.analysis.pt as mtpt
import click

ogr.UseExceptions()

class PTShapeFile(object):
    """
    write shape file for GIS plotting programs

    ======================== ==================================================
    key words/attributes      Description
    ======================== ==================================================
    edi_list                 list of edi files, full paths
    ellipse_size             size of normalized ellipse in map scale
                             *default* is .01
    mt_obj_list              list of mt.MT objects
                             *default* is None, filled if edi_list is given
    plot_period              list or value of period to convert to shape file
                             *default* is None, which will write a file for
                             every period in the edi files
    ptol                     tolerance to look for given periods
                             *default* is .05
    pt_dict                  dictionary with keys of plot_period.  Each
                             dictionary key is a structured array containing
                             the important information for the phase tensor.
    projection               projection of coordinates see EPSG for all options
                             *default* is WSG84 in lat and lon
    save_path                path to save files to
                             *default* is current working directory.
    ======================== ==================================================


    ======================== ==================================================
    Methods                   Description
    ======================== ==================================================
    _get_plot_period         get a list of all frequencies possible from
                             input files
    _get_pt_array            get phase tensors from input files and put the
                             information into a structured array
    write_shape_files        write shape files based on attributes of class
    ======================== ==================================================

    * This will project the data into UTM WSG84
    :Example: ::
        >>> edipath = r"/home/edi_files_rotated_to_geographic_north"
        >>> edilist = [os.path.join(edipath, edi) \
                      for edi in os.listdir(edipath)\
                      if edi.find('.edi')>0]
        >>> pts = PTShapeFile(edilist, save_path=r"/home/gis")
        >>> pts.write_shape_files()

    * To project into another datum, set the projection attribute
    :Example: ::
        >>> pts = PTShapeFile(edilist, save_path=r"/home/gis")
        >>> pts.projection = 'NAD27'
        >>> pts.write_shape_files()


    """

    def __init__(self, edi_list=None, proj='WGS84', esize=0.03, **kwargs):
        self.edi_list = edi_list
        self.projection = proj
        #self.projection = None
        # UTM zone 50 {'init': u'epsg:32750'} UTM zone 51 32751
        # WGS84: 'epsg:4326'  GDA94:  EPSG:4283 See  http://epsg.io/4283​
        # http://spatialreference.org/ref/epsg/4283/
        # self.utm_cs = None

        self.plot_period = None
        self.save_path = None  # os.getcwd()
        # self.ellipse_size = 500.0  # maximum ellipse major axis size in
        # metres
        self.ellipse_size = esize  # 0.002  # maximum ellipse major axis size in metres
        #self._theta = np.arange(0, 2 * np.pi, np.pi / 180.)
        self._theta = np.arange(0, 2 * np.pi, np.pi / 30.)  # FZ: adjusted number of points in array
        self.ptol = .05  # period value tolerance to be considered as equal

        self.mt_obj_list = None
        self.pt_dict = None

        if self.edi_list is not None:
            self.mt_obj_list = [mt.MT(edi) for edi in self.edi_list]
        # else:
        #     raise Exception("EDI files List is None")

        for key in list(kwargs.keys()):
            setattr(self, key, kwargs[key])

        if self.mt_obj_list is not None:
            self._get_plot_period()
            self._get_pt_array()

            print((self.plot_period))

        self._proj_dict = {'WGS84': 4326, 'NAD27': 4267, 'GDA94': 4283}

        self._rotation_angle = 0.0

    def _set_rotation_angle(self, rotation_angle):
        """
        rotate all mt_objs to rotation angle
        """

        self._rotation_angle = float(rotation_angle)

        for mt_obj in self.mt_obj_list:
            mt_obj.rotation_angle = float(self._rotation_angle)

    def _get_rotation_angle(self):
        return self._rotation_angle

    rotation_angle = property(_get_rotation_angle, _set_rotation_angle,
                              doc="rotation angle of Z and Tipper")

    def _get_plot_period(self):
        """
        from the list of edi's get a frequency list from all possible
        frequencies.

        """

        if self.plot_period is None:
            # get all frequencies from all edi files
            all_freqs = []
            for mt_obj in self.mt_obj_list:
                all_freqs.extend(list(mt_obj.Z.freq))

            # sort all frequencies so that they are in descending order,
            # use set to remove repeats and make an array
            self.plot_period = 1. / np.array(sorted(list(set(all_freqs)),
                                                    reverse=True))
        else:
            if isinstance(self.plot_period, list):
                pass
            if isinstance(self.plot_period, int) or isinstance(
                    self.plot_period, float):
                self.plot_period = [self.plot_period]

    def _get_pt_array(self, periods=None):
        """
        get the phase tensor information into a form that is more structured
        to manipulate easier later.

        make a dictionary with keys being the plot period values and each
        key has a structured array that contains all the important information
        collected from each station.

        Note: this only supports 2 coordinate system 1) lat-long; 2) UTM-WGS84 default

        :param periods: a list of the periods as a subset of the all periods.
        :return:
        """

        utm_cs_list=[]  # used to detect multi-UTM zones

        self.pt_dict = {}
        if self.plot_period is None:
            self._get_plot_period()

        if periods is None:
            periods=self.plot_period

        for plot_per in periods:
            self.pt_dict[plot_per] = []

            for mt_obj in self.mt_obj_list:
                p_index = [ff for ff, f2 in enumerate(1. / mt_obj.Z.freq)
                           if (f2 > plot_per * (1 - self.ptol)) and
                           (f2 < plot_per * (1 + self.ptol))]

                if len(p_index) >= 1:
                    p_index = p_index[0]

                    if self.projection is None:  # geographic-coord lat lon
                        east, north, elev = (mt_obj.lon, mt_obj.lat, 0)
                        self.utm_cs = osr.SpatialReference()
                        # Set geographic coordinate system to handle lat/lon
                        # self.utm_cs.SetWellKnownGeogCS(self.projection)
                        self.utm_cs.ImportFromEPSG(4326)
                        # create the spatial reference, WGS84=4326
                        # GDA94 = EPSG:4283 See  http://epsg.io/4283
                    elif self.projection == 'WGS84':  # UTM zones coordinate system
                        edi_proj = 'WGS84'
                        east, north, _ = project_point_ll2utm(mt_obj.lat, mt_obj.lon, edi_proj)
                        zone_number, is_northern, _ = get_utm_zone(mt_obj.lat, mt_obj.lon)
                        self.utm_cs = osr.SpatialReference()
                        self.utm_cs.SetUTM(zone_number, is_northern)
                        utm_cs_list.append(self.utm_cs.GetAttrValue('projcs'))
                    else:
                        raise Exception(
                            "%s is NOT supported" %
                            self.projection)

                    pt_tuple = (mt_obj.station, east, north,
                                mt_obj.pt.phimin[p_index],
                                mt_obj.pt.phimax[p_index],
                                mt_obj.pt.azimuth[p_index],
                                mt_obj.pt.beta[p_index],
                                2 * mt_obj.pt.beta[p_index],
                                mt_obj.pt.ellipticity[p_index])  # FZ: get ellipticity begin here

                    self.pt_dict[plot_per].append(pt_tuple)
                else:
                    print(("The period %s is NOT found for this station %s" %(plot_per, mt_obj.station)))

            self.pt_dict[plot_per] = np.array(self.pt_dict[plot_per],
                                              dtype=[('station', '|S15'),
                                                     ('east', np.float),
                                                     ('north', np.float),
                                                     ('phimin', np.float),
                                                     ('phimax', np.float),
                                                     ('azimuth', np.float),
                                                     ('skew', np.float),
                                                     ('n_skew', np.float),
                                                     ('ellipticity', np.float)])
        unique_utm_cs = sorted(list(set(utm_cs_list)))
        if len(unique_utm_cs) >1:
            print(("Warning: Multi-UTM-Zones found in the EDI files", unique_utm_cs))


    def write_shape_files(self, periods=None):
        """
        write shape file from given attributes
        https://pcjericks.github.io/py-gdalogr-cookbook/vector_layers.html
        #create-a-new-shapefile-and-add-data
        """

        # Why call again:
        # self._get_pt_array() # already called in __init__?

        if periods is None:
            periods = self.plot_period

        #for plot_per in self.plot_period:
        for plot_per in periods:
            # shape file path
            shape_fn = os.path.join(self.save_path,
                                    'PT_{0:.5g}s_{1}.shp'.format(plot_per,
                                                                 self.projection))

            # remove the shape file if it already exists, has trouble over
            # writing
            if os.path.isfile(shape_fn) == True:
                os.remove(shape_fn)

            # need to tell ogr which driver to use
            driver = ogr.GetDriverByName('ESRI Shapefile')

            if os.path.isfile(shape_fn) == True:
                driver.DeleteDataSource(shape_fn)

            # create shape file
            data_source = driver.CreateDataSource(shape_fn)

            #            ##if you read from a raster get the georeference point otherwise create one
            #            spatial_ref = osr.SpatialReference()
            #            #this puts it in the wsg84 reference frame.
            #            spatial_ref.ImportFromEPSG(self._proj_dict[self.projection])

            # create a layer to put the ellipses onto
            layer = data_source.CreateLayer('PT', self.utm_cs, ogr.wkbPolygon)

            # make field names
            field_name = ogr.FieldDefn("Name", ogr.OFTString)
            layer.CreateField(field_name)

            field_phimin = ogr.FieldDefn('phi_min', ogr.OFTReal)
            layer.CreateField(field_phimin)

            field_phimax = ogr.FieldDefn('phi_max', ogr.OFTReal)
            layer.CreateField(field_phimax)

            field_skew = ogr.FieldDefn('skew', ogr.OFTReal)
            layer.CreateField(field_skew)

            field_normalized_skew = ogr.FieldDefn('n_skew', ogr.OFTReal)
            layer.CreateField(field_normalized_skew)

            # FZ added azimuth
            field_azimuth = ogr.FieldDefn('azimuth', ogr.OFTReal)
            layer.CreateField(field_azimuth)

            field_ellipticity = ogr.FieldDefn('ellipt', ogr.OFTReal)
            # FZ: note osgeo gdal does not like name 'ellipticity'
            layer.CreateField(field_ellipticity)

            poly_list = []

#            print("period=", plot_per)
#            print(self.pt_dict.keys()) # (self.pt_dict[plot_per])['phimax'].size)
            phi_max_val = self.pt_dict[plot_per]['phimax'].max()


            for isite, pt_array in enumerate(self.pt_dict[plot_per]):

                # need to make an ellipse first using the parametric
                # equation
                azimuth = -np.deg2rad(pt_array['azimuth'])
                width = self.ellipse_size * (pt_array['phimax'] / phi_max_val)
                height = self.ellipse_size * (pt_array['phimin'] / phi_max_val)

                x0 = pt_array['east']
                y0 = pt_array['north']


                # apply formula to generate ellipses
                x = x0 + height * np.cos(self._theta) * np.cos(azimuth) - \
                    width * np.sin(self._theta) * np.sin(azimuth)
                y = y0 + height * np.cos(self._theta) * np.sin(azimuth) + \
                    width * np.sin(self._theta) * np.cos(azimuth)

                # 1) make a geometry shape of the ellipse
                ellipse = ogr.Geometry(ogr.wkbLinearRing)

                for ii, jj in zip(x, y):
                    ellipse.AddPoint(np.round(ii, 6), np.round(jj, 6))

                ellipse.CloseRings()

                # 2) make a polygon
                poly = ogr.Geometry(ogr.wkbPolygon)
                poly.AddGeometry(ellipse)

                poly_list.append(poly)

                # 4) this part is confusing but we need to create a feature that has the
                # same definition as the layer that we created.
                # get the layer definition
                feature_def = layer.GetLayerDefn()

                # create a new feature
                new_feature = ogr.Feature(feature_def)
                # set the geometry of that feature to be the ellipse
                new_feature.SetGeometry(poly)
                # create the feature in the layer.
                layer.CreateFeature(new_feature)

                #
                # 5) create a field to color by
                val_sta = pt_array['station']
                print(" the type of pt_array['station'] = ", type(val_sta))
                # decode the string cal_sta (numpy_bytes_ )  back into ASCII
                new_feature.SetField("Name", pt_array['station'].decode('UTF-8'))
                new_feature.SetField("phi_min", pt_array['phimin'])
                new_feature.SetField("phi_max", pt_array['phimax'])
                new_feature.SetField("skew", pt_array['skew'])
                new_feature.SetField("n_skew", pt_array['n_skew'])

                new_feature.SetField(
                    "azimuth", pt_array['azimuth'])  # FZ added
                new_feature.SetField(
                    "ellipt", pt_array['ellipticity'])  # FZ added
                # new_feature.SetField("ellipticity", pt_array['azimuth'])
                # # FZ added

                # add the new feature to the layer.
                layer.SetFeature(new_feature)

                # apparently need to destroy the feature
                new_feature.Destroy()


            # Need to be sure that all the new info is saved to
            data_source.SyncToDisk()

            # write a projection file
            #            spatial_ref.MorphToESRI()
            self.utm_cs.MorphToESRI()
            prj_file = open('{0}prj'.format(shape_fn[:-3]), 'w')
            #            prj_file.write(spatial_ref.ExportToWkt())
            prj_file.write(self.utm_cs.ExportToWkt())
            prj_file.close()

            data_source.Destroy()

            print('Wrote shape file to {0}'.format(shape_fn))

# ===========================
    def write_data_pt_shape_files_modem(self, modem_data_fn,
                                        rotation_angle=0.0):
        """
        write pt files from a modem data file.

        """

        modem_obj = mtpy.modeling.modem.Data()
        modem_obj.read_data_file(modem_data_fn)

        self.plot_period = modem_obj.period_list.copy()
        self.mt_obj_list = [modem_obj.mt_dict[key]
                            for key in list(modem_obj.mt_dict.keys())]

        self._get_pt_array()
        self._set_rotation_angle(rotation_angle)

        self.write_shape_files()

    def write_resp_pt_shape_files_modem(self, modem_data_fn, modem_resp_fn,
                                        rotation_angle=0.0):
        """
        write pt files from a modem response file where ellipses are normalized
        by the data file.

        """

        # first get the data and response and place them in array for later use
        modem_data_obj = mtpy.modeling.modem.Data()
        modem_data_obj.read_data_file(modem_data_fn)

        self.plot_period = modem_data_obj.period_list.copy()
        self.mt_obj_list = [modem_data_obj.mt_dict[key]
                            for key in list(modem_data_obj.mt_dict.keys())]
        self._get_pt_array()

        self._set_rotation_angle(rotation_angle)

        modem_resp_obj = mtpy.modeling.modem.Data()
        modem_resp_obj.read_data_file(modem_resp_fn)

        # rotate model response
        for r_key in list(modem_resp_obj.mt_dict.keys()):
            modem_resp_obj.mt_dict[
                r_key].rotation_angle = float(rotation_angle)

        resp_pt_dict = {}
        for p_index, plot_per in enumerate(self.plot_period):
            resp_pt_dict[plot_per] = []
            for key in list(modem_data_obj.mt_dict.keys()):
                mt_obj = modem_data_obj.mt_dict[key]
                if self.projection is None:
                    east, north, elev = (mt_obj.lon, mt_obj.lat, 0)
                    self.utm_cs = osr.SpatialReference()
                    # Set geographic coordinate system to handle lat/lon
                    self.utm_cs.SetWellKnownGeogCS(self.projection)
                else:
                    self.utm_cs, utm_point = project_point_ll2utm(mt_obj.lon,
                                                                 mt_obj.lat,
                                                                 self.projection)
                    east, north, elev = utm_point

                # get pt objects from data and model response
                try:
                    mpt = modem_resp_obj.mt_dict[key].pt

                    pt_tuple = (mt_obj.station, east, north,
                                mpt.phimin[p_index],
                                mpt.phimax[p_index],
                                mpt.azimuth[p_index],
                                mpt.beta[p_index],
                                2 * mpt.beta[p_index])
                except KeyError:
                    pt_tuple = (mt_obj.station, east, north,
                                0,
                                0,
                                0,
                                0,
                                0)
                resp_pt_dict[plot_per].append(pt_tuple)

            # now make each period an array for writing to file
            resp_pt_dict[plot_per] = np.array(resp_pt_dict[plot_per],
                                              dtype=[('station', '|S15'),
                                                     ('east', np.float),
                                                     ('north', np.float),
                                                     ('phimin', np.float),
                                                     ('phimax', np.float),
                                                     ('azimuth', np.float),
                                                     ('skew', np.float),
                                                     ('n_skew', np.float)])

        # write files
        for plot_per in self.plot_period:
            # shape file path
            shape_fn = os.path.join(self.save_path,
                                    'Resp_PT_{0:.5g}s_{1}.shp'.format(plot_per,
                                                                      self.projection))

            # remove the shape file if it already exists, has trouble over
            # writing
            if os.path.isfile(shape_fn) == True:
                os.remove(shape_fn)

            # need to tell ogr which driver to use
            driver = ogr.GetDriverByName('ESRI Shapefile')

            if os.path.isfile(shape_fn) == True:
                driver.DeleteDataSource(shape_fn)

            # create shape file
            data_source = driver.CreateDataSource(shape_fn)

            # create a layer to put the ellipses onto
            layer = data_source.CreateLayer('RPT', self.utm_cs, ogr.wkbPolygon)

            # make field names
            field_name = ogr.FieldDefn("Name", ogr.OFTString)
            layer.CreateField(field_name)

            field_phimin = ogr.FieldDefn('phi_min', ogr.OFTReal)
            layer.CreateField(field_phimin)

            field_phimax = ogr.FieldDefn('phi_max', ogr.OFTReal)
            layer.CreateField(field_phimax)

            field_skew = ogr.FieldDefn('skew', ogr.OFTReal)
            layer.CreateField(field_skew)

            field_normalized_skew = ogr.FieldDefn('n_skew', ogr.OFTReal)
            layer.CreateField(field_normalized_skew)

            poly_list = []
            phimax = self.pt_dict[plot_per]['phimax'].max()
            for pt_array in resp_pt_dict[plot_per]:

                # need to make an ellipse first using the parametric equation
                azimuth = -np.deg2rad(pt_array['azimuth'])
                width = self.ellipse_size * (pt_array['phimax'] / phimax)
                height = self.ellipse_size * (pt_array['phimin'] / phimax)
                x0 = pt_array['east']
                y0 = pt_array['north']

                x = x0 + height * np.cos(self._theta) * np.cos(azimuth) - \
                    width * np.sin(self._theta) * np.sin(azimuth)
                y = y0 + height * np.cos(self._theta) * np.sin(azimuth) + \
                    width * np.sin(self._theta) * np.cos(azimuth)

                # 1) make a geometry shape of the ellipse
                ellipse = ogr.Geometry(ogr.wkbLinearRing)

                for ii, jj in zip(x, y):
                    ellipse.AddPoint(np.round(ii, 6), np.round(jj, 6))

                ellipse.CloseRings()

                # 2) make a polygon
                poly = ogr.Geometry(ogr.wkbPolygon)
                poly.AddGeometry(ellipse)

                poly_list.append(poly)

                # 4) this part is confusing but we need to create a feature that has the
                # same definition as the layer that we created.
                # get the layer definition
                feature_def = layer.GetLayerDefn()

                # create a new feature
                new_feature = ogr.Feature(feature_def)
                # set the geometry of that feature to be the ellipse
                new_feature.SetGeometry(poly)
                # create the feature in the layer.
                layer.CreateFeature(new_feature)

                #
                # 5) create a field to color by
                new_feature.SetField("Name", pt_array['station'])
                new_feature.SetField("phi_min", pt_array['phimin'])
                new_feature.SetField("phi_max", pt_array['phimax'])
                new_feature.SetField("skew", pt_array['skew'])
                new_feature.SetField("n_skew", pt_array['n_skew'])

                # add the new feature to the layer.
                layer.SetFeature(new_feature)

                # apparently need to destroy the feature
                new_feature.Destroy()

            # Need to be sure that all the new info is saved to
            data_source.SyncToDisk()

            # write a projection file
            self.utm_cs.MorphToESRI()
            prj_file = open('{0}prj'.format(shape_fn[:-3]), 'w')
            prj_file.write(self.utm_cs.ExportToWkt())
            prj_file.close()

            data_source.Destroy()

            print('Wrote shape file to {0}'.format(shape_fn))

    def write_residual_pt_shape_files_modem(self, modem_data_fn, modem_resp_fn,
                                            rotation_angle=0.0, normalize='1'):
        """
        write residual pt shape files from ModEM output

        normalize [ '1' | 'all' ]
                   * '1' to normalize the ellipse by itself, all ellipses are
                         normalized to phimax, thus one axis is of length
                         1*ellipse_size
                   * 'all' to normalize each period by the largest phimax

        """

        # first get the data and response and place them in array for later use
        modem_data_obj = mtpy.modeling.modem.Data()
        modem_data_obj.read_data_file(modem_data_fn)

        self.plot_period = modem_data_obj.period_list.copy()
        self.mt_obj_list = [modem_data_obj.mt_dict[key]
                            for key in list(modem_data_obj.mt_dict.keys())]
        self._get_pt_array()

        self._set_rotation_angle(rotation_angle)

        modem_resp_obj = mtpy.modeling.modem.Data()
        modem_resp_obj.read_data_file(modem_resp_fn)

        # rotate model response
        for r_key in list(modem_resp_obj.mt_dict.keys()):
            modem_resp_obj.mt_dict[
                r_key].rotation_angle = float(rotation_angle)

        residual_pt_dict = {}
        for p_index, plot_per in enumerate(self.plot_period):
            residual_pt_dict[plot_per] = []
            for key in list(modem_data_obj.mt_dict.keys()):
                mt_obj = modem_data_obj.mt_dict[key]
                if self.projection is None:
                    east, north, elev = (mt_obj.lon, mt_obj.lat, 0)
                    self.utm_cs = osr.SpatialReference()
                    # Set geographic coordinate system to handle lat/lon
                    self.utm_cs.SetWellKnownGeogCS(self.projection)
                else:
                    self.utm_cs, utm_point = project_point_ll2utm(mt_obj.lon,
                                                                 mt_obj.lat,
                                                                 self.projection)
                    east, north, elev = utm_point

                # get pt objects from data and model response
                try:
                    dpt = modem_data_obj.mt_dict[key].pt
                    mpt = modem_resp_obj.mt_dict[key].pt
                except KeyError:
                    print('No information found for {0} in {1}').format(key, modem_resp_fn)
                    continue

                # calculate the residual pt
                try:
                    rpt = mtpt.ResidualPhaseTensor(pt_object1=dpt,
                                                   pt_object2=mpt)
                    rpt = rpt.residual_pt
                    rpt_mean = .25 * np.linalg.norm(rpt.pt[p_index], ord='fro')
                    #                    rpt_mean = .25*np.sqrt(abs(rpt.pt[p_index, 0, 0])**2+
                    #                                          abs(rpt.pt[p_index, 0, 1])**2+
                    #                                          abs(rpt.pt[p_index, 1, 0])**2+
                    # abs(rpt.pt[p_index, 1, 1])**2)
                    pt_tuple = (mt_obj.station, east, north,
                                rpt.phimin[p_index],
                                rpt.phimax[p_index],
                                rpt.azimuth[p_index],
                                rpt.beta[p_index],
                                rpt_mean)
                    #                                np.sqrt(abs(rpt.phimin[0][p_index]*
                    #                                            rpt.phimax[0][p_index])))
                    residual_pt_dict[plot_per].append(pt_tuple)
                except mtpt.MTex.MTpyError_PT:
                    print(key, dpt.pt.shape, mpt.pt.shape)
                    pt_tuple = (mt_obj.station, east, north,
                                0.0,
                                0.0,
                                0.0,
                                0.0,
                                0.0)
                    residual_pt_dict[plot_per].append(pt_tuple)
            # now make each period an array for writing to file
            residual_pt_dict[plot_per] = np.array(residual_pt_dict[plot_per],
                                                  dtype=[('station', '|S15'),
                                                         ('east', np.float),
                                                         ('north', np.float),
                                                         ('phimin', np.float),
                                                         ('phimax', np.float),
                                                         ('azimuth', np.float),
                                                         ('skew', np.float),
                                                         ('geometric_mean', np.float)])
        # return residual_pt_dict
        # write files
        for plot_per in self.plot_period:
            # shape file path
            shape_fn = os.path.join(self.save_path,
                                    'ResidualPT_{0:.5g}s_{1}.shp'.format(plot_per,
                                                                         self.projection))

            # remove the shape file if it already exists, has trouble over
            # writing
            if os.path.isfile(shape_fn) == True:
                os.remove(shape_fn)

            # need to tell ogr which driver to use
            driver = ogr.GetDriverByName('ESRI Shapefile')

            if os.path.isfile(shape_fn) == True:
                driver.DeleteDataSource(shape_fn)

            # create shape file
            data_source = driver.CreateDataSource(shape_fn)

            # create a layer to put the ellipses onto
            layer = data_source.CreateLayer('RPT', self.utm_cs, ogr.wkbPolygon)

            # make field names
            field_name = ogr.FieldDefn("Name", ogr.OFTString)
            layer.CreateField(field_name)

            field_phimin = ogr.FieldDefn('phi_min', ogr.OFTReal)
            layer.CreateField(field_phimin)

            field_phimax = ogr.FieldDefn('phi_max', ogr.OFTReal)
            layer.CreateField(field_phimax)

            field_skew = ogr.FieldDefn('skew', ogr.OFTReal)
            layer.CreateField(field_skew)

            field_geometric_mean = ogr.FieldDefn('mean', ogr.OFTReal)
            layer.CreateField(field_geometric_mean)

            poly_list = []
            phimax = self.pt_dict[plot_per]['phimax'].max()
            for pt_array in residual_pt_dict[plot_per]:
                # need to make an ellipse first using the parametric equation
                azimuth = -np.deg2rad(pt_array['azimuth'])
                if normalize == '1':
                    width = self.ellipse_size * \
                        (pt_array['phimax'] / pt_array['phimax'])
                    height = self.ellipse_size * \
                        (pt_array['phimin'] / pt_array['phimax'])
                elif normalize == 'all':
                    width = self.ellipse_size * (pt_array['phimax'] / phimax)
                    height = self.ellipse_size * (pt_array['phimin'] / phimax)
                x0 = pt_array['east']
                y0 = pt_array['north']

                x = x0 + height * np.cos(self._theta) * np.cos(azimuth) - \
                    width * np.sin(self._theta) * np.sin(azimuth)
                y = y0 + height * np.cos(self._theta) * np.sin(azimuth) + \
                    width * np.sin(self._theta) * np.cos(azimuth)

                # 1) make a geometry shape of the ellipse
                ellipse = ogr.Geometry(ogr.wkbLinearRing)

                for ii, jj in zip(x, y):
                    ellipse.AddPoint(np.round(ii, 6), np.round(jj, 6))

                ellipse.CloseRings()

                # 2) make a polygon
                poly = ogr.Geometry(ogr.wkbPolygon)
                poly.AddGeometry(ellipse)

                poly_list.append(poly)

                # 4) this part is confusing but we need to create a feature that has the
                # same definition as the layer that we created.
                # get the layer definition
                feature_def = layer.GetLayerDefn()

                # create a new feature
                new_feature = ogr.Feature(feature_def)
                # set the geometry of that feature to be the ellipse
                new_feature.SetGeometry(poly)
                # create the feature in the layer.
                layer.CreateFeature(new_feature)

                #
                # 5) create a field to color by
                new_feature.SetField('Name', pt_array['station'])
                new_feature.SetField('phi_min', pt_array['phimin'])
                new_feature.SetField('phi_max', pt_array['phimax'])
                new_feature.SetField('skew', pt_array['skew'])
                new_feature.SetField('mean', pt_array['geometric_mean'])

                # add the new feature to the layer.
                layer.SetFeature(new_feature)

                # apparently need to destroy the feature
                new_feature.Destroy()

            # Need to be sure that all the new info is saved to
            data_source.SyncToDisk()

            # write a projection file
            self.utm_cs.MorphToESRI()
            prj_file = open('{0}prj'.format(shape_fn[:-3]), 'w')
            prj_file.write(self.utm_cs.ExportToWkt())
            prj_file.close()

            data_source.Destroy()

            print('Wrote shape file to {0}'.format(shape_fn))


# ==============================================================================
# Tipper arrows
# ==============================================================================
class TipperShapeFile(object):
    """
    write shape file for GIS plotting programs.

    currently only writes the real induction vectors.

    ======================== ==================================================
    key words/attributes      Description
    ======================== ==================================================
    arrow_direction          [ 1 | -1 ] 1 for Weise convention --> point
                             toward conductors. *default* is 1
                             (-1 is not supported yet)
    arrow_head_height        height of arrow head in map units
                             *default* is .002
    arrow_head_width         width of arrow head in map units
                             *default* is .001
    arrow_lw                 width of arrow in map units
                             *default* is .0005

    arrow_size               size of normalized arrow length in map units
                             *default* is .01

    edi_list                 list of edi files, full paths
    mt_obj_list              list of mt.MT objects
                             *default* is None, filled if edi_list is given
    plot_period              list or value of period to convert to shape file
                             *default* is None, which will write a file for
                             every period in the edi files
    ptol                     tolerance to look for given periods
                             *default* is .05
    pt_dict                  dictionary with keys of plot_period.  Each
                             dictionary key is a structured array containing
                             the important information for the phase tensor.
    projection               projection of coordinates see EPSG for all options
                             *default* is WSG84
    save_path                path to save files to
                             *default* is current working directory.
    ======================== ==================================================

    ======================== ==================================================
    Methods                   Description
    ======================== ==================================================
    _get_plot_period         get a list of all possible frequencies from data
    _get_tip_array           get Tipper information from data and put into
                             a structured array for easy manipulation
    write_real_shape_files   write real induction arrow shape files
    write_imag_shape_files   write imaginary induction arrow shape files
    ======================== ==================================================

    :Example: ::
        >>> edipath = r"/home/edi_files_rotated_to_geographic_north"
        >>> edilist = [os.path.join(edipath, edi) \
                      for edi in os.listdir(edipath)\
                      if edi.find('.edi')>0]
        >>> tipshp = TipperShapeFile(edilist, save_path=r"/home/gis")
        >>> tipshp.arrow_head_height = .005
        >>> tipshp.arrow_lw = .0001
        >>> tipshp.arrow_size = .05
        >>> tipshp.write_shape_files()

    """

    def __init__(self, edi_list=None, **kwargs):
        self.edi_list = edi_list
        self.projection = 'WGS84'
        self.plot_period = None
        self.save_path = os.getcwd()
        self.arrow_size = 1000
        self.arrow_direction = 1
        self.ptol = .05
        self.arrow_head_width = 50
        self.arrow_head_height = 100
        self.arrow_lw = 20
        self.utm_cs = None

        self.mt_obj_list = None
        self.tip_dict = None

        if self.edi_list is not None:
            self.mt_obj_list = [mt.MT(edi) for edi in self.edi_list]

        for key in list(kwargs.keys()):
            setattr(self, key, kwargs[key])

        if self.mt_obj_list is not None:
            self._get_plot_period()
            self._get_tip_array()

        self._proj_dict = {'WGS84': 4326, 'NAD27': 4267}
        self._rotation_angle = 0.0

    def _set_rotation_angle(self, rotation_angle):
        """
        rotate all mt_objs to rotation angle
        """

        self._rotation_angle = float(rotation_angle)

        for mt_obj in self.mt_obj_list:
            mt_obj.rotation_angle = float(self._rotation_angle)

    def _get_rotation_angle(self):
        return self._rotation_angle

    rotation_angle = property(_get_rotation_angle, _set_rotation_angle,
                              doc="rotation angle of Z and Tipper")

    def _get_plot_period(self):
        """
        from the list of edi's get a frequency list to invert for.

        """

        if self.plot_period is None:
            # get all frequencies from all edi files
            all_freqs = []
            for mt_obj in self.mt_obj_list:
                all_freqs.extend(list(mt_obj.Z.freq))

            # sort all frequencies so that they are in descending order,
            # use set to remove repeats and make an array
            self.plot_period = 1. / np.array(sorted(list(set(all_freqs)),
                                                    reverse=True))
        else:
            if isinstance(self.plot_period, list):
                pass
            if isinstance(self.plot_period, int) or isinstance(
                    self.plot_period, float):
                self.plot_period = [self.plot_period]

    def _get_tip_array(self):
        """
        get the phase tensor information into a form that is more structured
        to manipulate easier later.

        make a dictionary with keys being the plot period values and each
        key has a structured array that contains all the important information
        collected from each station.
        """
        utm_cs_list=[]

        self.tip_dict = {}
        for plot_per in self.plot_period:
            self.tip_dict[plot_per] = []
            for mt_obj in self.mt_obj_list:
                mt_obj.Tipper.compute_mag_direction()
                try:
                    p_index = [ff for ff, f2 in enumerate(1. / mt_obj.Z.freq)
                               if (f2 > plot_per * (1 - self.ptol)) and
                               (f2 < plot_per * (1 + self.ptol))][0]
                    if self.projection is None:
                        east, north, elev = (mt_obj.lon, mt_obj.lat, 0)
                        self.utm_cs = osr.SpatialReference()
                        # Set geographic coordinate system to handle lat/lon
                        self.utm_cs.SetWellKnownGeogCS(self.projection)
                    else:
                        self.utm_cs, utm_point = project_point_ll2utm(mt_obj.lon,
                                                                     mt_obj.lat,
                                                                     self.projection)
                        east, north, elev = utm_point

                        utm_cs_list.append(self.utm_cs.GetAttrValue('projcs'))

                    if mt_obj.Tipper.tipper is not None:

                        if mt_obj.Tipper.tipper[p_index].all() != 0.0:
                            tp_tuple = (mt_obj.station,
                                        east,
                                        north,
                                        mt_obj.Tipper.mag_real[p_index],
                                        mt_obj.Tipper.mag_imag[p_index],
                                        mt_obj.Tipper.angle_real[p_index],
                                        mt_obj.Tipper.angle_imag[p_index])
                            self.tip_dict[plot_per].append(tp_tuple)
                        else:
                            tp_tuple = (mt_obj.station,
                                        east,
                                        north,
                                        0,
                                        0,
                                        0,
                                        0)
                            self.tip_dict[plot_per].append(tp_tuple)

                except IndexError:
                    pass

            self.tip_dict[plot_per] = np.array(self.tip_dict[plot_per],
                                               dtype=[('station', '|S15'),
                                                      ('east', np.float),
                                                      ('north', np.float),
                                                      ('mag_real', np.float),
                                                      ('mag_imag', np.float),
                                                      ('ang_real', np.float),
                                                      ('ang_imag', np.float)])

        unique_utm_cs = sorted(list(set(utm_cs_list)))
        if len(unique_utm_cs) >1:
            print(("Warning: Multi-UTM-Zones found in the EDI files", unique_utm_cs))

    def write_real_shape_files(self):
        """
        write shape file from given attributes
        """

        self._get_tip_array()

        for plot_per in self.plot_period:
            # shape file path
            shape_fn = os.path.join(self.save_path,
                                    'Tip_{0:.5g}s_{1}_real.shp'.format(plot_per,
                                                                       self.projection))

            # remove the shape file if it already exists, has trouble over
            # writing
            if os.path.isfile(shape_fn) == True:
                os.remove(shape_fn)

            # need to tell ogr which driver to use
            driver = ogr.GetDriverByName('ESRI Shapefile')

            if os.path.isfile(shape_fn) == True:
                driver.DeleteDataSource(shape_fn)

            # create shape file
            data_source = driver.CreateDataSource(shape_fn)

            # if you read from a raster get the georeference point otherwise create one
            # spatial_ref = osr.SpatialReference()
            # this puts it in the wsg84 reference frame.
            # spatial_ref.ImportFromEPSG(self._proj_dict[self.projection])

            # create a layer to put the ellipses onto
            layer = data_source.CreateLayer('TIPPER', self.utm_cs,
                                            ogr.wkbPolygon)

            # make field names
            field_name = ogr.FieldDefn("Name", ogr.OFTString)
            layer.CreateField(field_name)

            field_mag_real = ogr.FieldDefn('mag_real', ogr.OFTReal)
            layer.CreateField(field_mag_real)

            field_ang_real = ogr.FieldDefn('ang_real', ogr.OFTReal)
            layer.CreateField(field_ang_real)

            for tp_arr in self.tip_dict[plot_per]:
                cos_t = np.cos(-np.deg2rad(tp_arr['ang_real']))
                sin_t = np.sin(-np.deg2rad(tp_arr['ang_real']))
                # calculate the points to make the line
                txr = 0
                tyr = tp_arr['mag_real'] * self.arrow_size

                # make an arrow by drawing an outline.  have the arrow point
                # north to start and then rotate later with the rotation
                # matrix to properly orient it.
                x0 = 0
                y0 = 0

                x1 = x0 + self.arrow_lw
                y1 = y0

                x2 = x0 + self.arrow_lw
                y2 = y0 + tyr - self.arrow_head_height

                x3 = x0 + self.arrow_lw + self.arrow_head_width
                y3 = y2

                x4 = x0 + txr
                y4 = y0 + tyr

                x7 = x0 - self.arrow_lw
                y7 = y0

                x6 = x0 - self.arrow_lw
                y6 = y0 + tyr - self.arrow_head_height

                x5 = x0 - self.arrow_lw - self.arrow_head_width
                y5 = y6

                x = np.array([x0, x1, x2, x3, x4, x5, x6, x7])
                y = np.array([y0, y1, y2, y3, y4, y5, y6, y7])

                rot_matrix = np.array([[cos_t, -sin_t], [sin_t, cos_t]])

                # rotate the arrow to be properly oriented
                xy = np.array([x, y])
                rot_xy = np.dot(rot_matrix, xy)

                # shift the arrow to be centered on the station.
                x = tp_arr['east'] + rot_xy[0]
                y = tp_arr['north'] + rot_xy[1]

                # 1) make a geometry shape line
                arrow = ogr.Geometry(ogr.wkbLinearRing)
                for ii, jj in zip(x, y):
                    arrow.AddPoint(np.round(ii, 6), np.round(jj, 6))

                arrow.CloseRings()

                poly = ogr.Geometry(ogr.wkbPolygon)
                poly.AddGeometry(arrow)
                # 4) this part is confusing but we need to create a
                # feature that has the
                # same definition as the layer that we created.
                #    get the layer definition
                feature_def = layer.GetLayerDefn()

                # create a new feature
                new_feature = ogr.Feature(feature_def)
                # set the geometry of that feature to be the ellipse
                new_feature.SetGeometry(poly)
                # create the feature in the layer.
                layer.CreateFeature(new_feature)

                #
                # 5) create a field to color by
                new_feature.SetField("Name", tp_arr['station'])
                new_feature.SetField("mag_real", tp_arr['mag_real'])
                new_feature.SetField("ang_real", tp_arr['ang_real'])

                # add the new feature to the layer.
                layer.SetFeature(new_feature)

                # apparently need to destroy the feature
                new_feature.Destroy()

            # Need to be sure that all the new info is saved to
            data_source.SyncToDisk()

            # write a projection file
            # spatial_ref.MorphFromESRI()
            self.utm_cs.MorphFromESRI()
            prj_file = open('{0}prj'.format(shape_fn[:-3]), 'w')
            prj_file.write(self.utm_cs.ExportToWkt())
            prj_file.close()

            data_source.Destroy()

            print('Wrote shape file to {0}'.format(shape_fn))

    def write_imag_shape_files(self):
        """
        write shape file from given attributes
        """

        self._get_tip_array()

        for plot_per in self.plot_period:
            # shape file path
            shape_fn = os.path.join(self.save_path,
                                    'Tip_{0:.5g}s_{1}_imag.shp'.format(plot_per,
                                                                       self.projection))

            # remove the shape file if it already exists, has trouble over
            # writing
            if os.path.isfile(shape_fn) == True:
                os.remove(shape_fn)

            # need to tell ogr which driver to use
            driver = ogr.GetDriverByName('ESRI Shapefile')

            if os.path.isfile(shape_fn) == True:
                driver.DeleteDataSource(shape_fn)

            # create shape file
            data_source = driver.CreateDataSource(shape_fn)

            # if you read from a raster get the georeference point otherwise create one
            # spatial_ref = osr.SpatialReference()
            # this puts it in the wsg84 reference frame.
            # spatial_ref.ImportFromEPSG(self._proj_dict[self.projection])

            # create a layer to put the ellipses onto
            layer = data_source.CreateLayer('TIPPER', self.utm_cs,
                                            ogr.wkbPolygon)

            # make field names
            field_name = ogr.FieldDefn("Name", ogr.OFTString)
            layer.CreateField(field_name)

            field_mag_imag = ogr.FieldDefn('mag_imag', ogr.OFTReal)
            layer.CreateField(field_mag_imag)

            field_ang_imag = ogr.FieldDefn('ang_imag', ogr.OFTReal)
            layer.CreateField(field_ang_imag)

            for tp_arr in self.tip_dict[plot_per]:
                cos_t = np.cos(-np.deg2rad(tp_arr['ang_imag']))
                sin_t = np.sin(-np.deg2rad(tp_arr['ang_imag']))
                # calculate the points to make the line
                txr = 0
                tyr = tp_arr['mag_imag'] * self.arrow_size

                # make an arrow by drawing an outline.  have the arrow point
                # north to start and then rotate later with the rotation
                # matrix to properly orient it.
                x0 = 0
                y0 = 0

                x1 = x0 + self.arrow_lw
                y1 = y0

                x2 = x0 + self.arrow_lw
                y2 = y0 + tyr - self.arrow_head_height

                x3 = x0 + self.arrow_lw + self.arrow_head_width
                y3 = y2

                x4 = x0 + txr
                y4 = y0 + tyr

                x7 = x0 - self.arrow_lw
                y7 = y0

                x6 = x0 - self.arrow_lw
                y6 = y0 + tyr - self.arrow_head_height

                x5 = x0 - self.arrow_lw - self.arrow_head_width
                y5 = y6

                x = np.array([x0, x1, x2, x3, x4, x5, x6, x7])
                y = np.array([y0, y1, y2, y3, y4, y5, y6, y7])

                rot_matrix = np.array([[cos_t, -sin_t], [sin_t, cos_t]])

                # rotate the arrow to be properly oriented
                xy = np.array([x, y])
                rot_xy = np.dot(rot_matrix, xy)

                # shift the arrow to be centered on the station
                x = tp_arr['east'] + rot_xy[0]
                y = tp_arr['north'] + rot_xy[1]

                # 1) make a geometry shape line
                arrow = ogr.Geometry(ogr.wkbLinearRing)
                for ii, jj in zip(x, y):
                    arrow.AddPoint(np.round(ii, 6), np.round(jj, 6))

                arrow.CloseRings()

                poly = ogr.Geometry(ogr.wkbPolygon)
                poly.AddGeometry(arrow)
                # 4) this part is confusing but we need to create a
                # feature that has the
                # same definition as the layer that we created.
                #    get the layer definition
                feature_def = layer.GetLayerDefn()

                # create a new feature
                new_feature = ogr.Feature(feature_def)
                # set the geometry of that feature to be the ellipse
                new_feature.SetGeometry(poly)
                # create the feature in the layer.
                layer.CreateFeature(new_feature)

                #
                # 5) create a field to color by
                new_feature.SetField("Name", tp_arr['station'])
                new_feature.SetField("mag_imag", tp_arr['mag_imag'])
                new_feature.SetField("ang_imag", tp_arr['ang_imag'])

                # add the new feature to the layer.
                layer.SetFeature(new_feature)

                # apparently need to destroy the feature
                new_feature.Destroy()

            # Need to be sure that all the new info is saved to
            data_source.SyncToDisk()

            # write a projection file
            # spatial_ref.MorphFromESRI()
            self.utm_cs.MorphFromESRI()
            prj_file = open('{0}prj'.format(shape_fn[:-3]), 'w')
            prj_file.write(self.utm_cs.ExportToWkt())
            prj_file.close()

            data_source.Destroy()

            print('Wrote shape file to {0}'.format(shape_fn))

    def write_tip_shape_files_modem(self, modem_data_fn, rotation_angle=0.0):
        """
        write tip files from a modem data file.

        """

        modem_obj = mtpy.modeling.modem.Data()
        modem_obj.read_data_file(modem_data_fn)

        self.plot_period = modem_obj.period_list.copy()
        self.mt_obj_list = [modem_obj.mt_dict[key]
                            for key in list(modem_obj.mt_dict.keys())]

        self._set_rotation_angle(rotation_angle)

        self.write_imag_shape_files()
        self.write_real_shape_files()

    def write_tip_shape_files_modem_residual(self, modem_data_fn,
                                             modem_resp_fn,
                                             rotation_angle):
        """
        write residual tipper files for modem

        """
        modem_data_obj = mtpy.modeling.modem.Data()
        modem_data_obj.read_data_file(modem_data_fn)

        modem_resp_obj = mtpy.modeling.modem.Data()
        modem_resp_obj.read_data_file(modem_resp_fn)

        self.plot_period = modem_data_obj.period_list.copy()
        mt_keys = sorted(modem_data_obj.mt_dict.keys())
        self.mt_obj_list = [modem_data_obj.mt_dict[key]
                            for key in mt_keys]

        self._set_rotation_angle(rotation_angle)

        for mt_obj, key in zip(self.mt_obj_list, mt_keys):
            try:
                resp_tipper = modem_resp_obj.mt_dict[key].Tipper.tipper
                mt_obj.Tipper.tipper[:, :, :] -= resp_tipper[:, :, :]
            except KeyError:
                continue

        self.write_imag_shape_files()
        self.write_real_shape_files()


# ==============================================================================
# reproject a layer DOESNT WORK YET
# ==============================================================================
def reproject_layer(in_shape_file, out_shape_file=None, out_proj='WGS84'):
    """
    reproject coordinates into a different coordinate system

    """

    driver = ogr.GetDriverByName('ESRI Shapefile')

    # input SpatialReference
    inSpatialRef = osr.SpatialReference()
    inSpatialRef.ImportFromEPSG(2927)

    # output SpatialReference
    outSpatialRef = osr.SpatialReference()
    outSpatialRef.ImportFromEPSG(4326)

    # create the CoordinateTransformation
    coordTrans = osr.CoordinateTransformation(inSpatialRef, outSpatialRef)

    # get the input layer
    inDataSet = driver.Open(r'c:\data\spatial\basemap.shp')
    inLayer = inDataSet.GetLayer()

    # create the output layer
    outputShapefile = r'c:\data\spatial\basemap_4326.shp'
    if os.path.exists(outputShapefile):
        driver.DeleteDataSource(outputShapefile)
    outDataSet = driver.CreateDataSource(outputShapefile)
    outLayer = outDataSet.CreateLayer(
        "basemap_4326", geom_type=ogr.wkbMultiPolygon)

    # add fields
    inLayerDefn = inLayer.GetLayerDefn()
    for i in range(0, inLayerDefn.GetFieldCount()):
        fieldDefn = inLayerDefn.GetFieldDefn(i)
        outLayer.CreateField(fieldDefn)

    # get the output layer's feature definition
    outLayerDefn = outLayer.GetLayerDefn()

    # loop through the input features
    inFeature = inLayer.GetNextFeature()
    while inFeature:
        # get the input geometry
        geom = inFeature.GetGeometryRef()
        # reproject the geometry
        geom.Transform(coordTrans)
        # create a new feature
        outFeature = ogr.Feature(outLayerDefn)
        # set the geometry and attribute
        outFeature.SetGeometry(geom)
        for i in range(0, outLayerDefn.GetFieldCount()):
            outFeature.SetField(
                outLayerDefn.GetFieldDefn(i).GetNameRef(),
                inFeature.GetField(i))
        # add the feature to the shapefile
        outLayer.CreateFeature(outFeature)
        # destroy the features and get the next input feature
        outFeature.Destroy()
        inFeature.Destroy()
        inFeature = inLayer.GetNextFeature()

    # close the shapefiles
    inDataSet.Destroy()
    outDataSet.Destroy()


# ==============================================================================
# create a raster from an array
# ==============================================================================

def array2raster(newRasterfn, rasterOrigin, pixelWidth, pixelHeight, array):
    cols = array.shape[1]
    rows = array.shape[0]
    originX = rasterOrigin[0]
    originY = rasterOrigin[1]

    driver = gdal.GetDriverByName('GTiff')
    outRaster = driver.Create(newRasterfn, cols, rows, 1, gdal.GDT_Byte)
    outRaster.SetGeoTransform(
        (originX, pixelWidth, 0, originY, 0, pixelHeight))
    outband = outRaster.GetRasterBand(1)
    outband.WriteArray(array)
    outRasterSRS = osr.SpatialReference()
    outRasterSRS.ImportFromEPSG(4326)
    outRaster.SetProjection(outRasterSRS.ExportToWkt())
    outband.FlushCache()


#==============================================================================
# test
#==============================================================================
##edipath = r"c:\Users\jrpeacock\Documents\Mendenhall\MonoBasin\EDI_Files\GeographicNorth"
##edilst = [os.path.join(edipath, edi) for edi in os.listdir(edipath)
##          if edi.find('.edi') > 0]
##edilst.remove(os.path.join(edipath, 'mb035.edi'))
##
##pts = PTShapeFile(edilst, save_path=r"c:\Users\jrpeacock")
##pts.projection = 'NAD27'
##pts.ellipse_size = 1200
##pts.write_shape_files()
#
#
##tps = TipperShapeFile(edilst, save_path=r"c:\Users\jrpeacock")
##tps.projection = 'NAD27'
##tps.arrow_lw = 30
##tps.arrow_head_height = 100
##tps.arrow_head_width = 70
##tps.write_real_shape_files()
##tps.write_imag_shape_files()
#
#mfn = r"c:\Users\jrpeacock\Google Drive\Mono_Basin\Models\Modular_NLCG_110.dat"
#sv_path = r"c:\Users\jrpeacock\Google Drive\Mono_Basin\Models\GIS_Tip_Response"
##sv_path = r"c:\Users\jrpeacock\Google Drive\Mono_Basin\Models\GIS_PT_Response"
##pts = PTShapeFile(save_path=sv_path)
##pts.projection = 'NAD27'
##pts.ellipse_size = 1200
##pts.write_pt_shape_files_modem(mfn)
#
#tps = TipperShapeFile(save_path=sv_path)
#tps.projection = 'NAD27'
#tps.arrow_lw = 30
#tps.arrow_head_height = 100
#tps.arrow_head_width = 70
#tps.write_tip_shape_files_modem(mfn)


def modem_to_shapefiles(mfndat, save_dir):
    """
    create shape file representaiotn for ModEM model
    :param mfndat: \path2\Modular_NLCG_110.dat
    :param save_dir: \path2\outshp
    :return:
    """

    pts = PTShapeFile(save_path=save_dir)
    pts.projection = None  #'WGS84'  # ''NAD27'
    pts.ellipse_size = 0.1
    pts.write_data_pt_shape_files_modem(mfndat)

    # tipper shape run OK, need check results
    tipshp = TipperShapeFile(save_path=save_dir)
    tipshp.projection = 'WGS84'   #'NAD27'
    tipshp.arrow_lw = 30
    tipshp.arrow_head_height = 100
    tipshp.arrow_head_width = 70
    tipshp.write_tip_shape_files_modem(mfndat)

    return


def create_phase_tensor_shpfiles(
        edi_dir, save_dir, proj='WGS84', ellipse_size=1000, every_site=1, period_list=None):
    """
    generate shape file for a folder of edi files, and save the shape files a dir.
    :param edi_dir:
    :param save_dir:
    :param proj: defult is WGS84-UTM, with ellipse_size=1000 meters
    :param ellipse_size: the size of ellipse: 100-5000, try them out to suit your needs
    :param every_site: by default every MT station will be output, but user can sample down with 2, 3,..
    :return:
    """

    edipath = edi_dir

    edilst = [os.path.join(edipath, edi) for edi in os.listdir(edipath)
              if edi.find('.edi') > 0]
    # edilst.remove(os.path.join(edipath, 'mb035.edi'))

    #subset of the edilst:
    edilst2 = edilst[::every_site]
    pts = PTShapeFile(edilst2, save_path=save_dir, proj=proj)

    pts.ellipse_size = ellipse_size

    pts.write_shape_files( periods=period_list )


def create_tipper_shpfiles(edipath, save_dir):
    """
    Create Tipper (induction arrows real and imaginary) shape files
    :param edipath:
    :param save_dir:
    :return:
    """

    edilist = [os.path.join(edipath, edi) for edi in os.listdir(edipath)
               if edi.find('.edi') > 0]

    tipshp = TipperShapeFile(edilist, save_path=save_dir)

    #tipshp.projection = 'NAD27'
    tipshp.arrow_lw = 30
    tipshp.arrow_head_height = 100
    tipshp.arrow_head_width = 70

    tipshp.write_real_shape_files()
    tipshp.write_imag_shape_files()

    return

def create_modem_data_shapefiles():
    # modem: provide dat file and save_path below:
    # mfn = r"E:/Githubz/mtpy/examples/data/ModEM_files/VicSynthetic07/Modular_MPI_NLCG_016.dat"
    mfn = r"E:\Data\Modeling\Isa\100hs_flat_BB\Isa_run3_NLCG_048.dat"
    save_path = r"C:/tmp"
    modem_to_shapefiles(mfn, save_path)


# ===================================================
#  main test
# edi files-dir as input
#  python mtpy/utils/shapefiles.py  data/edifiles c:/temp/
#  python mtpy/utils/shapefiles.py examples/data/edi_files c:/temp/
# modem dat file as input
#  python mtpy/utils/shapefiles.py /e/Data/Modeling/Isa/100hs_flat_BB/Isa_run3_NLCG_048.dat c:/temp2/
#  python mtpy/utils/shapefiles.py examples/data/ModEM_files/VicSynthetic07/Modular_MPI_NLCG_016.dat c:/temp2/
# ----------------------------------------------------

# example codes of how to use this module
if __name__ == "__main__d":
    import sys
    if len(sys.argv) < 3:
        print((
            "USAGE: %s input_edifile_dir output_shape_file_dir" %
            sys.argv[0]))
        sys.exit(1)

    src_file_dir = sys.argv[1] # A modem data file  OR edi-folder

    dest_dir = sys.argv[2]
    if not os.path.isdir(dest_dir):
        os.mkdir(dest_dir)

    if src_file_dir[-4:].lower() == '.dat':      # modem dat file
        modem_to_shapefiles(src_file_dir, dest_dir)
    elif os.path.isdir(src_file_dir):            # input from edi folder
        create_phase_tensor_shpfiles(
            src_file_dir,
            dest_dir,
            proj='WGS84', ellipse_size=8000, # UTM and size in meters. 1deg=100KM
            #proj=None,  ellipse_size=0.01, # Lat-Long geographic coord and size in degree
            every_site=2,
            #period_list=[ 218.43599825251204, 218.43599825 ] # KeyError 218.43599825
            #period_list=[0.0128, 0.016]  # must get very accurate periods
        )

        # # create_phase_tensor_shpfiles(sys.argv[1], sys.argv[2], ,
        # # ellipse_size=3000, every_site=2) # projected into UTM coordinate

        #create_tipper_shpfiles(sys.argv[1],sys.argv[2])
    else:
        print("Nothing to do !!!")

# ===================================================
# Command Wrapper for shape files generation
# ===================================================
@click.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.option('-i','--input',type=str,
              default='examples/data/edi_files', \
              help='input edsi files dir  or Modem dat file examples/data/MoDEM_files/Modular_MPI_NLCG_028.dat')
@click.option('-o','--output',type=str,default="temp",help='Output directory')
def generate_shape_files(input,output):
    print("=======================================================================")
    print("Generating Shapes File requires following inputs edsi files directory")
    print("                                              or MoDEM input file")
    print("Default output is in temp directory                                     ")
    print("=======================================================================")
    if not os.path.isdir(output):
        os.mkdir(output)
    print(("input = {}".format(input)))
    print(("input = {}".format(input[-4:].lower())))
    if input[-4:].lower() == '.dat':
        modem_to_shapefiles(input, output)
    elif os.path.isdir(input):
        create_phase_tensor_shpfiles(
            input,
            output,
            proj='WGS84', ellipse_size=8000, # UTM and size in meters. 1deg=100KM
                                             # proj=None,  ellipse_size=0.01,
                                             # Lat-Long geographic coord and size in degree
            every_site=2,                    #period_list=[ 218.43599825251204, 218.43599825 ]
                                             #  KeyError 218.43599825
                                             #period_list=[0.0128, 0.016]  # must get very accurate periods
        )
    else:
        print ("Nothing to do !")


if __name__ == "__main__":
    generate_shape_files()
