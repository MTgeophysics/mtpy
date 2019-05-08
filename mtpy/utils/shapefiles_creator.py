#! /usr/bin/env python
"""
Description:
    Create shape files for Phase Tensor Ellipses, Tipper Real/Imag.
    export the phase tensor map and tippers into jpeg/png images

CreationDate:   2017-03-06
Developer:      fei.zhang@ga.gov.au

Revision History:
    LastUpdate:     10/11/2017   FZ fix bugs after the big merge
    LastUpdate:     20/11/2017   change from freq to period filenames, allow to specify a period
    LastUpdate:     30/10/2018   combine ellipses and tippers together, refactorings

"""



import glob
import logging
import os
import sys

import click
import geopandas as gpd
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from mpl_toolkits.axes_grid1 import make_axes_locatable
from shapely.geometry import Point, Polygon, LineString, LinearRing

from mtpy.core.edi_collection import EdiCollection
from mtpy.utils.decorator import deprecated
from mtpy.utils.mtpylog import MtPyLog
from mtpy.utils.edi_folders import recursive_glob

mpl.rcParams['lines.linewidth'] = 2
# mpl.rcParams['lines.color'] = 'r'
mpl.rcParams['figure.figsize'] = [10, 6]

_logger = MtPyLog.get_mtpy_logger(__name__)  # logger inside this file/module
_logger.setLevel(logging.DEBUG)  # set your logger level


class ShapeFilesCreator(EdiCollection):
    """ Extend the EdiCollection parent class,
    create phase tensor and tipper shapefiles for a list of edifiles

    :param edifile_list: [path2edi,...]
    :param outdir: path2output dir, where the shp file will be written.
    :param orig_crs = {'init': 'epsg:4283'}  # GDA94
    """

    def __init__(self, edifile_list, outdir, orig_crs={'init': 'epsg:4283'}):
        """
        loop through a list of edi files, create required shapefiles
        :param edifile_list: [path2edi,...]
        :param outdir: path2output dir, where the shp file will be written.
        :param orig_crs = {'init': 'epsg:4283'}  # GDA94
        # {'init': 'epsg:4326'}  # WGS84
        """

        self.orig_crs = orig_crs

        # ensure that outdir is specified, and be created if not there.
        if outdir is None:
            raise Exception("Error: OutputDir is not specified!!!")
        elif not os.path.exists(outdir):
            os.mkdir(outdir)

        self.outdir = outdir

        # call the super constructor
        super(ShapeFilesCreator, self).__init__(edilist=edifile_list, outdir=outdir)
        # python-3 syntax: super().__init__(edilist=edifile_list, outdir=outdir)

        self.stations_distances = self.get_stations_distances_stats()

        # These attributes below are defined in the parent class.
        # self.all_periods = self._get_all_periods()
        # self.ptol = 0.05  # this param controls how near-equal freqs/periods are grouped together:
        # 10% may result more double countings of freq/periods than 5%.
        # eg: E:\Data\MT_Datasets\WenPingJiang_EDI 18528 rows vs 14654 rows

        return

    def create_phase_tensor_shp(self, period, ellipsize=None, target_epsg_code=4283, export_fig=False):
        """
        create phase tensor ellipses shape file correspond to a MT period
        :return: (geopdf_obj, path_to_shapefile)
        """

        if ellipsize is None:  # automatically decide suitable ellipse size.
            ellipsize = self.stations_distances.get("Q1PERCENT") / 2  # a half or a third of the min_distance?
            self._logger.info("Automatically Selected Max-Ellispse Size = %s", ellipsize)

        pt = self.get_phase_tensor_tippers(period)

        self._logger.debug("phase tensor values =: %s", pt)

        if len(pt) < 1:
            self._logger.warn("No phase tensor for the period %s for any MT station", period)
            return None

        pdf = pd.DataFrame(pt)

        self._logger.debug(pdf['period'])

        mt_locations = [Point(xy) for xy in zip(pdf['lon'], pdf['lat'])]
        # OR pdf['geometry'] = pdf.apply(lambda z: Point(z.lon, z.lat), axis=1)
        # if you want to df = df.drop(['Lon', 'Lat'], axis=1)
        # orig_crs = {'init': 'epsg:4326'}  # initial crs WGS84
        # orig_crs = {'init': 'epsg:4283'}  # initial crs GDA94

        geopdf = gpd.GeoDataFrame(pdf, crs=self.orig_crs, geometry=mt_locations)

        # make  pt_ellispes using polygons
        phi_max_v = geopdf['phi_max'].max()  # the max of this group of ellipse

        print(phi_max_v)

        # points to trace out the polygon-ellipse

        theta = np.arange(0, 2 * np.pi, np.pi / 30.)

        azimuth = -np.deg2rad(geopdf['azimuth'])
        width = ellipsize * (geopdf['phi_max'] / phi_max_v)
        height = ellipsize * (geopdf['phi_min'] / phi_max_v)
        x0 = geopdf['lon']
        y0 = geopdf['lat']

        # apply formula to generate ellipses

        ellipse_list = []
        for i in range(0, len(azimuth)):
            x = x0[i] + height[i] * np.cos(theta) * np.cos(azimuth[i]) - \
                width[i] * np.sin(theta) * np.sin(azimuth[i])
            y = y0[i] + height[i] * np.cos(theta) * np.sin(azimuth[i]) + \
                width[i] * np.sin(theta) * np.cos(azimuth[i])

            polyg = Polygon(LinearRing([xy for xy in zip(x, y)]))

            # print polyg  # an ellispe

            ellipse_list.append(polyg)

        geopdf = gpd.GeoDataFrame(geopdf, crs=self.orig_crs, geometry=ellipse_list)

        if target_epsg_code is None:
            self._logger.info("The orginal Geopandas Dataframe CRS: %s", geopdf.crs)
            # {'init': 'epsg:4283', 'no_defs': True}
            # raise Exception("Must provide a target_epsg_code")
            target_epsg_code = geopdf.crs['init'][5:]
        else:
            geopdf.to_crs(epsg=target_epsg_code, inplace=True)
            # world = world.to_crs({'init': 'epsg:3395'})
            # world.to_crs(epsg=3395) would also work

        # to shape file
        shp_fname = 'Phase_Tensor_EPSG_%s_Period_%ss.shp' % (target_epsg_code, period)
        path2shp = os.path.join(self.outdir, shp_fname)
        self._logger.debug("To write to ESRI shp file %s", path2shp)
        geopdf.to_file(path2shp, driver='ESRI Shapefile')

        self._logger.info("Geopandas Dataframe CRS: %s", geopdf.crs)

        if export_fig is True:
            bbox_dict = self.get_bounding_box(epsgcode=target_epsg_code)
            # this bbox ensures that the whole MT-stations area is covered independent of periods
            print(bbox_dict)
            path2jpg = path2shp.replace(".shp", ".jpg")
            export_geopdf_to_image(geopdf, bbox_dict, path2jpg, colorby='phi_max',
                                   colormap='nipy_spectral_r')  # showfig=True)

        return (geopdf, path2shp)

    def create_tipper_real_shp(self, period, line_length=None, target_epsg_code=4283, export_fig=False):
        """
        create real tipper lines shapefile from a csv file
        The shapefile consists of lines without arrow.
        User can use GIS software such as ArcGIS to display and add an arrow at each line's end
        line_length is how long will be the line, auto-calculatable
        """

        if line_length is None:  # auto-calculate the tipper arrow length
            line_length = self.stations_distances.get("Q1PERCENT")
            self._logger.info("Automatically Selected Max Tipper Length  = %s", line_length)

        pt = self.get_phase_tensor_tippers(period)
        self._logger.debug("phase tensor values =: %s", pt)

        if len(pt) < 1:
            self._logger.warn("No phase tensor for the period %s for any MT station", period)
            return None

        pdf = pd.DataFrame(pt)

        tip_mag_re_maxval = pdf['tip_mag_re'].max()

        if (tip_mag_re_maxval > 0.00000001):
            line_length_normalized = line_length / tip_mag_re_maxval
        else:
            line_length_normalized = line_length

        self._logger.debug(pdf['period'])

        pdf['tip_re'] = pdf.apply(lambda x:
                                  LineString([(float(x.lon), float(x.lat)),
                                              (float(x.lon) + line_length_normalized * x.tip_mag_re * np.cos(
                                                  -np.deg2rad(x.tip_ang_re)),
                                               float(x.lat) + line_length_normalized * x.tip_mag_re * np.sin(
                                                   -np.deg2rad(x.tip_ang_re)))]), axis=1)

        geopdf = gpd.GeoDataFrame(pdf, crs=self.orig_crs, geometry='tip_re')

        if target_epsg_code is None:
            self._logger.info("Geopandas Datframe CRS: %s", geopdf.crs)
            # {'init': 'epsg:4283', 'no_defs': True}
            # raise Exception("Must provide a target_epsg_code")
            target_epsg_code = geopdf.crs['init'][5:]
        else:
            geopdf.to_crs(epsg=target_epsg_code, inplace=True)
            # world = world.to_crs({'init': 'epsg:3395'})
            # world.to_crs(epsg=3395) would also work

        # to shape file
        shp_fname = 'Tipper_Real_EPSG_%s_Period_%ss.shp' % (target_epsg_code, period)
        path2shp = os.path.join(self.outdir, shp_fname)
        self._logger.debug("To write to ESRI shp file %s", path2shp)
        geopdf.to_file(path2shp, driver='ESRI Shapefile')

        self._logger.info("Geopandas Dataframe CRS: %s", geopdf.crs)

        if export_fig is True:
            bbox_dict = self.get_bounding_box(epsgcode=target_epsg_code)
            # this bbox_dict ensures that we can set a consistent display area cover all ground stations,
            # not just this period-dependent geopdf
            self._logger.debug("All MT stations area bounding box %s", bbox_dict)
            path2jpg = path2shp.replace(".shp", ".jpg")
            export_geopdf_to_image(geopdf, bbox_dict, path2jpg, colorby='phi_max',
                                   colormap='nipy_spectral_r')  # showfig=True)

        return (geopdf, path2shp)

    def create_tipper_imag_shp(self, period, line_length=None, target_epsg_code=4283, export_fig=False):
        """
        create imagery tipper lines shapefile from a csv file
        The shapefile consists of lines without arrow.
        User can use GIS software such as ArcGIS to display and add an arrow at each line's end
        line_length is how long will be the line, auto-calculatable
        :return:(geopdf_obj, path_to_shapefile)
        """

        if line_length is None:  # auto-calculate the tipper arrow length
            line_length = self.stations_distances.get("Q1PERCENT")
            self._logger.info("Automatically Selected Max-Tipper Length =: %s", line_length)

        pt = self.get_phase_tensor_tippers(period)
        self._logger.debug("phase tensor values =: %s", pt)

        if len(pt) < 1:
            self._logger.warn("No phase tensor for the period %s for any MT station", period)
            return None

        pdf = pd.DataFrame(pt)

        tip_mag_im_maxval = pdf['tip_mag_im'].max()

        if (tip_mag_im_maxval > 0.00000001):
            line_length_normalized = line_length / tip_mag_im_maxval
        else:
            line_length_normalized = line_length

        self._logger.debug(pdf['period'])

        pdf['tip_im'] = pdf.apply(lambda x: LineString([(float(x.lon), float(x.lat)),
                                                        (float(x.lon) + line_length_normalized * x.tip_mag_im * np.cos(
                                                            -np.deg2rad(x.tip_ang_im)),
                                                         float(x.lat) + line_length_normalized * x.tip_mag_im * np.sin(
                                                             -np.deg2rad(x.tip_ang_im)))]),
                                  axis=1)

        geopdf = gpd.GeoDataFrame(pdf, crs=self.orig_crs, geometry='tip_im')

        if target_epsg_code is None:
            self._logger.info("Keep the Default/Original Geopandas Dataframe CRS: %s", geopdf.crs)
            # {'init': 'epsg:4283', 'no_defs': True}
            # raise Exception("Must provide a target_epsg_code")
            target_epsg_code = geopdf.crs['init'][5:]
        else:
            geopdf.to_crs(epsg=target_epsg_code, inplace=True)
            # world = world.to_crs({'init': 'epsg:3395'})
            # world.to_crs(epsg=3395) would also work

        # to shape file
        shp_fname = 'Tipper_Imag_EPSG_%s_Period_%ss.shp' % (target_epsg_code, period)
        path2shp = os.path.join(self.outdir, shp_fname)
        self._logger.debug("To write to ESRI shp file %s", path2shp)
        geopdf.to_file(path2shp, driver='ESRI Shapefile')

        self._logger.info("Geopandas Dataframe CRS: %s", geopdf.crs)

        if export_fig is True:
            bbox_dict = self.get_bounding_box(epsgcode=target_epsg_code)
            # this bbox_dict ensures that we can set a consistent display area cover all ground stations,
            # not just this period-dependent geopdf

            self._logger.debug("All MT stations area bounding box %s", bbox_dict)

            path2jpg = path2shp.replace(".shp", ".jpg")
            export_geopdf_to_image(geopdf, bbox_dict, path2jpg, colorby='phi_max',
                                   colormap='nipy_spectral_r')  # showfig=True)

        return (geopdf, path2shp)


def plot_phase_tensor_ellipses_and_tippers(edi_dir, outfile=None, iperiod=0):
    """
    plot phase tensor ellipses and tipers into one figure.
    :param edi_dir: edi directory
    :param outfile: save figure to output file
    :param iperiod: the index of periods
    :return: saved figure file
    """

    edifiles = recursive_glob(edi_dir)

    print("Number of EDI files found = %s" % len(edifiles))

    myobj = ShapeFilesCreator(edifiles, "c:/temp")

    allper = myobj.all_unique_periods

    gpd_phtensor = myobj.create_phase_tensor_shp(allper[iperiod], export_fig=False)[0]

    gpd_retip = myobj.create_tipper_real_shp(allper[iperiod], export_fig=False)[0]

    gpd_imtip = myobj.create_tipper_imag_shp(allper[iperiod], export_fig=False)[0]

    world = gpd.read_file(gpd.datasets.get_path('naturalearth_lowres'))

    # composing two layers in a map
    f, ax = plt.subplots(1, figsize=(20, 12))

    # ax.set_xlim([140.5,141])
    # ax.set_ylim([-21,-20])

    # Add layer of polygons on the axis

    # world.plot(ax=ax, alpha=0.5)  # background map
    gpd_phtensor.plot(ax=ax, linewidth=2, facecolor='grey', edgecolor='black')
    gpd_retip.plot(ax=ax, color='red', linewidth=4)
    gpd_imtip.plot(ax=ax, color='blue', linewidth=4)

    if outfile is not None:
        plt.savefig(outfile)  # ( 'C:/temp/phase_tensor_tippers.png')

    # Display
    # plt.show()

    return outfile

####################################################################
# Using geopandas to convert CSV files into shape files
# Refs:
#   http://toblerity.org/shapely/manual.html#polygons
#   https://geohackweek.github.io/vector/04-geopandas-intro/
# ===================================================================
def create_ellipse_shp_from_csv(csvfile, esize=0.03, target_epsg_code=4283):
    """
    create phase tensor ellipse geometry from a csv file. This function needs csv file as its input.
    :param csvfile: a csvfile with full path
    :param esize: ellipse size, defaut 0.03 is about 3KM in the max ellipse rad
    :return: a geopandas dataframe
    """
    # crs = {'init': 'epsg:4326'}  # if assume initial crs WGS84
    crs = {'init': 'epsg:4283'}  # if assume initial crs GDA94

    pdf = pd.read_csv(csvfile)
    mt_locations = [Point(xy) for xy in zip(pdf['lon'], pdf['lat'])]
    # OR pdf['geometry'] = pdf.apply(lambda z: Point(z.lon, z.lat), axis=1)
    # if you want to df = df.drop(['Lon', 'Lat'], axis=1)

    pdf = gpd.GeoDataFrame(pdf, crs=crs, geometry=mt_locations)

    # make  pt_ellispes using polygons
    phi_max_v = pdf['phi_max'].max()  # the max of this group of ellipse

    print(phi_max_v)

    # points to trace out the polygon-ellipse
    theta = np.arange(0, 2 * np.pi, np.pi / 30.)

    azimuth = -np.deg2rad(pdf['azimuth'])
    width = esize * (pdf['phi_max'] / phi_max_v)
    height = esize * (pdf['phi_min'] / phi_max_v)
    x0 = pdf['lon']
    y0 = pdf['lat']

    # apply formula to generate ellipses

    ellipse_list = []
    for i in range(0, len(azimuth)):
        x = x0[i] + height[i] * np.cos(theta) * np.cos(azimuth[i]) - \
            width[i] * np.sin(theta) * np.sin(azimuth[i])
        y = y0[i] + height[i] * np.cos(theta) * np.sin(azimuth[i]) + \
            width[i] * np.sin(theta) * np.cos(azimuth[i])

        polyg = Polygon(LinearRing([xy for xy in zip(x, y)]))

        # print polyg  # an ellispe

        ellipse_list.append(polyg)

    pdf = gpd.GeoDataFrame(pdf, crs=crs, geometry=ellipse_list)

    if target_epsg_code is None:
        raise Exception("Must provide a target_epsg_code")
    else:
        pdf.to_crs(epsg=target_epsg_code, inplace=True)
        # world = world.to_crs({'init': 'epsg:3395'})
        # world.to_crs(epsg=3395) would also work

    # to shape file
    shp_fname = csvfile.replace('.csv', '_ellip_epsg%s.shp' % target_epsg_code)
    pdf.to_file(shp_fname, driver='ESRI Shapefile')

    return pdf


def create_tipper_real_shp_from_csv(csvfile, line_length=0.03, target_epsg_code=4283):
    """ create tipper lines shape from a csv file. This function needs csv file as its input.
    The shape is a line without arrow.
    Must use a GIS software such as ArcGIS to display and add an arrow at each line's end
    line_length=4  how long will be the line (arrow)
    return: a geopandas dataframe object for further processing.
    """

    # crs = {'init': 'epsg:4326'}  # if assume initial crs WGS84
    crs = {'init': 'epsg:4283'}  # if assume initial crs GDA94

    pdf = pd.read_csv(csvfile)
    # mt_locations = [Point(xy) for xy in zip(pdf.lon, pdf.lat)]
    # OR pdf['geometry'] = pdf.apply(lambda z: Point(z.lon, z.lat), axis=1)
    # if you want to df = df.drop(['Lon', 'Lat'], axis=1)

    # geo_df = gpd.GeoDataFrame(pdf, crs=crs, geometry=mt_locations)

    pdf['tip_re'] = pdf.apply(lambda x:
                              LineString([(float(x.lon), float(x.lat)),
                                          (float(x.lon) + line_length * x.tip_mag_re * np.cos(
                                              -np.deg2rad(x.tip_ang_re)),
                                           float(x.lat) + line_length * x.tip_mag_re * np.sin(
                                               -np.deg2rad(x.tip_ang_re)))]), axis=1)

    pdf = gpd.GeoDataFrame(pdf, crs=crs, geometry='tip_re')

    if target_epsg_code is None:
        raise Exception("Must provide a target_epsg_code")
    else:
        pdf.to_crs(epsg=target_epsg_code, inplace=True)
        # world = world.to_crs({'init': 'epsg:3395'})
        # world.to_crs(epsg=3395) would also work

    # to shape file
    shp_fname = csvfile.replace('.csv', '_real_epsg%s.shp' % target_epsg_code)
    pdf.to_file(shp_fname, driver='ESRI Shapefile')

    return pdf


def create_tipper_imag_shp_from_csv(csvfile, line_length=0.03, target_epsg_code=4283):
    """ create imagery tipper lines shape from a csv file. this function needs csv file as input.
    The shape is a line without arrow.
    Must use a GIS software such as ArcGIS to display and add an arrow at each line's end
    line_length=4  how long will be the line (arrow)
    return: a geopandas dataframe object for further processing.
    """

    # crs = {'init': 'epsg:4326'}  # if assume initial crs WGS84
    crs = {'init': 'epsg:4283'}  # if assume initial crs GDA94

    pdf = pd.read_csv(csvfile)
    # mt_locations = [Point(xy) for xy in zip(pdf.lon, pdf.lat)]
    # OR pdf['geometry'] = pdf.apply(lambda z: Point(z.lon, z.lat), axis=1)
    # if you want to df = df.drop(['Lon', 'Lat'], axis=1)

    # geo_df = gpd.GeoDataFrame(pdf, crs=crs, geometry=mt_locations)

    pdf['tip_im'] = pdf.apply(lambda x:
                              LineString([(float(x.lon), float(x.lat)),
                                          (float(x.lon) + line_length * x.tip_mag_im * np.cos(
                                              -np.deg2rad(x.tip_ang_im)),
                                           float(x.lat) + line_length * x.tip_mag_im * np.sin(
                                               -np.deg2rad(x.tip_ang_im)))]), axis=1)

    pdf = gpd.GeoDataFrame(pdf, crs=crs, geometry='tip_im')

    if target_epsg_code is None:
        raise Exception("Must provide a target_epsg_code")  # EDI original lat/lon epsg 4326 or GDA94
    else:
        pdf.to_crs(epsg=target_epsg_code, inplace=True)
        # world = world.to_crs({'init': 'epsg:3395'})
        # world.to_crs(epsg=3395) would also work

    # to shape file
    shp_fname = csvfile.replace('.csv', '_imag_epsg%s.shp' % target_epsg_code)
    pdf.to_file(shp_fname, driver='ESRI Shapefile')

    return pdf


def export_geopdf_to_image(geopdf, bbox, jpg_file_name, target_epsg_code=None, colorby=None, colormap=None,
                           showfig=False):
    """
    Export a geopandas dataframe to a jpe_file, with optionally a new epsg projection.
    :param geopdf: a geopandas dataframe
    :param bbox: This param ensures that we can set a consistent display area defined by a dict with 4 keys
                    [MinLat, MinLon, MaxLat, MaxLon], cover all ground stations, not just this period-dependent geopdf
    :param output jpg_file_name: path2jpeg
    :param target_epsg_code: 4326 etc
    :param showfig: If True, then display fig on screen.
    :return:
    """

    if target_epsg_code is None:
        p = geopdf
        # target_epsg_code = '4283'  # EDI orginal lat/lon epsg 4326=WGS84 or 4283=GDA94
        target_epsg_code = geopdf.crs['init'][5:]
    else:
        p = geopdf.to_crs(epsg=target_epsg_code)
        # world = world.to_crs({'init': 'epsg:3395'})
        # world.to_crs(epsg=3395) would also work

    # bounds = p.total_bounds  # lat-lon bounds for this csv dataframe

    # plot and save

    fig_title = os.path.basename(jpg_file_name)
    _logger.info('saving figure to file %s', jpg_file_name)

    if colorby is None:
        colorby = 'phi_min'
    else:
        colorby = colorby

    if colormap is None:
        my_colormap = mpl.cm.gist_ncar  # a default choice:  jet_r #'jet'
    else:
        my_colormap = colormap

    if int(target_epsg_code) == 4326 or int(target_epsg_code) == 4283:

        myax = p.plot(figsize=[10, 10], linewidth=2.0, column=colorby, cmap=my_colormap)  # , marker='o', markersize=10)

        # add colorbar
        divider = make_axes_locatable(myax)
        # pad = separation from figure to colorbar
        cax = divider.append_axes("right", size="3%", pad=0.2)

        fig = myax.get_figure()

        sm = plt.cm.ScalarMappable(cmap=my_colormap)  # , norm=plt.Normalize(vmin=vmin, vmax=vmax))
        # fake up the array of the scalar mappable. Urgh...
        sm._A = p[colorby]  # [1,2,3]

        cb = fig.colorbar(sm, cax=cax, orientation='vertical')
        cb.set_label(colorby, fontdict={'size': 15, 'weight': 'bold'})
        # myax = p.plot(figsize=[10, 8], linewidth=2.0, column='phi_max', cmap='jet')  # , vmin=vmin, vmax=vmax)

        # calculate and set xy limit:
        # myax.set_xlim([140.2, 141.2])  #LieJunWang
        # myax.set_ylim([-20.8, -19.9])

        # myax.set_xlim([140, 150]) # GA-Vic
        # myax.set_ylim([-39, -34])
        #
        # myax.set_xlim([136.7, 137.0])  # 3D_MT_data_
        # myax.set_ylim([-20.65, -20.35])
        #
        # myax.set_xlim([140.0, 144.5])  # WPJ
        # myax.set_ylim([-23.5, -19.0])

        # automatically adjust plot xy-scope
        margin = 0.02  # degree
        margin = 0.05 * (bbox['MaxLon'] - bbox['MinLon'] + bbox['MaxLat'] - bbox['MinLat'])
        myax.set_xlim((bbox['MinLon'] - margin, bbox['MaxLon'] + margin))
        myax.set_ylim((bbox['MinLat'] - margin, bbox['MaxLat'] + margin))

        myax.set_xlabel('Longitude')
        myax.set_ylabel('Latitude')
        myax.set_title(fig_title)
    else:  # UTM kilometer units
        myax = p.plot(figsize=[10, 8], linewidth=2.0, column=colorby,
                      cmap=my_colormap)  # simple plot need to have details added

        myax.set_xlabel('East-West (KM)')
        myax.set_ylabel('North-South (KM)')
        myax.set_title(fig_title)

        # myax.set_xlim([400000, 1300000])
        # myax.set_ylim([5700000, 6200000])
        #
        # myax.set_xlim([400000, 900000])
        # myax.set_ylim([7400000, 7900000])

        # automatically adjust plot xy-scope
        # margin = 2000  # meters
        margin = 0.05 * (bbox['MaxLon'] - bbox['MinLon'] + bbox['MaxLat'] - bbox['MinLat'])
        myax.set_xlim((bbox['MinLon'] - margin, bbox['MaxLon'] + margin))
        myax.set_ylim((bbox['MinLat'] - margin, bbox['MaxLat'] + margin))

        xticks = myax.get_xticks() / 1000
        myax.set_xticklabels(xticks)
        yticks = myax.get_yticks() / 1000
        myax.set_yticklabels(yticks)

        # add colorbar
        divider = make_axes_locatable(myax)
        # pad = separation from figure to colorbar
        cax = divider.append_axes("right", size="3%", pad=0.2)

        fig = myax.get_figure()

        sm = plt.cm.ScalarMappable(cmap=my_colormap)  # , norm=plt.Normalize(vmin=vmin, vmax=vmax))
        # fake up the array of the scalar mappable. Urgh...
        sm._A = p[colorby]  # [1,2,3]

        cb = fig.colorbar(sm, cax=cax, orientation='vertical')
        cb.set_label(colorby, fontdict={'size': 15, 'weight': 'bold'})

    fig = plt.gcf()
    fig.savefig(jpg_file_name, dpi=400)

    if showfig is True:
        plt.show()

    # cleanup memory now
    plt.close()  # this will make prog faster and not too many plot obj kept.
    del (p)
    del (geopdf)
    del (fig)


def process_csv_folder(csv_folder, bbox_dict, target_epsg_code=4283):
    """
    process all *.csv files in a dir, ude target_epsg_code=4283 GDA94 as default.
    This function uses csv-files folder as its input.
    :param csv_folder:
    :return:
    """

    if csv_folder is None:
        _logger.critical("Must provide a csv folder")

    csvfiles = glob.glob(csv_folder + '/*Hz.csv')  # phase_tensor_tipper_0.004578Hz.csv

    # filter the csv files if you do not want to plot all of them
    print(len(csvfiles))

    # for acsv in csvfiles[:2]:
    for acsv in csvfiles:
        tip_re_gdf = create_tipper_real_shp_from_csv(acsv, line_length=0.02, target_epsg_code=target_epsg_code)
        my_gdf = tip_re_gdf
        jpg_file_name = acsv.replace('.csv', '_tip_re_epsg%s.jpg' % target_epsg_code)
        export_geopdf_to_image(my_gdf, bbox_dict, jpg_file_name, target_epsg_code)

        tip_im_gdf = create_tipper_imag_shp_from_csv(acsv, line_length=0.02, target_epsg_code=target_epsg_code)
        my_gdf = tip_im_gdf
        jpg_file_name = acsv.replace('.csv', '_tip_im_epsg%s.jpg' % target_epsg_code)
        export_geopdf_to_image(my_gdf, bbox_dict, jpg_file_name, target_epsg_code)

        ellip_gdf = create_ellipse_shp_from_csv(acsv, esize=0.01, target_epsg_code=target_epsg_code)
        # Now, visualize and output to image file from the geopandas dataframe
        my_gdf = ellip_gdf
        jpg_file_name = acsv.replace('.csv', '_ellips_epsg%s.jpg' % target_epsg_code)
        export_geopdf_to_image(my_gdf, bbox_dict, jpg_file_name, target_epsg_code)

    return



#############################################################################
# ==================================================================
# python mtpy/utils/shapefiles_creator.py data/edifiles /e/tmp
# ==================================================================
#############################################################################

if __name__ == "__main__OLD_V0":

    edidir = sys.argv[1]

    edifiles = glob.glob(os.path.join(edidir, "*.edi"))

    if len(sys.argv) > 2:
        path2out = sys.argv[2]
    else:
        path2out = None

    # filter the edi files here if desired, to get a subset:
    # edifiles2 = edifiles[0:-1:2]
    shp_maker = ShapeFilesCreator(edifiles, path2out)
    # ptdic = shp_maker.create_csv_files_deprecated()  # dest_dir=path2out)    #  create csv files E:/temp1
    ptdic = shp_maker.create_phase_tensor_csv(path2out)  # compare csv in E:/temp2

    # print ptdic
    # print ptdic[ptdic.keys()[0]]

    # edisobj = mtpy.core.edi_collection.EdiCollection(edifiles)
    edisobj = EdiCollection(edifiles)
    bbox_dict = edisobj.bound_box_dict
    print(bbox_dict)

    bbox_dict2 = shp_maker.bound_box_dict
    print(bbox_dict2)
    if bbox_dict != bbox_dict2:
        raise Exception("parent-child's attribute bbo_dic not equal!!!")

    # create shapefiles and plots
    # epsg projection 4283 - gda94
    # process_csv_folder(path2out, bbox_dict)

    # Further testing epsg codes:
    # epsg projection 28354 - gda94 / mga zone 54
    # epsg projection 32754 - wgs84 / utm zone 54s
    #  GDA94/GALCC =3112
    for my_epsgcode in [3112, ]:  # [4326, 4283, 3112, 32755]:   # 32754, 28355]:

        bbox_dict = edisobj.get_bounding_box(epsgcode=my_epsgcode)
        print(bbox_dict)
        process_csv_folder(path2out, bbox_dict, target_epsg_code=my_epsgcode)

###################################################################
# Example codes to use the ShapeFilesCreator class - new version

if __name__ == "__main__d":

    edidir = sys.argv[1]

    edifiles = glob.glob(os.path.join(edidir, "*.edi"))

    if len(sys.argv) > 2:
        path2out = sys.argv[2]
    else:
        path2out = None

    # esize=0.08 # specify ellipse size ?

    # Filter the edi files to get a subset:
    everysite = 1  # every 1,2,3,4, 5
    edi_list = edifiles[::everysite]  # subset of the edi files
    shp_maker = ShapeFilesCreator(edi_list, path2out)

    station_distance_stats = shp_maker.get_stations_distances_stats()

    esize = None  # if None, auto selected default in the method
    tipsize = None  # if None, auto selected default in the method

    _logger.info("User-defined Max-Ellispse Size =:%s", esize)
    _logger.info("User-defined Max-Tipper Length/Size =:%s", tipsize)

    shp_maker.create_phase_tensor_shp(999.99, ellipsize=esize)  # nothing created for non-existent peri

    min_period = shp_maker.all_unique_periods[0]
    max_period = shp_maker.all_unique_periods[-1]

    # for aper in [min_period, max_period]:
    for aper in shp_maker.all_unique_periods[::5]:  # ascending order: from short to long periods
        # default projection as target output
        # shp_maker.create_phase_tensor_shp(2.85)
        # shp_maker.create_phase_tensor_shp(aper,  ellipsize=esize,export_fig=True)
        # shp_maker.create_tipper_real_shp(aper, line_length=tipsize, export_fig=True)
        # shp_maker.create_tipper_imag_shp(aper, line_length=tipsize, export_fig=True)

        for my_epsgcode in [3112]:  # [3112, 4326, 4283, 32754, 32755, 28353, 28354, 28355]:
            shp_maker.create_phase_tensor_shp(aper, target_epsg_code=my_epsgcode, ellipsize=esize, export_fig=True)
            shp_maker.create_tipper_real_shp(aper, line_length=tipsize, target_epsg_code=my_epsgcode, export_fig=True)
            shp_maker.create_tipper_imag_shp(aper, line_length=tipsize, target_epsg_code=my_epsgcode, export_fig=True)


# ===================================================
# Click Command Wrapper for shape files from edi
# ===================================================
@click.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.option('-i', '--input', type=str,
              default='examples/data/edi_files_2', \
              help='input edi files dir ')
@click.option('-c', '--code', type=int, default=3112,
              help='epsg code [3112, 4326, 4283, 32754, 32755, 28353, 28354, 28355]')
@click.option('-o', '--output', type=str, default="temp", help='Output directory')
def generate_shape_files(input, output, code):
    print("=======================================================================")
    print("Generating Shapes File requires following inputs edi files directory   ")
    print("Default epsg code 3112                                                 ")
    print("        epsg_code(4326, 4283, 32754, 32755, 28353, 28354, 28355)       ")
    print("Default output is in temp directory                                    ")
    print("=======================================================================")

    if not os.path.isdir(output):
        os.mkdir(output)
    edifiles = glob.glob(os.path.join(input, "*.edi"))
    # Filter the edi files to get a subset:
    everysite = 1  # every 1,2,3,4, 5
    edi_list = edifiles[::everysite]  # subset of the edi files
    shp_maker = ShapeFilesCreator(edi_list, output)

    # station_distance_stats= shp_maker.get_stations_distances_stats()

    esize = None  # if None, auto selected default in the method
    tipsize = None  # if None, auto selected default in the method

    shp_maker.create_phase_tensor_shp(999.99, ellipsize=esize)  # nothing created for non-existent peri

    # min_period = shp_maker.all_unique_periods[0]
    # max_period = shp_maker.all_unique_periods[-1]

    # for aper in [min_period, max_period]:
    for aper in shp_maker.all_unique_periods[::5]:  # ascending order: from short to long periods
        # default projection as target output
        # shp_maker.create_phase_tensor_shp(2.85)
        # shp_maker.create_phase_tensor_shp(aper,  ellipsize=esize,export_fig=True)
        # shp_maker.create_tipper_real_shp(aper, line_length=tipsize, export_fig=True)
        # shp_maker.create_tipper_imag_shp(aper, line_length=tipsize, export_fig=True)

        for my_epsgcode in [code]:  # [3112, 4326, 4283, 32754, 32755, 28353, 28354, 28355]:
            shp_maker.create_phase_tensor_shp(aper, target_epsg_code=my_epsgcode, ellipsize=esize, export_fig=True)
            shp_maker.create_tipper_real_shp(aper, line_length=tipsize, target_epsg_code=my_epsgcode, export_fig=True)
            shp_maker.create_tipper_imag_shp(aper, line_length=tipsize, target_epsg_code=my_epsgcode, export_fig=True)


if __name__ == "__main__":

    print("Please see examples/scripts/create_pt_shapefiles.py")

  # generate_shape_files()  # click CLI interface
