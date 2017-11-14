#! /usr/bin/env python
"""
Description:

    This script creates shape files for MT datasets.
    Phase Tensor, Tipper Real/Imag, MT-site locations, etc


CreationDate:   2017-03-06
Developer:      fei.zhang@ga.gov.au

Revision History:
    LastUpdate:     10/11/2017   FZ fix bug after the big merge

"""

from __future__ import print_function

import csv
import glob
import logging
import os
import sys

import geopandas as gpd
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from mpl_toolkits.axes_grid1 import make_axes_locatable
from shapely.geometry import Point, Polygon, LineString, LinearRing

from mtpy.core.edi_collection import EdiCollection
import mtpy.core.mt as mt
from mtpy.utils.decorator import deprecated
from mtpy.utils.mtpylog import MtPyLog

mpl.rcParams['lines.linewidth'] = 2
# mpl.rcParams['lines.color'] = 'r'
mpl.rcParams['figure.figsize'] = [10, 6]

logger = MtPyLog().get_mtpy_logger(__name__)
logger.setLevel(logging.DEBUG)


class ShapeFilesCreator(EdiCollection):
    """ Extend the functions of the EdiCollection class,
    create phase tensor and tipper shapefiles for a list of edifiles
    """

    def __init__(self, edifile_list, outdir):
        """
        loop through a list of edi files, create required shapefiles
        :param edifile_list: [path2edi,...]
        :param outdir: path2output dir, where the shpe file weill be written.
        """

        super(ShapeFilesCreator,self).__init__(edilist=edifile_list, outdir=outdir)
        #python-3 syntax: super().__init__(edilist=edifile_list, outdir=outdir)

        self.outputdir = outdir
        self.all_periods = self._get_all_periods()
        self.ptol = 0.05  # this param controls what freqs/periods are grouped together:
        # 10% may result more double counting of freq/period data than 5%.
        # eg: E:\Data\MT_Datasets\WenPingJiang_EDI 18528 rows vs 14654 rows

        return

    def _get_all_periods(self):
        """
        from the list of edi files get a list of all possible frequencies.

        """
        if self.all_frequencies is not None:  # already initialized
            return

        # get all frequencies from all edi files
        all_freqs = []
        for mt_obj in self.mt_obj_list:
            all_freqs.extend(list(mt_obj.Z.freq))

        # sort all frequencies so that they are in ascending order,
        # use set to remove repeats and make an array
        self.all_frequencies = sorted(list(set(all_freqs)))

        logger.info("Number of MT Frequencies: %s", len(self.all_frequencies))

        all_periods = 1.0 / np.array(sorted(self.all_frequencies, reverse=True))

        logger.debug("Type of all_periods %s", type(all_periods))
        logger.info("Number of MT Periods: %s", len(all_periods))
        logger.debug(all_periods)

        # else:
        #     if isinstance(self.all_periods, list):
        #         pass
        #     if isinstance(self.all_periods, int) or isinstance(
        #             self.all_periods, float):
        #         self.all_periods = [self.all_periods]

        return all_periods

    @deprecated("This function is replaced by mtpy.core.edi_collection.EdiCollection.create_phase_tensor_csv()")
    def create_csv_files(self, dest_dir=None):
        """
        create csv file
        :return:
        """
        if dest_dir is None:
            dest_dir = self.outputdir

        # summary csv file
        csvfname = os.path.join(dest_dir, "phase_tensor.csv")

        pt_dict = {}

        csv_header = ['station', 'freq', 'lon', 'lat', 'phi_min', 'phi_max', 'azimuth', 'skew', 'n_skew', 'elliptic',
                      'tip_mag_re', 'tip_mag_im', 'tip_ang_re', 'tip_ang_im']

        with open(csvfname, "wb") as csvf:
            writer = csv.writer(csvf)
            writer.writerow(csv_header)

        for freq in self.all_frequencies:
            ptlist = []
            for mt_obj in self.mt_obj_list:

                f_index_list = [ff for ff, f2 in enumerate(mt_obj.Z.freq)
                                if (f2 > freq * (1 - self.ptol)) and
                                (f2 < freq * (1 + self.ptol))]
                if len(f_index_list) > 1:
                    logger.warn("more than one fre found %s", f_index_list)

                if len(f_index_list) >= 1:
                    p_index = f_index_list[0]
                    # geographic coord lat long and elevation
                    # long, lat, elev = (mt_obj.lon, mt_obj.lat, 0)
                    station, lon, lat = (mt_obj.station, mt_obj.lon, mt_obj.lat)

                    # print(type(mt_obj.pt.phimin))
                    # print(mt_obj.pt.phimin[p_index])

                    pt_stat = [station, freq, lon, lat,
                               mt_obj.pt.phimin[p_index],
                               mt_obj.pt.phimax[p_index],
                               mt_obj.pt.azimuth[p_index],
                               mt_obj.pt.beta[p_index],
                               2 * mt_obj.pt.beta[p_index],
                               mt_obj.pt.ellipticity[p_index],  # FZ: get ellipticity begin here
                               mt_obj.Tipper.mag_real[p_index],
                               mt_obj.Tipper.mag_imag[p_index],
                               mt_obj.Tipper.angle_real[p_index],
                               mt_obj.Tipper.angle_imag[p_index]]

                    ptlist.append(pt_stat)
                else:
                    logger.warn('Freq %s NOT found for this station %s', freq, mt_obj.station)

            with open(csvfname, "ab") as csvf:  # summary csv for all freqs
                writer = csv.writer(csvf)
                writer.writerows(ptlist)

            csvfile2 = csvfname.replace('.csv', '_%sHz.csv' % str(freq))

            with open(csvfile2, "wb") as csvf:  # individual csvfile for each freq
                writer = csv.writer(csvf)

                writer.writerow(csv_header)
                writer.writerows(ptlist)

            pt_dict[freq] = ptlist

        return pt_dict


    def create_phase_tensor_shp(self):
        """
        create phase tensor ellipses shape files
        :return:
        """

        pass

    def create_tipper_shp(self):
        """
        create tipper shape files
        :return:
        """

        pass

    def create_mt_sites_shp(self):
        """
        create MT site location points shape files
        :return:
        """

        pass


# http://toblerity.org/shapely/manual.html#polygons
# https://geohackweek.github.io/vector/04-geopandas-intro/
def create_ellipse_shp(csvfile, esize=0.03, target_epsg_code=None):
    """
    create phase tensor ellipse geometry from a csv file
    :param csvfile: a csvfile with full path
    :param esize: ellipse size, defaut 0.03 is about 3KM in the max ellipse rad
    :return: a geopandas dataframe
    """

    pdf = pd.read_csv(csvfile)
    mt_locations = [Point(xy) for xy in zip(pdf['lon'], pdf['lat'])]
    # OR pdf['geometry'] = pdf.apply(lambda z: Point(z.lon, z.lat), axis=1)
    # if you want to df = df.drop(['Lon', 'Lat'], axis=1)
    crs = {'init': 'epsg:4326'}  # initial crs WGS84

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
    for i in xrange(0, len(azimuth)):
        x = x0[i] + height[i] * np.cos(theta) * np.cos(azimuth[i]) - \
            width[i] * np.sin(theta) * np.sin(azimuth[i])
        y = y0[i] + height[i] * np.cos(theta) * np.sin(azimuth[i]) + \
            width[i] * np.sin(theta) * np.cos(azimuth[i])

        polyg = Polygon(LinearRing([xy for xy in zip(x, y)]))

        # print polyg  # an ellispe

        ellipse_list.append(polyg)

    pdf = gpd.GeoDataFrame(pdf, crs=crs, geometry=ellipse_list)

    if target_epsg_code is None:
        target_epsg_code = '4326'  # EDI original lat/lon epsg 4326 or GDA94
    else:
        pdf.to_crs(epsg=target_epsg_code, inplace=True)
        # world = world.to_crs({'init': 'epsg:3395'})
        # world.to_crs(epsg=3395) would also work

    # to shape file
    shp_fname = csvfile.replace('.csv', '_ellip_epsg%s.shp' % target_epsg_code)
    pdf.to_file(shp_fname, driver='ESRI Shapefile')

    return pdf


def plot_geopdf(pdf, bbox, out_file_name, target_epsg_code, showfig=False):
    if target_epsg_code is None:
        p = pdf
        target_epsg_code = '4326'  # EDI orginal lat/lon epsg 4326 or GDA94
    else:
        p = pdf.to_crs(epsg=target_epsg_code)
        # world = world.to_crs({'init': 'epsg:3395'})
        # world.to_crs(epsg=3395) would also work

    # bounds = p.total_bounds  # lat-lon bounds for this csv dataframe

    # plot and save
    jpg_fname = out_file_name.replace('.csv', '_epsg%s.jpg' % target_epsg_code)
    fig_title = os.path.basename(jpg_fname)
    logger.info('saving figure to file %s', jpg_fname)

    colorby = 'phi_min'
    my_cmap_r = 'jet'

    if int(target_epsg_code) == 4326:

        myax = p.plot(figsize=[10, 10], linewidth=2.0, column=colorby, cmap=my_cmap_r)  # , marker='o', markersize=10)

        # add colorbar
        divider = make_axes_locatable(myax)
        # pad = separation from figure to colorbar
        cax = divider.append_axes("right", size="3%", pad=0.2)

        fig = myax.get_figure()

        sm = plt.cm.ScalarMappable(cmap=my_cmap_r)  # , norm=plt.Normalize(vmin=vmin, vmax=vmax))
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
        myax.set_xlim((bbox['MinLon'] - margin, bbox['MaxLon'] + margin))
        myax.set_ylim((bbox['MinLat'] - margin, bbox['MaxLat'] + margin))

        myax.set_xlabel('Longitude')
        myax.set_ylabel('Latitude')
        myax.set_title(fig_title)
    else:
        myax = p.plot(figsize=[10, 8], linewidth=2.0, column=colorby,
                      cmap=my_cmap_r)  # simple plot need to have details added

        myax.set_xlabel('East-West (KM)')
        myax.set_ylabel('North-South (KM)')
        myax.set_title(fig_title)

        # myax.set_xlim([400000, 1300000])
        # myax.set_ylim([5700000, 6200000])
        #
        # myax.set_xlim([400000, 900000])
        # myax.set_ylim([7400000, 7900000])

        # automatically adjust plot xy-scope
        margin = 2000  # meters
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

        sm = plt.cm.ScalarMappable(cmap=my_cmap_r)  # , norm=plt.Normalize(vmin=vmin, vmax=vmax))
        # fake up the array of the scalar mappable. Urgh...
        sm._A = p[colorby]  # [1,2,3]

        cb = fig.colorbar(sm, cax=cax, orientation='vertical')
        cb.set_label(colorby, fontdict={'size': 15, 'weight': 'bold'})

    fig = plt.gcf()
    fig.savefig(jpg_fname, dpi=400)

    if showfig is True:
        plt.show()

    # cleanup memory now
    plt.close()  # this will make prog faster and not too many plot obj kept.
    del (p)
    del (pdf)
    del (fig)


def create_tipper_real_shp(csvfile, arr_size=0.03, target_epsg_code=None):
    """ create tipper lines shape from a csv file
    The shape is a line without arrow.
    Must use a GIS software such as ArcGIS to display and add an arrow at each line's end
    arr_size=4  how long will be the line (arrow)
    return: a geopandas dataframe object for further processing.
    """

    pdf = pd.read_csv(csvfile)
    # mt_locations = [Point(xy) for xy in zip(pdf.lon, pdf.lat)]
    # OR pdf['geometry'] = pdf.apply(lambda z: Point(z.lon, z.lat), axis=1)
    # if you want to df = df.drop(['Lon', 'Lat'], axis=1)

    crs = {'init': 'epsg:4326'}  # WGS84

    # geo_df = gpd.GeoDataFrame(pdf, crs=crs, geometry=mt_locations)

    pdf['tip_re'] = pdf.apply(lambda x:
                              LineString([(float(x.lon), float(x.lat)),
                                          (float(x.lon) + arr_size * x.tip_mag_re * np.cos(
                                              -np.deg2rad(x.tip_ang_re)),
                                           float(x.lat) + arr_size * x.tip_mag_re * np.sin(
                                               -np.deg2rad(x.tip_ang_re)))]), axis=1)

    pdf = gpd.GeoDataFrame(pdf, crs=crs, geometry='tip_re')

    if target_epsg_code is None:
        target_epsg_code = '4326'  # EDI original lat/lon epsg 4326 or GDA94
    else:
        pdf.to_crs(epsg=target_epsg_code, inplace=True)
        # world = world.to_crs({'init': 'epsg:3395'})
        # world.to_crs(epsg=3395) would also work

    # to shape file
    shp_fname = csvfile.replace('.csv', '_real_epsg%s.shp' % target_epsg_code)
    pdf.to_file(shp_fname, driver='ESRI Shapefile')

    return pdf


def create_tipper_imag_shp(csvfile, arr_size=0.03, target_epsg_code=None):
    """ create imagery tipper lines shape from a csv file
    The shape is a line without arrow.
    Must use a GIS software such as ArcGIS to display and add an arrow at each line's end
    arr_size=4  how long will be the line (arrow)
    return: a geopandas dataframe object for further processing.
    """

    pdf = pd.read_csv(csvfile)
    # mt_locations = [Point(xy) for xy in zip(pdf.lon, pdf.lat)]
    # OR pdf['geometry'] = pdf.apply(lambda z: Point(z.lon, z.lat), axis=1)
    # if you want to df = df.drop(['Lon', 'Lat'], axis=1)

    crs = {'init': 'epsg:4326'}  # WGS84

    # geo_df = gpd.GeoDataFrame(pdf, crs=crs, geometry=mt_locations)

    pdf['tip_im'] = pdf.apply(lambda x:
                              LineString([(float(x.lon), float(x.lat)),
                                          (float(x.lon) + arr_size * x.tip_mag_im * np.cos(
                                              -np.deg2rad(x.tip_ang_im)),
                                           float(x.lat) + arr_size * x.tip_mag_im * np.sin(
                                               -np.deg2rad(x.tip_ang_im)))]), axis=1)

    pdf = gpd.GeoDataFrame(pdf, crs=crs, geometry='tip_im')

    if target_epsg_code is None:
        target_epsg_code = '4326'  # EDI original lat/lon epsg 4326 or GDA94
    else:
        pdf.to_crs(epsg=target_epsg_code, inplace=True)
        # world = world.to_crs({'init': 'epsg:3395'})
        # world.to_crs(epsg=3395) would also work

    # to shape file
    shp_fname = csvfile.replace('.csv', '_imag_epsg%s.shp' % target_epsg_code)
    pdf.to_file(shp_fname, driver='ESRI Shapefile')

    return pdf


def process_csv_folder(csv_folder, bbox_dict, target_epsg_code=None):
    """
    process all *.csv files in a dir
    :param csv_folder:
    :return:
    """

    if csv_folder is None:
        logger.critical("Must provide a csv folder")

    csvfiles = glob.glob(csv_folder + '/*Hz.csv')  # phase_tensor_tipper_0.004578Hz.csv

    # filter the csv files if you do not want to plot all of them
    print(len(csvfiles))

    # for acsv in csvfiles[:2]:
    for acsv in csvfiles:
        tip_re_gdf = create_tipper_real_shp(acsv, arr_size=0.02, target_epsg_code=target_epsg_code)

        tip_im_gdf = create_tipper_imag_shp(acsv, arr_size=0.02, target_epsg_code=target_epsg_code)

        ellip_gdf = create_ellipse_shp(acsv, esize=0.01, target_epsg_code=target_epsg_code)

        # Now, visualize and output to image file from the geopandas dataframe
        my_gdf = ellip_gdf

        plot_geopdf(my_gdf, bbox_dict, acsv, target_epsg_code)

    return


# ==================================================================
# python mtpy/utils/shapefiles_creator.py data/edifiles /e/tmp
# ==================================================================
if __name__ == "__main__":

    edidir = sys.argv[1]

    edifiles = glob.glob(os.path.join(edidir, "*.edi"))

    if len(sys.argv) > 2:
        path2out = sys.argv[2]
    else:
        path2out = None

    # filter the edi files here if desired, to get a subset:
    # edifiles2 = edifiles[0:-1:2]
    shp_maker = ShapeFilesCreator(edifiles, path2out)
    ptdic = shp_maker.create_csv_files()  # dest_dir=path2out)    #  create csv files E:/temp1
    #use super class: ptdic =shp_maker.create_phase_tensor_csv(path2out)  # compare csv in E:/temp2

    # print ptdic
    # print ptdic[ptdic.keys()[0]]

    # edisobj = mtpy.core.edi_collection.EdiCollection(edifiles)
    edisobj = EdiCollection(edifiles)
    bbox_dict = edisobj.bound_box_dict
    print(bbox_dict)


    # create shapefiles and plots
    # epsg projection 4283 - gda94
    process_csv_folder(path2out, bbox_dict)  # , target_epsg_code will be default 4326?

    # epsg projection 28354 - gda94 / mga zone 54
    # epsg projection 32754 - wgs84 / utm zone 54s

    # for my_epsgcode in [32754,]:  #[3112, 32755]:   # 32754, 28355]:
    #     #my_epsgcode = 32755 # GDA94/GALCC =3112
    #     bbox_dict=edisobj.get_bounding_box(epsgcode=my_epsgcode)
    #     print(bbox_dict)
    #     process_csv_folder(path2out, bbox_dict, target_epsg_code=my_epsgcode)
