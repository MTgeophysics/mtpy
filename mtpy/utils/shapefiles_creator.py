"""
create shape files for MT datasets.
Phase Tensor, Tipper Real/Imag, MT-site locations,etc

fei.zhang@ga.gov.au
2017-03-06
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
import mtpy.core.mt as mt
import numpy as np
import pandas as pd
from mtpy.imaging.phase_tensor_maps import PlotPhaseTensorMaps
from mtpy.utils.mtpylog import MtPyLog
from shapely.geometry import Point, Polygon, LinearRing

from mpl_toolkits.axes_grid1 import make_axes_locatable

mpl.rcParams['lines.linewidth'] = 2
# mpl.rcParams['lines.color'] = 'r'
mpl.rcParams['figure.figsize'] = [10, 6]

logger = MtPyLog().get_mtpy_logger(__name__)
logger.setLevel(logging.DEBUG)


class ShapeFilesCreator(object):
    """
    create shape files for a list of MT edifiles
    """

    def __init__(self, edifile_list, outdir):
        """
        loop through a list of edi files, create required shapefiles
        :param edifile_list: [path2edi,...]
        :param outdir: path2output dir, where the shpe file weill be written.
        """

        self.edifiles = edifile_list
        logger.info("number of edi files to be processed: %s",
                    len(self.edifiles))
        assert len(self.edifiles) > 0

        self.outputdir = outdir

        if self.edifiles is not None:
            self.mt_obj_list = [mt.MT(edi) for edi in self.edifiles]

        # get all frequencies from all edi files
        self.all_frequencies = None
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

    def create_csv_files(self, dest_dir=None):
        """
        create csv from the shapefiles.py method
        :return:
        """
        if dest_dir is None:
            dest_dir = self.outputdir

        # summary csv file
        csvfname = os.path.join(dest_dir, "phase_tensor_tipper.csv")

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

                    pt_stat = [station, freq, lon, lat,
                               mt_obj.pt.phimin[0][p_index],
                               mt_obj.pt.phimax[0][p_index],
                               mt_obj.pt.azimuth[0][p_index],
                               mt_obj.pt.beta[0][p_index],
                               2 * mt_obj.pt.beta[0][p_index],
                               mt_obj.pt.ellipticity[0][p_index],  # FZ: get ellipticity begin here
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

    def create_csv_files2(self):
        """
        Using PlotPhaseTensorMaps class to generate csv file of phase tensor attributes, etc.
        Only for comparison. This method is more expensive because it will create plot object first.
        :return:
        """

        for freq in self.all_frequencies:
            ptm = PlotPhaseTensorMaps(fn_list=self.edifiles, plot_freq=freq)

            ptm.export_params_to_file(save_path=self.outputdir)

        return

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
def get_geopdf_from_csv(csvfile, esize=0.03):
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

    return pdf


def process_csv_folder(csv_folder, target_epsg_code=None):
    """
    process all *.csv files in a dir
    :param csv_folder:
    :return:
    """

    if csv_folder is None:
        logger.critical("Must provide a csv folder")

    csvfiles = glob.glob(csv_folder + '/*Hz.csv')  # phase_tensor_tipper_0.004578Hz.csv

    print(len(csvfiles))

    #for acsv in csvfiles[:2]:
    for acsv in csvfiles:
        p0 = get_geopdf_from_csv(acsv, esize=0.1)

        if target_epsg_code is None:
            p = p0
            target_epsg_code = '4326'  # EDI orginal lat/lon epsg 4326 or GDA94
        else:
            p = p0.to_crs(epsg=target_epsg_code)
            # world = world.to_crs({'init': 'epsg:3395'})
            # world.to_crs(epsg=3395) would also work

        bounds = p.total_bounds  # lat-lon bounds for this csv/pdf
        # plot and save
        jpg_fname = acsv.replace('.csv', '_epsg%s.jpg' % target_epsg_code)
        fig_title=os.path.basename(jpg_fname)
        logger.info('saving figure to file %s',jpg_fname)

        if int(target_epsg_code) == 4326:

            colorby = 'phi_min'
            my_cmap_r = 'jet'

            myax = p.plot(figsize=[10, 10], linewidth=2.0, column=colorby, cmap=my_cmap_r) #, marker='o', markersize=10)

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
            # myax.set_xlim([140.2, 141.2])
            # myax.set_ylim([-20.8, -19.9])

            # margin = 0.0
            #
            # xmin = (1-margin) * bounds[0]
            # xmax = (1+margin) * bounds[2]
            # ymin = (1-margin) * bounds[1]
            # ymax = (1+margin) * bounds[3]
            # myax.set_xlim([xmin, xmax])
            # myax.set_ylim([ymin, ymax])


            myax.set_xlim([140, 150])
            myax.set_ylim([-39, -34])

            myax.set_xlabel('Longitude')
            myax.set_ylabel('Latitude')
            myax.set_title(fig_title)
        else:
            myax = p.plot(figsize=[10, 8], linewidth=2.0, column='phi_max', cmap='jet')  # simple plot need to have details added
            myax.set_xlabel('East-West (KM)')
            myax.set_ylabel('North-South (KM)')
            myax.set_title(fig_title)
            myax.set_xlim([400000, 1300000])
            myax.set_ylim([5700000, 6200000])

            xticks = myax.get_xticks() / 1000
            myax.set_xticklabels(xticks)
            yticks = myax.get_yticks() / 1000
            myax.set_yticklabels(yticks)

        fig = plt.gcf()
        fig.savefig(jpg_fname, dpi=400)
        #plt.show()

        plt.close() # this will make prog faster and no too many plot obj kept.

        # to shape file
        shp_fname = acsv.replace('.csv', '_epsg%s.shp' % target_epsg_code)
        p.to_file(shp_fname, driver='ESRI Shapefile')

        # cleanup memory now
        del(p)
        del(p0)
        del(fig)


# ==================================================================
if __name__ == "__main__":

    edidir = sys.argv[1]

    edifiles = glob.glob(os.path.join(edidir, "*.edi"))

    if len(sys.argv) > 2:
        path2out = sys.argv[2]
    else:
        path2out = None

    # shp_maker = ShapeFilesCreator(edifiles, path2out)
    # # create csv files
    # ptdic = shp_maker.create_csv_files()  # dest_dir=path2out)

    # print ptdic
    # print ptdic[ptdic.keys()[0]]

    # shp_maker.create_mt_sites_shp()

    # create shapefiles and plots
    # epsg projection 4283 - gda94
    process_csv_folder(path2out)  # , target_epsg_code will be default 4326?

    # epsg projection 28354 - gda94 / mga zone 54
    # epsg projection 32754 - wgs84 / utm zone 54s
    #process_csv_folder(path2out, target_epsg_code=32754)
    # process_csv_folder(path2out, target_epsg_code=3112)
