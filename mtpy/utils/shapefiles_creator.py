"""
create shape files for MT datasets.
Phase Tensor, Tipper Real/Imag, MT-site locations,etc

fei.zhang@ga.gov.au
2017-03-06
"""
import glob
import os
import sys
import csv

import geopandas as gpd
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from shapely.geometry import Point, Polygon, LinearRing

import mtpy.core.mt as mt
from mtpy.imaging.phase_tensor_maps import PlotPhaseTensorMaps
from mtpy.utils.mtpylog import MtPyLog

mpl.rcParams['lines.linewidth']=2
# mpl.rcParams['lines.color'] = 'r'
mpl.rcParams['figure.figsize']=[10, 6]

logger=MtPyLog().get_mtpy_logger(__name__)


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

        self.edifiles=edifile_list
        logger.info("number of edi files to be processed: %s",
                    len(self.edifiles))
        assert len(self.edifiles) > 0

        self.outputdir=outdir

        if self.edifiles is not None:
            self.mt_obj_list=[mt.MT(edi) for edi in self.edifiles]

        # get all frequencies from all edi files
        self.all_periods=self._get_all_periods()

        self.ptol=0.05

        return

    def _get_all_periods(self):
        """
        from the list of edi files get a list of all possible frequencies.

        """

        # get all frequencies from all edi files
        all_freqs=[]
        for mt_obj in self.mt_obj_list:
            all_freqs.extend(list(mt_obj.Z.freq))

        # sort all frequencies so that they are in descending order,
        # use set to remove repeats and make an array
        self.all_frequencies=sorted(list(set(all_freqs)))

        all_periods=1.0 / \
                    np.array(sorted(list(set(all_freqs)), reverse=True))
        logger.debug(all_periods)

        logger.info("Number of MT Periods: %s", len(all_periods))

        # else:
        #     if isinstance(self.all_periods, list):
        #         pass
        #     if isinstance(self.all_periods, int) or isinstance(
        #             self.all_periods, float):
        #         self.all_periods = [self.all_periods]

        return all_periods

    def create_csv2(self, csvfname="E:/tmp/test.csv"):
        """
        create csv from the shapefiles.py method
        :return:
        """

        self.pt_dict={}

        csv_header=['station','freq', 'lon', 'lat', 'phi_min', 'phi_max','azimuth', 'skew', 'normskew','elliptic',
                    'tip_mag_re',  'tip_mag_im', 'tip_ang_re', 'tip_ang_im']
        with open(csvfname, "wb") as csvf:
            writer=csv.writer(csvf)
            writer.writerow(csv_header)

        for freq in self.all_frequencies:
            ptlist=[]
            for mt_obj in self.mt_obj_list:

                f_index_list=[ff for ff, f2 in enumerate( mt_obj.Z.freq)
                              if (f2 > freq * (1 - self.ptol)) and
                              (f2 < freq * (1 + self.ptol))]
                if len(f_index_list) > 1:
                    logger.warn("more than one fre found %s", f_index_list)

                p_index=f_index_list[0]
                # geographic coord lat long and elevation
                long, lat, elev=(mt_obj.lon, mt_obj.lat, 0)

                pt_stat=[mt_obj.station, freq, long, lat,
                          mt_obj.pt.phimin[0][p_index],
                          mt_obj.pt.phimax[0][p_index],
                          mt_obj.pt.azimuth[0][p_index],
                          mt_obj.pt.beta[0][p_index],
                          2 * mt_obj.pt.beta[0][p_index],
                          mt_obj.pt.ellipticity[0][p_index], # FZ: get ellipticity begin here
                         mt_obj.Tipper.mag_real[p_index],
                         mt_obj.Tipper.mag_imag[p_index],
                         mt_obj.Tipper.angle_real[p_index],
                         mt_obj.Tipper.angle_imag[p_index] ]

                ptlist.append(pt_stat)

            with open(csvfname, "ab") as csvf:
                writer = csv.writer(csvf)
                writer.writerows(ptlist)

            self.pt_dict[freq] = ptlist

        return self.pt_dict

    def create_csv(self):
        """
        generate a csv db for MT sites, phase tensor attributes, etc, from edi files
        :return:
        """

        for freq in self.all_frequencies:
            ptm=PlotPhaseTensorMaps(fn_list=self.edifiles, plot_freq=freq)

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
    """ create phase tensor ellipse
    esize is ellipse size, defaut 0.03 is about 3KM in the max ellipse rad
    """

    pdf=pd.read_csv(csvfile)
    mt_locations=[Point(xy) for xy in zip(pdf['lon'], pdf['lat'])]
    # OR pdf['geometry'] = pdf.apply(lambda z: Point(z.lon, z.lat), axis=1)
    # if you want to df = df.drop(['Lon', 'Lat'], axis=1)
    crs={'init': 'epsg:4326'}  # initial crs WGS84

    pdf=gpd.GeoDataFrame(pdf, crs=crs, geometry=mt_locations)

    # make  pt_ellispes using polygons
    phi_max_v=pdf['phi_max'].max()  # the max of this group of ellipse

    print phi_max_v

    # points to trace out the polygon-ellipse
    theta=np.arange(0, 2 * np.pi, np.pi / 30.)

    azimuth=-np.deg2rad(pdf['azimuth'])
    width=esize * (pdf['phi_max'] / phi_max_v)
    height=esize * (pdf['phi_min'] / phi_max_v)
    x0=pdf['lon']
    y0=pdf['lat']

    # apply formula to generate ellipses

    ellipse_list=[]
    for i in xrange(0, len(azimuth)):
        x=x0[i] + height[i] * np.cos(theta) * np.cos(azimuth[i]) - \
          width[i] * np.sin(theta) * np.sin(azimuth[i])
        y=y0[i] + height[i] * np.cos(theta) * np.sin(azimuth[i]) + \
          width[i] * np.sin(theta) * np.cos(azimuth[i])

        polyg=Polygon(LinearRing([xy for xy in zip(x, y)]))

        # print polyg  # an ellispe

        ellipse_list.append(polyg)

    pdf=gpd.GeoDataFrame(pdf, crs=crs, geometry=ellipse_list)

    return pdf


def process_csv_folder(csv_folder):
    """
    process all *.csv files in a dir
    :param csv_folder:
    :return:
    """

    csvfiles=glob.glob(csv_folder + '/*.csv')

    print (len(csvfiles))

    for acsv in csvfiles[:2]:
        # for acsv in csvfiles:
        p=get_geopdf_from_csv(acsv)

        # plot and save
        p.plot()
        fig=plt.gcf()
        jpg_fname=acsv.replace('.csv', '.jpg')

        print (jpg_fname)
        fig.savefig(jpg_fname, dpi=300)
        plt.show()

        # to shape file
        # shp_fname=acsv.replace('.csv', '.shp')
        # p.to_file(shp_fname, driver='ESRI Shapefile')


# ==================================================================
if __name__ == "__main__":
    edidir=sys.argv[1]

    edifiles=glob.glob(os.path.join(edidir, "*.edi"))

    path2out=sys.argv[2]

    shp_maker=ShapeFilesCreator(edifiles, path2out)

    ptdic=shp_maker.create_csv2()

    print ptdic
    #print ptdic[ptdic.keys()[0]]
    # shp_maker.create_mt_sites_shp()

    # process_csv_folder(r'E:/Data/MT_Datasets/WenPingJiang_SHP/')
