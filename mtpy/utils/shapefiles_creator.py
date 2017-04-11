"""
create shape files for MT datasets.
Phase Tensor, Tipper Real/Imag, MT-site locations,etc

fei.zhang@ga.gov.au
2017-03-06
"""
import os
import sys
import glob

import numpy as np
import geopandas

import mtpy.analysis.pt as mtpt
import mtpy.core.mt as mt
from mtpy.imaging.phase_tensor_maps import PlotPhaseTensorMaps

class ShapeFilesCreator:

    def __init__(self, edifile_list, outdir):
        """
        loop through a list of edi files, create required shapefiles
        :param edifile_list: [path2edi,...]
        :param outdir: path2output dir, where the shpe file weill be written.
        """

        self.edifiles = edifile_list
        print ("number of edi files to be processed:", len(self.edifiles))
        assert len(self.edifiles) > 0

        self.outputdir = outdir

        if self.edifiles is not None:
            self.mt_obj_list = [mt.MT(edi) for edi in self.edifiles]


        # get all frequencies from all edi files
        all_freqs=[]
        for mt_obj in self.mt_obj_list:
            all_freqs.extend(list(mt_obj.Z.freq))

        # sort all frequencies so that they are in descending order,
        # use set to remove repeats and make an array
        self.all_frequencies = sorted(list(set(all_freqs)))

        self.all_periods=1. / np.array(sorted(list(set(all_freqs)),
                                              reverse=True))
        print self.all_periods

        print ("Number of MT Periods:", len(self.all_periods))

        return


    def create_csv(self):
        """
        generate a csv db for MT sites, phase tensor attributes, etc, from edi files
        :return:
        """

        for per in self.all_frequencies:
            ptm=PlotPhaseTensorMaps(fn_list=self.edifiles, plot_freq=per)

            ptm.export_params_to_file(save_path=self.outputdir)

        pass


    def create_phase_tensor_shp(self):

        pass

    def create_tipper_shp(self):

        pass

    def create_MTsites_shp(self):

        pass

    # http://toblerity.org/shapely/manual.html#polygons
    # https://geohackweek.github.io/vector/04-geopandas-intro/

import os, sys,glob
import numpy as np
import pandas as pd
import geopandas as gpd
from shapely.geometry import Point, Polygon, LineString, LinearRing
import matplotlib.pyplot as plt

#
# import folium
# from IPython.display import display
#
# from shapely.geometry import mapping

def create_phase_tensor_ellipse_shp(csvfile, esize=0.03):
    """ create phase tensor ellipse
    esize is ellipse size, defaut 0.03 is about 3KM in the max ellipse rad
    """

    pdf=pd.read_csv(csvfile)
    mt_locations=[Point(xy) for xy in zip(pdf.lon, pdf.lat)]
    # OR pdf['geometry'] = pdf.apply(lambda z: Point(z.lon, z.lat), axis=1)
    # if you want to df = df.drop(['Lon', 'Lat'], axis=1)
    crs={'init': 'epsg:4326'}  #initial crs WGS84

    pdf=gpd.GeoDataFrame(pdf, crs=crs, geometry=mt_locations)

    # make  pt_ellispes using polygons
    PHIMAX=pdf['phi_max'].max()  # the max of this group of ellipse

    print PHIMAX

    theta=np.arange(0, 2 * np.pi, np.pi / 30.)  # points to trace out the polygon-ellipse

    azimuth=-np.deg2rad(pdf['azimuth'])
    width=esize * (pdf['phi_max'] / PHIMAX)
    height=esize * (pdf['phi_min'] / PHIMAX)
    x0=pdf['lon']
    y0=pdf['lat']

    # apply formula to generate ellipses

    ellipse_list=[]
    for i in xrange(0, len(azimuth)):
        x=x0[i] + height[i] * np.cos(theta) * np.cos(azimuth[i]) - width[i] * np.sin(theta) * np.sin(azimuth[i])
        y=y0[i] + height[i] * np.cos(theta) * np.sin(azimuth[i]) + width[i] * np.sin(theta) * np.cos(azimuth[i])

        polyg=Polygon(LinearRing([xy for xy in zip(x, y)]))

        # print polyg  # an ellispe

        ellipse_list.append(polyg)


        #     for xi, yi in zip(x, y):
        #         polyg.(np.round(xi, 6), np.round(yi, 6))


        #                     # 1) make a geometry shape of the ellipse
        #                     ellipse = ogr.Geometry(ogr.wkbLinearRing)
        #                     ellipse.CloseRings()

        #                     # 2) make a polygon
        #                     poly = ogr.Geometry(ogr.wkbPolygon)
        #                     poly.AddGeometry(ellipse)

        #                     poly_list.append(poly)

    pdf=gpd.GeoDataFrame(pdf, crs=crs, geometry=ellipse_list)

    return pdf


# ==================================================================
if __name__ == "__main__":

    edidir = sys.argv[1]

    edifiles = glob.glob(os.path.join(edidir, "*.edi"))

    path2out = sys.argv[2]

    shp_maker = ShapeFilesCreator(edifiles, path2out)

    #shp_maker.create_csv()

    #shp_maker.create_MTsites_shp()

    # process all *.csv files in a dir
    CSVDIR='E:/Data/MT_Datasets/WenPingJiang_SHP/'
    csvfiles=glob.glob(CSVDIR + '/*.csv')

    print (len(csvfiles))

    for acsv in csvfiles:
        p=create_phase_tensor_ellipse_shp(acsv)

        p.plot()
        plt.show()
        # shp_fname=acsv.replace('.csv', '.shp')
        # p.to_file(shp_fname, driver='ESRI Shapefile')
