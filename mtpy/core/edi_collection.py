"""
Description:
To compute and encapsulate the properties of a set of EDI files

Author: fei.zhang@ga.gov.au

InitDate: 2017-04-20
"""

from __future__ import print_function
import sys
import os
import glob
import logging
import csv
import numpy as np
import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt
from shapely.geometry import Point #, Polygon, LineString, LinearRing
# import matplotlib as mpl
# from mpl_toolkits.axes_grid1 import make_axes_locatable
import mtpy.core.mt as mt
import mtpy.imaging.mtplottools as mtplottools

from mtpy.utils.mtpylog import MtPyLog

logger = MtPyLog().get_mtpy_logger(__name__)
logger.setLevel(logging.DEBUG)


def is_num_in_seq(anum, aseq, tolerance=0.0000001):
    """
    check if anum is in a sequence by a small tolerance
    :param anum:
    :param aseq:
    :param tolerance:
    :return: True | False
    """

    for an_number in aseq:
        if abs(anum - an_number) < tolerance:
            return True
        else:
            pass

    return False


class EdiCollection(object):
    """
    A super class to encapsulate the properties pertinent to a set of EDI files
    """

    def __init__(self, edilist, ptol=0.05):
        """ constructor
        :param edilist: a list of edifiles with full path, for read-only
        :param ptol: period tolerance considered as equal, default 0.05 means 5 percent
        this param controls what freqs/periods are grouped together:
        10pct may result more double counting of freq/period data than 5pct.
        eg: E:/Data/MT_Datasets/WenPingJiang_EDI 18528 rows vs 14654 rows
        """

        self.edifiles = edilist
        logger.info("number of edi files in this collection: %s",
                    len(self.edifiles))
        assert len(self.edifiles) > 0

        self.num_of_edifiles = len(self.edifiles)  # number of stations
        print ("number of stations/edifiles = %s" % self.num_of_edifiles)

        self.ptol = ptol

        if self.edifiles is not None:
            self.mt_obj_list = [mt.MT(edi) for edi in self.edifiles]
        else:
            logger.error("None Edi file set")

        # get all frequencies from all edi files
        self.all_frequencies = None
        self.mt_periods=None
        self.all_unique_periods = self._get_all_periods()

        self.geopdf = self.create_mt_station_gdf()

        self.bound_box_dict = self.get_bounding_box()  # in orginal projection

        return

    def _get_all_periods(self):
        """
        from the list of edi files get a list of all unique periods from the frequencies.
        """
        if self.all_frequencies is not None:  # already initialized
            return

        # get all frequencies from all edi files
        all_freqs = []
        for mt_obj in self.mt_obj_list:
            all_freqs.extend(list(mt_obj.Z.freq))

        self.mt_periods = 1.0/np.array(all_freqs)

        # sort all frequencies so that they are in ascending order,
        # use set to remove repeats and make an array
        self.all_frequencies = sorted(list(set(all_freqs)))

        logger.debug("Number of MT Frequencies: %s", len(self.all_frequencies))

        all_periods = 1.0 / \
            np.array(sorted(self.all_frequencies, reverse=True))

        #logger.debug("Type of the all_periods %s", type(all_periods))
        logger.info("Number of MT Periods: %s", len(all_periods))
        logger.debug(all_periods)

        return all_periods

    def get_periods_by_stats(self, percentage=10.0):
        """
        check the presence of each period in all edi files, keep a list of periods which are at least percentage present
        :return: a list of periods which are present in at least percentage edi files
        """
        adict = {}
        for aper in self.all_unique_periods:
            station_list = []
            afreq = 1.0 / aper
            acount = 0
            for mt_obj in self.mt_obj_list:
                # if afreq in mt_obj.Z.freq:
                if is_num_in_seq(afreq, mt_obj.Z.freq):
                    acount = acount + 1
                    station_list.append(mt_obj.station)

            if (100.0 * acount) / self.num_of_edifiles >= percentage:
                adict.update({aper: acount})
                #print (aper, acount)
            else:
                logger.info("Period=%s is excluded. it is from stations: %s ", aper, station_list)

        mydict_ordered = sorted(
            adict.items(), key=lambda value: value[1], reverse=True)
        # for apair in mydict_ordered:
        #     print (apair)

        selected_periods = [pc[0] for pc in mydict_ordered]

        print("Selected periods %s out of the total %s:" % (len(selected_periods), len(self.all_unique_periods)))
        return selected_periods

    def create_mt_station_gdf(self, outshpfile=None):
        """
        create station location geopandas dataframe, and output to shape file outshpfile
        :return: gdf
        """

        mt_stations = []

        for mtobj in self.mt_obj_list:
            mt_stations.append(
                [mtobj.station, mtobj.lon, mtobj.lat, mtobj.elev, mtobj.utm_zone])

        pdf = pd.DataFrame(mt_stations, columns=['StationId', 'Lon', 'Lat', 'Elev', 'UtmZone'])

        #print (pdf.head())

        mt_points = [Point(xy) for xy in zip(pdf.Lon, pdf.Lat)]
        # OR pdf['geometry'] = pdf.apply(lambda z: Point(z.Lon, z.Lat), axis=1)
        # if you want to df = df.drop(['Lon', 'Lat'], axis=1)
        # crs0 = {'init': 'epsg:4326'}  # WGS84
        crs0 = {'init': 'epsg:4283'}  # GDA94
        gdf = gpd.GeoDataFrame(pdf, crs=crs0, geometry=mt_points)

        if outshpfile is not None:
            gdf.to_file(outshpfile, driver='ESRI Shapefile')

        return gdf

    def plot_stations(self, savefile=None, showfig=True):
        """
        visualise the geopandas df of MT stations
        :return:
        """

        gdf = self.geopdf
        gdf.plot(figsize=(10, 6), marker='o', color='blue', markersize=5)

        if savefile is not None:
            fig = plt.gcf()
            fig.savefig(savefile, dpi=300)

        if showfig is True:
            plt.show()

        return savefile

    def display_on_basemap(self):
        """

        :return:
        """

        world = gpd.read_file(gpd.datasets.get_path('naturalearth_lowres'))

        logger.debug(world.shape)

        myax = world.plot(alpha=0.5)
        # myax.set_xlim([149,150])
        # myax.set_xlim([110, 155])
        # myax.set_ylim([-40, -10])

        myax.set_xlim(
            (self.bound_box_dict['MinLon'], self.bound_box_dict['MaxLon']))
        myax.set_ylim(
            (self.bound_box_dict['MinLat'], self.bound_box_dict['MaxLat']))

        myax2 = self.geopdf.plot(ax=myax, figsize=(
            10, 6), marker='o', color='blue', markersize=8)

        plt.show()

        return myax2

    def display_on_image(self):
        """
        display/overlay the MT properties on a background geo-referenced map image
        :return:
        """
        import examples.sandpit.plot_geotiff_imshow as plotegoimg

        myax = plotegoimg.plot_geotiff(
            geofile='tests/data/PM_Gravity.tif', show=False)

        margin = 0.02  # degree
        myax.set_xlim(
            (self.bound_box_dict['MinLon'] - margin, self.bound_box_dict['MaxLon'] + margin))
        myax.set_ylim(
            (self.bound_box_dict['MinLat'] - margin, self.bound_box_dict['MaxLat'] + margin))

        myax2 = self.geopdf.plot(ax=myax, figsize=(
            10, 6), marker='o', color='r', markersize=10)

        plt.show()

        return myax2

    def create_measurement_csv(self, dest_dir='/e/tmp'):
        """
        create csv file from the data of EDI files: IMPEDANCE, APPARENT RESISTIVITIES AND PHASES
        see also utils/shapefiles_creator.py
        :return: csvfname
        """
        if dest_dir is None:
            raise Exception("output dir was not provided!!")
        else:
            logger.info("result will be in the dir %s", dest_dir)
            if not os.path.exists(dest_dir):
                os.mkdir(dest_dir)

        # summary csv file
        csv_basename = "edi_measurement"
        csvfname = os.path.join(dest_dir, "%s.csv"%csv_basename)

        pt_dict = {}

        csv_header = [
            'FREQ','STATION', 'LAT', 'LON',  'ZXXre','ZXXim',
            'ZXYre', 'ZXYim', 'ZYXre', 'ZYXim', 'ZYYre', 'ZYYim','TXre', 'TXim', 'TYre', 'TYim',
            'RHOxx', 'RHOxy', 'RHOyx', 'RHOyy', 'PHSxx', 'PHSxy', 'PHSyx', 'PHSyy'
                      ]

        with open(csvfname, "wb") as csvf:
            writer = csv.writer(csvf)
            writer.writerow(csv_header)

        for freq in self.all_frequencies:

            mtlist = []
            for mt_obj in self.mt_obj_list:

                f_index_list = [ff for ff, f2 in enumerate(mt_obj.Z.freq)
                                if (f2 > freq * (1 - self.ptol)) and
                                (f2 < freq * (1 + self.ptol))]
                if len(f_index_list) > 1:
                    logger.warn("more than one fre found %s", f_index_list)

                if len(f_index_list) >= 1:
                    p_index = f_index_list[0]

                    logger.debug("The freqs index %s", f_index_list)
                    # geographic coord lat long and elevation
                    # long, lat, elev = (mt_obj.lon, mt_obj.lat, 0)
                    station, lat, lon = (
                        mt_obj.station, mt_obj.lat, mt_obj.lon)

                    resist_phase =  mtplottools.ResPhase(z_object=mt_obj.Z)
                    # resist_phase.compute_res_phase()

                    mt_stat = [freq, station, lat, lon,
                               mt_obj.Z.z[p_index,0,0].real,
                               mt_obj.Z.z[p_index,0,0].imag,
                               mt_obj.Z.z[p_index,0,1].real,
                               mt_obj.Z.z[p_index,0,1].imag,
                               mt_obj.Z.z[p_index,1,0].real,
                               mt_obj.Z.z[p_index,1,0].imag,
                               mt_obj.Z.z[p_index,1,1].real,
                               mt_obj.Z.z[p_index,1,1].imag,
                               mt_obj.Tipper.tipper[p_index,0,0].real,
                               mt_obj.Tipper.tipper[p_index,0,0].imag,
                               mt_obj.Tipper.tipper[p_index, 0, 1].real,
                               mt_obj.Tipper.tipper[p_index, 0, 1].imag,
                               resist_phase.resxx[p_index], resist_phase.resxy[p_index],
                               resist_phase.resyx[p_index], resist_phase.resyy[p_index],
                               resist_phase.phasexx[p_index], resist_phase.phasexy[p_index],
                               resist_phase.phaseyx[p_index], resist_phase.phaseyy[p_index]
                               ]
                    mtlist.append(mt_stat)

                else:
                    logger.warn(
                        'Freq %s NOT found for this station %s', freq, mt_obj.station)

            with open(csvfname, "ab") as csvf:  # summary csv for all freqs
                writer = csv.writer(csvf)
                writer.writerows(mtlist)

            csv_basename2 = "%s_%sHz.csv" %(csv_basename, str(freq))
            csvfile2 = os.path.join(dest_dir, csv_basename2)

            with open(csvfile2, "wb") as csvf:  # individual csvfile for each freq
                writer = csv.writer(csvf)

                writer.writerow(csv_header)
                writer.writerows(mtlist)

            pt_dict[freq] = mtlist

        return csvfname

    def get_bounding_box(self, epsgcode=None):
        """ compute bounding box
        :return: bounding box in given proj coord system
        """

        if epsgcode is None:
            new_gdf = self.geopdf
        else:  # reproj
            new_gdf = self.geopdf.to_crs(epsg=epsgcode)

        tup = new_gdf.total_bounds

        bdict = {"MinLon": tup[0],
                 "MinLat": tup[1],
                 "MaxLon": tup[2],
                 "MaxLat": tup[3]}

        logger.debug(bdict)

        return bdict

    def show_prop(self):
        """
        show all properties
        :return:
        """
        print (len(self.all_unique_periods), 'unique periods (s)', self.all_unique_periods)
        print (len(self.all_frequencies),
               'unique frequencies (Hz)', self.all_frequencies)

        myper = obj.get_periods_by_stats(percentage=20)
        print(myper)


        print (self.bound_box_dict)

        self.plot_stations(savefile='/e/tmp/edi_collection_test.jpg')

        # self.display_on_basemap()

        self.display_on_image()

        # self.display_folium()

        return

    def create_phase_tensor_csv(self):
        """
        create phase tensor ellipse and tipper properties.
        :return:
        """
        return


if __name__ == "__main__":

# python mtpy/core/edi_collection.py examples/data/edi2/ /e/tmp3/edi2_csv
# python mtpy/core/edi_collection.py /e/Data/MT_Datasets/GA_UA_edited_10s-10000s /e/tmp3/GA_UA_edited_10s-10000s_csv

    if len(sys.argv) < 2:
        print ("\n  USAGE: %s edi_dir OR edi_list " % sys.argv[0])
        sys.exit(1)
    else:
        argv1 = sys.argv[1]
        if os.path.isdir(argv1):
            edis = glob.glob(argv1 + '/*.edi')
            assert len(edis) > 0  # must has edi files
            obj = EdiCollection(edis)

        elif os.path.isfile(argv1) and argv1.endswith('.edi'):
            # assume input is a list of EDI files
            obj = EdiCollection(sys.argv[1:])
        else:
            sys.exit(2)


        #obj.show_prop()

        #print(obj.get_bounding_box(epsgcode=28353))

        #obj.create_mt_station_gdf(outshpfile='/e/tmp/edi_collection_test.shp')

        obj.create_measurement_csv(dest_dir=sys.argv[2])
