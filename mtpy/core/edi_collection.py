# -*- coding: utf-8 -*-
"""
Description:
To compute and encapsulate the properties of a set of EDI files

Author: fei.zhang@ga.gov.au

CreateDate: 2017-04-20
"""

import csv
import glob
import os
import sys

from logging import INFO

import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from shapely.geometry import Point  # , Polygon, LineString, LinearRing

import mtpy.core.mt as mt
import mtpy.imaging.mtplottools as mtplottools
from mtpy.utils.mtpy_decorator import deprecated
from mtpy.utils.matplotlib_utils import gen_hist_bins
from mtpy.utils.mtpylog import MtPyLog
import mtpy.analysis.pt as MTpt
import mtpy.imaging.penetration

def is_num_in_seq(anum, aseq, atol=0.0001):
    """
    check if anum is in a sequence by a small tolerance

    :param anum: a number to be checked
    :param aseq: a sequence or a list of values
    :param atol: absolute tolerance

    :return: True | False
    """
    # print(np.isclose(anum, aseq, atol=atol))
    # return np.isclose(anum, aseq, atol=atol).any()
    for an_number in aseq:
        if abs(anum - an_number) < atol:
            return True
        else:
            pass

    return False


class EdiCollection(object):
    """
    A super class to encapsulate the properties pertinent to a set of EDI files

    :param edilist: a list of edifiles with full path, for read-only
    :param outdir:  computed result to be stored in outdir
    :param ptol: period tolerance considered as equal, default 0.05 means 5 percent

    The ptol parameter controls what freqs/periods are grouped together:
    10 percent may result more double counting of freq/period data than 5 pct.
    (eg: MT_Datasets/WPJ_EDI)
    """

    def __init__(self, edilist=None, mt_objs=None, outdir=None, ptol=0.05):
        """
        constructor
        """

        #self._logger = MtPyLog.get_mtpy_logger(self.__class__.__name__)  # will be EdiCollection
        self._logger = MtPyLog.get_mtpy_logger(__name__)  # __name__ will be  path.to.module OR __main__
        self._logger.setLevel(INFO)

        if edilist is not None:
            self.edifiles = edilist
            self._logger.info("number of edi files in this collection: %s",
                         len(self.edifiles))
        elif mt_objs is not None:
            self.edifiles = [mt_obj.fn for mt_obj in mt_objs]
        assert len(self.edifiles) > 0

        self.num_of_edifiles = len(self.edifiles)  # number of stations
        print("number of stations/edifiles = %s" % self.num_of_edifiles)

        self.ptol = ptol

        if edilist is not None:
            # if edilist is provided, always create MT objects from the list
            self._logger.debug("constructing MT objects from edi files")
            self.mt_obj_list = [mt.MT(edi) for edi in self.edifiles]
        elif mt_objs is not None:
            # use the supplied mt_objs
            self.mt_obj_list = list(mt_objs)
        else:
            self._logger.error("None Edi file set")

        # get all frequencies from all edi files
        self.all_frequencies = None
        self.mt_periods = None
        self.all_unique_periods = self._get_all_periods()

        self.geopdf = self.create_mt_station_gdf()

        self.bound_box_dict = self.get_bounding_box()  # in orginal projection

        # ensure that outdir is created if not exists.
        if outdir is None:
            #raise Exception("Error: OutputDir is not specified!!!")
            pass
        elif not os.path.exists(outdir):
            os.mkdir(outdir)

        self.outdir = outdir

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

        self.mt_periods = 1.0 / np.array(all_freqs)

        # sort all frequencies so that they are in ascending order,
        # use set to remove repeats and make an array
        self.all_frequencies = sorted(list(set(all_freqs)))

        self._logger.debug("Number of MT Frequencies: %s", len(self.all_frequencies))
        all_periods = 1.0 / np.array(sorted(self.all_frequencies, reverse=True))

        self._logger.debug("Type of all_periods %s", type(all_periods))
        self._logger.info("Number of MT Periods: %s", len(all_periods))
        self._logger.debug("Periods List: %s", str(all_periods))

        return sorted(all_periods)


    def get_period_occurance(self,aper):
        """
        For a given aperiod, compute its occurance frequencies among the stations/edi
        :param aper: a float value of the period
        :return:
        """
        station_list = []
        afreq = 1.0 / aper
        acount = 0
        for mt_obj in self.mt_obj_list:
            # if afreq in mt_obj.Z.freq:
            if is_num_in_seq(afreq, mt_obj.Z.freq):
                acount = acount + 1
                station_list.append(mt_obj.station)

        # print (station_list)

        occ_percentage = (100.0*acount)/self.num_of_edifiles

        return occ_percentage

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
                # print (aper, acount)
            else:
                self._logger.info("Period=%s is excluded. it is from stations: %s ", aper, station_list)

        mydict_ordered = sorted(
            list(adict.items()), key=lambda value: value[1], reverse=True)
        # for apair in mydict_ordered:
        #     print (apair)

        selected_periods = [pc[0] for pc in mydict_ordered]

        print("Selected periods %s out of the total %s:" % (len(selected_periods), len(self.all_unique_periods)))
        return selected_periods

    def select_periods(self, show=True, period_list=None, percentage=10.0):
        """
        Use edi_collection to analyse the whole set of EDI files

        :param show: True or false
        :param period_list:
        :param percentage:
        :return: select_period_list
        """

        uniq_period_list = self.all_unique_periods  # filtered list of periods ?
        print("Unique periods", len(uniq_period_list))

        if show:
            plt.figure()
            plt.clf()
            bins = gen_hist_bins(uniq_period_list)
            plt.hist(self.mt_periods, bins=bins)
            # plt.hist(self.mt_periods, bins=1000)
            plt.title("Histogram with uniq_periods bins")
            plt.xlabel("Periods")
            plt.ylabel("Occurance in number of MT stations")
            plt.show()

        if period_list:
            # 1 ASK user to input a Pmin and Pmax
            # assume uniq_period_list is sorted
            select_period_list = []
            index_start = 0
            for period in period_list:
                for index in range(index_start, len(uniq_period_list)):
                    if (isinstance(period, float) and np.isclose(uniq_period_list[index], period)) or \
                            (isinstance(period, tuple) and period[0] <= uniq_period_list[index] <= period[1]):
                        select_period_list.append(uniq_period_list[index])
                    elif (isinstance(period, float) and uniq_period_list[index] > period) or \
                            (isinstance(period, tuple) and period[1] < uniq_period_list[index]):
                        index_start = index
                        break
            select_period_list = np.array(select_period_list)
        else:
            # 2 percetage stats
            # select commonly occured frequencies from all stations.
            # This could miss some slightly varied frequencies in the middle range.
            select_period_list = np.array(self.get_periods_by_stats(percentage=percentage))

        print("Selected periods ", len(select_period_list))

        return select_period_list
    

    def create_mt_station_gdf(self, outshpfile=None):
        """
        create station location geopandas dataframe, and output to shape file

        :param outshpfile: output file

        :return: gdf
        """

        mt_stations = []

        for mtobj in self.mt_obj_list:
            mt_stations.append(
                [mtobj.station, mtobj.lon, mtobj.lat, mtobj.elev, mtobj.utm_zone])

        pdf = pd.DataFrame(mt_stations, columns=['StationId', 'Lon', 'Lat', 'Elev', 'UtmZone'])

        # print (pdf.head())

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
        Visualise the geopandas df of MT stations

        :param savefile:
        :param showfig:
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
        display MT stations which are in stored in geopandas dataframe in a base map.

        :return: plot object
        """

        world = gpd.read_file(gpd.datasets.get_path('naturalearth_lowres'))

        self._logger.debug(world.shape)

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

        :return: plot object
        """
        import mtpy.utils.plot_geotiff_imshow as plotegoimg

        myax = plotegoimg.plot_geotiff(
            geofile='data/PM_Gravity.tif', show=False)

        margin = 0.02  # degree
        myax.set_xlim(
            (self.bound_box_dict['MinLon'] - margin, self.bound_box_dict['MaxLon'] + margin))
        myax.set_ylim(
            (self.bound_box_dict['MinLat'] - margin, self.bound_box_dict['MaxLat'] + margin))

        myax2 = self.geopdf.plot(ax=myax, figsize=(
            10, 6), marker='o', color='r', markersize=10)

        plt.show()

        return myax2

    def get_phase_tensor_tippers(self, period, interpolate=True):
        """
        For a given MT period (s) value, compute the phase tensor and tippers etc.

        :param period: MT_period (s)
        :param interpolate: Boolean to indicate whether to interpolate on to the given period
        :return: dictionary pt_dict_list

        pt_dict keys ['station', 'freq', 'lon', 'lat', 'phi_min', 'phi_max', 'azimuth', 'skew', 'n_skew', 'elliptic',
                      'tip_mag_re', 'tip_mag_im', 'tip_ang_re', 'tip_ang_im']
        """

        pt_dict_list = []
        plot_per = period
        #plot_per = self.all_unique_periods[1]  # test first

        print("The plot period is ", plot_per)

        for mt_obj in self.mt_obj_list:
            pt_dict = {}
            pt = None
            ti = None

            if(interpolate == False):
                p_index = [ff for ff, f2 in enumerate(1.0/mt_obj.Z.freq)
                           if (f2 > plot_per * (1 - self.ptol)) and
                           (f2 < plot_per * (1 + self.ptol))]

                pt = mt_obj.pt
                ti = mt_obj.Tipper
            else:
                p_index = [0]
                newZ, newTipper = mt_obj.interpolate([1./plot_per], bounds_error=False)

                pt = MTpt.PhaseTensor(z_object=newZ)
                ti = newTipper
            # end if

            if len(p_index) >= 1:
                p_index = p_index[0]
                pt_dict['station']=mt_obj.station
                pt_dict['period'] =plot_per
                pt_dict['lon'] = mt_obj.lon
                pt_dict['lat'] = mt_obj.lat

                pt_dict['phi_min'] = pt.phimin[p_index]
                pt_dict['phi_max'] = pt.phimax[p_index]
                pt_dict['azimuth']= pt.azimuth[p_index]
                pt_dict['skew'] = pt.beta[p_index]
                pt_dict['n_skew'] = 2 * pt.beta[p_index]
                pt_dict['elliptic'] = pt.ellipticity[p_index]

                pt_dict['tip_mag_re']= ti.mag_real[p_index]
                pt_dict['tip_mag_im']= ti.mag_imag[p_index]
                pt_dict['tip_ang_re']= ti.angle_real[p_index]
                pt_dict['tip_ang_im']= ti.angle_imag[p_index]

                pt_dict_list.append(pt_dict)
            else:
                self._logger.warn(" the period %s is NOT found for this station %s. Skipping!!!" % (plot_per, mt_obj.station))

        return pt_dict_list


    def create_phase_tensor_csv(self, dest_dir, period_list=None,
                                interpolate=True,
                                file_name="phase_tensor.csv"):
        """
        create phase tensor ellipse and tipper properties.
        Implementation based on mtpy.utils.shapefiles_creator.ShapeFilesCreator.create_csv_files

        :param dest_dir: output directory
        :param period_list: list of periods; default=None, in which data for all available
                            frequencies are output
        :param interpolate: Boolean to indicate whether to interpolate data onto given period_list
        :param file_name: output file name

        :return: pt_dict
        """
        csvfname = os.path.join(dest_dir, file_name)

        pt_dict = {}

        csv_header = ['station', 'freq', 'lon', 'lat', 'phi_min', 'phi_max', 'azimuth', 'skew',
                      'n_skew', 'elliptic', 'tip_mag_re', 'tip_mag_im', 'tip_ang_re', 'tip_ang_im']

        freq_list = None
        if(period_list is None):
            freq_list = self.all_frequencies
        else:
            freq_list = 1./np.array(period_list)
        # end if

        #with open(csvfname, "wb") as csvf:
        with open(csvfname, "w",newline="") as csvf:
            writer = csv.writer(csvf)
            writer.writerow(csv_header)

            for freq in freq_list:
                ptlist = []
                for mt_obj in self.mt_obj_list:
                    f_index_list = None
                    pt = None
                    ti = None

                    if(interpolate):
                        f_index_list = [0]

                        newZ = None
                        newTipper = None
                        newZ, newTipper = mt_obj.interpolate([freq], bounds_error=False)

                        pt = MTpt.PhaseTensor(z_object=newZ)
                        ti = newTipper
                    else:
                        freq_min = freq * (1 - self.ptol)
                        freq_max = freq * (1 + self.ptol)

                        f_index_list = [ff for ff, f2 in enumerate(mt_obj.Z.freq)
                                        if (f2 > freq_min) and (f2 < freq_max)]
                        pt = mt_obj.pt
                        ti = mt_obj.Tipper
                    #end if

                    if len(f_index_list) > 1:
                        self._logger.warn("more than one freq found %s", f_index_list)
                    if len(f_index_list) >= 1:
                        p_index = f_index_list[0]
                        # geographic coord lat long and elevation
                        # long, lat, elev = (mt_obj.lon, mt_obj.lat, 0)
                        station, lon, lat = (mt_obj.station, mt_obj.lon, mt_obj.lat)

                        pt_stat = [station, freq, lon, lat,
                                   pt.phimin[p_index],
                                   pt.phimax[p_index],
                                   pt.azimuth[p_index],
                                   pt.beta[p_index],
                                   2 * pt.beta[p_index],
                                   pt.ellipticity[p_index],  # FZ: get ellipticity begin here
                                   ti.mag_real[p_index],
                                   ti.mag_imag[p_index],
                                   ti.angle_real[p_index],
                                   ti.angle_imag[p_index]]

                        ptlist.append(pt_stat)
                    else:
                        self._logger.warn("Freq %s NOT found for this station %s", freq, mt_obj.station)

                csv_freq_file = os.path.join(dest_dir,
                                             '{name[0]}_{freq}Hz{name[1]}'.format(
                                                 freq=str(freq), name=os.path.splitext(file_name)))
                #with open(csv_freq_file, "wb") as freq_csvf:
                with open(csv_freq_file, "w",newline="") as freq_csvf:
                    writer_freq = csv.writer(freq_csvf)
                    writer_freq.writerow(csv_header)
                    writer_freq.writerows(ptlist)

                writer.writerows(ptlist)

                pt_dict[freq] = ptlist

        return pt_dict

    # 2020-10: FZ started this new function on request of WPJ
    def create_penetration_depth_csv(self, dest_dir, period_list=None, interpolate=False, file_name="penetration_depth.csv"):
        """
        create penetration depth csv file for each frequency corresponding to the given input 1.0/period_list.
        of course subject to a tolerance.  Note that frequencies values are usually provided in MT EDI files.

        :param dest_dir: output directory
        :param period_list: list of periods; default=None all available periods will be output
        :param interpolate: Boolean to indicate whether to interpolate data onto given period_list
        :param file_name: output files basename

        :return: pt_dict
        """
        csvfname = os.path.join(dest_dir, file_name)

        p_dict = {}

        csv_header = ['station', 'freq', 'lon', 'lat', 'pen_depth_det', 'pen_depth_zxy', 'pen_depth_zyx']

        # convert the period_list into freq array
        freq_list = None
        if(period_list is None):
            freq_list = self.all_frequencies
        else:
            freq_list = 1./np.array(period_list)
        # end if

        with open(csvfname, "w",newline="") as csvf:
            writer = csv.writer(csvf)
            writer.writerow(csv_header)

            for freq in freq_list:
                pdlist = []

                stations, periods, pen_depth_det, latlons = mtpy.imaging.penetration.get_penetration_depth_by_period( self.mt_obj_list, 1.0/freq)  # whichrho='det')
                stations, periods, pen_depth_zxy, latlons = mtpy.imaging.penetration.get_penetration_depth_by_period( self.mt_obj_list, 1.0/freq,whichrho='zxy')
                stations, periods, pen_depth_zyx, latlons = mtpy.imaging.penetration.get_penetration_depth_by_period( self.mt_obj_list, 1.0/freq,whichrho='zyx')


                for iter in range(len(stations)):
                    pdlist.append([stations[iter], freq, latlons[iter][1], latlons[iter][0], pen_depth_det[iter], pen_depth_zxy[iter], pen_depth_zyx[iter]])

                csv_freq_file = os.path.join(dest_dir,
                                             '{name[0]}_{freq}Hz{name[1]}'.format(
                                                 freq=str(freq), name=os.path.splitext(file_name)))

                with open(csv_freq_file, "w",newline="") as freq_csvf:
                    writer_freq = csv.writer(freq_csvf)
                    writer_freq.writerow(csv_header)
                    writer_freq.writerows(pdlist)

                writer.writerows(pdlist)

                p_dict[freq] = pdlist

        return csvfname

    @deprecated("This function is more expensive compared with the method create_phase_tensor_csv(self,)")
    def create_phase_tensor_csv_with_image(self, dest_dir):
        """
        Using PlotPhaseTensorMaps class to generate csv file of phase tensor attributes, etc.
        Only for comparison. This method is more expensive because it will create plot object first.

        :return:
        """
        from mtpy.imaging.phase_tensor_maps import PlotPhaseTensorMaps
        for freq in self.all_frequencies:
            ptm = PlotPhaseTensorMaps(fn_list=self.edifiles, plot_freq=freq, fig_dpi=80, plot_yn='n')
            ptm.export_params_to_file(save_path=dest_dir)
        return

    def create_measurement_csv(self, dest_dir, period_list=None, interpolate=True):
        """
        create csv file from the data of EDI files: IMPEDANCE, APPARENT RESISTIVITIES AND PHASES
        see also utils/shapefiles_creator.py

        :param dest_dir: output directory
        :param period_list: list of periods; default=None, in which data for all available
                            frequencies are output
        :param interpolate: Boolean to indicate whether to interpolate data onto given period_list

        :return: csvfname
        """
        if dest_dir is None:
            dest_dir = self.outdir
        else:
            self._logger.info("result will be in the dir %s", dest_dir)
            if not os.path.exists(dest_dir):
                os.mkdir(dest_dir)

        # summary csv file
        csv_basename = "edi_measurement"
        csvfname = os.path.join(dest_dir, "%s.csv" % csv_basename)

        pt_dict = {}

        csv_header = [
            'FREQ', 'STATION', 'LAT', 'LON', 'ZXXre', 'ZXXim',
            'ZXYre', 'ZXYim', 'ZYXre', 'ZYXim', 'ZYYre', 'ZYYim', 'TXre', 'TXim', 'TYre', 'TYim',
            'RHOxx', 'RHOxy', 'RHOyx', 'RHOyy', 'PHSxx', 'PHSxy', 'PHSyx', 'PHSyy'
        ]

        freq_list = None
        if(period_list is None):
            freq_list = self.all_frequencies
        else:
            freq_list = 1./np.array(period_list)
        # end if

        with open(csvfname, "w",newline="") as csvf:
            writer = csv.writer(csvf)
            writer.writerow(csv_header)


        for freq in freq_list:
            mtlist = []
            for mt_obj in self.mt_obj_list:
                f_index_list = None
                pt = None
                ti = None
                zobj = None
                if (interpolate):
                    f_index_list = [0]

                    newZ = None
                    newTipper = None
                    self._logger.debug("Interpolating frequency ....%s", freq)
                    newZ, newTipper = mt_obj.interpolate([freq], bounds_error=False)

                    pt = MTpt.PhaseTensor(z_object=newZ)
                    ti = newTipper
                    zobj = newZ
                else:
                    freq_max = freq * (1 + self.ptol)
                    freq_min = freq * (1 - self.ptol)
                    f_index_list = np.where((mt_obj.Z.freq < freq_max) & (mt_obj.Z.freq > freq_min))

                    pt = mt_obj.pt
                    ti = mt_obj.Tipper
                    zobj = mt_obj.Z
                # end if

                if len(f_index_list) > 1:
                    self._logger.warn("more than one freq found %s", f_index_list)

                if len(f_index_list) >= 1:
                    p_index = f_index_list[0]

                    self._logger.debug("The freqs index %s", f_index_list)
                    # geographic coord lat long and elevation
                    # long, lat, elev = (mt_obj.lon, mt_obj.lat, 0)
                    station, lat, lon = (
                        mt_obj.station, mt_obj.lat, mt_obj.lon)

                    resist_phase = mtplottools.ResPhase(z_object=zobj)
                    # resist_phase.compute_res_phase()

                    mt_stat = [freq, station, lat, lon,
                               zobj.z[p_index, 0, 0].real,
                               zobj.z[p_index, 0, 0].imag,
                               zobj.z[p_index, 0, 1].real,
                               zobj.z[p_index, 0, 1].imag,
                               zobj.z[p_index, 1, 0].real,
                               zobj.z[p_index, 1, 0].imag,
                               zobj.z[p_index, 1, 1].real,
                               zobj.z[p_index, 1, 1].imag,
                               ti.tipper[p_index, 0, 0].real,
                               ti.tipper[p_index, 0, 0].imag,
                               ti.tipper[p_index, 0, 1].real,
                               ti.tipper[p_index, 0, 1].imag,
                               resist_phase.resxx[p_index], resist_phase.resxy[p_index],
                               resist_phase.resyx[p_index], resist_phase.resyy[p_index],
                               resist_phase.phasexx[p_index], resist_phase.phasexy[p_index],
                               resist_phase.phaseyx[p_index], resist_phase.phaseyy[p_index]
                               ]
                    mtlist.append(mt_stat)

                else:
                    self._logger.warn(
                        'Freq %s NOT found for this station %s', freq, mt_obj.station)

            with open(csvfname, "a",newline="") as csvf:  # summary csv for all freqs
                writer = csv.writer(csvf)
                writer.writerows(mtlist)

            csv_basename2 = "%s_%sHz.csv" % (csv_basename, str(freq))
            csvfile2 = os.path.join(dest_dir, csv_basename2)

            with open(csvfile2, "w", newline="") as csvf:  # individual csvfile for each freq
                writer = csv.writer(csvf)

                writer.writerow(csv_header)
                writer.writerows(mtlist)

            pt_dict[freq] = mtlist

        return csvfname

    def export_edi_files(self, dest_dir, period_list=None,
                                interpolate=True,period_buffer=None,longitude_format='LON'):
        """
        export edi files.
        :param dest_dir: output directory
        :param period_list: list of periods; default=None, in which data for all available
                            frequencies are output
        :param interpolate: Boolean to indicate whether to interpolate data onto given period_list; otherwise
                            a period_list is obtained from get_periods_by_stats()
        :param file_name: output file name
        :param period_buffer: buffer so that interpolation doesn't stretch too far
                              over periods. Provide a float or integer factor, 
                              greater than which interpolation will not stretch.
                              e.g. 1.5 means only interpolate to a maximum of
                              1.5 times each side of each frequency value

        :return:
        """

        if period_list is None:
            period_list = np.array(self.get_periods_by_stats())
        # end if

        for mt_obj in self.mt_obj_list:
            # interpolate each station onto the period list
            # check bounds of period list
            interp_periods = period_list[np.where(
                (period_list >= 1. / mt_obj.Z.freq.max()) &
                (period_list <= 1. / mt_obj.Z.freq.min()))]

            interp_periods = np.sort(interp_periods)
            
            # if specified, apply a buffer so that interpolation doesn't
            # stretch too far over periods
            if type(period_buffer) in [float, int]:
                interp_periods_new = []
                dperiods = 1. / mt_obj.Z.freq
                for iperiod in interp_periods:
                    # find nearest data period
                    difference = np.abs(iperiod - dperiods)
                    nearestdperiod = dperiods[difference == np.amin(difference)][0]
                    if max(nearestdperiod / iperiod, iperiod / nearestdperiod) < period_buffer:
                        interp_periods_new.append(iperiod)

                interp_periods = np.array(interp_periods_new)
                
                
            self._logger.debug("station_name and its original period: %s %s %s",
                               mt_obj.station, len(mt_obj.Z.freq), 1.0 / mt_obj.Z.freq)
            self._logger.debug("station_name and interpolation period: %s %s %s",
                               mt_obj.station, len(interp_periods), interp_periods)

            if len(interp_periods) > 0:  # not empty
                interp_z, interp_t = mt_obj.interpolate(1. / interp_periods)

                if dest_dir is not None and os.path.isdir(dest_dir):
                    mt_obj.write_mt_file(
                        save_dir=dest_dir,
                        fn_basename=mt_obj.station,
                        file_type='edi',
                        new_Z_obj=interp_z,
                        new_Tipper_obj=interp_t,
                        longitude_format=longitude_format)
            else:
                pass
        # end for
    # end func

        return

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

        self._logger.debug(bdict)

        return bdict

    def get_station_utmzones_stats(self):
        """
        A simple method to find what UTM zones these (edi files) MT stations belong to
        are they in a single UTM zone, which corresponds to a unique EPSG code?
        or do they belong to multiple UTM zones?

        :return: a_dict like {UTMZone:Number_of_MT_sites}
        """

        def get_utm_zone(latitude, longitude):
            zone_num = int(1 + (longitude + 180.0) / 6.0)

            if latitude >= 0:
                return "%s%s"%(zone_num,'N')
            else:
                return "%s%s"%(zone_num,'S')

        utm_zones={}
        for mt_obj in self.mt_obj_list:
            utmz = get_utm_zone(mt_obj.lat, mt_obj.lon)
            utm_zones[utmz] = utm_zones.get(utmz,0)+1

        return utm_zones

    def get_stations_distances_stats(self):
        """
        get the min max statistics of the distances between stations.
        useful for determining the ellipses tipper sizes etc

        :return: dict={}
        """
        import math

        mt_stations = []

        for mtobj in self.mt_obj_list:
            mt_stations.append( (mtobj.station, mtobj.lat, mtobj.lon,  mtobj.utm_zone) )

        pdf = pd.DataFrame(mt_stations, columns=['Station', 'Lat', 'Lon',  'UtmZone'])

        mt_distances = []
        for i in range(len(pdf)):
            xi=pdf.iloc[i]['Lat']
            yi=pdf.iloc[i]['Lon']
            for j in range(i+1, len(pdf)):
                xj = pdf.iloc[j]['Lat']
                yj = pdf.iloc[j]['Lon']
                dist = math.sqrt((xi-xj)**2 + (yi - yj)**2)
                mt_distances.append(dist)

                if(dist <0.004): # 0.004 is about 400 meters
                    self._logger.info("Small distances occurred between stations: %s %s", pdf.iloc[i].Station, pdf.iloc[j].Station)

        # print (mt_distances)

        anarray = pd.Series(mt_distances)
        print(anarray.describe())

        q01 = anarray.quantile(q=0.01)
        q02 = anarray.quantile(q=0.02)
        q03 = anarray.quantile(q=0.03)
        q04 = anarray.quantile(q=0.04)
        q05 = anarray.quantile(q=0.05)

        self._logger.info("1,2,3,4 5 Percentile distances: %s, %s, %s, %s, %s", q01, q02,q03,q04,q05)

        # anarray.plot()
        # plt.show()

        min_d = anarray.min()  # cold be very small due to two close stations, skew the result
        max_d = anarray.max()
        self._logger.debug("Minimum = %s", min_d )
        self._logger.debug("Maximum = %s", max_d )

        return {"MIN_DIST":min_d, "Q1PERCENT":q01, "Q2PERCENT":q02, "Q3PERCENT":q03, "MAX_DIST":max_d}

    def show_obj(self, dest_dir=None):
        """
        test call object's methods and show it's properties

        :return:
        """
        print(len(self.all_unique_periods), 'unique periods (s)', self.all_unique_periods)

        print(len(self.all_frequencies),
              'unique frequencies (Hz)', self.all_frequencies)

        myper = self.get_periods_by_stats(percentage=20)

        print(myper)

        print(self.bound_box_dict)

        print(self.get_bounding_box(epsgcode=28353))

        if dest_dir is None:
            self.plot_stations(savefile= os.path.join(self.outdir,'edi_collection_test.jpg'))
        else:
            self.plot_stations(savefile= os.path.join(dest_dir,'edi_collection_test.jpg'))

        # self.display_on_basemap()

        # self.display_on_image()

        # self.display_folium()

        utmzones=self.get_station_utmzones_stats()

        number_zones = len(list(utmzones.items()))
        self._logger.info("This Edi fileset has %s UTM Zone(s): %s ", number_zones, utmzones)

        return

    def get_min_max_distance(self):
        """
        get the min and max distance between all possible pairs of stations.

        :return: min_dist, max_dist
        """
        mt_distances = self.get_stations_distances_stats()
        min_dist = mt_distances.get("MIN_DIST")
        max_dist = mt_distances.get("MAX_DIST")

        return min_dist, max_dist

    def calculate_aver_impedance(self, component="det", rotation_angle=0, out_dir="/c/temp"):
        """
        calculate the average impedance tensor Z (related to apparent resistivity) of all edi (MT-stations) for each period.
        algorithm:
        -	1 make sure the stations all have the same period range, if not, interpolate onto common periods
        -	2 rotate to strike if necessary
        -	3 calculate: the determinant of the impedance tensor, or the geometric mean, if necessary
        -	4 get the median resistivity for each period
        -	5 get the median resistivity overall by taking the median of the above

        :param component: =det – default, returns average for determinant of impedance tensor
                          =geom_mean – returns average geometric mean of the off diagonals sqrt(ZxyXZyx)
                          =separate returns a 2x2 array containing average for each component of the impedance tensor.
        :param rotation_angle: angle to rotate the data by before calculating mean.
        :return: A_dictionary=: Period->Median_Resist_On_Stations, OVER_ALL-> Median_Resist
        """

        if not os.path.exists(out_dir):
            os.mkdir(out_dir)

        self._logger.info("result will be in the dir %s", out_dir)

        # summary csv file
        csv_basename = "z_average_impedance"
        csvfname = os.path.join(out_dir, "%s.csv" % csv_basename)

        pt_dict = {}

        csv_header = [
            'FREQ', 'STATION', 'LAT', 'LON', 'ZXXre', 'ZXXim',
            'ZXYre', 'ZXYim', 'ZYXre', 'ZYXim', 'ZYYre', 'ZYYim',
             "DETERM"
        ]

        freq_list = self.all_frequencies

        with open(csvfname, "w", newline="") as csvf:
            writer = csv.writer(csvf)
            writer.writerow(csv_header)

        for freq in freq_list:
            mtlist = []
            for mt_obj in self.mt_obj_list:
                f_index_list = None
                zobj = None

                #if (interpolate):
                if True:
                    f_index_list = [0]

                    newZ, newTipper = mt_obj.interpolate([freq], bounds_error=False)

                    zobj = newZ
                else:
                    freq_max = freq * (1 + self.ptol)
                    freq_min = freq * (1 - self.ptol)
                    f_index_list = np.where((mt_obj.Z.freq < freq_max) & (mt_obj.Z.freq > freq_min))

                    zobj = mt_obj.Z
                # end if

                # print("Debug type(zobj.det) ******", type(zobj.det), zobj.det.size, zobj.det, np.abs(zobj.det[0]))

                if len(f_index_list) > 1:
                    self._logger.warn("more than one freq found %s", f_index_list)

                if len(f_index_list) >= 1:
                    p_index = f_index_list[0]

                    self._logger.debug("The freqs index %s", f_index_list)
                    # geographic coord lat long and elevation
                    # long, lat, elev = (mt_obj.lon, mt_obj.lat, 0)
                    station, lat, lon = (
                        mt_obj.station, mt_obj.lat, mt_obj.lon)

                    mt_stat = [freq, station, lat, lon,
                               zobj.z[p_index, 0, 0].real,
                               zobj.z[p_index, 0, 0].imag,
                               zobj.z[p_index, 0, 1].real,
                               zobj.z[p_index, 0, 1].imag,
                               zobj.z[p_index, 1, 0].real,
                               zobj.z[p_index, 1, 0].imag,
                               zobj.z[p_index, 1, 1].real,
                               zobj.z[p_index, 1, 1].imag,
                                np.abs(zobj.det[0])
                               ]
                    mtlist.append(mt_stat)

                else:
                    self._logger.warn(
                        'Freq %s NOT found for this station %s', freq, mt_obj.station)

            with open(csvfname, "a", newline="") as csvf:  # summary csv for all freqs
                writer = csv.writer(csvf)
                writer.writerows(mtlist)

            csv_basename2 = "%s_%sHz.csv" % (csv_basename, str(freq))
            csvfile2 = os.path.join(out_dir, csv_basename2)

            with open(csvfile2, "w", newline="") as csvf:  # individual csvfile for each freq
                writer = csv.writer(csvf)

                writer.writerow(csv_header)
                writer.writerows(mtlist)

            pt_dict[freq] = mtlist

        return pt_dict

##################################################################
if __name__ == "__main__":

    # python mtpy/core/edi_collection.py examples/data/edi2/ /c/tmpedi2
    # python mtpy/core/edi_collection.py examples/data/edi_files/ /c/tmp1128
    # python mtpy/core/edi_collection.py examples/data/edi_files_2/ /c/tmp1127

    if len(sys.argv) < 3:
        print("\n  USAGE: %s edi_dir out_Dir" % sys.argv[0])
        sys.exit(1)
    else:
        argv1 = sys.argv[1]
        if os.path.isdir(argv1):
            edis = glob.glob(argv1 + '/*.edi')
            assert len(edis) > 0  # must has edi files
            obj = EdiCollection(edis)

        elif os.path.isfile(argv1) and argv1.endswith('.edi'):
            # assume input is a list of EDI files
            obj = EdiCollection(sys.argv[1:-2])
        else:
            sys.exit(2)

        outdir = sys.argv[2]

        #obj.show_obj(dest_dir = outdir)

        mt_distances = obj.get_stations_distances_stats()
        min_dist = mt_distances.get("MIN_DIST")
        max_dist = mt_distances.get("MAX_DIST")
        print( mt_distances )

        # obj.create_phase_tensor_csv(outdir)
        # obj.create_measurement_csv(dest_dir= outdir, interpolate=True)
        # this function has a bug in the case interpolate=False: check the output csv files you will see a lot of [].
        # obj.create_measurement_csv(dest_dir= outdir, interpolate=False)  # this function has a bug for interpolate=False
        #obj.calculate_aver_impedance(out_dir=outdir)

        # obj.create_mt_station_gdf(os.path.join(outdir, 'edi_collection_test.shp'))
        # obj.create_penetration_depth_csv(dest_dir= outdir, period_list=[0.1067,95.33], interpolate=False)
        obj.create_penetration_depth_csv(dest_dir= outdir, interpolate=False)

