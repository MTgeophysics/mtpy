# -*- coding: utf-8 -*-
"""
Description:
    Compute the average resistivity for a given list of EDI files
    
References:
    edi_collection.py
    penetration_depth3d.py

CreationDate:   26/11/2019
Developer:      fei.zhang@ga.gov.au

Revision History:
    LastUpdate:     26/11/2019   FZ
    LastUpdate:     dd/mm/yyyy  Who     Optional description
"""
import os
import sys
import csv
import glob

import numpy as np
import mtpy.core.mt as mt
import mtpy.imaging.mtplottools as mtplottools
from mtpy.utils.decorator import deprecated
from mtpy.utils.matplotlib_utils import gen_hist_bins

from logging import DEBUG, INFO, ERROR, WARNING
from mtpy.utils.mtpylog import MtPyLog

_logger = MtPyLog.get_mtpy_logger(
    __name__
)  # __name__ will be  path.to.module OR __main__
_logger.setLevel(INFO)


def get_resistivity_from_edi_file(edifile, rholist=["det"]):
    """Compute the resistivity values of an edi file
    :param edifile: input edifile
    :param rholist: flag the method to compute penetration depth: det zxy zyx
    :return: a tuple:(station_lat, statoin_lon, periods_list, pendepth_list)
    """
    _logger.debug("processing the edi file %s", edifile)

    mt_obj = mt.MT(edifile)
    zeta = mt_obj.Z  # the attribute Z represent the impedance tensor 2X2 matrix
    freqs = zeta.freq  # frequencies

    scale_param = np.sqrt(1.0 / (2.0 * np.pi * 4 * np.pi * 10 ** (-7)))

    _logger.debug("the scale parameter should be 355.88127 =?= %s", scale_param)

    # The periods array
    periods = 1.0 / freqs

    if "zxy" in rholist:
        # One of the 4-components: XY
        penetration_depth = scale_param * np.sqrt(zeta.resistivity[:, 0, 1] * periods)

    if "zyx" in rholist:
        penetration_depth = scale_param * np.sqrt(zeta.resistivity[:, 1, 0] * periods)

    if "det" in rholist:
        # determinant is |Zeta|**2
        det2 = np.abs(zeta.det[0])
        penetration_depth = scale_param * np.sqrt(0.2 * periods * det2 * periods)

    latlong_d = (mt_obj.lat, mt_obj.lon, periods, penetration_depth)
    return latlong_d


def create_penetration_depth_csv(edi_dir, outputcsv, zcomponent="det"):
    """ Loop over all edi files, and create a csv file with the columns:
    Header Lat, Lon, per0, per1,per2,.....
    lat, lon, pendepth0, pendepth1, ...
    :param edi_dir: path_to_edifiles_dir
    :param zcomponent: det | zxy  | zyx
    :param outputcsv: path2output.csv file
    :return:
    """
    import csv

    if not os.path.isdir(edi_dir):
        _logger.error("input edi directory not exists", edi_dir)
        raise Exception("MTPy Exception: EDI Dir not exist")

    edi_files = glob.glob(os.path.join(edi_dir, "*.edi"))

    _logger.debug(edi_files)

    # the first period list as a reference for checking other stations period
    periods_list0 = None
    latlon_dep = []  # CSV to be returned
    for afile in edi_files:
        # for efile in edi_files[:2]:
        _logger.debug("processing %s", afile)
        lat, lon, periods, depths = get_resistivity_from_edi_file(afile)
        if periods_list0 is None:
            periods_list0 = periods  # initial value assignment
            # depth_string = ','.join(['%.2f' % num for num in depths])
            # latlon_dep.append((lat, lon, depth_string))
            latlon_dep.append(["Lat", "Lon"] + list(periods))  # The first line header
            latlon_dep.append([lat, lon] + list(depths))

        # same length and same values.
        elif len(periods) == len(periods_list0) and (periods == periods_list0).all():
            # depth_string = ','.join(['%.2f' % num for num in depths])
            # latlon_dep.append((lat, lon, depth_string))
            latlon_dep.append([lat, lon] + list(depths))
        else:
            _logger.error("MT Periods Not Equal !! %s VS %s", periods, periods_list0)
            # raise Exception ("MTPy Exception: Periods Not Equal")
            # pass this edi, let's continue

    # logger.debug(latlon_dep)

    if outputcsv is None:
        _logger.error("Output CSV file must be provided", outputcsv)

    _logger.info("Saving to csv file: %s", outputcsv)
    with open(outputcsv, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerows(latlon_dep)

    return latlon_dep


def calculate_aver_resistivity(
    edifiles_list, component="det", rotation_angle=0, out_dir="/c/temp"
):
    """
    calculate the average apparent resistivity of the edifiles_list for each period.
    Algorithm:
    -	1 make sure the stations all have the same period range, if not, interpolate onto common periods
    -	2 rotate to strike angle if necessary
    -	3 calculate resistivity from the determinant of the impedance tensor, or the geometric mean, if necessary
    -	4 get the median resistivity for each period
    -	5 get the median resistivity overall by taking the median of the above

    :param component: =det – default, use determinant of impedance tensor
                      =geom_mean – use geometric mean of the off diagonals sqrt(ZxyXZyx)
                      =separate  a 2x2 array containing average for each component of the impedance tensor.
    :param rotation_angle: angle to rotate the data by before calculating mean.
    :return: A_dictionary=: Period->Median_Resist_On_Stations, OVER_ALL-> Median_Resist
    """

    if not os.path.exists(out_dir):
        os.mkdir(out_dir)

    _logger.info("result will be in the dir %s", out_dir)

    mt_obj_list = [mt.MT(edi) for edi in edifiles_list]

    # summary csv file
    csv_basename = "Average_Resistivity"
    csvfname = os.path.join(out_dir, "%s.csv" % csv_basename)

    pt_dict = {}

    csv_header = [
        "FREQ",
        "STATION",
        "LAT",
        "LON",
        "ZXXre",
        "ZXXim",
        "ZXYre",
        "ZXYim",
        "ZYXre",
        "ZYXim",
        "ZYYre",
        "ZYYim",
        "DETERM",
    ]

    freq_list = all_frequencies

    with open(csvfname, "w", newline="") as csvf:
        writer = csv.writer(csvf)
        writer.writerow(csv_header)

    for freq in freq_list:
        mtlist = []
        for mt_obj in mt_obj_list:
            f_index_list = None
            zobj = None

            # if (interpolate):
            if True:
                f_index_list = [0]

                newZ, newTipper = mt_obj.interpolate([freq], bounds_error=False)

                zobj = newZ
            else:
                freq_max = freq * (1 + ptol)
                freq_min = freq * (1 - ptol)
                f_index_list = np.where(
                    (mt_obj.Z.freq < freq_max) & (mt_obj.Z.freq > freq_min)
                )

                zobj = mt_obj.Z
            # end if

            # print("Debug type(zobj.det) ******", type(zobj.det), zobj.det.size, zobj.det, np.abs(zobj.det[0]))

            if len(f_index_list) > 1:
                _logger.warn("more than one freq found %s", f_index_list)

            if len(f_index_list) >= 1:
                p_index = f_index_list[0]

                _logger.debug("The freqs index %s", f_index_list)
                # geographic coord lat long and elevation
                # long, lat, elev = (mt_obj.lon, mt_obj.lat, 0)
                station, lat, lon = (mt_obj.station, mt_obj.lat, mt_obj.lon)

                mt_stat = [
                    freq,
                    station,
                    lat,
                    lon,
                    zobj.z[p_index, 0, 0].real,
                    zobj.z[p_index, 0, 0].imag,
                    zobj.z[p_index, 0, 1].real,
                    zobj.z[p_index, 0, 1].imag,
                    zobj.z[p_index, 1, 0].real,
                    zobj.z[p_index, 1, 0].imag,
                    zobj.z[p_index, 1, 1].real,
                    zobj.z[p_index, 1, 1].imag,
                    np.abs(zobj.det[0]),
                ]
                mtlist.append(mt_stat)

            else:
                _logger.warn(
                    "Freq %s NOT found for this station %s", freq, mt_obj.station
                )

        with open(csvfname, "a", newline="") as csvf:  # summary csv for all freqs
            writer = csv.writer(csvf)
            writer.writerows(mtlist)

        csv_basename2 = "%s_%sHz.csv" % (csv_basename, str(freq))
        csvfile2 = os.path.join(out_dir, csv_basename2)

        with open(
            csvfile2, "w", newline=""
        ) as csvf:  # individual csvfile for each freq
            writer = csv.writer(csvf)

            writer.writerow(csv_header)
            writer.writerows(mtlist)

        pt_dict[freq] = mtlist

    return pt_dict


# =============================================
# Section for quick test of this script
# ---------------------------------------------
if __name__ == "__main__":
    # call the main function
    calculate_aver_impedance()
