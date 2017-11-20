#!/bin/env python
# -*- coding: utf-8 -*-
"""
Description:
    Compute Phase Tensors from ModEM Dat File and output to CSV file

Usage Examples:
    python mtpy/modeling/modem_data_to_phase_tensor.py examples/data/ModEM_files/Modular_MPI_NLCG_028.dat [OutDir]
    python mtpy/modeling/modem_data_to_phase_tensor.py /e/tmp/GA_UA_edited_10s-10000s_16/ModEM_Data.dat [OutDir]

References:
    https://gajira.atlassian.net/browse/ALAMP-49

CreationDate:   06/09/2017
CreatedBy:      SysUser='u25656'

Developer:      fei.zhang@ga.gov.au
LastUpdate:     08/09/2017
"""

import csv
import os
import sys

import mtpy.analysis.pt as pt
from legacy.modeling import modem_data
from mtpy.core import z


class ModemDataToPhaseTensor(object):

    def __init__(self, mdfile_dat, outdir):

        self.datfile = mdfile_dat
        self.dest_dir= outdir

    def compute_phase_tensor(self):
        """
        Compute the phase tensors from a ModEM dat file
        :param datfile: path2/file.dat
        :return:
        """

        # Data file
        data_file = self.datfile
        dest_dir = self.dest_dir

        # Create a new ModEM data instance
        md = modem_data.Data()
        # Read the datafile
        md.read_data_file(data_fn=data_file)

        num_sites = md.data_array.shape[0]
        print ("ModEM data file number of sites:", num_sites)

        first_site_periods = md.data_array[0][9]  # (23L, 2L, 2L)
        print ("first_site_periods = %s" % str(first_site_periods.shape[0]))

        period_list = md.period_list
        freq_list = 1.0 / period_list
        num_periods = len(period_list)
        print ("ModEM data file number of periods:", num_periods)

        csv_basename ="modem_data_to_phase_tensor"
        csvfname = os.path.join(dest_dir, "%s.csv" % csv_basename)

        csv_header = [
            'Freq', 'Station', 'Lat', 'Long', 'Phimin', 'Phimax', 'Ellipticity', 'Azimuth']

        with open(csvfname, "wb") as csvf:
            writer = csv.writer(csvf)
            writer.writerow(csv_header)

        for period_num in xrange(num_periods):
            per= period_list[period_num]
            freq = freq_list[period_num]
            print ("Working on period",per , "frequency:", freq )
            csvrows = []
            for num_site in range(num_sites):
                # Obtain the site for this number
                this_site = md.data_array[num_site]

                # Station Label is first record in data array for a site
                site_label = this_site[0]
                # Latitude is the second record in data array for a site
                site_lat = this_site[1]
                # Longitude is the third record in data array for a site
                site_long = this_site[2]

                # Actual data is a list of 2x2 matrices (ordered by period), 10th record in data array for a site
                this_site_data = this_site[9]

                # Now get the data for this period only, going off the period number we're looping through currently
                this_period_data = this_site_data[period_num]
                # Create a Z object based on this 2x2 array of data
                this_period_data = z.Z(z_array=this_period_data, freq=freq_list)

                # Given the Z object we just created, give us a PhaseTensor object
                this_phase_tensor = pt.PhaseTensor(z_object=this_period_data)

                # Get the four parameters we care about
                this_phimin = this_phase_tensor.phimin[0][0]
                this_phimax = this_phase_tensor.phimax[0][0]
                this_ellipticity = this_phase_tensor.ellipticity[0][0]
                this_azimuth = this_phase_tensor.azimuth[0][0]

                # Print out comma delimited version of the parameters: label, lat, long, phimin, phimax, ellipticity, azimuth
                arow = [freq, site_label, site_lat, site_long, this_phimin, this_phimax, this_ellipticity, this_azimuth]
                # Done for this site

                csvrows.append(arow)

            with open(csvfname, "ab") as csvf:  # append to this summary csv file for all freqs
                writer = csv.writer(csvf)
                writer.writerows(csvrows)

            csv_basename2 = "%s_%sHz.csv" % (csv_basename, str(freq))
            csvfile2 = os.path.join(dest_dir, csv_basename2)

            with open(csvfile2, "wb") as csvf:  # csvfile  for eachindividual freq
                writer = csv.writer(csvf)
                writer.writerow(csv_header)
                writer.writerows(csvrows)

        # Done with all sites and periods

        return


if __name__ == "__main__":
    """
    How2Run Examples:
    python mtpy/modeling/modem_data_to_phase_tensor.py examples/data/ModEM_files/Modular_MPI_NLCG_028.dat [OutDir]
    python mtpy/modeling/modem_data_to_phase_tensor.py /e/tmp/GA_UA_edited_10s-10000s_16/ModEM_Data.dat [OutDir]
    """
    file_dat = sys.argv[1]
    if len(sys.argv)>2: outdir = sys.argv[2]
    else: outdir='E:/tmp'

    ob = ModemDataToPhaseTensor(file_dat, outdir)
    ob.compute_phase_tensor()
