"""
Description:
    Compute Phase Tensors from ModEM Dat File

LastUpdated:    2017-09-06
Author:         fei.zhang@ga.gov.au
"""
import os,sys

from mtpy.modeling import modem_data

# Import Z to deal with impedence data
from mtpy.core import z

# Phase tensor analysis toolkit
import mtpy.analysis.pt as pt

def compute_phase_tensor(datfile):
    """
    Compute the phase tensors from a ModEM dat file
    :param datfile: path2/file.dat
    :return:
    """

    # Data file
    data_file = datfile

    # Create a new ModEM data instance
    md = modem_data.Data()
    # Read the datafile
    md.read_data_file(data_fn=data_file)

    num_sites = md.data_array.shape[0]
    print "ModEM data file number of sites:", num_sites

    period_list = md.period_list
    freq_list = 1.0 / period_list
    num_periods = md.data_array[0][9].shape[0]
    num_periods = len(period_list)
    print "ModEM data file number of periods:", num_periods

    for period_num in range(num_periods):
        print "Working on period", period_list[period_num], "frequency:", freq_list[period_num]
        for num_site in range(num_sites):
            # Obtain the site for this number
            this_site = md.data_array[num_site]

            # Label is first record in data array for a site
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
            print str(site_label) + "," + str(site_lat) + "," + str(site_long) + "," + str(this_phimin) + "," + str(
                this_phimax) + "," + str(this_ellipticity) + "," + str(this_azimuth)
        # Done for this site
        # Done for this period
        # Print a blank line for neatness
        print ""
    # Done with all sites and periods

    return


if __name__ == "__main__":
    """
    How2Run:
    python examples/phase_tensor_from_model.py examples/data/ModEM_files/Modular_MPI_NLCG_028.dat
    """
    file_dat = sys.argv[1]
    compute_phase_tensor(file_dat)
