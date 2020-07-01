# -*- coding: utf-8 -*-
"""
Created on Mon Mar 02 13:22:02 2015

@author: a1655681
TODO: Original author credit
"""
import os
import csv
# os.chdir(r'C:\mtpywin\mtpy')

import scipy
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from mtpy.utils import mesh_tools, gis_tools, filehandling
import mtpy.modeling.occam2d as o2d

# path to work in and where files will be saved
wd = '/tmp/mare2dem'
# specify path to EDIs here if different to the working directory above
epath = '/home/bren/mtpy/examples/data/edi_files_2'
# whether you want the code to solve for static shift (at the moment only coded to have all sites 
# the same, could set independently though)
solve_statics = False
# BM: Topo?
surface_file = '/home/bren/mtpy/examples/data/AussieContinent_etopo1.asc'

# EDI to Occam2D
# Get MT station names from EDI file names
station_list = [os.path.splitext(os.path.basename(fn))[0]
                for fn in os.listdir(epath) if fn.endswith('.edi')]

gstrike = -72
data = o2d.Data(edi_path=epath, model_mode='1', station_list=station_list,
                interpolate_freq=False, geoelectric_strike=gstrike,
                res_te_err=20.,
                phase_te_err=10.,
                res_tm_err=10.,
                phase_tm_err=5.)
data.save_path = wd
data.write_data_file()
data_fn = os.path.join(data.save_path, 'OccamDataFile.dat')

# read in the Occam data file that was just created, to get the station list in the right order
# (different to the order of data.station_list)
new_o2d = o2d.Data()
new_o2d.read_data_file(data.data_fn)
ns, nfreq = len(new_o2d.station_list), len(new_o2d.freq)

# get the eastings and northings for the site locations when they are projected onto the profile
# has to be the edi_list from data not newData because newData doesn't have edi_list filled so can't get projected_east, projected_north
peasts = []
pnorths = []
for i in range(ns):
    peasts.append(data.edi_list[i].projected_east)
    pnorths.append(data.edi_list[i].projected_north)
peasts = np.array(peasts)
pnorths = np.array(pnorths)

# =============================================================================
# Extract information from the projected profile
# =============================================================================
[m, c1] = data.profile_line
# find the start of the profile, (x0, y0)
# assume that the start of the profile always has minimum eastings **don't think this is an accurate assumption, not sure how mtpy deals with profile orientations
x0 = peasts.min()
y0 = m * x0 + c1  # y coordinate at X minimum

# find the end of the profile, (x1,y1), y1 as the northing maximum **again, probably not true
x1 = peasts.max()
y1 = m * x1 + c1  # Y coordinate at X maximum
# this would be the origin for occam (the start of the line). The end of the line is x1, y1
data.profile_origin = [x0, y0]
# origin for Mare2d which is in the middle of the line
mare_origin_x, mare_origin_y = ((peasts.min() + (peasts.max()) / 2.), (pnorths.min() + (pnorths.max()) / 2.))
# calculate the length of the line
profile_length = np.sqrt((y1 - y0)**2 + (x1 - x0)**2)

# =============================================================================
# Elevation profile
# =============================================================================
# choose the points to sample the elevation on for loading into Mare2DEM, ensuring that the exact site locations are added in too,

elev_points_easts = np.linspace(x0, x1, 300, endpoint=False)  # make 300 points between the min and max of the profile for elevation profile, excluding x1 and y1
elev_points_norths = np.linspace(y0, y1, 300, endpoint=False)  # *******need to change 300 to be something sensible depending on profile length
# remove first value of each array (x0 and y0)
elev_points_easts = np.delete(elev_points_easts, [0])
elev_points_norths = np.delete(elev_points_norths, [0])

insert_array = []
for i, elev_e in enumerate(peasts):
    elev_points_easts = np.sort(np.concatenate((elev_points_easts, np.array([elev_e]))))
    # find index i where elev_e is now and put pnorths(elev_e) in that same location
    insert_index = np.where(elev_points_easts == elev_e)  # find the index where the peast has been inserted ******this wouldn't work for a North south line
    if insert_index[0][0] > len(elev_points_norths):
        elev_points_norths = np.append(elev_points_norths, pnorths[i])
    else:
        elev_points_norths = np.insert(elev_points_norths, insert_index[0][0], pnorths[i])  # insert the pnorths(elev_e) into elev_points_norths

# =============================================================================
# CONVERT THE ELEVATION TO PROFILE COORDINATES FOR ALL POINTS
# =============================================================================
loc_full_profile = []
# to convert to profile coordinates, the midpoint of the profile will be 0 in profile coordinates.
for i in range(len(elev_points_easts)):
    elev_profile_loc = np.sqrt((elev_points_norths[(i)] - elev_points_norths[(0)])**2 + (elev_points_easts[i] - elev_points_easts[0])**2)
    if elev_profile_loc < (profile_length / 2):
        loc_full_profile.append(-1 * (profile_length / 2 - elev_profile_loc))
    if elev_profile_loc == (profile_length / 2):
        loc_full_profile.append(0)
    if elev_profile_loc > (profile_length / 2):
        loc_full_profile.append(elev_profile_loc - profile_length / 2)

# read elevation file -
epsg = data.model_epsg
elevation_projected_cubic = mesh_tools.interpolate_elevation_to_grid(
    elev_points_easts, elev_points_norths, epsg=epsg, surfacefile=surface_file,
    method='cubic')
elevation_projected_cubic = -1 * elevation_projected_cubic

# then pull out the elevation from this profile at the locations of each of the sites, to add into the Z column of the data file
# =============================================================================
# GET THE ELEVATION AT EACH SITE AND CONVERT TO PROFILE COORDINATES FOR THE DATA FILE ELEVATIONS
# =============================================================================
# Okay, so this looks like it's enumerating the projected station locations, finding where
# these locations intersect with the 300 points generated for elevation profile and then
# extracting that point from the interpolated elevation
# But it's trying to pull from an X, Y grid using a single coordinate (where easting matches)
elevation_at_sites = []
for p, peast in enumerate(peasts):
    index = np.where(elev_points_easts == peast)
    if elev_points_norths[index] == pnorths[p]:
        elevation_at_sites.append(elevation_projected_cubic[index, index])
elevation_at_sites = np.squeeze(np.array(elevation_at_sites, dtype=np.float64))

# JUST NEED THIS PULLED FROM ELEV_PROJ_CUBIC WITH PNORTHS AND PEASTS, THEN CONVERT THOSE TO MODEL 
# COORDS FOR THE Z COLUMN OF DATA FILE
loc_at_sites = []  # is the elevation at sites but in Mare2D coordinates
for i in range(len(pnorths)):
    PROF_LOC = np.sqrt((pnorths[i] - pnorths[0])**2 + (peasts[i] - peasts[0])**2)
    if PROF_LOC < (profile_length / 2):
        loc_at_sites.append(-1 * (profile_length / 2 - PROF_LOC))
    if PROF_LOC == (profile_length / 2):
        loc_at_sites.append(0)
    if PROF_LOC > (profile_length / 2):
        loc_at_sites.append(PROF_LOC - profile_length / 2)
loc_at_sites = np.squeeze(np.array(loc_at_sites, dtype=np.float64))

# check the elevation along the profile in a plot
plt.figure()
plt.plot(loc_full_profile, elevation_projected_cubic, 'ko')
plt.plot(loc_at_sites, elevation_at_sites, color='r', marker='*')
plt.gca().invert_yaxis()
plt.show()
# =============================================================================
# Get the elev_points_easts and elev_points_norths in Mare2DEM reference system.
# =============================================================================
#elevationmodel = np.stack(([loc_full_profile, elevation_projected_cubic]), axis=1)
#elevationmodel = np.array(elevationmodel, dtype='float64')
#np.savetxt(os.path.join(wd, r'elevation_profile.txt'), elevationmodel)  # write the elevation file that can be imported straight into Mamba2D.m
#

# =============================================================================
# #rewrite the Occam2D data file into Mare2DEM format
# =============================================================================
data_fn = data_fn.strip()
with open(data_fn, 'r') as f:
    read_data = f.readlines()

data_list = []
new_data = []
flag = 0
mare_fn = os.path.join(wd, 'Mare2D_data.txt')
# reordering the columns of data from Occam to MARE2DEM style and replacing the data keys from Occam to Mare2D style:

with open(mare_fn, 'w') as output:
    for line in read_data:
        line.strip()
        line.split()
        if 'SITE ' in line:
            flag = 1
        else:
            if flag == 0:
                pass
            else:
                data_list.append(line)
    for i in range(0, ns + 1):  # nfreq):
        for j in data_list:
            if str(i) in j.split('  ')[1] and len(str(i)) == len(j.split('  ')[1]):
                reordered_j = '\n' + '\t' + j.split()[2] + '\t\t' + j.split()[1] + '\t\t' + j.split()[0] + '\t\t' + j.split()[0] + '\t\t' + j.split()[3] + '\t\t' + j.split()[4]
                new_data.append(reordered_j)
    new_data_replaced_keys = [datakeys.replace('\n\t2\t', '\n\t104\t') for datakeys in new_data]
    new_data_replaced_keys = [datakeys.replace('\n\t6\t', '\n\t106\t') for datakeys in new_data_replaced_keys]
    new_data_replaced_keys = [datakeys.replace('\n\t9\t', '\n\t103\t') for datakeys in new_data_replaced_keys]
    new_data_replaced_keys = [datakeys.replace('\n\t10\t', '\n\t105\t') for datakeys in new_data_replaced_keys]
    new_data_replaced_keys = [datakeys.replace('\n\t5\t', '\n\t125\t') for datakeys in new_data_replaced_keys]
    new_data_replaced_keys = [datakeys.replace('\n\t1\t', '\n\t123\t') for datakeys in new_data_replaced_keys]
    new_data_replaced_keys = [datakeys.replace('\n\t3\t', '\n\t133\t') for datakeys in new_data_replaced_keys]
    new_data_replaced_keys = [datakeys.replace('\n\t4\t', '\n\t134\t') for datakeys in new_data_replaced_keys]

    # 1. header
    fstring = 'Format:  EMData_2.2\n'
    fstring += 'UTM of x,y origin (UTM zone, N, E, 2D strike):'
    fstring += ' 54S{:>13.1f}{:>13.1f}\t{:d}\n'.format(mare_origin_x, mare_origin_y, gstrike)  # *********AT THE MOMENT UTM ZONE IS HARDCODED...

    # 2. frequencies
    fstring += '# MT Frequencies:    {}\n'.format(nfreq)
    fstring += '\n'.join([str(round(f, 8)) for f in data.freq])

    # 3. receiver info
    # Prepare the data for the MT reciever section
    # Zeros of shape (n_sites) for X (as float), Theta, Alpha, Beta and Length (ints) columns
    x_col = np.zeros(loc_at_sites.shape, dtype=np.float64)
    zero_ints = np.zeros(loc_at_sites.shape, dtype=np.int8)
    t_col, a_col, b_col, l_col = zero_ints, zero_ints, zero_ints, zero_ints
    # add 0.1 m (shift the sites 10 cm beneath subsurface as recommended)
    elevation_at_sites += 0.1
    site_names = [sn.split('_')[0] for sn in new_o2d.station_list]
    statics = np.ones(loc_at_sites.shape, dtype=np.int8) if solve_statics else zero_ints
    # Put into dataframe for easier stringifying
    header = ['X', 'Y', 'Z', 'Theta', 'Alpha', 'Beta', 'Length', 'SolveStatic', 'Name']
    recv_df = pd.DataFrame((x_col, loc_at_sites, elevation_at_sites, t_col, a_col, b_col, l_col,
                            statics, site_names)).T
    recv_df.columns = header
    recv_str = list(recv_df.to_string(index=False, float_format=lambda x: '%.6f' % x))
    recv_str[0] = '!'

    fstring += '\n# MT Receivers:      {}\n'.format(ns)
    fstring += "".join(recv_str)
    fstring += '\n'

    # 4. data
    fstring += '# Data:       {}\n'.format(len(new_data))
    fstring += '!  Type  Freq #    Tx #    Rx #           Data         StdErr'
    formats = ['%7i', '%7i', '%7i', '%7i', '%15.5f', '%15.5f']

    output.write(fstring)
    output.writelines(new_data_replaced_keys)
