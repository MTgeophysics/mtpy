# -*- coding: utf-8 -*-
"""
Created on Mon Mar 02 13:22:02 2015

@author: a1655681
"""
import os

os.chdir(r"C:\mtpywin\mtpy")
import mtpy.modeling.occam2d as o2d
import os
import os.path as op
import numpy as np
import csv
import pandas as pd
import dask.dataframe as dd
import mtpy.utils.gis_tools as gis
import scipy
import matplotlib.pyplot as plt

# path to work in and where files will be saved
wd = r"C:\Users\roberk20\Dropbox\Lithospheric_Architecture_Shared_Folder\Software\MARE2DEM_COMPILED\Inversions\Musgraves\mare2d"
epath = wd  # specify path to EDIs here if different to the working directory above
solve_statics = False  # whether you want the code to solve for static shift (at the moment only coded to have all sites the same, could set independently though)
csv_fn = r"C:\Users\roberk20\Dropbox\Lithospheric_Architecture_Shared_Folder\Software\MARE2DEM_COMPILED\DataFile\etopo1_bedrock_smaller_UTM__127_131_-32_-23.csv"
# =============================================================================


elst = [
    op.join(epath, ef) for ef in os.listdir(epath) if ef.endswith(".edi")
]  # make a list of EDIs
slst = [
    edi[0:-4] for edi in os.listdir(epath) if edi.find(".edi") > 0
]  # make a list of MT stations
gstrike = -72
data = o2d.Data(
    edi_path=epath,
    model_mode="1",
    station_list=slst,
    interpolate_freq=False,
    geoelectric_strike=gstrike,
    res_te_err=20.0,
    phase_te_err=10.0,
    res_tm_err=10.0,
    phase_tm_err=5.0,
)
data._fill_data()
data.save_path = wd
data.write_data_file()
data_fn = op.join(data.save_path, "OccamDataFile.dat")
# data_slst=data.station_list
# data_sloc=data.station_locations
# read in the Occam data file that was just created, to get the station list in the right order (different to the order of data.station_list)
newData = o2d.Data()
newData.read_data_file(data.data_fn)
newData_slst = newData.station_list
newData_sloc = newData.station_locations
ns, nfreq = len(newData.station_list), len(newData.freq)

## read elevation file -
df = dd.read_csv(csv_fn, header=None, sep="\t")

df = df.compute()
## extract numpy array from pandas dataframe
eastings, northings, elev = df.values[:, 0], df.values[:, 1], df.values[:, 2]

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
mare_origin_x, mare_origin_y = (
    (peasts.min() + (peasts.max()) / 2.0),
    (pnorths.min() + (pnorths.max()) / 2.0),
)
# calculate the length of the line
profile_length = np.sqrt((y1 - y0) ** 2 + (x1 - x0) ** 2)

# =============================================================================
# Elevation profile
# =============================================================================
# choose the points to sample the elevation on for loading into Mare2DEM, ensuring that the exact site locations are added in too,

elev_points_easts = np.linspace(
    x0, x1, 300, endpoint=False
)  # make 300 points between the min and max of the profile for elevation profile, excluding x1 and y1
elev_points_norths = np.linspace(
    y0, y1, 300, endpoint=False
)  # *******need to change 300 to be something sensible depending on profile length
# remove first value of each array (x0 and y0)
elev_points_easts = np.delete(elev_points_easts, [0])
elev_points_norths = np.delete(elev_points_norths, [0])

insert_array = []
for i, elev_e in enumerate(peasts):
    elev_points_easts = np.sort(np.concatenate((elev_points_easts, np.array([elev_e]))))
    # find index i where elev_e is now and put pnorths(elev_e) in that same location
    insert_index = np.where(
        elev_points_easts == elev_e
    )  # find the index where the peast has been inserted ******this wouldn't work for a North south line
    if insert_index[0][0] > len(elev_points_norths):
        elev_points_norths = np.append(elev_points_norths, pnorths[i])
    else:
        elev_points_norths = np.insert(
            elev_points_norths, insert_index[0][0], pnorths[i]
        )  # insert the pnorths(elev_e) into elev_points_norths

# =============================================================================
# CONVERT THE ELEVATION TO PROFILE COORDINATES FOR ALL POINTS
# =============================================================================
loc_full_profile = []
# to convert to profile coordinates, the midpoint of the profile will be 0 in profile coordinates.
for i in range(len(elev_points_easts)):
    elev_profile_loc = np.sqrt(
        (elev_points_norths[(i)] - elev_points_norths[(0)]) ** 2
        + (elev_points_easts[i] - elev_points_easts[0]) ** 2
    )
    if elev_profile_loc < (profile_length / 2):
        loc_full_profile.append(-1 * (profile_length / 2 - elev_profile_loc))
    if elev_profile_loc == (profile_length / 2):
        loc_full_profile.append(0)
    if elev_profile_loc > (profile_length / 2):
        loc_full_profile.append(elev_profile_loc - profile_length / 2)

elevation_projected_cubic = scipy.interpolate.griddata(
    (eastings, northings), elev, (elev_points_easts, elev_points_norths), method="cubic"
)
elevation_projected_cubic = -1 * elevation_projected_cubic
# then pull out the elevation from this profile at the locations of each of the sites, to add into the Z column of the data file
# =============================================================================
# GET THE ELEVATION AT EACH SITE AND CONVERT TO PROFILE COORDINATES FOR THE DATA FILE ELEVATIONS
# =============================================================================
elevation_at_sites = []
for p, peast in enumerate(peasts):
    index = np.where(elev_points_easts == peast)
    if elev_points_norths[index[0][0]] == pnorths[p]:
        elevation_at_sites.append(elevation_projected_cubic[index[0][0]])
    # JUST NEED THIS PULLED FROM ELEV_PROJ_CUBIC WITH PNORTHS AND PEASTS, THEN CONVERT THOSE TO MODEL COORDS FOR THE Z COLUMN OF DATA FILE
loc_at_sites = []  # is the elevation at sites but in Mare2D coordinates
for i in range(len(pnorths)):
    PROF_LOC = np.sqrt((pnorths[i] - pnorths[0]) ** 2 + (peasts[i] - peasts[0]) ** 2)
    if PROF_LOC < (profile_length / 2):
        loc_at_sites.append(-1 * (profile_length / 2 - PROF_LOC))
    if PROF_LOC == (profile_length / 2):
        loc_at_sites.append(0)
    if PROF_LOC > (profile_length / 2):
        loc_at_sites.append(PROF_LOC - profile_length / 2)

# check the elevation along the profile in a plot
plt.figure()
plt.plot(loc_full_profile, elevation_projected_cubic, "ko")
plt.plot(loc_at_sites, elevation_at_sites, color="r", marker="*")
plt.gca().invert_yaxis()
plt.show()
# =============================================================================
# Get the elev_points_easts and elev_points_norths in Mare2DEM reference system.
# =============================================================================

elevationmodel = np.stack(([loc_full_profile, elevation_projected_cubic]), axis=1)
elevationmodel = np.array(elevationmodel, dtype="float64")
np.savetxt(
    os.path.join(wd, r"elevation_profile.txt"), elevationmodel
)  # write the elevation file that can be imported straight into Mamba2D.m


# =============================================================================
# #rewrite the Occam2D data file into Mare2DEM format
# =============================================================================
data_fn = data_fn.strip()
with open(data_fn, "r") as f:
    read_data = f.readlines()

data_list = []
new_data = []
flag = 0
mare_fn = os.path.join(wd, "Mare2D_data.txt")
### reordering the columns of data from Occam to MARE2DEM style and replacing the data keys from Occam to Mare2D style:

with open(mare_fn, "w") as output:

    for line in read_data:
        line.strip()
        line.split()
        if "SITE " in line:
            flag = 1
        else:
            if flag == 0:
                pass
            else:
                data_list.append(line)
    for i in range(0, ns + 1):  # nfreq):
        for j in data_list:
            if str(i) in j.split("  ")[1] and len(str(i)) == len(j.split("  ")[1]):
                reordered_j = (
                    "\n"
                    + "\t"
                    + j.split()[2]
                    + "\t\t"
                    + j.split()[1]
                    + "\t\t"
                    + j.split()[0]
                    + "\t\t"
                    + j.split()[0]
                    + "\t\t"
                    + j.split()[3]
                    + "\t\t"
                    + j.split()[4]
                )
                new_data.append(reordered_j)
    new_data_replaced_keys = [
        datakeys.replace("\n\t2\t", "\n\t104\t") for datakeys in new_data
    ]
    new_data_replaced_keys = [
        datakeys.replace("\n\t6\t", "\n\t106\t") for datakeys in new_data_replaced_keys
    ]
    new_data_replaced_keys = [
        datakeys.replace("\n\t9\t", "\n\t103\t") for datakeys in new_data_replaced_keys
    ]
    new_data_replaced_keys = [
        datakeys.replace("\n\t10\t", "\n\t105\t") for datakeys in new_data_replaced_keys
    ]
    new_data_replaced_keys = [
        datakeys.replace("\n\t5\t", "\n\t125\t") for datakeys in new_data_replaced_keys
    ]
    new_data_replaced_keys = [
        datakeys.replace("\n\t1\t", "\n\t123\t") for datakeys in new_data_replaced_keys
    ]
    new_data_replaced_keys = [
        datakeys.replace("\n\t3\t", "\n\t133\t") for datakeys in new_data_replaced_keys
    ]
    new_data_replaced_keys = [
        datakeys.replace("\n\t4\t", "\n\t134\t") for datakeys in new_data_replaced_keys
    ]

    # 1. header
    fstring = "Format:  EMData_2.2\n"
    fstring += "UTM of x,y origin (UTM zone, N, E, 2D strike):"
    fstring += " 54S{:>13.1f}{:>13.1f}\t{:d}\n".format(
        mare_origin_x, mare_origin_y, gstrike
    )  ##*********AT THE MOMENT UTM ZONE IS HARDCODED...

    # 2. frequencies
    fstring += "# MT Frequencies:    {}\n".format(nfreq)
    fstring += "\n".join([str(round(f, 8)) for f in data.freq])

    # 3. receiver info
    fstring += "\n# MT Receivers:      {}\n".format(ns)
    fstring += "!            X            Y            Z   Theta   Alpha    Beta   Length SolveStatic  Name\n"
    zl = [0] * ns
    slst = [sn.split("_")[0] for sn in newData_slst]
    # add 0.1 m (shift the sites 10 cm beneath subsurface as recommended)
    elevation_at_sites = np.ndarray.tolist(np.array(elevation_at_sites) + 0.1)
    receiver_lst = np.vstack(
        [zl, loc_at_sites]
        + [elevation_at_sites]
        + [zl] * 4
        + [int(solve_statics)]
        + [slst]
    ).T
    formats = [
        "{:>14.2f}",
        "{:>13.2f}",
        "{:>13.2f}",
        "{:>8.2f}",
        "{:>8.2f}",
        "{:>8.2f}",
        "{:>8.2f}",
        "{:>12.0f}",
        "{:>6}\n",
    ]
    for i in range(ns):
        for j in range(9):
            if j < 8:
                value = float(receiver_lst[i, j])
            else:
                value = receiver_lst[i, j]
            fstring += formats[j].format(value)

    # 4. data
    fstring += "# Data:       {}\n".format(len(new_data))  #
    fstring += "!  Type  Freq #    Tx #    Rx #           Data         StdErr"
    formats = ["%7i", "%7i", "%7i", "%7i", "%15.5f", "%15.5f"]

    output.write(fstring)
    output.writelines(new_data_replaced_keys)
