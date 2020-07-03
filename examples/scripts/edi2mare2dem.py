# -*- coding: utf-8 -*-
"""
Created on Mon Mar 02 13:22:02 2015

@author: a1655681
TODO: Original author credit

Modified 03-07-2020 17:28:10 AEST 
brenainn.moushall@ga.gov.au
"""
import os
# os.chdir(r'C:\mtpywin\mtpy')

import mtpy.modeling.occam2d as o2d
import mtpy.modeling.mare2dem as m2d

# EDI directory
edi_dir = '/home/bren/mtpy/examples/data/edi_files_2'

# Full path to save Occam2D data file
o2d_path = '/tmp/o2d_data.dat'

# Full path to save Mare2D data file
m2d_path = '/tmp/mare2dem_test.txt'

# Whether to solve statics for all stations
solve_statics = False
# Alternatively, pass a list of station names
# Specified stations will have solve_statics==True, all others False
# solve_statics = ['Synth10', 'Synth11', 'Synth12']

# ASCII grid topo file for interpoalting elevation across the profile
surface_file = '/home/bren/mtpy/examples/data/AussieContinent_etopo1.asc'

# Generate an Occam2D data object from EDI data
gstrike = -72
station_list = m2d.station_list(edi_dir)
o2d_data = o2d.Data(edi_path=edi_dir, model_mode='1', station_list=station_list,
                    interpolate_freq=False, geoelectric_strike=gstrike, res_te_err=20.,
                    phase_te_err=10., res_tm_err=10., phase_tm_err=5.)

# Save the data file
o2d_data.save_path = o2d_path
o2d_data.write_data_file(data_fn=o2d_path)

# Convert the Occam2D profile to Mare2D
mare_origin_x, mare_origin_y, site_locations, site_elevations, m2d_profile = \
    m2d.occam2d_to_mare2dem(o2d_data, surface_file, elevation_sample_n=300)

m2d.write_mare2dem_data(o2d_path, site_locations, site_elevations,
                        (mare_origin_x, mare_origin_y, gstrike),
                        solve_statics=False, savepath=m2d_path)
