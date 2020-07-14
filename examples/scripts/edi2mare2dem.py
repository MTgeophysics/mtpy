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

mtpy_dir = '/path/to/mtpy'

# EDI directory
edi_dir = os.path.join(mtpy_dir,'examples','data','edi_files_2')


savepath = '/c/tmp'

# Full path to save Occam2D data file
o2d_path = os.path.join(savepath,'Occamdata.dat')

# Full path to save Mare2D data file
m2d_path = os.path.join(savepath,'MARE2Ddata.dat')

# Whether to solve statics for all stations
solve_statics = False
# Alternatively, pass a list of station names
# Specified stations will have solve_statics==True, all others False
# solve_statics = ['Synth10', 'Synth11', 'Synth12']

# ASCII grid topo file for interpoalting elevation across the profile
surface_file = os.path.join(mtpy_dir,'examples','data', 
                            'AussieContinent_etopo1.asc')

# Generate an Occam2D data object from EDI data
o2d_data = o2d.Data(edi_path=edi_dir, model_mode='1', optimize_line=True,
                    interpolate_freq=False, res_te_err=20.,
                    phase_te_err=10., res_tm_err=10., phase_tm_err=5.)

# Save the data file
o2d_data.save_path = o2d_path
o2d_data.write_data_file(data_fn=o2d_path)

# Convert the Occam2D profile to Mare2D
mare_origin_x, mare_origin_y, utm_zone, site_locations, site_elevations, m2d_profile, profile_elevation = \
    m2d.occam2d_to_mare2dem(o2d_data, surface_file, elevation_sample_n=300)

# Plot the profile and site locations against elevation
fig = m2d.plot(m2d_profile, profile_elevation, site_locations, site_elevations)
# fig.show()
fig.savefig(os.path.join(savepath,'m2d_plot.png'), dpi=400)

# Save the profile elevation to file that can be opened in Mamba2D
m2d.write_elevation_file(m2d_profile, profile_elevation,
                         os.path.join(savepath,'elevation.txt'))

gstrike = o2d_data.geoelectric_strike

m2d.write_mare2dem_data(o2d_path, site_locations, site_elevations,
                        (mare_origin_x, mare_origin_y, utm_zone),
                        gstrike=gstrike, solve_statics=False, savepath=m2d_path)
