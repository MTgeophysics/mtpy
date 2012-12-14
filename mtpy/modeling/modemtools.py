# -*- coding: utf-8 -*-
"""
Created on 12.12.2012

@author: LK
"""

import os
import os.path as op
import sys
import numpy as np
import glob

import mtpy.core.mttools as mtt
import mtpy.modeling.winglinktools as wlt
import mtpy.utils.latlongutmconversion as ll2utm
reload(wlt)

def winglinkmesh2modelfile(WLoutputfile, modelfilename= 'ModEM_initmodel', res_value=100):
    """return init3d file/start model

    mainly copied from ws3dtools...

    Inputs:
        WLoutputfile -  *.out-file from winglink
        modelfilename - name of the modelfile
        res_value - starting homogeneous half space in Ohm-meters (default: 100)


    Output:
        init file path

    """

    #create the output filename
    model_fn = op.abspath(modelfilename)

    #read widths for all blocks from WingLink output file
    dx,dy,dz=wlt.readWLOutFile(WLoutputfile,ncol=5)


    n_we_blocks=len(dx)
    n_ns_blocks=len(dy)
    nz=len(dz)

    init_modelFH = open(model_fn,'w')
    init_modelFH.write('#Initial halfspace model, based on WingLink generated mesh \n')
    init_modelFH.write('%i %i %i 0 \n'%(n_ns_blocks,n_we_blocks,nz))

    #write north block widths
    north_string=''
    count = 0
    for north_idx in range(n_ns_blocks):
        north_string += '%.3e '%(dy[north_idx])
        count +=1
        if count == 8:
            north_string +='\n'
            count = 0

    if n_ns_blocks%8 != 0:
        north_string +='\n'

    init_modelFH.write(north_string)

    #write east block widths
    east_string=''
    count=0
    for east_idx in range(n_we_blocks):
        east_string += '%.3e '%(dx[east_idx])
        count +=1
        if count == 8:
            east_string +='\n'
            count = 0

    if n_we_blocks%8 != 0:
        east_string +='\n'

    init_modelFH.write(east_string)

    #write z block widths/hights
    z_string=''
    count = 0
    for z_idx in range(nz):
        z_string += '%.3e '%(dz[z_idx])
        count +=1
        if count == 8:
            z_string +='\n'
            count = 0
    if nz%8 != 0:
        z_string +='\n'

    init_modelFH.write(z_string)


    #empty line required, if resistivity values are given instead of resistivity indices
    init_modelFH.write('\n')

    for idx_depth in range(nz):
        for idx_n in range(n_ns_blocks):
            we_profile_string =''
            for idx_e in range(n_we_blocks):
                we_profile_string +='%.1f '%(res_value)
            #linebreak after each west-east profile
            init_modelFH.write(we_profile_string+'\n')


    #define origin of model file ... just 0 at the moment
    #assumed to be at the lateral center of the model at the surface
    init_modelFH.write('%.1f %.1f %.1f \n'%(0.,0.,0.))

    #define rotation angle of model w.r.t. data set...just 0 at the moment
    init_modelFH.write('%.1f\n'%(0.))


    init_modelFH.close()

    print 'Wrote initial halfspace model to file: %s '%(model_fn)


    return model_fn



def latlon2xy(lat, lon, origin):
    """returns position in local rectangular coordinates.

    x is positive northwards, y positive eastwards, values given in meter

    origin must be given as a 2-tuple/-array in (lat,lon) form

    """

    #using existing tool for conversion from lat/lon to northing/easting, then just take planar approximation

    #comment: to be checked for overlapping UTM zones!!

    dummy1,utm_local_east,utm_local_north   = ll2utm.LLtoUTM(23, lat, lon)
    dummy2,utm_origin_east,utm_origin_north = ll2utm.LLtoUTM(23, origin[0], origin[1])

    x = utm_local_north - utm_origin_north
    y = utm_local_east  - utm_origin_east

    return x,y


def edis2datafile(winglink_outfile, edilist, sites_file,  comment='Generic datafile, generated from Python script'):

    datafilename = op.abspath('ModEM_inputdata')
    if len(comment)>100:
        sys.exit('comment string is too long (cannot exceed 100 characters)\n')


    lo_frequencies = []
    data_dict = {}

    counter = 0

    #define a little helper
    def closetoexisting(f,l):
        #check, if list l contains an element, which deviates less than 3% from value f
        a_l = np.array(l)
        dev = np.abs(f-a_l).min()/f*100.
        if dev < 3:
            return True
        else:
            return False


    #obtain x,y coordinateds using sites file:
    lo_mesh_xyz_coord_lists = wlt.getmeshblockcoordinates(WL_outfile)
    sites_dict              = wlt.readSitesFile2(sites_file)


    for idx_edi,edifilename in enumerate(edilist):
        #generate overall dictionary containing info from all files

        raw_dict = mtt.readedi(edifilename)
        stationname = raw_dict['station']
        #check, if the station has been used in the model setup:
        if not stationname in sites_dict.keys():
            continue

        data_dict[stationname] = raw_dict
        counter += 1


        raw_freqs = list(raw_dict['frequency'])


        for freq in raw_freqs:
            #check, if freq is already in list, ...
            if len(lo_frequencies) ==0:
                lo_frequencies.append(freq)
            elif freq in lo_frequencies:
                continue
            elif closetoexisting(freq,lo_frequencies):
                continue
            else:
                lo_frequencies.append(freq)


    so_frequencies = list (set(lo_frequencies))
    n_periods      = len(so_frequencies)
    n_stations     = counter

    print 'data from %i stations used for datafile and model setup'%(counter)

    #write header info
    F = open(datafilename,'w')
    F.write('# %s\n'%(comment))
    F.write('# Period Station Lat Lon X Y Z Component Real Imag Error\n')
    F.write('> Full_Impedance\n')
    F.write('> exp(-i\omega t)\n')
    F.write('> [mV/km]/[nT]\n')
    F.write('> 0.00\n')
    F.write('> 0 0 \n')
    F.write('> %i %i\n'%(n_periods, n_stations))

    #define components:
    z_components =['ZXX','ZXY','ZYX','ZYY']


    #iterate over general dictionary and write data file lines successively sorted by stations:
    for station in data_dict:

        station_dict = data_dict[station]
        station_frequencies = station_dict['frequency']
        lat = station_dict['lat']
        lon = station_dict['lon']
        Z   = station_dict['z']
        Z_var = station_dict['zvar']
        stationname = station_dict['station']
        print 'writing data for station %s'%(stationname)

        #no other choice so far...no depth given via EDI file:
        depth  = 0.

        stationdict = sites_dict[station]
        east_idx    = stationdict['idx_east']
        south_idx   = stationdict['idx_south']
        #WingLink indices go from North to South, but here we have South as reference
        north_coordinate = lo_mesh_xyz_coord_lists[0][-south_idx]
        east_coordinate  = lo_mesh_xyz_coord_lists[1][east_idx-1]
        #x,y = latlon2xy(lat, lon, origin)

        for idx_freq, tmp_freq in enumerate(station_frequencies):
            #take frequency from the frequency-list defined above
            correct_frequency = lo_frequencies[ np.abs(tmp_freq-np.array(lo_frequencies)).argmin() ]
            period = 1./correct_frequency

            for idx_comp,comp in enumerate(z_components):
                row = int(idx_comp/2.)
                column = idx_comp%2
                Z_value = Z[idx_freq,row,column]
                err = Z_var[idx_freq,row,column]

                current_data_line = '%f %s %f %f %.1f %.1f %.1f %s %f %f %f \n'%(period, stationname, lat, lon, north_coordinate, east_coordinate,depth,comp, np.real(Z_value), np.imag(Z_value), err)

                F.write(current_data_line)

    F.close()
    print 'wrote datafile %s'%(datafilename)
    return datafilename




def generate_edilist(edifolder):


    lo_edifiles = [op.abspath(i) for i in glob.glob(op.join(edifolder,'*.[eE][dD][iI]'))]

    return lo_edifiles




def winglink2modem(edifolder, winglinkoutput, sites_file, modelfilename='init_model', resistivity=100):


    #check input for consistency
    if not op.isdir(edifolder):
        sys.exit('cannot find EDI files: no such directory: \n%s'%(edifolder))

    edilist = generate_edilist(edifolder)


    if len(edilist)==0:
        sys.exit('cannot find EDI files in given directory: \n%s'%(edifolder))

    if not op.isfile(sites_file):
        sys.exit('cannot open sites information file: \n%s'%(sites_file))

    try:
        WL_outfile = op.abspath(winglinkoutput)
    except:
        sys.exit('cannot find specified WingLink output file: \n%s'%(winglinkoutput))

    try:
        HS_rho_value = float(resistivity)
    except:
        sys.exit('provided resistivity value is not a proper number')


    #set up model file
    modelfn = winglinkmesh2modelfile(WL_outfile, modelfilename=modelfilename, res_value=HS_rho_value)

    #set up data file
    datafn  = edis2datafile(WL_outfile, edilist, sites_file)


    return datafn, modelfn

