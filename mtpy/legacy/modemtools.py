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
import re

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



    for idx_depth in range(nz):
        #empty line required, if resistivity values are given instead of resistivity indices
        init_modelFH.write('\n')

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
    lo_mesh_xyz_coord_lists = wlt.getmeshblockcoordinates(winglink_outfile)
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
    """
    Conversion of WingLink output files into ModEM input.



    """

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




def wsinv2modem_data(wsinv_datafile, sites_file=None):
    """
    Convert an existing input data file from Weerachai's wsinv style into Egbert's ModEM type

    Just provide the wsinv data file and optional the 'sites' file
    The ModEM data file is then saved as 'ModEM_datafile' in the current working directory

    """



    if not op.isfile(wsinv_datafile):
        sys.exit('cannot find input data file:\n%s'%(wsinv_datafile))

    inFH  = open(wsinv_datafile,'r')

    outfn = op.abspath('ModEM_datafile')
    outFH = open(outfn,'w')

    #conversion factor from ohm (ws3dinv) into 'standard' mV/km/nT (ModEM):
    ohm2mvkmnt = 1./796

    #define Z components
    z_components =['ZXX','ZXY','ZYX','ZYY']


    outFH.write('# ModEM datafile, converted from Wsinv3D input\n')
    outFH.write('# Period Station 0 0 X Y 0 Component Real Imag Error\n')
    outFH.write('> Full_Impedance\n')
    outFH.write('> exp(-i\omega t)\n')
    outFH.write('> [mV/km]/[nT]\n')
    outFH.write('> 0.00\n')
    outFH.write('> 0 0 \n')


    indata_raw = inFH.readlines()
    inFH.close()
    indata_list =[]
    for row in indata_raw:
        indata_list.extend(row.strip().split() )

    n_stations = int(indata_list[0])
    n_periods  = int(indata_list[1])

    outFH.write('> %i %i\n'%(n_periods, n_stations))

    if not int(indata_list[2]) == 8:
        outFH.close()
        sys.exit('cannot handle input file - need 8 components, but %i given!'%int(indata_list[2]))

    #obtain ns coordinates:
    lo_coordinates_ns = []
    for n in range(n_stations):
        curr_north = float(indata_list[5+n])
        lo_coordinates_ns.append(curr_north)

    #obtain ew coordinates:
    lo_coordinates_ew = []
    for e in range(n_stations):
        curr_east = float(indata_list[5+n_stations+2+e])
        lo_coordinates_ew.append(curr_east)

    #setup list of station names
    lo_stations = []

    #check, if valid sites file is provided:
    try:
        if not op.isfile(sites_file):
            print 'ERROR -- could not find sites file:\n%s'%(sites_file)
            raise

        SF = open(sites_file,'r')
        SF_data = SF.readlines()
        if not len(SF_data) == n_stations:
            print 'ERROR -- sites file contains wrong number of stations (%i instead of %i)!'%(len(SF_data),n_stations)
            raise

        for i in range(n_stations):
            lo_stations.append(SF_data[i].strip().split()[0][:-4])


    except:
        print 'Could not find proper list of site names. Generic station names are used'
        for i in range(n_stations):
            stationname = 'Sta%02i'%(i+1)
            lo_stations.append(stationname)


    #find starting row for actual data block:
    repattern1 = re.compile(r'DATA_Period')
    r = None
    dummy1 = 0
    for i,line in enumerate(indata_raw):
        r = re.search(repattern1, line)
        if r:
            dummy1 = i
            break

    idx_first_data_row = dummy1


    #find starting row for error block
    repattern2 = re.compile(r'ERROR_Period')
    r = None
    dummy2 = 0
    for i,line in enumerate(indata_raw):
        r = re.search(repattern2, line)
        if r:
            dummy2 = i
            break

    idx_first_error_row = dummy2


    #loop for wrting data to output file; main looping over periods rather than stations for simplicity
    for idx_period in range(n_periods):
        current_starting_line_data = idx_first_data_row + (idx_period * (n_stations+1))
        print current_starting_line_data
        current_period = float(indata_raw[current_starting_line_data].strip().split(':')[1])

        current_starting_line_error = idx_first_error_row + (idx_period * (n_stations+1))

        print idx_period,current_period,current_starting_line_error

        for idx_station in range(n_stations):
            current_station = lo_stations[idx_station]
            north_coord = lo_coordinates_ns[idx_station]
            east_coord = lo_coordinates_ew[idx_station]


            current_data_line  = current_starting_line_data + 1 + idx_station
            current_error_line = current_starting_line_error + 1 + idx_station

            lo_current_data = indata_raw[current_data_line].strip().split()
            lo_current_error= indata_raw[current_error_line].strip().split()


            for idx_comp, comp in enumerate(z_components):

                real_value =  float(lo_current_data[idx_comp*2]) * ohm2mvkmnt
                imag_value = float(lo_current_data[idx_comp*2+1]) * ohm2mvkmnt

                real_error = float(lo_current_error[idx_comp*2]) * ohm2mvkmnt
                imag_error = float(lo_current_error[idx_comp*2+1]) * ohm2mvkmnt

                #only one error for ModEM, so taking the arithmetic mean:
                mean_error = 0.5 * (real_error + imag_error)


                current_data_line = '%.6E %s 0 0 %.1f %.1f 0 %s %.6E %.6E %.6E \n'%(current_period, current_station,  north_coord, east_coord,comp, real_value, imag_value, mean_error)

                outFH.write(current_data_line)

    outFH.close()
    print 'wrote datafile %s'%(outfn)


    return outfn



def wsinv2modem_model(wsinv_modelfile, modeltype='halfspace'):
    """
    Convert an existing input model file from Weerachai's wsinv style into Egbert's ModEM type

    Just provide the wsinv model file
    The ModEM model file is then saved as 'ModEM_modelfile' in the current working directorys

    So far, only halfspace models can be converted (Dec. 2012)

    """

    if not op.isfile(wsinv_modelfile):
        sys.exit('ERROR - could not find input model file:\n%s'%(wsinv_modelfile))

    outfilename = 'ModEM_modelfile'

    Fin = open(wsinv_modelfile,'r')
    modeldata_raw = Fin.readlines()
    Fin.close()

    blockline = modeldata_raw[1].strip().split()
    n_north_blocks = int(blockline[0])
    n_east_blocks  = int(blockline[1])
    n_z_blocks     = int(blockline[2])

    modeltype_index = int(blockline[3])

    if not modeltype_index == 1:
        sys.exit('ERROR - conversion of this model type not supported (yet)!')


    lo_blockwidths = []
    for row in modeldata_raw:
        lo_blockwidths.extend(row.strip().split())

    lo_blockwidths_north = [int(float(i)) for i in lo_blockwidths[6:6+n_north_blocks]]
    lo_blockwidths_east  = [int(float(j)) for j in lo_blockwidths[6+n_north_blocks:6+n_north_blocks+n_east_blocks]]
    lo_blockwidths_z     = [int(float(k)) for k in lo_blockwidths[6+n_north_blocks+n_east_blocks:6+n_north_blocks+n_east_blocks+n_z_blocks]]


    HS_resistivity_value = float(modeldata_raw[-1].strip().split()[0])


    #build new output file
    Fout = open(outfilename,'w')

    Fout.write('#Initial halfspace model, converted from wsinv input model file\n')
    Fout.write('%i %i %i 0 LOGE\n'%(n_north_blocks,n_east_blocks,n_z_blocks))

    #write north block widths
    north_string=''
    count = 0
    for north_idx in range(n_north_blocks):
        north_string += '%i '%(lo_blockwidths_north[north_idx])
        #count +=1
        #if count == 8:
            #north_string +='\n'
            #count = 0

    #if n_north_blocks%8 != 0:
    north_string +='\n'


    Fout.write( north_string)


    #write east block widths
    east_string=''
    count = 0
    for east_idx in range(n_east_blocks):
        east_string += '%i '%(lo_blockwidths_east[east_idx])
        count +=1
        #if count == 8:
            #east_string +='\n'
            #count = 0

    #if n_east_blocks%8 != 0:
    east_string +='\n'


    Fout.write(east_string)


    #write down block heights
    z_string=''
    count = 0
    for z_idx in range(n_z_blocks):
        z_string += '%i '%(lo_blockwidths_z[z_idx])
        #count +=1
        #if count == 8:
            #z_string +='\n'
            #count = 0

    #if n_z_blocks%8 != 0:
    z_string +='\n'

    Fout.write(z_string)

    blockcount = 0
    for idx_depth in range(n_z_blocks):
        #empty line required
        Fout.write('\n')

        for idx_e in range(n_east_blocks):
            ns_profile_string =''
            #data in one line for each north-south line:
            for idx_n in range(n_north_blocks):
                blockcount +=1
                ns_profile_string +='%.5E '%(np.log(HS_resistivity_value))
                #linebreak after each west-east profile
            Fout.write(ns_profile_string+'\n')

    #define origin of model file ... just 0 at the moment
    #assumed to be at the lateral center of the model at the surface
    Fout.write('%.1f %.1f %.1f \n'%(0.,0.,0.))

    #define rotation angle of model w.r.t. data set...just 0 at the moment
    Fout.write('%.1f\n'%(0.))


    Fout.close()
    print 'wrote modelfile %s'%(outfilename)

    return outfilename



def plotmodel3d(modem_modelfile, viewaxis='z',layer=0,savefile=None):
    """
    Plot routine for 3D model.
    Generates a 2D surface section plot of one layer (first layer as default) with the given orientation of viewing axis (default downwards).



    """


    import pylab as p


    pass

    return


def getmeshblockcoordinates(ModEM_modelfile):
    """
    Read a ModEM-style model file and return a list of 3 lists, which again contain the X/Y/Z coordinate of a mesh block

    Orientation is X-North, Y-East, Z-Down.
    Horizontal origin is in the center of the mesh,
    Indexing starts at the lower left (SouthWest) corner

    Referring to a block, which has the position (7 North,12 East,22 Depth), you get the coordinates as

    ( thislist[0][6], thislist[1][11],, thislist[2][21] )


    """

    try:
        ModEMmodelfn = os.path.abspath(os.path.realpath(ModEM_modelfile))

    except:
        sys.exit('ERROR - could not find file:\n%s'%(ModEM_modelfile))

    F = open(ModEMmodelfn, 'r')
    raw_data = F.readlines()
    F.close()

    dims = []
    modeldata_firstline = raw_data[1].strip().split()
    for n in range(3):
        dims.append(int(modeldata_firstline[n]))
    n_north_blocks = dims[0]
    n_east_blocks  = dims[1]
    n_depth_blocks = dims[2]

    north_blockwidths = [int(float(i)) for i in raw_data[2].strip().split()]
    east_blockwidths  = [int(float(j)) for j in raw_data[3].strip().split()]
    depth_blockwidths = [int(float(k)) for k in raw_data[4].strip().split()]

    coord_list_xyz =[]

    total_width_ew = np.sum(east_blockwidths)
    center_ew      = total_width_ew/2.

    total_width_ns = np.sum(north_blockwidths)
    center_ns      = total_width_ns/2.

    total_depth    = np.sum(depth_blockwidths)

    #depths
    lo_depths = []
    current_depth = depth_blockwidths[0]/2.
    lo_depths.append(current_depth)

    for idx_z in range(n_depth_blocks-1):
        current_depth += (depth_blockwidths[idx_z]/2. + depth_blockwidths[idx_z+1]/2.)
        lo_depths.append(current_depth)


    lo_norths = []
    current_north = north_blockwidths[0]/2.
    lo_norths.append(current_north)
    for idx_n in range(n_north_blocks-1):
        current_north += (north_blockwidths[idx_n]/2. + north_blockwidths[idx_n+1]/2.)
        lo_norths.append(current_north)

    lo_norths_centered = list(np.array(lo_norths)-center_ns)
    coord_list_xyz.append(lo_norths_centered)

    lo_easts = []
    current_east= east_blockwidths[0]/2.
    lo_easts.append(current_east)
    for idx_e in range(n_east_blocks-1):
        current_east+= (east_blockwidths[idx_e]/2. + east_blockwidths[idx_e+1]/2.)
        lo_easts.append(current_east)

    lo_easts_centered = list(np.array(lo_easts)-center_ew)
    coord_list_xyz.append(lo_easts_centered)

    coord_list_xyz.append(lo_depths)


    return coord_list_xyz
