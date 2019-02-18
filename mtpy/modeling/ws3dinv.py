# -*- coding: utf-8 -*-
"""
===============
ws3dinv
===============

    * Deals with input and output files for ws3dinv written by:
        Siripunvaraporn, W.; Egbert, G.; Lenbury, Y. & Uyeshima, M. 
        Three-dimensional magnetotelluric inversion: data-space method
        Physics of The Earth and Planetary Interiors, 2005, 150, 3-14
	* Dependencies: matplotlib 1.3.x, numpy 1.7.x, scipy 0.13               
                    and evtk if vtk files want to be written.
					
The intended use or workflow is something like this for getting started:

:Making input files: ::

    >>> import mtpy.modeling.ws3dinv as ws
    >>> import os
    >>> #1) make a list of all .edi files that will be inverted for 
    >>> edi_path = r"/home/EDI_Files"
    >>> edi_list = [os.path.join(edi_path, edi) for edi in edi_path 
    >>> ...         if edi.find('.edi') > 0]
    >>> #2) make a grid from the stations themselves with 200m cell spacing
    >>> wsmesh = ws.WSMesh(edi_list=edi_list, cell_size_east=200, 
    >>> ...                cell_size_north=200)
    >>> wsmesh.make_mesh()
    >>> # check to see if the mesh is what you think it should be
    >>> wsmesh.plot_mesh()
    >>> # all is good write the mesh file
    >>> wsmesh.write_initial_file(save_path=r"/home/ws3dinv/Inv1")
    >>> # note this will write a file with relative station locations
    >>> #change the starting model to be different than a halfspace
    >>> mm = ws.WS3DModelManipulator(initial_fn=wsmesh.initial_fn)
    >>> # an interactive gui will pop up to change the resistivity model
    >>> #once finished write a new initial file
    >>> mm.rewrite_initial_file()
    >>> #3) write data file
    >>> wsdata = ws.WSData(edi_list=edi_list, station_fn=wsmesh.station_fn)
    >>> wsdata.write_data_file()
    >>> #4) plot mt response to make sure everything looks ok
    >>> rp = ws.PlotResponse(data_fn=wsdata.data_fn)
    >>> #5) make startup file
    >>> sws = ws.WSStartup(data_fn=wsdata.data_fn, initial_fn=mm.new_initial_fn)
    
:checking the model and response: ::
    
    >>> mfn = r"/home/ws3dinv/Inv1/test_model.01"
    >>> dfn = r"/home/ws3dinv/Inv1/WSDataFile.dat"
    >>> rfn = r"/home/ws3dinv/Inv1/test_resp.01" 
    >>> sfn = r"/home/ws3dinv/Inv1/WS_Sation_Locations.txt"
    >>> # plot the data vs. model response
    >>> rp = ws.PlotResponse(data_fn=dfn, resp_fn=rfn, station_fn=sfn)
    >>> # plot model slices where you can interactively step through
    >>> ds = ws.PlotSlices(model_fn=mfn, station_fn=sfn)
    >>> # plot phase tensor ellipses on top of depth slices
    >>> ptm = ws.PlotPTMaps(data_fn=dfn, resp_fn=rfn, model_fn=mfn)
    >>> #write files for 3D visualization in Paraview or Mayavi
    >>> ws.write_vtk_files(mfn, sfn, r"/home/ParaviewFiles")
    
    
    
Created on Sun Aug 25 18:41:15 2013

@author: jpeacock-pr
"""

#==============================================================================

import os
import numpy as np

import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
from matplotlib.patches import Ellipse
from matplotlib.colors import Normalize
import matplotlib.colorbar as mcb
import matplotlib.gridspec as gridspec
import mtpy.core.z as mtz
import mtpy.core.mt as mt
import mtpy.imaging.mtplottools as mtplottools
import matplotlib.widgets as widgets
import matplotlib.colors as colors
import matplotlib.cm as cm
import mtpy.modeling.winglink as wl
import mtpy.utils.exceptions as mtex
import mtpy.analysis.pt as mtpt
import mtpy.imaging.mtcolors as mtcl

try:
    from evtk.hl import gridToVTK, pointsToVTK
except ImportError:
    print ('If you want to write a vtk file for 3d viewing, you need download '
           'and install evtk from https://bitbucket.org/pauloh/pyevtk')

    print ('Note: if you are using Windows you should build evtk first with'
           'either MinGW or cygwin using the command: \n'
           '    python setup.py build -compiler=mingw32  or \n'
           '    python setup.py build -compiler=cygwin')

#==============================================================================

#==============================================================================
# Data class
#==============================================================================
class WSData(object):
    """
    Includes tools for reading and writing data files intended to be used with
    ws3dinv.

    :Example: ::

        >>> import mtpy.modeling.ws3dinv as ws
        >>> import os
        >>> edi_path = r"/home/EDI_Files"
        >>> edi_list = [os.path.join(edi_path, edi) for edi in edi_path 
        >>> ...         if edi.find('.edi') > 0]
        >>> # create an evenly space period list in log space
        >>> p_list = np.logspace(np.log10(.001), np.log10(1000), 12)
        >>> wsdata = ws.WSData(edi_list=edi_list, period_list=p_list, 
        >>> ...                station_fn=r"/home/stations.txt")
        >>> wsdata.write_data_file()


    ====================== ====================================================
    Attributes              Description
    ====================== ====================================================
    data                   numpy structured array with keys:
                               * *station* --> station name
                               * *east*    --> relative eastern location in
                                               grid
                               * *north*   --> relative northern location in 
                                               grid
                               * *z_data*  --> impedance tensor array with 
                                               shape
                                         (n_stations, n_freq, 4, dtype=complex)
                               * *z_data_err--> impedance tensor error without
                                                error map applied
                               * *z_err_map --> error map from data file
    data_fn                full path to data file
    edi_list               list of edi files used to make data file
    n_z                    [ 4 | 8 ] number of impedance tensor elements
                           *default* is 8
    ncol                   number of columns in out file from winglink
                           *default* is 5
    period_list            list of periods to invert for
    ptol                   if periods in edi files don't match period_list
                           then program looks for periods within ptol
                           *defualt* is .15 or 15 percent
    rotation_angle         Angle to rotate the data relative to north.  Here 
                           the angle is measure clockwise from North, 
                           Assuming North is 0 and East is 90.  Rotating data,
                           and grid to align with regional geoelectric strike
                           can improve the inversion. *default* is None
    save_path              path to save the data file
    station_fn             full path to station file written by WSStation
    station_locations      numpy structured array for station locations keys:
                               * *station* --> station name
                               * *east*    --> relative eastern location in
                                               grid
                               * *north*   --> relative northern location in
                                               grid
                           if input a station file is written
    station_east           relative locations of station in east direction
    station_north          relative locations of station in north direction
    station_names          names of stations

    units                  [ 'mv' | 'else' ] units of Z, needs to be mv for 
                           ws3dinv. *default* is 'mv'
    wl_out_fn              Winglink .out file which describes a 3D grid
    wl_site_fn             Wingling .sites file which gives station locations
    z_data                 impedance tensors of data with shape:
                           (n_station, n_periods, 2, 2)
    z_data_err             error of data impedance tensors with error map 
                           applied, shape (n_stations, n_periods, 2, 2)
    z_err                  [ float | 'data' ] 
                           'data' to set errors as data errors or
                           give a percent error to impedance tensor elements
                           *default* is .05 or 5%  if given as percent, ie. 5%
                           then it is converted to .05.  
    z_err_floor            percent error floor, anything below this error will
                           be set to z_err_floor.  *default* is None
    z_err_map              [zxx, zxy, zyx, zyy] for n_z = 8
                           [zxy, zyx] for n_z = 4
                           Value in percent to multiply the error by, which 
                           give the user power to down weight bad data, so 
                           the resulting error will be z_err_map*z_err
    ====================== ====================================================

    ====================== ====================================================
    Methods                Description
    ====================== ====================================================
    build_data             builds the data from .edi files
    write_data_file        writes a data file from attribute data.  This way 
                           you can read in a data file, change some parameters
                           and rewrite.
    read_data_file         reads in a ws3dinv data file
    ====================== ====================================================

    """

    def __init__(self, **kwargs):

        self.save_path = kwargs.pop('save_path', None)
        self.data_basename = kwargs.pop('data_basename', 'WSDataFile.dat')
        self.units = kwargs.pop('units', 'mv')
        self.ncol = kwargs.pop('ncols', 5)
        self.ptol = kwargs.pop('ptol', 0.15)
        self.z_err = kwargs.pop('z_err', 0.05)
        self.z_err_floor = kwargs.pop('z_err_floor', None)
        self.z_err_map = kwargs.pop('z_err_map', [10,1,1,10])
        self.n_z = kwargs.pop('n_z', 8)
        self.period_list = kwargs.pop('period_list', None)
        self.edi_list = kwargs.pop('edi_list', None)
        self.station_locations = kwargs.pop('station_locations', None)
        self.rotation_angle = kwargs.pop('roatation_angle', None)

        self.station_east = None
        self.station_north = None
        self.station_names = None

        self.z_data = None
        self.z_data_err = None

        self.wl_site_fn = kwargs.pop('wl_site_fn', None)
        self.wl_out_fn = kwargs.pop('wl_out_fn', None)

        self.data_fn = kwargs.pop('data_fn', None)
        self.station_fn = kwargs.pop('station_fn', None)

        self.data = None

        # make sure the error given is a decimal percent
        if type(self.z_err) is not str and self.z_err > 1:
            self.z_err /= 100.

        # make sure the error floor given is a decimal percent
        if self.z_err_floor is not None and self.z_err_floor > 1:
            self.z_err_floor /= 100.

    def build_data(self):
        """
        Builds the data from .edi files to be written into a data file

        Need to call this if any parameters have been reset to write a 
        correct data file.

        """

        if self.edi_list is None:
            raise WSInputError('Need to input a list of .edi files to build '
                               'the data')

        if self.period_list is None:
            raise WSInputError('Need to input a list of periods to extract '
                               'from the .edi files.' )

        #get units correctly
        if self.units == 'mv':
            zconv = 1./796.

        self.period_list = np.array(self.period_list)

        #define some lengths
        n_stations = len(self.edi_list)
        n_periods = len(self.period_list)

        #make a structured array to keep things in for convenience
        z_shape = (n_periods, 2, 2)
        data_dtype = [('station', '|S10'),
                      ('east', np.float),
                      ('north', np.float),
                      ('z_data', (np.complex, z_shape)),
                      ('z_data_err', (np.complex, z_shape)),
                      ('z_err_map', (np.complex, z_shape))]
        self.data = np.zeros(n_stations, dtype=data_dtype)

        #------get station locations-------------------------------------------
        if self.wl_site_fn != None:
            if self.wl_out_fn is None:
                raise IOError('Need to input an .out file to get station'
                              'locations, this should be output by Winglink')

            #get x and y locations on a relative grid
            east_list, north_list, station_list = \
                            wl.get_station_locations(self.wl_site_fn,
                                                     self.wl_out_fn,
                                                     ncol=self.ncol)
            self.data['station'] = station_list
            self.data['east'] = east_list
            self.data['north'] = north_list

        #if a station location file is input
        if self.station_fn != None:
            stations = WSStation(self.station_fn)
            stations.read_station_file()
            self.data['station'] = stations.names
            self.data['east'] = stations.east
            self.data['north'] = stations.north

        #if the user made a grid in python or some other fashion
        if self.station_locations != None:
            try:
                for dd, sd in enumerate(self.station_locations):
                    self.data['east'][dd] = sd['east_c']
                    self.data['north'][dd] = sd['north_c']
                    self.data['station'][dd] = sd['station']

                    stations = WSStation()
                    stations.station_fn = os.path.join(self.save_path,
                                                    'WS_Station_locations.txt')
                    stations.east = self.data['east']
                    stations.north = self.data['north']
                    stations.names = self.data['station']
                    stations.write_station_file()

            except (KeyError, ValueError):
                self.data['east'] = self.station_locations[:, 0]
                self.data['north']= self.station_locations[:, 1]

        #--------find frequencies----------------------------------------------
        for ss, edi in enumerate(self.edi_list):
            if not os.path.isfile(edi):
                raise IOError('Could not find '+edi)

            mt_obj = mt.MT(edi)
            if self.rotation_angle is not None:
                mt_obj.rotation_angle = self.rotation_angle
            print('{0}{1}{0}'.format('-'*20, mt_obj.station))

            # get only those periods that are within the station data
            interp_periods = self.period_list[np.where(
                                (self.period_list >= 1./mt_obj.Z.freq.max()) &
                                (self.period_list <= 1./mt_obj.Z.freq.min()))]

            #interpolate over those periods
            interp_z, interp_t = mt_obj.interpolate(1./interp_periods)

            for kk, ff in enumerate(interp_periods):
                jj = np.where(self.period_list == ff)[0][0]
                print('    {0:.6g} (s)'.format(ff))

                self.data[ss]['z_data'][jj, :] = interp_z.z[kk, :, :]*zconv
                self.data[ss]['z_data_err'][jj, :] = interp_z.z_err[kk, :, :]*zconv


    def compute_errors(self):
        """
        compute the errors from the given attributes
        """

        for d_arr in self.data:
            if self.z_err == 'data':
                pass
            elif self.z_err_floor is None and type(self.z_err) is float:
                d_arr['z_data_err'][:] = d_arr['z_data'][:]*self.z_err
            elif self.z_err_floor is not None:
                ef_idx = np.where(d_arr['z_data_err'] < self.z_err_floor)
                d_arr['z_data_err'][ef_idx] = d_arr['z_data'][ef_idx]*self.z_err_floor


            d_arr['z_err_map'] = np.reshape(len(self.period_list)*self.z_err_map,
                                           (len(self.period_list), 2, 2))

    def write_data_file(self, **kwargs):
        """
        Writes a data file based on the attribute data

        Key Word Arguments:
        ---------------------
            **data_fn** : string
                          full path to data file name

            **save_path** : string
                            directory path to save data file, will be written
                            as save_path/data_basename
            **data_basename** : string
                                basename of data file to be saved as 
                                save_path/data_basename
                                *default* is WSDataFile.dat

        .. note:: if any of the data attributes have been reset, be sure
                  to call build_data() before write_data_file.  
        """
        if self.data is None:
            self.build_data()

        # compute errors, this helps when rewriting a data file
        self.compute_errors()
            
        for key in ['data_fn', 'save_path', 'data_basename']:
            try:
                setattr(self, key, kwargs[key])
            except KeyError:
                pass

        #create the output filename
        if self.save_path == None:
            if self.wl_out_fn is not None:
                self.save_path = os.path.dirname(self.wl_site_fn)
            else:
                self.save_path = os.getcwd()
            self.data_fn = os.path.join(self.save_path, self.data_basename)
        elif os.path.isdir(self.save_path) == True:
            self.data_fn = os.path.join(self.save_path, self.data_basename)
        else:
            self.data_fn = self.save_path

        #-----Write data file--------------------------------------------------
        n_stations = len(self.data)
        n_periods = self.data[0]['z_data'].shape[0]

        ofid = file(self.data_fn, 'w')
        ofid.write('{0:d} {1:d} {2:d}\n'.format(n_stations, n_periods,
                                                self.n_z))

        #write N-S locations
        ofid.write('Station_Location: N-S \n')
        for ii in range(n_stations/self.n_z+1):
            for ll in range(self.n_z):
                index = ii*self.n_z+ll
                try:
                    ofid.write('{0:+.4e} '.format(self.data['north'][index]))
                except IndexError:
                    pass
            ofid.write('\n')

        #write E-W locations
        ofid.write('Station_Location: E-W \n')
        for ii in range(n_stations/self.n_z+1):
            for ll in range(self.n_z):
                index = ii*self.n_z+ll
                try:
                    ofid.write('{0:+.4e} '.format(self.data['east'][index]))
                except IndexError:
                    pass
            ofid.write('\n')

        #write impedance tensor components
        for ii, p1 in enumerate(self.period_list):
            ofid.write('DATA_Period: {0:3.6f}\n'.format(p1))
            for ss in range(n_stations):
                zline = self.data[ss]['z_data'][ii].reshape(4,)
                for jj in range(self.n_z/2):
                    ofid.write('{0:+.4e} '.format(zline[jj].real))
                    ofid.write('{0:+.4e} '.format(-zline[jj].imag))
                ofid.write('\n')

        #write error as a percentage of Z
        for ii, p1 in enumerate(self.period_list):
            ofid.write('ERROR_Period: {0:3.6f}\n'.format(p1))
            for ss in range(n_stations):
                zline = self.data[ss]['z_data_err'][ii].reshape(4,)
                for jj in range(self.n_z/2):
                    ofid.write('{0:+.4e} '.format(zline[jj].real))
                    ofid.write('{0:+.4e} '.format(zline[jj].imag))
                ofid.write('\n')

        #write error maps
        for ii, p1 in enumerate(self.period_list):
            ofid.write('ERMAP_Period: {0:3.6f}\n'.format(p1))
            for ss in range(n_stations):
                zline = self.data[ss]['z_err_map'][ii].reshape(4,)
                for jj in range(self.n_z/2):
                    ofid.write('{0:.5e} '.format(self.z_err_map[jj]))
                    ofid.write('{0:.5e} '.format(self.z_err_map[jj]))
                ofid.write('\n')
        ofid.close()
        print('Wrote file to: {0}'.format(self.data_fn))

        self.station_east = self.data['east']
        self.station_north = self.data['north']
        self.station_names = self.data['station']
        self.z_data = self.data['z_data']
        self.z_data_err = self.data['z_data_err']*self.data['z_err_map']


    def read_data_file(self, data_fn=None, wl_sites_fn=None, station_fn=None):
        """
        read in data file

        Arguments:
        -----------
            **data_fn** : string
                          full path to data file
            **wl_sites_fn** : string
                              full path to sites file output by winglink.
                              This is to match the station name with station
                              number.
            **station_fn** : string
                             full path to station location file written by 
                             WSStation

        Fills Attributes:
        ------------------
            **data** : structure np.ndarray
                      fills the attribute WSData.data with values

            **period_list** : np.ndarray()
                             fills the period list with values.
        """

        if self.units == 'mv':
            zconv = 796.
        else:
            zconv = 1

        if data_fn is not None:
            self.data_fn = data_fn

        if self.data_fn is None:
            raise WSInputError('Need to input a data file')

        if os.path.isfile(self.data_fn) is False:
            raise WSInputError('Could not find {0}, check path'.format(
                                self.data_fn))

        self.save_path = os.path.dirname(self.data_fn)

        dfid = file(self.data_fn, 'r')
        dlines = dfid.readlines()

        #get size number of stations, number of frequencies,
        # number of Z components
        n_stations, n_periods, nz = np.array(dlines[0].strip().split(),
                                             dtype='int')
        nsstart = 2

        self.n_z = nz
        #make a structured array to keep things in for convenience
        z_shape = (n_periods, 2, 2)
        data_dtype = [('station', '|S10'),
                      ('east', np.float),
                      ('north', np.float),
                      ('z_data', (np.complex, z_shape)),
                      ('z_data_err', (np.complex, z_shape)),
                      ('z_err_map', (np.complex, z_shape))]
        self.data = np.zeros(n_stations, dtype=data_dtype)

        findlist = []
        for ii, dline in enumerate(dlines[1:50], 1):
            if dline.find('Station_Location: N-S') == 0:
                findlist.append(ii)
            elif dline.find('Station_Location: E-W') == 0:
                findlist.append(ii)
            elif dline.find('DATA_Period:') == 0:
                findlist.append(ii)

        ncol = len(dlines[nsstart].strip().split())

        #get site names if entered a sites file
        if wl_sites_fn != None:
            self.wl_site_fn = wl_sites_fn
            slist, station_list = wl.read_sites_file(self.wl_sites_fn)
            self.data['station'] = station_list

        elif station_fn != None:
            self.station_fn = station_fn
            stations = WSStation(self.station_fn)
            stations.read_station_file()
            self.data['station'] = stations.names
        else:
            self.data['station'] = np.arange(n_stations)


        #get N-S locations
        for ii, dline in enumerate(dlines[findlist[0]+1:findlist[1]],0):
            dline = dline.strip().split()
            for jj in range(ncol):
                try:
                    self.data['north'][ii*ncol+jj] = float(dline[jj])
                except IndexError:
                    pass
                except ValueError:
                    break

        #get E-W locations
        for ii, dline in enumerate(dlines[findlist[1]+1:findlist[2]],0):
            dline = dline.strip().split()
            for jj in range(self.n_z):
                try:
                    self.data['east'][ii*ncol+jj] = float(dline[jj])
                except IndexError:
                    pass
                except ValueError:
                    break
        #make some empty array to put stuff into
        self.period_list = np.zeros(n_periods)

        #get data
        per = 0
        error_find = False
        errmap_find = False
        for ii, dl in enumerate(dlines[findlist[2]:]):
            if dl.lower().find('period') > 0:
                st = 0

                if dl.lower().find('data') == 0:
                    dkey = 'z_data'
                    self.period_list[per] = float(dl.strip().split()[1])

                elif dl.lower().find('error') == 0:
                    dkey = 'z_data_err'
                    if not error_find:
                        error_find = True
                        per = 0

                elif dl.lower().find('ermap') == 0:
                    dkey = 'z_err_map'
                    if not errmap_find:
                        errmap_find = True
                        per = 0

                #print '-'*20+dkey+'-'*20
                per += 1

            else:
                if dkey == 'z_err_map':
                    zline = np.array(dl.strip().split(), dtype=np.float)
                    self.data[st][dkey][per-1,:] = np.array([[zline[0]-1j*zline[1],
                                                        zline[2]-1j*zline[3]],
                                                        [zline[4]-1j*zline[5],
                                                        zline[6]-1j*zline[7]]])
                else:
                    zline = np.array(dl.strip().split(), dtype=np.float)*zconv
                    self.data[st][dkey][per-1,:] = np.array([[zline[0]-1j*zline[1],
                                                        zline[2]-1j*zline[3]],
                                                        [zline[4]-1j*zline[5],
                                                        zline[6]-1j*zline[7]]])
                st += 1


        self.station_east = self.data['east']
        self.station_north = self.data['north']
        self.station_names = self.data['station']
        self.z_data = self.data['z_data']
        #need to be careful when multiplying complex numbers
        self.z_data_err = \
                self.data['z_data_err'].real*self.data['z_err_map'].real+1j*\
                self.data['z_data_err'].imag*self.data['z_err_map'].imag

        #make station_locations structure array
        self.station_locations = np.zeros(len(self.station_east),
                                          dtype=[('station','|S10'),
                                                 ('east', np.float),
                                                 ('north', np.float),
                                                 ('east_c', np.float),
                                                 ('north_c', np.float)])
        self.station_locations['east'] = self.data['east']
        self.station_locations['north'] = self.data['north']
        self.station_locations['station'] = self.data['station']

#==============================================================================
# stations
#==============================================================================
class WSStation(object):
    """
    read and write a station file where the locations are relative to the 
    3D mesh.

    ==================== ======================================================
    Attributes           Description
    ==================== ======================================================
    east                 array of relative locations in east direction    
    elev                 array of elevations for each station
    names                array of station names
    north                array of relative locations in north direction  
    station_fn           full path to station file
    save_path            path to save file to
    ==================== ======================================================

    ==================== ======================================================
    Methods              Description 
    ==================== ======================================================
    read_station_file    reads in a station file
    write_station_file   writes a station file  
    write_vtk_file       writes a vtk points file for station locations
    ==================== ======================================================
    """

    def __init__(self, station_fn=None, **kwargs):
        self.station_fn = station_fn
        self.east = kwargs.pop('east', None)
        self.north = kwargs.pop('north', None)
        self.elev = kwargs.pop('elev', None)
        self.names = kwargs.pop('names', None)
        self.save_path = kwargs.pop('save_path', None)

    def write_station_file(self, east=None, north=None, station_list=None,
                           save_path=None, elev=None):
        """
        write a station file to go with the data file.

        the locations are on a relative grid where (0, 0, 0) is the 
        center of the grid.  Also, the stations are assumed to be in the center
        of the cell.

        Arguments:
        -----------
            **east** : np.ndarray(n_stations)
                       relative station locations in east direction

            **north** : np.ndarray(n_stations)
                       relative station locations in north direction

            **elev** : np.ndarray(n_stations)
                       relative station locations in vertical direction

            **station_list** : list or np.ndarray(n_stations)
                               name of stations

            **save_path** : string
                            directory or full path to save station file to
                            if a directory  the file will be saved as
                            save_path/WS_Station_Locations.txt
                            if save_path is none the current working directory
                            is used as save_path

        Outputs:
        ---------
            **station_fn** : full path to station file

        """
        if east is not None:
            self.east = east
        if north is not None:
            self.north = north
        if station_list is not None:
            self.names = station_list
        if elev is not None:
            self.elev = elev
        else:
            if self.north is not None:
                self.elev = np.zeros_like(self.north)

        if save_path is not None:
            self.save_path = save_path
            if os.path.isdir(self.save_path):
                self.station_fn = os.path.join(self.save_path,
                                               'WS_Station_Locations.txt')
            else:
                self.station_fn = save_path
        elif self.save_path is None:
            self.save_path = os.getcwd()
            self.station_fn = os.path.join(self.save_path,
                                           'WS_Station_Locations.txt')
        elif os.path.isdir(self.save_path):
            self.station_fn = os.path.join(self.save_path,
                                           'WS_Station_Locations.txt')

        sfid = file(self.station_fn, 'w')
        sfid.write('{0:<14}{1:^14}{2:^14}{3:^14}\n'.format('station', 'east',
                                                    'north', 'elev'))
        for ee, nn, zz, ss in zip(self.east, self.north, self.elev, self.names):
            ee = '{0:+.4e}'.format(ee)
            nn = '{0:+.4e}'.format(nn)
            zz = '{0:+.4e}'.format(zz)
            sfid.write('{0:<14}{1:^14}{2:^14}{3:^14}\n'.format(ss, ee, nn, zz))
        sfid.close()

        print('Wrote station locations to {0}'.format(self.station_fn))

    def read_station_file(self, station_fn=None):
        """
        read in station file written by write_station_file

        Arguments:
        ----------
            **station_fn** : string
                             full path to station file

        Outputs:
        ---------
            **east** : np.ndarray(n_stations)
                       relative station locations in east direction

            **north** : np.ndarray(n_stations)
                       relative station locations in north direction
            **elev** : np.ndarray(n_stations)
                       relative station locations in vertical direction

            **station_list** : list or np.ndarray(n_stations)
                               name of stations

        """
        if station_fn is not None:
            self.station_fn = station_fn

        self.save_path = os.path.dirname(self.station_fn)

        self.station_locations = np.loadtxt(self.station_fn, skiprows=1,
                                            dtype=[('station', '|S10'),
                                                   ('east_c', np.float),
                                                  ('north_c', np.float),
                                                  ('elev', np.float)])

        self.east = self.station_locations['east_c']
        self.north = self.station_locations['north_c']
        self.names = self.station_locations['station']
        self.elev = self.station_locations['elev']


    def write_vtk_file(self, save_path, vtk_basename='VTKStations'):
        """
        write a vtk file to plot stations

        Arguments:
        ------------
            **save_path** : string
                            directory to save file to.  Will save as 
                            save_path/vtk_basename

            **vtk_basename** : string
                               base file name for vtk file, extension is 
                               automatically added.
        """
        if os.path.isdir(save_path) == True:
            save_fn = os.path.join(save_path, vtk_basename)

        if self.elev is None:
            self.elev = np.zeros_like(self.north)

        pointsToVTK(save_fn, self.north, self.east, self.elev,
                    data={'value':np.ones_like(self.north)})

        return save_fn

    def from_wl_write_station_file(self, sites_file, out_file, ncol=5):
        """
        write a ws station file from the outputs of winglink

        Arguments:
        -----------
            **sites_fn** : string
                           full path to sites file output from winglink

            **out_fn** : string
                         full path to .out file output from winglink

            **ncol** : int
                       number of columns the data is in
                       *default* is 5


        """

        wl_east, wl_north, wl_station_list = wl.get_station_locations(
                                                                    sites_file,
                                                                    out_file,
                                                                    ncol=ncol)
        self.write_station_file(east=wl_east, north=wl_north,
                                station_list=wl_station_list)

#==============================================================================
# mesh class
#==============================================================================
class WSMesh(object):
    """
    make and read a FE mesh grid

    The mesh assumes the coordinate system where:
        x == North
        y == East
        z == + down

    All dimensions are in meters.

    :Example: ::

        >>> import mtpy.modeling.ws3dinv as ws
        >>> import os
        >>> #1) make a list of all .edi files that will be inverted for 
        >>> edi_path = r"/home/EDI_Files"
        >>> edi_list = [os.path.join(edi_path, edi) for edi in edi_path 
        >>> ...         if edi.find('.edi') > 0]
        >>> #2) make a grid from the stations themselves with 200m cell spacing
        >>> wsmesh = ws.WSMesh(edi_list=edi_list, cell_size_east=200, 
        >>> ...                cell_size_north=200)
        >>> wsmesh.make_mesh()
        >>> # check to see if the mesh is what you think it should be
        >>> wsmesh.plot_mesh()
        >>> # all is good write the mesh file
        >>> wsmesh.write_initial_file(save_path=r"/home/ws3dinv/Inv1")

    ==================== ======================================================
    Attributes           Description    
    ==================== ======================================================
    cell_size_east       mesh block width in east direction
                         *default* is 500
    cell_size_north      mesh block width in north direction
                         *default* is 500
    edi_list             list of .edi files to invert for
    grid_east            overall distance of grid nodes in east direction 
    grid_north           overall distance of grid nodes in north direction 
    grid_z               overall distance of grid nodes in z direction 
    initial_fn           full path to initial file name
    n_layers             total number of vertical layers in model
    nodes_east           relative distance between nodes in east direction 
    nodes_north          relative distance between nodes in north direction 
    nodes_z              relative distance between nodes in east direction 
    pad_east             number of cells for padding on E and W sides
                         *default* is 5
    pad_north            number of cells for padding on S and N sides
                         *default* is 5
    pad_root_east        padding cells E & W will be pad_root_east**(x)
    pad_root_north       padding cells N & S will be pad_root_north**(x) 
    pad_z                number of cells for padding at bottom
                         *default* is 5
    res_list             list of resistivity values for starting model
    res_model            starting resistivity model
    rotation_angle       Angle to rotate the grid to. Angle is measured
                         positve clockwise assuming North is 0 and east is 90.
                         *default* is None
    save_path            path to save file to  
    station_fn           full path to station file
    station_locations    location of stations
    title                title in initial file
    z1_layer             first layer thickness
    z_bottom             absolute bottom of the model *default* is 300,000 
    z_target_depth       Depth of deepest target, *default* is 50,000
    ==================== ======================================================

    ==================== ======================================================
    Methods              Description
    ==================== ======================================================
    make_mesh            makes a mesh from the given specifications
    plot_mesh            plots mesh to make sure everything is good
    write_initial_file   writes an initial model file that includes the mesh
    ==================== ======================================================


    """

    def __init__(self, edi_list=None, **kwargs):

        self.edi_list = edi_list

        # size of cells within station area in meters
        self.cell_size_east = kwargs.pop('cell_size_east', 500)
        self.cell_size_north = kwargs.pop('cell_size_north', 500)

        #padding cells on either side
        self.pad_east = kwargs.pop('pad_east', 5)
        self.pad_north = kwargs.pop('pad_north', 5)
        self.pad_z = kwargs.pop('pad_z', 5)

        #root of padding cells
        self.pad_root_east = kwargs.pop('pad_root_east', 5)
        self.pad_root_north = kwargs.pop('pad_root_north', 5)

        self.z1_layer = kwargs.pop('z1_layer', 10)
        self.z_target_depth = kwargs.pop('z_target_depth', 50000)
        self.z_bottom = kwargs.pop('z_bottom', 300000)

        #number of vertical layers
        self.n_layers = kwargs.pop('n_layers', 30)

        #--> attributes to be calculated
        #station information
        self.station_locations = kwargs.pop('station_locations', None)

       #grid nodes
        self.nodes_east = None
        self.nodes_north = None
        self.nodes_z = None

        #grid locations
        self.grid_east = None
        self.grid_north = None
        self.grid_z = None

        #resistivity model
        self.res_model = None
        self.res_list = None
        self.res_model_int = None

        #rotation angle
        self.rotation_angle = kwargs.pop('rotation_angle', 0.0)
        
        #inital file stuff
        self.initial_fn = None
        self.station_fn = None
        self.save_path = kwargs.pop('save_path', None)
        self.title = 'Inital Model File made in MTpy'


    def make_mesh(self):
        """ 
        create finite element mesh according to parameters set.

        The mesh is built by first finding the center of the station area.  
        Then cells are added in the north and east direction with width
        cell_size_east and cell_size_north to the extremeties of the station 
        area.  Padding cells are then added to extend the model to reduce 
        edge effects.  The number of cells are pad_east and pad_north and the
        increase in size is by pad_root_east and pad_root_north.  The station
        locations are then computed as the center of the nearest cell as 
        required by the code.

        The vertical cells are built to increase in size exponentially with
        depth.  The first cell depth is first_layer_thickness and should be
        about 1/10th the shortest skin depth.  The layers then increase
        on a log scale to z_target_depth.  Then the model is
        padded with pad_z number of cells to extend the depth of the model.

        padding = np.round(cell_size_east*pad_root_east**np.arange(start=.5,
                           stop=3, step=3./pad_east))+west 

        """
        #if station locations are not input read from the edi files
        if self.station_locations is None:
            if self.edi_list is None:
                raise AttributeError('edi_list is None, need to input a list of '
                                     'edi files to read in.')

            n_stations = len(self.edi_list)

            #make a structured array to put station location information into
            self.station_locations = np.zeros(n_stations,
                                              dtype=[('station','|S10'),
                                                     ('east', np.float),
                                                     ('north', np.float),
                                                     ('east_c', np.float),
                                                     ('north_c', np.float),
                                                     ('elev', np.float)])
            #get station locations in meters
            for ii, edi in enumerate(self.edi_list):
                mt_obj = mt.MT(edi)
                self.station_locations[ii]['station'] = mt_obj.station
                self.station_locations[ii]['east'] = mt_obj.east
                self.station_locations[ii]['north'] = mt_obj.north
                self.station_locations[ii]['elev'] = mt_obj.elev




            #--> rotate grid if necessary
            #to do this rotate the station locations because ModEM assumes the
            #input mesh is a lateral grid.
            #needs to be 90 - because North is assumed to be 0 but the rotation
            #matrix assumes that E is 0.
            if self.rotation_angle != 0:
                cos_ang = np.cos(np.deg2rad(self.rotation_angle))
                sin_ang = np.sin(np.deg2rad(self.rotation_angle))
                rot_matrix = np.matrix(np.array([[cos_ang, sin_ang],
                                                 [-sin_ang, cos_ang]]))

                coords = np.array([self.station_locations['east'],
                                   self.station_locations['north']])

                #rotate the relative station locations
                new_coords = np.array(np.dot(rot_matrix, coords))

                self.station_locations['east'][:] = new_coords[0, :]
                self.station_locations['north'][:] = new_coords[1, :]

                print('Rotated stations by {0:.1f} deg clockwise from N'.format(
                                                        self.rotation_angle))
            #remove the average distance to get coordinates in a relative space
            self.station_locations['east'] -= self.station_locations['east'].mean()
            self.station_locations['north'] -= self.station_locations['north'].mean()

            #translate the stations so they are relative to 0,0
            east_center = (self.station_locations['east'].max()-
                            np.abs(self.station_locations['east'].min()))/2
            north_center = (self.station_locations['north'].max()-
                            np.abs(self.station_locations['north'].min()))/2

            #remove the average distance to get coordinates in a relative space
            self.station_locations['east'] -= east_center
            self.station_locations['north'] -= north_center

        #pickout the furtherst south and west locations 
        #and put that station as the bottom left corner of the main grid
        west = self.station_locations['east'].min()-(1.5*self.cell_size_east)
        east = self.station_locations['east'].max()+(1.5*self.cell_size_east)
        south = self.station_locations['north'].min()-(1.5*self.cell_size_north)
        north = self.station_locations['north'].max()+(1.5*self.cell_size_north)

        #make sure the variable n_stations is initialized
        try:
            n_stations
        except NameError:
            n_stations = self.station_locations.shape[0]

        #-------make a grid around the stations from the parameters above------
        #--> make grid in east-west direction
        #cells within station area
        midxgrid = np.arange(start=west,
                             stop=east+self.cell_size_east,
                             step=self.cell_size_east)

        #padding cells on the west side
        pad_west = np.round(-self.cell_size_east*\
                             self.pad_root_east**np.arange(start=.5, stop=3,
                             step=3./self.pad_east))+west

        #padding cells on east side
        pad_east = np.round(self.cell_size_east*\
                             self.pad_root_east**np.arange(start=.5, stop=3,
                             step=3./self.pad_east))+east

        #make the cells going west go in reverse order and append them to the
        #cells going east
        east_gridr = np.append(np.append(pad_west[::-1], midxgrid), pad_east)

        #--> make grid in north-south direction
        #N-S cells with in station area
        midygrid = np.arange(start=south,
                             stop=north+self.cell_size_north,
                             step=self.cell_size_north)

        #padding cells on south side
        south_pad = np.round(-self.cell_size_north*
                              self.pad_root_north**np.arange(start=.5,
                              stop=3, step=3./self.pad_north))+south

        #padding cells on north side
        north_pad = np.round(self.cell_size_north*
                              self.pad_root_north**np.arange(start=.5,
                              stop=3, step=3./self.pad_north))+north

        #make the cells going west go in reverse order and append them to the
        #cells going east
        north_gridr = np.append(np.append(south_pad[::-1], midygrid), north_pad)


        #--> make depth grid
        log_z = np.logspace(np.log10(self.z1_layer),
                            np.log10(self.z_target_depth-np.logspace(np.log10(self.z1_layer),
                            np.log10(self.z_target_depth),
                            num=self.n_layers)[-2]),
                            num=self.n_layers-self.pad_z)
        ztarget = np.array([zz-zz%10**np.floor(np.log10(zz)) for zz in
                           log_z])
        log_zpad = np.logspace(np.log10(self.z_target_depth),
                            np.log10(self.z_bottom-np.logspace(np.log10(self.z_target_depth),
                            np.log10(self.z_bottom),
                            num=self.pad_z)[-2]),
                            num=self.pad_z)
        zpadding = np.array([zz-zz%10**np.floor(np.log10(zz)) for zz in
                               log_zpad])

        z_nodes = np.append(ztarget, zpadding)
        z_grid = np.array([z_nodes[:ii+1].sum() for ii in range(z_nodes.shape[0])])

        #---Need to make an array of the individual cell dimensions for
        #   wsinv3d
        east_nodes = east_gridr.copy()
        nx = east_gridr.shape[0]
        east_nodes[:nx/2] = np.array([abs(east_gridr[ii]-east_gridr[ii+1])
                                          for ii in range(int(nx/2))])
        east_nodes[nx/2:] = np.array([abs(east_gridr[ii]-east_gridr[ii+1])
                                          for ii in range(int(nx/2)-1, nx-1)])

        north_nodes = north_gridr.copy()
        ny = north_gridr.shape[0]
        north_nodes[:ny/2] = np.array([abs(north_gridr[ii]-north_gridr[ii+1])
                                       for ii in range(int(ny/2))])
        north_nodes[ny/2:] = np.array([abs(north_gridr[ii]-north_gridr[ii+1])
                                       for ii in range(int(ny/2)-1, ny-1)])

        #--put the grids into coordinates relative to the center of the grid
        east_grid = east_nodes.copy()
        east_grid[:int(nx/2)] = -np.array([east_nodes[ii:int(nx/2)].sum()
                                           for ii in range(int(nx/2))])
        east_grid[int(nx/2):] = np.array([east_nodes[int(nx/2):ii+1].sum()
                                         for ii in range(int(nx/2), nx)])-\
                                         east_nodes[int(nx/2)]

        north_grid = north_nodes.copy()
        north_grid[:int(ny/2)] = -np.array([north_nodes[ii:int(ny/2)].sum()
                                            for ii in range(int(ny/2))])
        north_grid[int(ny/2):] = np.array([north_nodes[int(ny/2):ii+1].sum()
                                            for ii in range(int(ny/2),ny)])-\
                                            north_nodes[int(ny/2)]

        #make nodes attributes
        self.nodes_east = east_nodes
        self.nodes_north = north_nodes
        self.nodes_z = z_nodes
        self.grid_east = east_grid
        self.grid_north = north_grid
        self.grid_z = z_grid

        #make sure that the stations are in the center of the cell as requested
        #by the code.
        for ii in range(n_stations):
            #look for the closest grid line
            xx = [nn for nn, xf in enumerate(east_grid)
                if xf>(self.station_locations[ii]['east']-self.cell_size_east)
                and xf<(self.station_locations[ii]['east']+self.cell_size_east)]

            #shift the station to the center in the east-west direction
            if east_grid[xx[0]] < self.station_locations[ii]['east']:
                self.station_locations[ii]['east_c'] = \
                                        east_grid[xx[0]]+self.cell_size_east/2
            elif east_grid[xx[0]] > self.station_locations[ii]['east']:
                self.station_locations[ii]['east_c'] = \
                                        east_grid[xx[0]]-self.cell_size_east/2

            #look for closest grid line
            yy = [mm for mm, yf in enumerate(north_grid)
                 if yf>(self.station_locations[ii]['north']-self.cell_size_north)
                 and yf<(self.station_locations[ii]['north']+self.cell_size_north)]

            #shift station to center of cell in north-south direction
            if north_grid[yy[0]] < self.station_locations[ii]['north']:
                self.station_locations[ii]['north_c'] = \
                                    north_grid[yy[0]]+self.cell_size_north/2
            elif north_grid[yy[0]] > self.station_locations[ii]['north']:
                self.station_locations[ii]['north_c'] = \
                                    north_grid[yy[0]]-self.cell_size_north/2

        #--> print out useful information
        print('-'*15)
        print('   Number of stations = {0}'.format(len(self.station_locations)))
        print('   Dimensions: ')
        print('      e-w = {0}'.format(east_grid.shape[0]))
        print('      n-s = {0}'.format(north_grid.shape[0]))
        print('       z  = {0} (without 7 air layers)'.format(z_grid.shape[0]))
        print('   Extensions: ')
        print('      e-w = {0:.1f} (m)'.format(east_nodes.__abs__().sum()))
        print('      n-s = {0:.1f} (m)'.format(north_nodes.__abs__().sum()))
        print('      0-z = {0:.1f} (m)'.format(self.nodes_z.__abs__().sum()))
        print('-'*15)

        #write a station location file for later
        stations = WSStation()
        stations.write_station_file(east=self.station_locations['east_c'],
                                    north=self.station_locations['north_c'],
                                    elev=self.station_locations['elev'],
                                    station_list=self.station_locations['station'],
                                    save_path=self.save_path)
        self.station_fn = stations.station_fn


    def plot_mesh(self, east_limits=None, north_limits=None, z_limits=None,
                  **kwargs):
        """

        Arguments:
        ----------
            **east_limits** : tuple (xmin,xmax)
                             plot min and max distances in meters for the 
                             E-W direction.  If None, the east_limits
                             will be set to furthest stations east and west.
                             *default* is None

            **north_limits** : tuple (ymin,ymax)
                             plot min and max distances in meters for the 
                             N-S direction.  If None, the north_limits
                             will be set to furthest stations north and south.
                             *default* is None

            **z_limits** : tuple (zmin,zmax)
                            plot min and max distances in meters for the 
                            vertical direction.  If None, the z_limits is
                            set to the number of layers.  Z is positive down
                            *default* is None
        """

        fig_size = kwargs.pop('fig_size', [6, 6])
        fig_dpi = kwargs.pop('fig_dpi', 300)
        fig_num = kwargs.pop('fig_num', 1)

        station_marker = kwargs.pop('station_marker', 'v')
        marker_color = kwargs.pop('station_color', 'b')
        marker_size = kwargs.pop('marker_size', 2)

        line_color = kwargs.pop('line_color', 'k')
        line_width = kwargs.pop('line_width', .5)

        plt.rcParams['figure.subplot.hspace'] = .3
        plt.rcParams['figure.subplot.wspace'] = .3
        plt.rcParams['figure.subplot.left'] = .08
        plt.rcParams['font.size'] = 7

        fig = plt.figure(fig_num, figsize=fig_size, dpi=fig_dpi)
        plt.clf()

        #---plot map view
        ax1 = fig.add_subplot(1, 2, 1, aspect='equal')

        #make sure the station is in the center of the cell
        ax1.scatter(self.station_locations['east_c'],
                    self.station_locations['north_c'],
                    marker=station_marker,
                    c=marker_color,
                    s=marker_size)

        #plot the grid if desired
        east_line_xlist = []
        east_line_ylist = []
        for xx in self.grid_east:
            east_line_xlist.extend([xx, xx])
            east_line_xlist.append(None)
            east_line_ylist.extend([self.grid_north.min(),
                                    self.grid_north.max()])
            east_line_ylist.append(None)
        ax1.plot(east_line_xlist,
                      east_line_ylist,
                      lw=line_width,
                      color=line_color)

        north_line_xlist = []
        north_line_ylist = []
        for yy in self.grid_north:
            north_line_xlist.extend([self.grid_east.min(),
                                     self.grid_east.max()])
            north_line_xlist.append(None)
            north_line_ylist.extend([yy, yy])
            north_line_ylist.append(None)
        ax1.plot(north_line_xlist,
                      north_line_ylist,
                      lw=line_width,
                      color=line_color)

        if east_limits == None:
            ax1.set_xlim(self.station_locations['east'].min()-\
                            10*self.cell_size_east,
                         self.station_locations['east'].max()+\
                             10*self.cell_size_east)
        else:
            ax1.set_xlim(east_limits)

        if north_limits == None:
            ax1.set_ylim(self.station_locations['north'].min()-\
                            10*self.cell_size_north,
                         self.station_locations['north'].max()+\
                             10*self.cell_size_east)
        else:
            ax1.set_ylim(north_limits)

        ax1.set_ylabel('Northing (m)', fontdict={'size':9,'weight':'bold'})
        ax1.set_xlabel('Easting (m)', fontdict={'size':9,'weight':'bold'})

        ##----plot depth view
        ax2 = fig.add_subplot(1, 2, 2, aspect='auto', sharex=ax1)


        #plot the grid if desired
        east_line_xlist = []
        east_line_ylist = []
        for xx in self.grid_east:
            east_line_xlist.extend([xx, xx])
            east_line_xlist.append(None)
            east_line_ylist.extend([0,
                                    self.grid_z.max()])
            east_line_ylist.append(None)
        ax2.plot(east_line_xlist,
                 east_line_ylist,
                 lw=line_width,
                 color=line_color)

        z_line_xlist = []
        z_line_ylist = []
        for zz in self.grid_z:
            z_line_xlist.extend([self.grid_east.min(),
                                     self.grid_east.max()])
            z_line_xlist.append(None)
            z_line_ylist.extend([zz, zz])
            z_line_ylist.append(None)
        ax2.plot(z_line_xlist,
                 z_line_ylist,
                 lw=line_width,
                 color=line_color)


        #--> plot stations
        ax2.scatter(self.station_locations['east_c'],
                    [0]*self.station_locations.shape[0],
                    marker=station_marker,
                    c=marker_color,
                    s=marker_size)


        if z_limits == None:
            ax2.set_ylim(self.z_target_depth, -200)
        else:
            ax2.set_ylim(z_limits)

        if east_limits == None:
            ax1.set_xlim(self.station_locations['east'].min()-\
                            10*self.cell_size_east,
                         self.station_locations['east'].max()+\
                             10*self.cell_size_east)
        else:
            ax1.set_xlim(east_limits)

        ax2.set_ylabel('Depth (m)', fontdict={'size':9, 'weight':'bold'})
        ax2.set_xlabel('Easting (m)', fontdict={'size':9, 'weight':'bold'})

        plt.show()

    def convert_model_to_int(self):
        """
        convert the resistivity model that is in ohm-m to integer values
        corresponding to res_list

        """

        self.res_model_int = np.ones_like(self.res_model)
        #make a dictionary of values to write to file.
        self.res_dict = dict([(res, ii)
                              for ii, res in
                              enumerate(sorted(self.res_list), 1)])

        for ii, res in enumerate(self.res_list):
            indexes = np.where(self.res_model == res)
            self.res_model_int[indexes] = self.res_dict[res]
            if ii == 0:
                indexes = np.where(self.res_model <= res)
                self.res_model_int[indexes] = self.res_dict[res]
            elif ii == len(self.res_list)-1:
                indexes = np.where(self.res_model >= res)
                self.res_model_int[indexes] = self.res_dict[res]
            else:
                l_index = max([0, ii-1])
                h_index = min([len(self.res_list)-1, ii+1])
                indexes = np.where((self.res_model > self.res_list[l_index]) &
                                   (self.res_model < self.res_list[h_index]))
                self.res_model_int[indexes] = self.res_dict[res]

        print('Converted resistivity model to integers.')

    def write_initial_file(self, **kwargs):
        """
        will write an initial file for wsinv3d.  

        Note that x is assumed to be S --> N, y is assumed to be W --> E and
        z is positive downwards.  This means that index [0, 0, 0] is the 
        southwest corner of the first layer.  Therefore if you build a model
        by hand the layer block will look as it should in map view. 

        Also, the xgrid, ygrid and zgrid are assumed to be the relative 
        distance between neighboring nodes.  This is needed because wsinv3d 
        builds the  model from the bottom SW corner assuming the cell width
        from the init file.



        Key Word Arguments:
        ----------------------

            **nodes_north** : np.array(nx)
                        block dimensions (m) in the N-S direction. 
                        **Note** that the code reads the grid assuming that
                        index=0 is the southern most point.

            **nodes_east** : np.array(ny)
                        block dimensions (m) in the E-W direction.  
                        **Note** that the code reads in the grid assuming that
                        index=0 is the western most point.

            **nodes_z** : np.array(nz)
                        block dimensions (m) in the vertical direction.  
                        This is positive downwards.

            **save_path** : string
                          Path to where the initial file will be saved
                          to savepath/init3d

            **res_list** : float or list
                        The start resistivity as a float or a list of
                        resistivities that coorespond to the starting
                        resistivity model **res_model_int**.  
                        This must be input if you input **res_model_int**
                        If res_list is None, then Nr = 0 and the real 
                        resistivity values are used from **res_model**.
                        *default* is 100

            **title** : string
                        Title that goes into the first line of savepath/init3d

            **res_model** : np.array((nx,ny,nz))
                        Starting resistivity model.  Each cell is allocated an
                        integer value that cooresponds to the index value of
                        **res_list**.  .. note:: again that the modeling code 
                        assumes that the first row it reads in is the southern
                        most row and the first column it reads in is the 
                        western most column.  Similarly, the first plane it 
                        reads in is the Earth's surface.

            **res_model_int** : np.array((nx,ny,nz))
                        Starting resistivity model.  Each cell is allocated an
                        linear resistivity value.
                        .. note:: again that the modeling code 
                        assumes that the first row it reads in is the southern
                        most row and the first column it reads in is the 
                        western most column.  Similarly, the first plane it 
                        reads in is the Earth's surface.



        """
        keys = ['nodes_east', 'nodes_north', 'nodes_z', 'title', 'res_list',
                'res_model', 'res_model_int', 'save_path', 'initial_fn']
        for key in keys:
            try:
                setattr(self, key, kwargs[key])
            except KeyError:
                if self.__dict__[key] is None:
                    pass

        if self.initial_fn is None:
            if self.save_path is None:
                self.save_path = os.getcwd()
                self.initial_fn = os.path.join(self.save_path, "WSInitialModel")
            elif os.path.isdir(self.save_path) == True:
                self.initial_fn = os.path.join(self.save_path, "WSInitialModel")
            else:
                self.save_path = os.path.dirname(self.save_path)
                self.initial_fn= self.save_path

        #check to see what resistivity in input
        if self.res_list is None:
            nr = 0
        elif type(self.res_list) is not list and \
             type(self.res_list) is not np.ndarray:
            self.res_list = [self.res_list]
            nr = len(self.res_list)
        else:
            nr = len(self.res_list)

        #--> write file
        ifid = file(self.initial_fn, 'w')
        ifid.write('# {0}\n'.format(self.title.upper()))
        ifid.write('{0} {1} {2} {3}\n'.format(self.nodes_north.shape[0],
                                              self.nodes_east.shape[0],
                                              self.nodes_z.shape[0],
                                              nr))

        #write S --> N node block
        for ii, nnode in enumerate(self.nodes_north):
            ifid.write('{0:>12.1f}'.format(abs(nnode)))
            if ii != 0 and np.remainder(ii+1, 5) == 0:
                ifid.write('\n')
            elif ii == self.nodes_north.shape[0]-1:
                ifid.write('\n')

        #write W --> E node block
        for jj, enode in enumerate(self.nodes_east):
            ifid.write('{0:>12.1f}'.format(abs(enode)))
            if jj != 0 and np.remainder(jj+1, 5) == 0:
                ifid.write('\n')
            elif jj == self.nodes_east.shape[0]-1:
                ifid.write('\n')

        #write top --> bottom node block
        for kk, zz in enumerate(self.nodes_z):
            ifid.write('{0:>12.1f}'.format(abs(zz)))
            if kk != 0 and np.remainder(kk+1, 5) == 0:
                ifid.write('\n')
            elif kk == self.nodes_z.shape[0]-1:
                ifid.write('\n')

        #write the resistivity list
        if nr > 0:
            for ff in self.res_list:
                ifid.write('{0:.1f} '.format(ff))
            ifid.write('\n')
        else:
            pass

        if self.res_model == None:
            ifid.close()
        else:
            if nr > 0:
                if self.res_model_int is None:
                    self.convert_model_to_int()
                #need to flip the array such that the 1st index written is the
                #northern most value
                write_res_model = self.res_model_int[::-1, :, :]
                #get similar layers
            else:
                write_res_model = self.res_model[::-1, :, :]
            l1 = 0
            layers = []
            for zz in range(self.nodes_z.shape[0]-1):
                if (write_res_model[:, :, zz] ==
                    write_res_model[:, :, zz+1]).all() == False:
                    layers.append((l1, zz))
                    l1 = zz+1
            #need to add on the bottom layers
            layers.append((l1, self.nodes_z.shape[0]-1))

            #write out the layers from resmodel
            for ll in layers:
                ifid.write('{0} {1}\n'.format(ll[0]+1, ll[1]+1))
                for nn in range(self.nodes_north.shape[0]):
                    for ee in range(self.nodes_east.shape[0]):
                        if nr > 0:
                            ifid.write('{0:>3.0f}'.format(
                                          write_res_model[nn, ee, ll[0]]))
                        else:
                            ifid.write('{0:>8.1f}'.format(
                                          write_res_model[nn, ee, ll[0]]))
                    ifid.write('\n')
            ifid.close()

        print('Wrote file to: {0}'.format(self.initial_fn))

    def read_initial_file(self, initial_fn):
        """
        read an initial file and return the pertinent information including
        grid positions in coordinates relative to the center point (0,0) and 
        starting model.

        Arguments:
        ----------

            **initial_fn** : full path to initializing file.

        Outputs:
        --------

            **nodes_north** : np.array(nx)
                        array of nodes in S --> N direction

            **nodes_east** : np.array(ny) 
                        array of nodes in the W --> E direction

            **nodes_z** : np.array(nz)
                        array of nodes in vertical direction positive downwards

            **res_model** : dictionary
                        dictionary of the starting model with keys as layers

            **res_list** : list
                        list of resistivity values in the model

            **title** : string
                         title string

        """
        self.initial_fn = initial_fn
        ifid = file(self.initial_fn, 'r')
        ilines = ifid.readlines()
        ifid.close()

        self.title = ilines[0].strip()

        #get size of dimensions, remembering that x is N-S, y is E-W, z is + down
        nsize = ilines[1].strip().split()
        n_north = int(nsize[0])
        n_east = int(nsize[1])
        n_z = int(nsize[2])

        #initialize empy arrays to put things into
        self.nodes_north = np.zeros(n_north)
        self.nodes_east = np.zeros(n_east)
        self.nodes_z = np.zeros(n_z)
        self.res_model_int = np.zeros((n_north, n_east, n_z))
        self.res_model = np.zeros((n_north, n_east, n_z))

        #get the grid line locations
        line_index = 2       #line number in file
        count_n = 0  #number of north nodes found
        while count_n < n_north:
            iline = ilines[line_index].strip().split()
            for north_node in iline:
                self.nodes_north[count_n] = float(north_node)
                count_n += 1
            line_index += 1

        count_e = 0  #number of east nodes found
        while count_e < n_east:
            iline = ilines[line_index].strip().split()
            for east_node in iline:
                self.nodes_east[count_e] = float(east_node)
                count_e += 1
            line_index += 1

        count_z = 0  #number of vertical nodes
        while count_z < n_z:
            iline = ilines[line_index].strip().split()
            for z_node in iline:
                self.nodes_z[count_z] = float(z_node)
                count_z += 1
            line_index += 1

        #put the grids into coordinates relative to the center of the grid
        self.grid_north = self.nodes_north.copy()
        self.grid_north[:int(n_north/2)] =\
                        -np.array([self.nodes_north[ii:int(n_north/2)].sum()
                                   for ii in range(int(n_north/2))])
        self.grid_north[int(n_north/2):] = \
                        np.array([self.nodes_north[int(n_north/2):ii+1].sum()
                                 for ii in range(int(n_north/2), n_north)])-\
                                 self.nodes_north[int(n_north/2)]

        self.grid_east = self.nodes_east.copy()
        self.grid_east[:int(n_east/2)] = \
                            -np.array([self.nodes_east[ii:int(n_east/2)].sum()
                                       for ii in range(int(n_east/2))])
        self.grid_east[int(n_east/2):] = \
                            np.array([self.nodes_east[int(n_east/2):ii+1].sum()
                                     for ii in range(int(n_east/2),n_east)])-\
                                     self.nodes_east[int(n_east/2)]

        self.grid_z = np.array([self.nodes_z[:ii+1].sum() for ii in range(n_z)])

        #get the resistivity values
        self.res_list = [float(rr) for rr in ilines[line_index].strip().split()]
        line_index += 1

        #get model
        try:
            iline = ilines[line_index].strip().split()

        except IndexError:
            self.res_model[:, :, :] = self.res_list[0]
            self.res_model[:, :, :] = 1
            return

        if len(iline) == 0 or len(iline) == 1:
            self.res_model[:, :, :] = self.res_list[0]
            self.res_model_int[:, :, :] = 1
            return
        else:
            while line_index < len(ilines):
                iline = ilines[line_index].strip().split()
                if len(iline) == 2:
                    l1 = int(iline[0])-1
                    l2 = int(iline[1])
                    if l1 == l2:
                        l2 += 1
                    line_index += 1
                    count_n = 0
                elif len(iline) == 0:
                    break
                else:
                    count_e = 0
                    while count_e < n_east:
                        #be sure the indes of res list starts at 0 not 1 as
                        #in ws3dinv
                        self.res_model[count_n, count_e, l1:l2] =\
                                        self.res_list[int(iline[count_e])-1]
                        self.res_model_int[count_n, count_e, l1:l2] =\
                                            int(iline[count_e])
                        count_e += 1
                    count_n += 1
                    line_index += 1
            # Need to be sure that the resistivity array matches
            # with the grids, such that the first index is the
            # furthest south, even though ws3dinv outputs as first
            # index as furthest north.
            self.res_model = self.res_model[::-1, :, :]
            self.res_model_int = self.res_model_int[::-1, :, :]

#==============================================================================
# model class
#==============================================================================
class WSModel(object):
    """
    Reads in model file and fills necessary attributes.

    :Example: ::

        >>> mfn = r"/home/ws3dinv/test_model.00"
        >>> wsmodel = ws.WSModel(mfn)
        >>> wsmodel.write_vtk_file(r"/home/ParaviewFiles")

    ======================= ===================================================
    Attributes              Description
    ======================= ===================================================
    grid_east               overall distance of grid nodes in east direction 
    grid_north              overall distance of grid nodes in north direction 
    grid_z                  overall distance of grid nodes in z direction
    iteration_number        iteration number of the inversion
    lagrange                lagrange multiplier
    model_fn                full path to model file
    nodes_east              relative distance between nodes in east direction 
    nodes_north             relative distance between nodes in north direction 
    nodes_z                 relative distance between nodes in east direction 
    res_model               starting resistivity model
    rms                     root mean squared error of data and model
    ======================= ===================================================


    ======================= ===================================================
    Methods                 Description    
    ======================= ===================================================
    read_model_file         read model file and fill attributes
    write_vtk_file          write a vtk structured grid file for resistivity
                            model
    ======================= ===================================================

    """

    def __init__(self, model_fn=None):
        self.model_fn = model_fn
        self.iteration_number = None
        self.rms = None
        self.lagrange = None
        self.res_model = None

        self.nodes_north = None
        self.nodes_east = None
        self.nodes_z = None

        self.grid_north = None
        self.grid_east = None
        self.grid_z = None

        if self.model_fn is not None and os.path.isfile(self.model_fn) == True:
            self.read_model_file()

    def read_model_file(self):
        """
        read in a model file as x-north, y-east, z-positive down
        """

        mfid = file(self.model_fn, 'r')
        mlines = mfid.readlines()
        mfid.close()

        #get info at the beggining of file
        info = mlines[0].strip().split()
        self.iteration_number = int(info[2])
        self.rms = float(info[5])
        try:
            self.lagrange = float(info[8])
        except IndexError:
            print('Did not get Lagrange Multiplier')

        #get lengths of things
        n_north, n_east, n_z, n_res = np.array(mlines[1].strip().split(),
                                               dtype=np.int)

        #make empty arrays to put stuff into
        self.nodes_north = np.zeros(n_north)
        self.nodes_east = np.zeros(n_east)
        self.nodes_z = np.zeros(n_z)
        self.res_model = np.zeros((n_north, n_east, n_z))

        #get the grid line locations
        line_index = 2       #line number in file
        count_n = 0  #number of north nodes found
        while count_n < n_north:
            mline = mlines[line_index].strip().split()
            for north_node in mline:
                self.nodes_north[count_n] = float(north_node)
                count_n += 1
            line_index += 1

        count_e = 0  #number of east nodes found
        while count_e < n_east:
            mline = mlines[line_index].strip().split()
            for east_node in mline:
                self.nodes_east[count_e] = float(east_node)
                count_e += 1
            line_index += 1

        count_z = 0  #number of vertical nodes
        while count_z < n_z:
            mline = mlines[line_index].strip().split()
            for z_node in mline:
                self.nodes_z[count_z] = float(z_node)
                count_z += 1
            line_index += 1

        #put the grids into coordinates relative to the center of the grid
        self.grid_north = self.nodes_north.copy()
        self.grid_north[:int(n_north/2)] =\
                        -np.array([self.nodes_north[ii:int(n_north/2)].sum()
                                   for ii in range(int(n_north/2))])
        self.grid_north[int(n_north/2):] = \
                        np.array([self.nodes_north[int(n_north/2):ii+1].sum()
                                 for ii in range(int(n_north/2), n_north)])-\
                                 self.nodes_north[int(n_north/2)]

        self.grid_east = self.nodes_east.copy()
        self.grid_east[:int(n_east/2)] = \
                            -np.array([self.nodes_east[ii:int(n_east/2)].sum()
                                       for ii in range(int(n_east/2))])
        self.grid_east[int(n_east/2):] = \
                            np.array([self.nodes_east[int(n_east/2):ii+1].sum()
                                     for ii in range(int(n_east/2),n_east)])-\
                                     self.nodes_east[int(n_east/2)]

        self.grid_z = np.array([self.nodes_z[:ii+1].sum() for ii in range(n_z)])

        #--> get resistivity values
        #need to read in the north backwards so that the first index is
        #southern most point
        for kk in range(n_z):
            for jj in range(n_east):
                for ii in range(n_north):
                    self.res_model[(n_north-1)-ii, jj, kk] = \
                                             float(mlines[line_index].strip())
                    line_index += 1

    def write_vtk_file(self, save_fn):
        """

        """
        if os.path.isdir(save_fn) == True:
            save_fn = os.path.join(save_fn, 'VTKResistivity_Model')

        save_fn = gridToVTK(save_fn,
                            self.grid_north,
                            self.grid_east,
                            self.grid_z,
                            cellData={'resistivity':self.res_model})

        print('Wrote vtk file to {0}'.format(save_fn))


#==============================================================================
# Manipulate the model
#==============================================================================
class WSModelManipulator(object):
    """
    will plot a model from wsinv3d or init file so the user can manipulate the 
    resistivity values relatively easily.  At the moment only plotted
    in map view.


    :Example: ::
        >>> import mtpy.modeling.ws3dinv as ws
        >>> initial_fn = r"/home/MT/ws3dinv/Inv1/WSInitialFile"
        >>> mm = ws.WSModelManipulator(initial_fn=initial_fn)

    =================== =======================================================
    Buttons              Description    
    =================== =======================================================
    '='                 increase depth to next vertical node (deeper)
    '-'                 decrease depth to next vertical node (shallower)
    'q'                 quit the plot, rewrites initial file when pressed
    'a'                 copies the above horizontal layer to the present layer
    'b'                 copies the below horizonal layer to present layer
    'u'                 undo previous change
    =================== =======================================================


    =================== =======================================================
    Attributes          Description
    =================== =======================================================
    ax1                 matplotlib.axes instance for mesh plot of the model 
    ax2                 matplotlib.axes instance of colorbar
    cb                  matplotlib.colorbar instance for colorbar 
    cid_depth           matplotlib.canvas.connect for depth
    cmap                matplotlib.colormap instance
    cmax                maximum value of resistivity for colorbar. (linear)
    cmin                minimum value of resistivity for colorbar (linear)
    data_fn             full path fo data file
    depth_index         integer value of depth slice for plotting
    dpi                 resolution of figure in dots-per-inch
    dscale              depth scaling, computed internally
    east_line_xlist     list of east mesh lines for faster plotting
    east_line_ylist     list of east mesh lines for faster plotting
    fdict               dictionary of font properties
    fig                 matplotlib.figure instance
    fig_num              number of figure instance
    fig_size             size of figure in inches
    font_size           size of font in points
    grid_east           location of east nodes in relative coordinates
    grid_north          location of north nodes in relative coordinates
    grid_z              location of vertical nodes in relative coordinates
    initial_fn          full path to initial file
    m_height            mean height of horizontal cells
    m_width             mean width of horizontal cells
    map_scale            [ 'm' | 'km' ] scale of map
    mesh_east           np.meshgrid of east, north
    mesh_north          np.meshgrid of east, north
    mesh_plot           matplotlib.axes.pcolormesh instance
    model_fn            full path to model file
    new_initial_fn      full path to new initial file
    nodes_east          spacing between east nodes 
    nodes_north         spacing between north nodes 
    nodes_z             spacing between vertical nodes
    north_line_xlist    list of coordinates of north nodes for faster plotting
    north_line_ylist    list of coordinates of north nodes for faster plotting
    plot_yn             [ 'y' | 'n' ] plot on instantiation
    radio_res           matplotlib.widget.radio instance for change resistivity
    rect_selector       matplotlib.widget.rect_selector 
    res                 np.ndarray(nx, ny, nz) for model in linear resistivity
    res_copy            copy of res for undo
    res_dict            dictionary of segmented resistivity values 
    res_list            list of resistivity values for model linear scale
    res_model           np.ndarray(nx, ny, nz) of resistivity values from 
                        res_list (linear scale)
    res_model_int       np.ndarray(nx, ny, nz) of integer values corresponding
                        to res_list for initial model
    res_value           current resistivty value of radio_res
    save_path           path to save initial file to
    station_east        station locations in east direction
    station_north       station locations in north direction
    xlimits             limits of plot in e-w direction
    ylimits             limits of plot in n-s direction
    =================== =======================================================

    """

    def __init__(self, model_fn=None, initial_fn=None, data_fn=None, **kwargs):

        self.model_fn = model_fn
        self.initial_fn = initial_fn
        self.data_fn = data_fn
        self.new_initial_fn = None
        self.initial_fn_basename = kwargs.pop('initial_fn_basename',
                                              'WSInitialModel_mm')

        if self.model_fn is not None:
            self.save_path = os.path.dirname(self.model_fn)
        elif self.initial_fn is not None:
            self.save_path = os.path.dirname(self.initial_fn)
        elif self.data_fn is not None:
            self.save_path = os.path.dirname(self.data_fn)
        else:
            self.save_path = None

        #grid nodes
        self.nodes_east = None
        self.nodes_north = None
        self.nodes_z = None

        #grid locations
        self.grid_east = None
        self.grid_north = None
        self.grid_z = None

        #resistivity model
        self.res_model_int = None #model in ints
        self.res_model = None     #model in floats
        self.res = None

        #station locations in relative coordinates read from data file
        self.station_east = None
        self.station_north = None

        #--> set map scale
        self.map_scale = kwargs.pop('map_scale', 'km')

        self.m_width = 100
        self.m_height = 100

        #--> scale the map coordinates
        if self.map_scale=='km':
            self.dscale = 1000.
        if self.map_scale=='m':
            self.dscale = 1.

        #figure attributes
        self.fig = None
        self.ax1 = None
        self.ax2 = None
        self.cb = None
        self.east_line_xlist = None
        self.east_line_ylist = None
        self.north_line_xlist = None
        self.north_line_ylist = None

        #make a default resistivity list to change values
        self.res_dict = None
        self.res_list = kwargs.pop('res_list', None)
        if self.res_list is None:
            self.set_res_list(np.array([.3, 1, 10, 50, 100, 500, 1000, 5000],
                                      dtype=np.float))

        else:
            try:
                if len(self.res_list) > 10:
                    print ('!! Warning -- ws3dinv can only deal with 10 '
                           'resistivity values for the initial model')
            except TypeError:
                self.res_list = [self.res_list]
            self.set_res_list(self.res_list)


        #read in model or initial file
        self.read_file()

        #set initial resistivity value
        self.res_value = self.res_list[0]

        #--> set map limits
        self.xlimits = kwargs.pop('xlimits', None)
        self.ylimits = kwargs.pop('ylimits', None)

        self.font_size = kwargs.pop('font_size', 7)
        self.fig_dpi = kwargs.pop('fig_dpi', 300)
        self.fig_num = kwargs.pop('fig_num', 1)
        self.fig_size = kwargs.pop('fig_size', [6, 6])
        self.cmap = kwargs.pop('cmap', cm.jet_r)
        self.depth_index = kwargs.pop('depth_index', 0)

        self.fdict = {'size':self.font_size+2, 'weight':'bold'}

        #plot on initialization
        self.plot_yn = kwargs.pop('plot_yn', 'y')
        if self.plot_yn=='y':
            self.plot()

    def set_res_list(self, res_list):
        """
        on setting res_list also set the res_dict to correspond
        """
        self.res_list = res_list
        #make a dictionary of values to write to file.
        self.res_dict = dict([(res, ii)
                              for ii, res in enumerate(self.res_list,1)])
        if self.fig is not None:
            plt.close()
            self.plot()


    #---read files-------------------------------------------------------------
    def read_file(self):
        """
        reads in initial file or model file and set attributes:
            -resmodel
            -northrid
            -eastrid
            -zgrid
            -res_list if initial file

        """
        att_names = ['nodes_north', 'nodes_east', 'nodes_z', 'grid_east',
                     'grid_north', 'grid_z', 'res_model', 'res_list']

        #--> read model file
        if self.model_fn is not None and self.initial_fn is None:

            wsmodel = WSModel(self.model_fn)
            wsmodel.read_model_file()

            for name in att_names:
                if hasattr(wsmodel, name):
                    value = getattr(wsmodel, name)
                    setattr(self, name, value)

            #--> scale the resistivity values from the model into
            #    a segmented scale that cooresponds to res_list
            self.convert_res_to_model(self.res_model.copy())

        #--> read initial file
        elif self.initial_fn is not None and self.model_fn is None:
            wsmesh = WSMesh()
            wsmesh.read_initial_file(self.initial_fn)
            for name in att_names:
                if hasattr(wsmesh, name):
                    value = getattr(wsmesh, name)
                    setattr(self, name, value)

            self.res_model_int = wsmesh.res_model
            if len(wsmesh.res_list) == 1:
                self.set_res_list([.3, 1, 10, 100, 1000])
            else:
                self.set_res_list(wsmesh.res_list)

            #need to convert index values to resistivity values
            rdict = dict([(ii,res) for ii,res in enumerate(self.res_list,1)])

            for ii in range(len(self.res_list)):
                self.res_model[np.where(self.res_model_int==ii+1)] = rdict[ii+1]

        elif self.initial_fn is None and self.model_fn is None:
            print('Need to input either an initial file or model file to plot')
        else:
            print('Input just initial file or model file not both.')

        #--> read in data file if given
        if self.data_fn is not None:
            wsdata = WSData()
            wsdata.read_data_file(self.data_fn)

            #get station locations
            self.station_east = wsdata.data['east']
            self.station_north = wsdata.data['north']

        #get cell block sizes
        self.m_height = np.median(self.nodes_north[5:-5])/self.dscale
        self.m_width = np.median(self.nodes_east[5:-5])/self.dscale

        #make a copy of original in case there are unwanted changes
        self.res_copy = self.res_model.copy()



    #---plot model-------------------------------------------------------------
    def plot(self):
        """
        plots the model with:
            -a radio dial for depth slice 
            -radio dial for resistivity value

        """

        self.cmin = np.floor(np.log10(min(self.res_list)))
        self.cmax = np.ceil(np.log10(max(self.res_list)))

        #-->Plot properties
        plt.rcParams['font.size'] = self.font_size

        #need to add an extra row and column to east and north to make sure
        #all is plotted see pcolor for details.
        plot_east = np.append(self.grid_east, self.grid_east[-1]*1.25)/self.dscale
        plot_north = np.append(self.grid_north, self.grid_north[-1]*1.25)/self.dscale

        #make a mesh grid for plotting
        #the 'ij' makes sure the resulting grid is in east, north
        self.mesh_east, self.mesh_north = np.meshgrid(plot_east,
                                                      plot_north,
                                                      indexing='ij')

        self.fig = plt.figure(self.fig_num, self.fig_size, dpi=self.fig_dpi)
        plt.clf()
        self.ax1 = self.fig.add_subplot(1, 1, 1, aspect='equal')

        #transpose to make x--east and y--north
        plot_res = np.log10(self.res_model[:,:,self.depth_index].T)

        self.mesh_plot = self.ax1.pcolormesh(self.mesh_east,
                                             self.mesh_north,
                                             plot_res,
                                             cmap=self.cmap,
                                             vmin=self.cmin,
                                             vmax=self.cmax)

        #on plus or minus change depth slice
        self.cid_depth = \
                    self.mesh_plot.figure.canvas.mpl_connect('key_press_event',
                                                        self._on_key_callback)


        #plot the stations
        if self.station_east is not None:
            for ee, nn in zip(self.station_east, self.station_north):
                self.ax1.text(ee/self.dscale, nn/self.dscale,
                              '*',
                              verticalalignment='center',
                              horizontalalignment='center',
                              fontdict={'size':self.font_size-2,
                                        'weight':'bold'})

        #set axis properties
        if self.xlimits is not None:
            self.ax1.set_xlim(self.xlimits)
        else:
            self.ax1.set_xlim(xmin=self.grid_east.min()/self.dscale,
                              xmax=self.grid_east.max()/self.dscale)

        if self.ylimits is not None:
            self.ax1.set_ylim(self.ylimits)
        else:
            self.ax1.set_ylim(ymin=self.grid_north.min()/self.dscale,
                              ymax=self.grid_north.max()/self.dscale)

        #self.ax1.xaxis.set_minor_locator(MultipleLocator(100*1./dscale))
        #self.ax1.yaxis.set_minor_locator(MultipleLocator(100*1./dscale))

        self.ax1.set_ylabel('Northing ('+self.map_scale+')',
                            fontdict=self.fdict)
        self.ax1.set_xlabel('Easting ('+self.map_scale+')',
                            fontdict=self.fdict)

        depth_title = self.grid_z[self.depth_index]/self.dscale

        self.ax1.set_title('Depth = {:.3f} '.format(depth_title)+\
                           '('+self.map_scale+')',
                           fontdict=self.fdict)

        #plot the grid if desired
        self.east_line_xlist = []
        self.east_line_ylist = []
        for xx in self.grid_east:
            self.east_line_xlist.extend([xx/self.dscale, xx/self.dscale])
            self.east_line_xlist.append(None)
            self.east_line_ylist.extend([self.grid_north.min()/self.dscale,
                                         self.grid_north.max()/self.dscale])
            self.east_line_ylist.append(None)
        self.ax1.plot(self.east_line_xlist,
                      self.east_line_ylist,
                       lw=.25,
                       color='k')

        self.north_line_xlist = []
        self.north_line_ylist = []
        for yy in self.grid_north:
            self.north_line_xlist.extend([self.grid_east.min()/self.dscale,
                                          self.grid_east.max()/self.dscale])
            self.north_line_xlist.append(None)
            self.north_line_ylist.extend([yy/self.dscale, yy/self.dscale])
            self.north_line_ylist.append(None)
        self.ax1.plot(self.north_line_xlist,
                      self.north_line_ylist,
                      lw=.25,
                      color='k')

        #plot the colorbar
        self.ax2 = mcb.make_axes(self.ax1, orientation='vertical', shrink=.35)
        seg_cmap = cmap_discretize(self.cmap, len(self.res_list))
        self.cb = mcb.ColorbarBase(self.ax2[0],cmap=seg_cmap,
                                   norm=colors.Normalize(vmin=self.cmin,
                                                         vmax=self.cmax))


        self.cb.set_label('Resistivity ($\Omega \cdot$m)',
                          fontdict={'size':self.font_size})
        self.cb.set_ticks(np.arange(self.cmin, self.cmax+1))
        self.cb.set_ticklabels([mtplottools.labeldict[cc]
                                for cc in np.arange(self.cmin, self.cmax+1)])

        #make a resistivity radio button
        resrb = self.fig.add_axes([.85,.1,.1,.2])
        reslabels = ['{0:.4g}'.format(res) for res in self.res_list]
        self.radio_res = widgets.RadioButtons(resrb, reslabels,
                                        active=self.res_dict[self.res_value])

        #make a rectangular selector
        self.rect_selector = widgets.RectangleSelector(self.ax1,
                                                       self.rect_onselect,
                                                       drawtype='box',
                                                       useblit=True)


        plt.show()

        #needs to go after show()
        self.radio_res.on_clicked(self.set_res_value)


    def redraw_plot(self):
        """
        redraws the plot
        """

        current_xlimits = self.ax1.get_xlim()
        current_ylimits = self.ax1.get_ylim()

        self.ax1.cla()

        plot_res = np.log10(self.res_model[:,:,self.depth_index].T)

        self.mesh_plot = self.ax1.pcolormesh(self.mesh_east,
                                             self.mesh_north,
                                             plot_res,
                                             cmap=self.cmap,
                                             vmin=self.cmin,
                                             vmax=self.cmax)

         #plot the stations
        if self.station_east is not None:
            for ee,nn in zip(self.station_east, self.station_north):
                self.ax1.text(ee/self.dscale, nn/self.dscale,
                              '*',
                              verticalalignment='center',
                              horizontalalignment='center',
                              fontdict={'size':self.font_size-2,
                                        'weight':'bold'})

        #set axis properties
        if self.xlimits is not None:
            self.ax1.set_xlim(self.xlimits)
        else:
            self.ax1.set_xlim(current_xlimits)

        if self.ylimits is not None:
            self.ax1.set_ylim(self.ylimits)
        else:
            self.ax1.set_ylim(current_ylimits)

        self.ax1.set_ylabel('Northing ('+self.map_scale+')',
                            fontdict=self.fdict)
        self.ax1.set_xlabel('Easting ('+self.map_scale+')',
                            fontdict=self.fdict)

        depth_title = self.grid_z[self.depth_index]/self.dscale

        self.ax1.set_title('Depth = {:.3f} '.format(depth_title)+\
                           '('+self.map_scale+')',
                           fontdict=self.fdict)

        #plot finite element mesh
        self.ax1.plot(self.east_line_xlist,
                      self.east_line_ylist,
                      lw=.25,
                      color='k')


        self.ax1.plot(self.north_line_xlist,
                      self.north_line_ylist,
                      lw=.25,
                      color='k')

        #be sure to redraw the canvas
        self.fig.canvas.draw()

    def set_res_value(self, label):
        self.res_value = float(label)
        print('set resistivity to ', label)
        print(self.res_value)


    def _on_key_callback(self,event):
        """
        on pressing a key do something

        """

        self.event_change_depth = event

        #go down a layer on push of +/= keys
        if self.event_change_depth.key == '=':
            self.depth_index += 1

            if self.depth_index>len(self.grid_z)-1:
                self.depth_index = len(self.grid_z)-1
                print('already at deepest depth')

            print('Plotting Depth {0:.3f}'.format(self.grid_z[self.depth_index]/\
                    self.dscale)+'('+self.map_scale+')')

            self.redraw_plot()
        #go up a layer on push of - key
        elif self.event_change_depth.key == '-':
            self.depth_index -= 1

            if self.depth_index < 0:
                self.depth_index = 0

            print('Plotting Depth {0:.3f} '.format(self.grid_z[self.depth_index]/\
                    self.dscale)+'('+self.map_scale+')')

            self.redraw_plot()

        #exit plot on press of q
        elif self.event_change_depth.key == 'q':
            self.event_change_depth.canvas.mpl_disconnect(self.cid_depth)
            plt.close(self.event_change_depth.canvas.figure)
            self.rewrite_initial_file()

        #copy the layer above
        elif self.event_change_depth.key == 'a':
            try:
                if self.depth_index == 0:
                    print('No layers above')
                else:
                    self.res_model[:, :, self.depth_index] = \
                                       self.res_model[:, :, self.depth_index-1]
            except IndexError:
                print('No layers above')

            self.redraw_plot()

        #copy the layer below
        elif self.event_change_depth.key == 'b':
            try:
                self.res_model[:, :, self.depth_index] = \
                                    self.res_model[:, :, self.depth_index+1]
            except IndexError:
                print('No more layers below')

            self.redraw_plot()

        #undo
        elif self.event_change_depth.key == 'u':
            if type(self.xchange) is int and type(self.ychange) is int:
                self.res_model[self.ychange, self.xchange, self.depth_index] =\
                self.res_copy[self.ychange, self.xchange, self.depth_index]
            else:
                for xx in self.xchange:
                    for yy in self.ychange:
                        self.res_model[yy, xx, self.depth_index] = \
                        self.res_copy[yy, xx, self.depth_index]

            self.redraw_plot()

    def change_model_res(self, xchange, ychange):
        """
        change resistivity values of resistivity model

        """
        if type(xchange) is int and type(ychange) is int:
            self.res_model[ychange, xchange, self.depth_index] = self.res_value
        else:
            for xx in xchange:
                for yy in ychange:
                    self.res_model[yy, xx, self.depth_index] = self.res_value

        self.redraw_plot()

    def rect_onselect(self, eclick, erelease):
        """
        on selecting a rectangle change the colors to the resistivity values
        """
        x1, y1 = eclick.xdata, eclick.ydata
        x2, y2 = erelease.xdata, erelease.ydata

        self.xchange = self._get_east_index(x1, x2)
        self.ychange = self._get_north_index(y1, y2)

        #reset values of resistivity
        self.change_model_res(self.xchange, self.ychange)


    def _get_east_index(self, x1, x2):
        """
        get the index value of the points to be changed

        """
        if x1 < x2:
            xchange = np.where((self.grid_east/self.dscale >= x1) & \
                               (self.grid_east/self.dscale <= x2))[0]
            if len(xchange) == 0:
                xchange = np.where(self.grid_east/self.dscale >= x1)[0][0]-1
                return [xchange]

        if x1 > x2:
            xchange = np.where((self.grid_east/self.dscale <= x1) & \
                               (self.grid_east/self.dscale >= x2))[0]
            if len(xchange) == 0:
                xchange = np.where(self.grid_east/self.dscale >= x2)[0][0]-1
                return [xchange]


        #check the edges to see if the selection should include the square
        xchange = np.append(xchange, xchange[0]-1)
        xchange.sort()

        return xchange

    def _get_north_index(self, y1, y2):
        """
        get the index value of the points to be changed in north direction

        need to flip the index because the plot is flipped

        """

        if y1 < y2:
            ychange = np.where((self.grid_north/self.dscale > y1) & \
                               (self.grid_north/self.dscale < y2))[0]
            if len(ychange) == 0:
                ychange = np.where(self.grid_north/self.dscale >= y1)[0][0]-1
                return [ychange]

        elif y1 > y2:
            ychange = np.where((self.grid_north/self.dscale < y1) & \
                               (self.grid_north/self.dscale > y2))[0]
            if len(ychange) == 0:
                ychange = np.where(self.grid_north/self.dscale >= y2)[0][0]-1
                return [ychange]

        ychange -= 1
        ychange = np.append(ychange, ychange[-1]+1)

        return ychange


    def convert_model_to_int(self):
        """
        convert the resistivity model that is in ohm-m to integer values
        corresponding to res_list

        """

        self.res_model_int = np.ones_like(self.res_model)

        for ii, res in enumerate(self.res_list):
            indexes = np.where(self.res_model == res)
            self.res_model_int[indexes] = self.res_dict[res]
            if ii == 0:
                indexes = np.where(self.res_model <= res)
                self.res_model_int[indexes] = self.res_dict[res]
            elif ii == len(self.res_list)-1:
                indexes = np.where(self.res_model >= res)
                self.res_model_int[indexes] = self.res_dict[res]
            else:
                l_index = max([0, ii-1])
                h_index = min([len(self.res_list)-1, ii+1])
                indexes = np.where((self.res_model > self.res_list[l_index]) &
                                   (self.res_model < self.res_list[h_index]))
                self.res_model_int[indexes] = self.res_dict[res]

    def convert_res_to_model(self, res_array):
        """
        converts an output model into an array of segmented valued according
        to res_list. 

        output is an array of segemented resistivity values in ohm-m (linear)        

        """

        #make values in model resistivity array a value in res_list
        self.res_model = np.zeros_like(res_array)

        for ii, res in enumerate(self.res_list):
            indexes = np.where(res_array == res)
            self.res_model[indexes] = res
            if ii == 0:
                indexes = np.where(res_array <= res)
                self.res_model[indexes] = res
            elif ii == len(self.res_list)-1:
                indexes = np.where(res_array >= res)
                self.res_model[indexes] = res
            else:
                l_index = max([0, ii-1])
                h_index = min([len(self.res_list)-1, ii+1])
                indexes = np.where((res_array > self.res_list[l_index]) &
                                   (res_array < self.res_list[h_index]))
                self.res_model[indexes] = res

    def rewrite_initial_file(self, save_path=None):
        """
        write an initial file for wsinv3d from the model created.
        """
        #need to flip the resistivity model so that the first index is the
        #northern most block in N-S
        #self.res_model = self.res_model[::-1, :, :]

        if save_path is not None:
            self.save_path = save_path

        self.new_initial_fn = os.path.join(self.save_path,
                                           self.initial_fn_basename)
        wsmesh = WSMesh()
        #pass attribute to wsmesh
        att_names = ['nodes_north', 'nodes_east', 'nodes_z', 'grid_east',
                     'grid_north', 'grid_z', 'res_model', 'res_list',
                     'res_dict', 'res_model' ]
        for name in att_names:
            if hasattr(self, name):
                value = getattr(self, name)
                setattr(wsmesh, name, value)

        wsmesh.write_initial_file(initial_fn=self.new_initial_fn)


def cmap_discretize(cmap, N):
    """Return a discrete colormap from the continuous colormap cmap.

         cmap: colormap instance, eg. cm.jet. 
         N: number of colors.

     Example
         x = resize(arange(100), (5,100))
         djet = cmap_discretize(cm.jet, 5)
         imshow(x, cmap=djet)
    """

    colors_i = np.concatenate((np.linspace(0, 1., N), (0.,0.,0.,0.)))
    colors_rgba = cmap(colors_i)
    indices = np.linspace(0, 1., N+1)
    cdict = {}
    for ki,key in enumerate(('red','green','blue')):
        cdict[key] = [(indices[i], colors_rgba[i-1,ki], colors_rgba[i,ki])
                       for i in range(N+1)]
    # Return colormap object.
    return colors.LinearSegmentedColormap(cmap.name + "_%d"%N, cdict, 1024)
#==============================================================================
# response
#==============================================================================
class WSResponse(object):
    """
    class to deal with .resp file output by ws3dinv

    ====================== ====================================================
    Attributes             Description    
    ====================== ====================================================
    n_z                    number of vertical layers             
    period_list            list of periods inverted for
    resp                   np.ndarray structured with keys:
                              * *station* --> station name
                              * *east*    --> relative eastern location in
                                               grid
                              * *north*   --> relative northern location in 
                                               grid
                              * *z_resp*  --> impedance tensor array 
                                               of response with shape
                                         (n_stations, n_freq, 4, dtype=complex)
                              * *z_resp_err--> response impedance tensor error 
    resp_fn                full path to response file
    station_east           location of stations in east direction
    station_fn             full path to station file written by WSStation
    station_names          names of stations
    station_north          location of stations in north direction 
    units                  [ 'mv' | 'other' ] units of impedance tensor
    wl_sites_fn            full path to .sites file from Winglink
    z_resp                 impedance tensors of response with shape
                           (n_stations, n_periods, 2, 2)
    z_resp_err             impedance tensors errors of response with shape
                           (n_stations, n_periods, 2, 2) (zeros)
    ====================== ====================================================

    ====================== ====================================================
    Methods                Description     
    ====================== ====================================================
    read_resp_file         read response file and fill attributes 
    ====================== ====================================================
    """
    
    def __init__(self, resp_fn=None, station_fn=None, wl_station_fn=None):
        self.resp_fn = resp_fn
        self.station_fn = station_fn
        self.wl_sites_fn = wl_station_fn

        self.period_list = None
        self.resp = None

        self.n_z = None
        self.station_east = None
        self.station_north = None
        self.station_name = None
        self.z_resp = None
        self.z_resp_err = None

        self.units = 'mv'
        self._zconv = 796.
        
        if self.resp_fn is not None:
            self.read_resp_file()


    def read_resp_file(self, resp_fn=None, wl_sites_fn=None, station_fn=None):
        """
        read in data file

        Arguments:
        -----------
            **resp_fn** : string
                          full path to data file
            **sites_fn** : string
                           full path to sites file output by winglink.  This is
                           to match the station name with station number.
            **station_fn** : string
                             full path to station location file

        Outputs:
        --------
            **resp** : structure np.ndarray
                      fills the attribute WSData.data with values

            **period_list** : np.ndarray()
                             fills the period list with values.
        """

        if resp_fn is not None:
            self.resp_fn = resp_fn

        if wl_sites_fn is not None:
            self.wl_sites_fn = wl_sites_fn
        if station_fn is not None:
            self.station_fn = station_fn

        if not os.path.isfile(self.resp_fn):
            raise WSInputError('Cannot find {0}, check path'.format(self.resp_fn))

        dfid = file(self.resp_fn, 'r')
        dlines = dfid.readlines()

        #get size number of stations, number of frequencies,
        # number of Z components
        n_stations, n_periods, nz = np.array(dlines[0].strip().split(),
                                             dtype='int')
        nsstart = 2

        self.n_z = nz
        #make a structured array to keep things in for convenience
        z_shape = (n_periods, 2, 2)
        resp_dtype = [('station', '|S10'),
                      ('east', np.float),
                      ('north', np.float),
                      ('z_resp', (np.complex, z_shape)),
                      ('z_resp_err', (np.complex, z_shape))]
        self.resp = np.zeros(n_stations, dtype=resp_dtype)

        findlist = []
        for ii, dline in enumerate(dlines[1:50], 1):
            if dline.find('Station_Location: N-S') == 0:
                findlist.append(ii)
            elif dline.find('Station_Location: E-W') == 0:
                findlist.append(ii)
            elif dline.find('DATA_Period:') == 0:
                findlist.append(ii)

        ncol = len(dlines[nsstart].strip().split())

        #get site names if entered a sites file
        if self.wl_sites_fn != None:
            slist, station_list = wl.read_sites_file(self.wl_sites_fn)
            self.resp['station'] = station_list

        elif self.station_fn != None:
            stations = WSStation(self.station_fn)
            stations.read_station_file()
            self.resp['station'] = stations.names
        else:
            self.resp['station'] = np.arange(n_stations)


        #get N-S locations
        for ii, dline in enumerate(dlines[findlist[0]+1:findlist[1]],0):
            dline = dline.strip().split()
            for jj in range(ncol):
                try:
                    self.resp['north'][ii*ncol+jj] = float(dline[jj])
                except IndexError:
                    pass
                except ValueError:
                    break

        #get E-W locations
        for ii, dline in enumerate(dlines[findlist[1]+1:findlist[2]],0):
            dline = dline.strip().split()
            for jj in range(self.n_z):
                try:
                    self.resp['east'][ii*ncol+jj] = float(dline[jj])
                except IndexError:
                    pass
                except ValueError:
                    break
        #make some empty array to put stuff into
        self.period_list = np.zeros(n_periods)

        #get resp
        per = 0
        for ii, dl in enumerate(dlines[findlist[2]:]):
            if dl.lower().find('period') > 0:
                st = 0
                if dl.lower().find('data') == 0:
                    dkey = 'z_resp'
                    self.period_list[per] = float(dl.strip().split()[1])
                per += 1

            elif dl.lower().find('#iteration') >= 0:
                break
            else:
                zline = np.array(dl.strip().split(),dtype=np.float)*self._zconv
                self.resp[st][dkey][per-1,:] = np.array([[zline[0]-1j*zline[1],
                                                         zline[2]-1j*zline[3]],
                                                         [zline[4]-1j*zline[5],
                                                         zline[6]-1j*zline[7]]])
                st += 1

        self.station_east = self.resp['east']
        self.station_north = self.resp['north']
        self.station_name = self.resp['station']
        self.z_resp = self.resp['z_resp']
        self.z_resp_err = np.zeros_like(self.z_resp)

#==============================================================================
# WSError
#==============================================================================
class WSInputError(Exception):
    pass

#==============================================================================
# plot response
#==============================================================================
class PlotResponse(object):
    """
    plot data and response

    :Example: ::

        >>> import mtpy.modeling.ws3dinv as ws
        >>> dfn = r"/home/MT/ws3dinv/Inv1/WSDataFile.dat"
        >>> rfn = r"/home/MT/ws3dinv/Inv1/Test_resp.00"
        >>> sfn = r"/home/MT/ws3dinv/Inv1/WSStationLocations.txt"
        >>> wsrp = ws.PlotResponse(data_fn=dfn, resp_fn=rfn, station_fn=sfn)
        >>> # plot only the TE and TM modes
        >>> wsrp.plot_component = 2
        >>> wsrp.redraw_plot()

    ======================== ==================================================
    Attributes               Description    
    ======================== ==================================================
    color_mode               [ 'color' | 'bw' ] color or black and white plots
    cted                     color for data TE mode
    ctem                     color for data TM mode
    ctmd                     color for model TE mode
    ctmm                     color for model TM mode
    data_fn                  full path to data file
    data_object              WSResponse instance
    e_capsize                cap size of error bars in points (*default* is .5)
    e_capthick               cap thickness of error bars in points (*default*
                             is 1)
    fig_dpi                  resolution of figure in dots-per-inch (300)
    fig_list                 list of matplotlib.figure instances for plots
    fig_size                 size of figure in inches (*default* is [6, 6])
    font_size                size of font for tick labels, axes labels are
                             font_size+2 (*default* is 7)
    legend_border_axes_pad   padding between legend box and axes 
    legend_border_pad        padding between border of legend and symbols
    legend_handle_text_pad   padding between text labels and symbols of legend
    legend_label_spacing     padding between labels
    legend_loc               location of legend 
    legend_marker_scale      scale of symbols in legend
    lw                       line width response curves (*default* is .5)
    ms                       size of markers (*default* is 1.5)
    mted                     marker for data TE mode
    mtem                     marker for data TM mode
    mtmd                     marker for model TE mode
    mtmm                     marker for model TM mode 
    phase_limits             limits of phase
    plot_component           [ 2 | 4 ] 2 for TE and TM or 4 for all components
    plot_style               [ 1 | 2 ] 1 to plot each mode in a seperate
                             subplot and 2 to plot xx, xy and yx, yy in same 
                             plots
    plot_type                [ '1' | list of station name ] '1' to plot all 
                             stations in data file or input a list of station
                             names to plot if station_fn is input, otherwise
                             input a list of integers associated with the 
                             index with in the data file, ie 2 for 2nd station
    plot_z                   [ True | False ] *default* is True to plot 
                             impedance, False for plotting resistivity and 
                             phase
    plot_yn                  [ 'n' | 'y' ] to plot on instantiation
    res_limits               limits of resistivity in linear scale
    resp_fn                  full path to response file
    resp_object              WSResponse object for resp_fn, or list of 
                             WSResponse objects if resp_fn is a list of
                             response files
    station_fn               full path to station file written by WSStation
    subplot_bottom           space between axes and bottom of figure
    subplot_hspace           space between subplots in vertical direction
    subplot_left             space between axes and left of figure
    subplot_right            space between axes and right of figure
    subplot_top              space between axes and top of figure
    subplot_wspace           space between subplots in horizontal direction    
    ======================== ==================================================
    """

    def __init__(self, data_fn=None, resp_fn=None, station_fn=None, **kwargs):
        self.data_fn = data_fn
        self.resp_fn = resp_fn
        self.station_fn = station_fn

        self.data_object = None
        self.resp_object = []

        self.color_mode = kwargs.pop('color_mode', 'color')

        self.ms = kwargs.pop('ms', 1.5)
        self.lw = kwargs.pop('lw', .5)
        self.e_capthick = kwargs.pop('e_capthick', .5)
        self.e_capsize = kwargs.pop('e_capsize', 2)

        #color mode
        if self.color_mode == 'color':
            #color for data
            self.cted = kwargs.pop('cted', (0, 0, 1))
            self.ctmd = kwargs.pop('ctmd', (1, 0, 0))
            self.mted = kwargs.pop('mted', 's')
            self.mtmd = kwargs.pop('mtmd', 'o')

            #color for occam2d model
            self.ctem = kwargs.pop('ctem', (0, .6, .3))
            self.ctmm = kwargs.pop('ctmm', (.9, 0, .8))
            self.mtem = kwargs.pop('mtem', '+')
            self.mtmm = kwargs.pop('mtmm', '+')

        #black and white mode
        elif self.color_mode == 'bw':
            #color for data
            self.cted = kwargs.pop('cted', (0, 0, 0))
            self.ctmd = kwargs.pop('ctmd', (0, 0, 0))
            self.mted = kwargs.pop('mted', '*')
            self.mtmd = kwargs.pop('mtmd', 'v')

            #color for occam2d model
            self.ctem = kwargs.pop('ctem', (0.6, 0.6, 0.6))
            self.ctmm = kwargs.pop('ctmm', (0.6, 0.6, 0.6))
            self.mtem = kwargs.pop('mtem', '+')
            self.mtmm = kwargs.pop('mtmm', 'x')

        self.phase_limits = kwargs.pop('phase_limits', None)
        self.res_limits = kwargs.pop('res_limits', None)

        self.fig_num = kwargs.pop('fig_num', 1)
        self.fig_size = kwargs.pop('fig_size', [6, 6])
        self.fig_dpi = kwargs.pop('dpi', 300)

        self.subplot_wspace = .2
        self.subplot_hspace = .0
        self.subplot_right = .98
        self.subplot_left = .08
        self.subplot_top = .93
        self.subplot_bottom = .1

        self.legend_loc = 'upper left'
        self.legend_marker_scale = 1
        self.legend_border_axes_pad = .01
        self.legend_label_spacing = 0.07
        self.legend_handle_text_pad = .2
        self.legend_border_pad = .15

        self.font_size = kwargs.pop('font_size', 6)

        self.plot_type = kwargs.pop('plot_type', '1')
        self.plot_style = kwargs.pop('plot_style', 1)
        self.plot_component = kwargs.pop('plot_component', 4)
        self.plot_yn = kwargs.pop('plot_yn', 'y')
        self.plot_z = kwargs.pop('plot_z', True)
        self.ylabel_pad = kwargs.pop('ylabel_pad', 1.25)

        self.fig_list = []

        if self.plot_yn == 'y':
            self.plot()

    def plot_errorbar(self, ax, period, data, error, color, marker):
        """
        convinience function to make an error bar instance
        """

        errorbar_object = ax.errorbar(period,
                                      data,
                                      marker=marker,
                                      ms=self.ms,
                                      mfc='None',
                                      mec=color,
                                      ls=':',
                                      yerr=error,
                                      ecolor=color,
                                      color=color,
                                      picker=2,
                                      lw=self.lw,
                                      elinewidth=self.lw,
                                      capsize=self.e_capsize,
                                      capthick=self.e_capthick)
        return errorbar_object

    def plot(self):
        """
        plot
        """

        self.data_object = WSData()
        self.data_object.read_data_file(self.data_fn,
                                        station_fn=self.station_fn)

        #get shape of impedance tensors
        ns = self.data_object.data['station'].shape[0]
        nf = len(self.data_object.period_list)

        #read in response files
        if self.resp_fn != None:
            self.resp_object = []
            if type(self.resp_fn) is not list:
                self.resp_object = [WSResponse(self.resp_fn,
                                               station_fn=self.station_fn)]
            else:
                for rfile in self.resp_fn:
                    self.resp_object.append(WSResponse(rfile,
                                                   station_fn=self.station_fn))

        #get number of response files
        nr = len(self.resp_object)

        if type(self.plot_type) is list:
            ns = len(self.plot_type)

        #--> set default font size
        plt.rcParams['font.size'] = self.font_size

        fontdict = {'size':self.font_size+2, 'weight':'bold'}
        if self.plot_z == True:
            h_ratio = [1,1]
        elif self.plot_z == False:
            h_ratio = [2, 1.5]
        gs = gridspec.GridSpec(2, 2, height_ratios=h_ratio, hspace=.1)

        ax_list = []
        line_list = []
        label_list = []


        if self.plot_type != '1':
            pstation_list = []
            if type(self.plot_type) is not list:
                self.plot_type = [self.plot_type]
            for ii, station in enumerate(self.data_object.data['station']):
                if type(station) is not int:
                    for pstation in self.plot_type:
                        if station.find(str(pstation)) >= 0:
                            pstation_list.append(ii)
                else:
                    for pstation in self.plot_type:
                        if station == int(pstation):
                            pstation_list.append(ii)
        else:
            pstation_list = np.arange(ns)

        for jj in pstation_list:
            data_z = self.data_object.z_data[jj]
            data_z_err = self.data_object.z_data_err[jj]
            period = self.data_object.period_list
            station = self.data_object.station_names[jj]
            print('Plotting: {0}'.format(station))

            #check for masked points
            data_z[np.where(data_z == 7.95204E5-7.95204E5j)] = 0.0+0.0j
            data_z_err[np.where(data_z_err == 7.95204E5-7.95204E5j)] =\
                                                                1.0+1.0j

            #convert to apparent resistivity and phase
            z_object =  mtz.Z(z_array=data_z, z_err_array=data_z_err,
                              freq=1./period)

            rp = mtplottools.ResPhase(z_object)

            #find locations where points have been masked
            nzxx = np.where(rp.resxx!=0)[0]
            nzxy = np.where(rp.resxy!=0)[0]
            nzyx = np.where(rp.resyx!=0)[0]
            nzyy = np.where(rp.resyy!=0)[0]

            if self.resp_fn != None:
                plotr = True
            else:
                plotr = False

            #make figure
            fig = plt.figure(station, self.fig_size, dpi=self.fig_dpi)
            plt.clf()
            fig.suptitle(str(station), fontdict=fontdict)

            #set the grid of subplots
            gs = gridspec.GridSpec(2, 4,
                                   wspace=self.subplot_wspace,
                                   left=self.subplot_left,
                                   top=self.subplot_top,
                                   bottom=self.subplot_bottom,
                                   right=self.subplot_right,
                                   hspace=self.subplot_hspace,
                                   height_ratios=h_ratio)
            #---------plot the apparent resistivity-----------------------------------
            #plot each component in its own subplot
            if self.plot_style == 1:
                if self.plot_component == 2:
                    axrxy = fig.add_subplot(gs[0, 0:2])
                    axryx = fig.add_subplot(gs[0, 2:], sharex=axrxy)

                    axpxy = fig.add_subplot(gs[1, 0:2])
                    axpyx = fig.add_subplot(gs[1, 2:], sharex=axrxy)

                    if self.plot_z == False:
                        #plot resistivity
                        erxy = self.plot_errorbar(axrxy,
                                                  period[nzxy],
                                                  rp.resxy[nzxy],
                                                  rp.resxy_err[nzxy],
                                                  self.cted, self.mted)
                        eryx = self.plot_errorbar(axryx,
                                                  period[nzyx],
                                                  rp.resyx[nzyx],
                                                  rp.resyx_err[nzyx],
                                                  self.ctmd, self.mtmd)
                        #plot phase
                        erxy = self.plot_errorbar(axpxy,
                                                  period[nzxy],
                                                  rp.phasexy[nzxy],
                                                  rp.phasexy_err[nzxy],
                                                  self.cted, self.mted)
                        eryx = self.plot_errorbar(axpyx,
                                                  period[nzyx],
                                                  rp.phaseyx[nzyx],
                                                  rp.phaseyx_err[nzyx],
                                                  self.ctmd, self.mtmd)
                    elif self.plot_z == True:
                        #plot real
                        erxy = self.plot_errorbar(axrxy,
                                                  period[nzxy],
                                                  z_object.z[nzxy,0,1].real,
                                                  z_object.z_err[nzxy,0,1].real,
                                                  self.cted, self.mted)
                        eryx = self.plot_errorbar(axryx,
                                                  period[nzyx],
                                                  z_object.z[nzyx,1,0].real,
                                                  z_object.z_err[nzyx,1,0].real,
                                                  self.ctmd, self.mtmd)
                        #plot phase
                        erxy = self.plot_errorbar(axpxy,
                                                  period[nzxy],
                                                  z_object.z[nzxy,0,1].imag,
                                                  z_object.z_err[nzxy,0,1].imag,
                                                  self.cted, self.mted)
                        eryx = self.plot_errorbar(axpyx,
                                                  period[nzyx],
                                                  z_object.z[nzyx,1,0].imag,
                                                  z_object.z_err[nzyx,1,0].imag,
                                                  self.ctmd, self.mtmd)

                    ax_list = [axrxy, axryx, axpxy, axpyx]
                    line_list = [[erxy[0]], [eryx[0]]]
                    label_list = [['$Z_{xy}$'], ['$Z_{yx}$']]

                elif self.plot_component == 4:
                    axrxx = fig.add_subplot(gs[0, 0])
                    axrxy = fig.add_subplot(gs[0, 1], sharex=axrxx)
                    axryx = fig.add_subplot(gs[0, 2], sharex=axrxx)
                    axryy = fig.add_subplot(gs[0, 3], sharex=axrxx)

                    axpxx = fig.add_subplot(gs[1, 0])
                    axpxy = fig.add_subplot(gs[1, 1], sharex=axrxx)
                    axpyx = fig.add_subplot(gs[1, 2], sharex=axrxx)
                    axpyy = fig.add_subplot(gs[1, 3], sharex=axrxx)

                    if self.plot_z == False:
                        #plot resistivity
                        erxx= self.plot_errorbar(axrxx,
                                                  period[nzxx],
                                                  rp.resxx[nzxx],
                                                  rp.resxx_err[nzxx],
                                                  self.cted, self.mted)
                        erxy = self.plot_errorbar(axrxy,
                                                  period[nzxy],
                                                  rp.resxy[nzxy],
                                                  rp.resxy_err[nzxy],
                                                  self.cted, self.mted)
                        eryx = self.plot_errorbar(axryx,
                                                  period[nzyx],
                                                  rp.resyx[nzyx],
                                                  rp.resyx_err[nzyx],
                                                  self.ctmd, self.mtmd)
                        eryy = self.plot_errorbar(axryy,
                                                  period[nzyy],
                                                  rp.resyy[nzyy],
                                                  rp.resyy_err[nzyy],
                                                  self.ctmd, self.mtmd)
                        #plot phase
                        erxx= self.plot_errorbar(axpxx,
                                                  period[nzxx],
                                                  rp.phasexx[nzxx],
                                                  rp.phasexx_err[nzxx],
                                                  self.cted, self.mted)
                        erxy = self.plot_errorbar(axpxy,
                                                  period[nzxy],
                                                  rp.phasexy[nzxy],
                                                  rp.phasexy_err[nzxy],
                                                  self.cted, self.mted)
                        eryx = self.plot_errorbar(axpyx,
                                                  period[nzyx],
                                                  rp.phaseyx[nzyx],
                                                  rp.phaseyx_err[nzyx],
                                                  self.ctmd, self.mtmd)
                        eryy = self.plot_errorbar(axpyy,
                                                  period[nzyy],
                                                  rp.phaseyy[nzyy],
                                                  rp.phaseyy_err[nzyy],
                                                  self.ctmd, self.mtmd)
                    elif self.plot_z == True:
                        #plot real
                        erxx = self.plot_errorbar(axrxx,
                                                  period[nzxx],
                                                  z_object.z[nzxx,0,0].real,
                                                  z_object.z_err[nzxx,0,0].real,
                                                  self.cted, self.mted)
                        erxy = self.plot_errorbar(axrxy,
                                                  period[nzxy],
                                                  z_object.z[nzxy,0,1].real,
                                                  z_object.z_err[nzxy,0,1].real,
                                                  self.cted, self.mted)
                        eryx = self.plot_errorbar(axryx,
                                                  period[nzyx],
                                                  z_object.z[nzyx,1,0].real,
                                                  z_object.z_err[nzyx,1,0].real,
                                                  self.ctmd, self.mtmd)
                        eryy = self.plot_errorbar(axryy,
                                                  period[nzyy],
                                                  z_object.z[nzyy,1,1].real,
                                                  z_object.z_err[nzyy,1,1].real,
                                                  self.ctmd, self.mtmd)
                        #plot phase
                        erxx = self.plot_errorbar(axpxx,
                                                  period[nzxx],
                                                  z_object.z[nzxx,0,0].imag,
                                                  z_object.z_err[nzxx,0,0].imag,
                                                  self.cted, self.mted)
                        erxy = self.plot_errorbar(axpxy,
                                                  period[nzxy],
                                                  z_object.z[nzxy,0,1].imag,
                                                  z_object.z_err[nzxy,0,1].imag,
                                                  self.cted, self.mted)
                        eryx = self.plot_errorbar(axpyx,
                                                  period[nzyx],
                                                  z_object.z[nzyx,1,0].imag,
                                                  z_object.z_err[nzyx,1,0].imag,
                                                  self.ctmd, self.mtmd)
                        eryy = self.plot_errorbar(axpyy,
                                                  period[nzyy],
                                                  z_object.z[nzyy,1,1].imag,
                                                  z_object.z_err[nzyy,1,1].imag,
                                                  self.ctmd, self.mtmd)

                    ax_list = [axrxx, axrxy, axryx, axryy,
                               axpxx, axpxy, axpyx, axpyy]
                    line_list = [[erxx[0]], [erxy[0]], [eryx[0]], [eryy[0]]]
                    label_list = [['$Z_{xx}$'], ['$Z_{xy}$'],
                                  ['$Z_{yx}$'], ['$Z_{yy}$']]

                #set axis properties
                for aa, ax in enumerate(ax_list):
                    ax.tick_params(axis='y', pad=self.ylabel_pad)
                    if len(ax_list) == 4:
                        if aa < 2:
                            plt.setp(ax.get_xticklabels(), visible=False)
                            if self.plot_z == False:
                                ax.set_yscale('log', nonposy='clip')
                            if self.res_limits is not None:
                                ax.set_ylim(self.res_limits)
                        else:
                            ax.set_ylim(self.phase_limits)
                            ax.set_xlabel('Period (s)', fontdict=fontdict)
                        #set axes labels
                        if aa == 0:
                            if self.plot_z == False:
                                ax.set_ylabel('App. Res. ($\mathbf{\Omega \cdot m}$)',
                                              fontdict=fontdict)
                            elif self.plot_z == True:
                                ax.set_ylabel('Re[Z (m/s)]',
                                              fontdict=fontdict)
                        elif aa == 2:
                            if self.plot_z == False:
                                ax.set_ylabel('Phase (deg)',
                                              fontdict=fontdict)
                            elif self.plot_z == True:
                                ax.set_ylabel('Im[Z (m/s)]',
                                              fontdict=fontdict)
#                        else:
#                            plt.setp(ax.yaxis.get_ticklabels(), visible=False)

                    elif len(ax_list) == 8:
                        if aa < 4:
                            plt.setp(ax.get_xticklabels(), visible=False)
                            if self.plot_z == False:
                                ax.set_yscale('log', nonposy='clip')
                            if self.res_limits is not None:
                                ax.set_ylim(self.res_limits)
                        else:
                            ax.set_ylim(self.phase_limits)
                            ax.set_xlabel('Period (s)', fontdict=fontdict)
                        #set axes labels
                        if aa == 0:
                            if self.plot_z == False:
                                ax.set_ylabel('App. Res. ($\mathbf{\Omega \cdot m}$)',
                                              fontdict=fontdict)
                            elif self.plot_z == True:
                                ax.set_ylabel('Re[Z (m/s)]',
                                               fontdict=fontdict)
                        elif aa == 4:
                            if self.plot_z == False:
                                ax.set_ylabel('Phase (deg)',
                                              fontdict=fontdict)
                            elif self.plot_z == True:
                                ax.set_ylabel('Im[Z (m/s)]',
                                              fontdict=fontdict)
#                        else:
#                            plt.setp(ax.yaxis.get_ticklabels(), visible=False)

                    ax.set_xscale('log', nonposx='clip')
                    ax.set_xlim(xmin=10**(np.floor(np.log10(period[0]))) * 1.01,
                                xmax=10**(np.ceil(np.log10(period[-1]))) * .99)
                    ax.grid(True, alpha=.25)

            # plot xy and yx together and xx, yy together
            elif self.plot_style == 2:
                if self.plot_component == 2:
                    axrxy = fig.add_subplot(gs[0, 0:])
                    axpxy = fig.add_subplot(gs[1, 0:], sharex=axrxy)
                    if self.plot_z == False:
                        #plot resistivity
                        erxy = self.plot_errorbar(axrxy,
                                                  period[nzxy],
                                                  rp.resxy[nzxy],
                                                  rp.resxy_err[nzxy],
                                                  self.cted, self.mted)
                        eryx = self.plot_errorbar(axrxy,
                                                  period[nzyx],
                                                  rp.resyx[nzyx],
                                                  rp.resyx_err[nzyx],
                                                  self.ctmd, self.mtmd)
                        #plot phase
                        erxy = self.plot_errorbar(axpxy,
                                                  period[nzxy],
                                                  rp.phasexy[nzxy],
                                                  rp.phasexy_err[nzxy],
                                                  self.cted, self.mted)
                        eryx = self.plot_errorbar(axpxy,
                                                  period[nzyx],
                                                  rp.phaseyx[nzyx],
                                                  rp.phaseyx_err[nzyx],
                                                  self.ctmd, self.mtmd)
                    elif self.plot_z == True:
                        #plot real
                        erxy = self.plot_errorbar(axrxy,
                                                  period[nzxy],
                                                  z_object.z[nzxy,0,1].real,
                                                  z_object.z_err[nzxy,0,1].real,
                                                  self.cted, self.mted)
                        eryx = self.plot_errorbar(axrxy,
                                                  period[nzxy],
                                                  z_object.z[nzxy,1,0].real,
                                                  z_object.z_err[nzxy,1,0].real,
                                                  self.ctmd, self.mtmd)
                        #plot phase
                        erxy = self.plot_errorbar(axpxy,
                                                  period[nzxy],
                                                  z_object.z[nzxy,0,1].imag,
                                                  z_object.z_err[nzxy,0,1].imag,
                                                  self.cted, self.mted)
                        eryx = self.plot_errorbar(axpxy,
                                                  period[nzyx],
                                                  z_object.z[nzyx,1,0].imag,
                                                  z_object.z_err[nzyx,1,0].imag,
                                                  self.ctmd, self.mtmd)

                    ax_list = [axrxy, axpxy]
                    line_list = [erxy[0], eryx[0]]
                    label_list = ['$Z_{xy}$', '$Z_{yx}$']

                elif self.plot_component == 4:
                    axrxy = fig.add_subplot(gs[0, 0:2])
                    axpxy = fig.add_subplot(gs[1, 0:2], sharex=axrxy)

                    axrxx = fig.add_subplot(gs[0, 2:], sharex=axrxy)
                    axpxx = fig.add_subplot(gs[1, 2:], sharex=axrxy)
                    if self.plot_z == False:
                        #plot resistivity
                        erxx= self.plot_errorbar(axrxx,
                                                  period[nzxx],
                                                  rp.resxx[nzxx],
                                                  rp.resxx_err[nzxx],
                                                  self.cted, self.mted)
                        erxy = self.plot_errorbar(axrxy,
                                                  period[nzxy],
                                                  rp.resxy[nzxy],
                                                  rp.resxy_err[nzxy],
                                                  self.cted, self.mted)
                        eryx = self.plot_errorbar(axrxy,
                                                  period[nzyx],
                                                  rp.resyx[nzyx],
                                                  rp.resyx_err[nzyx],
                                                  self.ctmd, self.mtmd)
                        eryy = self.plot_errorbar(axrxx,
                                                  period[nzyy],
                                                  rp.resyy[nzyy],
                                                  rp.resyy_err[nzyy],
                                                  self.ctmd, self.mtmd)
                        #plot phase
                        erxx= self.plot_errorbar(axpxx,
                                                  period[nzxx],
                                                  rp.phasexx[nzxx],
                                                  rp.phasexx_err[nzxx],
                                                  self.cted, self.mted)
                        erxy = self.plot_errorbar(axpxy,
                                                  period[nzxy],
                                                  rp.phasexy[nzxy],
                                                  rp.phasexy_err[nzxy],
                                                  self.cted, self.mted)
                        eryx = self.plot_errorbar(axpxy,
                                                  period[nzyx],
                                                  rp.phaseyx[nzyx],
                                                  rp.phaseyx_err[nzyx],
                                                  self.ctmd, self.mtmd)
                        eryy = self.plot_errorbar(axpxx,
                                                  period[nzyy],
                                                  rp.phaseyy[nzyy],
                                                  rp.phaseyy_err[nzyy],
                                                  self.ctmd, self.mtmd)
                    elif self.plot_z == True:
                         #plot real
                        erxx = self.plot_errorbar(axrxx,
                                                  period[nzxx],
                                                  z_object.z[nzxx,0,0].real,
                                                  z_object.z_err[nzxx,0,0].real,
                                                  self.cted, self.mted)
                        erxy = self.plot_errorbar(axrxy,
                                                  period[nzxy],
                                                  z_object.z[nzxy,0,1].real,
                                                  z_object.z_err[nzxy,0,1].real,
                                                  self.cted, self.mted)
                        eryx = self.plot_errorbar(axrxy,
                                                  period[nzyx],
                                                  z_object.z[nzyx,1,0].real,
                                                  z_object.z_err[nzyx,1,0].real,
                                                  self.ctmd, self.mtmd)
                        eryy = self.plot_errorbar(axrxx,
                                                  period[nzyy],
                                                  z_object.z[nzyy,1,1].real,
                                                  z_object.z_err[nzyy,1,1].real,
                                                  self.ctmd, self.mtmd)
                        #plot phase
                        erxx = self.plot_errorbar(axpxx,
                                                  period[nzxx],
                                                  z_object.z[nzxx,0,0].imag,
                                                  z_object.z_err[nzxx,0,0].imag,
                                                  self.cted, self.mted)
                        erxy = self.plot_errorbar(axpxy,
                                                  period[nzxy],
                                                  z_object.z[nzxy,0,1].imag,
                                                  z_object.z_err[nzxy,0,1].imag,
                                                  self.cted, self.mted)
                        eryx = self.plot_errorbar(axpxy,
                                                  period[nzyx],
                                                  z_object.z[nzyx,1,0].imag,
                                                  z_object.z_err[nzyx,1,0].imag,
                                                  self.ctmd, self.mtmd)
                        eryy = self.plot_errorbar(axpxx,
                                                  period[nzyy],
                                                  z_object.z[nzyy,1,1].imag,
                                                  z_object.z_err[nzyy,1,1].imag,
                                                  self.ctmd, self.mtmd)

                    ax_list = [axrxy, axrxx, axpxy, axpxx]
                    line_list = [[erxy[0], eryx[0]], [erxx[0], eryy[0]]]
                    label_list = [['$Z_{xy}$', '$Z_{yx}$'],
                                  ['$Z_{xx}$', '$Z_{yy}$']]
                #set axis properties
                for aa, ax in enumerate(ax_list):
                    ax.tick_params(axis='y', pad=self.ylabel_pad)
                    if len(ax_list) == 2:
                        ax.set_xlabel('Period (s)', fontdict=fontdict)
                        if aa == 0:
                            plt.setp(ax.get_xticklabels(), visible=False)
                            if self.plot_z == False:
                                ax.set_yscale('log', nonposy='clip')
                                ax.set_ylabel('App. Res. ($\mathbf{\Omega \cdot m}$)',
                                              fontdict=fontdict)
                            elif self.plot_z == True:
                                ax.set_ylabel('Re[Impedance (m/s)]',
                                               fontdict=fontdict)
                            if self.res_limits is not None:
                                ax.set_ylim(self.res_limits)
                        else:
                            ax.set_ylim(self.phase_limits)
                            if self.plot_z == False:
                                ax.set_ylabel('Phase (deg)',
                                              fontdict=fontdict)
                            elif self.plot_z == True:
                                ax.set_ylabel('Im[Impedance (m/s)]',
                                              fontdict=fontdict)
                    elif len(ax_list) == 4:
                        if aa < 2:
                            plt.setp(ax.get_xticklabels(), visible=False)
                            if self.plot_z == False:
                                ax.set_yscale('log', nonposy='clip')
                            if self.res_limits is not None:
                                ax.set_ylim(self.res_limits)
                        else:
                            ax.set_ylim(self.phase_limits)
                            ax.set_xlabel('Period (s)', fontdict=fontdict)
                        if aa == 0:
                            if self.plot_z == False:
                                ax.set_ylabel('App. Res. ($\mathbf{\Omega \cdot m}$)',
                                              fontdict=fontdict)
                            elif self.plot_z == True:
                                ax.set_ylabel('Re[Impedance (m/s)]',
                                               fontdict=fontdict)
                        elif aa == 2:
                            if self.plot_z == False:
                                ax.set_ylabel('Phase (deg)',
                                              fontdict=fontdict)
                            elif self.plot_z == True:
                                ax.set_ylabel('Im[Impedance (m/s)]',
                                              fontdict=fontdict)
#                        else:
#                            plt.setp(ax.yaxis.get_ticklabels(), visible=False)

                    ax.set_xscale('log', nonposx='clip')
                    ax.set_xlim(xmin=10**(np.floor(np.log10(period[0]))) * 1.01,
                                xmax=10**(np.ceil(np.log10(period[-1]))) * .99)
                    ax.grid(True, alpha=.25)

            if plotr == True:
                for rr in range(nr):
                    if self.color_mode == 'color':
                        cxy = (0,.4+float(rr)/(3*nr),0)
                        cyx = (.7+float(rr)/(4*nr),.13,.63-float(rr)/(4*nr))
                    elif self.color_mode == 'bw':
                        cxy = (1-1.25/(rr+2.),1-1.25/(rr+2.),1-1.25/(rr+2.))
                        cyx = (1-1.25/(rr+2.),1-1.25/(rr+2.),1-1.25/(rr+2.))

                    resp_z = self.resp_object[rr].z_resp[jj]
                    resp_z_err = (data_z-resp_z)/(data_z_err)
                    resp_z_object =  mtz.Z(z_array=resp_z,
                                           z_err_array=resp_z_err,
                                           freq=1./period)

                    rrp = mtplottools.ResPhase(resp_z_object)

                    rms = resp_z_err.std()
                    rms_xx = resp_z_err[:, 0, 0].std()
                    rms_xy = resp_z_err[:, 0, 1].std()
                    rms_yx = resp_z_err[:, 1, 0].std()
                    rms_yy = resp_z_err[:, 1, 1].std()
                    print(' --- response {0} ---'.format(rr))
                    print('  RMS = {:.2f}'.format(rms))
                    print('      RMS_xx = {:.2f}'.format(rms_xx))
                    print('      RMS_xy = {:.2f}'.format(rms_xy))
                    print('      RMS_yx = {:.2f}'.format(rms_yx))
                    print('      RMS_yy = {:.2f}'.format(rms_yy))

                    if self.plot_style == 1:
                        if self.plot_component == 2:
                            if self.plot_z == False:
                                #plot resistivity
                                rerxy = self.plot_errorbar(axrxy,
                                                          period[nzxy],
                                                          rrp.resxy[nzxy],
                                                          None,
                                                          cxy, self.mted)
                                reryx = self.plot_errorbar(axryx,
                                                          period[nzyx],
                                                          rrp.resyx[nzyx],
                                                          None,
                                                          cyx, self.mtmd)
                                #plot phase
                                rerxy = self.plot_errorbar(axpxy,
                                                          period[nzxy],
                                                          rrp.phasexy[nzxy],
                                                          None,
                                                          cxy, self.mted)
                                reryx = self.plot_errorbar(axpyx,
                                                          period[nzyx],
                                                          rrp.phaseyx[nzyx],
                                                          None,
                                                          cyx, self.mtmd)
                            elif self.plot_z == True:
                                #plot real
                                rerxy = self.plot_errorbar(axrxy,
                                                          period[nzxy],
                                                          resp_z[nzxy,0,1].real,
                                                          None,
                                                          cxy, self.mted)
                                reryx = self.plot_errorbar(axryx,
                                                          period[nzyx],
                                                          resp_z[nzyx,1,0].real,
                                                          None,
                                                          cyx, self.mtmd)
                                #plot phase
                                rerxy = self.plot_errorbar(axpxy,
                                                          period[nzxy],
                                                          resp_z[nzxy,0,1].imag,
                                                          None,
                                                          cxy, self.mted)
                                reryx = self.plot_errorbar(axpyx,
                                                          period[nzyx],
                                                          resp_z[nzyx,1,0].imag,
                                                          None,
                                                          cyx, self.mtmd)

                            line_list[0] += [rerxy[0]]
                            line_list[1] += [reryx[0]]
                            label_list[0] += ['$Z^m_{xy}$ '+
                                               'rms={0:.2f}'.format(rms_xy)]
                            label_list[1] += ['$Z^m_{yx}$ '+
                                           'rms={0:.2f}'.format(rms_yx)]
                        elif self.plot_component == 4:
                            if self.plot_z == False:
                                #plot resistivity
                                rerxx= self.plot_errorbar(axrxx,
                                                          period[nzxx],
                                                          rrp.resxx[nzxx],
                                                          None,
                                                          cxy, self.mted)
                                rerxy = self.plot_errorbar(axrxy,
                                                          period[nzxy],
                                                          rrp.resxy[nzxy],
                                                          None,
                                                          cxy, self.mted)
                                reryx = self.plot_errorbar(axryx,
                                                          period[nzyx],
                                                          rrp.resyx[nzyx],
                                                          None,
                                                          cyx, self.mtmd)
                                reryy = self.plot_errorbar(axryy,
                                                          period[nzyy],
                                                          rrp.resyy[nzyy],
                                                          None,
                                                          cyx, self.mtmd)
                                #plot phase
                                rerxx= self.plot_errorbar(axpxx,
                                                          period[nzxx],
                                                          rrp.phasexx[nzxx],
                                                          None,
                                                          cxy, self.mted)
                                rerxy = self.plot_errorbar(axpxy,
                                                          period[nzxy],
                                                          rrp.phasexy[nzxy],
                                                          None,
                                                          cxy, self.mted)
                                reryx = self.plot_errorbar(axpyx,
                                                          period[nzyx],
                                                          rrp.phaseyx[nzyx],
                                                          None,
                                                          cyx, self.mtmd)
                                reryy = self.plot_errorbar(axpyy,
                                                          period[nzyy],
                                                          rrp.phaseyy[nzyy],
                                                          None,
                                                          cyx, self.mtmd)
                            elif self.plot_z == True:
                                #plot real
                                rerxx = self.plot_errorbar(axrxx,
                                                          period[nzxx],
                                                          resp_z[nzxx,0,0].real,
                                                          None,
                                                          cxy, self.mted)
                                rerxy = self.plot_errorbar(axrxy,
                                                          period[nzxy],
                                                          resp_z[nzxy,0,1].real,
                                                          None,
                                                          cxy, self.mted)
                                reryx = self.plot_errorbar(axryx,
                                                          period[nzyx],
                                                          resp_z[nzyx,1,0].real,
                                                          None,
                                                          cyx, self.mtmd)
                                reryy = self.plot_errorbar(axryy,
                                                          period[nzyy],
                                                          resp_z[nzyy,1,1].real,
                                                          None,
                                                          cyx, self.mtmd)
                                #plot phase
                                rerxx = self.plot_errorbar(axpxx,
                                                          period[nzxx],
                                                          resp_z[nzxx,0,0].imag,
                                                          None,
                                                          cxy, self.mted)
                                rerxy = self.plot_errorbar(axpxy,
                                                          period[nzxy],
                                                          resp_z[nzxy,0,1].imag,
                                                          None,
                                                          cxy, self.mted)
                                reryx = self.plot_errorbar(axpyx,
                                                          period[nzyx],
                                                          resp_z[nzyx,1,0].imag,
                                                          None,
                                                          cyx, self.mtmd)
                                reryy = self.plot_errorbar(axpyy,
                                                          period[nzyy],
                                                          resp_z[nzyy,1,1].imag,
                                                          None,
                                                          cyx, self.mtmd)

                            line_list[0] += [rerxx[0]]
                            line_list[1] += [rerxy[0]]
                            line_list[2] += [reryx[0]]
                            line_list[3] += [reryy[0]]
                            label_list[0] += ['$Z^m_{xx}$ '+
                                               'rms={0:.2f}'.format(rms_xx)]
                            label_list[1] += ['$Z^m_{xy}$ '+
                                           'rms={0:.2f}'.format(rms_xy)]
                            label_list[2] += ['$Z^m_{yx}$ '+
                                           'rms={0:.2f}'.format(rms_yx)]
                            label_list[3] += ['$Z^m_{yy}$ '+
                                           'rms={0:.2f}'.format(rms_yy)]
                    elif self.plot_style == 2:
                        if self.plot_component == 2:
                            if self.plot_z == False:
                                #plot resistivity
                                rerxy = self.plot_errorbar(axrxy,
                                                          period[nzxy],
                                                          rrp.resxy[nzxy],
                                                          None,
                                                          cxy, self.mted)
                                reryx = self.plot_errorbar(axrxy,
                                                          period[nzyx],
                                                          rrp.resyx[nzyx],
                                                          None,
                                                          cyx, self.mtmd)
                                #plot phase
                                rerxy = self.plot_errorbar(axpxy,
                                                          period[nzxy],
                                                          rrp.phasexy[nzxy],
                                                          None,
                                                          cxy, self.mted)
                                reryx = self.plot_errorbar(axpxy,
                                                          period[nzyx],
                                                          rrp.phaseyx[nzyx],
                                                          None,
                                                          cyx, self.mtmd)
                            elif self.plot_z == True:
                                #plot real
                                rerxy = self.plot_errorbar(axrxy,
                                                          period[nzxy],
                                                          resp_z[nzxy,0,1].real,
                                                          None,
                                                          cxy, self.mted)
                                reryx = self.plot_errorbar(axrxy,
                                                          period[nzyx],
                                                          resp_z[nzyx,1,0].real,
                                                          None,
                                                          cyx, self.mtmd)
                                #plot phase
                                rerxy = self.plot_errorbar(axpxy,
                                                          period[nzxy],
                                                          resp_z[nzxy,0,1].imag,
                                                          None,
                                                          cyx, self.mted)
                                reryx = self.plot_errorbar(axpxy,
                                                          period[nzyx],
                                                          resp_z[nzyx,1,0].imag,
                                                          None,
                                                          cyx, self.mtmd)
                            line_list += [rerxy[0], reryx[0]]
                            label_list += ['$Z^m_{xy}$ '+
                                           'rms={0:.2f}'.format(rms_xy),
                                           '$Z^m_{yx}$ '+
                                           'rms={0:.2f}'.format(rms_yx)]
                        elif self.plot_component == 4:
                            if self.plot_z == False:
                                #plot resistivity
                                rerxx= self.plot_errorbar(axrxx,
                                                          period[nzxx],
                                                          rrp.resxx[nzxx],
                                                          None,
                                                          cxy, self.mted)
                                rerxy = self.plot_errorbar(axrxy,
                                                          period[nzxy],
                                                          rrp.resxy[nzxy],
                                                          None,
                                                          cxy, self.mted)
                                reryx = self.plot_errorbar(axrxy,
                                                          period[nzyx],
                                                          rrp.resyx[nzyx],
                                                          None,
                                                          cyx, self.mtmd)
                                reryy = self.plot_errorbar(axrxx,
                                                          period[nzyy],
                                                          rrp.resyy[nzyy],
                                                          None,
                                                          cyx, self.mtmd)
                                #plot phase
                                rerxx= self.plot_errorbar(axpxx,
                                                          period[nzxx],
                                                          rrp.phasexx[nzxx],
                                                          None,
                                                          cxy, self.mted)
                                rerxy = self.plot_errorbar(axpxy,
                                                          period[nzxy],
                                                          rrp.phasexy[nzxy],
                                                          None,
                                                          cxy, self.mted)
                                reryx = self.plot_errorbar(axpxy,
                                                          period[nzyx],
                                                          rrp.phaseyx[nzyx],
                                                          None,
                                                          cyx, self.mtmd)
                                reryy = self.plot_errorbar(axpxx,
                                                          period[nzyy],
                                                          rrp.phaseyy[nzyy],
                                                          None,
                                                          cyx, self.mtmd)
                            elif self.plot_z == True:
                                #plot real
                                rerxx = self.plot_errorbar(axrxx,
                                                          period[nzxx],
                                                          resp_z[nzxx,0,0].real,
                                                          None,
                                                          cxy, self.mted)
                                rerxy = self.plot_errorbar(axrxy,
                                                          period[nzxy],
                                                          resp_z[nzxy,0,1].real,
                                                          None,
                                                          cxy, self.mted)
                                reryx = self.plot_errorbar(axrxy,
                                                          period[nzyx],
                                                          resp_z[nzyx,1,0].real,
                                                          None,
                                                          cyx, self.mtmd)
                                reryy = self.plot_errorbar(axrxx,
                                                          period[nzyy],
                                                          resp_z[nzyy,1,1].real,
                                                          None,
                                                          cyx, self.mtmd)
                                #plot phase
                                rerxx = self.plot_errorbar(axpxx,
                                                          period[nzxx],
                                                          resp_z[nzxx,0,0].imag,
                                                          None,
                                                          cxy, self.mted)
                                rerxy = self.plot_errorbar(axpxy,
                                                          period[nzxy],
                                                          resp_z[nzxy,0,1].imag,
                                                          None,
                                                          cxy, self.mted)
                                reryx = self.plot_errorbar(axpxy,
                                                          period[nzyx],
                                                          resp_z[nzyx,1,0].imag,
                                                          None,
                                                          cyx, self.mtmd)
                                reryy = self.plot_errorbar(axpxx,
                                                          period[nzyy],
                                                          resp_z[nzyy,1,1].imag,
                                                          None,
                                                          cyx, self.mtmd)
                            line_list[0] += [rerxy[0], reryx[0]]
                            line_list[1] += [rerxx[0], reryy[0]]
                            label_list[0] += ['$Z^m_{xy}$ '+
                                               'rms={0:.2f}'.format(rms_xy),
                                              '$Z^m_{yx}$ '+
                                              'rms={0:.2f}'.format(rms_yx)]
                            label_list[1] += ['$Z^m_{xx}$ '+
                                               'rms={0:.2f}'.format(rms_xx),
                                              '$Z^m_{yy}$ '+
                                              'rms={0:.2f}'.format(rms_yy)]

                #make legends
                if self.plot_style == 1:
                    for aa, ax in enumerate(ax_list[0:self.plot_component]):
                        ax.legend(line_list[aa],
                                  label_list[aa],
                                  loc=self.legend_loc,
                                  markerscale=self.legend_marker_scale,
                                  borderaxespad=self.legend_border_axes_pad,
                                  labelspacing=self.legend_label_spacing,
                                  handletextpad=self.legend_handle_text_pad,
                                  borderpad=self.legend_border_pad,
                                  prop={'size':max([self.font_size/nr, 5])})
                if self.plot_style == 2:
                    if self.plot_component == 2:
                        axrxy.legend(line_list,
                                      label_list,
                                      loc=self.legend_loc,
                                      markerscale=self.legend_marker_scale,
                                      borderaxespad=self.legend_border_axes_pad,
                                      labelspacing=self.legend_label_spacing,
                                      handletextpad=self.legend_handle_text_pad,
                                      borderpad=self.legend_border_pad,
                                      prop={'size':max([self.font_size/nr, 5])})
                    else:
                        for aa, ax in enumerate(ax_list[0:self.plot_component/2]):
                            ax.legend(line_list[aa],
                                      label_list[aa],
                                      loc=self.legend_loc,
                                      markerscale=self.legend_marker_scale,
                                      borderaxespad=self.legend_border_axes_pad,
                                      labelspacing=self.legend_label_spacing,
                                      handletextpad=self.legend_handle_text_pad,
                                      borderpad=self.legend_border_pad,
                                      prop={'size':max([self.font_size/nr, 5])})
    def redraw_plot(self):
        """
        redraw plot if parameters were changed

        use this function if you updated some attributes and want to re-plot.

        :Example: ::

            >>> # change the color and marker of the xy components
            >>> import mtpy.modeling.occam2d as occam2d
            >>> ocd = occam2d.Occam2DData(r"/home/occam2d/Data.dat")
            >>> p1 = ocd.plotAllResponses()
            >>> #change line width
            >>> p1.lw = 2
            >>> p1.redraw_plot()
        """
        for fig in self.fig_list:
            plt.close(fig)
        self.plot()

    def save_figure(self, save_fn, file_format='pdf', orientation='portrait',
                  fig_dpi=None, close_fig='y'):
        """
        save_plot will save the figure to save_fn.

        Arguments:
        -----------

            **save_fn** : string
                          full path to save figure to, can be input as
                          * directory path -> the directory path to save to
                            in which the file will be saved as 
                            save_fn/station_name_PhaseTensor.file_format

                          * full path -> file will be save to the given 
                            path.  If you use this option then the format
                            will be assumed to be provided by the path

            **file_format** : [ pdf | eps | jpg | png | svg ]
                              file type of saved figure pdf,svg,eps... 

            **orientation** : [ landscape | portrait ]
                              orientation in which the file will be saved
                              *default* is portrait

            **fig_dpi** : int
                          The resolution in dots-per-inch the file will be
                          saved.  If None then the dpi will be that at 
                          which the figure was made.  I don't think that 
                          it can be larger than dpi of the figure.

            **close_plot** : [ y | n ]
                             * 'y' will close the plot after saving.
                             * 'n' will leave plot open

        :Example: ::

            >>> # to save plot as jpg
            >>> import mtpy.modeling.occam2d as occam2d
            >>> dfn = r"/home/occam2d/Inv1/data.dat"
            >>> ocd = occam2d.Occam2DData(dfn)
            >>> ps1 = ocd.plotPseudoSection()
            >>> ps1.save_plot(r'/home/MT/figures', file_format='jpg')

        """

        if fig_dpi == None:
            fig_dpi = self.fig_dpi

        if os.path.isdir(save_fn) == False:
            file_format = save_fn[-3:]
            self.fig.savefig(save_fn, dpi=fig_dpi, format=file_format,
                             orientation=orientation, bbox_inches='tight')

        else:
            save_fn = os.path.join(save_fn, '_L2.'+
                                    file_format)
            self.fig.savefig(save_fn, dpi=fig_dpi, format=file_format,
                        orientation=orientation, bbox_inches='tight')

        if close_fig == 'y':
            plt.clf()
            plt.close(self.fig)

        else:
            pass

        self.fig_fn = save_fn
        print('Saved figure to: '+self.fig_fn)

    def update_plot(self):
        """
        update any parameters that where changed using the built-in draw from
        canvas.  

        Use this if you change an of the .fig or axes properties

        :Example: ::

            >>> # to change the grid lines to only be on the major ticks
            >>> import mtpy.modeling.occam2d as occam2d
            >>> dfn = r"/home/occam2d/Inv1/data.dat"
            >>> ocd = occam2d.Occam2DData(dfn)
            >>> ps1 = ocd.plotAllResponses()
            >>> [ax.grid(True, which='major') for ax in [ps1.axrte,ps1.axtep]]
            >>> ps1.update_plot()

        """

        self.fig.canvas.draw()

    def __str__(self):
        """
        rewrite the string builtin to give a useful message
        """

        return ("Plots data vs model response computed by WS3DINV")

#==============================================================================
# plot depth slices
#==============================================================================
class PlotDepthSlice(object):
    """
    Plots depth slices of resistivity model

    :Example: ::

        >>> import mtpy.modeling.ws3dinv as ws
        >>> mfn = r"/home/MT/ws3dinv/Inv1/Test_model.00"
        >>> sfn = r"/home/MT/ws3dinv/Inv1/WSStationLocations.txt"
        >>> # plot just first layer to check the formating        
        >>> pds = ws.PlotDepthSlice(model_fn=mfn, station_fn=sfn, 
        >>> ...                     depth_index=0, save_plots='n')
        >>> #move color bar up 
        >>> pds.cb_location
        >>> (0.64500000000000002, 0.14999999999999997, 0.3, 0.025)
        >>> pds.cb_location = (.645, .175, .3, .025)
        >>> pds.redraw_plot()
        >>> #looks good now plot all depth slices and save them to a folder
        >>> pds.save_path = r"/home/MT/ws3dinv/Inv1/DepthSlices"
        >>> pds.depth_index = None
        >>> pds.save_plots = 'y'
        >>> pds.redraw_plot()

    ======================= ===================================================
    Attributes              Description    
    ======================= ===================================================
    cb_location             location of color bar (x, y, width, height)
                            *default* is None, automatically locates
    cb_orientation          [ 'vertical' | 'horizontal' ] 
                            *default* is horizontal 
    cb_pad                  padding between axes and colorbar
                            *default* is None
    cb_shrink               percentage to shrink colorbar by
                            *default* is None
    climits                 (min, max) of resistivity color on log scale
                            *default* is (0, 4)
    cmap                    name of color map *default* is 'jet_r'
    data_fn                 full path to data file
    depth_index             integer value of depth slice index, shallowest
                            layer is 0
    dscale                  scaling parameter depending on map_scale 
    ew_limits               (min, max) plot limits in e-w direction in 
                            map_scale units. *default* is None, sets viewing
                            area to the station area
    fig_aspect              aspect ratio of plot. *default* is 1
    fig_dpi                 resolution of figure in dots-per-inch. *default* is
                            300
    fig_list                list of matplotlib.figure instances for each 
                            depth slice                 
    fig_size                [width, height] in inches of figure size
                            *default* is [6, 6]
    font_size               size of ticklabel font in points, labels are 
                            font_size+2. *default* is 7
    grid_east               relative location of grid nodes in e-w direction
                            in map_scale units
    grid_north              relative location of grid nodes in n-s direction
                            in map_scale units
    grid_z                  relative location of grid nodes in z direction
                            in map_scale units
    initial_fn              full path to initial file
    map_scale               [ 'km' | 'm' ] distance units of map. *default* is 
                            km
    mesh_east               np.meshgrid(grid_east, grid_north, indexing='ij')
    mesh_north              np.meshgrid(grid_east, grid_north, indexing='ij')
    model_fn                full path to model file
    nodes_east              relative distance betwen nodes in e-w direction
                            in map_scale units
    nodes_north             relative distance betwen nodes in n-s direction
                            in map_scale units
    nodes_z                 relative distance betwen nodes in z direction
                            in map_scale units
    ns_limits               (min, max) plot limits in n-s direction in 
                            map_scale units. *default* is None, sets viewing
                            area to the station area
    plot_grid               [ 'y' | 'n' ] 'y' to plot mesh grid lines. 
                            *default* is 'n'
    plot_yn                 [ 'y' | 'n' ] 'y' to plot on instantiation
    res_model               np.ndarray(n_north, n_east, n_vertical) of 
                            model resistivity values in linear scale
    save_path               path to save figures to
    save_plots              [ 'y' | 'n' ] 'y' to save depth slices to save_path
    station_east            location of stations in east direction in 
                            map_scale units  
    station_fn              full path to station locations file
    station_names           station names
    station_north           location of station in north direction in 
                            map_scale units
    subplot_bottom          distance between axes and bottom of figure window
    subplot_left            distance between axes and left of figure window  
    subplot_right           distance between axes and right of figure window
    subplot_top             distance between axes and top of figure window
    title                   titiel of plot *default* is depth of slice
    xminorticks             location of xminorticks
    yminorticks             location of yminorticks
    ======================= ===================================================
    """

    def __init__(self, model_fn=None, data_fn=None, station_fn=None,
                 initial_fn=None, **kwargs):
        self.model_fn = model_fn
        self.data_fn = data_fn
        self.station_fn = station_fn
        self.initial_fn = initial_fn

        self.save_path = kwargs.pop('save_path', None)
        if self.model_fn is not None and self.save_path is None:
            self.save_path = os.path.dirname(self.model_fn)
        elif self.initial_fn is not None and self.save_path is None:
            self.save_path = os.path.dirname(self.initial_fn)

        if self.save_path is not None:
            if not os.path.exists(self.save_path):
                os.mkdir(self.save_path)

        self.save_plots = kwargs.pop('save_plots', 'y')

        self.depth_index = kwargs.pop('depth_index', None)
        self.map_scale = kwargs.pop('map_scale', 'km')
        #make map scale
        if self.map_scale=='km':
            self.dscale=1000.
        elif self.map_scale=='m':
            self.dscale=1.
        self.ew_limits = kwargs.pop('ew_limits', None)
        self.ns_limits = kwargs.pop('ns_limits', None)

        self.plot_grid = kwargs.pop('plot_grid', 'n')

        self.fig_size = kwargs.pop('fig_size', [6, 6])
        self.fig_dpi = kwargs.pop('dpi', 300)
        self.fig_aspect = kwargs.pop('fig_aspect', 1)
        self.title = kwargs.pop('title', 'on')
        self.fig_list = []

        self.xminorticks = kwargs.pop('xminorticks', 1000)
        self.yminorticks = kwargs.pop('yminorticks', 1000)

        self.climits = kwargs.pop('climits', (0,4))
        self.cmap = kwargs.pop('cmap', 'jet_r')
        self.font_size = kwargs.pop('font_size', 8)

        self.cb_shrink = kwargs.pop('cb_shrink', .8)
        self.cb_pad = kwargs.pop('cb_pad', .01)
        self.cb_orientation = kwargs.pop('cb_orientation', 'horizontal')
        self.cb_location = kwargs.pop('cb_location', None)

        self.subplot_right = .99
        self.subplot_left = .085
        self.subplot_top = .92
        self.subplot_bottom = .1

        self.res_model = None
        self.grid_east = None
        self.grid_north = None
        self.grid_z  = None

        self.nodes_east = None
        self.nodes_north = None
        self.nodes_z = None

        self.mesh_east = None
        self.mesh_north = None

        self.station_east = None
        self.station_north = None
        self.station_names = None

        self.plot_yn = kwargs.pop('plot_yn', 'y')
        if self.plot_yn == 'y':
            self.plot()

    def read_files(self):
        """
        read in the files to get appropriate information
        """
        #--> read in model file
        if self.model_fn is not None:
            if os.path.isfile(self.model_fn) == True:
                wsmodel = WSModel(self.model_fn)
                self.res_model = wsmodel.res_model
                self.grid_east = wsmodel.grid_east/self.dscale
                self.grid_north = wsmodel.grid_north/self.dscale
                self.grid_z = wsmodel.grid_z/self.dscale
                self.nodes_east = wsmodel.nodes_east/self.dscale
                self.nodes_north = wsmodel.nodes_north/self.dscale
                self.nodes_z = wsmodel.nodes_z/self.dscale
            else:
                raise mtex.MTpyError_file_handling(
                        '{0} does not exist, check path'.format(self.model_fn))

        #--> read in data file to get station locations
        if self.data_fn is not None:
            if os.path.isfile(self.data_fn) == True:
                wsdata = WSData()
                wsdata.read_data_file(self.data_fn)
                self.station_east = wsdata.data['east']/self.dscale
                self.station_north = wsdata.data['north']/self.dscale
                self.station_names = wsdata.data['station']
            else:
                print('Could not find data file {0}'.format(self.data_fn))

        #--> read in station file
        if self.station_fn is not None:
            if os.path.isfile(self.station_fn) == True:
                wsstations = WSStation(self.station_fn)
                wsstations.read_station_file()
                self.station_east = wsstations.east/self.dscale
                self.station_north = wsstations.north/self.dscale
                self.station_names = wsstations.names
            else:
                print('Could not find station file {0}'.format(self.station_fn))

        #--> read in initial file
        if self.initial_fn is not None:
            if os.path.isfile(self.initial_fn) == True:
                wsmesh = WSMesh()
                wsmesh.read_initial_file(self.initial_fn)
                self.grid_east = wsmesh.grid_east/self.dscale
                self.grid_north = wsmesh.grid_north/self.dscale
                self.grid_z = wsmesh.grid_z/self.dscale
                self.nodes_east = wsmesh.nodes_east/self.dscale
                self.nodes_north = wsmesh.nodes_north/self.dscale
                self.nodes_z = wsmesh.nodes_z/self.dscale

                #need to convert index values to resistivity values
                rdict = dict([(ii,res) for ii,res in enumerate(wsmesh.res_list,1)])

                for ii in range(len(wsmesh.res_list)):
                    self.res_model[np.where(wsmesh.res_model==ii+1)] = \
                                                                    rdict[ii+1]
            else:
                raise mtex.MTpyError_file_handling(
                     '{0} does not exist, check path'.format(self.initial_fn))

        if self.initial_fn is None and self.model_fn is None:
            raise mtex.MTpyError_inputarguments('Need to input either a model'
                                                ' file or initial file.')

    def plot(self):
        """
        plot depth slices
        """
        #--> get information from files
        self.read_files()

        fdict = {'size':self.font_size+2, 'weight':'bold'}

        cblabeldict={-2:'$10^{-3}$',-1:'$10^{-1}$',0:'$10^{0}$',1:'$10^{1}$',
                     2:'$10^{2}$',3:'$10^{3}$',4:'$10^{4}$',5:'$10^{5}$',
                     6:'$10^{6}$',7:'$10^{7}$',8:'$10^{8}$'}

        #create an list of depth slices to plot
        if self.depth_index == None:
            zrange = list(range(self.grid_z.shape[0]))
        elif type(self.depth_index) is int:
            zrange = [self.depth_index]
        elif type(self.depth_index) is list or \
             type(self.depth_index) is np.ndarray:
            zrange = self.depth_index

        #set the limits of the plot
        if self.ew_limits == None:
            if self.station_east is not None:
                xlimits = (np.floor(self.station_east.min()),
                           np.ceil(self.station_east.max()))
            else:
                xlimits = (self.grid_east[5], self.grid_east[-5])
        else:
            xlimits = self.ew_limits

        if self.ns_limits == None:
            if self.station_north is not None:
                ylimits = (np.floor(self.station_north.min()),
                           np.ceil(self.station_north.max()))
            else:
                ylimits = (self.grid_north[5], self.grid_north[-5])
        else:
            ylimits = self.ns_limits


        #make a mesh grid of north and east
        self.mesh_east, self.mesh_north = np.meshgrid(self.grid_east,
                                                      self.grid_north,
                                                      indexing='ij')

        plt.rcParams['font.size'] = self.font_size

        #--> plot depths into individual figures
        for ii in zrange:
            depth = '{0:.3f} ({1})'.format(self.grid_z[ii],
                                     self.map_scale)
            fig = plt.figure(depth, figsize=self.fig_size, dpi=self.fig_dpi)
            plt.clf()
            ax1 = fig.add_subplot(1, 1, 1, aspect=self.fig_aspect)
            plot_res = np.log10(self.res_model[:, :, ii].T)
            mesh_plot = ax1.pcolormesh(self.mesh_east,
                                       self.mesh_north,
                                       plot_res,
                                       cmap=self.cmap,
                                       vmin=self.climits[0],
                                       vmax=self.climits[1])

            #plot the stations
            if self.station_east is not None:
                for ee, nn in zip(self.station_east, self.station_north):
                    ax1.text(ee, nn, '*',
                             verticalalignment='center',
                             horizontalalignment='center',
                             fontdict={'size':5, 'weight':'bold'})

            #set axis properties
            ax1.set_xlim(xlimits)
            ax1.set_ylim(ylimits)
            ax1.xaxis.set_minor_locator(MultipleLocator(self.xminorticks/self.dscale))
            ax1.yaxis.set_minor_locator(MultipleLocator(self.yminorticks/self.dscale))
            ax1.set_ylabel('Northing ('+self.map_scale+')',fontdict=fdict)
            ax1.set_xlabel('Easting ('+self.map_scale+')',fontdict=fdict)
            ax1.set_title('Depth = {0}'.format(depth), fontdict=fdict)

            #plot the grid if desired
            if self.plot_grid == 'y':
                east_line_xlist = []
                east_line_ylist = []
                for xx in self.grid_east:
                    east_line_xlist.extend([xx, xx])
                    east_line_xlist.append(None)
                    east_line_ylist.extend([self.grid_north.min(),
                                            self.grid_north.max()])
                    east_line_ylist.append(None)
                ax1.plot(east_line_xlist,
                              east_line_ylist,
                              lw=.25,
                              color='k')

                north_line_xlist = []
                north_line_ylist = []
                for yy in self.grid_north:
                    north_line_xlist.extend([self.grid_east.min(),
                                             self.grid_east.max()])
                    north_line_xlist.append(None)
                    north_line_ylist.extend([yy, yy])
                    north_line_ylist.append(None)
                ax1.plot(north_line_xlist,
                              north_line_ylist,
                              lw=.25,
                              color='k')


            #plot the colorbar
            if self.cb_location is None:
                if self.cb_orientation == 'horizontal':
                    self.cb_location = (ax1.axes.figbox.bounds[3]-.225,
                                        ax1.axes.figbox.bounds[1]+.05,.3,.025)

                elif self.cb_orientation == 'vertical':
                    self.cb_location = ((ax1.axes.figbox.bounds[2]-.15,
                                        ax1.axes.figbox.bounds[3]-.21,.025,.3))

            ax2 = fig.add_axes(self.cb_location)

            cb = mcb.ColorbarBase(ax2,
                                  cmap=self.cmap,
                                  norm=Normalize(vmin=self.climits[0],
                                                 vmax=self.climits[1]),
                                  orientation=self.cb_orientation)

            if self.cb_orientation == 'horizontal':
                cb.ax.xaxis.set_label_position('top')
                cb.ax.xaxis.set_label_coords(.5,1.3)


            elif self.cb_orientation == 'vertical':
                cb.ax.yaxis.set_label_position('right')
                cb.ax.yaxis.set_label_coords(1.25,.5)
                cb.ax.yaxis.tick_left()
                cb.ax.tick_params(axis='y',direction='in')

            cb.set_label('Resistivity ($\Omega \cdot$m)',
                         fontdict={'size':self.font_size+1})
            cb.set_ticks(np.arange(self.climits[0],self.climits[1]+1))
            cb.set_ticklabels([cblabeldict[cc]
                                for cc in np.arange(self.climits[0],
                                                    self.climits[1]+1)])

            self.fig_list.append(fig)

            #--> save plots to a common folder
            if self.save_plots == 'y':

                fig.savefig(os.path.join(self.save_path,
                            "Depth_{}_{:.4f}.png".format(ii, self.grid_z[ii])),
                            dpi=self.fig_dpi, bbox_inches='tight')
                fig.clear()
                plt.close()

            else:
                pass

    def redraw_plot(self):
        """
        redraw plot if parameters were changed

        use this function if you updated some attributes and want to re-plot.

        :Example: ::

            >>> # change the color and marker of the xy components
            >>> import mtpy.modeling.occam2d as occam2d
            >>> ocd = occam2d.Occam2DData(r"/home/occam2d/Data.dat")
            >>> p1 = ocd.plotAllResponses()
            >>> #change line width
            >>> p1.lw = 2
            >>> p1.redraw_plot()
        """
        for fig in self.fig_list:
            plt.close(fig)
        self.plot()

    def update_plot(self, fig):
        """
        update any parameters that where changed using the built-in draw from
        canvas.  

        Use this if you change an of the .fig or axes properties

        :Example: ::

            >>> # to change the grid lines to only be on the major ticks
            >>> import mtpy.modeling.occam2d as occam2d
            >>> dfn = r"/home/occam2d/Inv1/data.dat"
            >>> ocd = occam2d.Occam2DData(dfn)
            >>> ps1 = ocd.plotAllResponses()
            >>> [ax.grid(True, which='major') for ax in [ps1.axrte,ps1.axtep]]
            >>> ps1.update_plot()

        """

        fig.canvas.draw()

    def __str__(self):
        """
        rewrite the string builtin to give a useful message
        """

        return ("Plots depth slices of model from WS3DINV")

#==============================================================================
# plot phase tensors
#==============================================================================
class PlotPTMaps(mtplottools.MTEllipse):
    """
    Plot phase tensor maps including residual pt if response file is input.

    :Plot only data for one period: ::

        >>> import mtpy.modeling.ws3dinv as ws
        >>> dfn = r"/home/MT/ws3dinv/Inv1/WSDataFile.dat"
        >>> ptm = ws.PlotPTMaps(data_fn=dfn, plot_period_list=[0])

    :Plot data and model response: ::

        >>> import mtpy.modeling.ws3dinv as ws
        >>> dfn = r"/home/MT/ws3dinv/Inv1/WSDataFile.dat"
        >>> rfn = r"/home/MT/ws3dinv/Inv1/Test_resp.00"
        >>> mfn = r"/home/MT/ws3dinv/Inv1/Test_model.00"
        >>> ptm = ws.PlotPTMaps(data_fn=dfn, resp_fn=rfn, model_fn=mfn,
        >>> ...                 plot_period_list=[0])
        >>> # adjust colorbar
        >>> ptm.cb_res_pad = 1.25
        >>> ptm.redraw_plot()


    ========================== ================================================
    Attributes                 Description    
    ========================== ================================================
    cb_pt_pad                  percentage from top of axes to place pt 
                               color bar. *default* is .90
    cb_res_pad                 percentage from bottom of axes to place
                               resistivity color bar. *default* is 1.2
    cb_residual_tick_step      tick step for residual pt. *default* is 3
    cb_tick_step               tick step for phase tensor color bar, 
                               *default* is 45
    data                       np.ndarray(n_station, n_periods, 2, 2)
                               impedance tensors for station data                     
    data_fn                    full path to data fle               
    dscale                     scaling parameter depending on map_scale
    ellipse_cmap               color map for pt ellipses. *default* is
                               mt_bl2gr2rd
    ellipse_colorby            [ 'skew' | 'skew_seg' | 'phimin' | 'phimax'|
                                 'phidet' | 'ellipticity' ] parameter to color
                                 ellipses by. *default* is 'phimin'
    ellipse_range              (min, max, step) min and max of colormap, need
                               to input step if plotting skew_seg
    ellipse_size               relative size of ellipses in map_scale
    ew_limits                  limits of plot in e-w direction in map_scale
                               units.  *default* is None, scales to station 
                               area
    fig_aspect                 aspect of figure. *default* is 1
    fig_dpi                    resolution in dots-per-inch. *default* is 300
    fig_list                   list of matplotlib.figure instances for each
                               figure plotted.
    fig_size                   [width, height] in inches of figure window
                               *default* is [6, 6]
    font_size                  font size of ticklabels, axes labels are 
                               font_size+2. *default* is 7
    grid_east                  relative location of grid nodes in e-w direction
                               in map_scale units
    grid_north                 relative location of grid nodes in n-s direction
                               in map_scale units
    grid_z                     relative location of grid nodes in z direction
                               in map_scale units
    initial_fn                 full path to initial file
    map_scale                  [ 'km' | 'm' ] distance units of map. 
                               *default* is km
    mesh_east                  np.meshgrid(grid_east, grid_north, indexing='ij')
    mesh_north                 np.meshgrid(grid_east, grid_north, indexing='ij')
    model_fn                   full path to model file
    nodes_east                 relative distance betwen nodes in e-w direction
                               in map_scale units
    nodes_north                relative distance betwen nodes in n-s direction
                               in map_scale units
    nodes_z                    relative distance betwen nodes in z direction
                               in map_scale units
    ns_limits                  (min, max) limits of plot in n-s direction
                               *default* is None, viewing area is station area
    pad_east                   padding from extreme stations in east direction
    pad_north                  padding from extreme stations in north direction
    period_list                list of periods from data
    plot_grid                  [ 'y' | 'n' ] 'y' to plot grid lines
                               *default* is 'n'
    plot_period_list           list of period index values to plot
                               *default* is None 
    plot_yn                    ['y' | 'n' ] 'y' to plot on instantiation
                               *default* is 'y'
    res_cmap                   colormap for resisitivity values. 
                               *default* is 'jet_r'
    res_limits                 (min, max) resistivity limits in log scale
                               *default* is (0, 4)
    res_model                  np.ndarray(n_north, n_east, n_vertical) of 
                               model resistivity values in linear scale
    residual_cmap              color map for pt residuals. 
                               *default* is 'mt_wh2or' 
    resp                       np.ndarray(n_stations, n_periods, 2, 2)
                               impedance tensors for model response  
    resp_fn                    full path to response file
    save_path                  directory to save figures to
    save_plots                 [ 'y' | 'n' ] 'y' to save plots to save_path
    station_east               location of stations in east direction in 
                               map_scale units  
    station_fn                 full path to station locations file
    station_names              station names
    station_north              location of station in north direction in 
                               map_scale units
    subplot_bottom             distance between axes and bottom of figure window
    subplot_left               distance between axes and left of figure window  
    subplot_right              distance between axes and right of figure window
    subplot_top                distance between axes and top of figure window
    title                      titiel of plot *default* is depth of slice
    xminorticks                location of xminorticks
    yminorticks                location of yminorticks
    ========================== ================================================
    """

    def __init__(self, data_fn=None, resp_fn=None, station_fn=None,
                 model_fn=None, initial_fn=None, **kwargs):

        self.model_fn = model_fn
        self.data_fn = data_fn
        self.station_fn = station_fn
        self.resp_fn = resp_fn
        self.initial_fn = initial_fn

        self.save_path = kwargs.pop('save_path', None)
        if self.model_fn is not None and self.save_path is None:
            self.save_path = os.path.dirname(self.model_fn)
        elif self.initial_fn is not None and self.save_path is None:
            self.save_path = os.path.dirname(self.initial_fn)

        if self.save_path is not None:
            if not os.path.exists(self.save_path):
                os.mkdir(self.save_path)

        self.save_plots = kwargs.pop('save_plots', 'y')
        self.plot_period_list = kwargs.pop('plot_period_list', None)

        self.map_scale = kwargs.pop('map_scale', 'km')
        #make map scale
        if self.map_scale=='km':
            self.dscale=1000.
        elif self.map_scale=='m':
            self.dscale=1.
        self.ew_limits = kwargs.pop('ew_limits', None)
        self.ns_limits = kwargs.pop('ns_limits', None)

        self.pad_east = kwargs.pop('pad_east', 2)
        self.pad_north = kwargs.pop('pad_north', 2)

        self.plot_grid = kwargs.pop('plot_grid', 'n')

        self.fig_num = kwargs.pop('fig_num', 1)
        self.fig_size = kwargs.pop('fig_size', [6, 6])
        self.fig_dpi = kwargs.pop('dpi', 300)
        self.fig_aspect = kwargs.pop('fig_aspect', 1)
        self.title = kwargs.pop('title', 'on')
        self.fig_list = []

        self.xminorticks = kwargs.pop('xminorticks', 1000)
        self.yminorticks = kwargs.pop('yminorticks', 1000)

        self.residual_cmap = kwargs.pop('residual_cmap', 'mt_wh2or')
        self.font_size = kwargs.pop('font_size', 7)

        self.cb_tick_step = kwargs.pop('cb_tick_step', 45)
        self.cb_residual_tick_step = kwargs.pop('cb_residual_tick_step', 3)
        self.cb_pt_pad = kwargs.pop('cb_pt_pad', .90)
        self.cb_res_pad = kwargs.pop('cb_res_pad', 1.22)


        self.res_limits = kwargs.pop('res_limits', (0,4))
        self.res_cmap = kwargs.pop('res_cmap', 'jet_r')

        #--> set the ellipse properties -------------------
        self._ellipse_dict = kwargs.pop('ellipse_dict', {})
        self._read_ellipse_dict()

        self.subplot_right = .99
        self.subplot_left = .085
        self.subplot_top = .92
        self.subplot_bottom = .1
        self.subplot_hspace = .2
        self.subplot_wspace = .05

        self.res_model = None
        self.grid_east = None
        self.grid_north = None
        self.grid_z  = None

        self.nodes_east = None
        self.nodes_north = None
        self.nodes_z = None

        self.mesh_east = None
        self.mesh_north = None

        self.station_east = None
        self.station_north = None
        self.station_names = None

        self.data = None
        self.resp = None
        self.period_list = None

        self.plot_yn = kwargs.pop('plot_yn', 'y')
        if self.plot_yn == 'y':
            self.plot()

    def _get_pt(self):
        """
        get phase tensors
        """

        #--> read in data file
        if self.data_fn is None:
            raise mtex.MTpyError_inputarguments('Need to input a data file')
        wsdata = WSData()
        wsdata.read_data_file(self.data_fn, station_fn=self.station_fn)
        self.data = wsdata.z_data
        self.period_list = wsdata.period_list
        self.station_east = wsdata.station_east/self.dscale
        self.station_north = wsdata.station_north/self.dscale
        self.station_names = wsdata.station_names

        if self.plot_period_list is None:
            self.plot_period_list = self.period_list
        else:
            if type(self.plot_period_list) is list:
                #check if entries are index values or actual periods
                if type(self.plot_period_list[0]) is int:
                    self.plot_period_list = [self.period_list[ii]
                                             for ii in self.plot_period_list]
                else:
                    pass
            elif type(self.plot_period_list) is int:
                self.plot_period_list = self.period_list[self.plot_period_list]

        #--> read model file
        if self.model_fn is not None:
            wsmodel = WSModel(self.model_fn)
            self.res_model = wsmodel.res_model
            self.grid_east = wsmodel.grid_east/self.dscale
            self.grid_north = wsmodel.grid_north/self.dscale
            self.grid_z = wsmodel.grid_z/self.dscale
            self.mesh_east, self.mesh_north = np.meshgrid(self.grid_east,
                                                          self.grid_north,
                                                          indexing='ij')

        #--> read response file
        if self.resp_fn is not None:
            wsresp = WSResponse(self.resp_fn)
            self.resp = wsresp.z_resp




    def plot(self):
        """
        plot phase tensor maps for data and or response, each figure is of a
        different period.  If response is input a third column is added which is 
        the residual phase tensor showing where the model is not fitting the data 
        well.  The data is plotted in km in units of ohm-m.

        Inputs:
            data_fn = full path to data file
            resp_fn = full path to response file, if none just plots data
            sites_fn = full path to sites file
            periodlst = indicies of periods you want to plot
            esize = size of ellipses as:
                    0 = phase tensor ellipse
                    1 = phase tensor residual
                    2 = resistivity tensor ellipse
                    3 = resistivity tensor residual
            ecolor = 'phimin' for coloring with phimin or 'beta' for beta coloring
            colormm = list of min and max coloring for plot, list as follows:
                    0 = phase tensor min and max for ecolor in degrees
                    1 = phase tensor residual min and max [0,1]
                    2 = resistivity tensor coloring as resistivity on log scale
                    3 = resistivity tensor residual coloring as resistivity on 
                        linear scale
            xpad = padding of map from stations at extremities (km)
            units = 'mv' to convert to Ohm-m 
            dpi = dots per inch of figure
        """

        self._get_pt()
        plt.rcParams['font.size'] = self.font_size
        plt.rcParams['figure.subplot.left'] = self.subplot_left
        plt.rcParams['figure.subplot.right'] = self.subplot_right
        plt.rcParams['figure.subplot.bottom'] = self.subplot_bottom
        plt.rcParams['figure.subplot.top'] = self.subplot_top

        gs = gridspec.GridSpec(1, 3, hspace=self.subplot_hspace,
                               wspace=self.subplot_wspace)

        font_dict = {'size':self.font_size+2, 'weight':'bold'}
        n_stations = self.data.shape[0]
        #set some local parameters
        ckmin = float(self.ellipse_range[0])
        ckmax = float(self.ellipse_range[1])
        try:
            ckstep = float(self.ellipse_range[2])
        except IndexError:
            if self.ellipse_cmap == 'mt_seg_bl2wh2rd':
                raise ValueError('Need to input range as (min, max, step)')
            else:
                ckstep = 3
        nseg = float((ckmax-ckmin)/(2*ckstep))

        if self.ew_limits == None:
            if self.station_east is not None:
                self.ew_limits = (np.floor(self.station_east.min())-
                                                      self.pad_east,
                                  np.ceil(self.station_east.max())+
                                                      self.pad_east)
            else:
                self.ew_limits = (self.grid_east[5], self.grid_east[-5])

        if self.ns_limits == None:
            if self.station_north is not None:
                self.ns_limits = (np.floor(self.station_north.min())-
                                                        self.pad_north,
                                  np.ceil(self.station_north.max())+
                                                        self.pad_north)
            else:
                self.ns_limits = (self.grid_north[5], self.grid_north[-5])

        for ff, per in enumerate(self.plot_period_list):
            print('Plotting Period: {0:.5g}'.format(per))
            fig = plt.figure('{0:.5g}'.format(per), figsize=self.fig_size,
                             dpi=self.fig_dpi)
            fig.clf()

            if self.resp_fn is not None:
                axd = fig.add_subplot(gs[0, 0], aspect='equal')
                axm = fig.add_subplot(gs[0, 1], aspect='equal')
                axr = fig.add_subplot(gs[0, 2], aspect='equal')
                ax_list = [axd, axm, axr]

            else:
                axd = fig.add_subplot(gs[0, :], aspect='equal')
                ax_list = [axd]

            #plot model below the phase tensors
            if self.model_fn is not None:
                approx_depth, d_index = estimate_skin_depth(self.res_model,
                                                            self.grid_z,
                                                            per,
                                                            dscale=self.dscale)
                for ax in ax_list:
                    plot_res = np.log10(self.res_model[:, :, d_index].T)
                    ax.pcolormesh(self.mesh_east,
                                   self.mesh_north,
                                   plot_res,
                                   cmap=self.res_cmap,
                                   vmin=self.res_limits[0],
                                   vmax=self.res_limits[1])


            #--> get phase tensors
            pt = mtpt.PhaseTensor(z_array=self.data[:, ff],
                                  freq=np.repeat(per, n_stations))
            if self.resp is not None:
                mpt = mtpt.PhaseTensor(z_array=self.resp[:, ff],
                                       freq=np.repeat(per, n_stations))
                rpt = mtpt.ResidualPhaseTensor(pt_object1=pt, pt_object2=mpt)
                rpt = rpt.residual_pt
                rcarray = np.sqrt(abs(rpt.phimin[0]*rpt.phimax[0]))
                rcmin = np.floor(rcarray.min())
                rcmax = np.floor(rcarray.max())

            #--> get color array
            if self.ellipse_cmap == 'mt_seg_bl2wh2rd':
                bounds = np.arange(ckmin, ckmax+ckstep, ckstep)
                nseg = float((ckmax-ckmin)/(2*ckstep))

            #get the properties to color the ellipses by
            if self.ellipse_colorby == 'phiminang' or \
               self.ellipse_colorby == 'phimin':
                colorarray = pt.phimin[0]
                if self.resp is not None:
                    mcarray = mpt.phimin[0]

            elif self.ellipse_colorby == 'phidet':
                 colorarray = np.sqrt(abs(pt.det[0]))*(180/np.pi)
                 if self.resp is not None:
                    mcarray = np.sqrt(abs(mpt.det[0]))*(180/np.pi)


            elif self.ellipse_colorby == 'skew' or\
                 self.ellipse_colorby == 'skew_seg':
                colorarray = pt.beta[0]
                if self.resp is not None:
                    mcarray = mpt.beta[0]

            elif self.ellipse_colorby == 'ellipticity':
                colorarray = pt.ellipticity[0]
                if self.resp is not None:
                    mcarray = mpt.ellipticity[0]

            else:
                raise NameError(self.ellipse_colorby+' is not supported')


            #--> plot phase tensor ellipses for each stations
            for jj in range(n_stations):
                #-----------plot data phase tensors---------------
                eheight = pt.phimin[0][jj]/pt.phimax[0].max()*self.ellipse_size
                ewidth = pt.phimax[0][jj]/pt.phimax[0].max()*self.ellipse_size

                ellipse = Ellipse((self.station_east[jj],
                                   self.station_north[jj]),
                                   width=ewidth,
                                   height=eheight,
                                   angle=90-pt.azimuth[0][jj])

                #get ellipse color
                if self.ellipse_cmap.find('seg')>0:
                    ellipse.set_facecolor(mtcl.get_plot_color(colorarray[jj],
                                                         self.ellipse_colorby,
                                                         self.ellipse_cmap,
                                                         ckmin,
                                                         ckmax,
                                                         bounds=bounds))
                else:
                    ellipse.set_facecolor(mtcl.get_plot_color(colorarray[jj],
                                                         self.ellipse_colorby,
                                                         self.ellipse_cmap,
                                                         ckmin,
                                                         ckmax))

                axd.add_artist(ellipse)
                if self.resp is not None:
                    #-----------plot response phase tensors---------------
                    eheight = mpt.phimin[0][jj]/mpt.phimax[0].max()*\
                              self.ellipse_size
                    ewidth = mpt.phimax[0][jj]/mpt.phimax[0].max()*\
                              self.ellipse_size

                    ellipsem = Ellipse((self.station_east[jj],
                                       self.station_north[jj]),
                                       width=ewidth,
                                       height=eheight,
                                       angle=90-mpt.azimuth[0][jj])

                    #get ellipse color
                    if self.ellipse_cmap.find('seg')>0:
                        ellipsem.set_facecolor(mtcl.get_plot_color(mcarray[jj],
                                                         self.ellipse_colorby,
                                                         self.ellipse_cmap,
                                                         ckmin,
                                                         ckmax,
                                                         bounds=bounds))
                    else:
                        ellipsem.set_facecolor(mtcl.get_plot_color(mcarray[jj],
                                                         self.ellipse_colorby,
                                                         self.ellipse_cmap,
                                                         ckmin,
                                                         ckmax))

                    axm.add_artist(ellipsem)
                    #-----------plot residual phase tensors---------------
                    eheight = rpt.phimin[0][jj]/rpt.phimax[0].max()*\
                                self.ellipse_size
                    ewidth = rpt.phimax[0][jj]/rpt.phimax[0].max()*\
                                self.ellipse_size

                    ellipser = Ellipse((self.station_east[jj],
                                       self.station_north[jj]),
                                       width=ewidth,
                                       height=eheight,
                                       angle=rpt.azimuth[0][jj])

                    #get ellipse color
                    if self.ellipse_cmap.find('seg')>0:
                        ellipser.set_facecolor(mtcl.get_plot_color(rcarray[jj],
                                                     self.ellipse_colorby,
                                                     self.residual_cmap,
                                                     rcmin,
                                                     rcmax,
                                                     bounds=bounds))
                    else:
                        ellipser.set_facecolor(mtcl.get_plot_color(rcarray[jj],
                                                     self.ellipse_colorby,
                                                     self.residual_cmap,
                                                     rcmin,
                                                     rcmax))

                    axr.add_artist(ellipser)

            #--> set axes properties
            # data
            axd.set_xlim(self.ew_limits)
            axd.set_ylim(self.ns_limits)
            axd.set_xlabel('Easting ({0})'.format(self.map_scale),
                           fontdict=font_dict)
            axd.set_ylabel('Northing ({0})'.format(self.map_scale),
                           fontdict=font_dict)
            #make a colorbar for phase tensors
            #bb = axd.axes.get_position().bounds
            bb = axd.get_position().bounds
            y1 = .25*(2+(self.ns_limits[1]-self.ns_limits[0])/
                     (self.ew_limits[1]-self.ew_limits[0]))
            cb_location = (3.35*bb[2]/5+bb[0],
                            y1*self.cb_pt_pad, .295*bb[2], .02)
            cbaxd = fig.add_axes(cb_location)
            cbd = mcb.ColorbarBase(cbaxd,
                                   cmap=mtcl.cmapdict[self.ellipse_cmap],
                                   norm=Normalize(vmin=ckmin,
                                                  vmax=ckmax),
                                   orientation='horizontal')
            cbd.ax.xaxis.set_label_position('top')
            cbd.ax.xaxis.set_label_coords(.5, 1.75)
            cbd.set_label(mtplottools.ckdict[self.ellipse_colorby])
            cbd.set_ticks(np.arange(ckmin, ckmax+self.cb_tick_step,
                                    self.cb_tick_step))

            axd.text(self.ew_limits[0]*.95,
                     self.ns_limits[1]*.95,
                     'Data',
                     horizontalalignment='left',
                     verticalalignment='top',
                     bbox={'facecolor':'white'},
                     fontdict={'size':self.font_size+1})

            #Model and residual
            if self.resp is not None:
                for aa, ax in enumerate([axm, axr]):
                    ax.set_xlim(self.ew_limits)
                    ax.set_ylim(self.ns_limits)
                    ax.set_xlabel('Easting ({0})'.format(self.map_scale),
                                   fontdict=font_dict)
                    plt.setp(ax.yaxis.get_ticklabels(), visible=False)
                    #make a colorbar ontop of axis
                    bb = ax.axes.get_position().bounds
                    y1 = .25*(2+(self.ns_limits[1]-self.ns_limits[0])/
                     (self.ew_limits[1]-self.ew_limits[0]))
                    cb_location = (3.35*bb[2]/5+bb[0],
                                   y1*self.cb_pt_pad, .295*bb[2], .02)
                    cbax = fig.add_axes(cb_location)
                    if aa == 0:
                        cb = mcb.ColorbarBase(cbax,
                                              cmap=mtcl.cmapdict[self.ellipse_cmap],
                                              norm=Normalize(vmin=ckmin,
                                                             vmax=ckmax),
                                               orientation='horizontal')
                        cb.ax.xaxis.set_label_position('top')
                        cb.ax.xaxis.set_label_coords(.5, 1.75)
                        cb.set_label(mtplottools.ckdict[self.ellipse_colorby])
                        cb.set_ticks(np.arange(ckmin, ckmax+self.cb_tick_step,
                                    self.cb_tick_step))
                        ax.text(self.ew_limits[0]*.95,
                                self.ns_limits[1]*.95,
                                'Model',
                                horizontalalignment='left',
                                verticalalignment='top',
                                bbox={'facecolor':'white'},
                                 fontdict={'size':self.font_size+1})
                    else:
                        cb = mcb.ColorbarBase(cbax,
                                              cmap=mtcl.cmapdict[self.residual_cmap],
                                               norm=Normalize(vmin=rcmin,
                                                              vmax=rcmax),
                                               orientation='horizontal')
                        cb.ax.xaxis.set_label_position('top')
                        cb.ax.xaxis.set_label_coords(.5, 1.75)
                        cb.set_label(r"$\sqrt{\Phi_{min} \Phi_{max}}$")
                        cb_ticks = np.arange(rcmin,
                                             rcmax+self.cb_residual_tick_step,
                                             self.cb_residual_tick_step)
                        cb.set_ticks(cb_ticks)
                        ax.text(self.ew_limits[0]*.95,
                                self.ns_limits[1]*.95,
                                'Residual',
                                horizontalalignment='left',
                                verticalalignment='top',
                                bbox={'facecolor':'white'},
                                fontdict={'size':self.font_size+1})

            if self.model_fn is not None:
                for ax in ax_list:
                    ax.tick_params(direction='out')
                    bb = ax.axes.get_position().bounds
                    y1 = .25*(2-(self.ns_limits[1]-self.ns_limits[0])/
                             (self.ew_limits[1]-self.ew_limits[0]))
                    cb_position = (3.0*bb[2]/5+bb[0],
                                   y1*self.cb_res_pad, .35*bb[2], .02)
                    cbax = fig.add_axes(cb_position)
                    cb = mcb.ColorbarBase(cbax,
                                          cmap=self.res_cmap,
                                          norm=Normalize(vmin=self.res_limits[0],
                                                         vmax=self.res_limits[1]),
                                          orientation='horizontal')
                    cb.ax.xaxis.set_label_position('top')
                    cb.ax.xaxis.set_label_coords(.5, 1.5)
                    cb.set_label('Resistivity ($\Omega \cdot$m)')
                    cb_ticks = np.arange(np.floor(self.res_limits[0]),
                                         np.ceil(self.res_limits[1]+1), 1)
                    cb.set_ticks(cb_ticks)
                    cb.set_ticklabels([mtplottools.labeldict[ctk] for ctk in cb_ticks])



            plt.show()
            self.fig_list.append(fig)

    def redraw_plot(self):
        """
        redraw plot if parameters were changed

        use this function if you updated some attributes and want to re-plot.

        :Example: ::

            >>> # change the color and marker of the xy components
            >>> import mtpy.modeling.occam2d as occam2d
            >>> ocd = occam2d.Occam2DData(r"/home/occam2d/Data.dat")
            >>> p1 = ocd.plotAllResponses()
            >>> #change line width
            >>> p1.lw = 2
            >>> p1.redraw_plot()
        """
        for fig in self.fig_list:
            plt.close(fig)
        self.plot()

    def save_figure(self, save_path=None, fig_dpi=None, file_format='pdf',
                    orientation='landscape', close_fig='y'):
        """
        save_figure will save the figure to save_fn.

        Arguments:
        -----------

            **save_fn** : string
                          full path to save figure to, can be input as
                          * directory path -> the directory path to save to
                            in which the file will be saved as 
                            save_fn/station_name_PhaseTensor.file_format

                          * full path -> file will be save to the given 
                            path.  If you use this option then the format
                            will be assumed to be provided by the path

            **file_format** : [ pdf | eps | jpg | png | svg ]
                              file type of saved figure pdf,svg,eps... 

            **orientation** : [ landscape | portrait ]
                              orientation in which the file will be saved
                              *default* is portrait

            **fig_dpi** : int
                          The resolution in dots-per-inch the file will be
                          saved.  If None then the dpi will be that at 
                          which the figure was made.  I don't think that 
                          it can be larger than dpi of the figure.

            **close_plot** : [ y | n ]
                             * 'y' will close the plot after saving.
                             * 'n' will leave plot open

        :Example: ::

            >>> # to save plot as jpg
            >>> import mtpy.modeling.occam2d as occam2d
            >>> dfn = r"/home/occam2d/Inv1/data.dat"
            >>> ocd = occam2d.Occam2DData(dfn)
            >>> ps1 = ocd.plotPseudoSection()
            >>> ps1.save_plot(r'/home/MT/figures', file_format='jpg')

        """

        if fig_dpi == None:
            fig_dpi = self.fig_dpi

        if os.path.isdir(save_path) == False:
            try:
                os.mkdir(save_path)
            except:
                raise IOError('Need to input a correct directory path')

        for fig in self.fig_list:
            per = fig.canvas.get_window_title()
            save_fn = os.path.join(save_path, 'PT_DepthSlice_{0}s.{1}'.format(
                                    per, file_format))
            fig.savefig(save_fn, dpi=fig_dpi, format=file_format,
                        orientation=orientation, bbox_inches='tight')

            if close_fig == 'y':
                plt.close(fig)

            else:
                pass

            self.fig_fn = save_fn
            print('Saved figure to: '+self.fig_fn)

#==============================================================================
# ESTIMATE SKIN DEPTH FOR MODEL
#==============================================================================
def estimate_skin_depth(res_model, grid_z, period, dscale=1000):
    """
    estimate the skin depth from the resistivity model assuming that

        delta_skin ~ 500 * sqrt(rho_a*T)

    Arguments:
    -----------
        **resmodel** : np.ndarray (n_north, n_east, n_z)
                       array of resistivity values for model grid

        **grid_z** : np.ndarray (n_z)
                     array of depth layers in m or km, be sure to change
                     dscale accordingly

        **period** : float
                     period in seconds to estimate a skin depth for

        **dscale** : [1000 | 1]
                     scaling value to scale depth estimation to meters (1) or
                     kilometers (1000)

    Outputs:
    ---------
        **depth** : float
                    estimated skin depth in units according to dscale

        **depth_index** : int
                          index value of grid_z that corresponds to the 
                          estimated skin depth.
    """
    if dscale == 1000:
        ms = 'km'
        ds = .5
    if dscale == 1:
        ms = 'm'
        ds = 500.
    #find the apparent resisitivity of each depth slice within the station area
    apparent_res_xy = np.array([res_model[6:-6, 6:-6, 0:ii+1].mean()
                                        for ii in range(grid_z.shape[0])])

    #calculate the period for each skin depth
    skin_depth_period = np.array([(zz/ds)**2*(1/rho_a)
                                for zz, rho_a in zip(grid_z, apparent_res_xy)])

    #match the period
    try:
        period_index = np.where(skin_depth_period >= period)[0][0]
    except IndexError:
        period_index = len(skin_depth_period)-1

    #get the depth slice
    depth = grid_z[period_index]

    print('-'*60)
    print(' input period                   {0:.6g} (s)'.format(period))
    print(' estimated skin depth period    {0:.6g} (s)'.format(
                                               skin_depth_period[period_index]))
    print(' estimate apparent resisitivity {0:.0f} (Ohm-m)'.format(
           apparent_res_xy[period_index].mean()))
    print(' estimated depth                {0:.6g} ({1})'.format(depth, ms))
    print(' index                          {0}'.format(period_index))
    print('-'*60)


    return depth, period_index



#==============================================================================
# plot slices
#==============================================================================
class PlotSlices(object):
    """
    plot all slices and be able to scroll through the model

    :Example: ::

        >>> import mtpy.modeling.ws3dinv as ws
        >>> mfn = r"/home/MT/ws3dinv/Inv1/Test_model.00"
        >>> sfn = r"/home/MT/ws3dinv/Inv1/WSStationLocations.txt"
        >>> # plot just first layer to check the formating        
        >>> pds = ws.PlotSlices(model_fn=mfn, station_fn=sfn)

    ======================= ===================================================
    Buttons                  Description    
    ======================= ===================================================
    'e'                     moves n-s slice east by one model block
    'w'                     moves n-s slice west by one model block
    'n'                     moves e-w slice north by one model block
    'm'                     moves e-w slice south by one model block
    'd'                     moves depth slice down by one model block
    'u'                     moves depth slice up by one model block
    ======================= ===================================================


    ======================= ===================================================
    Attributes              Description    
    ======================= ===================================================
    ax_en                   matplotlib.axes instance for depth slice  map view 
    ax_ez                   matplotlib.axes instance for e-w slice
    ax_map                  matplotlib.axes instance for location map
    ax_nz                   matplotlib.axes instance for n-s slice
    climits                 (min , max) color limits on resistivity in log 
                            scale. *default* is (0, 4)
    cmap                    name of color map for resisitiviy.
                            *default* is 'jet_r'
    data_fn                 full path to data file name
    dscale                  scaling parameter depending on map_scale
    east_line_xlist         list of line nodes of east grid for faster plotting
    east_line_ylist         list of line nodes of east grid for faster plotting
    ew_limits               (min, max) limits of e-w in map_scale units
                            *default* is None and scales to station area
    fig                     matplotlib.figure instance for figure
    fig_aspect              aspect ratio of plots. *default* is 1
    fig_dpi                 resolution of figure in dots-per-inch
                            *default* is 300
    fig_num                 figure instance number
    fig_size                [width, height] of figure window. 
                            *default* is [6,6]
    font_dict               dictionary of font keywords, internally created
    font_size               size of ticklables in points, axes labes are 
                            font_size+2. *default* is 7
    grid_east               relative location of grid nodes in e-w direction
                            in map_scale units
    grid_north              relative location of grid nodes in n-s direction
                            in map_scale units
    grid_z                  relative location of grid nodes in z direction
                            in map_scale units
    index_east              index value of grid_east being plotted
    index_north             index value of grid_north being plotted
    index_vertical          index value of grid_z being plotted
    initial_fn              full path to initial file
    key_press               matplotlib.canvas.connect instance
    map_scale               [ 'm' | 'km' ] scale of map. *default* is km
    mesh_east               np.meshgrid(grid_east, grid_north)[0]
    mesh_en_east            np.meshgrid(grid_east, grid_north)[0]
    mesh_en_north           np.meshgrid(grid_east, grid_north)[1]
    mesh_ez_east            np.meshgrid(grid_east, grid_z)[0]
    mesh_ez_vertical        np.meshgrid(grid_east, grid_z)[1]
    mesh_north              np.meshgrid(grid_east, grid_north)[1]
    mesh_nz_north           np.meshgrid(grid_north, grid_z)[0]
    mesh_nz_vertical        np.meshgrid(grid_north, grid_z)[1]
    model_fn                full path to model file
    ms                      size of station markers in points. *default* is 2
    nodes_east              relative distance betwen nodes in e-w direction
                            in map_scale units
    nodes_north             relative distance betwen nodes in n-s direction
                            in map_scale units
    nodes_z                 relative distance betwen nodes in z direction
                            in map_scale units
    north_line_xlist        list of line nodes north grid for faster plotting  
    north_line_ylist        list of line nodes north grid for faster plotting
    ns_limits               (min, max) limits of plots in n-s direction
                            *default* is None, set veiwing area to station area 
    plot_yn                 [ 'y' | 'n' ] 'y' to plot on instantiation
                            *default* is 'y'
    res_model               np.ndarray(n_north, n_east, n_vertical) of 
                            model resistivity values in linear scale           
    station_color           color of station marker. *default* is black
    station_dict_east       location of stations for each east grid row
    station_dict_north      location of stations for each north grid row
    station_east            location of stations in east direction
    station_fn              full path to station file 
    station_font_color      color of station label 
    station_font_pad        padding between station marker and label
    station_font_rotation   angle of station label
    station_font_size       font size of station label
    station_font_weight     weight of font for station label
    station_id              [min, max] index values for station labels
    station_marker          station marker
    station_names           name of stations
    station_north           location of stations in north direction
    subplot_bottom          distance between axes and bottom of figure window
    subplot_hspace          distance between subplots in vertical direction
    subplot_left            distance between axes and left of figure window  
    subplot_right           distance between axes and right of figure window
    subplot_top             distance between axes and top of figure window
    subplot_wspace          distance between subplots in horizontal direction
    title                   title of plot 
    z_limits                (min, max) limits in vertical direction,
    ======================= ===================================================

    """

    def __init__(self, model_fn, data_fn=None, station_fn=None,
                 initial_fn=None, **kwargs):
        self.model_fn = model_fn
        self.data_fn = data_fn
        self.station_fn = station_fn
        self.initial_fn = initial_fn

        self.fig_num = kwargs.pop('fig_num', 1)
        self.fig_size = kwargs.pop('fig_size', [6, 6])
        self.fig_dpi = kwargs.pop('dpi', 300)
        self.fig_aspect = kwargs.pop('fig_aspect', 1)
        self.title = kwargs.pop('title', 'on')
        self.font_size = kwargs.pop('font_size', 7)

        self.subplot_wspace = .20
        self.subplot_hspace = .30
        self.subplot_right = .98
        self.subplot_left = .08
        self.subplot_top = .97
        self.subplot_bottom = .1

        self.index_vertical = kwargs.pop('index_vertical', 0)
        self.index_east = kwargs.pop('index_east', 0)
        self.index_north = kwargs.pop('index_north', 0)

        self.cmap = kwargs.pop('cmap', 'jet_r')
        self.climits = kwargs.pop('climits', (0, 4))

        self.map_scale = kwargs.pop('map_scale', 'km')
        #make map scale
        if self.map_scale=='km':
            self.dscale=1000.
        elif self.map_scale=='m':
            self.dscale=1.
        self.ew_limits = kwargs.pop('ew_limits', None)
        self.ns_limits = kwargs.pop('ns_limits', None)
        self.z_limits = kwargs.pop('z_limits', None)

        self.res_model = None
        self.grid_east = None
        self.grid_north = None
        self.grid_z  = None

        self.nodes_east = None
        self.nodes_north = None
        self.nodes_z = None

        self.mesh_east = None
        self.mesh_north = None

        self.station_east = None
        self.station_north = None
        self.station_names = None

        self.station_id = kwargs.pop('station_id', None)
        self.station_font_size = kwargs.pop('station_font_size', 8)
        self.station_font_pad = kwargs.pop('station_font_pad', 1.0)
        self.station_font_weight = kwargs.pop('station_font_weight', 'bold')
        self.station_font_rotation = kwargs.pop('station_font_rotation', 60)
        self.station_font_color = kwargs.pop('station_font_color', 'k')
        self.station_marker = kwargs.pop('station_marker',
                                         r"$\blacktriangledown$")
        self.station_color = kwargs.pop('station_color', 'k')
        self.ms = kwargs.pop('ms', 10)

        self.plot_yn = kwargs.pop('plot_yn', 'y')
        if self.plot_yn == 'y':
            self.plot()


    def read_files(self):
        """
        read in the files to get appropriate information
        """
        #--> read in model file
        if self.model_fn is not None:
            if os.path.isfile(self.model_fn) == True:
                wsmodel = WSModel(self.model_fn)
                self.res_model = wsmodel.res_model
                self.grid_east = wsmodel.grid_east/self.dscale
                self.grid_north = wsmodel.grid_north/self.dscale
                self.grid_z = wsmodel.grid_z/self.dscale
                self.nodes_east = wsmodel.nodes_east/self.dscale
                self.nodes_north = wsmodel.nodes_north/self.dscale
                self.nodes_z = wsmodel.nodes_z/self.dscale
            else:
                raise mtex.MTpyError_file_handling(
                        '{0} does not exist, check path'.format(self.model_fn))

        #--> read in data file to get station locations
        if self.data_fn is not None:
            if os.path.isfile(self.data_fn) == True:
                wsdata = WSData()
                wsdata.read_data_file(self.data_fn)
                self.station_east = wsdata.data['east']/self.dscale
                self.station_north = wsdata.data['north']/self.dscale
                self.station_names = wsdata.data['station']
            else:
                print('Could not find data file {0}'.format(self.data_fn))

        #--> read in station file
        if self.station_fn is not None:
            if os.path.isfile(self.station_fn) == True:
                wsstations = WSStation(self.station_fn)
                wsstations.read_station_file()
                self.station_east = wsstations.east/self.dscale
                self.station_north = wsstations.north/self.dscale
                self.station_names = wsstations.names
            else:
                print('Could not find station file {0}'.format(self.station_fn))

        #--> read in initial file
        if self.initial_fn is not None:
            if os.path.isfile(self.initial_fn) == True:
                wsmesh = WSMesh()
                wsmesh.read_initial_file(self.initial_fn)
                self.grid_east = wsmesh.grid_east/self.dscale
                self.grid_north = wsmesh.grid_north/self.dscale
                self.grid_z = wsmesh.grid_z/self.dscale
                self.nodes_east = wsmesh.nodes_east/self.dscale
                self.nodes_north = wsmesh.nodes_north/self.dscale
                self.nodes_z = wsmesh.nodes_z/self.dscale

                #need to convert index values to resistivity values
                rdict = dict([(ii,res) for ii,res in enumerate(wsmesh.res_list,1)])

                for ii in range(len(wsmesh.res_list)):
                    self.res_model[np.where(wsmesh.res_model==ii+1)] = \
                                                                    rdict[ii+1]
            else:
                raise mtex.MTpyError_file_handling(
                     '{0} does not exist, check path'.format(self.initial_fn))

        if self.initial_fn is None and self.model_fn is None:
            raise mtex.MTpyError_inputarguments('Need to input either a model'
                                                ' file or initial file.')

    def plot(self):
        """
        plot:
            east vs. vertical,
            north vs. vertical,
            east vs. north


        """

        self.read_files()

        self.get_station_grid_locations()

        self.font_dict = {'size':self.font_size+2, 'weight':'bold'}
        #set the limits of the plot
        if self.ew_limits == None:
            if self.station_east is not None:
                self.ew_limits = (np.floor(self.station_east.min()),
                                  np.ceil(self.station_east.max()))
            else:
                self.ew_limits = (self.grid_east[5], self.grid_east[-5])

        if self.ns_limits == None:
            if self.station_north is not None:
                self.ns_limits = (np.floor(self.station_north.min()),
                                  np.ceil(self.station_north.max()))
            else:
                self.ns_limits = (self.grid_north[5], self.grid_north[-5])

        if self.z_limits == None:
                self.z_limits = (self.grid_z[0]-5000/self.dscale,
                                 self.grid_z[-5])


        self.fig = plt.figure(self.fig_num, figsize=self.fig_size,
                              dpi=self.fig_dpi)
        plt.clf()
        gs = gridspec.GridSpec(2, 2,
                               wspace=self.subplot_wspace,
                               left=self.subplot_left,
                               top=self.subplot_top,
                               bottom=self.subplot_bottom,
                               right=self.subplot_right,
                               hspace=self.subplot_hspace)

        #make subplots
        self.ax_ez = self.fig.add_subplot(gs[0, 0], aspect=self.fig_aspect)
        self.ax_nz = self.fig.add_subplot(gs[1, 1], aspect=self.fig_aspect)
        self.ax_en = self.fig.add_subplot(gs[1, 0], aspect=self.fig_aspect)
        self.ax_map = self.fig.add_subplot(gs[0, 1])

        #make grid meshes being sure the indexing is correct
        self.mesh_ez_east, self.mesh_ez_vertical = np.meshgrid(self.grid_east,
                                                               self.grid_z,
                                                               indexing='ij')
        self.mesh_nz_north, self.mesh_nz_vertical = np.meshgrid(self.grid_north,
                                                                self.grid_z,
                                                                indexing='ij')
        self.mesh_en_east, self.mesh_en_north = np.meshgrid(self.grid_east,
                                                            self.grid_north,
                                                            indexing='ij')

        #--> plot east vs vertical
        self._update_ax_ez()

        #--> plot north vs vertical
        self._update_ax_nz()

        #--> plot east vs north
        self._update_ax_en()

        #--> plot the grid as a map view
        self._update_map()

        #plot color bar
        cbx = mcb.make_axes(self.ax_map, fraction=.15, shrink=.75, pad = .1)
        cb = mcb.ColorbarBase(cbx[0],
                              cmap=self.cmap,
                              norm=Normalize(vmin=self.climits[0],
                                             vmax=self.climits[1]))


        cb.ax.yaxis.set_label_position('right')
        cb.ax.yaxis.set_label_coords(1.25,.5)
        cb.ax.yaxis.tick_left()
        cb.ax.tick_params(axis='y',direction='in')

        cb.set_label('Resistivity ($\Omega \cdot$m)',
                     fontdict={'size':self.font_size+1})

        cb.set_ticks(np.arange(np.ceil(self.climits[0]),
                               np.floor(self.climits[1]+1)))
        cblabeldict={-2:'$10^{-3}$',-1:'$10^{-1}$',0:'$10^{0}$',1:'$10^{1}$',
                     2:'$10^{2}$',3:'$10^{3}$',4:'$10^{4}$',5:'$10^{5}$',
                     6:'$10^{6}$',7:'$10^{7}$',8:'$10^{8}$'}
        cb.set_ticklabels([cblabeldict[cc]
                            for cc in np.arange(np.ceil(self.climits[0]),
                                                np.floor(self.climits[1]+1))])

        plt.show()

        self.key_press = self.fig.canvas.mpl_connect('key_press_event',
                                                     self.on_key_press)


    def on_key_press(self, event):
        """
        on a key press change the slices

        """

        key_press = event.key

        if key_press == 'n':
            if self.index_north == self.grid_north.shape[0]:
                print('Already at northern most grid cell')
            else:
                self.index_north += 1
                if self.index_north > self.grid_north.shape[0]:
                    self.index_north = self.grid_north.shape[0]
            self._update_ax_ez()
            self._update_map()

        if key_press == 'm':
            if self.index_north == 0:
                print('Already at southern most grid cell')
            else:
                self.index_north -= 1
                if self.index_north < 0:
                    self.index_north = 0
            self._update_ax_ez()
            self._update_map()

        if key_press == 'e':
            if self.index_east == self.grid_east.shape[0]:
                print('Already at eastern most grid cell')
            else:
                self.index_east += 1
                if self.index_east > self.grid_east.shape[0]:
                    self.index_east = self.grid_east.shape[0]
            self._update_ax_nz()
            self._update_map()

        if key_press == 'w':
            if self.index_east == 0:
                print('Already at western most grid cell')
            else:
                self.index_east -= 1
                if self.index_east < 0:
                    self.index_east = 0
            self._update_ax_nz()
            self._update_map()

        if key_press == 'd':
            if self.index_vertical == self.grid_z.shape[0]:
                print('Already at deepest grid cell')
            else:
                self.index_vertical += 1
                if self.index_vertical > self.grid_z.shape[0]:
                    self.index_vertical = self.grid_z.shape[0]
            self._update_ax_en()
            print('Depth = {0:.5g} ({1})'.format(self.grid_z[self.index_vertical],
                                                 self.map_scale))

        if key_press == 'u':
            if self.index_vertical == 0:
                print('Already at surface grid cell')
            else:
                self.index_vertical -= 1
                if self.index_vertical < 0:
                    self.index_vertical = 0
            self._update_ax_en()
            print('Depth = {0:.5gf} ({1})'.format(self.grid_z[self.index_vertical],
                                                 self.map_scale))

    def _update_ax_ez(self):
        """
        update east vs vertical plot
        """
        self.ax_ez.cla()
        plot_ez = np.log10(self.res_model[self.index_north, :, :])
        self.ax_ez.pcolormesh(self.mesh_ez_east,
                              self.mesh_ez_vertical,
                              plot_ez,
                              cmap=self.cmap,
                              vmin=self.climits[0],
                              vmax=self.climits[1])
        #plot stations
        for sx in self.station_dict_north[self.grid_north[self.index_north]]:
            self.ax_ez.text(sx,
                            0,
                            self.station_marker,
                            horizontalalignment='center',
                            verticalalignment='baseline',
                            fontdict={'size':self.ms,
                                      'color':self.station_color})

        self.ax_ez.set_xlim(self.ew_limits)
        self.ax_ez.set_ylim(self.z_limits[1], self.z_limits[0])
        self.ax_ez.set_ylabel('Depth ({0})'.format(self.map_scale),
                              fontdict=self.font_dict)
        self.ax_ez.set_xlabel('Easting ({0})'.format(self.map_scale),
                              fontdict=self.font_dict)
        self.fig.canvas.draw()
        self._update_map()

    def _update_ax_nz(self):
        """
        update east vs vertical plot
        """
        self.ax_nz.cla()
        plot_nz = np.log10(self.res_model[:, self.index_east, :])
        self.ax_nz.pcolormesh(self.mesh_nz_north,
                              self.mesh_nz_vertical,
                              plot_nz,
                              cmap=self.cmap,
                              vmin=self.climits[0],
                              vmax=self.climits[1])
        #plot stations
        for sy in self.station_dict_east[self.grid_east[self.index_east]]:
            self.ax_nz.text(sy,
                            0,
                            self.station_marker,
                            horizontalalignment='center',
                            verticalalignment='baseline',
                            fontdict={'size':self.ms,
                                      'color':self.station_color})
        self.ax_nz.set_xlim(self.ns_limits)
        self.ax_nz.set_ylim(self.z_limits[1], self.z_limits[0])
        self.ax_nz.set_xlabel('Northing ({0})'.format(self.map_scale),
                              fontdict=self.font_dict)
        self.ax_nz.set_ylabel('Depth ({0})'.format(self.map_scale),
                              fontdict=self.font_dict)
        self.fig.canvas.draw()
        self._update_map()

    def _update_ax_en(self):
        """
        update east vs vertical plot
        """

        self.ax_en.cla()
        plot_en = np.log10(self.res_model[:, :, self.index_vertical].T)
        self.ax_en.pcolormesh(self.mesh_en_east,
                              self.mesh_en_north,
                              plot_en,
                              cmap=self.cmap,
                              vmin=self.climits[0],
                              vmax=self.climits[1])
        self.ax_en.set_xlim(self.ew_limits)
        self.ax_en.set_ylim(self.ns_limits)
        self.ax_en.set_ylabel('Northing ({0})'.format(self.map_scale),
                              fontdict=self.font_dict)
        self.ax_en.set_xlabel('Easting ({0})'.format(self.map_scale),
                              fontdict=self.font_dict)
        #--> plot the stations
        if self.station_east is not None:
            for ee, nn in zip(self.station_east, self.station_north):
                self.ax_en.text(ee, nn, '*',
                                 verticalalignment='center',
                                 horizontalalignment='center',
                                 fontdict={'size':5, 'weight':'bold'})

        self.fig.canvas.draw()
        self._update_map()

    def _update_map(self):
        self.ax_map.cla()
        self.east_line_xlist = []
        self.east_line_ylist = []
        for xx in self.grid_east:
            self.east_line_xlist.extend([xx, xx])
            self.east_line_xlist.append(None)
            self.east_line_ylist.extend([self.grid_north.min(),
                                         self.grid_north.max()])
            self.east_line_ylist.append(None)
        self.ax_map.plot(self.east_line_xlist,
                         self.east_line_ylist,
                         lw=.25,
                         color='k')

        self.north_line_xlist = []
        self.north_line_ylist = []
        for yy in self.grid_north:
            self.north_line_xlist.extend([self.grid_east.min(),
                                          self.grid_east.max()])
            self.north_line_xlist.append(None)
            self.north_line_ylist.extend([yy, yy])
            self.north_line_ylist.append(None)
        self.ax_map.plot(self.north_line_xlist,
                         self.north_line_ylist,
                         lw=.25,
                         color='k')
        #--> e-w indication line
        self.ax_map.plot([self.grid_east.min(),
                          self.grid_east.max()],
                         [self.grid_north[self.index_north],
                          self.grid_north[self.index_north]],
                         lw=1,
                         color='g')

        #--> e-w indication line
        self.ax_map.plot([self.grid_east[self.index_east],
                          self.grid_east[self.index_east]],
                         [self.grid_north.min(),
                          self.grid_north.max()],
                         lw=1,
                         color='b')
         #--> plot the stations
        if self.station_east is not None:
            for ee, nn in zip(self.station_east, self.station_north):
                self.ax_map.text(ee, nn, '*',
                                 verticalalignment='center',
                                 horizontalalignment='center',
                                 fontdict={'size':5, 'weight':'bold'})

        self.ax_map.set_xlim(self.ew_limits)
        self.ax_map.set_ylim(self.ns_limits)
        self.ax_map.set_ylabel('Northing ({0})'.format(self.map_scale),
                              fontdict=self.font_dict)
        self.ax_map.set_xlabel('Easting ({0})'.format(self.map_scale),
                              fontdict=self.font_dict)

        #plot stations
        self.ax_map.text(self.ew_limits[0]*.95, self.ns_limits[1]*.95,
                 '{0:.5g} ({1})'.format(self.grid_z[self.index_vertical],
                                        self.map_scale),
                 horizontalalignment='left',
                 verticalalignment='top',
                 bbox={'facecolor': 'white'},
                 fontdict=self.font_dict)


        self.fig.canvas.draw()

    def get_station_grid_locations(self):
        """
        get the grid line on which a station resides for plotting

        """
        self.station_dict_east = dict([(gx, []) for gx in self.grid_east])
        self.station_dict_north = dict([(gy, []) for gy in self.grid_north])
        if self.station_east is not None:
            for ss, sx in enumerate(self.station_east):
                gx = np.where(self.grid_east <= sx)[0][-1]
                self.station_dict_east[self.grid_east[gx]].append(self.station_north[ss])

            for ss, sy in enumerate(self.station_north):
                gy = np.where(self.grid_north <= sy)[0][-1]
                self.station_dict_north[self.grid_north[gy]].append(self.station_east[ss])
        else:
            return

    def redraw_plot(self):
        """
        redraw plot if parameters were changed

        use this function if you updated some attributes and want to re-plot.

        :Example: ::

            >>> # change the color and marker of the xy components
            >>> import mtpy.modeling.occam2d as occam2d
            >>> ocd = occam2d.Occam2DData(r"/home/occam2d/Data.dat")
            >>> p1 = ocd.plotAllResponses()
            >>> #change line width
            >>> p1.lw = 2
            >>> p1.redraw_plot()
        """

        plt.close(self.fig)
        self.plot()

    def save_figure(self, save_fn=None, fig_dpi=None, file_format='pdf',
                    orientation='landscape', close_fig='y'):
        """
        save_figure will save the figure to save_fn.

        Arguments:
        -----------

            **save_fn** : string
                          full path to save figure to, can be input as
                          * directory path -> the directory path to save to
                            in which the file will be saved as 
                            save_fn/station_name_PhaseTensor.file_format

                          * full path -> file will be save to the given 
                            path.  If you use this option then the format
                            will be assumed to be provided by the path

            **file_format** : [ pdf | eps | jpg | png | svg ]
                              file type of saved figure pdf,svg,eps... 

            **orientation** : [ landscape | portrait ]
                              orientation in which the file will be saved
                              *default* is portrait

            **fig_dpi** : int
                          The resolution in dots-per-inch the file will be
                          saved.  If None then the dpi will be that at 
                          which the figure was made.  I don't think that 
                          it can be larger than dpi of the figure.

            **close_plot** : [ y | n ]
                             * 'y' will close the plot after saving.
                             * 'n' will leave plot open

        :Example: ::

            >>> # to save plot as jpg
            >>> import mtpy.modeling.occam2d as occam2d
            >>> dfn = r"/home/occam2d/Inv1/data.dat"
            >>> ocd = occam2d.Occam2DData(dfn)
            >>> ps1 = ocd.plotPseudoSection()
            >>> ps1.save_plot(r'/home/MT/figures', file_format='jpg')

        """

        if fig_dpi == None:
            fig_dpi = self.fig_dpi

        if os.path.isdir(save_fn) == False:
            file_format = save_fn[-3:]
            self.fig.savefig(save_fn, dpi=fig_dpi, format=file_format,
                             orientation=orientation, bbox_inches='tight')

        else:
            save_fn = os.path.join(save_fn, '_E{0}_N{1}_Z{2}.{3}'.format(
                                    self.index_east, self.index_north,
                                    self.index_vertical, file_format))
            self.fig.savefig(save_fn, dpi=fig_dpi, format=file_format,
                        orientation=orientation, bbox_inches='tight')

        if close_fig == 'y':
            plt.clf()
            plt.close(self.fig)

        else:
            pass

        self.fig_fn = save_fn
        print('Saved figure to: '+self.fig_fn)

#==============================================================================
# STARTUP FILES
#==============================================================================
class WSStartup(object):
    """
    read and write startup files

    :Example: ::

        >>> import mtpy.modeling.ws3dinv as ws
        >>> dfn = r"/home/MT/ws3dinv/Inv1/WSDataFile.dat"
        >>> ifn = r"/home/MT/ws3dinv/Inv1/init3d"
        >>> sws = ws.WSStartup(data_fn=dfn, initial_fn=ifn)


    =================== =======================================================
    Attributes          Description    
    =================== =======================================================
    apriori_fn          full path to *a priori* model file
                        *default* is 'default'
    control_fn          full path to model index control file
                        *default* is 'default'
    data_fn             full path to data file
    error_tol           error tolerance level
                        *default* is 'default'
    initial_fn          full path to initial model file
    lagrange            starting lagrange multiplier
                        *default* is 'default'
    max_iter            max number of iterations
                        *default* is 10
    model_ls            model length scale
                        *default* is 5 0.3 0.3 0.3
    output_stem         output file name stem
                        *default* is 'ws3dinv'
    save_path           directory to save file to
    startup_fn          full path to startup file
    static_fn           full path to statics file
                        *default* is 'default'
    target_rms          target rms
                        *default* is 1.0
    =================== =======================================================

    """

    def __init__(self, data_fn=None, initial_fn=None, **kwargs):
        self.data_fn = data_fn
        self.initial_fn = initial_fn

        self.output_stem = kwargs.pop('output_stem', 'ws3dinv')
        self.apriori_fn = kwargs.pop('apriori_fn', 'default')
        self.model_ls = kwargs.pop('model_ls', [5, 0.3, 0.3, 0.3])
        self.target_rms = kwargs.pop('target_rms', 1.0)
        self.control_fn = kwargs.pop('control_fn', 'default')
        self.max_iter = kwargs.pop('max_iter', 10)
        self.error_tol = kwargs.pop('error_tol', 'default')
        self.static_fn = kwargs.pop('static_fn', 'default')
        self.lagrange = kwargs.pop('lagrange', 'default')
        self.save_path = kwargs.pop('save_path', None)

        self.startup_fn = kwargs.pop('startup_fn', None)

        self._startup_keys = ['data_file',
                               'output_file',
                               'initial_model_file',
                               'prior_model_file',
                               'control_model_index',
                               'target_rms',
                               'max_no_iteration',
                               'model_length_scale',
                               'lagrange_info',
                               'error_tol_level',
                               'static_file']



    def write_startup_file(self):
        """
        makes a startup file for WSINV3D.  

        """
        if self.data_fn is None:
            raise IOError('Need to input data file name')

        if self.initial_fn is None:
            raise IOError('Need to input initial model file name')

        #create the output filename
        if self.save_path == None and self.data_fn != None:
            self.startup_fn = os.path.join(os.path.dirname(self.data_fn),
                                           'startup')
        elif os.path.isdir(self.save_path) == True:
            self.startup_fn = os.path.join(self.save_path, 'startup')
        else:
            self.startup_fn = self.save_path

        slines = []

        if os.path.dirname(self.startup_fn) == os.path.dirname(self.data_fn):
            slines.append('{0:<20}{1}\n'.format('DATA_FILE',
                                                os.path.basename(self.data_fn)))
            if len(os.path.basename(self.data_fn)) > 70:
                print('Data file is too long, going to get an error at runtime')
        else:
            slines.append('{0:<20}{1}\n'.format('DATA_FILE',self.data_fn))
            if len(self.data_fn) > 70:
                print('Data file is too long, going to get an error at runtime')

        slines.append('{0:<20}{1}\n'.format('OUTPUT_FILE', self.output_stem))

        if os.path.dirname(self.startup_fn) == os.path.dirname(self.initial_fn):
            slines.append('{0:<20}{1}\n'.format('INITIAL_MODEL_FILE',
                                             os.path.basename(self.initial_fn)))
        else:
            slines.append('{0:<20}{1}\n'.format('INITIAL_MODEL_FILE',
                                                self.initial_fn))
        slines.append('{0:<20}{1}\n'.format('PRIOR_MODEL_FILE',
                                            self.apriori_fn))
        slines.append('{0:<20}{1}\n'.format('CONTROL_MODEL_INDEX ',
                                            self.control_fn))
        slines.append('{0:<20}{1}\n'.format('TARGET_RMS', self.target_rms))
        slines.append('{0:<20}{1}\n'.format('MAX_NO_ITERATION',
                                            self.max_iter))
        slines.append('{0:<20}{1:.0f} {2:.1f} {3:.1f} {4:.1f}\n'.format(
                                                        'MODEL_LENGTH_SCALE',
                                                        self.model_ls[0],
                                                        self.model_ls[1],
                                                        self.model_ls[2],
                                                        self.model_ls[3]))

        slines.append('{0:<20}{1} \n'.format('LAGRANGE_INFO', self.lagrange))
        slines.append('{0:<20}{1} \n'.format('ERROR_TOL_LEVEL',
                                             self.error_tol))
        slines.append('{0:<20}{1} \n'.format('STATIC_FILE', self.static_fn))

        sfid = file(self.startup_fn, 'w')
        sfid.write(''.join(slines))
        sfid.close()

        print('Wrote startup file to: {0}'.format(self.startup_fn))

    def read_startup_file(self, startup_fn=None):
        """
        read startup file fills attributes

        """
        if startup_fn is not None:
            self.startup_fn = startup_fn
        if self.startup_fn is None:
            raise IOError('Need to input startup file name')

        self.save_path = os.path.dirname(self.startup_fn)

        sfid = file(self.startup_fn, 'r')
        slines = sfid.readlines()
        sfid.close()

        slines = [ss.strip().split()[1:] for ss in slines]

        self.data_fn = slines[0][0].strip()
        if self.data_fn.find(os.path.sep) == -1:
            self.data_fn = os.path.join(self.save_path, self.data_fn)
        self.output_stem = slines[1][0].strip()
        self.initial_fn = slines[2][0].strip()
        if self.initial_fn.find(os.path.sep) == -1:
            self.initial_fn = os.path.join(self.save_path, self.initial_fn)
        self.apriori_fn = slines[3][0].strip()
        self.control_fn = slines[4][0].strip()
        self.target_rms = float(slines[5][0].strip())
        self.max_iter = int(slines[6][0].strip())
        try:
            self.model_ls = [int(slines[7][0]), float(slines[7][1]),
                             float(slines[7][2]), float(slines[7][3])]
        except ValueError:
            self.model_ls = slines[7][0]

        self.lagrange = slines[8][0].strip()
        self.error_tol = slines[9][0].strip()
        try:
            self.static_fn = slines[10][0].strip()
        except IndexError:
            print('Did not find static_fn')






#==============================================================================
# WRITE A VTK FILE TO IMAGE IN PARAVIEW OR MAYAVI
#==============================================================================

def write_vtk_res_model(res_model, grid_north, grid_east, grid_z, save_fn):
    """
    Write a vtk file for resistivity as a structured grid 
    to be read into paraview or mayavi

    **Doesn't work properly under windows**

    adds extension automatically

    """

    if os.path.isdir(save_fn) == True:
        save_fn = os.path.join(save_fn, 'VTKResistivity_Model')

    save_fn = gridToVTK(save_fn, grid_north, grid_east, grid_z,
                        cellData={'resistivity':res_model})

    return save_fn

def write_vtk_stations(station_north, station_east, save_fn, station_z=None):
    """
    Write a vtk file as points to be read into paraview or mayavi

    **Doesn't work properly under windows**

    adds extension automatically

    """

    if os.path.isdir(save_fn) == True:
        save_fn = os.path.join(save_fn, 'VTKStations')

    if station_z is None:
        station_z = np.zeros_like(station_north)

    pointsToVTK(save_fn, station_north, station_east, station_z,
              cellData={'value':np.ones_like(station_north)})

    return save_fn

def write_vtk_files(model_fn, station_fn, save_path):
    """
    writes vtk files
    """

    wsstation = WSStation(station_fn=station_fn)
    wsstation.read_station_file()
    wsstation.write_vtk_file(save_path)

    wsmodel = WSModel(model_fn)
    wsmodel.write_vtk_file(save_path)


def computeMemoryUsage(nx, ny, nz, n_stations, n_zelements, n_period):
    """
    compute the memory usage of a model

    Arguments:
    ----------
        **nx** : int
                 number of cells in N-S direction

        **ny** : int
                 number of cells in E-W direction

        **nz** : int
                 number of cells in vertical direction including air layers (7)

        **n_stations** : int
                         number of stations

        **n_zelements** : int
                          number of impedence tensor elements either 4 or 8

        **n_period** : int
                       number of periods to invert for

    Returns:
    --------
        **mem_req** : float
                      approximate memory useage in GB
    """

    mem_req = 1.2*(8*(n_stations*n_period*n_zelements)**2+
                   8*(nx*ny*nz*n_stations*n_period*n_zelements))
    return mem_req*1E-9





