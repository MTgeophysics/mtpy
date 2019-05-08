# -*- coding: utf-8 -*-
"""
Created on Tue Aug 23 11:54:38 2011

@author: a1185872
"""

import numpy as np
import scipy as sp
import os
import subprocess
import shutil
import fnmatch
from operator import itemgetter
import time
import matplotlib.colorbar as mcb
from matplotlib.colors import Normalize
from matplotlib.ticker import MultipleLocator
import matplotlib.gridspec as gridspec
import mtpy.core.edi as mtedi
import mtpy.modeling.winglinktools as wlt
import matplotlib.pyplot as plt

import mtpy.utils.gis_tools


occamdict = {'1': 'resxy', '2': 'phasexy', '3': 'realtip', '4': 'imagtip', '5': 'resyx',
             '6': 'phaseyx'}


class Occam1D:
    """
    ==============================================
    This class will deal with everything occam 1D
    =============================================
    """

    def __init__(self, savepath=None):

        self.savepath = None
        self.modelfn = None
        self.inputfn = None
        self.datafn_te = None
        self.datafn_tm = None
        self.iternum = 0

    def make1DdataFile(self, station, edipath=None, savepath=None,
                       polarization='both', reserr='data', phaseerr='data',
                       string_fmt='%+.6e', ss=3 * ' ', thetar=0):
        """
        make1Ddatafile will write a data file for Occam1D

        Arguments:
        ---------    
            **station** : string
                          the station name and path if edipath=None

            **edipath** : string
                          path to the edi files to be written into a data file,
                          useful for multile data files

            **savepath** : string
                           path to save the file, if None set to dirname of 
                           station if edipath = None.  Otherwise set to 
                           dirname of edipath.

            **thetar** : float
                         rotation angle to rotate Z. Clockwise positive and N=0
                         *default* = 0

            **polarization** : [ 'both' | 'TE' | 'TM' | 'det']
                              polarization to model can be (*default*='both'):

                                - 'both' for TE and TM as separate files
                                - 'TE' for just TE mode
                                - 'TM' for just TM mode
                                - 'det' for the determinant of Z.
                                .. note:: 

                                    if polarization = 'det' two files 
                                    will be created stationDet_TE.dat and 
                                    stationDet_TM.dat.  These files both use
                                    the determinant however the code for the 
                                    Occam input is different so you can test 
                                    the difference between inverting the 
                                    determinant as TE or TM because there is
                                    no option in Occam for the Det.


            **reserr** : float
                        errorbar for resistivity values.  Can be set to (
                        *default* = 'data'): 

                        - 'data' for errorbars from the data
                        - percent number ex. 10 for ten percent

            **phaseerr** : float
                          errorbar for phase values.  Can be set to (
                          *default* = 'data'):

                            - 'data' for errorbars from the data
                            - percent number ex. 10 for ten percent

            **string_fmt** : format of the values written to the file. 
                      *default* = %+.6e

            **ss** : spacing between values in file.  *default* = ' '*3

        Returns:
        --------
            **Occam1D.datafn_te** : full path to data file for TE mode

            **Occam1D.datafn_tm** : full path to data file for TM mode

        :Example: ::

            >>> old = occam.Occam1D()
            >>> old.make1DdataFile('MT01',edipath=r"/home/Line1",
            >>>                    savepath=r"/home/Occam1D/Line1/Inv1_TE",
            >>>                    mode='TE')
            >>> Wrote Data File: /home/Occam1D/Line1/Inv1_TE/MT01TE.dat 
        """

        if os.path.dirname(station) == '':
            if edipath == None:
                raise IOError('Need to input a path for the file.')
            else:
                # find the edifile
                for fn in os.listdir(edipath):
                    if fn.lower().find(station.lower()) >= 0:
                        edifile = os.path.join(edipath, fn)
        else:
            edifile = station

        self.station = os.path.basename(edifile)[:-4]

        # raise an error if can't find the edifile
        if edifile == None:
            raise NameError('No edifile exists, check path and station name')

        #read in edifile
        impz = Z.Z(edifile)

        # make sure the savepath exists, if not create it
        if savepath == None:
            if not self.savepath:
                savepath = os.path.dirname(edifile)
                if not os.path.exists(savepath):
                    os.mkdir(savepath)
                savepath = self.savepath
            else:
                savepath = self.savepath
        elif os.path.basename(savepath).find('.') > 0:
            savepath = os.path.dirname(savepath)
            if not os.path.exists(savepath):
                os.mkdir(os.path.dirname(savepath))
            self.savepath = savepath
        else:
            if not os.path.exists(savepath):
                os.mkdir(savepath)
            self.savepath = savepath

        # load the edifile and get resistivity and phase
        rp = impz.getResPhase(thetar=thetar)
        freq = impz.frequency
        nf = len(freq)
        returnfn = []

        pdict = {'TE': ['xy'],
                 'TM': ['yx'],
                 'both': ['xy', 'yx'],
                 'det': ['det'],
                 'all': ['xy', 'yx', 'det']}
        if polarization == 'both' or polarization == 'det':
            for pol in ['xy', 'yx']:
                if pol == 'xy':
                    if polarization == 'det':
                        dfilesave = os.path.join(self.savepath,
                                                 impz.station + 'Det_TE.dat')
                    else:
                        dfilesave = os.path.join(self.savepath,
                                                 impz.station + 'TE.dat')
                elif pol == 'yx':
                    if polarization == 'det':
                        dfilesave = os.path.join(self.savepath,
                                                 impz.station + 'Det_TM.dat')
                    else:
                        dfilesave = os.path.join(self.savepath,
                                                 impz.station + 'TM.dat')

                datafid = open(dfilesave, 'w')

                datafid.write('Format:  EMData_1.1 \n')
                datafid.write('!Polarization:' + ss + pol + '\n')

                # needs a transmitter to work so put in a dummy one
                datafid.write('# Transmitters: 1\n')
                datafid.write('0 0 0 0 0 \n')

                # write frequencies
                datafid.write('# Frequencies:' + ss + str(nf) + '\n')
                for ff in freq:
                    datafid.write(ss + '%.6f' % ff + '\n')

                # needs a receiver to work so put in a dummy one
                datafid.write('# Receivers: 1 \n')
                datafid.write('0 0 0 0 0 0 \n')

                # write data
                datafid.write('# Data:' + 2 * ss + str(2 * nf) + '\n')
                datafid.write('!' + 2 * ss + 'Type' + 2 * ss + 'Freq#' + 2 * ss + 'Tx#' + 2 * ss +
                              'Rx#' + 2 * ss + 'Data' + 2 * ss + 'Std_Error' + '\n')

                # put the yx phase component in the first quadrant as
                # prescribed
                if pol == 'yx':
                    rp.phaseyx = rp.phaseyx + 180
                    # check if there are any negative phases
                    negphase = np.where(rp.phaseyx > 180)
                    if len(negphase) > 0:
                        rp.phaseyx[negphase[0]] = rp.phaseyx\
                            [negphase[0]] - 360

                # write the resistivity and phase components
                for ii in range(nf):
                    #------------write resistivity components------------------
                    if reserr == 'data':
                        if pol == 'xy':
                            if polarization == 'det':
                                datafid.write(2 * ss + 'RhoZ' + pol + 2 * ss + str(ii + 1) +
                                              2 * ss + '0' + 2 * ss + '0' + 2 * ss +
                                              string_fmt % rp.resdet[ii] + 2 * ss +
                                              string_fmt % rp.resdeterr[ii] + '\n')
                            else:
                                datafid.write(2 * ss + 'RhoZ' + pol + 2 * ss + str(ii + 1) +
                                              2 * ss + '0' + 2 * ss + '0' + 2 * ss +
                                              string_fmt % rp.resxy[ii] + 2 * ss +
                                              string_fmt % rp.resxyerr[ii] + '\n')
                        elif pol == 'yx':
                            if polarization == 'det':
                                datafid.write(2 * ss + 'RhoZ' + pol + 2 * ss + str(ii + 1) +
                                              2 * ss + '0' + 2 * ss + '0' + 2 * ss +
                                              string_fmt % rp.resdet[ii] + 2 * ss +
                                              string_fmt % rp.resdeterr[ii] + '\n')
                            else:
                                datafid.write(2 * ss + 'RhoZ' + pol + 2 * ss + str(ii + 1) +
                                              2 * ss + '0' + 2 * ss + '0' + 2 * ss +
                                              string_fmt % rp.resyx[ii] + 2 * ss +
                                              string_fmt % rp.resyxerr[ii] + '\n')
                    #-----------if percent error is given--------------------
                    else:
                        if pol == 'xy':
                            if polarization == 'det':
                                datafid.write(2 * ss + 'RhoZ' + pol + 2 * ss + str(ii + 1) +
                                              2 * ss + '0' + 2 * ss + '0' + 2 * ss +
                                              string_fmt % rp.resdet[ii] + 2 * ss +
                                              string_fmt % (rp.resdet[ii] * reserr / 100.) +
                                              '\n')
                            else:
                                datafid.write(2 * ss + 'RhoZ' + pol + 2 * ss + str(ii + 1) +
                                              2 * ss + '0' + 2 * ss + '0' + 2 * ss +
                                              string_fmt % rp.resxy[ii] + 2 * ss +
                                              string_fmt % (rp.resxy[ii] * reserr / 100.) +
                                              '\n')
                        elif pol == 'yx':
                            if polarization == 'det':
                                datafid.write(2 * ss + 'RhoZ' + pol + 2 * ss + str(ii + 1) +
                                              2 * ss + '0' + 2 * ss + '0' + 2 * ss +
                                              string_fmt % rp.resdet[ii] + 2 * ss +
                                              string_fmt % (rp.resdet[ii] * reserr / 100.) +
                                              '\n')
                            else:
                                datafid.write(2 * ss + 'RhoZ' + pol + 2 * ss + str(ii + 1) +
                                              2 * ss + '0' + 2 * ss + '0' + 2 * ss +
                                              string_fmt % rp.resyx[ii] + 2 * ss +
                                              string_fmt % (rp.resyx[ii] * reserr / 100.) +
                                              '\n')

                    #---------------write phase components--------------------
                    if phaseerr == 'data':
                        if pol == 'xy':
                            if polarization == 'det':
                                datafid.write(2 * ss + 'PhsZ' + pol + 2 * ss + str(ii + 1) +
                                              2 * ss + '0' + 2 * ss + '0' + 2 * ss +
                                              string_fmt % rp.phasedet[ii] + 2 * ss +
                                              string_fmt % rp.phasedeterr[ii] + '\n')
                            else:
                                datafid.write(2 * ss + 'PhsZ' + pol + 2 * ss + str(ii + 1) +
                                              2 * ss + '0' + 2 * ss + '0' + 2 * ss +
                                              string_fmt % rp.phasexy[ii] + 2 * ss +
                                              string_fmt % rp.phasexyerr[ii] + '\n')
                        if pol == 'yx':
                            if polarization == 'det':
                                datafid.write(2 * ss + 'PhsZ' + pol + 2 * ss + str(ii + 1) +
                                              2 * ss + '0' + 2 * ss + '0' + 2 * ss +
                                              string_fmt % rp.phasedet[ii] + 2 * ss +
                                              string_fmt % rp.phasedeterr[ii] + '\n')
                            else:
                                datafid.write(2 * ss + 'PhsZ' + pol + 2 * ss + str(ii + 1) +
                                              2 * ss + '0' + 2 * ss + '0' + 2 * ss +
                                              string_fmt % rp.phasexy[ii] + 2 * ss +
                                              string_fmt % rp.phasexyerr[ii] + '\n')
                    #-----------if percent error is given--------------------
                    else:
                        if pol == 'xy':
                            if polarization == 'det':
                                datafid.write(2 * ss + 'PhsZ' + pol + 2 * ss + str(ii + 1) +
                                              2 * ss + '0' + 2 * ss + '0' + 2 * ss +
                                              string_fmt % rp.phasedet[ii] + 2 * ss +
                                              string_fmt % (phaseerr / 100. * (180 / np.pi)) +
                                              '\n')
                            else:
                                datafid.write(2 * ss + 'PhsZ' + pol + 2 * ss + str(ii + 1) +
                                              2 * ss + '0' + 2 * ss + '0' + 2 * ss +
                                              string_fmt % rp.phasexy[ii] + 2 * ss +
                                              string_fmt % (phaseerr / 100. * (180 / np.pi)) +
                                              '\n')
                        if pol == 'yx':
                            if polarization == 'det':
                                datafid.write(2 * ss + 'PhsZ' + pol + 2 * ss + str(ii + 1) +
                                              2 * ss + '0' + 2 * ss + '0' + 2 * ss +
                                              string_fmt % rp.phasedet[ii] + 2 * ss +
                                              string_fmt % (phaseerr / 100. * (180 / np.pi)) +
                                              '\n')
                            else:
                                datafid.write(2 * ss + 'PhsZ' + pol + 2 * ss + str(ii + 1) +
                                              2 * ss + '0' + 2 * ss + '0' + 2 * ss +
                                              string_fmt % rp.phaseyx[ii] + 2 * ss +
                                              string_fmt % (phaseerr / 100. * (180 / np.pi)) +
                                              '\n')
                datafid.write('\n')
                datafid.close()
                print('Wrote Data File: ', dfilesave)
                returnfn.append(dfilesave)
            self.datafn_te = returnfn[0]
            self.datafn_tm = returnfn[1]
        else:
            if polarization == 'TE':
                pol = 'xy'
                dfilesave = os.path.join(savepath, impz.station + 'TE.dat')
                self.datafn_te = dfilesave
            elif polarization == 'TM':
                pol = 'yx'
                dfilesave = os.path.join(savepath, impz.station + 'TM.dat')
                self.datafn_te = dfilesave

            # open file to write to
            datafid = open(dfilesave, 'w')
            datafid.write('Format:  EMData_1.1 \n')
            datafid.write('!Polarization:' + ss + pol + '\n')

            # needs a transmitter to work so put in a dummy one
            datafid.write('# Transmitters: 1\n')
            datafid.write('0 0 0 0 0 \n')

            # write frequencies
            datafid.write('# Frequencies:' + ss + str(nf) + '\n')
            for ff in freq:
                datafid.write(ss + '%.6f' % ff + '\n')

            # needs a receiver to work so put in a dummy one
            datafid.write('# Receivers: 1 \n')
            datafid.write('0 0 0 0 0 0 \n')

            # write header line
            datafid.write('# Data:' + 2 * ss + str(2 * nf) + '\n')
            datafid.write('!' + 2 * ss + 'Type' + 2 * ss + 'Freq#' + 2 * ss + 'Tx#' + 2 * ss + 'Rx#' +
                          2 * ss + 'Data' + 2 * ss + 'Std_Error' + '\n')

            # put the yx phase component in the first quadrant as prescribed
            if pol == 'yx':
                rp.phaseyx = rp.phaseyx + 180
                # check if there are any negative phases
                negphase = np.where(rp.phaseyx > 180)
                if len(negphase) > 0:
                    rp.phaseyx[negphase[0]] = rp.phaseyx\
                        [negphase[0]] - 360

            # write the resistivity and phase components
            for ii in range(nf):
                # write resistivity components
                if reserr == 'data':
                    if pol == 'xy':
                        datafid.write(2 * ss + 'RhoZ' + pol + 2 * ss + str(ii + 1) + 2 * ss + '0' +
                                      2 * ss + '1' + 2 * ss + string_fmt % rp.resxy[ii] + 2 * ss +
                                      string_fmt % rp.resxyerr[ii] + '\n')
                    elif pol == 'yx':
                        datafid.write(2 * ss + 'RhoZ' + pol + 2 * ss + str(ii + 1) + 2 * ss + '0' +
                                      2 * ss + '1' + 2 * ss + string_fmt % rp.resyx[ii] + 2 * ss +
                                      string_fmt % rp.resyxerr[ii] + '\n')
                    elif pol == 'det':
                        datafid.write(2 * ss + 'RhoZ' + pol + 2 * ss + str(ii + 1) + 2 * ss + '0' +
                                      2 * ss + '1' + 2 * ss + string_fmt % rp.resyx[ii] + 2 * ss +
                                      string_fmt % rp.resyxerr[ii] + '\n')
                else:
                    if pol == 'xy':
                        datafid.write(2 * ss + 'RhoZ' + pol + 2 * ss + str(ii + 1) + 2 * ss + '0' +
                                      2 * ss + '1' + 2 * ss + string_fmt % rp.resxy[ii] + 2 * ss +
                                      string_fmt % (rp.resxy[ii] * reserr / 100.) + '\n')
                    elif pol == 'yx':
                        datafid.write(2 * ss + 'RhoZ' + pol + 2 * ss + str(ii + 1) + 2 * ss + '0' +
                                      2 * ss + '1' + 2 * ss + string_fmt % rp.resyx[ii] + 2 * ss +
                                      string_fmt % (rp.resyx[ii] * reserr / 100.) + '\n')

                # write phase components
                if phaseerr == 'data':
                    if pol == 'xy':
                        datafid.write(2 * ss + 'PhsZ' + pol + 2 * ss + str(ii + 1) + 2 * ss + '0' +
                                      2 * ss + '1' + 2 * ss + string_fmt % rp.phasexy[ii] + 2 * ss +
                                      string_fmt % rp.phasexyerr[ii] + '\n')
                    if pol == 'yx':
                        datafid.write(2 * ss + 'PhsZ' + pol + 2 * ss + str(ii + 1) + 2 * ss + '0' +
                                      2 * ss + '1' + 2 * ss + string_fmt % rp.phaseyx[ii] + 2 * ss +
                                      string_fmt % rp.phaseyxerr[ii] + '\n')
                else:
                    if pol == 'xy':
                        datafid.write(2 * ss + 'PhsZ' + pol + 2 * ss + str(ii + 1) + 2 * ss + '0' +
                                      2 * ss + '1' + 2 * ss + string_fmt % rp.phasexy[ii] + 2 * ss +
                                      string_fmt % (phaseerr / 100. * (180 / np.pi)) + '\n')
                    if pol == 'yx':
                        datafid.write(2 * ss + 'PhsZ' + pol + 2 * ss + str(ii + 1) + 2 * ss + '0' +
                                      2 * ss + '1' + 2 * ss + string_fmt % rp.phaseyx[ii] + 2 * ss +
                                      string_fmt % (phaseerr / 100. * (180 / np.pi)) + '\n')
            datafid.close()
            print('Wrote Data File: ', dfilesave)

    def make1DModelFile(self, savepath=None, nlayers=100, bottomlayer=10000,
                        basestep=10, z1layer=10, airlayerheight=10000):
        """
        Makes a 1D model file for Occam1D.  

        Arguments:
        ----------

            **savepath** :path to save file to, if just path saved as 
                          savepath\model.mod, if None defaults to dirpath

            **nlayers** : number of layers

            **bottomlayer** : depth of bottom layer in meters

            **basestep** : numerical base of logarithmic depth step 10 or 2 or 
                          1 for linear

            **z1layer** : depth of first layer in meters

            **airlayerheight** : height of air layers in meters

        Returns:
        --------

            **Occam1D.modelfn** = full path to model file

        ..Note: This needs to be redone.

        :Example: ::

            >>> old = occam.Occam1D()
            >>> old.make1DModelFile(savepath=r"/home/Occam1D/Line1/Inv1_TE",
            >>>                     nlayers=50,bottomlayer=10000,z1layer=50)
            >>> Wrote Model file: /home/Occam1D/Line1/Inv1_TE/Model1D 
        """

        ss = '   '
        # if the savepath was not entered test to see if there is one
        if savepath == None:
            if not self.savepath:
                raise IOError('No savepath found.  Please input one.')
            self.modelfn = os.path.join(self.savepath, 'Model1D')

        # if the save path was entered as just a path
        elif os.path.basename(savepath).find('.') == -1:
            if not os.path.exists(savepath):
                os.mkdir(savepath)
            self.savepath = savepath
            self.modelfn = os.path.join(self.savepath, 'Model1D')

        # if the save path was entered as a file with full path
        else:
            self.modelfn = savepath

        #---------need to refine this--------------------

        layers = np.logspace(
            np.log10(z1layer), np.log10(bottomlayer), num=nlayers)

        # make the model file
        modfid = open(self.modelfn, 'w')
        modfid.write('Format: Resistivity1DMod_1.0' + '\n')
        modfid.write('#LAYERS:    ' + str(nlayers + 2) + '\n')
        modfid.write('!Set free values to -1 or ? \n')
        modfid.write('!penalize between 1 and 0,' +
                     '0 allowing jump between layers and 1 smooth. \n')
        modfid.write(
            '!preference is the assumed resistivity on linear scale. \n')
        modfid.write(
            '!pref_penalty needs to be put if preference is not 0 [0,1]. \n')
        modfid.write('! top_depth' + ss + 'resistivity' + ss + 'penalty' + ss + 'preference' + ss +
                     'pref_penalty \n')
        modfid.write(ss + '-10000' + ss + '1d12' + ss + '0' +
                     ss + '0' + ss + '0' + ss + '!air layer \n')
        modfid.write(ss + '0' + ss + '-1' + ss + '0' + ss +
                     '0' + ss + '0' + ss + '!first ground layer \n')
        for ll in layers:
            modfid.write(ss + '{0:.2f}'.format(ll) + ss + '-1' + ss + '1' + ss + '0' + ss + '0' +
                         '\n')

        modfid.close()

        print('Wrote Model file: ', self.modelfn)

    def make1DInputFile(self, savepath=None, imode='TE', roughtype=1,
                        maxiter=20, targetrms=1.0, rhostart=100,
                        description='1dInv', lagrange=5.0, roughness=1.0E7,
                        debuglevel=1, iteration=0, misfit=100.0):
        """
        Make a 1D input file for Occam 1D

        Arguments:
        ---------
            **savepath** : full path to save input file to, if just path then 
                           saved as savepath/input

            **modelfile** : full path to model file, if None then assumed to be in 
                            savepath/model.mod

            **datafile** : full path to data file, if None then assumed to be 
                            in savepath/TE.dat or TM.dat

            **roughtype** : roughness type. *default* = 0

            **maxiter** : maximum number of iterations. *default* = 20 

            **targetrms** : target rms value. *default* = 1.0

            **rhostart** : starting resistivity value on linear scale. 
                            *default* = 100

            **description** : description of the inversion. 

            **lagrange** : starting Lagrange multiplier for smoothness.
                           *default* = 5

            **roughness** : starting roughness value. *default* = 1E7

            **debuglevel** : something to do with how Fortran debuggs the code
                             Almost always leave at *default* = 1

            **iteration** : the starting iteration number, handy if the
                            starting model is from a previous run.
                            *default* = 0

            **misfit** : starting misfit value. *default* = 100

        Returns:
        --------
            **Occam1D.inputfn** : full path to input file. 

        :Example: ::

            >>> old = occam.Occam1D()
            >>> old.make1DdataFile('MT01',edipath=r"/home/Line1",
            >>>                    savepath=r"/home/Occam1D/Line1/Inv1_TE",
            >>>                    mode='TE')
            >>> Wrote Data File: /home/Occam1D/Line1/Inv1_TE/MT01TE.dat
            >>>            
            >>> old.make1DModelFile(savepath=r"/home/Occam1D/Line1/Inv1_TE",
            >>>                     nlayers=50,bottomlayer=10000,z1layer=50)
            >>> Wrote Model file: /home/Occam1D/Line1/Inv1_TE/Model1D
            >>>
            >>> old.make1DInputFile(rhostart=10,targetrms=1.5,maxiter=15)
            >>> Wrote Input File: /home/Occam1D/Line1/Inv1_TE/Input1D  
        """

        ss = '   '

        # make input data file name
        # if no savepath is input, test if there is already one
        if savepath == None:
            if not self.savepath:
                raise IOError('No savepath.  Please input one.')
            self.inputfn = os.path.join(self.savepath, 'Input1D')

        # if the savepath was input as just a path
        elif os.path.basename(savepath).find('.') == -1:
            if not os.path.exists(savepath):
                os.mkdir(savepath)
            self.inputfn = os.path.join(savepath, 'Input1D')

        # if the savepath was input as a full path to file
        else:
            self.inputfn = savepath

        if not self.modelfn:
            if self.savepath:
                self.modelfn = os.path.join(self.savepath, 'Model1D')
            else:
                raise IOError('No savepath.  Please input one.')

        # try to get data file name
        if imode == 'TE':
            if not self.datafn_te:
                if self.savepath:
                    try:
                        self.datafn_te = [os.path.join(self.savepath, dd)
                                          for dd in os.listdir(self.savepath)
                                          if dd.find('TE.dat') > 0][0]
                    except IndexError:
                        raise IOError(
                            'No TE data file found. Please input one.')
                else:
                    raise IOError('No savepth found. Please input one.')
            else:
                pass
        if imode == 'TM':
            if not self.datafn_tm:
                if self.savepath:
                    try:
                        self.datafn_tm = [os.path.join(self.savepath, dd)
                                          for dd in os.listdir(self.savepath)
                                          if dd.find('TM.dat') > 0][0]
                    except IndexError:
                        raise IOError(
                            'No TM data file found. Please input one.')
            else:
                pass

        # read in the model and get number of parameters
        self.read1DModelFile()
        paramcount = self.mdict['nparam']

        # write input file
        infid = open(self.inputfn, 'w')
        infid.write(
            'Format:             OCCAMITER_FLEX      ! Flexible format \n')
        infid.write('Description:        ' + description +
                    '     !For your own notes. \n')
        infid.write('Model File:         ' + self.modelfn + '       \n')
        if imode == 'TE':
            infid.write('Data File:          ' + self.datafn_te + '        \n')
        if imode == 'TM':
            infid.write('Data File:          ' + self.datafn_tm + '        \n')
        infid.write('Date/Time:          ' + time.ctime() + '\n')
        infid.write('Max Iter:           ' + str(maxiter) + '\n')
        infid.write('Target Misfit:      ' + str(targetrms) + '\n')
        infid.write('Roughness Type:     ' + str(roughtype) + '\n')
        infid.write('!Model Bounds:      min,max             ! Optional, places bounds' +
                    ' on log10(rho) values. \n')
        infid.write('!Model Value Steps: stepsize            ! Optional, forces model' +
                    ' into discrete steps of stepsize. \n')
        infid.write('Debug Level:        ' + str(debuglevel) +
                    ' ' * 19 + '! Console output. ' +
                    '0: minimal, 1: default, 2: detailed \n')
        infid.write('Iteration:          ' + str(iteration) +
                    ' ' * 19 + '! Iteration number,' +
                    ' use 0 for starting from scratch. \n')
        infid.write('Lagrange Value:     ' + str(lagrange) +
                    ' ' * 17 + '! log10(largrance ' +
                    'multiplier), starting value.\n')
        infid.write('Roughness Value:    ' + str(roughness) +
                    ' ' * 10 + '! Roughness of last' +
                    ' model, ignored on startup. \n')
        infid.write('Misfit Value:       ' + str(misfit) +
                    ' ' * 15 + '! Misfit of model listed' +
                    'below. Ignored on startup.\n')
        infid.write('Misfit Reached:     0	                ! 0: not reached,' +
                    ' 1: reached.  Useful when restarting.\n')
        infid.write('Param Count:        ' + str(paramcount) +
                    ' ' * 17 + '! Number of free' +
                    ' inversion parameters. \n')
        for ii in range(paramcount):
            infid.write(ss + str(np.log10(rhostart)) + '\n')

        infid.close()
        print('Wrote Input File: ', self.inputfn)

    def read1DModelFile(self):
        """

        will read in model 1D file

        Arguments:
        ----------
            **modelfn** : full path to model file

        Returns:
        --------
            **Occam1D.mdict** : dictionary of values with keys: 

                *'depth'* : depth of model in meters

                *'res'* : value of resisitivity

                *'pen'* : penalty

                *'pre'* : preference

                *'prefpen'* : preference penalty

        :Example: ::

            >>> old = occam.Occam1d()
            >>> old.savepath = r"/home/Occam1D/Line1/Inv1_TE"
            >>> old.read1DModelFile()
        """
        if not self.modelfn:
            if not self.savepath:
                raise IOError('No model file found.  Please input one.')
            self.modelfn = os.path.join(self.savepath, 'Model1D')

        mfid = open(self.modelfn, 'r')
        mlines = mfid.readlines()
        mfid.close()
        try:
            self.mdict
        except AttributeError:
            mdict = {}
            mdict['nparam'] = 0
            for key in ['depth', 'res', 'pen', 'pref', 'prefpen']:
                mdict[key] = []

            for mm, mline in enumerate(mlines):
                if mline.find('!') == 0:
                    pass
                elif mline.find(':') >= 0:
                    mlst = mline.strip().split(':')
                    mdict[mlst[0]] = mlst[1]
                else:
                    mlst = mlst = mline.strip().split()
                    mdict['depth'].append(float(mlst[0]))
                    if mlst[1] == '?':
                        mdict['res'].append(-1)
                    elif mlst[1] == '1d12':
                        mdict['res'].append(1.0E12)
                    else:
                        try:
                            mdict['res'].append(float(mlst[1]))
                        except ValueError:
                            mdict['res'].append(-1)
                    mdict['pen'].append(float(mlst[2]))
                    mdict['pref'].append(float(mlst[3]))
                    mdict['prefpen'].append(float(mlst[4]))
                    if mlst[1] == '-1' or mlst[1] == '?':
                        mdict['nparam'] += 1

            # make everything an array
            for key in ['depth', 'res', 'pen', 'pref', 'prefpen']:
                mdict[key] = np.array(mdict[key])

            # create an array with empty columns to put the TE and TM models
            # into
            mres = np.zeros((len(mdict['res']), 3))
            mres[:, 0] = mdict['res']
            mdict['res'] = mres
            # make dictionary an attribute of Occam1D class
            self.mdict = mdict

    def read1DInputFile(self):
        """
        reads in a 1D input file

        Arguments:
        ---------
            **inputfn** : full path to input file

        Returns:
        --------
            **Occam1D.indict** : dictionary with keys following the header and

                *'res'* : an array of resistivity values

        :Example: ::

            >>> old = occam.Occam1d()
            >>> old.savepath = r"/home/Occam1D/Line1/Inv1_TE"
            >>> old.read1DInputFile()
        """
        if not self.inputfn:
            if not self.savepath:
                raise IOError('No input file found.  Please input one.')
            self.inputfn = os.path.join(self.savepath, 'Input1D')

        infid = open(self.inputfn, 'r')
        ilines = infid.readlines()
        infid.close()

        self.indict = {}
        res = []

        # split the keys and values from the header information
        for iline in ilines:
            if iline.find(':') >= 0:
                ikey = iline[0:20].strip()
                ivalue = iline[20:].split('!')[0].strip()
                self.indict[ikey[:-1]] = ivalue
            else:
                try:
                    res.append(float(iline.strip()))
                except ValueError:
                    pass

        # make the resistivity array ready for models to be input
        self.indict['res'] = np.zeros((len(res), 3))
        self.indict['res'][:, 0] = res

        # get data file
        if self.indict['Data File'].find('TE') > 0:
            self.datafn_te = self.indict['Data File']

        elif self.indict['Data File'].find('TM') > 0:
            self.datafn_tm = self.indict['Data File']

    def read1DdataFile(self, imode='TE'):
        """
        reads a 1D data file

        Arguments:
        ----------
            **datafile** : full path to data file

            **imode** : mode to read from can be TE or TM

        Returns:
        --------
            **Occam1D.rpdict** : dictionary with keys:

                *'freq'* : an array of frequencies with length nf

                *'resxy'* : TE resistivity array with shape (nf,4) for (0) data,
                          (1) dataerr, (2) model, (3) modelerr

                *'resyx'* : TM resistivity array with shape (nf,4) for (0) data,
                          (1) dataerr, (2) model, (3) modelerr

                *'phasexy'* : TE phase array with shape (nf,4) for (0) data,
                            (1) dataerr, (2) model, (3) modelerr

                *'phaseyx'* : TM phase array with shape (nf,4) for (0) data,
                            (1) dataerr, (2) model, (3) modelerr

        :Example: ::

            >>> old = occam.Occam1d()
            >>> old.datafn_te = r"/home/Occam1D/Line1/Inv1_TE/MT01TE.dat"
            >>> old.read1DdataFile()
        """

        # get the data file for the correct mode
        if imode == 'TE':
            if not self.datafn_te:
                raise IOError('No TE data file found.  Please input one.')

            dfid = open(self.datafn_te, 'r')

        elif imode == 'TM':
            if not self.datafn_tm:
                raise IOError('No TM data file found.  Please input one.')

            dfid = open(self.datafn_te, 'r')

        #read in lines
        dlines = dfid.readlines()
        dfid.close()

        # make a dictionary of all the fields found so can put them into arrays
        finddict = {}
        for ii, dline in enumerate(dlines):
            if dline.find('#') <= 3:
                fkey = dline[2:].strip().split(':')[0]
                fvalue = ii
                finddict[fkey] = fvalue

        # get number of frequencies
        nfreq = int(dlines[finddict['Frequencies']][
                    2:].strip().split(':')[1].strip())

        # frequency list
        freq = np.array([float(ff) for ff in dlines[finddict['Frequencies'] + 1:
                                                    finddict['Receivers']]])

        # data dictionary to put things into
        # check to see if there is alread one, if not make a new one
        try:
            self.rpdict
        except NameError:
            self.rpdict = {'freq': freq,
                           'resxy': np.zeros((4, nfreq)),
                           'resyx': np.zeros((4, nfreq)),
                           'phasexy': np.zeros((4, nfreq)),
                           'phaseyx': np.zeros((4, nfreq))
                           }

        # get data
        for dline in dlines[finddict['Data'] + 1:]:
            if dline.find('!') == 0:
                pass
            else:
                dlst = dline.strip().split()
                if len(dlst) > 4:
                    jj = int(dlst[1]) - 1
                    dvalue = float(dlst[4])
                    derr = float(dlst[5])
                    if dlst[0] == 'RhoZxy' or dlst[0] == '103':
                        self.rpdict['resxy'][0, jj] = dvalue
                        self.rpdict['resxy'][1, jj] = derr
                    if dlst[0] == 'PhsZxy' or dlst[0] == '104':
                        self.rpdict['phasexy'][0, jj] = dvalue
                        self.rpdict['phasexy'][1, jj] = derr
                    if dlst[0] == 'RhoZyx' or dlst[0] == '105':
                        self.rpdict['resyx'][0, jj] = dvalue
                        self.rpdict['resyx'][1, jj] = derr
                    if dlst[0] == 'PhsZyx' or dlst[0] == '106':
                        self.rpdict['phaseyx'][0, jj] = dvalue
                        self.rpdict['phaseyx'][1, jj] = derr

    def read1DIterFile(self, iterfn, imode='TE'):
        """
        read an 1D iteration file

        Arguments:
        ----------
            **imode** : mode to read from 

        Returns:
        --------
            **Occam1D.itdict** : dictionary with keys of the header:

            **Occam1D.mdict['res']** : fills this array with the appropriate 
                                        values (0) for data, (1) TE, (2) TM

        :Example: ::

            >>> old = occam.Occam1d()
            >>> old.read1DIterFile(r"/home/Occam1D/Inv1_TE/M01TE_15.iter")

        """

        if not self.savepath:
            self.savepath = os.path.dirname(iterfn)

        self.read1DModelFile()

        freeparams = np.where(self.mdict['res'] == -1)[0]

        if imode == 'TE':
            self.iterfn_te = iterfn
            ifid = open(self.iterfn_te, 'r')
        elif imode == 'TM':
            self.iterfn_tm = iterfn
            ifid = open(self.iterfn_tm, 'r')

        ilines = ifid.readlines()
        ifid.close()

        self.itdict = {}
        model = []
        for ii, iline in enumerate(ilines):
            if iline.find(':') >= 0:
                ikey = iline[0:20].strip()
                ivalue = iline[20:].split('!')[0].strip()
                self.itdict[ikey[:-1]] = ivalue
            else:
                try:
                    ilst = iline.strip().split()
                    for kk in ilst:
                        model.append(float(kk))
                except ValueError:
                    pass

        # put the model values into the model dictionary into the res array
        # for easy manipulation and access.  Also so you can compare TE and TM
        model = np.array(model)
        if imode == 'TE':
            self.mdict['res'][:, 1] = self.mdict['res'][:, 0]
            self.mdict['res'][freeparams, 1] = model
        if imode == 'TM':
            self.mdict['res'][:, 2] = self.mdict['res'][:, 0]
            self.mdict['res'][freeparams, 2] = model

    def read1DRespFile(self, respfn, imode='TE'):
        """
        read response file

        Arguments:
        ---------
            **repsfn** : full path to response file

        Returns:
        --------
            *Occam1D.*rpdict** : dictionary with keys:

                *freq* : an array of frequencies with length nf

                *resxy* : TE resistivity array with shape (nf,4) for (0) data,
                          (1) dataerr, (2) model, (3) modelerr

                *resyx* : TM resistivity array with shape (nf,4) for (0) data,
                          (1) dataerr, (2) model, (3) modelerr

                *phasexy* : TE phase array with shape (nf,4) for (0) data,
                            (1) dataerr, (2) model, (3) modelerr

                *phaseyx* : TM phase array with shape (nf,4) for (0) data,
                            (1) dataerr, (2) model, (3) modelerr

       :Example: ::

            >>> old = occam.Occam1d()
            >>> old.read1DRespFile(r"/home/Occam1D/Inv1_TE/M01TE_15.resp")
        """

        if imode == 'TE':
            self.respfn_te = respfn
        elif imode == 'TM':
            self.respfn_tm = respfn

        if not self.savepath:
            self.savepath = os.path.dirname(respfn)

        dfid = open(respfn, 'r')

        dlines = dfid.readlines()
        dfid.close()

        finddict = {}
        for ii, dline in enumerate(dlines):
            if dline.find('#') <= 3:
                fkey = dline[2:].strip().split(':')[0]
                fvalue = ii
                finddict[fkey] = fvalue
        nfreq = int(dlines[finddict['Frequencies']][
                    2:].strip().split(':')[1].strip())

        # frequency list
        freq = np.array([float(ff) for ff in dlines[finddict['Frequencies'] + 1:
                                                    finddict['Receivers']]])

        # data dictionary
        try:
            self.rpdict
        except AttributeError:
            self.rpdict = {'freq': freq,
                           'resxy': np.zeros((4, nfreq)),
                           'resyx': np.zeros((4, nfreq)),
                           'phasexy': np.zeros((4, nfreq)),
                           'phaseyx': np.zeros((4, nfreq))
                           }

        for dline in dlines[finddict['Data'] + 1:]:
            if dline.find('!') == 0:
                pass
            else:
                dlst = dline.strip().split()
                if len(dlst) > 4:
                    jj = int(dlst[1]) - 1
                    dvalue = float(dlst[4])
                    derr = float(dlst[5])
                    rvalue = float(dlst[6])
                    rerr = float(dlst[7])
                    if dlst[0] == 'RhoZxy' or dlst[0] == '103':
                        self.rpdict['resxy'][0, jj] = dvalue
                        self.rpdict['resxy'][1, jj] = derr
                        self.rpdict['resxy'][2, jj] = rvalue
                        self.rpdict['resxy'][3, jj] = rerr
                    if dlst[0] == 'PhsZxy' or dlst[0] == '104':
                        self.rpdict['phasexy'][0, jj] = dvalue
                        self.rpdict['phasexy'][1, jj] = derr
                        self.rpdict['phasexy'][2, jj] = rvalue
                        self.rpdict['phasexy'][3, jj] = rerr
                    if dlst[0] == 'RhoZyx' or dlst[0] == '105':
                        self.rpdict['resyx'][0, jj] = dvalue
                        self.rpdict['resyx'][1, jj] = derr
                        self.rpdict['resyx'][2, jj] = rvalue
                        self.rpdict['resyx'][3, jj] = rerr
                    if dlst[0] == 'PhsZyx' or dlst[0] == '106':
                        self.rpdict['phaseyx'][0, jj] = dvalue
                        self.rpdict['phaseyx'][1, jj] = derr
                        self.rpdict['phaseyx'][2, jj] = rvalue
                        self.rpdict['phaseyx'][3, jj] = rerr

    def plot1D(self, iternum=10, savepath=None, iterfn=None, respfn=None,
               imode='TE', fignum=1, ms=4, dpi=150, fs=10, lw=2, dlimits=None):
        """
        Plots the results of a 1D inversion.  The left plot is the response
        and the right hand plot is the model as a function of depth.

        Arguments:
        ----------
            **iternum** : iteration number to plot. *Default* is 10 

            **savepath** : path to the 1D inversion response and iteration 
                           files. *Default* is None.

            **iterfn** : full path to iteration file if savepath is not entered
                         *Default* is None.

            **respfn** : full path to response file if savepath is not entered.
                         *Default* is None.

            **imode** : mode to plot. can be input as:
                * 'TE' for TE mode
                * 'TM' for TM mode
                * 'both' for both TE and TM modes
                * Default* is 'TE'

            **fignum** : figure number that the plot will be. *Default* is 1

            **ms** : marker size. *Default* is 4

            **dpi** : dots per inch resolution of the plot. *Default* is 150.

            **fs** : font size of the labels. *Default* is 10

            **dlimits** : limits on the depth axes. Input as a tuple 
                          (dmin,dmax).  *Default* is None.

        :Example: ::

            >>> old = occam.Occam1d()
            >>> old.savepath = r"/home/Occam1D/Line1/Inv1_TE"
            >>> #Look at only the interval between 5 and 10 kilometers
            >>> old.plot1D(dlimits=(5,10))
        """

        self.iternum = iternum
        # get files
        try:
            self.modelfn
        except AttributeError:
            if not self.dirpath:
                self.dirpath = os.path.dirname(respfn)

            self.modelfn = os.path.join(self.dirpath, 'Model1D')
            if os.path.isfile(self.modelfn) == False:
                raise IOError('Could not find ' + self.modelfn)

        #-------------read in response files---------------------
        if respfn == None:
            if imode == 'TE':
                try:
                    self.respfn_te = [os.path.join(self.savepath, rr)
                                      for rr in os.listdir(self.savepath)
                                      if rr.find('TE_{0}.resp'.format(self.iternum)) > 0][0]

                    self.read1DRespFile(self.respfn_te, imode='TE')

                except IndexError:
                    raise IOError('Could not find response TE file.')
            elif imode == 'TM':
                try:
                    self.respfn_tm = [os.path.join(self.savepath, rr)
                                      for rr in os.listdir(self.savepath)
                                      if rr.find('TM_{0}.resp'.format(self.iternum)) > 0][0]

                    self.read1DRespFile(self.respfn_tm, imode='TM')
                except IndexError:
                    raise IOError('Could not find response TM file.')
            elif imode == 'both':
                try:
                    self.respfn_te = [os.path.join(self.savepath, rr)
                                      for rr in os.listdir(self.savepath)
                                      if rr.find('TE_{0}.resp'.format(self.iternum)) > 0][0]
                    self.read1DRespFile(self.respfn_te, imode='TE')
                except IndexError:
                    raise IOError('Could not find response TE file.')
                try:
                    self.respfn_tm = [os.path.join(self.savepath, rr)
                                      for rr in os.listdir(self.savepath)
                                      if rr.find('TM_{0}.resp'.format(self.iternum)) > 0][0]
                    self.read1DRespFile(self.respfn_tm, imode='TM')
                except IndexError:
                    raise IOError('Could not find response TM file.')

        # if the response files are input read them in
        else:
            if imode == 'TE':
                self.respfn_te = respfn
                self.read1DRespFile(self.respfn_te, imode='TE')

            elif imode == 'TM':
                self.respfn_tm = respfn
                self.read1DRespFile(self.respfn_tm, imode='TM')

            elif imode == 'both':
                if type(iterfn) is not list or type(iterfn) is not tuple:
                    raise IOError(
                        'Please enter iteration files as a list or tuple.')
                self.respfn_te = respfn[0]
                self.read1DRespFile(self.respfn_te, imode='TE')

                self.respfn_tm = respfn[1]
                self.read1DRespFile(self.respfn_tm, imode='TM')

        #------Read in iteration files--------------------
        if iterfn == None:
            if imode == 'TE':
                try:
                    self.iterfn_te = [os.path.join(self.savepath, rr)
                                      for rr in os.listdir(self.savepath)
                                      if rr.find('TE_{0}.iter'.format(self.iternum)) > 0][0]
                    print(self.iterfn_te)
                    self.read1DIterFile(self.iterfn_te, imode='TE')

                except IndexError:
                    raise IOError('Could not find iteration TE file.')
            elif imode == 'TM':
                try:
                    self.iterfn_tm = [os.path.join(self.savepath, rr)
                                      for rr in os.listdir(self.savepath)
                                      if rr.find('TM_{0}.iter'.format(self.iternum)) > 0][0]

                    self.read1DIterFile(self.iterfn_tm, imode='TM')
                except IndexError:
                    raise IOError('Could not find iteration TM file.')
            elif imode == 'both':
                try:
                    self.iterfn_te = [os.path.join(self.savepath, rr)
                                      for rr in os.listdir(self.savepath)
                                      if rr.find('TE_{0}.iter'.format(self.iternum)) > 0][0]
                    self.read1DIterFile(self.iterfn_te, imode='TE')
                except IndexError:
                    raise IOError('Could not find iteration TE file.')
                try:
                    self.iterfn_tm = [os.path.join(self.savepath, rr)
                                      for rr in os.listdir(self.savepath)
                                      if rr.find('TM_{0}.iter'.format(self.iternum)) > 0][0]
                    self.read1DIterFile(self.iterfn_tm, imode='TM')
                except IndexError:
                    raise IOError('Could not find iteration TM file.')
        else:
            if imode == 'TE':
                self.iterfn_te = iterfn
                self.read1DIterFile(self.iterfn_te, imode='TE')

            elif imode == 'TM':
                self.iterfn_tm = iterfn
                self.read1DIterFile(self.iterfn_tm, imode='TM')

            elif imode == 'both':
                if type(iterfn) is not list or type(iterfn) is not tuple:
                    raise IOError(
                        'Please enter iteration files as a list or tuple.')
                self.iterfn_te = iterfn[0]
                self.read1DIterFile(self.iterfn_te, imode='TE')

                self.iterfn_tm = iterfn[1]
                self.read1DIterFile(self.iterfn_tm, imode='TM')

        period = 1 / self.rpdict['freq']

        # make a grid of subplots
        gs = gridspec.GridSpec(6, 5, hspace=.25, wspace=.75)

        # make a figure
        fig = plt.figure(fignum, [8, 8], dpi=dpi)
        plt.clf()

        # set some plot parameters
        plt.rcParams['font.size'] = fs - 2
        plt.rcParams['figure.subplot.left'] = .1
        plt.rcParams['figure.subplot.right'] = .93
        plt.rcParams['figure.subplot.bottom'] = .1
        plt.rcParams['figure.subplot.top'] = .90

        # subplot resistivity
        axr = fig.add_subplot(gs[:4, :4])

        # subplot for phase
        axp = fig.add_subplot(gs[4:, :4], sharex=axr)

        # check for data in resistivity
        rxy = np.where(self.rpdict['resxy'][0] != 0)[0]
        ryx = np.where(self.rpdict['resyx'][0] != 0)[0]

        pxy = np.where(self.rpdict['phasexy'][0] != 0)[0]
        pyx = np.where(self.rpdict['phaseyx'][0] != 0)[0]

        # check to make sure a model was read in for resistivity
        rxym = np.where(self.rpdict['resxy'][2] != 0)[0]
        ryxm = np.where(self.rpdict['resyx'][2] != 0)[0]

        pxym = np.where(self.rpdict['phasexy'][2] != 0)[0]
        pyxm = np.where(self.rpdict['phaseyx'][2] != 0)[0]

        #----------Plot TE mode-------------------
        if imode == 'TE':
            titlestr = '$Z_{TE}$'
            # plot data resistivity
            if len(rxy) != 0:
                r1 = axr.loglog(period[rxy], self.rpdict['resxy'][0][rxy],
                                ls='None', marker='o', color='k', mfc='k', ms=ms)

            # plot data phase
            if len(pxy) != 0:
                p1 = axp.semilogx(period[pxy], self.rpdict['phasexy'][0][pxy],
                                  ls='None', marker='o', color='k', mfc='k', ms=ms)

            # plot model resistivity
            if len(rxym) != 0:
                r2 = axr.loglog(period[rxym], self.rpdict['resxy'][2][rxym],
                                ls=':', color='b', lw=lw)
            # plot model phase
            if len(pxym) != 0:
                p2 = axp.semilogx(period[pxym], self.rpdict['phasexy'][2][pxym],
                                  ls=':', color='b', lw=lw)

            # add legend
            axr.legend([r1[0], r2[0]], ['Data', 'Model'], loc='upper left',
                       markerscale=1,
                       borderaxespad=.15,
                       labelspacing=.18,
                       handletextpad=.15, borderpad=.15)

        #--------Plot TM mode-----------------------
        elif imode == 'TM':
            titlestr = '$Z_{TM}$'
            # plot data resistivity
            if len(ryx) != 0:
                r1 = axr.loglog(period[ryx], self.rpdict['resyx'][0][ryx],
                                ls='None', marker='o', color='k', mfc='k', ms=ms)
            # plot data phase
            if len(pyx) != 0:
                p1 = axp.semilogx(period[pyx], self.rpdict['phaseyx'][0][pyx],
                                  ls='None', marker='o', color='k', mfc='k', ms=ms)

            # plot model resistivity
            if len(ryxm) != 0:
                r2 = axr.loglog(period[ryxm], self.rpdict['resyx'][2][ryxm],
                                ls=':', color='b', lw=lw)
            # plot model phase
            if len(pyxm) != 0:
                p2 = axp.semilogx(period[pyxm], self.rpdict['phaseyx'][2][pyxm],
                                  ls=':', color='b', lw=lw)

            axr.legend([r1[0], r2[0]], ['Data', 'Model'],
                       loc='upper left', markerscale=1,
                       borderaxespad=.15,
                       labelspacing=.18,
                       handletextpad=.15, borderpad=.15)

        #-------------Plot Both Modes--------------------------------
        elif imode == 'both':
            titlestr = '$Z_{TE}$ and $Z_{TM}$'
            # plot data resistivity
            if len(rxy) != 0:
                r1te = axr.loglog(period[rxy], self.rpdict['resxy'][0][rxy],
                                  ls='None', marker='s', color='k', mfc='k', ms=ms)
            if len(ryx) != 0:
                r1tm = axr.loglog(period[ryx], self.rpdict['resyx'][0][ryx],
                                  ls='None', marker='o', color='k', mfc='k', ms=ms)

            # plot data phase
            if len(pxy) != 0:
                p1te = axp.semilogx(period[pxy], self.rpdict['phasexy'][0][pxy],
                                    ls='None', marker='s', color='k', mfc='k', ms=ms)

            if len(pyx) != 0:
                p1tm = axp.semilogx(period[pyx], self.rpdict['phaseyx'][0][pyx],
                                    ls='None', marker='o', color='k', mfc='k', ms=ms)

            # plot model resistivity
            if len(rxym) != 0:
                r2te = axr.loglog(period[rxym], self.rpdict['resxy'][2][rxym],
                                  ls=':', color='b', lw=lw)

            if len(ryxm) != 0:
                r2tm = axr.loglog(period[ryxm], self.rpdict['resyx'][2][ryxm],
                                  ls=':', color='r', lw=lw)
            # plot model phase
            if len(pxym) != 0:
                p2 = axp.semilogx(period[pxym], self.rpdict['phasexy'][2][pxym],
                                  ls=':', color='b', lw=lw)
            if len(pyxm) != 0:
                p2 = axp.semilogx(period[pyxm], self.rpdict['phaseyx'][2][pyxm],
                                  ls=':', color='r', lw=lw)

            # add legend
            axr.legend([r1te[0], r2te[0], r1tm[0], r2tm[0]],
                       ['Data$_{TE}$', 'Model$_{TE}$',
                        'Data$_{TM}$', 'Model$_{TM}$'],
                       loc='upper left', markerscale=1,
                       borderaxespad=.15,
                       labelspacing=.18,
                       handletextpad=.15, borderpad=.15)

        axr.grid(True, alpha=.4, which='both')
        plt.setp(axr.xaxis.get_ticklabels(), visible=False)
        axp.grid(True, alpha=.4, which='both')
        axp.yaxis.set_major_locator(MultipleLocator(10))
        axp.yaxis.set_minor_locator(MultipleLocator(1))

        axr.set_ylabel('App. Res. ($\Omega \cdot m$)',
                       fontdict={'size': fs, 'weight': 'bold'})
        axp.set_ylabel('Phase (deg)',
                       fontdict={'size': fs, 'weight': 'bold'})
        axp.set_xlabel('Period (s)', fontdict={'size': fs, 'weight': 'bold'})
        plt.suptitle(titlestr, fontsize=fs + 2, fontweight='bold')

        #------------------plot 1D inversion---------------------------------
        axm = fig.add_subplot(gs[:, 4])
        depthp = self.mdict['depth'][1:]
        if imode == 'TE':
            modelresp = abs(10**self.mdict['res'][1:, 1])
            axm.loglog(modelresp[::-1], depthp[::-1], ls='steps-', color='b',
                       lw=lw)
        elif imode == 'TM':
            modelresp = abs(10**self.mdict['res'][1:, 2])
            axm.loglog(modelresp[::-1], depthp[::-1], ls='steps-', color='b',
                       lw=lw)
        elif imode == 'both':
            modelrespte = abs(10**self.mdict['res'][1:, 1])
            axm.loglog(modelrespte[::-1], depthp[::-1], ls='steps-', color='b',
                       lw=lw)
            modelresptm = abs(10**self.mdict['res'][1:, 2])
            axm.loglog(modelresptm[::-1], depthp[::-1], ls='steps-', color='r',
                       lw=lw)

        if dlimits == None:
            axm.set_ylim(ymin=depthp[-1], ymax=depthp[0])
        else:
            axm.set_ylim(dlimits)
        axm.set_ylabel('Depth (m)', fontdict={'size': fs, 'weight': 'bold'})
        axm.set_xlabel('Resistivity ($\Omega \cdot m$)',
                       fontdict={'size': fs, 'weight': 'bold'})
        axm.grid(True, which='both', alpha=.4)

        plt.show()

    def plotL2Curve(self, savepath=None, imode='TE', fignum=1, dpi=150, fs=10):
        """
        Plot the L curve for RMS vs Iteration and RMS vs Roughness.

        Arguments:
        ----------

            **savepath** : path to iteration files

            **imode** : mode to plot.  Can be:
                    * 'TE' for TE mode. *Default*
                    * 'TM' for TM mode
                    * 'both' for both TE ant TM modes

            **fignum** : figure number for the plot

            **dpi** : dots per inch resolution of the plot

            **fs** : font size of labels

        :Example: ::

            >>> old = occam.Occam1d()
            >>> old.savepath = r"/home/Occam1D/Line1/Inv1_TE"
            >>> old.plotL2Curve()
        """

        if savepath == None:
            if not self.savepath:
                raise IOError('No savepath found, please enter one.')

        else:
            self.savepath = savepath

        self.rms_te = []
        self.rms_tm = []
        self.roughness_te = []
        self.roughness_tm = []

        # get rms and roughness from each iteration for the different modes
        if imode == 'TE':
            # get all iteration files for TE mode
            iterlstte = [os.path.join(self.savepath, itfn)
                         for itfn in os.listdir(self.savepath)
                         if itfn.find('TE') > 0 and itfn.find('iter') > 0]

            self.rms_te = np.zeros(len(iterlstte))
            self.roughness_te = np.zeros(len(iterlstte))

            # get rms and roughness
            for itfn in iterlstte:
                self.read1DIterFile(itfn, imode='TE')

                # get iteration number to make sure the items are in sequence
                ii = int(self.itdict['Iteration'])

                # put the values in appropriate place
                self.rms_te[ii] = float(self.itdict['Misfit Value'])
                self.roughness_te[ii] = float(self.itdict['Roughness Value'])

        elif imode == 'TM':
            # get all iteration files for TM mode
            iterlsttm = [os.path.join(self.savepath, itfn)
                         for itfn in os.listdir(self.savepath)
                         if itfn.find('TM') > 0 and itfn.find('iter') > 0]

            self.rms_tm = np.zeros(len(iterlsttm))
            self.roughness_tm = np.zeros(len(iterlsttm))

            # get rms and roughness
            for itfn in iterlsttm:
                self.read1DIterFile(itfn, imode='TM')

                # get iteration number to make sure the items are in sequence
                ii = int(self.itdict['Iteration'])

                # put the values in appropriate place
                self.rms_tm[ii] = float(self.itdict['Misfit Value'])
                self.roughness_tm[ii] = float(self.itdict['Roughness Value'])

        elif imode == 'both':
            # get all iteration files for TE mode
            iterlstte = [os.path.join(self.savepath, itfn)
                         for itfn in os.listdir(self.savepath)
                         if itfn.find('TE') > 0 and itfn.find('iter') > 0]

            self.rms_te = np.zeros(len(iterlstte))
            self.roughness_te = np.zeros(len(iterlstte))

            # get rms and roughness
            for itfn in iterlstte:
                self.read1DIterFile(itfn, imode='TE')

                # get iteration number to make sure the items are in sequence
                ii = int(self.itdict['Iteration'])

                # put the values in appropriate place
                self.rms_te[ii] = float(self.itdict['Misfit Value'])
                self.roughness_te[ii] = float(self.itdict['Roughness Value'])

            # get all iteration files for TM mode
            iterlsttm = [os.path.join(self.savepath, itfn)
                         for itfn in os.listdir(self.savepath)
                         if itfn.find('TM') > 0 and itfn.find('iter') > 0]

            self.rms_tm = np.zeros(len(iterlsttm))
            self.roughness_tm = np.zeros(len(iterlsttm))

            # get rms and roughness
            for itfn in iterlsttm:
                self.read1DIterFile(itfn, imode='TM')

                # get iteration number to make sure the items are in sequence
                ii = int(self.itdict['Iteration'])

                # put the values in appropriate place
                self.rms_tm[ii] = float(self.itdict['Misfit Value'])
                self.roughness_tm[ii] = float(self.itdict['Roughness Value'])

        # plot the rms vs iteration, roughness vs rms
        #---------plot TE mode-------------------
        if imode == 'TE':
            fig = plt.figure(fignum, dpi=dpi)
            plt.clf()
            ax1 = fig.add_subplot(1, 1, 1)

            nr = len(self.rms_te)
            # plot the rms vs iteration
            l1, = ax1.plot(np.arange(1, nr, 1), self.rms_te[1:], '-k', lw=1,
                           marker='d', ms=5)

            # plot the median of the RMS
            medte = np.median(self.rms_te[1:])
            m1, = ax1.plot(np.arange(0, nr, 1),
                           np.repeat(medte, nr),
                           '--r', lw=.75)

            # make subplot for RMS vs Roughness Plot
            ax2 = ax1.twiny()

            # plot the rms vs roughness
            l2, = ax2.plot(self.roughness_te[1:], self.rms_te[1:],
                           '--b', lw=.75, marker='o', ms=7, mfc='white')
            for ii, rms in enumerate(self.rms_te[1:], 1):
                ax2.text(self.roughness_te[ii], rms, '{0}'.format(ii),
                         horizontalalignment='center',
                         verticalalignment='center',
                         fontdict={'size': 6, 'weight': 'bold', 'color': 'blue'})

            # make a legend
            ax1.legend([l1, l2, m1], ['RMS_TE', 'Roughness_TE',
                                      'Median_RMS={0:.2f}'.format(medte)],
                       ncol=4, loc='upper center', columnspacing=.25,
                       markerscale=.75, handletextpad=.15)

            ax1.set_ylim(medte - 1, medte + 1)
            ax1.set_ylabel('RMS', fontdict={'size': fs, 'weight': 'bold'})
            ax1.set_xlabel('Iteration', fontdict={
                           'size': fs, 'weight': 'bold'})
            ax1.grid(alpha=.25, which='both')
            ax2.set_xlabel('Roughness', fontdict={'size': fs, 'weight': 'bold',
                                                  'color': 'blue'})
            for t2 in ax2.get_xticklabels():
                t2.set_color('blue')

        #-------Plot TM mode-------------------
        elif imode == 'TM':
            fig = plt.figure(fignum, dpi=dpi)
            plt.clf()
            ax1 = fig.add_subplot(1, 1, 1)

            nr = len(self.rms_tm)
            # plot the rms vs iteration
            l1, = ax1.plot(np.arange(1, nr, 1), self.rms_tm[1:], '-k', lw=1,
                           marker='d', ms=5)

            # plot the median of the RMS
            medtm = np.median(self.rms_tm[1:])
            m1, = ax1.plot(np.arange(0, nr, 1),
                           np.repeat(medtm, nr),
                           '--r', lw=.75)

            # make subplot for RMS vs Roughness Plot
            ax2 = ax1.twiny()

            # plot the rms vs roughness
            l2, = ax2.plot(self.roughness_tm[1:], self.rms_tm[1:],
                           '--b', lw=.75, marker='o', ms=7, mfc='white')
            for ii, rms in enumerate(self.rms_tm[1:], 1):
                ax2.text(self.roughness_tm[ii], rms, '{0}'.format(ii),
                         horizontalalignment='center',
                         verticalalignment='center',
                         fontdict={'size': fs - 2, 'weight': 'bold', 'color': 'blue'})

            # make a legend
            ax1.legend([l1, l2, m1], ['RMS_TM', 'Roughness_TM',
                                      'Median_RMS={0:.2f}'.format(medtm)],
                       ncol=4, loc='upper center', columnspacing=.25,
                       markerscale=.75, handletextpad=.15)

            ax1.set_ylim(medtm - 1, medtm + 1)
            ax1.set_ylabel('RMS', fontdict={'size': fs, 'weight': 'bold'})
            ax1.set_xlabel('Iteration', fontdict={
                           'size': fs, 'weight': 'bold'})
            ax1.grid(alpha=.25, which='both')
            ax2.set_xlabel('Roughness', fontdict={'size': fs, 'weight': 'bold',
                                                  'color': 'blue'})
            for t2 in ax2.get_xticklabels():
                t2.set_color('blue')

        elif imode == 'both':
            fig = plt.figure(fignum, dpi=dpi)
            plt.clf()
            ax1 = fig.add_subplot(2, 1, 1)
            ax3 = fig.add_subplot(2, 1, 2, sharex=ax1)

            plt.rcParams['figure.subplot.hspace'] = .4
            plt.rcParams['figure.subplot.left'] = .1
            plt.rcParams['figure.subplot.right'] = .97
            plt.rcParams['figure.subplot.bottom'] = .1
            plt.rcParams['figure.subplot.top'] = .92

            nr = len(self.rms_te)
            # plot the rms vs iteration
            l1, = ax1.plot(np.arange(1, nr, 1), self.rms_te[1:], '-k', lw=1,
                           marker='d', ms=5)

            # plot the median of the RMS
            medte = np.median(self.rms_te[1:])
            m1, = ax1.plot(np.arange(0, nr, 1),
                           np.repeat(medte, nr),
                           '--r', lw=.75)

            # make subplot for RMS vs Roughness Plot
            ax2 = ax1.twiny()

            # plot the rms vs roughness
            l2, = ax2.plot(self.roughness_te[1:], self.rms_te[1:],
                           '--b', lw=.75, marker='o', ms=7, mfc='white')
            for ii, rms in enumerate(self.rms_te[1:], 1):
                ax2.text(self.roughness_te[ii], rms, '{0}'.format(ii),
                         horizontalalignment='center',
                         verticalalignment='center',
                         fontdict={'size': fs - 2, 'weight': 'bold', 'color': 'blue'})

            # make a legend
            ax1.legend([l1, l2, m1], ['RMS_TE', 'Roughness_TE',
                                      'Median_RMS={0:.2f}'.format(medte)],
                       ncol=4, loc='upper center', columnspacing=.25,
                       markerscale=.75, handletextpad=.15)

            ax1.set_ylim(medte - 1, medte + 1)
            ax1.set_ylabel('RMS', fontdict={'size': fs, 'weight': 'bold'})
            # ax1.set_xlabel('Iteration',fontdict={'size':8,'weight':'bold'})
            ax1.grid(alpha=.25, which='both')
            ax2.set_xlabel('Roughness', fontdict={'size': fs, 'weight': 'bold',
                                                  'color': 'blue'})
            for t2 in ax2.get_xticklabels():
                t2.set_color('blue')

            # plot TM
            nr = len(self.rms_te)
            # plot the rms vs iteration
            l3, = ax3.plot(np.arange(1, nr, 1), self.rms_tm[1:], '-k', lw=1,
                           marker='d', ms=5)

            # plot the median of the RMS
            medtm = np.median(self.rms_tm[1:])
            m3, = ax3.plot(np.arange(0, nr, 1),
                           np.repeat(medtm, nr),
                           '--r', lw=.75)

            # make subplot for RMS vs Roughness Plot
            ax4 = ax3.twiny()

            # plot the rms vs roughness
            l4, = ax4.plot(self.roughness_tm[1:], self.rms_tm[1:],
                           '--b', lw=.75, marker='o', ms=7, mfc='white')
            for ii, rms in enumerate(self.rms_tm[1:], 1):
                ax4.text(self.roughness_tm[ii], rms, '{0}'.format(ii),
                         horizontalalignment='center',
                         verticalalignment='center',
                         fontdict={'size': 6, 'weight': 'bold', 'color': 'blue'})

            # make a legend
            ax3.legend([l1, l2, m1], ['RMS_TM', 'Roughness_TM',
                                      'Median_RMS={0:.2f}'.format(medtm)],
                       ncol=4, loc='upper center', columnspacing=.25,
                       markerscale=.75, handletextpad=.15)

            ax3.set_ylim(medtm - 1, medtm + 1)
            ax3.set_ylabel('RMS', fontdict={'size': fs, 'weight': 'bold'})
            ax3.set_xlabel('Iteration', fontdict={
                           'size': fs, 'weight': 'bold'})
            ax3.grid(alpha=.25, which='both')
            ax4.set_xlabel('Roughness', fontdict={'size': fs, 'weight': 'bold',
                                                  'color': 'blue'})
            for t2 in ax4.get_xticklabels():
                t2.set_color('blue')

        plt.show()


def getdatetime():

    return time.asctime(time.gmtime())


def makestartfiles(parameter_dict):

    read_datafile(parameter_dict)

    parameter_dict['n_sideblockelements'] = 7
    parameter_dict['n_bottomlayerelements'] = 4

    parameter_dict['itform'] = 'not specified'
    parameter_dict['description'] = 'N/A'

    parameter_dict['datetime'] = getdatetime()

    parameter_dict['iruf'] = 1
    parameter_dict['idebug'] = 1
    parameter_dict['nit'] = 0
    parameter_dict['pmu'] = 5.0
    parameter_dict['rlast'] = 1.0E+07
    parameter_dict['tobt'] = 100.
    parameter_dict['ifftol'] = 0

    blocks_elements_setup(parameter_dict)

    get_model_setup(parameter_dict)

    writemeshfile(parameter_dict)
    writemodelfile(parameter_dict)
    writestartupfile(parameter_dict)

    MeshF = parameter_dict['meshfn']
    ModF = parameter_dict['inmodelfn']
    SF = parameter_dict['startupfn']

    return (MeshF, ModF, SF)


def writemeshfile(parameter_dict):

    mesh_positions_vert = parameter_dict['mesh_positions_vert']
    mesh_positions_hor = parameter_dict['mesh_positions_hor']
    n_nodes_hor = parameter_dict['n_nodes_hor']
    n_nodes_vert = parameter_dict['n_nodes_vert']

    fh_mesh = file(parameter_dict['meshfn'], 'w')
    mesh_outstring = ''

    temptext = "MESH FILE FROM MTpy\n"
    mesh_outstring += temptext

    temptext = "%i %i %i %i %i %i\n" % (0, n_nodes_hor, n_nodes_vert, 0, 0, 2)
    mesh_outstring += temptext

    temptext = ""
    for i in range(n_nodes_hor - 1):
        temptext += "%.1f " % (mesh_positions_hor[i])
    temptext += "\n"
    mesh_outstring += temptext

    temptext = ""
    for i in range(n_nodes_vert - 1):
        temptext += "%.1f " % (mesh_positions_vert[i])
    temptext += "\n"
    mesh_outstring += temptext

    mesh_outstring += "%i\n" % (0)

    for j in range(4 * (n_nodes_vert - 1)):
        tempstring = ''
        tempstring += (n_nodes_hor - 1) * "?"
        tempstring += '\n'
        mesh_outstring += tempstring

    fh_mesh.write(mesh_outstring)
    fh_mesh.close()


def writemodelfile(parameter_dict):
    "needed : filename,binding_offset,startcolumn, n_layers,layer_thickness,block_width"

    modelblockstrings = parameter_dict['modelblockstrings']
    nfev = parameter_dict['nfev']
    lo_colnumbers = parameter_dict['lo_colnumbers']
    boffset = float(parameter_dict['binding_offset'])
    n_layers = int(float(parameter_dict['n_layers']))

    fh_model = file(parameter_dict['inmodelfn'], 'w')
    model_outstring = ''

    temptext = "Format:           %s\n" % ("OCCAM2MTMOD_1.0")
    model_outstring += temptext
    temptext = "Model Name:       %s\n" % (parameter_dict['modelname'])
    model_outstring += temptext
    temptext = "Description:      %s\n" % ("Random Text")
    model_outstring += temptext
    temptext = "Mesh File:        %s\n" % (
        os.path.basename(parameter_dict['meshfn']))
    model_outstring += temptext
    temptext = "Mesh Type:        %s\n" % ("PW2D")
    model_outstring += temptext
    temptext = "Statics File:     %s\n" % ("none")
    model_outstring += temptext
    temptext = "Prejudice File:   %s\n" % ("none")
    model_outstring += temptext
    temptext = "Binding Offset:   %.1f\n" % (boffset)
    model_outstring += temptext
    temptext = "Num Layers:       %i\n" % (n_layers)
    model_outstring += temptext

    for k in range(n_layers):
        n_meshlayers = nfev[k]
        n_meshcolumns = lo_colnumbers[k]
        temptext = "%i %i\n" % (n_meshlayers, n_meshcolumns)
        model_outstring += temptext

        temptext = modelblockstrings[k]
        model_outstring += temptext
        #model_outstring += "\n"

    temptext = "Number Exceptions:%i\n" % (0)
    model_outstring += temptext

    fh_model.write(model_outstring)
    fh_model.close()


def writestartupfile(parameter_dict):

    fh_startup = file(parameter_dict['startupfn'], 'w')
    startup_outstring = ''

    temptext = "Format:           %s\n" % (parameter_dict['itform'])
    startup_outstring += temptext
    temptext = "Description:      %s\n" % (parameter_dict['description'])
    startup_outstring += temptext
    temptext = "Model File:       %s\n" % (
        os.path.basename(parameter_dict['inmodelfn']))
    startup_outstring += temptext
    temptext = "Data File:        %s\n" % (
        os.path.basename(parameter_dict['datafile']))
    startup_outstring += temptext
    temptext = "Date/Time:        %s\n" % (parameter_dict['datetime'])
    startup_outstring += temptext
    temptext = "Max Iter:         %i\n" % (
        int(float(parameter_dict['n_max_iterations'])))
    startup_outstring += temptext
    temptext = "Req Tol:          %.1g\n" % (
        float(parameter_dict['targetrms']))
    startup_outstring += temptext
    temptext = "IRUF:             %s\n" % (parameter_dict['iruf'])
    startup_outstring += temptext
    temptext = "Debug Level:      %s\n" % (parameter_dict['idebug'])
    startup_outstring += temptext
    temptext = "Iteration:        %i\n" % (
        int(float(parameter_dict['n_max_iterations'])))
    startup_outstring += temptext
    temptext = "PMU:              %s\n" % (parameter_dict['pmu'])
    startup_outstring += temptext
    temptext = "Rlast:            %s\n" % (parameter_dict['rlast'])
    startup_outstring += temptext
    temptext = "Tlast:            %s\n" % (parameter_dict['tobt'])
    startup_outstring += temptext
    temptext = "IffTol:           %s\n" % (parameter_dict['ifftol'])
    startup_outstring += temptext
    temptext = "No. Parms:        %i\n" % (
        int(float(parameter_dict['n_parameters'])))
    startup_outstring += temptext
    temptext = ""
    for l in range(int(float(parameter_dict['n_parameters']))):
        temptext += "%.1g  " % (2.0)
    temptext += "\n"
    startup_outstring += temptext

    fh_startup.write(startup_outstring)
    fh_startup.close()


def read_datafile(parameter_dict):

    df = parameter_dict['datafile']
    F = file(df, 'r')
    datafile_content = F.readlines()
    F.close()

    # RELYING ON A CONSTANT FORMAT, ACCESSING THE PARTS BY COUNTING OF
    # LINES!!!:

    n_sites = int(datafile_content[2].strip().split()[1])
    sitenames = []
    for i in range(n_sites):
        sitenames.append(datafile_content[3 + i].strip())

    sitelocations = []
    for i in range(n_sites):
        idx = 4 + n_sites + i
        sitelocations.append(float(datafile_content[idx].strip()))

    n_freqs = int(datafile_content[2 * n_sites + 4].strip().split()[1])
    freqs = []
    for i in range(n_freqs):
        idx = 2 * n_sites + 5 + i
        freqs.append(float(datafile_content[idx].strip()))

    n_data = int(datafile_content[2 * n_sites +
                                  5 + n_freqs].strip().split()[2])

    parameter_dict['lo_site_names'] = sitenames
    parameter_dict['lo_site_locations'] = sitelocations
    parameter_dict['n_sites'] = n_sites
    parameter_dict['n_datapoints'] = n_data
    parameter_dict['n_freqs'] = n_freqs
    parameter_dict['lo_freqs'] = freqs


def get_model_setup(parameter_dict):

    ncol0 = int(float(parameter_dict['ncol0']))
    n_layer = int(float(parameter_dict['n_layers']))
    nfe = parameter_dict['nfe']
    thickness = parameter_dict['thickness']
    width = parameter_dict['width']
    trigger = float(parameter_dict['trigger'])
    dlz = parameter_dict['dlz']

    modelblockstrings = []
    lo_colnumbers = []

    ncol = ncol0
    np = 0

    for layer_idx in range(n_layer):
        block_idx = 1

        # print layer_idx,len(thickness),len(width),ncol,len(dlz)

        while block_idx + 2 < ncol - 1:

            # PROBLEM : 'thickness' has only "n_layer'-1 entries!!
            if not dlz[layer_idx] > (trigger * (width[block_idx] + width[block_idx + 1])):
                block_idx += 1
                continue

            else:
                width[block_idx] += width[block_idx + 1]
                nfe[block_idx] += nfe[block_idx + 1]

                for m in range(block_idx + 2, ncol):
                    width[m - 1] = width[m]
                    nfe[m - 1] = nfe[m]

                ncol -= 1

        lo_colnumbers.append(ncol)

        tempstring = ""
        for j in range(ncol):
            tempstring += "%i " % (nfe[j])
        tempstring += "\n"
        modelblockstrings.append(tempstring)

        np = np + ncol

        # completely unnecessary!!! :
        if layer_idx == 0:
            mcol = ncol

    parameter_dict['modelblockstrings'] = modelblockstrings
    parameter_dict['lo_colnumbers'] = lo_colnumbers
    parameter_dict['n_parameters'] = np
    parameter_dict['n_cols_max'] = mcol


def blocks_elements_setup(parameter_dict):

    lo_sites = parameter_dict['lo_site_locations']
    n_sites = len(lo_sites)
    maxwidth = float(parameter_dict['max_blockwidth'])

    nbot = int(float(parameter_dict['n_bottomlayerelements']))
    nside = int(float(parameter_dict['n_sideblockelements']))

    # j: index for finite elements
    # k: index for regularisation bricks
    # Python style: start with 0 instead of 1

    sitlok = []
    sides = []
    width = []
    dly = []
    nfe = []
    thickness = []
    nfev = []
    dlz = []
    bot = []

    j = 0
    sitlok.append(lo_sites[0])

    for idx in range(1, n_sites - 1):

        spacing = lo_sites[idx] - lo_sites[idx - 1]
        n_localextrasites = int(spacing / maxwidth) + 1

        for idx2 in range(n_localextrasites):
            sitlok.append(lo_sites[idx - 1] + (idx2 + 1.) /
                          float(n_localextrasites) * spacing)
            j += 1

    # nrk: number of total dummy stations
    nrk = j
    print("%i dummy stations defined" % (nrk))

    spacing1 = (sitlok[1] - sitlok[0]) / 2.
    sides.append(3 * spacing1)

    for idx in range(1, nside):
        curr_side = 3 * sides[idx - 1]
        if curr_side > 1000000.:
            curr_side = 1000000.
        sides.append(curr_side)

    #-------------------------------------------

    j = 0
    k = 0

    firstblockwidth = 0.

    for idx in range(nside - 1, -1, -1):
        firstblockwidth += sides[idx]
        dly.append(sides[idx])
        j += 1

    width.append(firstblockwidth)
    nfe.append(nside)

    dly.append(spacing1)
    dly.append(spacing1)
    j += 2
    nfe.append(2)
    width.append(2 * spacing1)

    block_offset = width[1]

    k += 1

    dly.append(spacing1)
    dly.append(spacing1)
    j += 2
    nfe.append(2)
    width.append(2 * spacing1)

    block_offset += spacing1

    k += 1

    #------------------------

    for idx in range(1, nrk - 1):
        spacing2 = (sitlok[idx + 1] - sitlok[idx]) / 2.
        dly.append(spacing1)
        dly.append(spacing2)
        j += 2
        nfe.append(2)
        width.append(spacing1 + spacing2)
        k += 1
        spacing1 = spacing2

    dly.append(spacing2)
    dly.append(spacing2)

    j += 2
    nfe.append(2)
    width.append(2 * spacing2)
    k += 1

    dly.append(spacing2)
    dly.append(spacing2)

    j += 2
    nfe.append(2)
    width.append(2 * spacing2)
    k += 1

    width[-1] = 0.
    sides[0] = 3 * spacing2

    #------------------------------

    for idx in range(1, nside):
        curr_side = 3 * sides[idx - 1]
        if curr_side > 1000000.:
            curr_side = 1000000.
        sides[idx] = curr_side

    lastblockwidth = 0.
    for idx in range(nside):
        j += 1
        lastblockwidth += sides[idx]
        dly.append(sides[idx])

    width[-1] = lastblockwidth

    #---------------------------------

    k += 1
    nfe.append(nside)

    nodey = j + 1
    ncol0 = k

    block_offset = sitlok[0] - block_offset

    #----------------------------------

    layers_per_decade = float(parameter_dict['n_layersperdecade'])
    first_layer_thickness = float(parameter_dict['firstlayer_thickness'])

    t = 10.**(1. / layers_per_decade)
    t1 = first_layer_thickness
    thickness.append(t1)

    d1 = t1

    n_layers = int(float(parameter_dict['n_layers']))

    for idx in range(1, n_layers - 1):
        d2 = d1 * t
        curr_thickness = d2 - d1
        if curr_thickness < t1:
            curr_thickness = t1
        thickness.append(curr_thickness)
        d1 += curr_thickness

    bot.append(3 * thickness[n_layers - 2])

    for idx in range(1, nbot):
        bot.append(bot[idx - 1] * 3)

    #--------------------------------------------------

    k = 0

    dlz.append(thickness[0] / 2.)
    dlz.append(thickness[0] / 2.)
    nfev.append(2)

    k += 2

    dlz.append(thickness[1] / 2.)
    dlz.append(thickness[1] / 2.)
    nfev.append(2)

    k += 2

    for idx in range(2, n_layers - 1):
        k += 1
        nfev.append(1.)
        dlz.append(thickness[idx])

    for idx in range(nbot):
        k += 1
        dlz.append(bot[idx])

    nfev.append(nbot)

    nodez = k + 1

    parameter_dict['ncol0'] = ncol0
    parameter_dict['nfe'] = nfe
    parameter_dict['nfev'] = nfev
    parameter_dict['thickness'] = thickness
    parameter_dict['width'] = width
    parameter_dict['binding_offset'] = block_offset
    #parameter_dict['y_nodes']           = nodey
    #parameter_dict['z_nodes']           = nodez
    parameter_dict['dlz'] = dlz
    #parameter_dict['dly']               = dly
    parameter_dict['mesh_positions_vert'] = dlz
    parameter_dict['mesh_positions_hor'] = dly
    parameter_dict['n_nodes_hor'] = nodey
    parameter_dict['n_nodes_vert'] = nodez


class OccamPointPicker(object):
    """
    This class helps the user interactively pick points to mask and add 
    error bars. 

    Useage:
    -------
    To mask just a single point right click over the point and a gray point 
    will appear indicating it has been masked

    To mask both the apparent resistivity and phase left click over the point.
    Gray points will appear over both the apparent resistivity and phase.  
    Sometimes the points don't exactly matchup, haven't quite worked that bug
    out yet, but not to worry it picks out the correct points

    To add error bars to a point click the middle or scroll bar button.  This
    only adds error bars to the point and does not reduce them so start out
    with reasonable errorbars.  You can change the increment that the error
    bars are increased with reserrinc and phaseerrinc.

    Arguments:
    ----------
        **axlst** : list of the resistivity and phase axis that have been 
                    plotted as [axr_te,axr_tm,axp_te,axp_tm]

        **linelst** : list of lines used to plot the responses, not the error 
                      bars as [res_te,res_tm,phase_te,phase_tm]

        **errlst** : list of the errorcaps and errorbar lines as 
                   [[cap1,cap2,bar],...]

        **reserrinc** : increment to increase the errorbars for resistivity.
                        put .20 for 20 percent change. *Default* is .05

        **phaseerrinc** : increment to increase the errorbars for the phase
                          put .10 for 10 percent change. *Defualt* is .02 

        **marker** : marker type for masked points.  See matplotlib.pyplot.plot
                    for options of markers.  *Default* is h for hexagon.

    Attributes:
    -----------

        **axlst** : axes list used to plot the data

        **linelst** : line list used to plot the data

        **errlst** : error list used to plot the data

        **data** : list of data points that were not masked for each plot.

        **fdict** : dictionary of frequency arrays for each plot and data set.

        **fndict** : dictionary of figure numbers to corresponed with data.

        **cidlst** : list of event ids.

        **reserrinc** : increment to increase resistivity error bars

        **phaseinc** : increment to increase phase error bars

        **marker** : marker of masked points

        **fignum** : figure numbers

        **occamlines** : list of lines to write into the occam data file.

    :Example: ::

        >>> ocd = occam.Occam2DData()
        >>> ocd.datafn = r"/home/Occam2D/Line1/Inv1/Data.dat"
        >>> ocd.plotMaskPoints()
    """

    def __init__(self, axlst, linelst, errlst, reserrinc=.05, phaseerrinc=.02,
                 marker='h'):

        # give the class some attributes
        self.axlst = axlst
        self.linelst = linelst
        self.errlst = errlst
        self.data = []
        self.error = []
        self.fdict = []
        self.fndict = {}
        # see if just one figure is plotted or multiple figures are plotted
        self.ax = axlst[0][0]
        self.line = linelst[0][0]
        self.cidlst = []
        for nn in range(len(axlst)):
            self.data.append([])
            self.error.append([])
            self.fdict.append([])

            # get data from lines and make a dictionary of frequency points for
            # easy indexing
            for ii, line in enumerate(linelst[nn]):
                self.data[nn].append(line.get_data()[1])
                self.fdict[nn].append(dict([('{0:.5g}'.format(kk), ff) for ff, kk in
                                            enumerate(line.get_data()[0])]))
                self.fndict['{0}'.format(line.figure.number)] = nn

                # set some events
                if ii == 0:
                    cid1 = line.figure.canvas.mpl_connect('pick_event', self)
                    cid2 = line.figure.canvas.mpl_connect('axes_enter_event',
                                                          self.inAxes)
                    cid3 = line.figure.canvas.mpl_connect('key_press_event',
                                                          self.on_close)
                    cid4 = line.figure.canvas.mpl_connect('figure_enter_event',
                                                          self.inFigure)
                    self.cidlst.append([cid1, cid2, cid3, cid4])

            # read in the error in a useful way so that it can be translated to
            # the data file.  Make the error into an array
            for ee, err in enumerate(errlst[nn]):
                errpath = err[2].get_paths()
                errarr = np.zeros(len(list(self.fdict[nn][ee].keys())))
                for ff, epath in enumerate(errpath):
                    errv = epath.vertices
                    errarr[ff] = abs(errv[0, 1] - self.data[nn][ee][ff])
                self.error[nn].append(errarr)

        # set the error bar increment values
        self.reserrinc = reserrinc
        self.phaseerrinc = phaseerrinc

        # set the marker
        self.marker = marker

        # set the figure number
        self.fignum = self.line.figure.number

        # make a list of occam lines to write later
        self.occamlines = []

    def __call__(self, event):
        """
        When the function is called the mouse events will be recorder for 
        picking points to mask or change error bars.  The axes is redrawn with
        a gray marker to indicate a masked point and/or increased size in 
        errorbars.

        Arguments:
        ----------
            **event** : type mouse_click_event

        Useage:
        -------

            **Left mouse button** will mask both resistivity and phase point

            **Right mouse button** will mask just the point selected

            **Middle mouse button** will increase the error bars

            **q** will close the figure.
        """
        self.event = event
        # make a new point that is an PickEvent type
        npoint = event.artist
        # if the right button is clicked mask the point
        if event.mouseevent.button == 3:
            # get the point that was clicked on
            ii = event.ind
            xd = npoint.get_xdata()[ii]
            yd = npoint.get_ydata()[ii]

            # set the x index from the frequency dictionary
            ll = self.fdict[self.fignum][self.jj]['{0:.5g}'.format(xd[0])]

            # change the data to be a zero
            self.data[self.fignum][self.jj][ll] = 0

            # reset the point to be a gray x
            self.ax.plot(xd, yd, ls='None', color=(.7, .7, .7), marker=self.marker,
                         ms=4)

        # if the left button is clicked change both resistivity and phase
        # points
        elif event.mouseevent.button == 1:
            # get the point that was clicked on
            ii = event.ind
            xd = npoint.get_xdata()[ii]
            yd = npoint.get_ydata()[ii]

            # set the x index from the frequency dictionary
            ll = self.fdict[self.fignum][self.jj]['{0:.5g}'.format(xd[0])]

            # set the data point to zero
            self.data[self.fignum][self.jj][ll] = 0

            # reset the point to be a gray x
            self.ax.plot(xd, yd, ls='None', color=(.7, .7, .7), marker=self.marker,
                         ms=4)

            # check to make sure there is a corresponding res/phase point
            try:
                # get the corresponding y-value
                yd2 = self.data[self.fignum][self.kk][ll]

                # set that data point to 0 as well
                self.data[self.fignum][self.kk][ll] = 0

                # make that data point a gray x
                self.axlst[self.fignum][self.kk].plot(xd, yd2, ls='None',
                                                      color=(.7, .7, .7), marker=self.marker,
                                                      ms=4)
            except KeyError:
                print('Axis does not contain res/phase point')

        # if click the scroll button or middle button change increase the
        # errorbars by the given amount
        elif event.mouseevent.button == 2:
            ii = event.ind
            xd = npoint.get_xdata()[ii]
            yd = npoint.get_ydata()[ii]

            # get x index
            ll = self.fdict[self.fignum][self.jj]['{0:.5g}'.format(xd[0])]

            # make error bar array
            eb = self.errlst[self.fignum][self.jj][2].get_paths()[ll].vertices

            # make ecap array
            ecapl = self.errlst[self.fignum][self.jj][0].get_data()[1][ll]
            ecapu = self.errlst[self.fignum][self.jj][1].get_data()[1][ll]

            # change apparent resistivity error
            if self.jj == 0 or self.jj == 1:
                nebu = eb[0, 1] - self.reserrinc * eb[0, 1]
                nebl = eb[1, 1] + self.reserrinc * eb[1, 1]
                ecapl = ecapl - self.reserrinc * ecapl
                ecapu = ecapu + self.reserrinc * ecapu

            # change phase error
            elif self.jj == 2 or self.jj == 3:
                nebu = eb[0, 1] - eb[0, 1] * self.phaseerrinc
                nebl = eb[1, 1] + eb[1, 1] * self.phaseerrinc
                ecapl = ecapl - ecapl * self.phaseerrinc
                ecapu = ecapu + ecapu * self.phaseerrinc

            # put the new error into the error array
            self.error[self.fignum][self.jj][ll] = abs(nebu -
                                                       self.data[self.fignum][self.jj][ll])

            # set the new error bar values
            eb[0, 1] = nebu
            eb[1, 1] = nebl

            # reset the error bars and caps
            ncapl = self.errlst[self.fignum][self.jj][0].get_data()
            ncapu = self.errlst[self.fignum][self.jj][1].get_data()
            ncapl[1][ll] = ecapl
            ncapu[1][ll] = ecapu

            # set the values
            self.errlst[self.fignum][self.jj][0].set_data(ncapl)
            self.errlst[self.fignum][self.jj][1].set_data(ncapu)
            self.errlst[self.fignum][self.jj][2].get_paths()[ll].vertices = eb

        # redraw the canvas
        self.ax.figure.canvas.draw()

    # get the axis number that the mouse is in and change to that axis
    def inAxes(self, event):
        """
        gets the axes that the mouse is currently in.

        Arguments:
        ---------
            **event**: is a type axes_enter_event

        Returns:
        --------

            **OccamPointPicker.jj** : index of resistivity axes for axlst

            **OccamPointPicker.kk** : index of phase axes for axlst

        """

        self.event2 = event
        self.ax = event.inaxes
        for jj, axj in enumerate(self.axlst):
            for ll, axl in enumerate(axj):
                if self.ax == axl:
                    self.jj = ll

        # set complimentary resistivity and phase plots together
        if self.jj == 0:
            self.kk = 2
        if self.jj == 1:
            self.kk = 3
        if self.jj == 2:
            self.kk = 0
        if self.jj == 3:
            self.kk = 1

    # get the figure number that the mouse is in
    def inFigure(self, event):
        """
        gets the figure number that the mouse is in

        Arguments:
        ----------
            **event** : figure_enter_event

        Returns:
        --------
            **OccamPointPicker.fignum** : figure number that corresponds to the
                                          index in the axlst, datalst, errorlst
                                          and linelst.

        """
        self.event3 = event
        self.fignum = self.fndict['{0}'.format(event.canvas.figure.number)]
        self.line = self.linelst[self.fignum][0]

    # type the q key to quit the figure and disconnect event handling
    def on_close(self, event):
        """
        close the figure with a 'q' key event and disconnect the event ids

        Arguments:
        ----------
            **event** : key_press_event

        Returns:
        --------
            print statement saying the figure is closed
        """
        self.event3 = event
        if self.event3.key == 'q':
            for cid in self.cidlst[self.fignum]:
                event.canvas.mpl_disconnect(cid)
            plt.close(event.canvas.figure)
            print('Closed figure ', self.fignum)


class Occam2DData:
    """
    Occam2DData covers all aspects of dealing with data for an Occam 2D
    inversion using the code of Constable et al. [1987] and deGroot-Hedlin and 
    Constable [1990] from Scripps avaliable at 
    http://marineemlab.ucsd.edu/Projects/Occam/2DMT/index.html.


    """

    def __init__(self, datafn=None):
        self.datafn = datafn

    def make2DdataFile(self, edipath, mmode='both', savepath=None,
                       stationlst=None, title=None, thetar=0, resxyerr=10,
                       resyxerr=10, phasexyerr=5, phaseyxerr=5, ss=3 * ' ',
                       string_fmt='%+2.6f', freqstep=1, plotyn='y',
                       lineori='ew', proj_strike='yes', tipper_err=None,
                       ftol=.05):
        """
        Make a data file that Occam can read.  At the moment the inversion line
        is the best fit line through all the stations used for the inversion.

        Arguments:
        ----------
            **edipath** : path to edifiles

            **mmode** : modes to invert for.  Can be: 
                        * 'both' -> will model both TE and TM modes
                        * 'TM'   -> will model just TM mode
                        * 'TE'   -> will model just TE mode

            **savepath** : path to save the data file to, this can include the 
                           name of the data file, if not the file will be 
                           named: savepath/Data.dat or edipath/Data.dat if 
                           savepath=None

            **stationlst** : list of stations to put in the data file, doesn't 
                             need to be in order, the relative distance will be
                             calculated internally.  If stationlst:None, it
                             will be assumed all the files in edipath will be 
                             input into the data file

            **title** : title input into the data file. *Default* is None

            **thetar** : rotation angle (deg) of the MT response to align 
                         TE mode along strike which is assumed to be 
                         perpendicular to the inversion line and TM
                         parallel to the inversion line.  Angle is on 
                         the unit circle with an orientation that north is 0 
                         degree, east 90. Rotations are postive clockwise.
                         *Default* is 0 (North).  

                         If proj_strike = 'yes' then the stations will be 
                         projected onto a line that is perpendicular to 
                         the thetar, which is assumed to be the strike
                         direction. 

            **resxyerr** : percent error in the res_xy component (TE), 
                          can be entered as 'data' where the errors from the 
                          data are used otherwise enter as a percentage.
                          enter 10 for 10 percent. *Default* is 10

            **resyxerr** : percent error in the res_yx component (TM), 
                          can be entered as 'data' where the errors from the 
                          data are used otherwise enter as a percentage.
                          enter 10 for 10 percent. *Default* is 10

            **phasexyerr** : percent error in the phase_xy component (TE), 
                             can be entered as 'data' where the errors from the 
                             data are used otherwise enter as a percentage.
                             enter 10 for 10 percent. *Default* is 5  

            **phaseyxerr** : percent error in the phase_yx component (TM), 
                             can be entered as 'data' where the errors from the 
                             data are used otherwise enter as a percentage.
                             enter 10 for 10 percent. *Default* is 5  

            **ss** : string
                    is the spacing parameter for the data file. 
                     *Default* is ' '*3 (3 spaces) 

            **string_fmt** : format of the numbers for the data file, see string 
                      formats for a full description. *Default* is '%+2.6f

            **freqstep** : take frequencies at this step, so if you want to 
                           take every third frequency enter 3.  
                           Can input as a list of specific frequencies.  
                           Note that the frequencies must match the frequencies
                           in the EDI files within ftol, otherwise they will 
                           not be input.  

            **ftol** : tolerance level (decimal %) to match frequencies to 
                       freqstep if input as a list.  *Default* is .05

            **plotyn** : y or n to plot the stations that are pojected on to 
                         the best fitting line through the stations.

            **lineori** : predominant line orientation with respect to 
                          geographic north
                          ew for east-west line-> will orientate so first 
                                                  station is farthest to the 
                                                  west
                          ns for north-south line-> will orientate so first 
                                                    station is farthest to the
                                                    south

            **proj_strike** : 'yes' to project the line perpendicular to the 
                                    strike direction
                              'no' to project the line on the best fitting line
                                   through the stations.

            **tipper_err** : error for tipper in percent.  If this value is 
                            entered than the tipper will be included in the 
                            inversion, if the value is None than the tipper 
                            will not be included. 
                            Can be entered as 'data' where the errors from the 
                            data are used otherwise enter as a percentage.
                            enter 10 for 10 percent.

        Returns:
        --------
            **Occam2DData.datafn** : full path of data file

        :Example: ::

            >>> import mtpy.core.occamtools as occam
            >>> ocd = occam.Occam2DData()
            >>> #define a path to where the edifiles are
            >>> edipath = r"/home/EDIfiles"
            >>>
            >>> #create a save path to put the data file
            >>> svpath = r"/home/Occam2D/Line1/Inv1"
            >>>
            >>> #create a station list of stations to use in inversion
            >>> slst = ['MT0{0}'.format(ii) for ii in range(1,10)]
            >>> 
            >>> #write data file that is rotated 50 degrees east of north and
            >>> #projected on to the strike direction 50 degrees east of north
            >>> ocd.make2DdataFile(edipath,stationlst=slst,savepath=svpath,
            >>> ...                thetar=50,proj_strike='yes',lineori='ew')
            >>>                       
            >>> 'Wrote Occam2D data file to: /home/Occam2D/Line1/Inv1/Data.dat' 


        """

        if abs(thetar) > 2 * np.pi:
            thetar = thetar * (np.pi / 180)

        #-----------------------Station Locations------------------------------
        # create a list to put all the station dictionaries into
        surveylst = []
        eastlst = []
        northlst = []
        pstationlst = []
        freqlst = []

        # get edi files for all stations in edipath if stationlst is None
        if stationlst == None:
            stationlst = [edifile[:-4]
                          for edifile in os.listdir(edipath) if edifile.find('.edi')]

        for kk, station in enumerate(stationlst):
            # search for filenames in the given directory and match to station
            # name
            for filename in os.listdir(edipath):
                if fnmatch.fnmatch(filename, station + '*.edi'):
                    print('Found station edifile: ', filename)

                    # create a dictionary for the station data and info
                    surveydict = {}
                    edifile = os.path.join(edipath, filename)

                    # read in edi file
                    z1 = mtedi.Edi()
                    z1.readfile(edifile)
                    freq = z1.freq

                    # rotate data
                    if thetar != 0:
                        z1.Z.rotate(thetar)
                        if z1.Tipper.tipper is not None:
                            z1.Tipper.rotate(thetar)

                    # check to see if the frequency is in descending order
                    if freq[0] < freq[-1]:
                        freq = freq[::-1]
                        z = z1.Z.z[::-1, :, :]
                        zvar = z1.Z.z_err[::-1, :, :]
                        if z1.Tipper.tipper is not None:
                            tip = z1.Tipper.tipper[::-1, :, :]
                            tipvar = z1.Tipper.tipper_err[::-1, :]

                        print(('Flipped frequency to descending for station: '
                               '{0}'.format(station)))
                    else:
                        z = z1.Z.z
                        zvar = z1.Z.z_err
                        if z1.Tipper.tipper is not None:
                            tip = z1.Tipper.tipper
                            tipvar = z1.Tipper.tipper_err

                    # get eastings and northings so everything is in meters
                    zone, east, north = mtpy.utils.gis_tools.ll_to_utm(23, z1.lat, z1.lon)

                    # put things into a dictionary to sort out order of
                    # stations
                    surveydict['station'] = station
                    surveydict['east'] = east
                    surveydict['north'] = north
                    surveydict['zone'] = zone
                    surveydict['z'] = z
                    surveydict['zvar'] = zvar
                    surveydict['freq'] = freq
                    surveydict['tipper'] = tip
                    surveydict['tippervar'] = tipvar
                    surveydict['lat'] = z1.lat
                    surveydict['lon'] = z1.lon
                    freqlst.append(freq)
                    eastlst.append(east)
                    northlst.append(north)
                    pstationlst.append(station)
                    surveylst.append(surveydict)

        self.eastlst = np.array(eastlst)
        self.northlst = np.array(northlst)
        #-----------------------------------------------------------------
        # project stations onto a best fitting line taking into account the
        # strike direction to get relative MT distances correct
        #-----------------------------------------------------------------

        # get bestfitting line
        p = sp.polyfit(self.eastlst, self.northlst, 1)

        # now project this line onto the strike direction so that the distances
        # are relative to geoelectric strike.  Needs to be negative cause
        # the strike angles is measured clockwise, where as the line angle is
        # measured counterclockwise.
        if proj_strike == 'yes':
            p[0] = thetar
        else:
            pass

        # need to work on this

        # the angle of the line is now the angle of the best fitting line added
        # to the geoelectric strike direction, which gives the relative distance
        # along the strike direction.
        theta = np.arctan(p[0])
        print('Profile Line Angle is: {0:.4g}'.format(theta * 180 / np.pi))

        # plot stations on profile line
        if plotyn == 'y':
            lfig = plt.figure(4, dpi=200)
            plt.clf()
            ploty = sp.polyval(p, self.eastlst)
            lax = lfig.add_subplot(1, 1, 1, aspect='equal')
            lax.plot(self.eastlst, ploty, '-b', lw=2)
            lax.set_title('Projected Stations')
            lax.set_ylim(ploty.min() - 1000., ploty.max() + 1000.)
            lax.set_xlim(self.eastlst.min() - 1000, self.eastlst.max() + 1000.)
            lax.set_xlabel('Easting (m)',
                           fontdict={'size': 12, 'weight': 'bold'})
            lax.set_ylabel('Northing (m)',
                           fontdict={'size': 12, 'weight': 'bold'})
            plt.show()
        for ii in range(len(surveylst)):
            if surveylst[ii]['zone'] != surveylst[0]['zone']:
                print(surveylst[ii]['station'])
            d = (northlst[ii] - sp.polyval(p,
                                           self.eastlst[ii])) * np.cos(theta)
            x0 = self.eastlst[ii] + d * np.sin(theta)
            y0 = self.northlst[ii] - d * np.cos(theta)
            surveylst[ii]['east'] = x0
            surveylst[ii]['north'] = y0

            # need to figure out a way to account for zone changes

            if lineori == 'ew':
                if surveylst[0]['east'] < surveylst[ii]['east']:
                    surveylst[ii]['offset'] = np.sqrt((surveylst[0]['east'] -
                                                       surveylst[ii]['east'])**2 +
                                                      (surveylst[0]['north'] -
                                                       surveylst[ii]['north'])**2)
                elif surveylst[0]['east'] > surveylst[ii]['east']:
                    surveylst[ii]['offset'] = -1 * np.sqrt((surveylst[0]['east'] -
                                                            surveylst[ii]['east'])**2 +
                                                           (surveylst[0]['north'] -
                                                            surveylst[ii]['north'])**2)
                else:
                    surveylst[ii]['offset'] = 0
            elif lineori == 'ns':
                if surveylst[0]['north'] < surveylst[ii]['north']:
                    surveylst[ii]['offset'] = np.sqrt((surveylst[0]['east'] -
                                                       surveylst[ii]['east'])**2 +
                                                      (surveylst[0]['north'] -
                                                       surveylst[ii]['north'])**2)
                elif surveylst[0]['north'] > surveylst[ii]['north']:
                    surveylst[ii]['offset'] = -1 * np.sqrt((surveylst[0]['east'] -
                                                            surveylst[ii]['east'])**2 +
                                                           (surveylst[0]['north'] -
                                                            surveylst[ii]['north'])**2)
                else:
                    surveylst[ii]['offset'] = 0

            if plotyn == 'y':
                lax.plot(x0, y0, 'v', color='k', ms=8, mew=3)
                lax.text(x0, y0 + 100, pstationlst[ii],
                         horizontalalignment='center',
                         verticalalignment='baseline',
                         fontdict={'size': 12, 'weight': 'bold'})

        # sort by ascending order of distance from first station
        surveylst = sorted(surveylst, key=itemgetter('offset'))

        # number of stations read
        nstat = len(surveylst)

        #--------------------------Match Frequencies---------------------------
        # a dictionary is created with the frequency as the key and the value is
        # the frequency number in the list. Each edi file is iterated over
        # extracting only the matched frequencies.  This makes it necessary to
        # have the same frequency content in each edifile.  If the frequencies
        # do not match then you can specify a tolerance to look around for
        # each frequency.

        # make a list to iterate over frequencies
        if type(freqstep) is list or type(freqstep) is not int:
            if type(freqstep[0]) is int:
                # find the median frequency list
                maxflen = max([len(ff) for ff in freqlst])
                farray = np.zeros((nstat, maxflen))
                for ii in range(nstat):
                    farray[ii, 0:len(freqlst[ii])] = freqlst[ii]

                mfreq = np.median(farray, axis=0)
                print(len(mfreq), len(freqstep))
                fdict = dict([('%.6g' % mfreq[ff], ii)
                              for ii, ff in enumerate(freqstep, 1)
                              if mfreq[ff] != 0])
            else:
                fdict = dict([('%.6g' % ff, ii)
                              for ii, ff in enumerate(freqstep, 1)])
        else:
            # find the median frequency list
            maxflen = max([len(ff) for ff in freqlst])
            farray = np.zeros((nstat, maxflen))
            for ii in range(nstat):
                farray[ii, 0:len(freqlst[ii])] = freqlst[ii]

            mfreq = np.median(farray, axis=0)

            # make a dictionary of values
            fdict = dict([('%.6g' % ff, ii) for ii, ff in
                          enumerate(mfreq[list(range(0, maxflen, freqstep))], 1)
                          if ff != 0])

        # print the frequencies to look for to make sure its what the user wants
        # make a list of keys that is sorted in descending order
        klst = [float(dd) for dd in list(fdict.keys())]
        klst.sort(reverse=True)
        klst = ['%.6g' % dd for dd in klst]

        print('Frequencies to look for are: (# freq(Hz) Period(s)) ')
        for key in klst:
            print(fdict[key], key, 1. / float(key))

        # make lists of parameters to write to file
        reslst = []
        offsetlst = []
        stationlstsort = []
        for kk in range(nstat):
            z = surveylst[kk]['z']
            zvar = surveylst[kk]['zvar']
            freq = surveylst[kk]['freq']
            offsetlst.append(surveylst[kk]['offset'])
            stationlstsort.append(surveylst[kk]['station'])
            tip = surveylst[kk]['tipper']
            tipvar = surveylst[kk]['tippervar']
            # loop over frequencies to pick out the ones desired
            for jj, ff in enumerate(freq):
                # jj is the index of edi file frequency list, this index
                # corresponds to the impedance tensor component index
                # ff is the frequency from the edi file frequency list
                try:
                    # nn is the frequency number out of extracted frequency
                    # list
                    nn = fdict['%.6g' % ff]

                    # calculate apparent resistivity
                    wt = .2 / (ff)
                    resxy = wt * abs(z[jj, 0, 1])**2
                    resyx = wt * abs(z[jj, 1, 0])**2

                    # calculate the phase putting the yx in the 1st quadrant
                    phasexy = np.arctan2(z[jj, 0, 1].imag,
                                         z[jj, 0, 1].real) * (180 / np.pi)
                    phaseyx = np.arctan2(z[jj, 1, 0].imag,
                                         z[jj, 1, 0].real) * (180 / np.pi) + 180
                    # put phases in correct quadrant if should be negative
                    if phaseyx > 180:
                        phaseyx = phaseyx - 360
                        print('Found Negative Phase', surveylst[kk]['station'], ff)

                    # calculate errors
                    #res_xy (TE)
                    if resxyerr == 'data':
                        dresxyerr = wt * \
                            (abs(z[jj, 0, 1]) + zvar[jj, 0, 1])**2 - resxy
                        lresxyerr = (dresxyerr / resxy) / np.log(10)

                    else:
                        lresxyerr = (resxyerr / 100.) / np.log(10)

                    # Res_yx(TM)
                    if resyxerr == 'data':
                        dresyxerr = wt * \
                            (abs(z[jj, 1, 0]) + zvar[jj, 1, 0])**2 - resyx
                        lresyxerr = (dresyxerr / resyx) / np.log(10)
                    else:
                        lresyxerr = (resyxerr / 100.) / np.log(10)

                    # phase_xy(TE)
                    if phasexyerr == 'data':
                        dphasexyerr = np.arcsin(zvar[jj, 0, 1] / abs(z[jj, 0, 1])) *\
                            (180 / np.pi)
                    else:
                        dphasexyerr = (phasexyerr / 100.) * 57 / 2.

                    #phase_yx (TM)
                    if phaseyxerr == 'data':
                        dphaseyxerr = np.arcsin(zvar[jj, 1, 0] / abs(z[jj, 1, 0])) *\
                            (180 / np.pi)
                    else:
                        dphaseyxerr = (phaseyxerr / 100.) * 57 / 2.

                    # calculate log10 of resistivity as prescribed by OCCAM
                    lresyx = np.log10(resyx)
                    lresxy = np.log10(resxy)

                    # if include the tipper
                    if tipper_err != None:
                        if tip[jj, 0, 0].real == 0.0 or tip[jj, 0, 1] == 0.0:
                            tipyn = 'n'
                        else:
                            # calculate the projection angle for real and
                            # imaginary
                            tipphir = np.arctan(tip[jj, 0, 0].real / tip[jj, 0, 1].real) -\
                                theta
                            tipphii = np.arctan(tip[jj, 0, 0].imag / tip[jj, 0, 1].imag) -\
                                theta

                            # project the tipper onto the profile line
                            projtipr = np.sqrt(tip[jj, 0, 0].real**2 + tip[jj, 0, 1].real**2) *\
                                np.cos(tipphir)
                            projtipi = np.sqrt(tip[jj, 0, 0].imag**2 + tip[jj, 0, 1].imag**2) *\
                                np.cos(tipphii)

                            # error of tipper is a decimal percentage
                            projtiperr = tipper_err / 100.

                            tipyn = 'y'

                    # make a list of lines to write to the data file
                    if mmode == 'both':
                        reslst.append(ss + str(kk + 1) + ss + str(nn) + ss + '1' + ss +
                                      string_fmt % lresxy + ss + string_fmt % lresxyerr + '\n')
                        reslst.append(ss + str(kk + 1) + ss + str(nn) + ss + '2' + ss +
                                      string_fmt % phasexy + ss + string_fmt % dphasexyerr + '\n')
                        reslst.append(ss + str(kk + 1) + ss + str(nn) + ss + '5' + ss +
                                      string_fmt % lresyx + ss + string_fmt % lresyxerr + '\n')
                        reslst.append(ss + str(kk + 1) + ss + str(nn) + ss + '6' + ss +
                                      string_fmt % phaseyx + ss + string_fmt % dphaseyxerr + '\n')
                        if tipper_err != None and tipyn == 'y':
                            reslst.append(ss + str(kk + 1) + ss + str(nn) + ss + '3' + ss +
                                          string_fmt % projtipr + ss + string_fmt % projtiperr + '\n')
                            reslst.append(ss + str(kk + 1) + ss + str(nn) + ss + '4' + ss +
                                          string_fmt % projtipi + ss + string_fmt % projtiperr + '\n')
                    elif mmode == 'TM':
                        reslst.append(ss + str(kk + 1) + ss + str(nn) + ss + '5' + ss +
                                      string_fmt % lresyx + ss + string_fmt % lresyxerr + '\n')
                        reslst.append(ss + str(kk + 1) + ss + str(nn) + ss + '6' + ss +
                                      string_fmt % phaseyx + ss + string_fmt % dphaseyxerr + '\n')
                        if tipper_err != None and tipyn == 'y':
                            reslst.append(ss + str(kk + 1) + ss + str(nn) + ss + '3' + ss +
                                          string_fmt % projtipr + ss + string_fmt % projtiperr + '\n')
                            reslst.append(ss + str(kk + 1) + ss + str(nn) + ss + '4' + ss +
                                          string_fmt % projtipi + ss + string_fmt % projtiperr + '\n')
                    elif mmode == 'TE':
                        reslst.append(ss + str(kk + 1) + ss + str(nn) + ss + '1' + ss +
                                      string_fmt % lresxy + ss + string_fmt % lresxyerr + '\n')
                        reslst.append(ss + str(kk + 1) + ss + str(nn) + ss + '2' + ss +
                                      string_fmt % phasexy + ss + string_fmt % dphasexyerr + '\n')
                        if tipper_err != None and tipyn == 'y':
                            reslst.append(ss + str(kk + 1) + ss + str(nn) + ss + '3' + ss +
                                          string_fmt % projtipr + ss + string_fmt % projtiperr + '\n')
                            reslst.append(ss + str(kk + 1) + ss + str(nn) + ss + '4' + ss +
                                          string_fmt % projtipi + ss + string_fmt % projtiperr + '\n')
                    else:
                        raise NameError('mmode' + mmode + ' not defined')
                except KeyError:
                    # search around the frequency given by ftol
                    try:
                        for key in list(fdict.keys()):
                            if ff > float(key) * (1 - ftol) and ff < float(key) * (1 + ftol):
                                nn = fdict[key]
                                wt = .2 / (ff)
                                resxy = wt * abs(z[jj, 0, 1])**2
                                resyx = wt * abs(z[jj, 1, 0])**2

                                # calculate the phase putting the yx in the 1st
                                # quadrant
                                phasexy = np.arctan2(z[jj, 0, 1].imag, z[jj, 0, 1].real) *\
                                    (180 / np.pi)
                                phaseyx = np.arctan2(z[jj, 1, 0].imag, z[jj, 1, 0].real) *\
                                    (180 / np.pi) + 180
                                # put phases in correct quadrant if should be
                                # negative
                                if phaseyx > 180:
                                    phaseyx = phaseyx - 360
                                    print('Found Negative Phase', surveylst[kk]['station'], ff)

                                # calculate errors
                                #res_xy (TE)
                                if resxyerr == 'data':
                                    dresxyerr = wt * \
                                        (abs(z[jj, 0, 1]) +
                                         zvar[jj, 0, 1])**2 - resxy
                                    lresxyerr = (
                                        dresxyerr / resxy) / np.log(10)

                                else:
                                    lresxyerr = (resxyerr / 100.) / np.log(10)

                                # Res_yx(TM)
                                if resyxerr == 'data':
                                    dresyxerr = wt * \
                                        (abs(z[jj, 1, 0]) +
                                         zvar[jj, 1, 0])**2 - resyx
                                    lresyxerr = (
                                        dresyxerr / resyx) / np.log(10)
                                else:
                                    lresyxerr = (resyxerr / 100.) / np.log(10)

                                # phase_xy(TE)
                                if phasexyerr == 'data':
                                    dphasexyerr = np.arcsin(zvar[jj, 0, 1] / abs(z[jj, 0, 1])) *\
                                        (180 / np.pi)
                                else:
                                    dphasexyerr = (phasexyerr / 100.) * 57 / 2.

                                #phase_yx (TM)
                                if phaseyxerr == 'data':
                                    dphaseyxerr = np.arcsin(zvar[jj, 1, 0] / abs(z[jj, 1, 0])) *\
                                        (180 / np.pi)
                                else:
                                    dphaseyxerr = (phaseyxerr / 100.) * 57 / 2.

                                # calculate log10 of resistivity as prescribed
                                # by OCCAM
                                lresyx = np.log10(resyx)
                                lresxy = np.log10(resxy)

                                # if include the tipper
                                if tipper_err != None:
                                    if tip[jj, 0].real == 0.0 or tip[jj, 1] == 0.0:
                                        tipyn = 'n'
                                    else:
                                        # calculate the projection angle for
                                        # real and imaginary
                                        tipphir = np.arctan(
                                            tip[jj, 0, 0].real / tip[jj, 0, 1].real) - theta
                                        tipphii = np.arctan(
                                            tip[jj, 0, 0].imag / tip[jj, 0, 1].imag) - theta

                                        # project the tipper onto the profile
                                        # line
                                        projtipr = np.sqrt(tip[jj, 0, 0].real**2 + tip[jj, 0, 1].real**2) *\
                                            np.cos(tipphir)
                                        projtipi = np.sqrt(tip[jj, 0, 0].imag**2 + tip[jj, 0, 1].imag**2) *\
                                            np.cos(tipphii)

                                        # error of tipper is a decimal
                                        # percentage
                                        projtiperr = tipper_err / 100.

                                        tipyn = 'y'

                                # make a list of lines to write to the data
                                # file
                                if mmode == 'both':
                                    reslst.append(ss + str(kk + 1) + ss + str(nn) + ss + '1' + ss +
                                                  string_fmt % lresxy + ss + string_fmt % lresxyerr + '\n')
                                    reslst.append(ss + str(kk + 1) + ss + str(nn) + ss + '2' + ss +
                                                  string_fmt % phasexy + ss + string_fmt % dphasexyerr + '\n')
                                    reslst.append(ss + str(kk + 1) + ss + str(nn) + ss + '5' + ss +
                                                  string_fmt % lresyx + ss + string_fmt % lresyxerr + '\n')
                                    reslst.append(ss + str(kk + 1) + ss + str(nn) + ss + '6' + ss +
                                                  string_fmt % phaseyx + ss + string_fmt % dphaseyxerr + '\n')
                                    if tipper_err != None and tipyn == 'y':
                                        reslst.append(ss + str(kk + 1) + ss + str(nn) + ss + '3' + ss +
                                                      string_fmt % projtipr + ss + string_fmt % projtiperr + '\n')
                                        reslst.append(ss + str(kk + 1) + ss + str(nn) + ss + '4' + ss +
                                                      string_fmt % projtipi + ss + string_fmt % projtiperr + '\n')
                                elif mmode == 'TM':
                                    reslst.append(ss + str(kk + 1) + ss + str(nn) + ss + '5' + ss +
                                                  string_fmt % lresyx + ss + string_fmt % lresyxerr + '\n')
                                    reslst.append(ss + str(kk + 1) + ss + str(nn) + ss + '6' + ss +
                                                  string_fmt % phaseyx + ss + string_fmt % dphaseyxerr + '\n')
                                    if tipper_err != None and tipyn == 'y':
                                        reslst.append(ss + str(kk + 1) + ss + str(nn) + ss + '3' + ss +
                                                      string_fmt % projtipr + ss + string_fmt % projtiperr + '\n')
                                        reslst.append(ss + str(kk + 1) + ss + str(nn) + ss + '4' + ss +
                                                      string_fmt % projtipi + ss + string_fmt % projtiperr + '\n')
                                elif mmode == 'TE':
                                    reslst.append(ss + str(kk + 1) + ss + str(nn) + ss + '1' + ss +
                                                  string_fmt % lresxy + ss + string_fmt % lresxyerr + '\n')
                                    reslst.append(ss + str(kk + 1) + ss + str(nn) + ss + '2' + ss +
                                                  string_fmt % phasexy + ss + string_fmt % dphasexyerr + '\n')
                                    if tipper_err != None and tipyn == 'y':
                                        reslst.append(ss + str(kk + 1) + ss + str(nn) + ss + '3' + ss +
                                                      string_fmt % projtipr + ss + string_fmt % projtiperr + '\n')
                                        reslst.append(ss + str(kk + 1) + ss + str(nn) + ss + '4' + ss +
                                                      string_fmt % projtipi + ss + string_fmt % projtiperr + '\n')
                                else:
                                    raise NameError(
                                        'mmode' + mmode + ' not defined')

                                break
                            else:
                                pass
            #                           print 'Did not find frequency {0} for station {1}'.format(ff,surveylst[kk]['station'])
                                # calculate resistivity

                    except KeyError:
                        pass

        #======================================================================
        #                             write dat file
        #======================================================================
        if savepath != None:
            if os.path.basename(savepath).find('.') > 0:
                self.datafn = savepath
            else:
                if not os.path.exists(savepath):
                    os.mkdir(savepath)
                self.datafn = os.path.join(savepath, 'Data.dat')
        else:
            self.datafn = os.path.join(edipath, 'Data.dat')

        if title == None:
            title = 'Occam Inversion'

        datfid = open(self.datafn, 'w')
        datfid.write('FORMAT:' + ' ' * 11 + 'OCCAM2MTDATA_1.0' + '\n')
        datfid.write('TITLE:' + ' ' * 12 + '{0:.4g}--'.format(theta * 180 / np.pi) + ' ' +
                     title + '\n')

        # write station sites
        datfid.write('SITES:' + ' ' * 12 + str(nstat) + '\n')
        for station in stationlstsort:
            datfid.write(ss + station + '\n')

        # write offsets
        datfid.write('OFFSETS (M):' + '\n')
        for offset in offsetlst:
            datfid.write(ss + string_fmt % offset + '\n')

        # write frequencies
        #writefreq=[freq[ff] for ff in range(0,len(freq),freqstep)]
        datfid.write('FREQUENCIES:' + ' ' * 8 + str(len(fdict)) + '\n')
        for fkey in klst:
            datfid.write(ss + string_fmt % float(fkey) + '\n')

        # write data block
        datfid.write('DATA BLOCKS:' + ' ' * 10 + str(len(reslst)) + '\n')
        datfid.write('SITE' + ss + 'FREQ' + ss + 'TYPE' +
                     ss + 'DATUM' + ss + 'ERROR' + '\n')
        for ll, datline in enumerate(reslst):
            if datline.find('#IND') >= 0:
                print('Found #IND on line ', ll)
                ndline = datline.replace('#IND', '00')
                print('Replaced with 00')
                datfid.write(ndline)
            else:
                datfid.write(datline)
        datfid.close()

        print('Wrote Occam2D data file to: ', self.datafn)

    def read2DdataFile(self):
        """
            read2DdataFile will read in data from a 2D occam data file.  
            Only supports the first 6 data types of occam2D

        Arguments:
        ----------

            **OccamPointPicker.datafn** : full path to data file

        Returns:
        --------
            **OccamPointPicker.rplst** : list of dictionaries for each station 
                                         with keywords:

                *'station'* : string
                              station name

                *'offset'* : float
                            relative offset

                *'resxy'* : np.array(nf,4)
                          TE resistivity and error as row 0 and 1 ressectively

                *'resyx'* : np.array(fn,4)
                          TM resistivity and error as row 0 and 1 respectively

                *'phasexy'* : np.array(nf,4)
                            TE phase and error as row 0 and 1 respectively

                *'phaseyx'* : np.array(nf,4)
                            Tm phase and error as row 0 and 1 respectively

                *'realtip'* : np.array(nf,4)
                            Real Tipper and error as row 0 and 1 respectively

                *'imagtip'* : np.array(nf,4)
                            Imaginary Tipper and error as row 0 and 1 
                            respectively

                Note: that the resistivity will be in log10 space.  Also, there
                are 2 extra rows in the data arrays, this is to put the 
                response from the inversion. 

        :Example: ::

            >>> import mtpy.modeling.occamtools as occam
            >>> ocd = occam.Occam2DData()
            >>> ocd.datafn = r"/home/Occam2D/Line1/Inv1/Data.dat"
            >>> ocd.read2DdataFile()

        """

        dfid = open(self.datafn, 'r')

        dlines = dfid.readlines()
        # get format of input data
        self.occamfmt = dlines[0].strip().split(':')[1].strip()

        # get title
        self.titlestr = dlines[1].strip().split(':')[1].strip()

        if self.titlestr.find('--') > 0:
            tstr = self.titlestr.split('--')
            self.theta_profile = float(tstr[0])
            self.title = tstr[1]
        else:
            self.title = self.titlestr
            self.theta_profile = 0
            print('Need to figure out angle of profile line')
        # get number of sits
        nsites = int(dlines[2].strip().split(':')[1].strip())

        # get station names
        self.stationlst = [dlines[ii].strip() for ii in range(3, nsites + 3)]

        # get offsets in meters
        offsets = [float(dlines[ii].strip())
                   for ii in range(4 + nsites, 4 + 2 * nsites)]

        # get number of frequencies
        nfreq = int(dlines[4 + 2 * nsites].strip().split(':')[1].strip())

        # get frequencies
        self.freq = np.array([float(dlines[ii].strip())
                              for ii in range(5 + 2 * nsites, 5 + 2 * nsites + nfreq)])

        # get periods
        self.period = 1. / self.freq

        #-----------get data-------------------
        # set zero array size the first row will be the data and second the
        # error
        asize = (4, nfreq)
        # make a list of dictionaries for each station.
        self.rplst = [{'station': station, 'offset': offsets[ii],
                       'resxy':np.zeros(asize),
                       'resyx':np.zeros(asize),
                       'phasexy':np.zeros(asize),
                       'phaseyx':np.zeros(asize),
                       'realtip':np.zeros(asize),
                       'imagtip':np.zeros(asize),
                       } for ii, station in enumerate(self.stationlst)]
        for line in dlines[7 + 2 * nsites + nfreq:]:
            ls = line.split()
            # station index
            ss = int(float(ls[0])) - 1
            # component key
            comp = str(int(float(ls[2])))
            # frequency index
            ff = int(float(ls[1])) - 1
            # print ls,ss,comp,ff
            # put into array
            # input data
            self.rplst[ss][occamdict[comp]][0, ff] = float(ls[3])
            # error
            self.rplst[ss][occamdict[comp]][1, ff] = float(ls[4])

    def rewrite2DdataFile(self, edipath=None, thetar=0, resxyerr='prev',
                          resyxerr='prev', phasexyerr='prev', phaseyxerr='prev',
                          tipper_err=None, mmode='both', flst=None,
                          removestation=None, savepath=None):
        """
        rewrite2DDataFile will rewrite an existing data file so you can 
        redefine some of the parameters, such as rotation angle, or errors for 
        the different components or only invert for one mode or add one or add
        tipper or remove tipper.

        Arguments:
        ----------

            **datafn** : string
                         full path to data file to rewrite

            **thetar** : float
                         rotation angle with positive clockwise

            **resxyerr** : float
                           error for TE mode resistivity (percent) or 'data' 
                           for data or 'prev' to take errors from data file.

            **resyxerr** : float
                           error for TM mode resistivity (percent) or 'data' 
                           for data or 'prev' to take errors from data file.

            **phasexyerr** : float
                             error for TE mode phase (percent) or 'data' 
                             for data or 'prev' to take errors from data file.

            **phaseyxerr** : float
                             error for TM mode phase (percent) or 'data' 
                             for data or 'prev' to take errors from data file.

            **tipper_err** : float 
                            error for tipper (percent) input only if you want
                            to invert for the tipper or 'data' for data errors
                            or prev to take errors from data file.

            **mmodes** : string can be:
                            * 'both' for both TE and TM
                            * 'TE' for TE
                            * 'TM' for TM

            **flst** : list or np.array
                       frequency list in Hz to rewrite, needs to be similar to
                       the datafile, cannot add frequencies

            **removestation** : list of stations to remove if desired

            **savepath** : string
                           full path to save the file to, if None then file 
                           is saved to os.path.dirname(datafn,'DataRW.dat')

        Returns:
        --------

            **Occam2DData.ndatafn** : string
                                      full path to new data file

        :Example: ::

            >>> import mtpy.modeling.occamtools as occam
            >>> ocd = occam.Occam2DData()
            >>> ocd.datafn = r"/home/Occam2D/Line1/Inv1/Data.dat"
            >>>
            >>> #rotate by 10 degrees east of North, remove station MT03 and 
            >>> #increase the resistivity error to 30 percent and put into a new 
            >>> #folder using savepath
            >>> svpath = r"/home/Occam2D/Line1/Inv2"
            >>> edipath = r"/home/EDIfiles"
            >>> ocd.rewrite2DdataFile(edipath=edipath,thetar=10,resxyerr=30,
            >>>                       resyxerr=30,removstation='MT03',
            >>>                       savepath=svpath)
            >>> Rewrote the data file to: /home/Occam2D/Line1/Inv2/DataRW.dat
        """

        ss = 3 * ' '
        string_fmt = '%2.6f'

        # load the data for the data file
        self.read2DdataFile()
        # copy the information into local lists in case you want to keep the
        # original data
        rplst = list(self.rplst)
        stationlst = list(self.stationlst)
        # make a dictionary of rplst for easier extraction of data
        rpdict = dict([(station, rplst[ii]) for ii, station in
                       enumerate(stationlst)])

        # remove stations from rplst and stationlst if desired
        if removestation != None:
            # if removestation is not a list make it one
            if type(removestation) is not list:
                removestation = [removestation]

            # remove station from station list
            for rstation in removestation:
                try:
                    stationlst.remove(rstation)
                except ValueError:
                    print('Did not find ' + rstation)

        # if flst is not the same as freq make freq=flst
        if flst != None:
            freq = flst
        else:
            freq = self.freq

        # if the rotation angle is not 0 than need to read the original data in
        if thetar != 0:
            if edipath == None:
                raise IOError('Need to input the edipath to original edifiles to' +
                              ' get rotations correct')

            # get list of edifiles already in data file
            edilst = [os.path.join(edipath, edi) for stat in stationlst
                      for edi in os.listdir(edipath) if edi[0:len(stat)] == stat]
            reslst = []
            for kk, edifn in enumerate(edilst, 1):
                imp1 = Z.Z(edifn)
                rp = imp1.getResPhase(thetar=thetar)
                imptip = imp1.getTipper()
                tip = imptip.tipper
                station = stationlst[kk - 1]
                fdict = dict([('{0:.6g}'.format(fr), ii) for ii, fr in
                              enumerate(imp1.frequency)])
                # loop over frequencies to pick out the ones desired
                for jj, ff in enumerate(freq, 1):
                    # jj is the index of edi file frequency list, this index corresponds
                    # to the impedance tensor component index
                    # ff is the frequency from the edi file frequency list
                    try:
                        # nn is the frequency number out of extracted frequency
                        # list
                        nn = fdict['%.6g' % ff]

                        # calculate resistivity
                        resxy = rp.resxy[nn]
                        resyx = rp.resyx[nn]

                        # calculate the phase putting the yx in the 1st
                        # quadrant
                        phasexy = rp.phasexy[nn]
                        phaseyx = rp.phaseyx[nn] + 180
                        # put phases in correct quadrant if should be negative
                        if phaseyx > 180:
                            phaseyx = phaseyx - 360
                            print('Found Negative Phase at', imp1.station, kk, ff)

                        # calculate errors
                        #res_xy (TE)
                        if resxyerr == 'data':
                            lresxyerr = (rp.resxyerr[nn] / resxy) / np.log(10)
                        # take errors from data file
                        elif resxyerr == 'prev':
                            lresxyerr = rpdict[station]['resxy'][1, jj - 1]
                        else:
                            lresxyerr = (resxyerr / 100.) / np.log(10)

                        # Res_yx(TM)
                        if resyxerr == 'data':
                            lresxyerr = rpdict[station]['resyx'][1, jj - 1]
                        # take errors from data file
                        elif resyxerr == 'prev':
                            lresyxerr = rpdict[station]['resyx'][1, jj - 1]
                        else:
                            lresyxerr = (resyxerr / 100.) / np.log(10)

                        # phase_xy(TE)
                        if phasexyerr == 'data':
                            dphasexyerr = rp.phasexyerr[nn]
                            # take errors from data file
                        elif phasexyerr == 'prev':
                            dphasexyerr = rpdict[station]['phasexy'][1, jj - 1]
                        else:
                            dphasexyerr = (phasexyerr / 100.) * 57 / 2.

                        #phase_yx (TM)
                        if phaseyxerr == 'data':
                            dphaseyxerr = rp.phaseyxerr[nn]
                        elif phaseyxerr == 'prev':
                            dphaseyxerr = rpdict[station]['phaseyx'][1, jj - 1]
                        else:
                            dphaseyxerr = (phaseyxerr / 100.) * 57 / 2.

                        # calculate log10 of resistivity as prescribed by OCCAM
                        lresyx = np.log10(resyx)
                        lresxy = np.log10(resxy)

                        # if include the tipper
                        if tipper_err != None:
                            if tip[nn, 0] == 0.0 or tip[nn, 1] == 0.0:
                                tipyn = 'n'
                            else:
                                # calculate the projection angle for real and
                                # imaginary
                                tipphir = np.arctan(tip[nn, 0].real / tip[nn, 1].real) -\
                                    self.theta_profile
                                tipphii = np.arctan(tip[nn, 0].imag / tip[nn, 1].imag) -\
                                    self.theta_profile

                                # project the tipper onto the profile line
                                projtipr = np.sqrt(tip[nn, 0].real**2 + tip[nn, 1].real**2) *\
                                    np.cos(tipphir)
                                projtipi = np.sqrt(tip[nn, 0].imag**2 + tip[nn, 1].imag**2) *\
                                    np.cos(tipphii)

                                # error of tipper is a decimal percentage
                                projtiperr = tipper_err / 100.

                                tipyn = 'y'

                        # make a list of lines to write to the data file
                        if mmode == 'both':
                            if rpdict[station]['resxy'][0, jj - 1] != 0.0:
                                reslst.append(ss + str(kk) + ss + str(jj) + ss + '1' + ss +
                                              string_fmt % lresxy + ss + string_fmt % lresxyerr + '\n')
                            if rpdict[station]['phasexy'][0, jj - 1] != 0.0:
                                reslst.append(ss + str(kk) + ss + str(jj) + ss + '2' + ss +
                                              string_fmt % phasexy + ss + string_fmt % dphasexyerr + '\n')
                            if rpdict[station]['resyx'][0, jj - 1] != 0.0:
                                reslst.append(ss + str(kk) + ss + str(jj) + ss + '5' + ss +
                                              string_fmt % lresyx + ss + string_fmt % lresyxerr + '\n')
                            if rpdict[station]['phaseyx'][0, jj - 1] != 0.0:
                                reslst.append(ss + str(kk) + ss + str(jj) + ss + '6' + ss +
                                              string_fmt % phaseyx + ss + string_fmt % dphaseyxerr + '\n')
                            if tipper_err != None and tipyn == 'y':
                                if rpdict[station]['realtip'][0, jj - 1] != 0.0:
                                    reslst.append(ss + str(kk) + ss + str(jj) + ss + '3' + ss +
                                                  string_fmt % projtipr + ss + string_fmt % projtiperr +
                                                  '\n')
                                if rpdict[station]['imagtip'][0, jj - 1] != 0.0:
                                    reslst.append(ss + str(kk) + ss + str(jj) + ss + '4' + ss +
                                                  string_fmt % projtipi + ss + string_fmt % projtiperr +
                                                  '\n')
                        elif mmode == 'TM':
                            if rpdict[station]['resyx'][0, jj - 1] != 0.0:
                                reslst.append(ss + str(kk) + ss + str(jj) + ss + '5' + ss +
                                              string_fmt % lresyx + ss + string_fmt % lresyxerr + '\n')
                            if rpdict[station]['phaseyx'][0, jj - 1] != 0.0:
                                reslst.append(ss + str(kk) + ss + str(jj) + ss + '6' + ss +
                                              string_fmt % phaseyx + ss + string_fmt % dphaseyxerr + '\n')
                            if tipper_err != None and tipyn == 'y':
                                if rpdict[station]['realtip'][0, jj - 1] != 0.0:
                                    reslst.append(ss + str(kk) + ss + str(jj) + ss + '3' + ss +
                                                  string_fmt % projtipr + ss + string_fmt % projtiperr +
                                                  '\n')
                                if rpdict[station]['imagtip'][0, jj - 1] != 0.0:
                                    reslst.append(ss + str(kk) + ss + str(jj) + ss + '4' + ss +
                                                  string_fmt % projtipi + ss + string_fmt % projtiperr +
                                                  '\n')
                        elif mmode == 'TE':
                            if rpdict[station]['resxy'][0, jj - 1] != 0.0:
                                reslst.append(ss + str(kk) + ss + str(jj) + ss + '1' + ss +
                                              string_fmt % lresxy + ss + string_fmt % lresxyerr + '\n')
                            if rpdict[station]['phasexy'][0, jj - 1] != 0.0:
                                reslst.append(ss + str(kk) + ss + str(jj) + ss + '2' + ss +
                                              string_fmt % phasexy + ss + string_fmt % dphasexyerr + '\n')
                            if tipper_err != None and tipyn == 'y':
                                if rpdict[station]['realtip'][0, jj - 1] != 0.0:
                                    reslst.append(ss + str(kk) + ss + str(jj) + ss + '3' + ss +
                                                  string_fmt % projtipr + ss + string_fmt % projtiperr +
                                                  '\n')
                                if rpdict[station]['imagtip'][0, jj - 1] != 0.0:
                                    reslst.append(ss + str(kk) + ss + str(jj) + ss + '4' + ss +
                                                  string_fmt % projtipi + ss + string_fmt % projtiperr +
                                                  '\n')
                        else:
                            raise NameError('mmode' + mmode + ' not defined')
                    except KeyError:
                        pass

        # If no rotation is desired but error bars are than...
        else:
            reslst = []
            for kk, station in enumerate(stationlst, 1):
                srp = rpdict[station]
                nr = srp['resxy'].shape[1]
                # calculate errors and rewrite
                #res_xy (TE)
                if resxyerr != None:
                    if resxyerr == 'prev':
                        lresxyerr = rpdict[station]['resxy'][1, :]
                    else:
                        lresxyerr = np.repeat(
                            (resxyerr / 100.) / np.log(10), nr)
                    srp['resxy'][1, :] = lresxyerr

                # Res_yx(TM)
                if resyxerr != None:
                    if resyxerr == 'prev':
                        lresyxerr = rpdict[station]['resyx'][1, :]
                    else:
                        lresyxerr = np.repeat(
                            (resyxerr / 100.) / np.log(10), nr)
                    srp['resyx'][1, :] = lresyxerr

                # phase_xy(TE)
                if phasexyerr != None:
                    if phasexyerr == 'prev':
                        dphasexyerr = rpdict[station]['phasexy'][1, :]
                    else:
                        dphasexyerr = np.repeat(
                            (phasexyerr / 100.) * 57 / 2., nr)
                    srp['phasexy'][1, :] = dphasexyerr

                #phase_yx (TM)
                if phaseyxerr != None:
                    if phaseyxerr == 'prev':
                        dphaseyxerr = rpdict[station]['phaseyx'][1, :]
                    else:
                        dphaseyxerr = np.repeat(
                            (phaseyxerr / 100.) * 57 / 2., nr)
                    srp['phaseyx'][1, :] = dphaseyxerr

                if tipper_err != None:
                    # error of tipper is a decimal percentage
                    projtiperr = tipper_err / 100.
                    srp['realtip'][1, :] = np.repeat(projtiperr, nr)
                    srp['imagtip'][1, :] = np.repeat(projtiperr, nr)

                for jj, ff in enumerate(freq, 1):
                    # make a list of lines to write to the data file
                    if mmode == 'both':
                        if srp['resxy'][0, jj - 1] != 0.0:
                            reslst.append(ss + str(kk) + ss + str(jj) + ss + '1' + ss +
                                          string_fmt % srp['resxy'][0, jj - 1] + ss +
                                          string_fmt % srp['resxy'][1, jj - 1] + '\n')
                        if srp['phasexy'][0, jj - 1] != 0.0:
                            reslst.append(ss + str(kk) + ss + str(jj) + ss + '2' + ss +
                                          string_fmt % srp['phasexy'][0, jj - 1] + ss +
                                          string_fmt % srp['phasexy'][1, jj - 1] + '\n')
                        if srp['resyx'][0, jj - 1] != 0.0:
                            reslst.append(ss + str(kk) + ss + str(jj) + ss + '5' + ss +
                                          string_fmt % srp['resyx'][0, jj - 1] + ss +
                                          string_fmt % srp['resyx'][1, jj - 1] + '\n')
                        if srp['phaseyx'][0, jj - 1] != 0.0:
                            reslst.append(ss + str(kk) + ss + str(jj) + ss + '6' + ss +
                                          string_fmt % srp['phaseyx'][0, jj - 1] + ss +
                                          string_fmt % srp['phaseyx'][1, jj - 1] + '\n')
                        if tipper_err != None:
                            if srp['realtip'][0, jj - 1] != 0.0:
                                reslst.append(ss + str(kk) + ss + str(jj) + ss + '3' + ss +
                                              string_fmt % srp['realtip'][0, jj - 1] + ss +
                                              string_fmt % srp['realtip'][1, jj - 1] + '\n')
                            if srp['imagtip'][0, jj - 1] != 0.0:
                                reslst.append(ss + str(kk) + ss + str(jj) + ss + '4' + ss +
                                              string_fmt % srp['imagtip'][0, jj - 1] + ss +
                                              string_fmt % srp['imagtip'][1, jj - 1] + '\n')
                    elif mmode == 'TM':
                        if srp['resyx'][0, jj - 1] != 0.0:
                            reslst.append(ss + str(kk) + ss + str(jj) + ss + '5' + ss +
                                          string_fmt % srp['resyx'][0, jj - 1] + ss +
                                          string_fmt % srp['resyx'][1, jj - 1] + '\n')
                        if srp['phaseyx'][0, jj - 1] != 0.0:
                            reslst.append(ss + str(kk) + ss + str(jj) + ss + '6' + ss +
                                          string_fmt % srp['phaseyx'][0, jj - 1] + ss +
                                          string_fmt % srp['phaseyx'][1, jj - 1] + '\n')
                        if tipper_err != None:
                            if srp['realtip'][0, jj - 1] != 0.0:
                                reslst.append(ss + str(kk) + ss + str(jj) + ss + '3' + ss +
                                              string_fmt % srp['realtip'][0, jj - 1] + ss +
                                              string_fmt % srp['realtip'][1, jj - 1] + '\n')
                            if srp['imagtip'][0, jj - 1] != 0.0:
                                reslst.append(ss + str(kk) + ss + str(jj) + ss + '4' + ss +
                                              string_fmt % srp['imagtip'][0, jj - 1] + ss +
                                              string_fmt % srp['imagtip'][1, jj - 1] + '\n')
                    elif mmode == 'TE':
                        if srp['resxy'][0, jj - 1] != 0.0:
                            reslst.append(ss + str(kk) + ss + str(jj) + ss + '1' + ss +
                                          string_fmt % srp['resxy'][0, jj - 1] + ss +
                                          string_fmt % srp['resxy'][1, jj - 1] + '\n')
                        if srp['phasexy'][0, jj - 1] != 0.0:
                            reslst.append(ss + str(kk) + ss + str(jj) + ss + '2' + ss +
                                          string_fmt % srp['phasexy'][0, jj - 1] + ss +
                                          string_fmt % srp['phasexy'][1, jj - 1] + '\n')
                        if tipper_err != None:
                            if srp['realtip'][0, jj - 1] != 0.0:
                                reslst.append(ss + str(kk) + ss + str(jj) + ss + '3' + ss +
                                              string_fmt % srp['realtip'][0, jj - 1] + ss +
                                              string_fmt % srp['realtip'][1, jj - 1] + '\n')
                            if srp['imagtip'][0, jj - 1] != 0.0:
                                reslst.append(ss + str(kk) + ss + str(jj) + ss + '4' + ss +
                                              string_fmt % srp['imagtip'][0, jj - 1] + ss +
                                              string_fmt % srp['imagtip'][1, jj - 1] + '\n')

        #======================================================================
        #                             write dat file
        #======================================================================

        # make the file name of the data file
        if self.datafn.find('RW') > 0:
            if savepath == None:
                self.ndatafn = self.datafn
            elif os.path.isdir(savepath) == True:
                self.ndatafn = os.path.join(savepath, 'DataRW.dat')
            elif os.path.isfile(savepath) == True:
                self.ndatafn = savepath
            elif savepath.find('.dat') > 0:
                self.ndatafn = savepath
        else:
            if savepath == None:
                self.ndatafn = self.datafn[:-4] + 'RW.dat'
            elif os.path.isdir(savepath) == True:
                self.ndatafn = os.path.join(savepath, 'DataRW.dat')
            elif os.path.isfile(savepath) == True:
                self.ndatafn = savepath
            elif savepath.find('.dat') > 0:
                self.ndatafn = savepath

        nstat = len(stationlst)

        if self.titlestr == None:
            self.titlestr = 'Occam Inversion'

        datfid = open(self.ndatafn, 'w')
        datfid.write('FORMAT:' + ' ' * 11 + 'OCCAM2MTDATA_1.0' + '\n')
        datfid.write('TITLE:' + ' ' * 12 + self.titlestr + '\n')

        # write station sites
        datfid.write('SITES:' + ' ' * 12 + str(nstat) + '\n')
        for station in stationlst:
            datfid.write(ss + station + '\n')

        # write offsets
        datfid.write('OFFSETS (M):' + '\n')
        projangle = (thetar - self.theta_profile) * np.pi / 180.
        for station in stationlst:
            # need to project the stations on to the strike direction
            datfid.write(ss + string_fmt % (rpdict[station]['offset'] * np.cos(projangle)) +
                         '\n')

        # write frequencies
        #writefreq=[freq[ff] for ff in range(0,len(freq),freqstep)]
        datfid.write('FREQUENCIES:' + ' ' * 8 + str(len(freq)) + '\n')
        for ff in self.freq:
            datfid.write(ss + string_fmt % ff + '\n')

        # write data block
        datfid.write('DATA BLOCKS:' + ' ' * 10 + str(len(reslst)) + '\n')
        datfid.write('SITE' + ss + 'FREQ' + ss + 'TYPE' +
                     ss + 'DATUM' + ss + 'ERROR' + '\n')
        for ll, datline in enumerate(reslst):
            if datline.find('#IND') >= 0:
                print('Found #IND on line ', ll)
                ndline = datline.replace('#IND', '00')
                print('Replaced with 00')
                datfid.write(ndline)
            else:
                datfid.write(datline)
        datfid.close()

        print('Rewrote the data file to: ', self.ndatafn)

    def makeModelFiles(self, niter=20, targetrms=1.0, nlayers=100, nlperdec=30,
                       z1layer=50, bwidth=200, trigger=.75, savepath=None, rhostart=100,
                       occampath=r"c:\Peacock\PHD\OCCAM\MakeFiles"):
        """
        makeModel will make an the input files for occam using Steve Constable's
        MakeModel2DMT.f code.

        Inputs:
            datafn = full path to data file
            niter = maximum number of iterations
            targetrms = target root mean square error
            nlayers = total number of layers in mesh
            nlperdec = number of layers per decade
            z1layer = thickness of the first layer in meters
            bwidth = maximum block width for regularization grid in meters
            trigger = triger point to amalgamate blocks
            savepath = path to save files to
            rhostart = starting resistivity for homogeneous half space in ohm-m
            occampath = path to MakeModel2DMT.exe

        Outputs:
            meshfn = mesh file for finite element grid saved ats MESH
            inmodelfn = input model, starting model with rhostart as starting value
                        saved as INMODEL
            startupfn = start up filepath, saved as startup
        """
        # get the base name of data file
        dfnb = os.path.basename(self.datafn)

        # put data file into the same directory as MakeModel2DMT
        if os.path.dirname(self.datafn) != occampath:
            shutil.copy(self.datafn, os.path.join(occampath, dfnb))

        # write input file for MakeModel2DMT
        mmfid = open(os.path.join(occampath, 'inputMakeModel.txt'), 'w')
        mmfid.write(dfnb + '\n')
        mmfid.write(str(niter) + '\n')
        mmfid.write(str(targetrms) + '\n')
        mmfid.write(str(nlayers) + '\n')
        mmfid.write(str(nlperdec) + '\n')
        mmfid.write(str(z1layer) + '\n')
        mmfid.write(str(bwidth) + '\n')
        mmfid.write(str(trigger) + '\n')
        mmfid.write('\n')
        mmfid.close()

        # get current working directory
        cdir = os.getcwd()

        # change directory path to occam path
        os.chdir(occampath)

        #---call MakeModel2DMT---
        subprocess.os.system("MakeModel2DMT < inputMakeModel.txt")

        # change back to original working directory
        os.chdir(cdir)

        if savepath == None:
            savepath = os.path.dirname(self.datafn)

        if not os.path.exists(savepath):
            os.mkdir(savepath)

        meshfn = os.path.join(savepath, 'MESH')
        inmodelfn = os.path.join(savepath, 'INMODEL')
        startupfn = os.path.join(savepath, 'startup')

        # copy ouput files to savepath
        shutil.copy(os.path.join(occampath, 'MESH'), meshfn)
        shutil.copy(os.path.join(occampath, 'INMODEL'), inmodelfn)
        shutil.copy(os.path.join(occampath, 'startup'), startupfn)
        shutil.copy(os.path.join(occampath, 'inputMakeModel.txt'),
                    os.path.join(savepath, 'inputMakeModel.txt'))
        if not os.path.exists(os.path.join(savepath, dfnb)):
            shutil.copy(self.datafn, os.path.join(savepath, dfnb))
        if os.path.getctime(os.path.join(savepath, dfnb)) <\
                os.path.getctime(self.datafn):
            shutil.copy(self.datafn, os.path.join(savepath, dfnb))


#        #rewrite mesh so it contains the right number of columns and rows
#        rewriteMesh(meshfn)

        # write startup file to have the starting desired starting rho value
        ifid = open(startupfn, 'r')
        ilines = ifid.readlines()
        ifid.close()

        if rhostart != 100:
            # make startup model a homogeneous half space of rhostart
            rhostart = np.log10(rhostart)
            ifid = open(startupfn, 'w')
            for line in ilines:
                if line.find('2.000000') >= 0:
                    line = line.replace('2.000000', '%.6f' % rhostart)
                ifid.write(line)
        ifid.close()

        print('Be sure to check the INMODEL file for clumped numbers near the bottom.')
        print('Also, check the MESH and startup files to make sure they are correct.')

        self.meshfn = meshfn
        self.inmodelfn = inmodelfn
        self.startupfn = startupfn

    def plotMaskPoints(self, plottype=None, reserrinc=.20, phaseerrinc=.05,
                       marker='h', colormode='color', dpi=300, ms=2,
                       reslimits=None, phaselimits=(-5, 95)):
        """
        An interactive plotting tool to mask points an add errorbars

        Arguments:
        ----------
            **plottype** : string
                           describes the way the responses are plotted. Can be:
                               *list of stations to plot ['mt01','mt02']
                               *one station 'mt01'
                               *None plots all stations *Default*

            **reserrinc** : float
                            amount to increase the error bars. Input as a 
                            decimal percentage.  0.3 for 30 percent
                            *Default* is 0.2 (20 percent)

            **phaseerrinc** : float
                              amount to increase the error bars. Input as a 
                              decimal percentage.  0.3 for 30 percent
                              *Default* is 0.05 (5 percent)

            **marker** : string
                         marker that the masked points will be
                         *Default* is 'h' for hexagon

            **colormode** : string
                            defines the color mode to plot the responses in
                            *'color'* for color plots
                            *'bw'* for black and white plots
                            *Default* is 'color'

            **dpi** : int
                      dot-per-inch resolution of the plots.
                      *Default* is 300

            **ms** : float
                     size of the marker in the response plots
                     *Default* is 2

            **reslimits**: tuple (min,max)
                           min and max limits of the resistivity values in 
                           linear scale
                           *Default* is None

            **phaselimits**: tuple (min,max)
                            min and max phase limits 
                            *Default* is (-5,90)

        Returns:
        ---------
            data type **OccamPointPicker**  


        :Example: ::

            >>> import mtpy.modeling.occamtools as occam
            >>> ocd = occam.Occam2DData()
            >>> ocd.datafn = r"/home/Occam2D/Line1/Inv1/Data.dat"
            >>> ocd.plotMaskPoints()                   

        """

        if colormode == 'color':
            # color for data
            cted = (0, 0, 1)
            ctmd = (1, 0, 0)
            mted = 's'
            mtmd = 'o'

        elif colormode == 'bw':
            # color for data
            cted = (0, 0, 0)
            ctmd = (0, 0, 0)
            mted = 's'
            mtmd = 'o'

        # read in data file
        self.read2DdataFile()
        rplst = list(self.rplst)

        # get periods
        period = self.period

        # define some empty lists to put things into
        pstationlst = []
        axlst = []
        linelst = []
        errlst = []

        # get the stations to plot
        # if none plot all of them
        if plottype == None:
            pstationlst = list(range(len(self.stationlst)))

        # otherwise pick out the stations to plot along with their index number
        elif type(plottype) is not list:
            plottype = [plottype]
            for ii, station in enumerate(self.stationlst):
                for pstation in plottype:
                    if station.find(pstation) >= 0:
                        pstationlst.append(ii)

        # set the subplot grid
        gs = gridspec.GridSpec(6, 2, wspace=.1, left=.1, top=.93, bottom=.07)
        for jj, ii in enumerate(pstationlst):
            fig = plt.figure(ii + 1, dpi=dpi)
            plt.clf()

            # make subplots
            axrte = fig.add_subplot(gs[:4, 0])
            axrtm = fig.add_subplot(gs[:4, 1])
            axpte = fig.add_subplot(gs[-2:, 0], sharex=axrte)
            axptm = fig.add_subplot(gs[-2:, 1], sharex=axrtm)

            # plot resistivity TE Mode
            # cut out missing data points first
            rxy = np.where(rplst[ii]['resxy'][0] != 0)[0]
            rte = axrte.errorbar(period[rxy], 10**rplst[ii]['resxy'][0][rxy],
                                 ls=':', marker=mted, ms=ms, mfc=cted, mec=cted,
                                 color=cted,
                                 yerr=np.log(10) * rplst[ii]['resxy'][1][rxy] *
                                 10**rplst[ii]['resxy'][0][rxy],
                                 ecolor=cted, picker=2)

            # plot Phase TE Mode
            # cut out missing data points first
            pxy = [np.where(rplst[ii]['phasexy'][0] != 0)[0]]
            pte = axpte.errorbar(period[pxy], rplst[ii]['phasexy'][0][pxy],
                                 ls=':', marker=mted, ms=ms, mfc=cted, mec=cted,
                                 color=cted, yerr=rplst[ii]['phasexy'][1][pxy],
                                 ecolor=cted, picker=1)

            # plot resistivity TM Mode
            # cut out missing data points first
            ryx = np.where(rplst[ii]['resyx'][0] != 0)[0]
            rtm = axrtm.errorbar(period[ryx], 10**rplst[ii]['resyx'][0][ryx],
                                 ls=':', marker=mtmd, ms=ms, mfc=ctmd, mec=ctmd,
                                 color=ctmd,
                                 yerr=np.log(10) * rplst[ii]['resyx'][1][ryx] *
                                 10**rplst[ii]['resyx'][0][ryx],
                                 ecolor=ctmd, picker=2)
            # plot Phase TM Mode
            # cut out missing data points first
            pyx = [np.where(rplst[ii]['phaseyx'][0] != 0)[0]]
            ptm = axptm.errorbar(period[pyx], rplst[ii]['phaseyx'][0][pyx],
                                 ls=':', marker=mtmd, ms=ms, mfc=ctmd, mec=ctmd,
                                 color=ctmd, yerr=rplst[ii]['phaseyx'][1][pyx],
                                 ecolor=ctmd, picker=1)

            # make the axis presentable
            # set the apparent resistivity scales to log and x-axis to log
            axplst = [axrte, axrtm, axpte, axptm]
            llst = [rte[0], rtm[0], pte[0], ptm[0]]
            elst = [[rte[1][0], rte[1][1], rte[2][0]],
                    [rtm[1][0], rtm[1][1], rtm[2][0]],
                    [pte[1][0], pte[1][1], pte[2][0]],
                    [ptm[1][0], ptm[1][1], ptm[2][0]]]

            axlst.append(axplst)
            linelst.append(llst)
            errlst.append(elst)

            # set the axes properties for each subplot
            for nn, xx in enumerate(axplst):
                # set xscale to logarithmic in period
                xx.set_xscale('log', nonposx='clip')

                # if apparent resistivity
                if nn == 0 or nn == 1:
                    # set x-ticklabels to invisible
                    plt.setp(xx.xaxis.get_ticklabels(), visible=False)

                    # set apparent resistivity scale to logarithmic
                    try:
                        xx.set_yscale('log', nonposy='clip')
                    except ValueError:
                        pass

                    # if there are resistivity limits set those
                    if reslimits != None:
                        xx.set_ylim(reslimits)

                # Set the title of the TE plot
                if nn == 0:
                    xx.set_title(self.stationlst[ii] + ' Obs$_{xy}$ (TE-Mode)',
                                 fontdict={'size': 9, 'weight': 'bold'})
                    xx.yaxis.set_label_coords(-.075, .5)
                    xx.set_ylabel('App. Res. ($\Omega \cdot m$)',
                                  fontdict={'size': 9, 'weight': 'bold'})
                # set the title of the TM plot
                if nn == 1:
                    xx.set_title(self.stationlst[ii] + ' Obs$_{yx}$ (TM-Mode)',
                                 fontdict={'size': 9, 'weight': 'bold'})

                # set the phase axes properties
                if nn == 2 or nn == 3:
                    # set the phase limits
                    xx.set_ylim(phaselimits)

                    # set label coordinates
                    xx.yaxis.set_label_coords(-.075, .5)

                    # give the y-axis label to the bottom left plot
                    if nn == 2:
                        xx.set_ylabel('Phase (deg)',
                                      fontdict={'size': 9, 'weight': 'bold'})
                    # set the x-axis label
                    xx.set_xlabel('Period (s)',
                                  fontdict={'size': 9, 'weight': 'bold'})

                    # set tick marks of the y-axis
                    xx.yaxis.set_major_locator(MultipleLocator(10))
                    xx.yaxis.set_minor_locator(MultipleLocator(2))

                xx.grid(True, alpha=.4, which='both')

        # make points an attribute of self which is a data type
        # OccamPointPicker
        self.points = OccamPointPicker(axlst, linelst, errlst, reserrinc=reserrinc,
                                       phaseerrinc=phaseerrinc, marker=marker)

        # be sure to show the plot
        plt.show()

    def maskPoints(self):
        """
        maskPoints will take in points found from plotMaskPoints and rewrite 
        the data file to nameRW.dat.  **Be sure to run plotMaskPoints first**

        Arguments:
        ---------
            None

        Returns:
        ---------

            **OccamPointPicker.ndatafn** : full path to rewritten data file

        :Example: ::

            >>> import mtpy.modeling.occamtools as occam
            >>> ocd = occam.Occam2DData()
            >>> ocd.datafn = r"/home/Occam2D/Line1/Inv1/Data.dat"
            >>> ocd.plotMaskPoints()
            >>>
            >>> #after all the points are masked rewrite the data file
            >>> ocd.maskPoints()
            >>>
            >>> Rewrote Occam2D data file to: /home/Occam2D/Line1/Inv1/DataRW.dat 

        """

        self.read2DdataFile()
        rplst = list(self.rplst)
        # rewrite the data file
        # make a reverse dictionary for locating the masked points in the data
        # file
        rploc = dict([('{0}'.format(self.points.fndict[key]), int(key) - 1)
                      for key in list(self.points.fndict.keys())])

        # make a period dictionary to locate points changed
        frpdict = dict([('{0:.5g}'.format(fr), ff)
                        for ff, fr in enumerate(1. / self.freq)])

        # loop over the data list
        for dd, dat in enumerate(self.points.data):
            derror = self.points.error[dd]
            # loop over the 4 main entrie
            for ss, skey in enumerate(['resxy', 'resyx', 'phasexy', 'phaseyx']):
                # rewrite any coinciding points
                for frpkey in list(frpdict.keys()):
                    try:
                        ff = frpdict[frpkey]
                        floc = self.points.fdict[dd][ss][frpkey]

                        # CHANGE APPARENT RESISTIVITY
                        if ss == 0 or ss == 1:
                            # change the apparent resistivity value
                            if rplst[rploc[str(dd)]][skey][0][ff] !=\
                                    np.log10(dat[ss][floc]):
                                if dat[ss][floc] == 0:
                                    rplst[rploc[str(dd)]][skey][0][ff] = 0.0
                                else:
                                    rplst[rploc[str(dd)]][skey][0][ff] =\
                                        np.log10(dat[ss][floc])

                            # change the apparent resistivity error value
                            if dat[ss][floc] == 0.0:
                                rerr = 0.0
                            else:
                                rerr = derror[ss][floc] / \
                                    dat[ss][floc] / np.log(10)
                            if rplst[rploc[str(dd)]][skey][1][ff] != rerr:
                                rplst[rploc[str(dd)]][skey][1][ff] = rerr

                        # DHANGE PHASE
                        elif ss == 2 or ss == 3:
                            # change the phase value
                            if rplst[rploc[str(dd)]][skey][0][ff] !=\
                                    dat[ss][floc]:
                                if dat[ss][floc] == 0:
                                    rplst[rploc[str(dd)]][skey][0][ff] = 0.0
                                else:
                                    rplst[rploc[str(dd)]][skey][0][ff] =\
                                        dat[ss][floc]

                            # change the apparent resistivity error value
                            if dat[ss][floc] == 0.0:
                                rerr = 0.0
                            else:
                                rerr = derror[ss][floc]
                            if rplst[rploc[str(dd)]][skey][1][ff] != rerr:
                                rplst[rploc[str(dd)]][skey][1][ff] = rerr
                    except KeyError:
                        pass

        # rewrite the data file
        ss = 3 * ' '
        string_fmt = '%2.6f'
        reslst = []

        # make a dictionary of rplst for easier extraction of data
        rpdict = dict([(station, rplst[ii])
                       for ii, station in enumerate(self.stationlst)])

        # loop over stations in the data file
        for kk, station in enumerate(self.stationlst, 1):
            srp = rpdict[station]

            # loop over frequencies
            for jj, ff in enumerate(self.freq, 1):
                # make a list of lines to write to the data file
                if srp['resxy'][0, jj - 1] != 0.0:
                    reslst.append(ss + str(kk) + ss + str(jj) + ss + '1' + ss +
                                  string_fmt % srp['resxy'][0, jj - 1] + ss +
                                  string_fmt % srp['resxy'][1, jj - 1] + '\n')
                if srp['phasexy'][0, jj - 1] != 0.0:
                    reslst.append(ss + str(kk) + ss + str(jj) + ss + '2' + ss +
                                  string_fmt % srp['phasexy'][0, jj - 1] + ss +
                                  string_fmt % srp['phasexy'][1, jj - 1] + '\n')
                if srp['resyx'][0, jj - 1] != 0.0:
                    reslst.append(ss + str(kk) + ss + str(jj) + ss + '5' + ss +
                                  string_fmt % srp['resyx'][0, jj - 1] + ss +
                                  string_fmt % srp['resyx'][1, jj - 1] + '\n')
                if srp['phaseyx'][0, jj - 1] != 0.0:
                    reslst.append(ss + str(kk) + ss + str(jj) + ss + '6' + ss +
                                  string_fmt % srp['phaseyx'][0, jj - 1] + ss +
                                  string_fmt % srp['phaseyx'][1, jj - 1] + '\n')
                if srp['realtip'][0, jj - 1] != 0.0:
                    reslst.append(ss + str(kk) + ss + str(jj) + ss + '3' + ss +
                                  string_fmt % srp['realtip'][0, jj - 1] + ss +
                                  string_fmt % srp['realtip'][1, jj - 1] + '\n')
                if srp['imagtip'][0, jj - 1] != 0.0:
                    reslst.append(ss + str(kk) + ss + str(jj) + ss + '4' + ss +
                                  string_fmt % srp['imagtip'][0, jj - 1] + ss +
                                  string_fmt % srp['imagtip'][1, jj - 1] + '\n')

        #======================================================================
        #                             write dat file
        #======================================================================
        # make the file name of the data file
        if self.datafn.find('RW') > 0:
            self.ndatafn = self.datafn
        else:
            self.ndatafn = self.datafn[:-4] + 'RW.dat'

        # get number of stations
        nstat = len(self.stationlst)

        # set title string
        if self.titlestr == None:
            self.titlestr = 'Occam Inversion'

        datfid = open(self.ndatafn, 'w')
        datfid.write('FORMAT:' + ' ' * 11 + 'OCCAM2MTDATA_1.0' + '\n')
        datfid.write('TITLE:' + ' ' * 12 + self.titlestr + '\n')

        # write station sites
        datfid.write('SITES:' + ' ' * 12 + str(nstat) + '\n')
        for station in self.stationlst:
            datfid.write(ss + station + '\n')

        # write offsets
        datfid.write('OFFSETS (M):' + '\n')
        for station in self.stationlst:
            datfid.write(ss + string_fmt % rpdict[station]['offset'] + '\n')

        # write frequencies
        #writefreq=[freq[ff] for ff in range(0,len(freq),freqstep)]
        datfid.write('FREQUENCIES:' + ' ' * 8 + str(len(self.freq)) + '\n')
        for ff in self.freq:
            datfid.write(ss + string_fmt % ff + '\n')

        # write data block
        datfid.write('DATA BLOCKS:' + ' ' * 10 + str(len(reslst)) + '\n')
        datfid.write('SITE' + ss + 'FREQ' + ss + 'TYPE' +
                     ss + 'DATUM' + ss + 'ERROR' + '\n')
        for ll, datline in enumerate(reslst):
            if datline.find('#IND') >= 0:
                print('Found #IND on line ', ll)
                ndline = datline.replace('#IND', '00')
                print('Replaced with 00')
                datfid.write(ndline)
            else:
                datfid.write(datline)
        datfid.close()

        print('Wrote Occam2D data file to: ', self.ndatafn)

    def read2DRespFile(self, respfn):
        """
        read2DRespFile will read in a response file and combine the data with info 
        from the data file.

        Arguments:
        ----------
            **respfn** : full path to the response file

            **datafn** : full path to data file

        Returns:
        --------
            for each data array, the rows are ordered as:
                - 0 -> input data
                - 1 -> input error
                - 2 -> model output
                - 3 -> relative error (data-model)/(input error)

            **rplst** : list of dictionaries for each station with keywords:

                *station* : string
                            station name

                *offset* : float
                            relative offset

                *resxy* : np.array(nf,4)
                          TE resistivity and error as row 0 and 1 ressectively

                *resyx* : np.array(fn,4)
                          TM resistivity and error as row 0 and 1 respectively

                *phasexy* : np.array(nf,4)
                            TE phase and error as row 0 and 1 respectively

                *phaseyx* : np.array(nf,4)
                            Tm phase and error as row 0 and 1 respectively

                *realtip* : np.array(nf,4)
                            Real Tipper and error as row 0 and 1 respectively

                *imagtip* : np.array(nf,4)
                            Imaginary Tipper and error as row 0 and 1 
                            respectively

                Note: that the resistivity will be in log10 space.  Also, there
                are 2 extra rows in the data arrays, this is to put the 
                response from the inversion.  

        :Example: ::

            >>> import mtpy.modeling.occamtools as occam
            >>> ocd = occam.Occam2DData()
            >>> ocd.read2DRespFile(r"/home/Occam2D/Line1/Inv1/Test_15.resp")
        """
        # make the response file an attribute
        self.respfn = respfn

        # read in the current data file
        self.read2DdataFile()

        rfid = open(self.respfn, 'r')

        rlines = rfid.readlines()
        for line in rlines:
            ls = line.split()
            # station index
            ss = int(float(ls[0])) - 1
            # component key
            comp = str(int(float(ls[2])))
            # frequency index
            ff = int(float(ls[1])) - 1
            # put into array
            # model response
            self.rplst[ss][occamdict[comp]][2, ff] = float(ls[5])
            # relative error
            self.rplst[ss][occamdict[comp]][3, ff] = float(ls[6])

    def plot2DResponses(self, respfn=None, wlfn=None, maxcol=8, plottype='1',
                        ms=2, fs=10, phaselimits=(-5, 95), colormode='color',
                        reslimits=None, plotnum=2, **kwargs):
        """
        plotResponse will plot the responses modeled from winglink against the 
        observed data.

        Arguments:
        ----------

            **respfn** : string
                         full path to response file

            **wlfn** : string
                       full path to a winglink data file used for a similar
                       inversion.  This will be plotted on the response
                       plots for comparison of fits.
                       *Default* is None

            **maxcol** : int
                         maximum number of columns for the plot
                         *Default* is 8

            **plottype** : string
                          * 'all' to plot all on the same plot
                          * '1' to plot each respones in a different figure
                               station to plot a single station 
                          * or enter as a list of stations to plot a few 
                            stations [station1,station2]. Does not have to be 
                            verbatim but should have similar unique characters
                            input pb01 for pb01cs in outputfile

                          * *Default* is '1' to plot all stations

            **ms** : float
                     marker size
                     *Default* is 2

            **phaselimits** : tuple (min,max)
                              limits of phase in degrees (min,max)
                              *Default* is (-5,95)

            **colormode** : string
                            * 'color' for color plots
                            * 'bw' for black and white plots

            **reslimits** : tuple (min,max)
                            resistivity limits on a log scale 
                            (log10(min),log10(max))
                            *Default* is None

            **plotnum** : int
                            * 1 to plot both TE and TM in the same plot
                            * 2 to plot TE and TM in separate subplots
                            * *Default* is 2

        :Example: ::

            >>> import mtpy.modeling.occamtools as occam
            >>> ocd = occam.Occam2DData()
            >>> rfile = r"/home/Occam2D/Line1/Inv1/Test_15.resp"
            >>> ocd.datafn = r"/home/Occam2D/Line1/Inv1/DataRW.dat"
            >>>
            >>> #plot all responses in their own figure and modes in separate 
            >>> #suplots
            >>> ocd.plot2DResponses(respfn=rfile)

        """

        plt.rcParams['font.size'] = fs - 2

        try:
            dpi = kwargs['dpi']
        except KeyError:
            dpi = 200

        # color mode
        if colormode == 'color':
            # color for data
            cted = (0, 0, 1)
            ctmd = (1, 0, 0)
            mted = 's'
            mtmd = 'o'

            # color for occam model
            ctem = (0, .6, .3)
            ctmm = (.9, 0, .8)
            mtem = '+'
            mtmm = '+'

            # color for Winglink model
            ctewl = (0, .6, .8)
            ctmwl = (.8, .7, 0)
            mtewl = 'x'
            mtmwl = 'x'

        # black and white mode
        elif colormode == 'bw':
            # color for data
            cted = (0, 0, 0)
            ctmd = (0, 0, 0)
            mted = '*'
            mtmd = 'v'

            # color for occam model
            ctem = (0.6, .6, .6)
            ctmm = (.6, .6, .6)
            mtem = '+'
            mtmm = 'x'

            # color for Wingling model
            ctewl = (.3, .3, .3)
            ctmwl = (.3, .3, .3)
            mtewl = '|'
            mtmwl = '_'

        # if there is a response file to plot
        if respfn != None:
            # read in the data
            self.read2DRespFile(respfn)
            # boolean for plotting response
            plotresp = True
        else:
            # read in current data file
            self.read2DdataFile()
            # boolean for plotting response
            plotresp = False

        legend_keys = ['teo', 'tem', 'tmo', 'tmm', 'wlte', 'wltm']

        # make a local copy of the rplst
        rplst = list(self.rplst)

        # boolean for adding winglink output to the plots 0 for no, 1 for yes
        addwl = 0
        hspace = .15
        # read in winglink data file
        if wlfn != None:
            addwl = 1
            hspace = .25
            wld, wlrplst, wlplst, wlslst, wltlst = wlt.readOutputFile(wlfn)
            sdict = dict([(ostation, wlstation) for wlstation in wlslst
                          for ostation in self.stationlst
                          if wlstation.find(ostation) >= 0])

        # set a local parameter period for less typing
        period = self.period

        #---------------plot each respones in a different figure---------------
        if plottype == '1':

            # set the grid of subplots
            if plotnum == 1:
                gs = gridspec.GridSpec(6, 2, wspace=.1, left=.09, top=.93, bottom=.1,
                                       hspace=hspace)
            elif plotnum == 2:
                gs = gridspec.GridSpec(6, 2, wspace=.1, left=.07, top=.93, bottom=.1,
                                       hspace=hspace)
            # loop over each station
            for ii, station in enumerate(self.stationlst):

                rlst = []
                llst = []

#                rmslst=np.hstack((rplst[ii]['resxy'][3],
#                                           rplst[ii]['resyx'][3],
#                                            rplst[ii]['phasexy'][3],
#                                            rplst[ii]['phaseyx'][3]))
#                rms=np.sqrt(np.sum(ms**2 for ms in rmslst)/len(rmslst))
                # get the RMS values for each TE and TM modes separately
                rmslstte = np.hstack((rplst[ii]['resxy'][3],
                                      rplst[ii]['phasexy'][3]))
                rmslsttm = np.hstack((rplst[ii]['resyx'][3],
                                      rplst[ii]['phaseyx'][3]))
                rmste = np.sqrt(
                    np.sum(rms**2 for rms in rmslstte) / len(rmslstte))
                rmstm = np.sqrt(
                    np.sum(rms**2 for rms in rmslsttm) / len(rmslsttm))

                fig = plt.figure(ii + 1, [9, 10], dpi=dpi)
                plt.clf()

                # set subplot instances
                # plot both TE and TM in same subplot
                if plotnum == 1:
                    axrte = fig.add_subplot(gs[:4, :])
                    axrtm = axrte
                    axpte = fig.add_subplot(gs[-2:, :], sharex=axrte)
                    axptm = axpte

                # plot TE and TM in separate subplots
                elif plotnum == 2:
                    axrte = fig.add_subplot(gs[:4, 0])
                    axrtm = fig.add_subplot(gs[:4, 1])
                    axpte = fig.add_subplot(gs[-2:, 0], sharex=axrte)
                    axptm = fig.add_subplot(gs[-2:, 1], sharex=axrtm)

                # Plot Resistivity

                # cut out missing data points first
                rxy = np.where(rplst[ii]['resxy'][0] != 0)[0]
                ryx = np.where(rplst[ii]['resyx'][0] != 0)[0]

                # check to see if there is a xy component (TE Mode)
                if len(rxy) > 0:
                    rte = axrte.errorbar(period[rxy],
                                         10**rplst[ii]['resxy'][0][rxy],
                                         ls=':', marker=mted, ms=ms, mfc=cted,
                                         mec=cted, color=cted,
                                         yerr=np.log(10) * rplst[ii]['resxy'][1][rxy] *
                                         10**rplst[ii]['resxy'][0][rxy],
                                         ecolor=cted, picker=2)
                    rlst.append(rte[0])
                    llst.append('$Obs_{TE}$')
                else:
                    pass

                # check to see if there is a yx component (TM Mode)
                if len(ryx) > 0:
                    rtm = axrtm.errorbar(period[ryx],
                                         10**rplst[ii]['resyx'][0][ryx],
                                         ls=':', marker=mtmd, ms=ms, mfc=ctmd,
                                         mec=ctmd, color=ctmd,
                                         yerr=np.log(10) * rplst[ii]['resyx'][1][ryx] *
                                         10**rplst[ii]['resyx'][0][ryx],
                                         ecolor=ctmd, picker=2)
                    rlst.append(rtm[0])
                    llst.append('$Obs_{TM}$')
                else:
                    pass

                # plot phase
                # cut out missing data points first
                pxy = np.where(rplst[ii]['phasexy'][0] != 0)[0]
                pyx = np.where(rplst[ii]['phaseyx'][0] != 0)[0]

                # plot the xy component (TE Mode)
                if len(pxy) > 0:
                    axpte.errorbar(period[pxy], rplst[ii]['phasexy'][0][pxy],
                                   ls=':', marker=mted, ms=ms, mfc=cted, mec=cted,
                                   color=cted,
                                   yerr=rplst[ii]['phasexy'][1][pxy],
                                   ecolor=cted, picker=1)
                else:
                    pass

                # plot the yx component (TM Mode)
                if len(pyx) > 0:
                    axptm.errorbar(period[pyx], rplst[ii]['phaseyx'][0][pyx],
                                   ls=':', marker=mtmd, ms=ms, mfc=ctmd, mec=ctmd,
                                   color=ctmd,
                                   yerr=rplst[ii]['phaseyx'][1][pyx],
                                   ecolor=ctmd, picker=1)
                else:
                    pass

                # if there is a response file
                if plotresp == True:
                    mrxy = np.where(rplst[ii]['resxy'][2] != 0)[0]
                    mryx = np.where(rplst[ii]['resyx'][2] != 0)[0]

                    # plot the Model Resistivity
                    # check for the xy of model component
                    if len(mrxy) > 0:
                        r3 = axrte.errorbar(period[mrxy],
                                            10**rplst[ii]['resxy'][2][mrxy],
                                            ls='--', marker=mtem, ms=ms, mfc=ctem,
                                            mec=ctem, color=ctem,
                                            yerr=10**(rplst[ii]['resxy'][3][mrxy] *
                                                      rplst[ii]['resxy'][2][mrxy] / np.log(10)),
                                            ecolor=ctem)
                        rlst.append(r3[0])
                        llst.append('$Mod_{TE}$')
                    else:
                        pass

                    # check for the yx model component  of resisitivity
                    if len(mryx) > 0:
                        r4 = axrtm.errorbar(period[mryx],
                                            10**rplst[ii]['resyx'][2][mryx],
                                            ls='--', marker=mtmm, ms=ms, mfc=ctmm,
                                            mec=ctmm, color=ctmm,
                                            yerr=10**(rplst[ii]['resyx'][3][mryx] *
                                                      rplst[ii]['resyx'][2][mryx] / np.log(10)),
                                            ecolor=ctmm)
                        rlst.append(r4[0])
                        llst.append('$Mod_{TM}$')

                    # plot the model phase
                    # check for removed points
                    mpxy = np.where(rplst[ii]['phasexy'][2] != 0)[0]
                    mpyx = np.where(rplst[ii]['phaseyx'][2] != 0)[0]

                    # plot the xy component (TE Mode)
                    if len(mpxy) > 0:
                        axpte.errorbar(period[mpxy],
                                       rplst[ii]['phasexy'][2][mpxy],
                                       ls='--', marker=mtem, ms=ms, mfc=ctem,
                                       mec=ctem, color=ctem,
                                       yerr=rplst[ii]['phasexy'][3][mpxy],
                                       ecolor=ctem)
                    else:
                        pass

                    # plot the yx component (TM Mode)
                    if len(mpyx) > 0:
                        axptm.errorbar(period[mpyx],
                                       rplst[ii]['phaseyx'][2][mpyx],
                                       ls='--', marker=mtmm, ms=ms, mfc=ctmm,
                                       mec=ctmm, color=ctmm,
                                       yerr=rplst[ii]['phaseyx'][3][mpyx],
                                       ecolor=ctmm)
                    else:
                        pass

                # add in winglink responses
                if addwl == 1:
                    try:
                        wlrms = wld[sdict[station]]['rms']
                        axr.set_title(self.stationlst[ii] + '\n' +
                                      'rms_occ_TE={0:.2f}, rms_occ_TM={1:.2f}, rms_wl= {2:.2f}'.format(
                                          rmste, rmstm, wlrms),
                                      fontdict={'size': fs + 1, 'weight': 'bold'})
                        for ww, wlstation in enumerate(wlslst):
                            #                        print station,wlstation
                            if wlstation.find(station) == 0:
                                print(station, wlstation)
                                wlrpdict = wlrplst[ww]

                        zrxy = [np.where(wlrpdict['resxy'][0] != 0)[0]]
                        zryx = [np.where(wlrpdict['resyx'][0] != 0)[0]]

                        # plot winglink resistivity
                        r5 = axrte.loglog(wlplst[zrxy],
                                          wlrpdict['resxy'][1][zrxy],
                                          ls='-.', marker=mtewl, ms=5, color=ctewl,
                                          mfc=ctewl)
                        r6 = axrtm.loglog(wlplst[zryx],
                                          wlrpdict['resyx'][1][zryx],
                                          ls='-.', marker=mtmwl, ms=5, color=ctmwl,
                                          mfc=ctmwl)

                        # plot winglink phase
                        axpte.semilogx(wlplst[zrxy],
                                       wlrpdict['phasexy'][1][zrxy],
                                       ls='-.', marker=mtewl, ms=5, color=ctewl,
                                       mfc=ctewl)
                        axptm.semilogx(wlplst[zryx],
                                       wlrpdict['phaseyx'][1][zryx],
                                       ls='-.', marker=mtmwl, ms=5, color=ctmwl,
                                       mfc=ctmwl)

                        rlst.append(r5[0])
                        rlst.append(r6[0])
                        llst.append('$WLMod_{TE}$')
                        llst.append('$WLMod_{TM}$')
                    except IndexError:
                        print('Station not present')
                else:
                    if plotnum == 1:
                        axrte.set_title(self.stationlst[ii] +
                                        ' rms_TE={0:.2f}, rms_TM={1:.2f}'.format(
                                            rmste, rmstm),
                                        fontdict={'size': fs + 1, 'weight': 'bold'})
                    elif plotnum == 2:
                        axrte.set_title(self.stationlst[ii] +
                                        ' rms_TE={0:.2f}'.format(rmste),
                                        fontdict={'size': fs + 1, 'weight': 'bold'})
                        axrtm.set_title(self.stationlst[ii] +
                                        ' rms_TM={0:.2f}'.format(rmstm),
                                        fontdict={'size': fs + 1, 'weight': 'bold'})

                # set the axis properties
                for aa, axr in enumerate([axrte, axrtm]):
                    # set both axes to logarithmic scale
                    axr.set_xscale('log', nonposx='clip')
                    axr.set_yscale('log', nonposy='clip')

                    # put on a grid
                    axr.grid(True, alpha=.3, which='both')
                    axr.yaxis.set_label_coords(-.07, .5)

                    # set resistivity limits if desired
                    if reslimits != None:
                        axr.set_ylim(10**reslimits[0], 10**reslimits[1])

                    # set the tick labels to invisible
                    plt.setp(axr.xaxis.get_ticklabels(), visible=False)
                    if aa == 0:
                        axr.set_ylabel('App. Res. ($\Omega \cdot m$)',
                                       fontdict={'size': fs, 'weight': 'bold'})

                    # set legend based on the plot type
                    if plotnum == 1:
                        if aa == 0:
                            axr.legend(rlst, llst,
                                       loc=2, markerscale=1,
                                       borderaxespad=.05,
                                       labelspacing=.08,
                                       handletextpad=.15,
                                       borderpad=.05,
                                       prop={'size': fs})
                    elif plotnum == 2:
                        if aa == 0:
                            if plotresp == True:
                                try:
                                    axr.legend([rlst[0], rlst[2]],
                                               [llst[0], llst[2]],
                                               loc=2, markerscale=1,
                                               borderaxespad=.05,
                                               labelspacing=.08,
                                               handletextpad=.15,
                                               borderpad=.05,
                                               prop={'size': fs})
                                except IndexError:
                                    pass

                            else:
                                try:
                                    axr.legend([rlst[0]], [llst[0]],
                                               loc=2, markerscale=1,
                                               borderaxespad=.05,
                                               labelspacing=.08,
                                               handletextpad=.15,
                                               borderpad=.05,
                                               prop={'size': fs})
                                except IndexError:
                                    pass
                        if aa == 1:
                            if plotresp == True:
                                try:
                                    axr.legend([rlst[1], rlst[3]],
                                               [llst[1], llst[3]],
                                               loc=2, markerscale=1,
                                               borderaxespad=.05,
                                               labelspacing=.08,
                                               handletextpad=.15, borderpad=.05,
                                               prop={'size': fs})
                                except IndexError:
                                    pass
                            else:
                                try:
                                    axr.legend([rlst[1]], [llst[1]],
                                               loc=2, markerscale=1,
                                               borderaxespad=.05,
                                               labelspacing=.08,
                                               handletextpad=.15,
                                               borderpad=.05,
                                               prop={'size': fs})
                                except IndexError:
                                    pass

                # set Properties for the phase axes
                for aa, axp in enumerate([axpte, axptm]):
                    # set the x-axis to log scale
                    axp.set_xscale('log', nonposx='clip')

                    # set the phase limits
                    axp.set_ylim(phaselimits)

                    # put a grid on the subplot
                    axp.grid(True, alpha=.3, which='both')

                    # set the tick locations
                    axp.yaxis.set_major_locator(MultipleLocator(10))
                    axp.yaxis.set_minor_locator(MultipleLocator(2))

                    # set the x axis label
                    axp.set_xlabel('Period (s)',
                                   fontdict={'size': fs, 'weight': 'bold'})

                    # put the y label on the far left plot
                    axp.yaxis.set_label_coords(-.07, .5)
                    if aa == 0:
                        axp.set_ylabel('Phase (deg)',
                                       fontdict={'size': fs, 'weight': 'bold'})

            # set the plot to be full screen well at least try
            plt.show()

        #---Plot single or subset of stations----------------------------------
        else:
            pstationlst = []

            if type(plottype) is not list:
                plottype = [plottype]
            for ii, station in enumerate(self.stationlst):
                for pstation in plottype:
                    if station.find(pstation) >= 0:
                        #                    print 'plotting ',station
                        pstationlst.append(ii)
            if addwl == 1:
                pwlstationlst = []
                for ww, wlstation in enumerate(wlslst):
                    for pstation in plottype:
                        if wlstation.find(pstation) >= 0:
                            # print 'plotting ',wlstation
                            pwlstationlst.append(ww)
            if plotnum == 1:
                gs = gridspec.GridSpec(6, 2, wspace=.1, left=.09, top=.93, bottom=.1,
                                       hspace=hspace)
            elif plotnum == 2:
                gs = gridspec.GridSpec(6, 2, wspace=.1, left=.07, top=.93, bottom=.1,
                                       hspace=hspace)
            for jj, ii in enumerate(pstationlst):
                legend_dict = dict([(lkey, None) for lkey in legend_keys])
                legend_label = dict([(lkey, None) for lkey in legend_keys])

                # get RMS values for TE and TM separately
                rmslstte = np.hstack((rplst[ii]['resxy'][3],
                                      rplst[ii]['phasexy'][3]))
                rmslsttm = np.hstack((rplst[ii]['resyx'][3],
                                      rplst[ii]['phaseyx'][3]))
                rmste = np.sqrt(
                    np.sum(rms**2 for rms in rmslstte) / len(rmslstte))
                rmstm = np.sqrt(
                    np.sum(rms**2 for rms in rmslsttm) / len(rmslsttm))

                fig = plt.figure(ii + 1, [9, 10], dpi=dpi)
                plt.clf()

                # set subplot instances
                # plot both TE and TM in same subplot
                if plotnum == 1:
                    axrte = fig.add_subplot(gs[:4, :])
                    axrtm = axrte
                    axpte = fig.add_subplot(gs[-2:, :], sharex=axrte)
                    axptm = axpte

                # plot TE and TM in separate subplots
                elif plotnum == 2:
                    axrte = fig.add_subplot(gs[:4, 0])
                    axrtm = fig.add_subplot(gs[:4, 1])
                    axpte = fig.add_subplot(gs[-2:, 0], sharex=axrte)
                    axptm = fig.add_subplot(gs[-2:, 1], sharex=axrtm)

                # Plot Resistivity

                # cut out missing data points first
                rxy = np.where(rplst[ii]['resxy'][0] != 0)[0]
                ryx = np.where(rplst[ii]['resyx'][0] != 0)[0]

                # check to see if there is a xy component (TE Mode)
                if len(rxy) > 0:
                    rte = axrte.errorbar(period[rxy],
                                         10**rplst[ii]['resxy'][0][rxy],
                                         ls=':', marker=mted, ms=ms, mfc=cted,
                                         mec=cted, color=cted,
                                         yerr=np.log(10) * rplst[ii]['resxy'][1][rxy] *
                                         10**rplst[ii]['resxy'][0][rxy],
                                         ecolor=cted, picker=2)
                    legend_dict['teo'] = rte[0]
                    legend_label['teo'] = '$Obs_{TE}$'
                else:
                    pass

                # check to see if there is a yx component (TM Mode)
                if len(ryx) > 0:
                    rtm = axrtm.errorbar(period[ryx],
                                         10**rplst[ii]['resyx'][0][ryx],
                                         ls=':', marker=mtmd, ms=ms, mfc=ctmd,
                                         mec=ctmd, color=ctmd,
                                         yerr=np.log(10) * rplst[ii]['resyx'][1][ryx] *
                                         10**rplst[ii]['resyx'][0][ryx],
                                         ecolor=ctmd, picker=2)
                    legend_dict['tmo'] = rtm[0]
                    legend_label['tmo'] = '$Obs_{TM}$'

                else:
                    pass

                # plot phase
                # cut out missing data points first
                pxy = np.where(rplst[ii]['phasexy'][0] != 0)[0]
                pyx = np.where(rplst[ii]['phaseyx'][0] != 0)[0]

                # plot the xy component (TE Mode)
                if len(pxy) > 0:
                    axpte.errorbar(period[pxy], rplst[ii]['phasexy'][0][pxy],
                                   ls=':', marker=mted, ms=ms, mfc=cted, mec=cted,
                                   color=cted,
                                   yerr=rplst[ii]['phasexy'][1][pxy],
                                   ecolor=cted, picker=1)
                else:
                    pass

                # plot the yx component (TM Mode)
                if len(pyx) > 0:
                    axptm.errorbar(period[pyx], rplst[ii]['phaseyx'][0][pyx],
                                   ls=':', marker=mtmd, ms=ms, mfc=ctmd, mec=ctmd,
                                   color=ctmd,
                                   yerr=rplst[ii]['phaseyx'][1][pyx],
                                   ecolor=ctmd, picker=1)
                else:
                    pass

                # if there is a response file
                if plotresp == True:
                    mrxy = np.where(rplst[ii]['resxy'][2] != 0)[0]
                    mryx = np.where(rplst[ii]['resyx'][2] != 0)[0]

                    # plot the Model Resistivity
                    # check for the xy of model component
                    if len(mrxy) > 0:
                        r3 = axrte.errorbar(period[mrxy],
                                            10**rplst[ii]['resxy'][2][mrxy],
                                            ls='--', marker=mtem, ms=ms, mfc=ctem,
                                            mec=ctem, color=ctem,
                                            yerr=10**(rplst[ii]['resxy'][3][mrxy] *
                                                      rplst[ii]['resxy'][2][mrxy] / np.log(10)),
                                            ecolor=ctem)
                        legend_dict['tem'] = r3[0]
                        legend_label['tem'] = '$Mod_{TE}$'
                    else:
                        pass

                    # check for the yx model component  of resisitivity
                    if len(mryx) > 0:
                        r4 = axrtm.errorbar(period[mryx],
                                            10**rplst[ii]['resyx'][2][mryx],
                                            ls='--', marker=mtmm, ms=ms, mfc=ctmm,
                                            mec=ctmm, color=ctmm,
                                            yerr=10**(rplst[ii]['resyx'][3][mryx] *
                                                      rplst[ii]['resyx'][2][mryx] / np.log(10)),
                                            ecolor=ctmm)
                        legend_dict['tmm'] = r4[0]
                        legend_label['tmm'] = '$Mod_{TM}$'

                    # plot the model phase
                    # check for removed points
                    mpxy = np.where(rplst[ii]['phasexy'][2] != 0)[0]
                    mpyx = np.where(rplst[ii]['phaseyx'][2] != 0)[0]

                    # plot the xy component (TE Mode)
                    if len(mpxy) > 0:
                        axpte.errorbar(period[mpxy],
                                       rplst[ii]['phasexy'][2][mpxy],
                                       ls='--', marker=mtem, ms=ms, mfc=ctem,
                                       mec=ctem, color=ctem,
                                       yerr=rplst[ii]['phasexy'][3][mpxy],
                                       ecolor=ctem)
                    else:
                        pass

                    # plot the yx component (TM Mode)
                    if len(mpyx) > 0:
                        axptm.errorbar(period[mpyx],
                                       rplst[ii]['phaseyx'][2][mpyx],
                                       ls='--', marker=mtmm, ms=ms, mfc=ctmm,
                                       mec=ctmm, color=ctmm,
                                       yerr=rplst[ii]['phaseyx'][3][mpyx],
                                       ecolor=ctmm)
                    else:
                        pass

                # add in winglink responses
                if addwl == 1:
                    try:
                        wlrms = wld[sdict[station]]['rms']
                        axr.set_title(self.stationlst[ii] + '\n' +
                                      ' rms_occ_TE={0:.2f}, rms_occ_TM={1:.2f}, rms_wl= {2:.2f}'.format(
                                          rmste, rmstm, wlrms),
                                      fontdict={'size': fs + 1, 'weight': 'bold'})
                        for ww, wlstation in enumerate(wlslst):
                            #                        print station,wlstation
                            if wlstation.find(station) == 0:
                                print(station, wlstation)
                                wlrpdict = wlrplst[ww]

                        zrxy = [np.where(wlrpdict['resxy'][0] != 0)[0]]
                        zryx = [np.where(wlrpdict['resyx'][0] != 0)[0]]

                        # plot winglink resistivity
                        r5 = axrte.loglog(wlplst[zrxy],
                                          wlrpdict['resxy'][1][zrxy],
                                          ls='-.', marker=mtewl, ms=5, color=ctewl,
                                          mfc=ctewl)
                        r6 = axrtm.loglog(wlplst[zryx],
                                          wlrpdict['resyx'][1][zryx],
                                          ls='-.', marker=mtmwl, ms=5, color=ctmwl,
                                          mfc=ctmwl)

                        # plot winglink phase
                        axpte.semilogx(wlplst[zrxy],
                                       wlrpdict['phasexy'][1][zrxy],
                                       ls='-.', marker=mtewl, ms=5, color=ctewl,
                                       mfc=ctewl)
                        axptm.semilogx(wlplst[zryx],
                                       wlrpdict['phaseyx'][1][zryx],
                                       ls='-.', marker=mtmwl, ms=5, color=ctmwl,
                                       mfc=ctmwl)

                        legend_dict['wlte'] = r5[0]
                        legend_label['wlte'] = '$WLMod_{TE}$'
                        legend_dict['wltm'] = r6[0]
                        legend_label['wltm'] = '$WLMod_{TM}$'
                    except IndexError:
                        print('Station not present')
                else:
                    if plotnum == 1:
                        axrte.set_title(self.stationlst[ii] +
                                        ' rms_TE={0:.2f}, rms_TM={1:.2f}'.format(
                                            rmste, rmstm),
                                        fontdict={'size': fs + 1, 'weight': 'bold'})
                    elif plotnum == 2:
                        axrte.set_title(self.stationlst[ii] +
                                        ' rms_TE={0:.2f}'.format(rmste),
                                        fontdict={'size': fs + 1, 'weight': 'bold'})
                        axrtm.set_title(self.stationlst[ii] +
                                        ' rms_TM={0:.2f}'.format(rmstm),
                                        fontdict={'size': fs + 1, 'weight': 'bold'})

                # set the axis properties
                for aa, axr in enumerate([axrte, axrtm]):
                    # set both axes to logarithmic scale
                    axr.set_xscale('log', nonposx='clip')
                    axr.set_yscale('log', nonposy='clip')

                    # put on a grid
                    axr.grid(True, alpha=.3, which='both')
                    axr.yaxis.set_label_coords(-.07, .5)

                    # set resistivity limits if desired
                    if reslimits != None:
                        axr.set_ylim(10**reslimits[0], 10**reslimits[1])

                    # set the tick labels to invisible
                    plt.setp(axr.xaxis.get_ticklabels(), visible=False)
                    if aa == 0:
                        axr.set_ylabel('App. Res. ($\Omega \cdot m$)',
                                       fontdict={'size': fs, 'weight': 'bold'})
                    if plotnum == 1:
                        if aa == 0:
                            rlst = [legend_dict[lkey] for lkey in legend_keys
                                    if legend_dict[lkey] != None]
                            llst = [legend_label[lkey] for lkey in legend_keys
                                    if legend_label[lkey] != None]
                            axr.legend(rlst, llst,
                                       loc=2, markerscale=1, borderaxespad=.05,
                                       labelspacing=.08,
                                       handletextpad=.15, borderpad=.05, prop={'size': fs})
                    elif plotnum == 2:
                        if aa == 0:
                            rlst = [legend_dict[lkey] for lkey in legend_keys
                                    if legend_dict[lkey] != None and
                                    lkey.find('te') >= 0]
                            llst = [legend_label[lkey] for lkey in legend_keys
                                    if legend_label[lkey] != None and
                                    lkey.find('te') >= 0]
                            if plotresp == True:

                                axr.legend(rlst, llst,
                                           loc=2, markerscale=1, borderaxespad=.05,
                                           labelspacing=.08,
                                           handletextpad=.15, borderpad=.05,
                                           prop={'size': fs})
                            else:
                                axr.legend(rlst, llst,
                                           loc=2, markerscale=1, borderaxespad=.05,
                                           labelspacing=.08,
                                           handletextpad=.15, borderpad=.05,
                                           prop={'size': fs})
                        if aa == 1:
                            rlst = [legend_dict[lkey] for lkey in legend_keys
                                    if legend_dict[lkey] != None and
                                    lkey.find('tm') >= 0]
                            llst = [legend_label[lkey] for lkey in legend_keys
                                    if legend_label[lkey] != None and
                                    lkey.find('tm') >= 0]
                            if plotresp == True:
                                axr.legend(rlst, llst,
                                           loc=2, markerscale=1, borderaxespad=.05,
                                           labelspacing=.08,
                                           handletextpad=.15, borderpad=.05,
                                           prop={'size': fs})
                            else:
                                axr.legend(rlst, llst,
                                           loc=2, markerscale=1, borderaxespad=.05,
                                           labelspacing=.08,
                                           handletextpad=.15, borderpad=.05,
                                           prop={'size': fs})

                for aa, axp in enumerate([axpte, axptm]):
                    # set the x-axis to log scale
                    axp.set_xscale('log', nonposx='clip')

                    # set the phase limits
                    axp.set_ylim(phaselimits)

                    # put a grid on the subplot
                    axp.grid(True, alpha=.3, which='both')

                    # set the tick locations
                    axp.yaxis.set_major_locator(MultipleLocator(10))
                    axp.yaxis.set_minor_locator(MultipleLocator(2))

                    # set the x axis label
                    axp.set_xlabel('Period (s)',
                                   fontdict={'size': fs, 'weight': 'bold'})

                    # put the y label on the far left plot
                    axp.yaxis.set_label_coords(-.07, .5)
                    if aa == 0:
                        axp.set_ylabel('Phase (deg)',
                                       fontdict={'size': fs, 'weight': 'bold'})
           # set figure to full window size
            plt.show()

    def plotPseudoSection(self, respfn=None, fignum=1, rcmap='jet_r', pcmap='jet',
                          rlim=((0, 4), (0, 4)), plim=((0, 90), (0, 90)), ml=2,
                          stationid=(0, 4)):
        """
        plots a pseudo section of the data and response if input.

        Arguments:
        ----------
            **respfn** : string
                         full path to response file

            **fignum** : int
                         figure number to put the plot in

            **rcmap** : string
                        colormap to plot the resistivity in see matplotlib.cm
                        *Default* 'jet_r' 

            **pcmap** : string
                        colormap to plot the resistivity in see matplotlib.cm
                        *Default* 'jet'
            **rlim** : tuple ((te_min,te_max),(tm_min,tm_max))
                       min and max ranges for resistivity of TE and TM modes
                       in log10 scale
                       *Default* is ((0,4),(0,4))

            **plim** : tuple ((te_min,te_max),(tm_min,tm_max))
                       min and max ranges for phase of TE and TM modes
                       *Default* is ((0,90),(0,90))

            **ml** : int
                    tick markers between stations
                    *Default* is 2 (for one every second station)

            **stationid** : tuple (min_string_index,max_string_index)
                            indicating how long the station name is
                            *Default* is (0,4) for a string of length 4

       :Example: ::

            >>> import mtpy.modeling.occamtools as occam
            >>> ocd = occam.Occam2DData()
            >>> rfile = r"/home/Occam2D/Line1/Inv1/Test_15.resp"
            >>> ocd.datafn = r"/home/Occam2D/Line1/Inv1/DataRW.dat"
            >>> ocd.plot2PseudoSection(respfn=rfile) 

        """

        try:
            self.read2DRespFile(respfn)
            nr = 2
        except TypeError:
            nr = 1

        ns = len(self.stationlst)
        nf = len(self.freq)
        ylimits = (1. / self.freq.min(), 1. / self.freq.max())
    #    print ylimits

        # make a grid for pcolormesh so you can have a log scale
        # get things into arrays for plotting
        offsetlst = np.zeros(ns)
        resxyarr = np.zeros((nf, ns, nr))
        resyxarr = np.zeros((nf, ns, nr))
        phasexyarr = np.zeros((nf, ns, nr))
        phaseyxarr = np.zeros((nf, ns, nr))

        for ii, rpdict in enumerate(self.rplst):
            offsetlst[ii] = rpdict['offset']
            resxyarr[:, ii, 0] = rpdict['resxy'][0]
            resyxarr[:, ii, 0] = rpdict['resyx'][0]
            phasexyarr[:, ii, 0] = rpdict['phasexy'][0]
            phaseyxarr[:, ii, 0] = rpdict['phaseyx'][0]
            if respfn != None:
                resxyarr[:, ii, 1] = rpdict['resxy'][2]
                resyxarr[:, ii, 1] = rpdict['resyx'][2]
                phasexyarr[:, ii, 1] = rpdict['phasexy'][2]
                phaseyxarr[:, ii, 1] = rpdict['phaseyx'][2]

        # make a meshgrid for plotting
        # flip frequency so bottom corner is long period
        dgrid, fgrid = np.meshgrid(offsetlst, 1. / self.freq[::-1])

        # make list for station labels
        slabel = [self.stationlst[ss][stationid[0]:stationid[1]]
                  for ss in range(0, ns, ml)]
        labellst = ['$r_{TE-Data}$', '$r_{TE-Model}$',
                    '$r_{TM-Data}$', '$r_{TM-Model}$',
                    '$\phi_{TE-Data}$', '$\phi_{TE-Model}$',
                    '$\phi_{TM-Data}$', '$\phi_{TM-Model}$']
        xloc = offsetlst[0] + abs(offsetlst[0] - offsetlst[1]) / 5
        yloc = 1. / self.freq[1]

        if respfn != None:

            plt.rcParams['font.size'] = 7
            plt.rcParams['figure.subplot.bottom'] = .09
            plt.rcParams['figure.subplot.top'] = .96

            fig = plt.figure(fignum, dpi=200)
            plt.clf()

            # make subplot grids
            gs1 = gridspec.GridSpec(2, 2, left=0.06, right=.48, hspace=.1,
                                    wspace=.005)
            gs2 = gridspec.GridSpec(2, 2, left=0.52, right=.98, hspace=.1,
                                    wspace=.005)

            # plot TE resistivity data
            ax1r = fig.add_subplot(gs1[0, 0])
            ax1r.pcolormesh(dgrid, fgrid, np.flipud(resxyarr[:, :, 0]), cmap=rcmap,
                            vmin=rlim[0][0], vmax=rlim[0][1])

            # plot TE resistivity model
            ax2r = fig.add_subplot(gs1[0, 1])
            ax2r.pcolormesh(dgrid, fgrid, np.flipud(resxyarr[:, :, 1]), cmap=rcmap,
                            vmin=rlim[0][0], vmax=rlim[0][1])

            # plot TM resistivity data
            ax3r = fig.add_subplot(gs2[0, 0])
            ax3r.pcolormesh(dgrid, fgrid, np.flipud(resyxarr[:, :, 0]), cmap=rcmap,
                            vmin=rlim[1][0], vmax=rlim[1][1])

            # plot TM resistivity model
            ax4r = fig.add_subplot(gs2[0, 1])
            ax4r.pcolormesh(dgrid, fgrid, np.flipud(resyxarr[:, :, 1]), cmap=rcmap,
                            vmin=rlim[1][0], vmax=rlim[1][1])

            # plot TE phase data
            ax1p = fig.add_subplot(gs1[1, 0])
            ax1p.pcolormesh(dgrid, fgrid, np.flipud(phasexyarr[:, :, 0]),
                            cmap=pcmap, vmin=plim[0][0], vmax=plim[0][1])

            # plot TE phase model
            ax2p = fig.add_subplot(gs1[1, 1])
            ax2p.pcolormesh(dgrid, fgrid, np.flipud(phasexyarr[:, :, 1]),
                            cmap=pcmap, vmin=plim[0][0], vmax=plim[0][1])

            # plot TM phase data
            ax3p = fig.add_subplot(gs2[1, 0])
            ax3p.pcolormesh(dgrid, fgrid, np.flipud(phaseyxarr[:, :, 0]),
                            cmap=pcmap, vmin=plim[1][0], vmax=plim[1][1])

            # plot TM phase model
            ax4p = fig.add_subplot(gs2[1, 1])
            ax4p.pcolormesh(dgrid, fgrid, np.flipud(phaseyxarr[:, :, 1]),
                            cmap=pcmap, vmin=plim[1][0], vmax=plim[1][1])

            axlst = [ax1r, ax2r, ax3r, ax4r, ax1p, ax2p, ax3p, ax4p]

            # make everthing look tidy
            for xx, ax in enumerate(axlst):
                ax.semilogy()
                ax.set_ylim(ylimits)
                ax.xaxis.set_ticks(offsetlst[np.arange(0, ns, ml)])
                ax.xaxis.set_ticks(offsetlst, minor=True)
                ax.xaxis.set_ticklabels(slabel)
                ax.set_xlim(offsetlst.min(), offsetlst.max())
                if np.remainder(xx, 2.0) == 1:
                    plt.setp(ax.yaxis.get_ticklabels(), visible=False)
                    cbx = mcb.make_axes(ax, shrink=.7, pad=.015)
                    if xx < 4:
                        if xx == 1:
                            cb = mcb.ColorbarBase(cbx[0], cmap=rcmap,
                                                  norm=Normalize(vmin=rlim[0][0],
                                                                 vmax=rlim[0][1]))
                        if xx == 3:
                            cb = mcb.ColorbarBase(cbx[0], cmap=rcmap,
                                                  norm=Normalize(vmin=rlim[1][0],
                                                                 vmax=rlim[1][1]))
                            cb.set_label('App. Res. ($\Omega \cdot$m)',
                                         fontdict={'size': 9})
                    else:
                        if xx == 5:
                            cb = mcb.ColorbarBase(cbx[0], cmap=pcmap,
                                                  norm=Normalize(vmin=plim[0][0],
                                                                 vmax=plim[0][1]))
                        if xx == 7:
                            cb = mcb.ColorbarBase(cbx[0], cmap=pcmap,
                                                  norm=Normalize(vmin=plim[1][0],
                                                                 vmax=plim[1][1]))
                            cb.set_label('Phase (deg)', fontdict={'size': 9})
                ax.text(xloc, yloc, labellst[xx],
                        fontdict={'size': 10},
                        bbox={'facecolor': 'white'},
                        horizontalalignment='left',
                        verticalalignment='top')
                if xx == 0 or xx == 4:
                    ax.set_ylabel('Period (s)',
                                  fontdict={'size': 10, 'weight': 'bold'})
                if xx > 3:
                    ax.set_xlabel('Station', fontdict={'size': 10,
                                                       'weight': 'bold'})

            plt.show()

        else:

            plt.rcParams['font.size'] = 7
            plt.rcParams['figure.subplot.bottom'] = .09
            plt.rcParams['figure.subplot.top'] = .96

            fig = plt.figure(fignum, dpi=200)
            plt.clf()

            # make subplot grids
            gs1 = gridspec.GridSpec(2, 2, left=0.06, right=.48, hspace=.1,
                                    wspace=.005)
            gs2 = gridspec.GridSpec(2, 2, left=0.52, right=.98, hspace=.1,
                                    wspace=.005)

            # plot TE resistivity data
            ax1r = fig.add_subplot(gs1[0, :])
            ax1r.pcolormesh(dgrid, fgrid, np.flipud(resxyarr[:, :, 0]), cmap=rcmap,
                            vmin=rlim[0][0], vmax=rlim[0][1])

            # plot TM resistivity data
            ax3r = fig.add_subplot(gs2[0, :])
            ax3r.pcolormesh(dgrid, fgrid, np.flipud(resyxarr[:, :, 0]), cmap=rcmap,
                            vmin=rlim[1][0], vmax=rlim[1][1])

            # plot TE phase data
            ax1p = fig.add_subplot(gs1[1, :])
            ax1p.pcolormesh(dgrid, fgrid, np.flipud(phasexyarr[:, :, 0]), cmap=pcmap,
                            vmin=plim[0][0], vmax=plim[0][1])

            # plot TM phase data
            ax3p = fig.add_subplot(gs2[1, :])
            ax3p.pcolormesh(dgrid, fgrid, np.flipud(phaseyxarr[:, :, 0]), cmap=pcmap,
                            vmin=plim[1][0], vmax=plim[1][1])

            axlst = [ax1r, ax3r, ax1p, ax3p]

            # make everything look tidy
            for xx, ax in enumerate(axlst):
                ax.semilogy()
                ax.set_ylim(ylimits)
                ax.xaxis.set_ticks(offsetlst[np.arange(0, ns, ml)])
                ax.xaxis.set_ticks(offsetlst, minor=True)
                ax.xaxis.set_ticklabels(slabel)
                ax.set_xlim(offsetlst.min(), offsetlst.max())
                plt.setp(ax.yaxis.get_ticklabels(), visible=False)
                cbx = mcb.make_axes(ax, shrink=.7, pad=.015)
                if xx == 0:
                    cb = mcb.ColorbarBase(cbx[0], cmap=rcmap,
                                          norm=Normalize(vmin=rlim[0][0],
                                                         vmax=rlim[0][1]))
                elif xx == 1:
                    cb = mcb.ColorbarBase(cbx[0], cmap=rcmap,
                                          norm=Normalize(vmin=rlim[1][0],
                                                         vmax=rlim[1][1]))
                    cb.set_label('App. Res. ($\Omega \cdot$m)',
                                 fontdict={'size': 9})
                elif xx == 2:
                    cb = mcb.ColorbarBase(cbx[0], cmap=pcmap,
                                          norm=Normalize(vmin=plim[0][0],
                                                         vmax=plim[0][1]))
                elif xx == 3:
                    cb = mcb.ColorbarBase(cbx[0], cmap=pcmap,
                                          norm=Normalize(vmin=plim[1][0],
                                                         vmax=plim[1][1]))
                    cb.set_label('Phase (deg)', fontdict={'size': 9})
                ax.text(xloc, yloc, labellst[xx],
                        fontdict={'size': 10},
                        bbox={'facecolor': 'white'},
                        horizontalalignment='left',
                        verticalalignment='top')
                if xx == 0 or xx == 2:
                    ax.set_ylabel('Period (s)',
                                  fontdict={'size': 10, 'weight': 'bold'})
                if xx > 1:
                    ax.set_xlabel('Station', fontdict={'size': 10,
                                                       'weight': 'bold'})

            plt.show()

    def plotAllResponses(self, station, fignum=1):
        """
        Plot all the responses of occam inversion from data file.  This assumes
        the response curves are in the same folder as the datafile.

        Arguments:
        ----------
            **station** : string
                          station name to plot

            **fignum** : int
                         plot number to put figure into

        :Example: ::

            >>> import mtpy.modeling.occamtools as occam
            >>> ocd = occam.Occam2DData()
            >>> ocd.datafn = r"/home/Occam2D/Line1/Inv1/DataRW.dat"
            >>> ocd.plotAllResponses('MT01')

        """

        rpath = os.path.dirname(self.datafn)

        gs = gridspec.GridSpec(6, 2, wspace=.20)

        plt.rcParams['font.size'] = int(7)
        plt.rcParams['figure.subplot.left'] = .08
        plt.rcParams['figure.subplot.right'] = .98
        plt.rcParams['figure.subplot.bottom'] = .1
        plt.rcParams['figure.subplot.top'] = .92

        rlst = [os.path.join(rpath, rfile) for rfile in os.listdir(rpath)
                if rfile.find('.resp') > 0]

        nresp = len(rlst)

        colorlst = [(cc, 0, 1 - cc) for cc in np.arange(0, 1, 1. / nresp)]
        fig = plt.figure(fignum, [7, 8], dpi=200)
        plt.clf()
        axrte = fig.add_subplot(gs[:4, 0])
        axrtm = fig.add_subplot(gs[:4, 1])
        axpte = fig.add_subplot(gs[-2:, 0])
        axptm = fig.add_subplot(gs[-2:, 1])
        rmstelst = []
        rmstmlst = []
        rmstestr = []
        rmstmstr = []
        # read responses
        for jj, rfile in enumerate(rlst):
            respfn = os.path.join(rpath, rfile)
            self.read2DRespFile(respfn)

            ii = np.where(np.array(self.stationlst) == station)[0][0]

            period = 1. / self.freq

            rmslstte = np.hstack((self.rplst[ii]['resxy'][3],
                                  self.rplst[ii]['phasexy'][3]))
            rmslsttm = np.hstack((self.rplst[ii]['resyx'][3],
                                  self.rplst[ii]['phaseyx'][3]))
            rmste = np.sqrt(np.sum(ms**2 for ms in rmslstte) / len(rmslstte))
            rmstm = np.sqrt(np.sum(ms**2 for ms in rmslsttm) / len(rmslsttm))
            rmstelst.append('%d rms=%.3f ' % (jj, rmste))
            rmstmlst.append('%d rms=%.3f ' % (jj, rmstm))
            rmstestr.append(rmste)
            rmstmstr.append(rmstm)
            # plot resistivity

            if jj == 0:
                # cut out missing data points first
                rxy = np.where(self.rplst[ii]['resxy'][0] != 0)[0]
                ryx = np.where(self.rplst[ii]['resyx'][0] != 0)[0]
                r1, = axrte.loglog(period[rxy],
                                   10**self.rplst[ii]['resxy'][0][rxy],
                                   ls=':', marker='s', ms=4, color='k', mfc='k')
                r2, = axrtm.loglog(period[ryx],
                                   10**self.rplst[ii]['resyx'][0][ryx],
                                   ls=':', marker='o', ms=4, color='k', mfc='k')
                rlstte = [r1]
                rlsttm = [r2]

            mrxy = [np.where(self.rplst[ii]['resxy'][2] != 0)[0]]
            mryx = [np.where(self.rplst[ii]['resyx'][2] != 0)[0]]
            r3, = axrte.loglog(period[mrxy], 10**self.rplst[ii]['resxy'][2][mrxy],
                               ls='-', color=colorlst[jj])
            r4, = axrtm.loglog(period[mryx], 10**self.rplst[ii]['resyx'][2][mryx],
                               ls='-', color=colorlst[jj])

            rlstte.append(r3)
            rlsttm.append(r4)

            # plot phase
            # cut out missing data points first
            pxy = [np.where(self.rplst[ii]['phasexy'][0] != 0)[0]]
            pyx = [np.where(self.rplst[ii]['phaseyx'][0] != 0)[0]]

            if jj == 0:
                axpte.semilogx(period[pxy], self.rplst[ii]['phasexy'][0][pxy],
                               ls=':', marker='s', ms=4, color='k', mfc='k')
                axptm.semilogx(period[pyx], self.rplst[ii]['phaseyx'][0][pyx],
                               ls=':', marker='o', ms=4, color='k', mfc='k')

            mpxy = [np.where(self.rplst[ii]['phasexy'][2] != 0)[0]]
            mpyx = [np.where(self.rplst[ii]['phaseyx'][2] != 0)[0]]
            axpte.semilogx(period[mpxy], self.rplst[ii]['phasexy'][2][mpxy],
                           ls='-', color=colorlst[jj])
            axptm.semilogx(period[mpyx], self.rplst[ii]['phaseyx'][2][mpyx],
                           ls='-', color=colorlst[jj])

        axrte.grid(True, alpha=.4)
        axrtm.grid(True, alpha=.4)

        axrtm.set_xticklabels(['' for ii in range(10)])
        axrte.set_xticklabels(['' for ii in range(10)])

        rmstestr = np.median(np.array(rmstestr)[1:])
        rmstmstr = np.median(np.array(rmstmstr)[1:])
        axrte.set_title('TE rms={0:.2f}'.format(rmstestr),
                        fontdict={'size': 10, 'weight': 'bold'})
        axrtm.set_title('TM rms={0:.2f}'.format(rmstmstr),
                        fontdict={'size': 10, 'weight': 'bold'})

        axpte.grid(True, alpha=.4)
        axpte.yaxis.set_major_locator(MultipleLocator(10))
        axpte.yaxis.set_minor_locator(MultipleLocator(1))

        axrte.set_ylabel('App. Res. ($\Omega \cdot m$)',
                         fontdict={'size': 10, 'weight': 'bold'})
        axpte.set_ylabel('Phase (deg)',
                         fontdict={'size': 10, 'weight': 'bold'})
        axpte.set_xlabel('Period (s)', fontdict={'size': 10, 'weight': 'bold'})

        axrte.yaxis.set_label_coords(-.08, .5)
        axpte.yaxis.set_label_coords(-.08, .5)

        axrtm.set_xticklabels(['' for ii in range(10)])
        axptm.grid(True, alpha=.4)
        axptm.yaxis.set_major_locator(MultipleLocator(10))
        axptm.yaxis.set_minor_locator(MultipleLocator(1))

        axrtm.set_ylabel('App. Res. ($\Omega \cdot m$)',
                         fontdict={'size': 12, 'weight': 'bold'})
        axptm.set_ylabel('Phase (deg)',
                         fontdict={'size': 12, 'weight': 'bold'})
        axptm.set_xlabel('Period (s)', fontdict={'size': 12, 'weight': 'bold'})

        axrtm.yaxis.set_label_coords(-.08, .5)
        axptm.yaxis.set_label_coords(-.08, .5)
        plt.suptitle(station, fontsize=12, fontweight='bold')
        plt.show()


class Occam2DModel(Occam2DData):
    """
    This class deals with the model side of Occam inversions, including 
    plotting the model, the L-curve, depth profiles.  It will also be able to 
    build a mesh and regularization grid at some point.  

    It inherits Occam2DData and the data can be extracted from the method
    get2DData().  After this call you can use all the methods of Occam2DData,
    such as plotting the model responses and pseudo sections.


    """

    def __init__(self, iterfn, meshfn=None, inmodelfn=None):
        self.iterfn = iterfn

        self.invpath = os.path.dirname(self.iterfn)

        # get meshfile if none is provides assuming the mesh file is named
        # with mesh
        if self.invpath != None:
            self.meshfn = os.path.join(self.invpath, 'MESH')
            if os.path.isfile(self.meshfn) == False:
                for ff in os.listdir(self.invpath):
                    if ff.lower().find('mesh') >= 0:
                        self.meshfn = os.path.join(self.invpath, ff)
                if os.path.isfile(self.meshfn) == False:
                    raise NameError('Could not find a mesh file, ' +
                                    'input manually')

        # get inmodelfile if none is provides assuming the mesh file is
        # named with inmodel
        if inmodelfn == None:
            self.inmodelfn = os.path.join(self.invpath, 'INMODEL')
            if os.path.isfile(self.inmodelfn) == False:
                for ff in os.listdir(self.invpath):
                    if ff.lower().find('inmodel') >= 0:
                        self.inmodelfn = os.path.join(self.invpath, ff)
                if os.path.isfile(self.inmodelfn) == False:
                    raise NameError('Could not find a model file, ' +
                                    'input manually')

    def read2DIter(self):
        """
        read2DIter will read an iteration file and combine that info from the 
        datafn and return a dictionary of variables.

        Arguments:
        ----------
            **iterfn** : string
                        full path to iteration file if iterpath=None.  If 
                        iterpath is input then iterfn is just the name
                        of the file without the full path.

        Returns:
        --------
            **Occam2DModel.idict** : dictionary of parameters, 
                                     keys are verbatim from the file, 
                                     except for the key 'model' which is the 

        :Example: ::

            >>> import mtpy.modeling.occamtools as occam
            >>> itfn = r"/home/Occam2D/Line1/Inv1/Test_15.iter"
            >>> ocm = occam.Occam2DModel(itfn)
            >>> ocm.read2DIter()

        """

        # check to see if the file exists
        if os.path.exists(self.iterfn) == False:
            raise IOError('File: ' + self.iterfn +
                          ' does not exist, check path')

        # open file, read lines, close file
        ifid = file(self.iterfn, 'r')
        ilines = ifid.readlines()
        ifid.close()

        # create dictionary to put things
        self.idict = {}
        ii = 0
        # put header info into dictionary with similar keys
        while ilines[ii].lower().find('param') != 0:
            iline = ilines[ii].strip().split(':')
            self.idict[iline[0].lower()] = iline[1].strip()
            ii += 1

        # get number of parameters
        iline = ilines[ii].strip().split(':')
        nparam = int(iline[1].strip())
        self.idict[iline[0]] = nparam
        self.idict['model'] = np.zeros(nparam)
        kk = int(ii + 1)

        jj = 0
        while jj < len(ilines) - kk:
            iline = ilines[jj + kk].strip().split()
            for ll in range(4):
                try:
                    self.idict['model'][jj * 4 + ll] = float(iline[ll])
                except IndexError:
                    pass
            jj += 1

        # get the data file name from the iteration header
        self.datafn = self.idict['data file']
        if self.datafn.find(os.sep) == -1:
            self.datafn = os.path.join(self.invpath, self.datafn)
        if os.path.isfile(self.datafn) == False:
            for ff in os.listdir(self.invpath):
                if ff.lower().find('.dat') >= 0:
                    self.datafn = os.path.join(self.invpath, ff)
            if os.path.isfile(self.datafn) == False:
                raise NameError('Could not find a data file, input manually')

    def read2DInmodel(self):
        """
        read an INMODEL file for occam 2D

        Arguments:
        ----------
            **inmodelfn** : string
                            full path to INMODEL file


        Returns:
        --------
            **Occam2DModel.rows** : list of combined data blocks where first 
                                    number of each list represents the number 
                                    of combined mesh layers for this 
                                    regularization block.  The second number is
                                    the number of columns in the regularization
                                    block layer.

            **Occam2DModel.cols** : list of combined mesh columns for the 
                                    regularization layer. The sum of this list 
                                    must be equal to the number of mesh columns

            **Occam2DModel.headerdict** : dictionary of all the header 
                                          information including the binding 
                                          offset.

        :Example: ::

            >>> import mtpy.modeling.occamtools as occam
            >>> itfn = r"/home/Occam2D/Line1/Inv1/Test_15.iter"
            >>> ocm = occam.Occam2DModel(itfn)
            >>> ocm.read2DInmodel()
        """

        ifid = open(self.inmodelfn, 'r')

        headerdict = {}
        rows = []
        cols = []
        ncols = []

        ilines = ifid.readlines()

        for ii, iline in enumerate(ilines):
            if iline.find(':') > 0:
                iline = iline.strip().split(':')
                headerdict[iline[0].lower()] = iline[1]
                # append the last line
                if iline[0].lower().find('exception') > 0:
                    cols.append(ncols)
            else:
                iline = iline.strip().split()
                iline = [int(jj) for jj in iline]
                if len(iline) == 2:
                    if len(ncols) > 0:
                        cols.append(ncols)
                    rows.append(iline)
                    ncols = []
                elif len(iline) > 2:
                    ncols = ncols + iline

        self.rows = np.array(rows)
        self.cols = cols
        self.inmodel_headerdict = headerdict

    def read2DMesh(self):
        """
        reads an Occam 2D mesh file

        Arguments:
        ----------
            **Occam2DModel.meshfn** : string 
                                      full path to mesh file

        Returns:
        --------
            **Occam2DModel.hnodes**: array of horizontal nodes 
                                    (column locations (m))

            **Occam2DModel.vnodes** : array of vertical nodes 
                                      (row locations(m))

            **Occam2DModel.mdata** : np.array of free parameters

        To do:
        ------
            incorporate fixed values

        :Example: ::

            >>> import mtpy.modeling.occamtools as occam
            >>> itfn = r"/home/Occam2D/Line1/Inv1/Test_15.iter"
            >>> ocm = occam.Occam2DModel(itfn)
            >>> ocm.read2DMesh()
        """

        mfid = file(self.meshfn, 'r')

        mlines = mfid.readlines()

        nh = int(mlines[1].strip().split()[1]) - 1
        nv = int(mlines[1].strip().split()[2]) - 1

        hnodes = np.zeros(nh + 1)
        vnodes = np.zeros(nv + 1)
        mdata = np.zeros((nh, nv, 4), dtype=str)

        # get horizontal nodes
        jj = 2
        ii = 0
        while ii < nh:
            hline = mlines[jj].strip().split()
            for mm in hline:
                hnodes[ii] = float(mm)
                print(ii, hnodes)
                ii += 1
            jj += 1

        # get vertical nodes
        ii = 0
        while ii < nv:
            vline = mlines[jj].strip().split()
            for mm in vline:
                vnodes[ii] = float(mm)
                ii += 1
            jj += 1

        # get free parameters
        for ii, mm in enumerate(mlines[jj + 1:]):
            kk = 0
            while kk < 4:
                mline = mm.rstrip()
                if mline.lower().find('exception') > 0:
                    break
                for jj in range(nh):
                    try:
                        mdata[jj, ii, kk] = mline[jj]
                    except IndexError:
                        pass
                kk += 1

        # make the node information an attributes of the occamModel class
        self.hnodes = hnodes
        self.vnodes = vnodes
        self.meshdata = mdata

    def get2DData(self):
        """
        get data from data file using the inherited :func:'read2DdataFile' 
        """
        try:
            self.read2DdataFile()
        except AttributeError:
            print('No Data file defined')

    def get2DModel(self):
        """
        get2DModel will create an array based on the FE mesh and fill the 
        values found from the regularization grid.  This way the array can 
        be manipulated as a 2D object and plotted as an image or a mesh.

        Returns:
        --------

            **Occam2DModel.resmodel** : model array with log resistivity values

            **Occam2DModel.plotx** : np.array
                                    horizontal distance of FE mesh (m) blocks

            **Occam2DModel.ploty** : np.array
                                    depth of vertical nodes of FE mesh (m)
        """

        # read iteration file to get model and data file
        self.read2DIter()

        # read in data file as an OccamData type
        print('Reading data from: ', self.datafn)
        self.get2DData()

        # read in MESH file
        print('Reading mesh from: ', self.meshfn)
        self.read2DMesh()

        #read in INMODEL
        print('Reading model from: ', self.inmodelfn)
        self.read2DInmodel()
        # get the binding offset which is the right side of the furthest left
        # block, this helps locate the model in relative space
        bndgoff = float(self.inmodel_headerdict['binding offset'])

        # make sure that the number of rows and number of columns are the same
        assert len(self.rows) == len(self.cols)

        # initiate the resistivity model to the shape of the FE mesh
        resmodel = np.zeros((self.vnodes.shape[0], self.hnodes.shape[0]))

        # read in the model and set the regularization block values to map onto
        # the FE mesh so that the model can be plotted as an image or regular
        # mesh.
        mm = 0
        for ii in range(len(self.rows)):
            # get the number of layers to combine
            # this index will be the first index in the vertical direction
            ny1 = self.rows[:ii, 0].sum()
            # the second index  in the vertical direction
            ny2 = ny1 + self.rows[ii][0]
            # make the list of amalgamated columns an array for ease
            lc = np.array(self.cols[ii])
            # loop over the number of amalgamated blocks
            for jj in range(len(self.cols[ii])):
                # get first in index in the horizontal direction
                nx1 = lc[:jj].sum()
                # get second index in horizontal direction
                nx2 = nx1 + lc[jj]
                # put the apporpriate resistivity value into all the amalgamated
                # model blocks of the regularization grid into the forward model
                # grid
                resmodel[ny1:ny2, nx1:nx2] = self.idict['model'][mm]
                mm += 1

        # make some arrays for plotting the model
        plotx = np.array([self.hnodes[:ii + 1].sum()
                          for ii in range(len(self.hnodes))])
        ploty = np.array([self.vnodes[:ii + 1].sum()
                          for ii in range(len(self.vnodes))])

        # center the grid onto the station coordinates
        x0 = bndgoff - plotx[self.cols[0][0] - 1]
        plotx = plotx + x0

        # flip the arrays around for plotting purposes
        # plotx=plotx[::-1] and make the first layer start at zero
        ploty = ploty[::-1] - ploty[0]

        # make a mesh grid to plot in the model coordinates
        self.meshx, self.meshy = np.meshgrid(plotx, ploty)

        # flip the resmodel upside down so that the top is the stations
        resmodel = np.flipud(resmodel)

        # make attributes of the class
        self.resmodel = resmodel
        self.plotx = plotx
        self.ploty = ploty

        # set the offsets of the stations and station list.
        self.offsetlst = []
        for rpdict in self.rplst:
            self.offsetlst.append(rpdict['offset'])

    def plot2DModel(self, datafn=None,
                    xpad=1.0, ypad=1.0, spad=1.0, ms=10, stationid=None,
                    fdict={'size': 8, 'rotation': 60, 'weight': 'normal'},
                    dpi=300, ylimits=None, xminorticks=5, yminorticks=1,
                    climits=(0, 4), cmap='jet_r', fs=8, femesh='off',
                    regmesh='off', aspect='auto', title='on', meshnum='off',
                    blocknum='off', blkfdict={'size': 3}, fignum=1,
                    plotdimensions=(10, 10), grid='off', yscale='km',
                    xlimits=None):
        """
        plotModel will plot the model output by occam in the iteration file.

        Arguments:
        ----------

            **datafn** : string 
                        full path to data file.  If none is input it will use 
                        the data file found in the iteration file.

            **xpad** : float (units of km) 
                       padding in the horizontal direction of model.

            **ypad** : float (units according to **yscale**)
                       padding in the vertical direction of the top of the 
                       model to fit the station names and markers.

            **spad** : float (units according to **yscale**)
                       padding of station names away from the top of the model,
                       this is kind of awkward at the moment especially if you
                       zoom into the model, it usually looks retarded and 
                       doesn't fit

            **ms** : float  
                     marker size in ambiguous points

            **stationid** : tuple (min,max)
                           index of station names to plot -> ex. pb01sdr would 
                           be stationid=(0,4) to plot pb01

            **fdict** : font dictionary for the station names, can have keys:

                    *'size'* : font size

                    *'rotation'* : angle of rotation (deg) of font

                    *'weight'* : weight of font 

                    *'color'* : color of font

                    *'style'* : style of font ex. 'italics'

            **plotdimensions** : tuple (x,y)
                               x-y dimensions of the figure in inches

            **dpi** : int 
                      dot per inch of figure, should be 300 for publications

            **ylimits** : tuple (min,max)
                          limits of depth scale (km). ex, ylimits=(0,30)

            **xlimits** : tuple (min,max)
                          limits of horizontal scale (km). ex, ylimits=(0,30)


            **xminorticks** : int or float
                              location of minor tick marks for the horizontal 
                              axis

            **yminorticks** : int or float
                              location of minor tick marks for vertical axis

            **climits** : tuple (min,max)
                          limits of log10(resistivity). ex. climits=(0,4)

            **cmap** : string
                       color map to plot the model image
                       see matplotlib.cm for all options

            **fs** : float
                     font size of axis labels

            **femesh** : string ('on','off')
                        'on' to plot finite element forward modeling mesh 
                        (black)

            **regmesh** : string ('on','off')
                         'on' to plot regularization mesh (blue)

            **aspect** : tuple (width,height)
                        aspect ratio of the figure, depends on your line 
                        length and the depth you want to investigate

            **title** : string ('on,'off',input,None)

                        * 'on' to put the RMS and Roughness as the title
                        * input a string that will be added to the RMS and
                          roughness put 
                        * None to not put a title on the plot and print out RMS
                          and roughness

            **meshnum** : string ('on','off')
                         'on' to plot FE mesh block numbers

            **fignum** : int
                        figure number to plot to

            **blocknum** : tuple ('on','off')
                          'on' to plot numbers on the regularization blocks

            **blkfdict** : font dictionary for the numbering of regularization 
                           blocks with keys:

                               *'size'* : float font size

                               *'weight'* : font weight

            **grid** : string ('major','minor','both')
                        * 'major' for major ticks grid
                        * 'minor' for a grid of the minor ticks
                        * 'both' for a grid with major and minor ticks

            **yscale** : string ('km','m')
                        * 'km' for depth in km 
                        * 'm' for depth in meters

        :Example: ::

            >>> import mtpy.modeling.occamtools as occam
            >>> itfn = r"/home/Occam2D/Line1/Inv1/Test_15.iter"
            >>> ocm = occam.Occam2DModel(itfn)
            >>> ocm.plot2DModel(ms=20,ylimits=(0,.350),yscale='m',spad=.10,
            >>>                 ypad=.125,xpad=.025,climits=(0,2.5),
            >>>                 aspect='equal')
        """

        # set the scale of the plot
        if yscale == 'km':
            dfactor = 1000.
            pfactor = 1.0
        elif yscale == 'm':
            dfactor = 1.
            pfactor = 1000.
        else:
            dfactor = 1000.
            pfactor = 1.0

        # get the model
        self.get2DModel()

        # set some figure properties to use the maiximum space
        plt.rcParams['font.size'] = int(dpi / 40.)
        plt.rcParams['figure.subplot.left'] = .08
        plt.rcParams['figure.subplot.right'] = .99
        plt.rcParams['figure.subplot.bottom'] = .1
        plt.rcParams['figure.subplot.top'] = .92
        plt.rcParams['figure.subplot.wspace'] = .01
#        plt.rcParams['text.usetex']=True

        # plot the model as a mesh
        fig = plt.figure(fignum, plotdimensions, dpi=dpi)
        plt.clf()

        # add a subplot to the figure with the specified aspect ratio
        ax = fig.add_subplot(1, 1, 1, aspect=aspect)

        # plot the model as a pcolormesh so the extents are constrained to
        # the model coordinates
        self.mesh_plot = ax.pcolormesh(self.meshx / dfactor, self.meshy / dfactor,
                                       self.resmodel, cmap=cmap, vmin=climits[
                                           0],
                                       vmax=climits[1])

        # make a colorbar for the resistivity
        cbx = mcb.make_axes(ax, shrink=.8, pad=.01)
        cb = mcb.ColorbarBase(cbx[0], cmap=cmap, norm=Normalize(vmin=climits[0],
                                                                vmax=climits[1]))
        cb.set_label('Resistivity ($\Omega \cdot$m)',
                     fontdict={'size': fs, 'weight': 'bold'})
        cb.set_ticks(np.arange(int(climits[0]), int(climits[1]) + 1))
        cb.set_ticklabels(['10$^{0}$'.format(nn) for nn in
                           np.arange(int(climits[0]), int(climits[1]) + 1)])

        # set the offsets of the stations and plot the stations
        # need to figure out a way to set the marker at the surface in all
        # views.
        for rpdict in self.rplst:
            # plot the station marker
            # plots a V for the station cause when you use scatter the spacing
            # is variable if you change the limits of the y axis, this way it
            # always plots at the surface.
            ax.text(rpdict['offset'] / dfactor, self.ploty.min(), 'V',
                    horizontalalignment='center',
                    verticalalignment='baseline',
                    fontdict={'size': ms, 'weight': 'bold', 'color': 'black'})

            # put station id onto station marker
            # if there is a station id index
            if stationid != None:
                ax.text(rpdict['offset'] / dfactor, -spad * pfactor,
                        rpdict['station'][stationid[0]:stationid[1]],
                        horizontalalignment='center',
                        verticalalignment='baseline',
                        fontdict=fdict)
            # otherwise put on the full station name found form data file
            else:
                ax.text(rpdict['offset'] / dfactor, -spad * pfactor,
                        rpdict['station'],
                        horizontalalignment='center',
                        verticalalignment='baseline',
                        fontdict=fdict)

        # set the initial limits of the plot to be square about the profile
        # line
        if ylimits == None:
            ax.set_ylim(abs(max(self.offsetlst) - min(self.offsetlst)) / dfactor,
                        -ypad * pfactor)
        else:
            ax.set_ylim(ylimits[1] * pfactor, (ylimits[0] - ypad) * pfactor)
        ax.set_xlim(min(self.offsetlst) / dfactor - (xpad * pfactor),
                    (max(self.offsetlst) / dfactor + (xpad * pfactor)))
        # set the axis properties
        ax.xaxis.set_minor_locator(MultipleLocator(xminorticks * pfactor))
        ax.yaxis.set_minor_locator(MultipleLocator(yminorticks * pfactor))
        if yscale == 'km':
            ax.set_xlabel('Horizontal Distance (km)',
                          fontdict={'size': fs, 'weight': 'bold'})
            ax.set_ylabel('Depth (km)', fontdict={
                          'size': fs, 'weight': 'bold'})
        elif yscale == 'm':
            ax.set_xlabel('Horizontal Distance (m)',
                          fontdict={'size': fs, 'weight': 'bold'})
            ax.set_ylabel('Depth (m)', fontdict={'size': fs, 'weight': 'bold'})

        # put a grid on if one is desired
        if grid == 'major':
            ax.grid(alpha=.3, which='major')
        if grid == 'minor':
            ax.grid(alpha=.3, which='minor')
        if grid == 'both':
            ax.grid(alpha=.3, which='both')
        else:
            pass

        # set title as rms and roughness
        if type(title) is str:
            if title == 'on':
                titlestr = os.path.join(os.path.basename(os.path.dirname(self.iterfn)),
                                        os.path.basename(self.iterfn))
                ax.set_title(titlestr +
                             ': RMS {0:.2f}, Roughness={1:.0f}'.format(
                                 float(self.idict['misfit value']),
                                 float(self.idict['roughness value'])),
                             fontdict={'size': fs + 1, 'weight': 'bold'})
            else:
                ax.set_title(title + '; RMS {0:.2f}, Roughness={1:.0f}'.format(
                    float(self.idict['misfit value']),
                    float(self.idict['roughness value'])),
                    fontdict={'size': fs + 1, 'weight': 'bold'})
        else:
            print('RMS {0:.2f}, Roughness={1:.0f}'.format(
                float(self.idict['misfit value']),
                float(self.idict['roughness value'])))

        # plot forward model mesh
        if femesh == 'on':
            for xx in self.plotx / dfactor:
                ax.plot([xx, xx], [0, self.ploty[0] / dfactor],
                        color='k', lw=.5)
            for yy in self.ploty / dfactor:
                ax.plot([self.plotx[0] / dfactor, self.plotx[-1] / dfactor],
                        [yy, yy], color='k', lw=.5)

        # plot the regularization mesh
        if regmesh == 'on':
            linelst = []
            for ii in range(len(self.rows)):
                # get the number of layers to combine
                # this index will be the first index in the vertical direction
                ny1 = self.rows[:ii, 0].sum()
                # the second index  in the vertical direction
                ny2 = ny1 + self.rows[ii][0]
                # make the list of amalgamated columns an array for ease
                lc = np.array(self.cols[ii])
                yline = ax.plot([self.plotx[0] / dfactor, self.plotx[-1] / dfactor],
                                [self.ploty[-ny1] / dfactor,
                                 self.ploty[-ny1] / dfactor],
                                color='b', lw=.5)
                linelst.append(yline)
                # loop over the number of amalgamated blocks
                for jj in range(len(self.cols[ii])):
                    # get first in index in the horizontal direction
                    nx1 = lc[:jj].sum()
                    # get second index in horizontal direction
                    nx2 = nx1 + lc[jj]
                    try:
                        if ny1 == 0:
                            ny1 = 1
                        xline = ax.plot([self.plotx[nx1] / dfactor,
                                         self.plotx[nx1] / dfactor],
                                        [self.ploty[-ny1] / dfactor,
                                         self.ploty[-ny2] / dfactor],
                                        color='b', lw=.5)
                        linelst.append(xline)
                    except IndexError:
                        pass

        # plot the mesh block numbers
        if meshnum == 'on':
            kk = 1
            for yy in self.ploty[::-1] / dfactor:
                for xx in self.plotx / dfactor:
                    ax.text(xx, yy, '{0}'.format(kk), fontdict={'size': 3})
                    kk += 1

        # plot regularization block numbers
        if blocknum == 'on':
            kk = 1
            for ii in range(len(self.rows)):
                # get the number of layers to combine
                # this index will be the first index in the vertical direction
                ny1 = self.rows[:ii, 0].sum()
                # the second index  in the vertical direction
                ny2 = ny1 + self.rows[ii][0]
                # make the list of amalgamated columns an array for ease
                lc = np.array(self.cols[ii])
                # loop over the number of amalgamated blocks
                for jj in range(len(self.cols[ii])):
                    # get first in index in the horizontal direction
                    nx1 = lc[:jj].sum()
                    # get second index in horizontal direction
                    nx2 = nx1 + lc[jj]
                    try:
                        if ny1 == 0:
                            ny1 = 1
                        # get center points of the blocks
                        yy = self.ploty[-ny1] - (self.ploty[-ny1] -
                                                 self.ploty[-ny2]) / 2
                        xx = self.plotx[nx1] - \
                            (self.plotx[nx1] - self.plotx[nx2]) / 2
                        # put the number
                        ax.text(xx / dfactor, yy / dfactor, '{0}'.format(kk),
                                fontdict=blkfdict,
                                horizontalalignment='center',
                                verticalalignment='center')
                    except IndexError:
                        pass
                    kk += 1
        self.model_axes = ax
        plt.show()

    def plotL2Curve(self, fnstem=None, fignum=1, dpi=300):
        """
        PlotL2Curve will plot the RMS vs iteration number for the given 
        inversion folder and roughness vs iteration number

        Arguments:
        ----------
            **fnstem** : string
                         filename stem to look for in case multiple inversions 
                         were run in the same folder.  If none then searches 
                         for anything ending in .iter

            **fignum** : int
                         figure number

            dpi** : int
                    dots per inch resolution of the figure


        :Example: ::

            >>> import mtpy.modeling.occamtools as occam
            >>> itfn = r"/home/Occam2D/Line1/Inv1/Test_15.iter"
            >>> ocm = occam.Occam2DModel(itfn)
            >>> ocm.plotL2Curve(fignum=2)
        """

        invpath = os.path.dirname(self.iterfn)

        if fnstem == None:
            iterlst = [os.path.join(invpath, itfile)
                       for itfile in os.listdir(invpath) if itfile.find('.iter') > 0]
        else:
            iterlst = [os.path.join(invpath, itfile)
                       for itfile in os.listdir(invpath) if itfile.find('.iter') > 0 and
                       itfile.find(fnstem) > 0]

        nr = len(iterlst)

        rmsarr = np.zeros((nr, 2))

        for itfile in iterlst:
            self.iterfn = itfile
            self.read2DIter()
            ii = int(self.idict['iteration'])
            rmsarr[ii, 0] = float(self.idict['misfit value'])
            rmsarr[ii, 1] = float(self.idict['roughness value'])

        # set the dimesions of the figure
        plt.rcParams['font.size'] = int(dpi / 40.)
        plt.rcParams['figure.subplot.left'] = .08
        plt.rcParams['figure.subplot.right'] = .90
        plt.rcParams['figure.subplot.bottom'] = .1
        plt.rcParams['figure.subplot.top'] = .90
        plt.rcParams['figure.subplot.wspace'] = .01

        # make figure instance
        fig = plt.figure(fignum, [6, 5], dpi=dpi)
        plt.clf()

        # make a subplot for RMS vs Iteration
        ax1 = fig.add_subplot(1, 1, 1)

        # plot the rms vs iteration
        l1, = ax1.plot(np.arange(1, nr, 1), rmsarr[
                       1:, 0], '-k', lw=1, marker='d', ms=5)

        # plot the median of the RMS
        m1, = ax1.plot(np.arange(0, nr, 1), np.repeat(np.median(rmsarr[1:, 0]), nr),
                       '--r', lw=.75)

        # plot the mean of the RMS
        m2, = ax1.plot(np.arange(0, nr, 1), np.repeat(np.mean(rmsarr[1:, 0]), nr),
                       ls='--', color='orange', lw=.75)

        # make subplot for RMS vs Roughness Plot
        ax2 = ax1.twiny()

        # plot the rms vs roughness
        l2, = ax2.plot(rmsarr[1:, 1], rmsarr[1:, 0], '--b', lw=.75, marker='o', ms=7,
                       mfc='white')
        for ii, rms in enumerate(rmsarr[1:, 0], 1):
            ax2.text(rmsarr[ii, 1], rms, '{0}'.format(ii),
                     horizontalalignment='center',
                     verticalalignment='center',
                     fontdict={'size': 6, 'weight': 'bold', 'color': 'blue'})

        # make a legend
        ax1.legend([l1, l2, m1, m2], ['RMS', 'Roughness',
                                      'Median_RMS={0:.2f}'.format(
                                          np.median(rmsarr[1:, 0])),
                                      'Mean_RMS={0:.2f}'.format(np.mean(rmsarr[1:, 0]))],
                   ncol=4, loc='upper center', columnspacing=.25, markerscale=.75,
                   handletextpad=.15)

        # set the axis properties for RMS vs iteration
        ax1.yaxis.set_minor_locator(MultipleLocator(.1))
        ax1.xaxis.set_minor_locator(MultipleLocator(1))
        ax1.set_ylabel('RMS', fontdict={'size': 8, 'weight': 'bold'})
        ax1.set_xlabel('Iteration', fontdict={'size': 8, 'weight': 'bold'})
        ax1.grid(alpha=.25, which='both')
        ax2.set_xlabel('Roughness', fontdict={'size': 8, 'weight': 'bold',
                                              'color': 'blue'})
        for t2 in ax2.get_xticklabels():
            t2.set_color('blue')

        plt.show()

    def plotDepthModel(self, dpi=300, depthmm=(1, 10000), plottype='1',
                       yscale='log', plotdimensions=(3, 6), plotnum=1, fignum=1):
        """
        Plots a depth section profile for a given set of stations.

        Arguments:
        ----------

            **dpi** : int
                      dots-per-inch resolution of figure

            **depthmm** : tuple (min,max)
                          minimum and maximum depth to plot in meters

            **plottype** : input as:
                            * '1' to plot all stations found
                            * [station list] to plot only a few stations

            **plotnum** : input as:
                          * 1 to plot in different figures
                          * 'all' to plot in all into one figure.

            **yscale** : 'log' for logarithmic or 'linear' for linear

            **plotdimensions** : tuple (width,height)
                                 figure dimensions in inches

            **fignum** : int
                         figure number

        :Example: ::

            >>> import mtpy.modeling.occamtools as occam
            >>> itfn = r"/home/Occam2D/Line1/Inv1/Test_15.iter"
            >>> ocm = occam.Occam2DModel(itfn)
            >>> #plot just a few stations depth profile in one figure
            >>> ocm.plotDepthModel(plottype=['MT01','MT05'],plotnum='all')


        """

        try:
            self.offsetlst
        except AttributeError:
            self.get2DModel()
        # get stations to plot
        if plottype == '1':
            pstationlst = np.arange(len(self.stationlst))
        else:
            pstationlst = []
            if type(plottype) is not list:
                plottype = [plottype]
            for ps in plottype:
                for ii, ss in enumerate(self.stationlst):
                    if ss.find(ps) == 0:
                        pstationlst.append(ii)

        # get the average x-spacing within the station region, occam pads by
        # 7 cells by default
        xavg = np.floor(np.mean([abs(self.plotx[ii] - self.plotx[ii + 1])
                                 for ii in range(7, len(self.plotx) - 7)]))

        # get the station indices to extract from the model
        slst = []
        for ff in pstationlst:
            offset = self.offsetlst[ff]
            for ii, xx in enumerate(self.plotx):
                if offset >= xx - xavg / 2. and offset <= xx + xavg / 2.:
                    slst.append(ii)

        # get depth limits
        if depthmm == None:
            depthmm = (self.ploty.min(), self.ploty.max())
        if depthmm[0] == 0:
            depthmm[0] = 1

        # set the dimesions of the figure
        plt.rcParams['font.size'] = int(dpi / 40.)
        plt.rcParams['figure.subplot.left'] = .15
        plt.rcParams['figure.subplot.right'] = .95
        plt.rcParams['figure.subplot.bottom'] = .15
        plt.rcParams['figure.subplot.top'] = .90
        plt.rcParams['figure.subplot.wspace'] = .05

        if plotnum == 'all':
            # set the dimesions of the figure
            plt.rcParams['font.size'] = int(dpi / 60.)
            plt.rcParams['figure.subplot.left'] = .09
            plt.rcParams['figure.subplot.right'] = .95
            plt.rcParams['figure.subplot.bottom'] = .15
            plt.rcParams['figure.subplot.top'] = .90
            plt.rcParams['figure.subplot.wspace'] = .1

            fig = plt.figure(fignum, plotdimensions, dpi=dpi)
            plt.clf()
            ns = len(slst)
            # plot the depth section for each station
            for ii, ss in enumerate(slst):
                ax = fig.add_subplot(1, ns, ii + 1)

                # plot resistivity vs depth
                if yscale == 'linear':
                    p1, = ax.semilogx(10**self.resmodel[:, ss], self.ploty,
                                      ls='steps-')
                elif yscale == 'log':
                    if self.ploty[-1] == 0.0:
                        self.ploty[-1] = 1
                    p1, = ax.loglog(10**self.resmodel[:, ss], self.ploty,
                                    ls='steps-')
                ax.set_ylim(depthmm[1], depthmm[0])

                ax.set_title(self.data.stationlst[pstationlst[ii]],
                             fontdict={'size': 10, 'weight': 'bold'})
                if ii == 0:
                    ax.set_ylabel('Depth (m)',
                                  fontdict={'size': 8, 'weight': 'bold'})
                else:
                    plt.setp(ax.yaxis.get_ticklabels(), visible=False)
                if ii == np.round(ns / 2.):
                    ax.set_xlabel('Resistivity ($\Omega \cdot$m)',
                                  fontdict={'size': 8, 'weight': 'bold'})
                ax.grid(True, alpha=.3, which='both')
                ax.set_xlim(10**self.resmodel.min(), 10**self.resmodel.max())
        else:
            # plot the depth section for each station
            for ii, ss in enumerate(slst):
                fig = plt.figure(ii + 1, plotdimensions, dpi=dpi)
                plt.clf()
                ax = fig.add_subplot(1, 1, 1)

                # plot resistivity vs depth
                if yscale == 'linear':
                    p1, = ax.semilogx(10**self.resmodel[:, ss], self.ploty,
                                      ls='steps-')
                elif yscale == 'log':
                    if self.ploty[-1] == 0.0:
                        self.ploty[-1] = 1
                    p1, = ax.loglog(10**self.resmodel[:, ss], self.ploty,
                                    ls='steps-')
                ax.set_ylim(depthmm[1], depthmm[0])

                ax.set_title(self.stationlst[pstationlst[ii]],
                             fontdict={'size': 10, 'weight': 'bold'})
                ax.set_ylabel('Depth (m)', fontdict={
                              'size': 8, 'weight': 'bold'})
                ax.set_xlabel('Resistivity ($\Omega \cdot$m)',
                              fontdict={'size': 8, 'weight': 'bold'})
                ax.grid(True, alpha=.3, which='both')
