# -*- coding: utf-8 -*-
"""
==================
Occam1D
==================

    * Wrapper class to interact with Occam1D written by Kerry Keys at Scripps 
      adapted from the method of Constable et al., [1987].

    * This class only deals with the MT functionality of the Fortran code, so
      it can make the input files for computing the 1D MT response of an input
      model and or data.  It can also read the output and plot them in a
      useful way.

    * Note that when you run the inversion code, the convergence is quite
      quick, within the first few iterations, so have a look at the L2 cure
      to decide which iteration to plot, otherwise if you look at iterations
      long after convergence the models will be unreliable.



    * Key, K., 2009, 1D inversion of multicomponent, multi-frequency marine
      CSEM data: Methodology and synthetic studies for resolving thin
      resistive layers: Geophysics, 74, F9–F20.

    * The original paper describing the Occam's inversion approach is:

    * Constable, S. C., R. L. Parker, and C. G. Constable, 1987,
      Occam’s inversion –– A practical algorithm for generating smooth
      models from electromagnetic sounding data, Geophysics, 52 (03), 289–300.


    :Intended Use: ::

        >>> import mtpy.modeling.occam1d as occam1d
        >>> #--> make a data file
        >>> d1 = occam1d.Data()
        >>> d1.write_data_file(edi_file=r'/home/MT/mt01.edi', res_err=10, phase_err=2.5,
        >>> ...                save_path=r"/home/occam1d/mt01/TE", mode='TE')
        >>> #--> make a model file
        >>> m1 = occam1d.Model()
        >>> m1.write_model_file(save_path=d1.save_path, target_depth=15000)
        >>> #--> make a startup file
        >>> s1 = occam1d.Startup()
        >>> s1.data_fn = d1.data_fn
        >>> s1.model_fn = m1.model_fn
        >>> s1.save_path = m1.save_path
        >>> s1.write_startup_file()
        >>> #--> run occam1d from python
        >>> occam_path = r"/home/occam1d/Occam1D_executable"
        >>> occam1d.Run(s1.startup_fn, occam_path, mode='TE')
        >>> #--plot the L2 curve
        >>> l2 = occam1d.PlotL2(d1.save_path, m1.model_fn)
        >>> #--> see that iteration 7 is the optimum model to plot
        >>> p1 = occam1d.Plot1DResponse()
        >>> p1.data_te_fn = d1.data_fn
        >>> p1.model_fn = m1.model_fn
        >>> p1.iter_te_fn = r"/home/occam1d/mt01/TE/TE_7.iter"
        >>> p1.resp_te_fn = r"/home/occam1d/mt01/TE/TE_7.resp"
        >>> p1.plot()

@author: J. Peacock (Oct. 2013)
"""
# ------------------------------------------------------------------------------
import numpy as np
import os
import os.path as op
import time
from matplotlib.ticker import MultipleLocator
import matplotlib.gridspec as gridspec
import mtpy.core.mt as mt
import mtpy.utils.calculator as mtcc
import mtpy.analysis.geometry as mtg
import matplotlib.pyplot as plt
import subprocess
import string


# ------------------------------------------------------------------------------

class Data(object):
    """
    reads and writes occam 1D data files

    ===================== =====================================================
    Attributes             Description
    ===================== =====================================================
    _data_fn              basename of data file *default* is Occam1DDataFile
    _header_line          header line for description of data columns
    _ss                   string spacing *default* is 6*' '
    _string_fmt           format of data *default* is '+.6e'
    data                  array of data
    data_fn               full path to data file
    freq                  frequency array of data
    mode                  mode to invert for [ 'TE' | 'TM' | 'det' ]
    phase_te              array of TE phase
    phase_tm              array of TM phase
    res_te                array of TE apparent resistivity
    res_tm                array of TM apparent resistivity
    resp_fn               full path to response file
    save_path             path to save files to
    ===================== =====================================================


    ===================== =====================================================
    Methods               Description
    ===================== =====================================================
    write_data_file       write an Occam1D data file
    read_data_file        read an Occam1D data file
    read_resp_file        read a .resp file output by Occam1D
    ===================== =====================================================

    :Example: ::

        >>> import mtpy.modeling.occam1d as occam1d
        >>> #--> make a data file for TE mode
        >>> d1 = occam1d.Data()
        >>> d1.write_data_file(edi_file=r'/home/MT/mt01.edi', res_err=10, phase_err=2.5,
        >>> ...                save_path=r"/home/occam1d/mt01/TE", mode='TE')

    """

    def __init__(self, data_fn=None, **kwargs):
        self.data_fn = data_fn

        if self.data_fn is not None:
            self.save_path = os.path.dirname(self.data_fn)
        else:
            self.save_path = os.getcwd()

        self._string_fmt = '+.6e'
        self._ss = 6 * ' '
        self._data_fn = 'Occam1d_DataFile'
        self._header_line = '!{0}\n'.format('      '.join(['Type', 'Freq#',
                                                           'TX#', 'Rx#', 'Data',
                                                           'Std_Error']))
        self.mode = 'det'
        self.data = None

        self.freq = None
        self.res_te = None
        self.res_tm = None
        self.phase_te = None
        self.phase_tm = None
        self.resp_fn = None

        for key in list(kwargs.keys()):
            setattr(self, key, kwargs[key])

    def write_data_file(self, rp_tuple=None, edi_file=None, save_path=None,
                        mode='det', res_err='data', phase_err='data', thetar=0,
                        res_errorfloor=0., phase_errorfloor=0., z_errorfloor=0.,
                        remove_outofquadrant=False):
        """
        make1Ddatafile will write a data file for Occam1D

        Arguments:
        ---------
            **rp_tuple** : np.ndarray (freq, res, res_err, phase, phase_err)
                            with res, phase having shape (num_freq, 2, 2).

            **edi_file** : string
                          full path to edi file to be modeled.

            **save_path** : string
                           path to save the file, if None set to dirname of
                           station if edipath = None.  Otherwise set to
                           dirname of edipath.

            **thetar** : float
                         rotation angle to rotate Z. Clockwise positive and N=0
                         *default* = 0

            **mode** : [ 'te' | 'tm' | 'det']
                              mode to model can be (*default*='both'):
                                - 'te' for just TE mode (res/phase)
                                - 'tm' for just TM mode (res/phase)
                                - 'det' for the determinant of Z (converted to
                                        res/phase)
                              add 'z' to any of these options to model
                              impedance tensor values instead of res/phase


            **res_err** : float
                        errorbar for resistivity values.  Can be set to (
                        *default* = 'data'):

                        - 'data' for errorbars from the data
                        - percent number ex. 10 for ten percent

            **phase_err** : float
                          errorbar for phase values.  Can be set to (
                          *default* = 'data'):

                            - 'data' for errorbars from the data
                            - percent number ex. 10 for ten percent
            **res_errorfloor**: float
                                error floor for resistivity values
                                in percent
            **phase_errorfloor**: float
                                  error floor for phase in degrees
            **remove_outofquadrant**: True/False; option to remove the resistivity and
                                      phase values for points with phases out
                                      of the 1st/3rd quadrant (occam requires
                                      0 < phase < 90 degrees; phases in the 3rd
                                      quadrant are shifted to the first by
                                      adding 180 degrees)

        :Example: ::

            >>> import mtpy.modeling.occam1d as occam1d
            >>> #--> make a data file
            >>> d1 = occam1d.Data()
            >>> d1.write_data_file(edi_file=r'/home/MT/mt01.edi', res_err=10,
            >>> ...                phase_err=2.5, mode='TE',
            >>> ...                save_path=r"/home/occam1d/mt01/TE")
        """
        # be sure that the input mode is not case sensitive
        self.mode = mode.lower()

        if self.mode == 'te':
            d1_str = 'RhoZxy'
            d2_str = 'PhsZxy'
        elif self.mode == 'tm':
            d1_str = 'RhoZyx'
            d2_str = 'PhsZyx'
        elif self.mode == 'det':
            d1_str = 'RhoZxy'
            d2_str = 'PhsZxy'
        elif self.mode == 'detz':
            d1_str = 'RealZxy'
            d2_str = 'ImagZxy'
        elif self.mode == 'tez':
            d1_str = 'RealZxy'
            d2_str = 'ImagZxy'
        elif self.mode == 'tmz':
            d1_str = 'RealZyx'
            d2_str = 'ImagZyx'

        # read in data as a tuple or .edi file
        if edi_file is None and rp_tuple is None:
            raise IOError('Need to input either an edi file or rp_array')

        if edi_file is not None:

            # read in edifile
            mt_obj = mt.MT(edi_file)
            z_obj = mt_obj.Z
            z_obj.compute_resistivity_phase()

            # get frequencies to invert
            freq = z_obj.freq
            nf = len(freq)

            # rotate if necessary
            if thetar != 0:
                z_obj.rotate(thetar)

            # get the data requested by the given mode
            if self.mode == 'te':
                data_1 = z_obj.resistivity[:, 0, 1]
                data_1_err = z_obj.resistivity_err[:, 0, 1]

                data_2 = z_obj.phase[:, 0, 1]
                data_2_err = z_obj.phase_err[:, 0, 1]

            elif self.mode == 'tm':
                data_1 = z_obj.resistivity[:, 1, 0]
                data_1_err = z_obj.resistivity_err[:, 1, 0]

                # need to put the angle in the right quadrant
                data_2 = z_obj.phase[:, 1, 0] % 180
                data_2_err = z_obj.phase_err[:, 1, 0]

            elif self.mode.startswith('det'):

                # take square root as determinant is similar to z squared
                zdetreal = (z_obj.det ** 0.5).real
                zdetimag = (z_obj.det ** 0.5).imag

                # error propagation - new error is 0.5 * relative error in zdet
                # then convert back to absolute error
                det_err = mtcc.compute_determinant_error(z_obj.z,z_obj.z_err)
                
#                # relative errors of real and imaginary components of sqrt determinant
                zereal = zdetreal * det_err * 0.5 / z_obj.det.real
                zeimag = zdetimag * det_err * 0.5 / z_obj.det.imag


                if self.mode.endswith('z'):
                    # convert to si units if we are modelling impedance tensor
                    
                    data_1 = zdetreal * np.pi * 4e-4
                    data_1_err = zereal * np.pi * 4e-4
                    data_2 = zdetimag * np.pi * 4e-4
                    data_2_err = zeimag * np.pi * 4e-4
                else:
                    # convert to res/phase
                    # data_1 is resistivity
                    # need to take absolute of square root before squaring it
                    data_1 = .2 / freq * np.abs(z_obj.det**0.5)**2
                    # data_2 is phase
                    data_2 = np.rad2deg(np.arctan2(zdetimag, zdetreal))

                    # initialise error arrays
                    data_1_err = np.zeros_like(data_1, dtype=np.float)
                    data_2_err = np.zeros_like(data_2, dtype=np.float)

                    # assign errors, use error based on sqrt of z
                    for zdr, zdi, zer, zei, ii in zip(zdetreal, zdetimag,
                                                      zereal, zeimag,
                                                      list(range(len(z_obj.det)))):
                        # now we can convert errors to polar coordinates
                        de1, de2 = mtcc.z_error2r_phi_error(zdr, zdi, (zei+zer)/2.)
                        # convert relative resistivity error to absolute
                        de1 *= data_1[ii]
                        data_1_err[ii] = de1
                        data_2_err[ii] = de2
                        
            elif self.mode == 'tez':
                # convert to si units
                data_1 = z_obj.z[:, 0, 1].real * np.pi * 4e-4
                data_1_err = z_obj.z_err[:, 0, 1] * np.pi * 4e-4

                data_2 = z_obj.z[:, 0, 1].imag * np.pi * 4e-4
                data_2_err = z_obj.z_err[:, 0, 1] * np.pi * 4e-4

            elif self.mode == 'tmz':
                # convert to si units
                data_1 = z_obj.z[:, 1, 0].real * np.pi * 4e-4
                data_1_err = z_obj.z_err[:, 1, 0] * np.pi * 4e-4

                data_2 = z_obj.z[:, 1, 0].imag * np.pi * 4e-4
                data_2_err = z_obj.z_err[:, 1, 0] * np.pi * 4e-4

            else:
                raise IOError('Mode {0} is not supported.'.format(self.mode))

        if rp_tuple is not None:
            if len(rp_tuple) != 5:
                raise IOError('Be sure rp_array is correctly formated\n'
                              'should be freq, res, res_err, phase, phase_err')
            freq, rho, rho_err, phi, phi_err = rp_tuple
            nf = len(freq)

            if self.mode in 'te':
                data_1 = rho[:, 0, 1]
                data_1_err = rho_err[:, 0, 1]

                data_2 = phi[:, 0, 1]
                data_2_err = phi_err[:, 0, 1]

            elif self.mode in 'tm':
                data_1 = rho[:, 1, 0]
                data_1_err = rho_err[:, 1, 0]

                data_2 = phi[:, 1, 0] % 180
                data_2_err = phi_err[:, 1, 0]

            if 'det' in mode.lower():
                data_1 = rho[:, 0, 1]
                data_1_err = rho_err[:, 0, 1]

                data_2 = phi[:, 0, 1]
                data_2_err = phi_err[:, 0, 1]

        if remove_outofquadrant:
            freq, data_1, data_1_err, data_2, data_2_err = self._remove_outofquadrant_phase(
                freq,
                data_1,
                data_1_err,
                data_2,
                data_2_err)
            nf = len(freq)

        # ---> get errors--------------------------------------
        # set error floors
        if 'z' in self.mode:
            if z_errorfloor > 0:
                data_1_err = np.abs(data_1_err)
                test = data_1_err / np.abs(data_1 + 1j * data_2) < z_errorfloor / 100.
                data_1_err[test] = np.abs(data_1 + 1j * data_2)[test] * z_errorfloor / 100.
                data_2_err = data_1_err.copy()
        else:
            if res_errorfloor > 0:
                test = data_1_err / data_1 < res_errorfloor / 100.
                data_1_err[test] = data_1[test] * res_errorfloor / 100.
            if phase_errorfloor > 0:
                data_2_err[data_2_err < phase_errorfloor] = phase_errorfloor

            if res_err != 'data':
                data_1_err = data_1 * res_err / 100.
            if phase_err != 'data':
                data_2_err = np.repeat(phase_err / 100. * (180 / np.pi), nf)

        # --> write file
        # make sure the savepath exists, if not create it
        if save_path is not None:
            self.save_path = save_path
        if self.save_path == None:
            try:
                self.save_path = os.path.dirname(edi_file)
            except TypeError:
                pass
        elif os.path.basename(self.save_path).find('.') > 0:
            self.save_path = os.path.dirname(self.save_path)
            self._data_fn = os.path.basename(self.save_path)
        if not os.path.exists(self.save_path):
            os.mkdir(self.save_path)

        if self.data_fn is None:
            self.data_fn = os.path.join(self.save_path,
                                    '{0}_{1}.dat'.format(self._data_fn, mode.upper()))

        # --> write file as a list of lines
        dlines = []

        dlines.append('Format:  EMData_1.1 \n')
        dlines.append('!mode:   {0}\n'.format(mode.upper()))
        dlines.append('!rotation_angle = {0:.2f}\n'.format(thetar))

        # needs a transmitter to work so put in a dummy one
        dlines.append('# Transmitters: 1\n')
        dlines.append('0 0 0 0 0 \n')

        # write frequencies
        dlines.append('# Frequencies:   {0}\n'.format(nf))
        if freq[0] < freq[-1]:
            freq = freq[::-1]
            data_1 = data_1[::-1]
            data_2 = data_2[::-1]
            data_1_err = data_1_err[::-1]
            data_2_err = data_2_err[::-1]
        for ff in freq:
            dlines.append('   {0:{1}}\n'.format(ff, self._string_fmt))

        # needs a receiver to work so put in a dummy one
        dlines.append('# Receivers: 1 \n')
        dlines.append('0 0 0 0 0 0 \n')

        # write data
        dlines.append('# Data:{0}{1}\n'.format(self._ss, 2 * nf))
        num_data_line = len(dlines)

        dlines.append(self._header_line)
        data_count = 0

        #        data1 = np.abs(data1)
        #        data2 = np.abs(data2)

        for ii in range(nf):
            # write lines
            if data_1[ii] != 0.0:
                dlines.append(self._ss.join([d1_str, str(ii + 1), '0', '1',
                                             '{0:{1}}'.format(data_1[ii], self._string_fmt),
                                             '{0:{1}}\n'.format(data_1_err[ii], self._string_fmt)]))
                data_count += 1
            if data_2[ii] != 0.0:
                dlines.append(self._ss.join([d2_str, str(ii + 1), '0', '1',
                                             '{0:{1}}'.format(data_2[ii], self._string_fmt),
                                             '{0:{1}}\n'.format(data_2_err[ii], self._string_fmt)]))
                data_count += 1

        # --> write file
        dlines[num_data_line - 1] = '# Data:{0}{1}\n'.format(self._ss, data_count)

        with open(self.data_fn, 'w') as dfid:
            dfid.writelines(dlines)

        print('Wrote Data File to : {0}'.format(self.data_fn))

        # --> set attributes

        if 'z' in mode.lower():
            self.z = data_1 + 1j * data_2
            self.z_err = data_1_err
        else:
            if 'det' in mode.lower():
                self.res_det = data_1
                self.phase_det = data_2
            elif self.mode == 'te':
                self.res_te = data_1
                self.phase_te = data_2
            elif self.mode == 'tm':
                self.res_tm = data_1
                self.phase_tm = data_2
            self.res_err = data_1_err
            self.phase_err = data_2_err

        self.freq = freq

    def _remove_outofquadrant_phase(self, freq, d1, d1_err, d2, d2_err):
        """
        remove out of quadrant phase from data
        """
        # remove data points with phase out of quadrant
        if 'z' in self.mode:
            include = (d1 / d2 > 0) & (d1 / d2 > 0)
        elif self.mode in ['det', 'te', 'tm']:
            include = (d2 % 180 <= 90) & (d2 % 180 >= 0) & (d2 % 180 <= 90) & (d2 % 180 >= 0)

        newfreq, nd1, nd1_err, nd2, nd2_err = [arr[include] for arr in [freq,
                                                                        d1, d1_err, d2, d2_err]]
        # fix any zero errors to 100% of the res value or 90 degrees for phase
        nd1_err[nd1_err == 0] = nd1[nd1_err == 0]
        if 'z' in self.mode:
            nd2_err[nd2_err == 0] = nd2[nd2_err == 0]
        else:
            nd2_err[nd2_err == 0] = 90

        return newfreq, nd1, nd1_err, nd2, nd2_err

    def read_data_file(self, data_fn=None):
        """
        reads a 1D data file

        Arguments:
        ----------
            **data_fn** : full path to data file

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

            >>> old = occam1d.Data()
            >>> old.data_fn = r"/home/Occam1D/Line1/Inv1_TE/MT01TE.dat"
            >>> old.read_data_file()
        """

        if data_fn is not None:
            self.data_fn = data_fn
        if self.data_fn is None:
            raise IOError('Need to input a data file')
        elif os.path.isfile(self.data_fn) == False:
            raise IOError('Could not find {0}, check path'.format(self.data_fn))

        self._data_fn = os.path.basename(self.data_fn)
        self.save_path = os.path.dirname(self.data_fn)

        dfid = open(self.data_fn, 'r')

        # read in lines
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
        nfreq = int(dlines[finddict['Frequencies']][2:].strip().split(':')[1].strip())

        # frequency list
        freq = np.array([float(ff) for ff in dlines[finddict['Frequencies'] + 1:
        finddict['Receivers']]])

        # data dictionary to put things into
        # check to see if there is alread one, if not make a new one
        if self.data is None:
            self.data = {'freq': freq,
                         'zxy': np.zeros((4, nfreq), dtype=complex),
                         'zyx': np.zeros((4, nfreq), dtype=complex),
                         'resxy': np.zeros((4, nfreq)),
                         'resyx': np.zeros((4, nfreq)),
                         'phasexy': np.zeros((4, nfreq)),
                         'phaseyx': np.zeros((4, nfreq))}

        # get data
        for dline in dlines[finddict['Data'] + 1:]:
            if dline.find('!') == 0:
                pass
            else:
                dlst = dline.strip().split()
                dlst = [dd.strip() for dd in dlst]
                if len(dlst) > 4:
                    jj = int(dlst[1]) - 1
                    dvalue = float(dlst[4])
                    derr = float(dlst[5])
                    if dlst[0] == 'RhoZxy' or dlst[0] == '103':
                        self.mode = 'TE'
                        self.data['resxy'][0, jj] = dvalue
                        self.data['resxy'][1, jj] = derr
                    if dlst[0] == 'PhsZxy' or dlst[0] == '104':
                        self.mode = 'TE'
                        self.data['phasexy'][0, jj] = dvalue
                        self.data['phasexy'][1, jj] = derr
                    if dlst[0] == 'RhoZyx' or dlst[0] == '105':
                        self.mode = 'TM'
                        self.data['resyx'][0, jj] = dvalue
                        self.data['resyx'][1, jj] = derr
                    if dlst[0] == 'PhsZyx' or dlst[0] == '106':
                        self.mode = 'TM'
                        self.data['phaseyx'][0, jj] = dvalue
                        self.data['phaseyx'][1, jj] = derr
                    if dlst[0] == 'RealZxy' or dlst[0] == '113':
                        self.mode = 'TEz'
                        self.data['zxy'][0, jj] = dvalue / (np.pi * 4e-4)
                        self.data['zxy'][1, jj] = derr / (np.pi * 4e-4)
                    if dlst[0] == 'ImagZxy' or dlst[0] == '114':
                        self.mode = 'TEz'
                        self.data['zxy'][0, jj] += 1j * dvalue / (np.pi * 4e-4)
                        self.data['zxy'][1, jj] = derr / (np.pi * 4e-4)
                    if dlst[0] == 'RealZyx' or dlst[0] == '115':
                        self.mode = 'TMz'
                        self.data['zyx'][0, jj] = dvalue / (np.pi * 4e-4)
                        self.data['zyx'][1, jj] = derr / (np.pi * 4e-4)
                    if dlst[0] == 'ImagZyx' or dlst[0] == '116':
                        self.mode = 'TMz'
                        self.data['zyx'][0, jj] += 1j * dvalue / (np.pi * 4e-4)
                        self.data['zyx'][1, jj] = derr / (np.pi * 4e-4)

        if 'z' in self.mode:
            if 'TE' in self.mode:
                pol = 'xy'
            elif 'TM' in self.mode:
                pol = 'yx'

            self.data['res' + pol][0] = 0.2 * np.abs(self.data['z' + pol][0]) ** 2. / freq
            self.data['phase' + pol][0] = np.rad2deg(
                np.arctan(self.data['res' + pol][0].imag / self.data['res' + pol][0].real))
            for jjj in range(len(freq)):
                res_rel_err, phase_err = \
                    mtcc.z_error2r_phi_error(self.data['z' + pol][0, jjj].real,
                                             self.data['z' + pol][0, jjj].imag, 
                                             self.data['z' + pol][1, jjj])
                    
                self.data['res' + pol][1, jjj], self.data['phase' + pol][1, jjj] = \
                    res_rel_err*self.data['res' + pol][0, jjj], phase_err

            self.data['resyx'][0] = 0.2 * np.abs(self.data['zxy'][0]) ** 2. / freq

        self.freq = freq
        self.res_te = self.data['resxy']
        self.res_tm = self.data['resyx']
        self.phase_te = self.data['phasexy']
        self.phase_tm = self.data['phaseyx']

    def read_resp_file(self, resp_fn=None, data_fn=None):
        """
        read response file

        Arguments:
        ---------
            **resp_fn** : full path to response file

            **data_fn** : full path to data file

        Fills:
        --------

            *freq* : an array of frequencies with length nf

            *res_te* : TE resistivity array with shape (nf,4) for (0) data,
                      (1) dataerr, (2) model, (3) modelerr

            *res_tm* : TM resistivity array with shape (nf,4) for (0) data,
                      (1) dataerr, (2) model, (3) modelerr

            *phase_te* : TE phase array with shape (nf,4) for (0) data,
                        (1) dataerr, (2) model, (3) modelerr

            *phase_tm* : TM phase array with shape (nf,4) for (0) data,
                        (1) dataerr, (2) model, (3) modelerr

       :Example: ::
            >>> o1d = occam1d.Data()
            >>> o1d.data_fn = r"/home/occam1d/mt01/TE/Occam1D_DataFile_TE.dat"
            >>> o1d.read_resp_file(r"/home/occam1d/mt01/TE/TE_7.resp")

        """

        if resp_fn is not None:
            self.resp_fn = resp_fn
        if self.resp_fn is None:
            raise IOError('Need to input response file')

        if data_fn is not None:
            self.data_fn = data_fn
        if self.data_fn is None:
            raise IOError('Need to input data file')
        # --> read in data file
        self.read_data_file()

        # --> read response file
        dfid = open(self.resp_fn, 'r')

        dlines = dfid.readlines()
        dfid.close()

        finddict = {}
        for ii, dline in enumerate(dlines):
            if dline.find('#') <= 3:
                fkey = dline[2:].strip().split(':')[0]
                fvalue = ii
                finddict[fkey] = fvalue

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
                    try:
                        rerr = float(dlst[7])
                    except ValueError:
                        rerr = 1000.
                    if dlst[0] == 'RhoZxy' or dlst[0] == '103':
                        self.res_te[0, jj] = dvalue
                        self.res_te[1, jj] = derr
                        self.res_te[2, jj] = rvalue
                        self.res_te[3, jj] = rerr
                    if dlst[0] == 'PhsZxy' or dlst[0] == '104':
                        self.phase_te[0, jj] = dvalue
                        self.phase_te[1, jj] = derr
                        self.phase_te[2, jj] = rvalue
                        self.phase_te[3, jj] = rerr
                    if dlst[0] == 'RhoZyx' or dlst[0] == '105':
                        self.res_tm[0, jj] = dvalue
                        self.res_tm[1, jj] = derr
                        self.res_tm[2, jj] = rvalue
                        self.res_tm[3, jj] = rerr
                    if dlst[0] == 'PhsZyx' or dlst[0] == '106':
                        self.phase_tm[0, jj] = dvalue
                        self.phase_tm[1, jj] = derr
                        self.phase_tm[2, jj] = rvalue
                        self.phase_tm[3, jj] = rerr
                    if dlst[0] == 'RealZxy' or dlst[0] == '113':
                        self.mode = 'TEz'
                        self.data['zxy'][0, jj] = dvalue / (np.pi * 4e-4)
                        self.data['zxy'][1, jj] = derr / (np.pi * 4e-4)
                        self.data['zxy'][2, jj] = rvalue / (np.pi * 4e-4)
                        self.data['zxy'][3, jj] = rerr
                    if dlst[0] == 'ImagZxy' or dlst[0] == '114':
                        self.mode = 'TEz'
                        self.data['zxy'][0, jj] += 1j * dvalue / (np.pi * 4e-4)
                        self.data['zxy'][1, jj] = derr / (np.pi * 4e-4)
                        self.data['zxy'][2, jj] += 1j * rvalue / (np.pi * 4e-4)
                        self.data['zxy'][3, jj] = rerr
                    if dlst[0] == 'RealZyx' or dlst[0] == '115':
                        self.mode = 'TMz'
                        self.data['zyx'][0, jj] = dvalue / (np.pi * 4e-4)
                        self.data['zyx'][1, jj] = derr / (np.pi * 4e-4)
                        self.data['zyx'][2, jj] = rvalue / (np.pi * 4e-4)
                        self.data['zyx'][3, jj] = rerr
                    if dlst[0] == 'ImagZyx' or dlst[0] == '116':
                        self.mode = 'TMz'
                        self.data['zyx'][0, jj] += 1j * dvalue / (np.pi * 4e-4)
                        self.data['zyx'][1, jj] = derr / (np.pi * 4e-4)
                        self.data['zyx'][2, jj] += 1j * rvalue / (np.pi * 4e-4)
                        self.data['zyx'][3, jj] = rerr
        if 'z' in self.mode:
            if 'TE' in self.mode:
                pol = 'xy'
            elif 'TM' in self.mode:
                pol = 'yx'
            for ii in [0, 2]:
                self.data['res' + pol][0 + ii] = 0.2 * np.abs(self.data['z' + pol][0 + ii]) ** 2. / self.freq
                self.data['phase' + pol][0 + ii] = np.rad2deg(
                    np.arctan(self.data['z' + pol][0 + ii].imag / self.data['z' + pol][0 + ii].real))

                self.data['res' + pol][1 + ii] = self.data['res' + pol][0 + ii] * self.data['z' + pol][
                    1 + ii].real / np.abs(self.data['z' + pol][0 + ii])

                for jjj in range(len(self.freq)):
                    self.data['phase' + pol][1 + ii, jjj] = \
                        mtcc.z_error2r_phi_error(self.data['z' + pol][0 + ii, jjj].real,
                                                 self.data['z' + pol][0 + ii, jjj].imag,
                                                 self.data['z' + pol][1 + ii, jjj].real)[1]
            if pol == 'xy':
                self.res_te = self.data['resxy']
                self.phase_te = self.data['phasexy']
            elif pol == 'yx':
                self.res_tm = self.data['resyx']
                self.phase_tm = self.data['phaseyx']


class Model(object):
    """
    read and write the model file fo Occam1D

    All depth measurements are in meters.

    ======================== ==================================================
    Attributes               Description
    ======================== ==================================================
    _model_fn                basename for model file *default* is Model1D
    _ss                      string spacing in model file *default* is 3*' '
    _string_fmt              format of model layers *default* is '.0f'
    air_layer_height         height of air layer *default* is 10000
    bottom_layer             bottom of the model *default* is 50000
    itdict                   dictionary of values from iteration file
    iter_fn                  full path to iteration file
    model_depth              array of model depths
    model_fn                 full path to model file
    model_penalty            array of penalties for each model layer
    model_preference_penalty array of model preference penalties for each layer
    model_prefernce          array of preferences for each layer
    model_res                array of resistivities for each layer
    n_layers                 number of layers in the model
    num_params               number of parameters to invert for (n_layers+2)
    pad_z                    padding of model at depth *default* is 5 blocks
    save_path                path to save files
    target_depth             depth of target to investigate
    z1_layer                 depth of first layer *default* is 10
    ======================== ==================================================

    ======================== ==================================================
    Methods                  Description
    ======================== ==================================================
    write_model_file         write an Occam1D model file, where depth increases
                             on a logarithmic scale
    read_model_file          read an Occam1D model file
    read_iter_file           read an .iter file output by Occam1D
    ======================== ==================================================

    :Example: ::

        >>> #--> make a model file
        >>> m1 = occam1d.Model()
        >>> m1.write_model_file(save_path=r"/home/occam1d/mt01/TE")
   """

    def __init__(self, model_fn=None, **kwargs):
        self.model_fn = model_fn
        self.iter_fn = None

        self.n_layers = kwargs.pop('n_layers', 100)
        self.bottom_layer = kwargs.pop('bottom_layer', None)
        self.target_depth = kwargs.pop('target_depth', None)
        self.pad_z = kwargs.pop('pad_z', 5)
        self.z1_layer = kwargs.pop('z1_layer', 10)
        self.air_layer_height = kwargs.pop('zir_layer_height', 10000)
        self._set_layerdepth_defaults()

        self.save_path = kwargs.pop('save_path', None)
        if self.model_fn is not None and self.save_path is None:
            self.save_path = os.path.dirname(self.model_fn)

        self._ss = ' ' * 3
        self._string_fmt = '.0f'
        self._model_fn = 'Model1D'
        self.model_res = None
        self.model_depth = None
        self.model_penalty = None
        self.model_prefernce = None
        self.model_preference_penalty = None
        self.num_params = None

    def _set_layerdepth_defaults(self, z1_threshold=3., bottomlayer_threshold=2.):
        """
        set target depth, bottom layer and z1 layer, making sure all the layers
        are consistent with each other and will work in the inversion
        (e.g. check target depth is not deeper than bottom layer)
        """

        if self.target_depth is None:
            if self.bottom_layer is None:
                # if neither target_depth nor bottom_layer are set, set defaults
                self.target_depth = 10000.
            else:
                self.target_depth = mtcc.roundsf(self.bottom_layer / 5., 1.)

        if self.bottom_layer is None:
            self.bottom_layer = 5. * self.target_depth
        # if bottom layer less than a factor of 2 greater than target depth then adjust deeper
        elif float(self.bottom_layer) / self.target_depth < bottomlayer_threshold:
            self.bottom_layer = bottomlayer_threshold * self.target_depth
            print("bottom layer not deep enough for target depth, set to {} m".format(self.bottom_layer))

        if self.z1_layer is None:
            self.z1_layer = mtcc.roundsf(self.target_depth / 1000., 0)
        elif self.target_depth / self.z1_layer < z1_threshold:
            self.z1_layer = self.target_depth / z1_threshold
            print("z1 layer not deep enough for target depth, set to {} m".format(self.z1_layer))

    def write_model_file(self, save_path=None, **kwargs):
        """
        Makes a 1D model file for Occam1D.

        Arguments:
        ----------

            **save_path** :path to save file to, if just path saved as
                          savepath\model.mod, if None defaults to dirpath

            **n_layers** : number of layers

            **bottom_layer** : depth of bottom layer in meters

            **target_depth** : depth to target under investigation

            **pad_z** : padding on bottom of model past target_depth

            **z1_layer** : depth of first layer in meters

            **air_layer_height** : height of air layers in meters

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
        if save_path is not None:
            self.save_path = save_path
            if os.path.isdir == False:
                os.mkdir(self.save_path)

        self.model_fn = os.path.join(self.save_path, self._model_fn)

        for key in list(kwargs.keys()):
            setattr(self, key, kwargs[key])

        if self.model_depth is None:
            # ---------create depth layers--------------------
            log_z = np.logspace(np.log10(self.z1_layer),
                                np.log10(self.target_depth -
                                         np.logspace(np.log10(self.z1_layer),
                                                     np.log10(self.target_depth),
                                                     num=self.n_layers)[-2]),
                                num=self.n_layers - self.pad_z)
            ztarget = np.array([zz - zz % 10 ** np.floor(np.log10(zz)) for zz in
                                log_z])
            log_zpad = np.logspace(np.log10(self.target_depth),
                                   np.log10(self.bottom_layer -
                                            np.logspace(np.log10(self.target_depth),
                                                        np.log10(self.bottom_layer),
                                                        num=self.pad_z)[-2]),
                                   num=self.pad_z)
            zpadding = np.array([zz - zz % 10 ** np.floor(np.log10(zz)) for zz in
                                 log_zpad])
            z_nodes = np.append(ztarget, zpadding)
            self.model_depth = np.array([z_nodes[:ii + 1].sum()
                                         for ii in range(z_nodes.shape[0])])
        else:
            self.n_layers = len(self.model_depth)

        self.num_params = self.n_layers + 2
        # make the model file
        modfid = open(self.model_fn, 'w')
        modfid.write('Format: Resistivity1DMod_1.0' + '\n')
        modfid.write('#LAYERS:    {0}\n'.format(self.num_params))
        modfid.write('!Set free values to -1 or ? \n')
        modfid.write('!penalize between 1 and 0,' +
                     '0 allowing jump between layers and 1 smooth. \n')
        modfid.write('!preference is the assumed resistivity on linear scale. \n')
        modfid.write('!pref_penalty needs to be put if preference is not 0 [0,1]. \n')
        modfid.write('! {0}\n'.format(self._ss.join(['top_depth', 'resistivity',
                                                     'penalty', 'preference',
                                                     'pref_penalty'])))
        modfid.write(self._ss.join([str(-self.air_layer_height),
                                    '1d12', '0', '0', '0', '!air layer', '\n']))
        modfid.write(self._ss.join(['0', '-1', '0', '0', '0',
                                    '!first ground layer', '\n']))
        for ll in self.model_depth:
            modfid.write(self._ss.join(['{0:{1}}'.format(np.ceil(ll),
                                                         self._string_fmt),
                                        '-1', '1', '0', '0', '\n']))

        modfid.close()

        print('Wrote Model file: {0}'.format(self.model_fn))

    def read_model_file(self, model_fn=None):
        """

        will read in model 1D file

        Arguments:
        ----------
            **modelfn** : full path to model file

        Fills attributes:
        --------

            * model_depth' : depth of model in meters

            * model_res : value of resisitivity

            * model_penalty : penalty

            * model_preference : preference

            * model_penalty_preference : preference penalty

        :Example: ::

            >>> m1 = occam1d.Model()
            >>> m1.savepath = r"/home/Occam1D/Line1/Inv1_TE"
            >>> m1.read_model_file()
        """
        if model_fn is not None:
            self.model_fn = model_fn
        if self.model_fn is None:
            raise IOError('Need to input a model file')
        elif os.path.isfile(self.model_fn) == False:
            raise IOError('Could not find{0}, check path'.format(self.model_fn))

        self._model_fn = os.path.basename(self.model_fn)
        self.save_path = os.path.dirname(self.model_fn)
        mfid = open(self.model_fn, 'r')
        mlines = mfid.readlines()
        mfid.close()
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
                mlst = mline.strip().split()
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

        # create an array with empty columns to put the TE and TM models into
        mres = np.zeros((len(mdict['res']), 2))
        mres[:, 0] = mdict['res']
        mdict['res'] = mres

        # make attributes
        self.model_res = mdict['res']
        self.model_depth = mdict['depth']
        self.model_penalty = mdict['pen']
        self.model_prefernce = mdict['pref']
        self.model_preference_penalty = mdict['prefpen']
        self.num_params = mdict['nparam']

    def read_iter_file(self, iter_fn=None, model_fn=None):
        """
        read an 1D iteration file

        Arguments:
        ----------
            **imode** : mode to read from

        Returns:
        --------
            **Occam1D.itdict** : dictionary with keys of the header:

            **model_res** : fills this array with the appropriate
                            values (0) for data, (1) for model

        :Example: ::

            >>> m1 = occam1d.Model()
            >>> m1.model_fn = r"/home/occam1d/mt01/TE/Model1D"
            >>> m1.read_iter_file(r"/home/Occam1D/Inv1_TE/M01TE_15.iter")

        """

        if iter_fn is not None:
            self.iter_fn = iter_fn

        if self.iter_fn is None:
            raise IOError('Need to input iteration file')

        if model_fn is not None:
            self.model_fn = model_fn
        if self.model_fn is None:
            raise IOError('Need to input a model file')
        else:
            self.read_model_file()

        freeparams = np.where(self.model_res == -1)[0]

        with open(self.iter_fn, 'r') as ifid:
            ilines = ifid.readlines()

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
        # for easy manipulation and access.
        model = np.array(model)
        self.model_res[freeparams, 1] = model


class Startup(object):
    """
    read and write input files for Occam1D

    ====================== ====================================================
    Attributes             Description
    ====================== ====================================================
    _ss                    string spacing
    _startup_fn            basename of startup file *default* is OccamStartup1D
    data_fn                full path to data file
    debug_level            debug level *default* is 1
    description            description of inversion for your self
                           *default* is 1D_Occam_Inv
    max_iter               maximum number of iterations *default* is 20
    model_fn               full path to model file
    rough_type             roughness type *default* is 1
    save_path              full path to save files to
    start_iter             first iteration number *default* is 0
    start_lagrange         starting lagrange number on log scale
                           *default* is 5
    start_misfit           starting misfit value *default* is 100
    start_rho              starting resistivity value (halfspace) in log scale
                           *default* is 100
    start_rough            starting roughness (ignored by Occam1D)
                           *default* is 1E7
    startup_fn             full path to startup file
    target_rms             target rms *default* is 1.0
    ====================== ====================================================
    """

    def __init__(self, data_fn=None, model_fn=None, **kwargs):
        self.data_fn = data_fn
        self.model_fn = model_fn

        if self.data_fn is not None:
            self.save_path = os.path.dirname(self.data_fn)
        elif self.model_fn is not None:
            self.save_path = os.path.dirname(self.model_fn)

        self.startup_fn = None
        self.rough_type = kwargs.pop('rough_type', 1)
        self.max_iter = kwargs.pop('max_iter', 20)
        self.target_rms = kwargs.pop('target_rms', 1)
        self.start_rho = kwargs.pop('start_rho', 100)
        self.description = kwargs.pop('description', '1D_Occam_Inv')
        self.start_lagrange = kwargs.pop('start_lagrange', 5.0)
        self.start_rough = kwargs.pop('start_rough', 1.0E7)
        self.debug_level = kwargs.pop('debug_level', 1)
        self.start_iter = kwargs.pop('start_iter', 0)
        self.start_misfit = kwargs.pop('start_misfit', 100)
        self.min_max_bounds = kwargs.pop('min_max_bounds', None)
        self.model_step = kwargs.pop('model_step', None)
        self._startup_fn = 'OccamStartup1D'
        self._ss = ' ' * 3
        
        
    def write_startup_file(self, save_path=None, **kwargs):
        """
        Make a 1D input file for Occam 1D

        Arguments:
        ---------
            **savepath** : full path to save input file to, if just path then
                           saved as savepath/input

            **model_fn** : full path to model file, if None then assumed to be in
                            savepath/model.mod

            **data_fn** : full path to data file, if None then assumed to be
                            in savepath/TE.dat or TM.dat

            **rough_type** : roughness type. *default* = 0

            **max_iter** : maximum number of iterations. *default* = 20

            **target_rms** : target rms value. *default* = 1.0

            **start_rho** : starting resistivity value on linear scale.
                            *default* = 100

            **description** : description of the inversion.

            **start_lagrange** : starting Lagrange multiplier for smoothness.
                           *default* = 5

            **start_rough** : starting roughness value. *default* = 1E7

            **debuglevel** : something to do with how Fortran debuggs the code
                             Almost always leave at *default* = 1

            **start_iter** : the starting iteration number, handy if the
                            starting model is from a previous run.
                            *default* = 0

            **start_misfit** : starting misfit value. *default* = 100

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

        if save_path is not None:
            self.save_path = save_path
        if os.path.isdir(self.save_path) == False:
            os.mkdir(self.save_path)

        self.startup_fn = os.path.join(self.save_path, self._startup_fn)

        # --> read data file
        if self.data_fn is None:
            raise IOError('Need to input data file name.')
        else:
            data = Data()
            data.read_data_file(self.data_fn)

        # --> read model file
        if self.model_fn is None:
            raise IOError('Need to input model file name.')
        else:
            model = Model()
            model.read_model_file(self.model_fn)

            # --> get any keywords
        for key in list(kwargs.keys()):
            setattr(self, key, kwargs[key])

            # --> write input file
        infid = open(self.startup_fn, 'w')
        infid.write('{0:<21}{1}\n'.format('Format:', 'OCCAMITER_FLEX'))
        infid.write('{0:<21}{1}\n'.format('Description:', self.description))
        infid.write('{0:<21}{1}\n'.format('Model File:',
                                          os.path.basename(self.model_fn)))
        infid.write('{0:<21}{1}\n'.format('Data File:',
                                          os.path.basename(self.data_fn)))
        infid.write('{0:<21}{1}\n'.format('Date/Time:', time.ctime()))
        infid.write('{0:<21}{1}\n'.format('Max Iter:', self.max_iter))
        infid.write('{0:<21}{1}\n'.format('Target Misfit:', self.target_rms))
        infid.write('{0:<21}{1}\n'.format('Roughness Type:', self.rough_type))
        if self.min_max_bounds == None:
            infid.write('{0:<21}{1}\n'.format('!Model Bounds:', 'min,max'))
        else:
            infid.write('{0:<21}{1},{2}\n'.format('Model Bounds:',
                                                  self.min_max_bounds[0],
                                                  self.min_max_bounds[1]))
        if self.model_step == None:
            infid.write('{0:<21}{1}\n'.format('!Model Value Steps:',
                                              'stepsize'))
        else:
            infid.write('{0:<21}{1}\n'.format('Model Value Steps:',
                                              self.model_step))
        infid.write('{0:<21}{1}\n'.format('Debug Level:', self.debug_level))
        infid.write('{0:<21}{1}\n'.format('Iteration:', self.start_iter))
        infid.write('{0:<21}{1}\n'.format('Lagrange Value:', self.start_lagrange))
        infid.write('{0:<21}{1}\n'.format('Roughness Value:', self.start_rough))
        infid.write('{0:<21}{1}\n'.format('Misfit Value:', self.start_misfit))
        infid.write('{0:<21}{1}\n'.format('Misfit Reached:', 0))
        infid.write('{0:<21}{1}\n'.format('Param Count:', model.num_params))

        for ii in range(model.num_params):
            infid.write('{0}{1:.2f}\n'.format(self._ss,
                                              np.log10(self.start_rho)))

        infid.close()
        print('Wrote Input File: {0}'.format(self.startup_fn))

    def read_startup_file(self, startup_fn):
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
        if startup_fn is not None:
            self.startup_fn = startup_fn

        if self.startup_fn is None:
            raise IOError('Need to input a startup file.')

        self._startup_fn = os.path.basename(self.startup_fn)
        self.save_path = os.path.dirname(self.startup_fn)

        infid = open(self.startup_fn, 'r')
        ilines = infid.readlines()
        infid.close()

        self.indict = {}
        res = []

        # split the keys and values from the header information
        for iline in ilines:
            if iline.find(':') >= 0:
                ikey = iline[0:20].strip()[:-1]
                ivalue = iline[20:].split('!')[0].strip()
                if ikey.find('!') == 0:
                    pass
                else:
                    setattr(self, ikey.lower().replace(' ', '_'), ivalue)
                self.indict[ikey[:-1]] = ivalue
            else:
                try:
                    res.append(float(iline.strip()))
                except ValueError:
                    pass

        # make the resistivity array ready for models to be input
        self.indict['res'] = np.zeros((len(res), 3))
        self.indict['res'][:, 0] = res


class Plot1DResponse(object):
    """
    plot the 1D response and model.  Plots apparent resisitivity and phase
    in different subplots with the model on the far right.  You can plot both
    TE and TM modes together along with different iterations of the model.
    These will be plotted in different colors or shades of gray depneng on
    color_scale.

    :Example: ::

        >>> import mtpy.modeling.occam1d as occam1d
        >>> p1 = occam1d.Plot1DResponse(plot_yn='n')
        >>> p1.data_te_fn = r"/home/occam1d/mt01/TE/Occam_DataFile_TE.dat"
        >>> p1.data_tm_fn = r"/home/occam1d/mt01/TM/Occam_DataFile_TM.dat"
        >>> p1.model_fn = r"/home/occam1d/mt01/TE/Model1D"
        >>> p1.iter_te_fn = [r"/home/occam1d/mt01/TE/TE_{0}.iter".format(ii)
        >>> ...              for ii in range(5,10)]
        >>> p1.iter_tm_fn = [r"/home/occam1d/mt01/TM/TM_{0}.iter".format(ii)
        >>> ...              for ii in range(5,10)]
        >>> p1.resp_te_fn = [r"/home/occam1d/mt01/TE/TE_{0}.resp".format(ii)
        >>> ...              for ii in range(5,10)]
        >>> p1.resp_tm_fn = [r"/home/occam1d/mt01/TM/TM_{0}.resp".format(ii)
        >>> ...              for ii in range(5,10)]
        >>> p1.plot()

    ==================== ======================================================
    Attributes           Description
    ==================== ======================================================
    axm                  matplotlib.axes instance for model subplot
    axp                  matplotlib.axes instance for phase subplot
    axr                  matplotlib.axes instance for app. res subplot
    color_mode           [ 'color' | 'bw' ]
    cted                 color of TE data markers
    ctem                 color of TM data markers
    ctmd                 color of TE model markers
    ctmm                 color of TM model markers
    data_te_fn           full path to data file for TE mode
    data_tm_fn           full path to data file for TM mode
    depth_limits         (min, max) limits for depth plot in depth_units
    depth_scale          [ 'log' | 'linear' ] *default* is linear
    depth_units          [ 'm' | 'km' ] *default is 'km'
    e_capsize            capsize of error bars
    e_capthick           cap thickness of error bars
    fig                  matplotlib.figure instance for plot
    fig_dpi              resolution in dots-per-inch for figure
    fig_num              number of figure instance
    fig_size             size of figure in inches [width, height]
    font_size            size of axes tick labels, axes labels are +2
    grid_alpha           transparency of grid
    grid_color           color of grid
    iter_te_fn           full path or list of .iter files for TE mode
    iter_tm_fn           full path or list of .iter files for TM mode
    lw                   width of lines for model
    model_fn             full path to model file
    ms                   marker size
    mted                 marker for TE data
    mtem                 marker for TM data
    mtmd                 marker for TE model
    mtmm                 marker for TM model
    phase_limits         (min, max) limits on phase in degrees
    phase_major_ticks    spacing for major ticks in phase
    phase_minor_ticks    spacing for minor ticks in phase
    plot_yn              [ 'y' | 'n' ] plot on instantiation
    res_limits           limits of resistivity in linear scale
    resp_te_fn           full path or list of .resp files for TE mode
    resp_tm_fn           full path or list of .iter files for TM mode
    subplot_bottom       spacing of subplots from bottom of figure
    subplot_hspace       height spacing between subplots
    subplot_left         spacing of subplots from left of figure
    subplot_right        spacing of subplots from right of figure
    subplot_top          spacing of subplots from top of figure
    subplot_wspace       width spacing between subplots
    title_str            title of plot
    ==================== ======================================================

    """

    def __init__(self, data_te_fn=None, data_tm_fn=None, model_fn=None,
                 resp_te_fn=None, resp_tm_fn=None, iter_te_fn=None,
                 iter_tm_fn=None, **kwargs):
        self.data_te_fn = data_te_fn
        self.data_tm_fn = data_tm_fn

        self.model_fn = model_fn

        self.override_legend_subscript = kwargs.pop('override_legend_subscript',None)
        self.resp_te_fn = resp_te_fn
        if type(self.resp_te_fn) is not list:
            self.resp_te_fn = [self.resp_te_fn]

        self.resp_tm_fn = resp_tm_fn
        if type(self.resp_tm_fn) is not list:
            self.resp_tm_fn = [self.resp_tm_fn]

        self.iter_te_fn = iter_te_fn
        if type(self.iter_te_fn) is not list:
            self.iter_te_fn = [self.iter_te_fn]

        self.iter_tm_fn = iter_tm_fn
        if type(self.iter_tm_fn) is not list:
            self.iter_tm_fn = [self.iter_tm_fn]

        self.color_mode = kwargs.pop('color_mode', 'color')

        self.ms = kwargs.pop('ms', 1.5)
        self.lw = kwargs.pop('lw', .5)
        self.ls = kwargs.pop('ls', ':')
        self.e_capthick = kwargs.pop('e_capthick', .5)
        self.e_capsize = kwargs.pop('e_capsize', 2)

        self.phase_major_ticks = kwargs.pop('phase_major_ticks', 10)
        self.phase_minor_ticks = kwargs.pop('phase_minor_ticks', 5)

        self.grid_color = kwargs.pop('grid_color', (.25, .25, .25))
        self.grid_alpha = kwargs.pop('grid_alpha', .3)

        # color mode
        if self.color_mode == 'color':
            # color for data
            self.cted = kwargs.pop('cted', (0, 0, 1))
            self.ctmd = kwargs.pop('ctmd', (1, 0, 0))
            self.mted = kwargs.pop('mted', 's')
            self.mtmd = kwargs.pop('mtmd', 'o')

            # color for occam2d model
            self.ctem = kwargs.pop('ctem', (0, .6, .3))
            self.ctmm = kwargs.pop('ctmm', (.9, 0, .8))
            self.mtem = kwargs.pop('mtem', '+')
            self.mtmm = kwargs.pop('mtmm', '+')

        # black and white mode
        elif self.color_mode == 'bw':
            # color for data
            self.cted = kwargs.pop('cted', (0, 0, 0))
            self.ctmd = kwargs.pop('ctmd', (0, 0, 0))
            self.mted = kwargs.pop('mted', '*')
            self.mtmd = kwargs.pop('mtmd', 'v')

            # color for occam2d model
            self.ctem = kwargs.pop('ctem', (0.6, 0.6, 0.6))
            self.ctmm = kwargs.pop('ctmm', (0.6, 0.6, 0.6))
            self.mtem = kwargs.pop('mtem', '+')
            self.mtmm = kwargs.pop('mtmm', 'x')

        self.phase_limits = kwargs.pop('phase_limits', (-5, 95))
        self.res_limits = kwargs.pop('res_limits', None)
        self.depth_limits = kwargs.pop('depth_limits', None)
        self.depth_scale = kwargs.pop('depth_scale', 'linear')
        self.depth_units = kwargs.pop('depth_units', 'km')

        self.fig_num = kwargs.pop('fig_num', 1)
        self.fig_size = kwargs.pop('fig_size', [6, 6])
        self.fig_dpi = kwargs.pop('dpi', 300)
        self.fig = None
        self.axr = None
        self.axp = None
        self.axm = None

        self.subplot_wspace = .25
        self.subplot_hspace = .15
        self.subplot_right = .92
        self.subplot_left = .085
        self.subplot_top = .93
        self.subplot_bottom = .1

        self.font_size = kwargs.pop('font_size', 6)

        self.title_str = kwargs.pop('title_str', '')
        self.plot_yn = kwargs.pop('plot_yn', 'y')

        if self.plot_yn == 'y':
            self.plot()

    def plot(self):
        """
        plot data, response and model

        """
        if type(self.resp_te_fn) is not list:
            self.resp_te_fn = [self.resp_te_fn]

        if type(self.resp_tm_fn) is not list:
            self.resp_tm_fn = [self.resp_tm_fn]

        if type(self.iter_te_fn) is not list:
            self.iter_te_fn = [self.iter_te_fn]

        if type(self.iter_tm_fn) is not list:
            self.iter_tm_fn = [self.iter_tm_fn]

        # make a grid of subplots
        gs = gridspec.GridSpec(6, 5, hspace=self.subplot_hspace,
                               wspace=self.subplot_wspace)

        # make a figure
        self.fig = plt.figure(self.fig_num, self.fig_size, dpi=self.fig_dpi)
        plt.clf()

        # set some plot parameters
        plt.rcParams['font.size'] = self.font_size
        plt.rcParams['figure.subplot.left'] = self.subplot_left
        plt.rcParams['figure.subplot.right'] = self.subplot_right
        plt.rcParams['figure.subplot.bottom'] = self.subplot_bottom
        plt.rcParams['figure.subplot.top'] = self.subplot_top

        # subplot resistivity
        self.axr = self.fig.add_subplot(gs[:4, :4])

        # subplot for phase
        self.axp = self.fig.add_subplot(gs[4:, :4], sharex=self.axr)

        # subplot for model
        self.axm = self.fig.add_subplot(gs[:, 4])

        legend_marker_list_te = []
        legend_label_list_te = []
        legend_marker_list_tm = []
        legend_label_list_tm = []
        # --> plot data apparent resistivity and phase-------------------------
        if self.data_te_fn is not None:
            d1 = Data()
            d1.read_data_file(self.data_te_fn)

            # --> cut out missing data
            rxy = np.where(d1.res_te[0] != 0)[0]

            # --> TE mode Data
            if len(rxy) > 0:
                rte = self.axr.errorbar(1. / d1.freq[rxy],
                                        d1.res_te[0][rxy],
                                        ls=self.ls,
                                        marker=self.mted,
                                        ms=self.ms,
                                        mfc=self.cted,
                                        mec=self.cted,
                                        color=self.cted,
                                        yerr=d1.res_te[1][rxy],
                                        ecolor=self.cted,
                                        picker=2,
                                        lw=self.lw,
                                        elinewidth=self.lw,
                                        capsize=self.e_capsize,
                                        capthick=self.e_capthick)
                legend_marker_list_te.append(rte[0])
                if self.override_legend_subscript is not None:
                    legend_label_list_tm.append('$Obs_{'+str.upper(self.override_legend_subscript)+'}$')
                else:
                    legend_label_list_te.append('$Obs_{TM}$')
            else:
                pass
            # --------------------plot phase--------------------------------
            # cut out missing data points first
            pxy = np.where(d1.phase_te[0] != 0)[0]

            # --> TE mode data
            if len(pxy) > 0:
                self.axp.errorbar(1. / d1.freq[pxy],
                                  d1.phase_te[0][pxy],
                                  ls=self.ls,
                                  marker=self.mted,
                                  ms=self.ms,
                                  mfc=self.cted,
                                  mec=self.cted,
                                  color=self.cted,
                                  yerr=d1.phase_te[1][pxy],
                                  ecolor=self.cted,
                                  picker=1,
                                  lw=self.lw,
                                  elinewidth=self.lw,
                                  capsize=self.e_capsize,
                                  capthick=self.e_capthick)
            else:
                pass
        # --> plot tm data------------------------------------------------------
        if self.data_tm_fn is not None:
            d1 = Data()
            d1.read_data_file(self.data_tm_fn)

            ryx = np.where(d1.res_tm[0] != 0)[0]

            # --> TM mode data
            if len(ryx) > 0:
                rtm = self.axr.errorbar(1. / d1.freq[ryx],
                                        d1.res_tm[0][ryx],
                                        ls=self.ls,
                                        marker=self.mtmd,
                                        ms=self.ms,
                                        mfc=self.ctmd,
                                        mec=self.ctmd,
                                        color=self.ctmd,
                                        yerr=d1.res_tm[1][ryx],
                                        ecolor=self.ctmd,
                                        picker=2,
                                        lw=self.lw,
                                        elinewidth=self.lw,
                                        capsize=self.e_capsize,
                                        capthick=self.e_capthick)
                legend_marker_list_tm.append(rtm[0])
                if self.override_legend_subscript is not None:
                    legend_label_list_tm.append('$Obs_{'+str.upper(self.override_legend_subscript)+'}$')
                else:
                    legend_label_list_te.append('$Obs_{TM}$')
            else:
                pass

                # --------------------plot phase--------------------------------
            # cut out missing data points first
            pyx = np.where(d1.phase_tm[0] != 0)[0]

            # --> TM mode data
            if len(pyx) > 0:
                self.axp.errorbar(1. / d1.freq[pyx],
                                  d1.phase_tm[0][pyx],
                                  ls=self.ls,
                                  marker=self.mtmd,
                                  ms=self.ms,
                                  mfc=self.ctmd,
                                  mec=self.ctmd,
                                  color=self.ctmd,
                                  yerr=d1.phase_tm[1][pyx],
                                  ecolor=self.ctmd,
                                  picker=1,
                                  lw=self.lw,
                                  elinewidth=self.lw,
                                  capsize=self.e_capsize,
                                  capthick=self.e_capthick)
            else:
                pass

        # --> plot model apparent resistivity and phase-------------------------
        nr = len(self.resp_te_fn)
        for rr, rfn in enumerate(self.resp_te_fn):
            if rfn is None:
                break
            # accommodate larger number of iterations that might have > 2 digits
            itnum = rfn[-8:-5]
            while not str.isdigit(itnum[0]):
                itnum = itnum[1:]
                if itnum == '':
                    break
            if self.color_mode == 'color':
                cxy = (0, .4 + float(rr) / (3 * nr), 0)
            elif self.color_mode == 'bw':
                cxy = (1 - 1.25 / (rr + 2.), 1 - 1.25 / (rr + 2.), 1 - 1.25 / (rr + 2.))

            d1 = Data()

            d1.read_resp_file(rfn, data_fn=self.data_te_fn)

            # get non zero data
            rxy = np.where(d1.res_te[2] != 0)[0]

            # --> TE mode Data
            if len(rxy) > 0:
                rte = self.axr.errorbar(1. / d1.freq[rxy],
                                        d1.res_te[2][rxy],
                                        ls=self.ls,
                                        marker=self.mtem,
                                        ms=self.ms,
                                        mfc=cxy,
                                        mec=cxy,
                                        color=cxy,
                                        yerr=None,
                                        ecolor=cxy,
                                        picker=2,
                                        lw=self.lw,
                                        elinewidth=self.lw,
                                        capsize=self.e_capsize,
                                        capthick=self.e_capthick)
                legend_marker_list_te.append(rte[0])
                if self.override_legend_subscript is not None:
                    legend_label_list_tm.append('$Mod_{'+str.upper(self.override_legend_subscript)+'}$' + itnum)
                else:
                    legend_label_list_te.append('$Mod_{TE}$' + itnum)
            else:
                pass

            # --------------------plot phase--------------------------------
            # cut out missing data points first
            # --> data
            pxy = np.where(d1.phase_te[2] != 0)[0]

            # --> TE mode phase
            if len(pxy) > 0:
                self.axp.errorbar(1. / d1.freq[pxy],
                                  d1.phase_te[2][pxy],
                                  ls=self.ls,
                                  marker=self.mtem,
                                  ms=self.ms,
                                  mfc=cxy,
                                  mec=cxy,
                                  color=cxy,
                                  yerr=None,
                                  ecolor=cxy,
                                  picker=1,
                                  lw=self.lw,
                                  elinewidth=self.lw,
                                  capsize=self.e_capsize,
                                  capthick=self.e_capthick)
            else:
                pass
        # ---------------plot TM model response---------------------------------
        nr = len(self.resp_tm_fn)
        for rr, rfn in enumerate(self.resp_tm_fn):
            if rfn is None:
                break
            # accommodate larger number of iterations that might have > 2 digits
            itnum = rfn[-8:-5]
            while not str.isdigit(itnum[0]):
                itnum = itnum[1:]
                if itnum == '':
                    break
            if self.color_mode == 'color':
                cyx = (.7 + float(rr) / (4 * nr), .13, .63 - float(rr) / (4 * nr))
            elif self.color_mode == 'bw':
                cyx = (1 - 1.25 / (rr + 2.), 1 - 1.25 / (rr + 2.), 1 - 1.25 / (rr + 2.))
            d1 = Data()

            d1.read_resp_file(rfn, data_fn=self.data_tm_fn)
            ryx = np.where(d1.res_tm[2] != 0)[0]
            # --> TM mode model
            if len(ryx) > 0:
                rtm = self.axr.errorbar(1. / d1.freq[ryx],
                                        d1.res_tm[2][ryx],
                                        ls=self.ls,
                                        marker=self.mtmm,
                                        ms=self.ms,
                                        mfc=cyx,
                                        mec=cyx,
                                        color=cyx,
                                        yerr=None,
                                        ecolor=cyx,
                                        picker=2,
                                        lw=self.lw,
                                        elinewidth=self.lw,
                                        capsize=self.e_capsize,
                                        capthick=self.e_capthick)
                legend_marker_list_tm.append(rtm[0])
                if self.override_legend_subscript is not None:
                    legend_label_list_tm.append('$Mod_{'+str.upper(self.override_legend_subscript)+'}$' + itnum)
                else:
                    legend_label_list_te.append('$Mod_{TM}$' + itnum)
            else:
                pass

            pyx = np.where(d1.phase_tm[2] != 0)[0]

            # --> TM mode model
            if len(pyx) > 0:
                self.axp.errorbar(1. / d1.freq[pyx],
                                  d1.phase_tm[0][pyx],
                                  ls=self.ls,
                                  marker=self.mtmm,
                                  ms=self.ms,
                                  mfc=cyx,
                                  mec=cyx,
                                  color=cyx,
                                  yerr=None,
                                  ecolor=cyx,
                                  picker=1,
                                  lw=self.lw,
                                  elinewidth=self.lw,
                                  capsize=self.e_capsize,
                                  capthick=self.e_capthick)
            else:
                pass

        # --> set axis properties-----------------------------------------------
        self.axr.set_xscale('log', nonposx='clip')
        self.axp.set_xscale('log', nonposx='clip')
        self.axr.set_yscale('log', nonposy='clip')
        self.axr.grid(True, alpha=self.grid_alpha, which='both',
                      color=self.grid_color)
        plt.setp(self.axr.xaxis.get_ticklabels(), visible=False)
        self.axp.grid(True, alpha=self.grid_alpha, which='both',
                      color=self.grid_color)
        self.axp.yaxis.set_major_locator(MultipleLocator(self.phase_major_ticks))
        self.axp.yaxis.set_minor_locator(MultipleLocator(self.phase_minor_ticks))

        if self.res_limits is not None:
            self.axr.set_ylim(self.res_limits)

        self.axp.set_ylim(self.phase_limits)
        self.axr.set_ylabel('App. Res. ($\Omega \cdot m$)',
                            fontdict={'size': self.font_size, 'weight': 'bold'})
        self.axp.set_ylabel('Phase (deg)',
                            fontdict={'size': self.font_size, 'weight': 'bold'})
        self.axp.set_xlabel('Period (s)',
                            fontdict={'size': self.font_size, 'weight': 'bold'})
        plt.suptitle(self.title_str, fontsize=self.font_size + 2, fontweight='bold')
        if legend_marker_list_te == [] or legend_marker_list_tm == []:
            num_col = 1
        else:
            num_col = 2
        self.axr.legend(legend_marker_list_te + legend_marker_list_tm,
                        legend_label_list_te + legend_label_list_tm,
                        loc=2, markerscale=1,
                        borderaxespad=.05,
                        labelspacing=.08,
                        handletextpad=.15,
                        borderpad=.05,
                        ncol=num_col,
                        prop={'size': self.font_size + 1})

        # --> plot depth model--------------------------------------------------
        if self.model_fn is not None:
            # put axis labels on the right side for clarity
            self.axm.yaxis.set_label_position('right')
            self.axm.yaxis.set_tick_params(left='off', right='on',
                                           labelright='on')
            self.axm.yaxis.tick_right()

            if self.depth_units == 'km':
                dscale = 1000.
            else:
                dscale = 1.

            # --> plot te models
            nr = len(self.iter_te_fn)
            for ii, ifn in enumerate(self.iter_te_fn):
                if ifn is None:
                    break
                if self.color_mode == 'color':
                    cxy = (0, .4 + float(ii) / (3 * nr), 0)
                elif self.color_mode == 'bw':
                    cxy = (1 - 1.25 / (ii + 2.), 1 - 1.25 / (ii + 2.), 1 - 1.25 / (ii + 2.))
                m1 = Model()
                m1.read_iter_file(ifn, self.model_fn)
                plot_depth = m1.model_depth[1:] / dscale
                plot_model = abs(10 ** m1.model_res[1:, 1])
                self.axm.semilogx(plot_model[::-1],
                                  plot_depth[::-1],
                                  ls='-',
                                  color=cxy,
                                  lw=self.lw)

            # --> plot TM models
            nr = len(self.iter_tm_fn)
            for ii, ifn in enumerate(self.iter_tm_fn):
                if ifn is None:
                    break
                if self.color_mode == 'color':
                    cyx = (.7 + float(ii) / (4 * nr), .13, .63 - float(ii) / (4 * nr))
                elif self.color_mode == 'bw':
                    cyx = (1 - 1.25 / (ii + 2.), 1 - 1.25 / (ii + 2.), 1 - 1.25 / (ii + 2.))
                m1 = Model()
                m1.read_iter_file(ifn, self.model_fn)
                plot_depth = m1.model_depth[1:] / dscale
                plot_model = abs(10 ** m1.model_res[1:, 1])
                self.axm.semilogx(plot_model[::-1],
                                  plot_depth[::-1],
                                  ls='steps-',
                                  color=cyx,
                                  lw=self.lw)

            m1 = Model()
            m1.read_model_file(self.model_fn)
            if self.depth_limits is None:
                dmin = min(plot_depth)
                if dmin == 0:
                    dmin = 1
                dmax = max(plot_depth)
                self.depth_limits = (dmin, dmax)

            self.axm.set_ylim(ymin=max(self.depth_limits),
                              ymax=min(self.depth_limits))
            if self.depth_scale == 'log':
                self.axm.set_yscale('log', nonposy='clip')
            self.axm.set_ylabel('Depth ({0})'.format(self.depth_units),
                                fontdict={'size': self.font_size, 'weight': 'bold'})
            self.axm.set_xlabel('Resistivity ($\Omega \cdot m$)',
                                fontdict={'size': self.font_size, 'weight': 'bold'})
            self.axm.grid(True, which='both', alpha=.25)

        plt.show()

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

    def save_figure(self, save_fn, file_format='pdf', orientation='portrait',
                    fig_dpi=None, close_plot='y'):
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

        if fig_dpi is None:
            fig_dpi = self.fig_dpi

        if not os.path.isdir(save_fn):
            file_format = save_fn[-3:]
            self.fig.savefig(save_fn, dpi=fig_dpi, format=file_format,
                             orientation=orientation, bbox_inches='tight')

        else:
            save_fn = os.path.join(save_fn, 'Occam1d.' +
                                   file_format)
            self.fig.savefig(save_fn, dpi=fig_dpi, format=file_format,
                             orientation=orientation, bbox_inches='tight')

        if close_plot == 'y':
            plt.clf()
            plt.close(self.fig)

        else:
            pass

        self.fig_fn = save_fn
        print('Saved figure to: ' + self.fig_fn)

    def __str__(self):
        """
        rewrite the string builtin to give a useful message
        """

        return "Plots model responses and model for 1D occam inversion"


class Run(object):
    """
    run occam 1d from python given the correct files and location of occam1d
    executable

    """

    def __init__(self, startup_fn=None, occam_path=None, **kwargs):
        self.startup_fn = startup_fn
        self.occam_path = occam_path
        self.mode = kwargs.pop('mode', 'TE')

        self.run_occam1d()

    def run_occam1d(self):

        if self.startup_fn is None:
            raise IOError('Need to input startup file')
        if self.occam_path is None:
            raise IOError('Need to input path to occam1d executable')

        os.chdir(os.path.dirname(self.startup_fn))
        test = subprocess.call([self.occam_path,
                                os.path.basename(self.startup_fn),
                                self.mode])
        if test == 0:
            print('=========== Ran Inversion ==========')
            print('  check {0} for files'.format(os.path.dirname(self.startup_fn)))


class PlotL2(object):
    """
    plot L2 curve of iteration vs rms and roughness

    Arguments:
    ----------
        **rms_arr** : structured array with keys:
                      * 'iteration' --> for iteration number (int)
                      * 'rms' --> for rms (float)
                      * 'roughness' --> for roughness (float)

    ======================= ===================================================
    Keywords/attributes     Description
    ======================= ===================================================
    ax1                     matplotlib.axes instance for rms vs iteration
    ax2                     matplotlib.axes instance for roughness vs rms
    fig                     matplotlib.figure instance
    fig_dpi                 resolution of figure in dots-per-inch
    fig_num                 number of figure instance
    fig_size                size of figure in inches (width, height)
    font_size               size of axes tick labels, axes labels is +2
    plot_yn                 [ 'y' | 'n']
                            'y' --> to plot on instantiation
                            'n' --> to not plot on instantiation
    rms_arr                 structure np.array as described above
    rms_color               color of rms marker and line
    rms_lw                  line width of rms line
    rms_marker              marker for rms values
    rms_marker_size         size of marker for rms values
    rms_mean_color          color of mean line
    rms_median_color        color of median line
    rough_color             color of roughness line and marker
    rough_font_size         font size for iteration number inside roughness
                            marker
    rough_lw                line width for roughness line
    rough_marker            marker for roughness
    rough_marker_size       size of marker for roughness
    subplot_bottom          subplot spacing from bottom
    subplot_left            subplot spacing from left
    subplot_right           subplot spacing from right
    subplot_top             subplot spacing from top
    ======================= ===================================================

    =================== =======================================================
    Methods             Description
    =================== =======================================================
    plot                plots L2 curve.
    redraw_plot         call redraw_plot to redraw the figures,
                        if one of the attributes has been changed
    save_figure         saves the matplotlib.figure instance to desired
                        location and format
    =================== ======================================================

    """

    def __init__(self, dir_path, model_fn, **kwargs):
        self.dir_path = dir_path
        self.model_fn = model_fn
        self._get_iter_list()

        self.subplot_right = .98
        self.subplot_left = .085
        self.subplot_top = .91
        self.subplot_bottom = .1

        self.fig_num = kwargs.pop('fig_num', 1)
        self.fig_size = kwargs.pop('fig_size', [6, 6])
        self.fig_dpi = kwargs.pop('dpi', 300)
        self.font_size = kwargs.pop('font_size', 8)

        self.rms_lw = kwargs.pop('rms_lw', 1)
        self.rms_marker = kwargs.pop('rms_marker', 'd')
        self.rms_color = kwargs.pop('rms_color', 'k')
        self.rms_marker_size = kwargs.pop('rms_marker_size', 5)
        self.rms_median_color = kwargs.pop('rms_median_color', 'red')
        self.rms_mean_color = kwargs.pop('rms_mean_color', 'orange')

        self.rough_lw = kwargs.pop('rough_lw', .75)
        self.rough_marker = kwargs.pop('rough_marker', 'o')
        self.rough_color = kwargs.pop('rough_color', 'b')
        self.rough_marker_size = kwargs.pop('rough_marker_size', 7)
        self.rough_font_size = kwargs.pop('rough_font_size', 6)

        self.plot_yn = kwargs.pop('plot_yn', 'y')
        if self.plot_yn == 'y':
            self.plot()

    def _get_iter_list(self):
        """
        get all iteration files in dir_path
        """

        if os.path.isdir(self.dir_path) == False:
            raise IOError('Could not find {0}'.format(self.dir_path))

        iter_list = [os.path.join(self.dir_path, fn)
                     for fn in os.listdir(self.dir_path)
                     if fn.find('.iter') > 0]

        self.rms_arr = np.zeros(len(iter_list),
                                dtype=np.dtype([('iteration', np.int),
                                                ('rms', np.float),
                                                ('roughness', np.float)]))
        for ii, fn in enumerate(iter_list):
            m1 = Model()
            m1.read_iter_file(fn, self.model_fn)
            self.rms_arr[ii]['iteration'] = int(m1.itdict['Iteration'])
            self.rms_arr[ii]['rms'] = float(m1.itdict['Misfit Value'])
            self.rms_arr[ii]['roughness'] = float(m1.itdict['Roughness Value'])

        self.rms_arr.sort(order='iteration')

    def plot(self):
        """
        plot L2 curve
        """

        nr = self.rms_arr.shape[0]
        med_rms = np.median(self.rms_arr['rms'])
        mean_rms = np.mean(self.rms_arr['rms'])

        # set the dimesions of the figure
        plt.rcParams['font.size'] = self.font_size
        plt.rcParams['figure.subplot.left'] = self.subplot_left
        plt.rcParams['figure.subplot.right'] = self.subplot_right
        plt.rcParams['figure.subplot.bottom'] = self.subplot_bottom
        plt.rcParams['figure.subplot.top'] = self.subplot_top

        # make figure instance
        self.fig = plt.figure(self.fig_num, self.fig_size, dpi=self.fig_dpi)
        plt.clf()

        # make a subplot for RMS vs Iteration
        self.ax1 = self.fig.add_subplot(1, 1, 1)

        # plot the rms vs iteration
        l1, = self.ax1.plot(self.rms_arr['iteration'],
                            self.rms_arr['rms'],
                            '-k',
                            lw=1,
                            marker='d',
                            ms=5)

        # plot the median of the RMS
        m1, = self.ax1.plot(self.rms_arr['iteration'],
                            np.repeat(med_rms, nr),
                            ls='--',
                            color=self.rms_median_color,
                            lw=self.rms_lw * .75)

        # plot the mean of the RMS
        m2, = self.ax1.plot(self.rms_arr['iteration'],
                            np.repeat(mean_rms, nr),
                            ls='--',
                            color=self.rms_mean_color,
                            lw=self.rms_lw * .75)

        # make subplot for RMS vs Roughness Plot
        self.ax2 = self.ax1.twiny()

        self.ax2.set_xlim(self.rms_arr['roughness'][1:].min(),
                          self.rms_arr['roughness'][1:].max())

        self.ax1.set_ylim(0, self.rms_arr['rms'][1])

        # plot the rms vs roughness
        l2, = self.ax2.plot(self.rms_arr['roughness'],
                            self.rms_arr['rms'],
                            ls='--',
                            color=self.rough_color,
                            lw=self.rough_lw,
                            marker=self.rough_marker,
                            ms=self.rough_marker_size,
                            mfc='white')

        # plot the iteration number inside the roughness marker
        for rms, ii, rough in zip(self.rms_arr['rms'], self.rms_arr['iteration'],
                                  self.rms_arr['roughness']):
            # need this because if the roughness is larger than this number
            # matplotlib puts the text out of bounds and a draw_text_image
            # error is raised and file cannot be saved, also the other
            # numbers are not put in.
            if rough > 1e8:
                pass
            else:
                self.ax2.text(rough,
                              rms,
                              '{0}'.format(ii),
                              horizontalalignment='center',
                              verticalalignment='center',
                              fontdict={'size': self.rough_font_size,
                                        'weight': 'bold',
                                        'color': self.rough_color})

        # make a legend
        self.ax1.legend([l1, l2, m1, m2],
                        ['RMS', 'Roughness',
                         'Median_RMS={0:.2f}'.format(med_rms),
                         'Mean_RMS={0:.2f}'.format(mean_rms)],
                        ncol=1,
                        loc='upper right',
                        columnspacing=.25,
                        markerscale=.75,
                        handletextpad=.15)

        # set the axis properties for RMS vs iteration
        self.ax1.yaxis.set_minor_locator(MultipleLocator(.1))
        self.ax1.xaxis.set_minor_locator(MultipleLocator(1))
        self.ax1.set_ylabel('RMS',
                            fontdict={'size': self.font_size + 2,
                                      'weight': 'bold'})
        self.ax1.set_xlabel('Iteration',
                            fontdict={'size': self.font_size + 2,
                                      'weight': 'bold'})
        self.ax1.grid(alpha=.25, which='both', lw=self.rough_lw)
        self.ax2.set_xlabel('Roughness',
                            fontdict={'size': self.font_size + 2,
                                      'weight': 'bold',
                                      'color': self.rough_color})

        for t2 in self.ax2.get_xticklabels():
            t2.set_color(self.rough_color)

        plt.show()

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
            save_fn = os.path.join(save_fn, '_L2.' +
                                   file_format)
            self.fig.savefig(save_fn, dpi=fig_dpi, format=file_format,
                             orientation=orientation, bbox_inches='tight')

        if close_fig == 'y':
            plt.clf()
            plt.close(self.fig)

        else:
            pass

        self.fig_fn = save_fn
        print('Saved figure to: ' + self.fig_fn)

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

        return "Plots RMS vs Iteration computed by Occam2D"


def parse_arguments(arguments):
    """
    takes list of command line arguments obtained by passing in sys.argv
    reads these and returns a parser object

    author: Alison Kirkby (2016)
    """

    import argparse

    parser = argparse.ArgumentParser(description='Set up and run a set of isotropic occam1d model runs')

    parser.add_argument('edipath',
                        help='folder containing edi files to use, full path or relative to working directory',
                        type=str)
    parser.add_argument('-l', '--program_location',
                        help='path to the inversion program',
                        type=str, default=r'/home/547/alk547/occam1d/OCCAM1DCSEM')
    parser.add_argument('-efr', '--resistivity_errorfloor',
                        help='error floor in resistivity, percent',
                        type=float, default=0)
    parser.add_argument('-efp', '--phase_errorfloor',
                        help='error floor in phase, degrees',
                        type=float, default=0)
    parser.add_argument('-efz', '--z_errorfloor',
                        help='error floor in z, percent',
                        type=float, default=0)
    parser.add_argument('-wd', '--working_directory',
                        help='working directory',
                        type=str, default='.')
    parser.add_argument('-m', '--modes', nargs='*',
                        help='modes to run, any or all of TE, TM, det (determinant)',
                        type=str, default=['TE'])
    parser.add_argument('-r', '--rotation_angle',
                        help='angle to rotate the data by, in degrees or can define option "strike" to rotate to strike, or "file" to get rotation angle from file',
                        type=str, default='0')
    parser.add_argument('-rfile', '--rotation_angle_file',
                        help='file containing rotation angles, first column is station name (must match edis) second column is rotation angle',
                        type=str, default=None)
    parser.add_argument('-spr', '--strike_period_range', nargs=2,
                        help='period range to use for calculation of strike if rotating to strike, two floats',
                        type=float, default=[1e-3, 1e3])
    parser.add_argument('-sapp', '--strike_approx',
                        help='approximate strike angle, the strike closest to this value is chosen',
                        type=float, default=0.)
    parser.add_argument('-q', '--remove_outofquadrant',
                        help='whether or not to remove points outside of the first or third quadrant, True or False',
                        type=bool, default=True)
    parser.add_argument('-itermax', '--iteration_max',
                        help='maximum number of iterations',
                        type=int, default=100)
    parser.add_argument('-rf', '--rms_factor',
                        help='factor to multiply the minimum possible rms by to get the target rms for the second run',
                        type=float, default=1.05)
    parser.add_argument('-rmsmin','--rms_min',
                        help='minimum target rms to assign, e.g. set a value of 1.0 to prevent overfitting data',
                        type=float, default=1.0)
    parser.add_argument('-nl', '--n_layers',
                        help='number of layers in the inversion',
                        type=int, default=80)
    parser.add_argument('-z1', '--z1_layer',
                        help='thickness of z1 layer',
                        type=float, default=10)
    parser.add_argument('-td', '--target_depth',
                        help='target depth for the inversion in metres',
                        type=int, default=10000)
    parser.add_argument('-rho0', '--start_rho',
                        help='starting resistivity value for the inversion',
                        type=float, default=100)
    parser.add_argument('-s', '--master_savepath',
                        help='master directory to save suite of runs into',
                        default='inversion_suite')

    args = parser.parse_args(arguments)
    args.working_directory = os.path.abspath(args.working_directory)
    if args.rotation_angle not in ['file', 'strike']:
        try:
            args.rotation_angle = float(args.rotation_angle)
        except:
            args.rotation_angle = 0.

    return args


def update_inputs():
    """
    update input parameters from command line

    author: Alison Kirkby (2016)
    """
    from sys import argv

    args = parse_arguments(argv[1:])
    cline_inputs = {}
    cline_keys = [i for i in dir(args) if i[0] != '_']

    for key in cline_keys:
        cline_inputs[key] = getattr(args, key)

    return cline_inputs


def get_strike(mt_object, fmin, fmax, strike_approx=0):
    """
    get the strike from the z array, choosing the strike angle that is closest
    to the azimuth of the PT ellipse (PT strike).

    if there is not strike available from the z array use the PT strike.

    """
    fselect = (mt_object.Z.freq > fmin) & (mt_object.Z.freq < fmax)

    # get median strike angles for frequencies needed (two strike angles due to 90 degree ambiguity)
    zstrike = mtg.strike_angle(z_object=mt_object.Z)[fselect]
    # put both strikes in the same quadrant for averaging
    zstrike = zstrike % 90
    zstrike = np.median(zstrike[np.isfinite(zstrike[:, 0])], axis=0)
    # add 90 to put one back in the other quadrant
    zstrike[1] += 90
    # choose closest value to approx_strike
    zstrike = zstrike[np.abs(zstrike - strike_approx) - np.amin(np.abs(zstrike - strike_approx)) < 1e-3]

    if len(zstrike) > 0:
        strike = zstrike[0]
    else:
        # if the data are 1d set strike to 90 degrees (i.e. no rotation)
        strike = 90.

    return strike


def generate_inputfiles(**input_parameters):
    """
    generate all the input files to run occam1d, return the path and the
    startup files to run.

    author: Alison Kirkby (2016)
    """
    edipath = op.join(input_parameters['working_directory'], input_parameters['edipath'])
    edilist = [ff for ff in os.listdir(edipath) if ff.endswith('.edi')]

    wkdir_master = op.join(input_parameters['working_directory'],
                           input_parameters['master_savepath'])
    if not os.path.exists(wkdir_master):
        os.mkdir(wkdir_master)

    rundirs = {}

    for edifile in edilist:
        # read the edi file to get the station name
        eo = mt.MT(op.join(edipath, edifile))
        #print(input_parameters['rotation_angle'], input_parameters['working_directory'], input_parameters[
        #    'rotation_angle_file'])
        if input_parameters['rotation_angle'] == 'strike':
            spr = input_parameters['strike_period_range']
            fmax, fmin = [1. / np.amin(spr), 1. / np.amax(spr)]
            rotangle = (get_strike(eo, fmin, fmax,
                                   strike_approx=input_parameters['strike_approx']) - 90.) % 180
        elif input_parameters['rotation_angle'] == 'file':
            with open(op.join(input_parameters['working_directory'], input_parameters['rotation_angle_file'])) as f:
                line = f.readline().strip().split()

                while string.upper(line[0]) != string.upper(eo.station):
                    line = f.readline().strip().split()
                    if len(line) == 0:
                        line = ['', '0.0']
                        break
            rotangle = float(line[1])
        else:
            rotangle = input_parameters['rotation_angle']
            
        # create a working directory to store the inversion files in
        svpath = 'station' + eo.station
        wd = op.join(wkdir_master, svpath)
        if not os.path.exists(wd):
            os.mkdir(wd)
        rundirs[svpath] = []

        # create the model file
        ocm = Model(n_layers=input_parameters['n_layers'], save_path=wd,
                    target_depth=input_parameters['target_depth'],
                    z1_layer=input_parameters['z1_layer'])
        ocm.write_model_file()

        for mode in input_parameters['modes']:
            # create a data file for each mode
            ocd = Data()
            ocd._data_fn = 'Occam1d_DataFile_rot%03i' % rotangle
            ocd.write_data_file(
                res_errorfloor=input_parameters['resistivity_errorfloor'],
                phase_errorfloor=input_parameters['phase_errorfloor'],
                z_errorfloor=input_parameters['z_errorfloor'],
                remove_outofquadrant=input_parameters['remove_outofquadrant'],
                mode=mode,
                edi_file=op.join(edipath, edifile),
                thetar=rotangle,
                save_path=wd)

            ocs = Startup(data_fn=ocd.data_fn,
                          model_fn=ocm.model_fn,
                          start_rho=input_parameters['start_rho'])
            startup_fn = 'OccamStartup1D' + mode
            ocs.write_startup_file(save_path=wd,
                                   startup_fn=op.join(wd, startup_fn),
                                   max_iter=input_parameters['iteration_max'],
                                   target_rms=input_parameters['rms_min']/input_parameters['rms_factor'])
            rundirs[svpath].append(startup_fn)

    return wkdir_master, rundirs


def divide_inputs(work_to_do, size):
    """
    divide list of inputs into chunks to send to each processor

    """
    chunks = [[] for _ in range(size)]
    for i, d in enumerate(work_to_do):
        chunks[i % size].append(d)

    return chunks


def build_run():
    """
    build input files and run a suite of models in series (pretty quick so won't bother parallelise)

    run Occam1d on each set of inputs.
    Occam is run twice. First to get the lowest possible misfit.
    we then set the target rms to a factor (default 1.05) times the minimum rms achieved
    and run to get the smoothest model.

    author: Alison Kirkby (2016)
    """
    # from mpi4py import MPI

    # get command line arguments as a dictionary
    input_parameters = update_inputs()

    # create the inputs and get the run directories
    master_wkdir, run_directories = generate_inputfiles(**input_parameters)

    # run Occam1d on each set of inputs.
    # Occam is run twice. First to get the lowest possible misfit.
    # we then set the target rms to a factor (default 1.05) times the minimum rms achieved
    # and run to get the smoothest model.
    for rundir in list(run_directories.keys()):
        wd = op.join(master_wkdir, rundir)
        os.chdir(wd)
        for startupfile in run_directories[rundir]:
            # define some parameters
            mode = startupfile[14:]
            iterstring = 'RMSmin' + mode
            # run for minimum rms
            subprocess.call([input_parameters['program_location'],
                             startupfile,
                             iterstring])
            # read the iter file to get minimum rms
            iterfilelist = [ff for ff in os.listdir(wd) if (ff.startswith(iterstring) and ff.endswith('.iter'))]
            # only run a second lot of inversions if the first produced outputs
            if len(iterfilelist) > 0:
                iterfile = max(iterfilelist)
                startup = Startup()
                startup.read_startup_file(op.join(wd, iterfile))
                # create a new startup file the same as the previous one but target rms is factor*minimum_rms
                target_rms = float(startup.misfit_value) * input_parameters['rms_factor']
                if target_rms < input_parameters['rms_min']:
                    target_rms = input_parameters['rms_min']
                startupnew = Startup(data_fn=op.join(wd, startup.data_file),
                                     model_fn=op.join(wd, startup.model_file),
                                     max_iter=input_parameters['iteration_max'],
                                     start_rho=input_parameters['start_rho'],
                                     target_rms=target_rms)
                startupnew.write_startup_file(startup_fn=op.join(wd, startupfile), save_path=wd)
                # run occam again
                subprocess.call([input_parameters['program_location'],
                                 startupfile,
                                 'Smooth' + mode])


if __name__ == '__main__':
    build_run()
