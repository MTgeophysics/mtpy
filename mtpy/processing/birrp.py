# -*- coding: utf-8 -*-
"""
BIRRP
===========
    * deals with inputs and outputs from BIRRP

Created on Tue Sep 20 14:33:20 2016

@author: jrpeacock
"""

#==============================================================================
import numpy as np
import os
import subprocess
from datetime import datetime

import mtpy.utils.configfile as mtcfg
import mtpy.utils.filehandling as mtfh
import mtpy.utils.exceptions as mtex
import mtpy.core.mt as mt

#==============================================================================
class BIRRPParameters(object):
    """
    class to hold and produce the appropriate parameters given the input
    parameters.
    """

    def __init__(self, ilev=0, **kwargs):
        # set the input level
        self.ilev = ilev

        self.ninp = 2
        self._nout = 3
        self._nref = 2
        self.nr2 = 0
        self.nr3 = 0
        self.tbw = 2.0
        self.nfft = 2**18
        self.nsctmax = 14
        self.ofil = 'mt'
        self.nlev = 0
        self.nar = 5
        self.imode = 0
        self.jmode = 0
        self.nfil = 0
        self.nrr = 1
        self.nsctinc = 2
        self.nsctmax = int(np.floor(np.log2(self.nfft))-4)
        self.nf1 = int(self.tbw+2)
        self.nfinc = int(self.tbw)
        self.nfsect = 2
        self.mfft = 2
        self.uin = 0
        self.ainlin = 0.0001
        self.ainuin = 0.9999
        self.c2threshb = 0.0
        self.c2threshe = 0.0
        self.nz = 0
        self.perlo = 1000
        self.perhi = .0001
        self.nprej = 0
        self.prej = None
        self.c2threshe1 = 0

        self.thetae = [0, 90, 0]
        self.thetab = [0, 90, 0]
        self.thetaf = [0, 90, 0]

        # set any attributes that are given
        for key in list(kwargs.keys()):
            setattr(self, key, kwargs[key])

        self._validate_parameters()

    def to_dict(self):
        """
        get appropriate parameters
        """

        param_dict = {}
        param_dict['ninp'] = self.ninp
        param_dict['nout'] = self._nout
        param_dict['nref'] = self._nref
        param_dict['tbw'] = self.tbw
        param_dict['nfft'] = self.nfft
        param_dict['nsctmax'] = self.nsctmax
        param_dict['ofil'] = self.ofil
        param_dict['nlev'] = self.nlev
        param_dict['nar'] = self.nar
        param_dict['imode'] = self.imode
        param_dict['jmode'] = self.jmode
        param_dict['nfil'] = self.nfil
        param_dict['thetae'] = self.thetae
        param_dict['thetab'] = self.thetab
        param_dict['thetaf'] = self.thetaf

        if self.ilev == 0:

            param_dict['uin'] = self.uin
            param_dict['ainuin'] = self.ainuin
            param_dict['c2threshe'] = self.c2threshe
            param_dict['nz'] = self.nz
            param_dict['c2thresh1'] = self.c2threshe1

        elif self.ilev == 1:
            if self._nref > 3:
                param_dict['nref2'] = self._nref_2
                param_dict['nref3'] = self._nref_3
            param_dict['nrr'] = self.nrr
            param_dict['nsctinc'] = self.nsctinc
            param_dict['nsctmax'] = self.nsctmax
            param_dict['nf1'] = self.nf1
            param_dict['nfinc'] = self.nfinc
            param_dict['nfsect'] = self.nfsect
            param_dict['uin'] = self.uin
            param_dict['ainlin'] = self.ainlin
            param_dict['ainuin'] = self.ainuin
            if self.nrr == 1:
                param_dict['c2thresh'] = self.c2threshb
                param_dict['c2threse'] = self.c2threshe
                if self.c2threshe == 0 and self.c2threshb == 0:
                    param_dict['nz'] = self.nz
                else:
                    param_dict['nz'] = self.nz
                    param_dict['perlo'] = self.perlo
                    param_dict['perhi'] = self.perhi
            elif self.nrr == 0:
                param_dict['c2threshb'] = self.c2threshb
                param_dict['c2threse'] = self.c2threshe
            param_dict['nprej'] = self.nprej
            param_dict['prej'] = self.prej
            param_dict['c2thresh1'] = self.c2threshe1

        return param_dict

    def from_dict(self, birrp_dict):
        """
        set birrp parameters from dict
        """

        for key, value in birrp_dict.items():
            try:
                setattr(self, key, value)
            except AttributeError:
                print('WARNING: cannot set {0}, skipping'.format(key))

        self._validate_parameters()

    def _validate_parameters(self):
        """
        check to make sure the parameters are legit.
        """

        # be sure the
        if self.ninp not in [1, 2, 3]:
            print('WARNING: Number of inputs {0} not allowed.'.format(self.ninp))
            self.ninp = 2
            print('  --> setting ninp to {0}'.format(self.ninp))

        if self._nout not in [2, 3]:
            print('WARNING: Number of outputs {0} not allowed.'.format(self._nout))
            self._nout = 2
            print('  --> setting nout to {0}'.format(self._nout))

        if self._nref > 3:
            print('WARNING: nref > 3, setting ilev to 1')
            self.ilev = 1

        if self.tbw < 0 or self.tbw > 4:
            print('WARNING: Total bandwidth of slepian window {0} not allowed.'.format(self.tbw))
            self.tbw = 2
            print('  --> setting tbw to {0}'.format(self.tbw))

        if (np.log2(self.nfft) % 1) != 0.0:
            print('WARNING: Window length nfft should be a power of 2 not {0}'.format(self.nfft))
            self.nfft = 2**np.floor(np.log2(self.nfft))
            print('  -- > setting nfft to {0}, (2**{1:.0f})'.format(self.nfft,
                                                              np.log2(self.nfft)))

        if np.log2(self.nfft)-self.nsctmax < 4:
            print('WARNING: Maximum number of windows {0} is too high'.format(self.nsctmax))
            self.nsctmax = np.log2(self.nfft)-4
            print('  --> setting nsctmax to {0}'.format(self.nsctmax))

        if self.uin != 0:
            print('WARNING: You\'re playing with fire if uin is not 0.')
            self.uin = 0
            print('  --> setting uin to 0, if you don\'t want that change it back')

        if self.imode not in [0, 1, 2, 3]:
            raise BIRRPParameterError('Invalid number for time series mode,'
                                       'imode, {0}, should be 0, 1, 2, or 3'.format(self.imode))

        if self.jmode not in [0, 1]:
            raise BIRRPParameterError('Invalid number for time mode,'
                                       'imode, {0}, should be 0, or 1'.format(self.imode))
        if self.ilev == 1:
            if self.nsctinc != 2:
                print('WARNING: Check decimation increment nsctinc, should be 2 not {0}'.format(self.nsctinc))

            if self.nfsect != 2:
                print('WARNING: Will get an error from BIRRP if nfsect is not 2.')
                print('number of frequencies per section is {0}'.format(self.nfsect))
                self.nfsect = 2
                print('  --> setting nfsect to 2')

            if self.nf1 != self.tbw + 2:
                print('WARNING: First frequency should be around tbw+2.')
                print('nf1 currently set to {0}'.format(self.nf1))

            if self.nfinc != self.tbw:
                print('WARNING: sequence of frequencies per window should be around tbw.')
                print('nfinc currently set to {0}'.format(self.nfinc))

            if self.nprej != 0:
                if self.prej is None or type(self.prej) is not list:
                    raise BIRRPParameterError('Need to input a prejudice list if nprej != 0'+
                                     '\nInput as a list of frequencies' )

            if self.nrr not in [0, 1]:
                print(('WARNING: Value for picking remote reference or '+
                     'two stage processing, nrr, '+
                     'should be 0 or 1 not {0}'.format(self.nrr)))
                self.nrr = 0
                print('  --> setting nrr to {0}'.format(self.nrr))


            if self.c2threshe != 0 or self.c2threshb != 0:
                if not self.perhi:
                    raise BIRRPParameterError('Need to input a high period (s) threshold as perhi')

                if not self.perlo:
                    raise BIRRPParameterError('Need to input a low period (s) threshold as perlo')

        if len(self.thetae) != 3:
            print('WARNING: Electric rotation angles not input properly {0}'.format(self.thetae))
            print('input as north, east, orthogonal rotation')
            self.thetae = [0, 90, 0]
            print('  --> setting thetae to {0}'.format(self.thetae))

        if len(self.thetab) != 3:
            print('WARNING: Magnetic rotation angles not input properly {0}'.format(self.thetab))
            print('input as north, east, orthogonal rotation')
            self.thetab = [0, 90, 0]
            print('  --> setting thetab to {0}'.format(self.thetab))


        if len(self.thetaf) != 3:
            print('WARNING: Field rotation angles not input properly {0}'.format(self.thetaf))
            print('WARNING: input as north, east, orthogonal rotation')
            self.thetaf = [0, 90, 0]
            print('  --> setting thetaf to {0}'.format(self.thetaf))


    def read_config_file(self, birrp_config_fn):
        """
        read in a configuration file and fill in the appropriate parameters
        """

        birrp_dict = mtcfg.read_configfile(birrp_config_fn)

        for birrp_key in list(birrp_dict.keys()):
            try:
                b_value = float(birrp_dict[birrp_key])
                if np.remainder(b_value, 1) == 0:
                    b_value = int(b_value)
            except ValueError:
                b_value = birrp_dict[birrp_key]

            setattr(self, birrp_key, b_value)

    def write_config_file(self, save_fn):
        """
        write a config file for birrp parameters
        """

        cfg_fn = mtfh.make_unique_filename('{0}_birrp_params.cfg'.format(save_fn))

        station = os.path.basename(save_fn)
        birrp_dict = self.to_dict()
        mtcfg.write_dict_to_configfile({station:birrp_dict}, cfg_fn)
        print('INFO: Wrote BIRRP config file for edi file to {0}'.format(cfg_fn))
#==============================================================================
# Error classes
#==============================================================================
class BIRRPParameterError(Exception):
    pass

class ScriptFileError(Exception):
    pass

#==============================================================================
# write script file
#==============================================================================
class ScriptFile(BIRRPParameters):
    """
    class to read and write script file

    **Arguments**
    --------------------

        **fn_arr** : numpy.ndarray
                     numpy.ndarray([[block 1], [block 2]])

    .. note::  [block n] is a numpy structured array with data type

        =================== ================================= =================
        Name                Description                       Type
        =================== ================================= =================
        fn                  file path/name                    string
        nread               number of points to read          int
        nskip               number of points to skip          int
        comp                component                         [ ex | ey | hx hy | hz ]
        calibration_fn      calibration file path/name        string
        rr                  a remote reference channel        [ True | False ]
        rr_num              remote reference pair number      int
        start               start time iso format             Timestamp
        stop                stop time iso format              Timestamp
        station             station name                      string
        sampling_rate       sampling rate                     int
        =================== ================================= =================


    **BIRRP Parameters**
    -------------------------

    ================== ========================================================
    parameter          description
    ================== ========================================================
    ilev               processing mode 0 for basic and 1 for advanced RR-2
                       stage
    nout               Number of Output time series (2 or 3-> for BZ)
    ninp               Number of input time series for E-field (1,2,3)
    nref               Number of reference channels (2 for MT)
    nrr                bounded remote reference (0) or 2 stage bounded
                       influence (1)
    tbw                Time bandwidth for Sepian sequence
    deltat             Sampling rate (+) for (s), (-) for (Hz)
    nfft               Length of FFT (should be even)
    nsctinc            section increment divisor (2 to divide by half)
    nsctmax            Number of windows used in FFT
    nf1                1st frequency to extract from FFT window (>=3)
    nfinc              frequency extraction increment
    nfsect             number of frequencies to extract
    mfft               AR filter factor, window divisor (2 for half)
    uin                Quantile factor determination
    ainlin             Residual rejection factor low end (usually 0)
    ainuin             Residual rejection factor high end (.95-.99)
    c2threshb          Coherence threshold for magnetics (0 if undesired)
    c2threshe          Coherence threshold for electrics (0 if undesired)
    nz                 Threshold for Bz (0=separate from E, 1=E threshold,
                                         2=E and B)
                       Input if 3 B components else None
    c2thresh1          Squared coherence for Bz, input if NZ=0, Nout=3
    perlo              longest period to apply coherence threshold over
    perhi              shortes period to apply coherence threshold over
    ofil               Output file root(usually three letters, can add full
                                        path)
    nlev               Output files (0=Z; 1=Z,qq; 2=Z,qq,w; 3=Z,qq,w,d)
    nprej              number of frequencies to reject
    prej               frequencies to reject (+) for period, (-) for frequency
    npcs               Number of independent data to be processed (1 for one
                       segement)
    nar                Prewhitening Filter (3< >15) or 0 if not desired',
    imode              Output file mode (0=ascii; 1=binary; 2=headerless ascii;
                       3=ascii in TS mode',
    jmode              input file mode (0=user defined; 1=sconvert2tart time
                                        YYYY-MM-DD HH:MM:SS)',
    nread              Number of points to be read for each data set
                       (if segments>1 -> npts1,npts2...)',
    nfil               Filter parameters (0=none; >0=input parameters;
                                          <0=filename)
    nskip              Skip number of points in time series (0) if no skip,
                        (if segements >1 -> input1,input2...)',
    nskipr             Number of points to skip over (0) if none,
                       (if segements >1 -> input1,input2...)',
    thetae             Rotation angles for electrics (relative to geomagnetic
                       North)(N,E,rot)',
    thetab             Rotation angles for magnetics (relative to geomagnetic
                       North)(N,E,rot)',
    thetar             Rotation angles for calculation (relative to geomagnetic
                       North)(N,E,rot)'
    ================== ========================================================

    .. note:: Currently only supports jmode = 0 and imode = 0

    .. seealso:: BIRRP Manual and publications by Chave and Thomson
                for more details on the parameters found at:

                http://www.whoi.edu/science/AOPE/people/achave/Site/Next1.html

    """

    def __init__(self, script_fn=None, fn_arr=None, **kwargs):
        super(ScriptFile, self).__init__(fn_arr=None, **kwargs)

        self.fn_arr = fn_arr
        self.script_fn = None

        self._npcs = 0
        self._nref = 2
        self._nref_2 = 0
        self._nref_3 = 0
        self._comp_list = None

        self._fn_dtype = np.dtype([('fn', 'U100'),
                                   ('nread', np.int),
                                   ('nskip', np.int),
                                   ('comp', 'U2'),
                                   ('calibration_fn', 'U100'),
                                   ('rr', np.bool),
                                   ('rr_num', np.int),
                                   ('start', 'U19'),
                                   ('stop', 'U19'),
                                   ('sampling_rate', np.int),
                                   ('station', 'U10')])

        if self.fn_arr is not None:
            self._validate_fn_arr()

        for key in list(kwargs.keys()):
            setattr(self, key, kwargs[key])

    def _validate_fn_arr(self):
        """
        make sure fn_arr is an np.array
        """

        for aa in self.fn_arr:
            if type(aa) not in [np.ndarray, np.recarray]:
                raise ScriptFileError('Input fn_arr elements should be numpy'+
                                      ' arrays or recarray' +
                                      ' with dtype {0}'.format(self._fn_dtype))

        names = list(self.fn_arr[0].dtype.names)
        if 'index' in names:
            names.remove('index')
        if not sorted(names) == sorted(self._fn_dtype.names):
            raise ScriptFileError('fn_arr.dtype needs to be {0}'.format(self._fn_dtype))

        # make sure the shapes are the same
        shapes = [aa.shape[0] for aa in self.fn_arr]
        if min(shapes) != max(shapes):
            raise ScriptFileError('fn_arr does not have all the same shapes.'+
                                  '{0}'.format(shapes))

        # make sure that rr is bool
        for aa in self.fn_arr:
            if aa['rr'].dtype != np.dtype(bool):
                for index, test in enumerate(aa['rr']):
                    if test in ['true', 'True', 'TRUE', True]:
                        aa['rr'][index] = True
                    else:
                        aa['rr'][index] = False
                aa = aa['rr'].astype(bool)

        # make sure all the same sampling rate
        sr = np.unique(self.fn_arr['sampling_rate'])
        if len(sr) != 1:
            raise ScriptFileError('Samping rates are not the same, found'+
                                  ' {0}'.format(sr))

        # make sure components are named correctly
        for aa in self.fn_arr:
            for element in aa:
                if element['rr']:
                    if element['comp'].find('_') < 0:
                        element['comp'] = 'rr{0}_{1:02}'.format(element['comp'],
                                                                element['rr_num'])


    @property
    def nout(self):
        if self.fn_arr is not None:
            self._nout = len(np.where(self.fn_arr[0]['rr'] == False)[0])-2
        else:
            print('WARNING: fn_arr is None, set nout to 0')
            self._nout = 0
        return self._nout

    @property
    def npcs(self):
        if self.fn_arr is not None:
            self._npcs = len(self.fn_arr)
        else:
            print('WARNING: fn_arr is None, set npcs to 0')
            self._npcs = 0
        return self._npcs

    @property
    def nref(self):
        if self.fn_arr is not None:
            num_ref = np.where(self.fn_arr[0]['rr'] == True)[0]
            self._nref = len(num_ref)
        else:
            print('WARNING: fn_arr is None, set nref to 0')
            self._nref = 0

        if self._nref > 3:
            self.nr2 = self.fn_arr[0]['rr_num'].max()

        return self._nref

    @property
    def deltat(self):
        return -1 * np.unique(self.fn_arr['sampling_rate'])[0]

    @property
    def comp_list(self):
        num_comp = self.ninp + self.nout
        if num_comp == 4:
            self._comp_list = ['ex', 'ey', 'hx', 'hy']
        elif num_comp == 5:
            self._comp_list = ['ex', 'ey', 'hz', 'hx', 'hy']
        else:
            raise ValueError('Number of components {0} invalid, check inputs'.format(num_comp))

        if self.nref == 0:
            self._comp_list += ['hx', 'hy']

        else:
            for ii in range(int(self.nref/2)):
                self._comp_list += ['rrhx_{0:02}'.format(ii),
                                    'rrhy_{0:02}'.format(ii)]

        return self._comp_list


    def write_script_file(self, script_fn=None, ofil=None):
        if ofil is not None:
            self.ofil = ofil

        if script_fn is not None:
            self.script_fn = script_fn

        # be sure all the parameters are correct according to BIRRP
        self.nout
        self.nref
        self.npcs
        self.comp_list
        self._validate_parameters()

        # begin writing script file
        s_lines = []
        s_lines += ['{0:0.0f}'.format(self.ilev)]
        s_lines += ['{0:0.0f}'.format(self.nout)]
        s_lines += ['{0:0.0f}'.format(self.ninp)]

        if self.ilev == 0:

            s_lines += ['{0:.3f}'.format(self.tbw)]
            s_lines += ['{0:.3f}'.format(self.deltat)]
            s_lines += ['{0:0.0f},{1:0.0f}'.format(self.nfft, self.nsctmax)]
            s_lines += ['y']
            s_lines += ['{0:.5f},{1:.5f}'.format(self.uin, self.ainuin)]
            s_lines += ['{0:.3f}'.format(self.c2threshe)]
            #parameters for bz component if ninp=3
            if self.nout == 3:
                if self.c2threshe == 0:
                    s_lines += ['{0:0.0f}'.format(0)]
                    s_lines += ['{0:.3f}'.format(self.c2threshe1)]
                else:
                    s_lines += ['{0:0.0f}'.format(self.nz)]
                    s_lines += ['{0:.3f}'.format(self.c2threshe1)]
            else:
                pass
            s_lines += [self.ofil]
            s_lines += ['{0:0.0f}'.format(self.nlev)]

        elif self.ilev == 1:
            print('INFO: Writing Advanced mode')
            s_lines += ['{0:0.0f}'.format(self.nref)]
            if self.nref > 3:
                s_lines += ['{0:0.0f},{1:0.0f}'.format(self.nr3, self.nr2)]
            s_lines += ['{0:0.0f}'.format(self.nrr)]
            s_lines += ['{0:.3f}'.format(self.tbw)]
            s_lines += ['{0:.3f}'.format(self.deltat)]
            s_lines += ['{0:0.0f},{1:.2g},{2:0.0f}'.format(self.nfft,
                                                           self.nsctinc,
                                                           self.nsctmax)]
            s_lines += ['{0:0.0f},{1:.2g},{2:0.0f}'.format(self.nf1,
                                                           self.nfinc,
                                                           self.nfsect)]
            s_lines += ['y']
            s_lines += ['{0:.2g}'.format(self.mfft)]
            s_lines += ['{0:.5g},{1:.5g},{2:.5g}'.format(self.uin,
                                                            self.ainlin,
                                                            self.ainuin)]
            #if remote referencing
            if int(self.nrr) == 0:

                s_lines += ['{0:.3f}'.format(self.c2threshe)]
                #parameters for bz component if ninp=3
                if self.nout == 3:
                    if self.c2threshe != 0:
                        s_lines += ['{0:0.0f}'.format(self.nz)]
                        s_lines += ['{0:.3f}'.format(self.c2threshe1)]
                    else:
                        s_lines += ['{0:0.0f}'.format(0)]
                        s_lines += ['{0:.3f}'.format(self.c2threshe1)]
                    if self.c2threshe1 != 0.0 or self.c2threshe != 0.0:
                        s_lines += ['{0:.6g},{1:.6g}'.format(self.perlo,
                                                                self.perhi)]
                else:
                    if self.c2threshe != 0.0:
                        s_lines += ['{0:.6g},{1:.6g}'.format(self.perlo,
                                                                self.perhi)]
            #if 2 stage processing
            elif int(self.nrr) == 1:
                s_lines += ['{0:.3f}'.format(self.c2threshb)]
                s_lines += ['{0:.3f}'.format(self.c2threshe)]
                if self.nout == 3:
                    if self.c2threshb != 0 or self.c2threshe != 0:
                        s_lines += ['{0:0.0f}'.format(self.nz)]
                        s_lines += ['{0:.3f}'.format(self.c2threshe1)]
                    elif self.c2threshb == 0 and self.c2threshe == 0:
                        s_lines += ['{0:0.0f}'.format(0)]
                        s_lines += ['{0:.3f}'.format(0)]
                if self.c2threshb != 0.0 or self.c2threshe != 0.0:
                    s_lines += ['{0:.6g},{1:.6g}'.format(self.perlo,
                                                            self.perhi)]
            s_lines += [self.ofil]
            s_lines += ['{0:0.0f}'.format(self.nlev)]
            s_lines += ['{0:0.0f}'.format(self.nprej)]
            if self.nprej != 0:
                if type(self.prej) is not list:
                    self.prej = [self.prej]
                s_lines += ['{0:.5g}'.format(nn) for nn in self.prej]

        s_lines += ['{0:0.0f}'.format(self.npcs)]
        s_lines += ['{0:0.0f}'.format(self.nar)]
        s_lines += ['{0:0.0f}'.format(self.imode)]
        s_lines += ['{0:0.0f}'.format(self.jmode)]

       #write in filenames
        if self.jmode == 0:
            # loop over each data block
            for ff, fn_arr in enumerate(self.fn_arr):

                # get the least amount of data points to read
                s_lines += ['{0:0.0f}'.format(fn_arr['nread'].min())]

                for cc in self.comp_list:
                    fn_index = np.where(fn_arr['comp'] == cc)[0][0]

                    if ff == 0:
                        fn_lines = self.make_fn_lines_block_00(fn_arr[fn_index])
                    else:
                        fn_lines = self.make_fn_lines_block_n(fn_arr[fn_index])
                    s_lines += fn_lines

        #write rotation angles
        s_lines += [' '.join(['{0:.2f}'.format(theta)
                              for theta in self.thetae])]
        s_lines += [' '.join(['{0:.2f}'.format(theta)
                              for theta in self.thetab])]
        s_lines += [' '.join(['{0:.2f}'.format(theta)
                              for theta in self.thetaf])]

        if self.nref > 3:
            for kk in range(self.nref):
                s_lines += [' '.join(['{0:.2f}'.format(theta)
                                      for theta in self.thetab])]

        with open(self.script_fn, 'w') as fid:
            try:
                fid.write('\n'.join(s_lines))
            except TypeError:
                print(s_lines)

        print('INFO: Wrote script file to {0}'.format(self.script_fn))

    def make_fn_lines_block_00(self, fn_arr):
        """
        make lines for file in script file which includes

            - nread
            - filter_fn
            - fn
            - nskip

        """
        lines = []
        if fn_arr['calibration_fn'] in ['', 0, '0']:
            lines += ['0']
        else:
            lines += ['-2']
            lines += [str(fn_arr['calibration_fn'])]
        lines += [str(fn_arr['fn'])]
        lines += ['{0:d}'.format(fn_arr['nskip'])]

        return lines

    def make_fn_lines_block_n(self, fn_arr):
        """
        make lines for file in script file which includes

            - nread
            - filter_fn
            - fn
            - nskip

        """
        lines = []
        lines += [str(fn_arr['fn'])]
        lines += ['{0:d}'.format(fn_arr['nskip'])]

        return lines

# =============================================================================
# run birrp
# =============================================================================
def run(birrp_exe, script_file):
    """
    run a birrp script file from command line via python subprocess.

    Arguments
    --------------

        **birrp_exe** : string
                        full path to the compiled birrp executable

        **script_file** : string
                          full path to input script file following the
                          guidelines of the BIRRP documentation.

    Outputs
    ---------------

        **log_file.log** : a log file of how BIRRP ran


    .. seealso:: BIRRP Manual and publications by Chave and Thomson
                for more details on the parameters found at:

                http://www.whoi.edu/science/AOPE/people/achave/Site/Next1.html

    """
    # check to make sure the given executable is legit
    if not os.path.isfile(birrp_exe):
        raise mtex.MTpyError_inputarguments('birrp executable not found:'+
                                            '{0}'.format(birrp_exe))

    # get the current working directory so we can go back to it later
    current_dir = os.path.abspath(os.curdir)

    #change directory to directory of the script file
    os.chdir(os.path.dirname(script_file))
    local_script_fn = os.path.basename(script_file)
    st = datetime.now()

    print('*'*10)
    print('INFO: Processing {0} with {1}'.format(script_file, birrp_exe))
    print('INFO: Starting Birrp processing at {0}...'.format(st))

    birrp_process = subprocess.Popen(birrp_exe+'< {0}'.format(local_script_fn),
                                     stdin=subprocess.PIPE,
                                     shell=True)
    birrp_process.wait()
    et = datetime.now()
    print('_'*60)
    print('INFO: Starting Birrp processing at {0}...'.format(st))
    print('INFO: Ended Birrp processing at   {0}...'.format(et))

    #go back to initial directory
    os.chdir(current_dir)
    t_diff = (et - st).total_seconds()
    print('\n{0} DONE !!! {0}'.format('='*20))
    print('\tTook {0:02}:{1:02} minutes:seconds'.format(int(t_diff // 60),
                                                        int(t_diff % 60)))
    
    return birrp_process.communicate()[0]

#==============================================================================
# Write edi file from birrp outputs
#==============================================================================
class J2Edi(object):
    """
    Read in BIRRP out puts, in this case the .j file and convert that into
    an .edi file using the survey_config_fn parameters.

    Key Word Arguments
    ----------------------

        **birrp_dir** : string
                        full path to directory where birrp outputs are

        **station** : string
                      station name

        **survey_config_fn** : string
                               full path to survey configuration file with
                               information on location and site setup
                               must have a key that is the same as station.

        **birrp_config_fn** : string
                              full path to configuration file that was used to
                              process with (all the birrp parameters used).  If
                              None is input, the file is searched for, if it
                              is not found, the processing parameters are
                              used from the .j file.

        **j_fn** : string
                   full path to j file.  If none is input the .j file is
                   searched for in birrp_dir.


    ============================ ==============================================
    Methods                      Description
    ============================ ==============================================
    read_survey_config_fn        read in survey configuration file
    get_birrp_config_fn          get the birrp_config_fn in birrp_dir
    read_birrp_config_fn         read in birrp_config_fn
    get_j_file                   find .j file in birrp_dir
    write_edi_file               write an .edi file fro all the provided
                                 information.
    ============================ ==============================================

    Example
    -------------

        >>> import mtpy.proceessing.birrp as birrp
        >>> j2edi_obj = birrp.J_To_Edi()
        >>> j2edi_obj.birrp_dir = r"/home/data/mt01/BF/256"
        >>> j2edi_obj.station = 'mt01'
        >>> j2edi_obj.survey_config_fn = r"/home/data/2016_survey.cfg"
        >>> j2edi_obj.write_edi_file()

    """

    def __init__(self, **kwargs):

        self.birrp_dir = None
        self.birrp_config_fn = None
        self.birrp_dict = None

        self.survey_config_fn = None
        self.survey_config_dict = None

        self.station = None
        self.j_fn = None
        self.mt_obj = mt.MT()

        for key in list(kwargs.keys()):
            setattr(self, key, kwargs[key])

    def read_survey_config_fn(self, survey_config_fn=None):
        """
        read in survey configuration file and output into a useful dictionary
        """
        if survey_config_fn is not None:
            self.surve_config_fn = survey_config_fn

        if not os.path.isfile(self.survey_config_fn):
            raise mtex.MTpyError_inputarguments('Could not find {0}, check path'.format(survey_config_fn))

        # read in survey information
        self.survey_config_dict = mtcfg.read_survey_configfile(self.survey_config_fn)[self.station]

    def get_birrp_config_fn(self):
        """
        get birrp configuration file from birrp directory
        """

        if self.birrp_dir is None:
            print('WARNING: Could not get birrp_config_fn because no birrp directory specified')
            self.birrp_config_fn = None
            return

        try:
            self.birrp_config_fn = [os.path.join(self.birrp_dir, fn)
                                    for fn in os.listdir(self.birrp_dir)
                                    if fn.find('birrp_params') > 0][-1]
            print('INFO: Using {0}'.format(self.birrp_config_fn))

        except IndexError:
            print('WARNING: Could not find a birrp_params config file in {0}'.format(self.birrp_dir))
            self.birrp_config_fn = None
            return

    def read_birrp_config_fn(self, birrp_config_fn=None):
        """
        read in birrp configuration file
        """

        if birrp_config_fn is not None:
            self.birrp_config_fn = birrp_config_fn

        if self.birrp_config_fn is None:
            self.get_birrp_config_fn()
            if self.birrp_config_fn is None:
                self.birrp_dict = None
                return

        self.birrp_dict = mtcfg.read_configfile(self.birrp_config_fn)

    def get_j_file(self, birrp_dir=None):
        """
        get .j file output by birrp

        """

        if birrp_dir is not None:
            self.birrp_dir = birrp_dir

        if self.birrp_dir is None or not os.path.exists(self.birrp_dir):
            raise mtex.MTpyError_inputarguments('No birrp directory input,'
                                                'check path {0}'.format(self.birrp_dir))
        try:
            self.j_fn = [os.path.join(self.birrp_dir, fn)
                         for fn in os.listdir(self.birrp_dir)
                         if fn.endswith('.j')][0]
        except IndexError:
            print('ERROR: Could not find a .j file in {0}, check path.'.format(self.birrp_dir))
            self.j_fn = None

    def _fill_site(self):
        """
        fill header data
        """
        self.mt_obj.lat = self.survey_config_dict['latitude']
        self.mt_obj.lon = self.survey_config_dict['longitude']
        self.mt_obj.Site.start_date = self.survey_config_dict['date']
        self.mt_obj.Site.survey = self.survey_config_dict['location']
        self.mt_obj.station = self.survey_config_dict['station']
        self.mt_obj.elev = self.survey_config_dict['elevation']
        self.mt_obj.Site.acquired_by = self.survey_config_dict['network']

    def _fill_info(self):
        """
        fill information section
        """
        self.mt_obj.Notes.info_dict = {}

        for key in sorted(self.birrp_dict.keys()):
            setattr(self.mt_obj.Processing, 'birrp_'+key, self.birrp_dict[key])


    def _fill_field_notes(self):
        """
        get define measurement blocks
        """

        # --> hx
        self.mt_obj.FieldNotes.Magnetometer_hx.id = 1
        self.mt_obj.FieldNotes.Magnetometer_hx.chtype = 'hx'
        self.mt_obj.FieldNotes.Magnetometer_hx.x = 0
        self.mt_obj.FieldNotes.Magnetometer_hx.y = 0
        self.mt_obj.FieldNotes.Magnetometer_hx.azm = float(self.survey_config_dict['b_xaxis_azimuth'])
        self.mt_obj.FieldNotes.Magnetometer_hx.acqchan = self.survey_config_dict['hx']


        #--> hy
        self.mt_obj.FieldNotes.Magnetometer_hy.id = 2
        self.mt_obj.FieldNotes.Magnetometer_hy.chtype = 'hy'
        self.mt_obj.FieldNotes.Magnetometer_hy.x = 0
        self.mt_obj.FieldNotes.Magnetometer_hy.y = 0
        self.mt_obj.FieldNotes.Magnetometer_hy.azm = float(self.survey_config_dict['b_yaxis_azimuth'])
        self.mt_obj.FieldNotes.Magnetometer_hy.acqchan = self.survey_config_dict['hy']

        ch_count = 2
        #--> hz
        try:
            int(self.survey_config_dict['hz'])
            self.mt_obj.FieldNotes.Magnetometer_hz.id = 3
            self.mt_obj.FieldNotes.Magnetometer_hz.chtype = 'hz'
            self.mt_obj.FieldNotes.Magnetometer_hz.x = 0
            self.mt_obj.FieldNotes.Magnetometer_hz.y = 0
            self.mt_obj.FieldNotes.Magnetometer_hz.azm = 90
            self.mt_obj.FieldNotes.Magnetometer_hz.acqchan = self.survey_config_dict['hz']
            ch_count += 1
        except ValueError:
            pass

        #--> ex
        self.mt_obj.FieldNotes.Electrode_ex.id = ch_count+1
        self.mt_obj.FieldNotes.Electrode_ex.chtype = 'ex'
        self.mt_obj.FieldNotes.Electrode_ex.x = 0
        self.mt_obj.FieldNotes.Electrode_ex.y = 0
        self.mt_obj.FieldNotes.Electrode_ex.x2 = float(self.survey_config_dict['e_xaxis_length'])
        self.mt_obj.FieldNotes.Electrode_ex.y2 = 0
        self.mt_obj.FieldNotes.Electrode_ex.azm = float(self.survey_config_dict['e_xaxis_azimuth'])
        self.mt_obj.FieldNotes.Electrode_ex.acqchan = ch_count+1

        #--> ex
        self.mt_obj.FieldNotes.Electrode_ey.id = ch_count+2
        self.mt_obj.FieldNotes.Electrode_ey.chtype = 'ey'
        self.mt_obj.FieldNotes.Electrode_ey.x = 0
        self.mt_obj.FieldNotes.Electrode_ey.y = 0
        self.mt_obj.FieldNotes.Electrode_ey.x2 = 0
        self.mt_obj.FieldNotes.Electrode_ey.y2 = float(self.survey_config_dict['e_yaxis_length'])
        self.mt_obj.FieldNotes.Electrode_ey.azm = float(self.survey_config_dict['e_yaxis_azimuth'])
        self.mt_obj.FieldNotes.Electrode_ey.acqchan = ch_count+2

#        #--> rhx
#        ch_count += 2
#        try:
#            self.mt_obj.Define_measurement.meas_rhx = mtedi.HMeasurement()
#            self.mt_obj.Define_measurement.meas_rhx.id = ch_count+1
#            self.mt_obj.Define_measurement.meas_rhx.chtype = 'rhx'
#            self.mt_obj.Define_measurement.meas_rhx.x = 0
#            self.mt_obj.Define_measurement.meas_rhx.y = 0
#            self.mt_obj.Define_measurement.meas_rhx.azm = 0
#            self.mt_obj.Define_measurement.meas_rhx.acqchan = self.survey_config_dict['rr_hx']
#            ch_count += 1
#        except KeyError:
#            pass
#
#        #--> rhy
#        try:
#            self.mt_obj.Define_measurement.meas_rhy = mtedi.HMeasurement()
#            self.mt_obj.Define_measurement.meas_rhy.id = ch_count+1
#            self.mt_obj.Define_measurement.meas_rhy.chtype = 'rhy'
#            self.mt_obj.Define_measurement.meas_rhy.x = 0
#            self.mt_obj.Define_measurement.meas_rhy.y = 0
#            self.mt_obj.Define_measurement.meas_rhy.azm = 90
#            self.mt_obj.Define_measurement.meas_rhy.acqchan = self.survey_config_dict['rr_hy']
#            ch_count += 1
#        except KeyError:
#            pass
#
#        self.mt_obj.Define_measurement.maxchan = ch_count

    def write_edi_file(self, station=None, birrp_dir=None,
                       survey_config_fn=None, birrp_config_fn=None,
                       copy_path=None):

        """
        Read in BIRRP out puts, in this case the .j file and convert that into
        an .edi file using the survey_config_fn parameters.

        Arguments
        -------------

            **station** : string
                          name of station

            **birrp_dir** : string
                            full path to output directory for BIRRP

            **survey_config_fn** : string
                                  full path to survey configuration file

            **birrp_config_fn** : string
                                  full path to birrp configuration file
                                  *default* is none and is looked for in the
                                  birrp_dir

            **copy_path** : string
                            full path to directory to copy the edi file to

        Outputs
        -------------

            **edi_fn** : string
                         full path to edi file

        .. note::

        The survey_config_fn is a file that has the structure:
            [station]
                b_instrument_amplification = 1
                b_instrument_type = coil
                b_logger_gain = 1
                b_logger_type = zen
                b_xaxis_azimuth = 0
                b_yaxis_azimuth = 90
                box = 26
                date = 2015/06/09
                e_instrument_amplification = 1
                e_instrument_type = Ag-Agcl electrodes
                e_logger_gain = 1
                e_logger_type = zen
                e_xaxis_azimuth = 0
                e_xaxis_length = 100
                e_yaxis_azimuth = 90
                e_yaxis_length = 100
                elevation = 2113.2
                hx = 2274
                hy = 2284
                hz = 2254
                lat = 37.7074236995
                location = Earth
                lon = -118.999542099
                network = USGS
                notes = Generic config file
                rr_box = 25
                rr_date = 2015/06/09
                rr_hx = 2334
                rr_hy = 2324
                rr_lat = 37.6909139779
                rr_lon = -119.028707542
                rr_station = 302
                sampling_interval = all
                save_path = \home\mtdata\survey_01\mt_01
                station = 300
                station_type = mt

        This file can be written using mtpy.utils.configfile::

            >>> import mtpy.utils.configfile as mtcfg
            >>> station_dict = {}
            >>> station_dict['lat'] = 21.346
            >>> station_dict['lon'] = 122.45654
            >>> station_dict['elev'] = 123.43
            >>> cfg_fn = r"\home\mtdata\survey_01"
            >>> mtcfg.write_dict_to_configfile({station: station_dict}, cfg_fn)


        """
        # check to make sure all the files exist
        if station is not None:
            self.station = station

        if self.station is None:
            raise mtex.MTpyError_inputarguments('Need to input the station name')

        # birrp directory
        if birrp_dir is not None:
            self.birrp_dir = birrp_dir

        if not os.path.isdir(self.birrp_dir):
            raise mtex.MTpyError_inputarguments('Could not find {0}, check path'.format(birrp_dir))

        # survey configuratrion
        if survey_config_fn is not None:
            self.survey_config_fn = survey_config_fn
        self.read_survey_config_fn()

        # birrp configuration file
        if birrp_config_fn is not None:
            self.birrp_config_fn = birrp_config_fn
        self.read_birrp_config_fn()

        # get .j file first\
        self.get_j_file()

        # read in .j file
        self.mt_obj = mt.MT()
        self.mt_obj.read_mt_file(self.j_fn)

        # get birrp parameters from .j file if birrp dictionary is None
        if self.birrp_dict is None:
            self.birrp_dict = self.j_obj.header_dict
            for b_key in list(self.birrp_dict.keys()):
                if 'filnam' in b_key:
                    self.birrp_dict.pop(b_key)


        # fill in different blocks of the edi file
        self._fill_site()
        self._fill_info()
        self._fill_field_notes()

        # write edi file
        edi_fn = mtfh.make_unique_filename(os.path.join(self.birrp_dir,
                                                        '{0}.edi'.format(self.station)))

        edi_fn = self.mt_obj._write_edi_file(new_edi_fn=edi_fn)

        return edi_fn
