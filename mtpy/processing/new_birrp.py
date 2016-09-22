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
import time


import mtpy.core.z as mtz
import mtpy.utils.configfile as mtcfg
import mtpy.utils.filehandling as mtfh
import mtpy.utils.exceptions as mtex

#==============================================================================
class BIRRP_Parameters(object):
    """
    class to hold and produce the appropriate parameters given the input 
    parameters.
    """
    
    def __init__(self, ilev=0, **kwargs):
        # set the input level
        self.ilev = ilev

        # get default parameters        
        self._get_parameters()
        
        # set any attributes that are given
        for key in kwargs.keys():
            setattr(self, key, kwargs[key])
            
        self._validate_parameters()
        
    def _get_parameters(self):
        """
        get appropriate parameters
        """
        
        if self.ilev not in [0, 1]:
            raise ValueError('Input level, ilev, must be 0 for basic or 1 d'+ 
                             'for advanced, not {0}'.format(self.ilev))
        
        if self.ilev == 0:
            self.nout = 3
            self.ninp = 2
            self.tbw = 2.0
            self.nfft = 2**18
            self.nsctmax = 14
            self.uin = 0
            self.ainuin = .9999
            self.c2threshe = 0
            self.nz = 0
            self.c2threshe1 = 0
            self.ofil = 'mt'
            self.nlev = 0
            self.nar = 5
            self.imode = 0
            self.jmode = 0
            self.nfil = 0
            self.thetae = [0, 90, 0]
            self.thetab = [0, 90, 0]
            self.thetaf = [0, 90, 0]
            
        elif self.ilev == 1:
            self.nout = 3
            self.ninp = 2
            self.nref = 2
            self.nrr = 1
            self.tbw = 2
            self.nfft = 2**18
            self.nsctinc = 2
            self.nsctmax = np.floor(np.log2(self.nfft))-4
            self.nf1 = self.tbw+2
            self.nfinc = self.tbw
            self.nfsect = 2
            self.uin = 0
            self.ainlin = 0.0001
            self.ainuin = 0.9999
            if self.nrr == 1:
                self.c2threshb = 0.0
                self.c2threshe = 0.0
                if self.c2threshe == 0 and self.c2threshb == 0:
                    self.nz = 0
                else:
                    self.nz = 0
                    self.perlo = 1000
                    self.perhi = .0001
            elif self.nrr == 0:
                self.c2threshb = 0
                self.c2threshe = 0
            self.nprej = 0
            self.prej = None
            self.c2threshe1 = 0
            self.ofil = 'mt'
            self.nlev = 0
            self.nar = 5
            self.imode = 0
            self.jmode = 0
            self.nfil = 0
            self.thetae = [0, 90, 0]
            self.thetab = [0, 90, 0]
            self.thetaf = [0, 90, 0]
            
    def _validate_parameters(self):
        """
        check to make sure the parameters are legit.
        """
        
        # be sure the 
        if self.ninp not in [1, 2, 3]:
            print 'Number of inputs {0} not allowed.'.format(self.ninp)
            self.ninp = 2
            print '  --> setting ninp to {0}'.format(self.ninp)            
        
        if self.nout not in [2, 3]:
            print 'Number of outputs {0} not allowed.'.format(self.nout)
            self.nout = 2
            print '  --> setting nout to {0}'.format(self.nout)
        
        if self.tbw < 0 or self.tbw > 4:
            print 'Total bandwidth of slepian window {0} not allowed.'.format(self.tbw)
            self.tbw = 2
            print '  --> setting tbw to {0}'.format(self.tbw)
            
        if np.remainder(np.log2(self.nfft), 2) != 0:
            print 'Window length nfft should be a power of 2 not {0}'.format(self.nfft)
            self.nfft = 2**np.floor(np.log2(self.nfft))
            print '  -- > setting nfft to {0}, (2**{1:.0f})'.format(self.nfft,
                                                              np.log2(self.nfft))
                                                              
        if np.log2(self.nfft)-self.nsctmax < 4:
            print 'Maximum number of windows {0} is too high'.format(self.nsctmax)
            self.nsctmax = np.log2(self.nfft)-4
            print '  --> setting nsctmax to {0}'.format(self.nsctmax)
            
        if self.uin != 0:
            print 'You\'re playing with fire if uin is not 0.'
            self.uin = 0
            print '  --> setting uin to 0, if you don\'t want that change it back'
        
        if self.imode not in [0, 1, 2, 3]:
            raise BIRRP_Parameter_Error('Invalid number for time series mode,'
                                       'imode, {0}, should be 0, 1, 2, or 3'.format(self.imode))
        
        if self.jmode not in [0, 1]:
            raise BIRRP_Parameter_Error('Invalid number for time mode,'
                                       'imode, {0}, should be 0, or 1'.format(self.imode))
        if self.ilev == 1:
            if self.nsctinc != 2:
                print '!WARNING! Check decimation increment nsctinc, should be 2 not {0}'.format(self.nsctinc)
                
            if self.nfsect != 2:
                print 'Will get an error from BIRRP if nfsect is not 2.'
                print 'number of frequencies per section is {0}'.format(self.nfsect)
                self.nfsect = 2
                print '  --> setting nfsect to 2'
                
            if self.nf1 != self.tbw+2:
                print '!WARNING! First frequency should be around tbw+2.'
                print 'nf1 currently set to {0}'.format(self.nf1)
                
            if self.nfinc != self.tbw:
                print '!WARNING! sequence of frequencies per window should be around tbw.'
                print 'nfinc currently set to {0}'.format(self.nfinc) 
                
            if self.nprej != 0:
                if self.prej is None or type(self.prej) is not list:
                    raise BIRRP_Parameter_Error('Need to input a prejudice list if nprej != 0'+
                                     '\nInput as a list of frequencies' )
                                     
            if self.nrr not in [0, 1]:
                print('!WARNING! Value for picking remote reference or '+ 
                     'two stage processing, nrr, '+ 
                     'should be 0 or 1 not {0}'.format(self.nrr))
                self.nrr = 0
                print '  --> setting nrr to {0}'.format(self.nrr)
                
                
            if self.c2threshe != 0 or self.c2threshb != 0:
                if not self.perhi:
                    raise BIRRP_Parameter_Error('Need to input a high period (s) threshold as perhi')
        
                if not self.perlo:
                    raise BIRRP_Parameter_Error('Need to input a low period (s) threshold as perlo')
         
        if len(self.thetae) != 3:
            print 'Electric rotation angles not input properly {0}'.format(self.thetae)
            print 'input as north, east, orthogonal rotation'            
            self.thetae = [0, 90, 0]
            print '  --> setting thetae to {0}'.format(self.thetae)
        
        if len(self.thetab) != 3:
            print 'Magnetic rotation angles not input properly {0}'.format(self.thetab)
            print 'input as north, east, orthogonal rotation'            
            self.thetab = [0, 90, 0]
            print '  --> setting thetab to {0}'.format(self.thetab)
            
                
        if len(self.thetaf) != 3:
            print 'Fiedl rotation angles not input properly {0}'.format(self.thetaf)
            print 'input as north, east, orthogonal rotation'            
            self.thetaf = [0, 90, 0]
            print '  --> setting thetaf to {0}'.format(self.thetaf)


    def read_config_file(self, birrp_config_fn):
        """
        read in a configuration file and fill in the appropriate parameters
        """
        
        birrp_dict = mtcfg.read_configfile(birrp_config_fn)
        
        for birrp_key in birrp_dict.keys():
            try:
                b_value = float(birrp_dict[birrp_key])
                if np.remainder(b_value, 1) == 0:
                    b_value = int(b_value)
            except ValueError:
                b_value = birrp_dict[birrp_key]
                
            setattr(self, birrp_key, b_value)
            
    def write_config_file(self, birrp_dict, save_fn):
        """
        write a config file for birrp parameters
        """
        
        cfg_fn = mtfh.make_unique_filename('{0}_birrp_params.cfg'.format(save_fn))
                                                                 
        mtcfg.write_dict_to_configfile(birrp_dict, cfg_fn)
        print 'Wrote BIRRP config file for edi file to {0}'.format(cfg_fn)
#==============================================================================
# Error classes            
#==============================================================================
class BIRRP_Parameter_Error(Exception):
    pass

#==============================================================================
# write script file
#==============================================================================
def write_script_file(processing_dict, save_path=None):
    """
    writescript_fn(processingdict will write a script file for BIRRP using 
    info in processingdict which is a dictionary with keys:
    
    ================== ========================================================
    parameter          description
    ================== ======================================================== 
    station            station name
    fn_list             list of file names to be processed, this must be in 
                       the correct order [EX, EY, HZ, HX, HY] and if multiple
                       sections are to be processed at the same time then 
                       must be input as a nested loop 
                       [[EX1, EY1, HZ1, HX1, HY1], 
                       [EX2, EY2, HZ2, HX2, HY2], ...]
    rrfn_list           list of remote reference file names, similar to the 
                       fn_list [[HX1, HY1], [HX2, HY2], ...]
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
    hx_cal             full path to calibration file for hx 
    hy_cal             full path to calibration file for hy 
    hz_cal             full path to calibration file for hz 
    rrhx_cal           full path to calibration file for remote reference hx 
    rrhy_cal           full path to calibration file for remote reference hy 
                       Note that BIRRP assumes the calibrations are the same
                       as in the first segment, so the remote reference must
                       be the same for all segments, or at least the same
                       instruments have to be used so the calibration is 
                       the same.
    ================== ========================================================              
    
    .. seealso:: BIRRP Manual and publications by Chave and Thomson
                for more details on the parameters found at:
                
                http://www.whoi.edu/science/AOPE/people/achave/Site/Next1.html
            
    Arguments
    -----------
        **processing_dict** : dictionary with keys as above
        
        **save_path** : string (full path to directory to save script file)
                        if none saves as:
                            os.path.join(os.path.dirname(fn_list[0]),'BF')
        
    Outputs
    ----------
        **script_file** : full path to script file to guide birrp
        
        **birrp_dict** : dictionary of birrp parameters input into script file
        
    
    """
     
    #===================================================================
    # Write a script file for BIRRP, Chave et al. [2004]            
    #===================================================================
    pdict = dict(processing_dict)
    
    try:
        fn_array = np.array(pdict['fn_list'])
    except KeyError:
        raise KeyError('fn_list --> Need to input a list of files to process')
    
    try:
        nds, ndf = fn_array.shape
    except ValueError:
        ndf = fn_array.shape[0]
        nds = 0
    if save_path is None:
        if nds == 0:
            bfpath = os.path.join(os.path.dirname(pdict['fn_list'][0]),
                                 'BF')
        else:
            bfpath = os.path.join(os.path.dirname(pdict['fn_list'][0][0]),
                                  'BF')
    else:
        bfpath = save_path
    
    
    if nds == 0:
        npcs = 1
    elif nds == 1:
        nds = 0
        npcs = 1
        pdict['fn_list'] = pdict['fn_list'][0]
        try:
            pdict['rrfn_list'] = pdict['rrfn_list'][0]  
        except KeyError:
            pass

    else:
        npcs = int(nds)
        
    #make a directory to put BIRRP Files (BF)
    if not os.path.exists(bfpath):
        os.mkdir(bfpath)
        print 'Made directory: ', bfpath

    #output file stem, full path
    ofil = os.path.join(bfpath,pdict['station']) 
    
    #mode to process: default is basic
    ilev = int(pdict.pop('ilev', 0))
    
    #number of output channels
    try:
        nout = int(pdict['nout'])
    except KeyError:
        if nds != 0:
            if ndf == 5:
                nout = 3
            else:
                nout = 2
        elif ndf == 5:
            nout = 3
        else:
            nout = 2
    
    #number of input channels default is 2
    ninp = int(pdict.pop('ninp', 2))
        
    #time bandwidth window size            
    tbw = int(pdict.pop('tbw', 2))
        
    #------------Options for Advanced mode-------------------------
    if ilev == 1:
        #Advanced: number of remote reference channels
        nref = int(pdict.pop('nref', 2))
            
        #Advanced: remote reference type processing
        nrr = int(pdict.pop('nrr', 1))          
            
        #Advanced: magnetic coherence threshold
        c2threshb = float(pdict.pop('c2threshb', 0))
            
        #Advanced: window increment divisor
        nsctinc = int(pdict.pop('nsctinc', 2))
            
        #Advanced: first frequency to extract
        nf1 = int(pdict.pop('nf1', tbw+2))
        
        #Advanced: frequency increment
        nfinc = int(pdict.pop('nfinc', tbw))
            
        #number of frequencies to extract
        nfsect = int(pdict.pop('nfsec', 2))
            
        #number AR filter is divided by
        mfft = int(pdict.pop('mfft', 2))
        
        #Advanced: lower bound of leverage point rejection
        ainlin = float(pdict.pop('ainlin', .0001))
        
        #Advanced: coherence threshold low period
        perlo = int(pdict.pop('perlo', 1000))
            
        #Advanced: coherenct threshold high period
        perhi = float(pdict.pop('perhi', .0001))
            
        #Advanced:  number of frequencies to reject
        nprej = int(pdict.pop('nprej', 0))
        
        #Advanced
        try:
            prej = pdict['prej'].split(',')
            if type(prej) is list:
                prej = [float(ff) for ff in prej]
            if nprej != len(prej):
                nprej = len(prej)
        except KeyError:
            prej = []
            
    #---------------------Options for Basic Mode--------------------------
    
    #time series sampling rate
    deltat = int(pdict.pop('deltat', -100))
    
    #max length of fft window 
    nfft = int(pdict.pop('nfft', 2**16))

    #maximum number of sections
    nsctmax = int(pdict.pop('nsctmax', 12))
    
    #quantile factor
    uin = int(pdict.pop('uin', 0))            
    
    #upper bound of leverage point rejection
    ainuin = float(pdict.pop('ainuin', .9999))
 
    #electric channel coherence threshold
    c2threshe = int(pdict.pop('c2threshe', 0))

    #Bz coherency threshold mode
    try:
        nz = int(pdict['nz'])
    except KeyError:
        if nout == 3:
            nz = 0
        else:
            nz = None
    
    #Bz coherence threshold
    try:
        c2threshe1 = float(pdict['c2threshe1'])
    except KeyError:
        if nout == 3:
            c2threshe1 = 0
        else:
            c2threshe1 = None
    
    #output level
    nlev = int(pdict.pop('nlev', 0))
    
    #order of prewhitening auto regressive filter
    nar = int(pdict.pop('nar', 5))
        
    #input mode
    imode = int(pdict.pop('imode', 0))
    
    #output mode
    jmode = int(pdict.pop('jmode', 0))       
    
    #name of filter file
    nfil = int(pdict.pop('nfil', 0))
    
    #calibration for coils
    hx_cal = pdict.pop('hx_cal', None) 
    hy_cal = pdict.pop('hy_cal', None) 
    hz_cal = pdict.pop('hz_cal', None)
    
    rrhx_cal = pdict.pop('rrhx_cal', None) 
    rrhy_cal = pdict.pop('rrhy_cal', None) 
    
    if jmode == 0:
        #number of points to read
        nread = pdict.pop('nread', 1440000)
        if type(nread) is not list and type(nread) is not np.ndarray:
            nread = int(nread)
        #number of points to skip in time series
        try:
            nskip = pdict['nskip']
            if nds != 0 and type(nskip) is not list and \
               type(nskip) is not np.ndarray:
                nskip = [nskip for ii in range(nds)]
        except KeyError:
            if nds != 0:
                nskip = [0 for ii in range(nds)]
            else:
                nskip = 0
        
        #number of point to skip from remote reference time series    
        try:
            nskipr = pdict['nskipr']
            if nds != 0 and type(nskipr) is not list and \
               type(nskipr) is not np.ndarray:
                nskipr = [nskipr for ii in range(nds)]
        except KeyError:
            if nds==0:
                nskipr = 0
            else:
                nskipr = [0 for ii in range(nds)]
    
    if jmode == 1:
        #start time of data
        dstim = pdict.pop('dstim', '1970-01-01 00:00:00')
        
        #window start time
        wstim = pdict.pop('wstim', '1970-01-01 00:00:00')
        
        #window end time
        wetim = pdict.pop('wetim', '1970-01-02 00:00:00')
        
    
    #rotation angle of electric channels
    thetae = pdict.pop('thetae', '0,90,0')
    
    #rotation angle of magnetic channels
    thetab = pdict.pop('thetab', '0,90,0')
    
    #rotation angle of final impedance tensor
    thetaf = pdict.pop('thetaf', '0,90,0')

    #===================================================================
    # Write values to a .script file
    #===================================================================
    #print '+++ ',nskipr
    #print ndf,nds
    #write to a file
    script_fn = ofil+'.script'
    s_lines = []
    fid = file(script_fn,'w')
    if ilev==0: 
        s_lines.append('{0:d}'.format(ilev))
        s_lines.append('{0:d}'.format(nout))
        s_lines.append('{0:d}'.format(ninp))
        s_lines.append('{0:.3f}'.format(tbw))
        s_lines.append('{0:.3f}'.format(deltat))
        s_lines.append('{0:d},{1:d}'.format(nfft,nsctmax))
        s_lines.append('y')
        s_lines.append('{0:.5f},{1:.5f}'.format(uin,ainuin))
        s_lines.append('{0:.3f}'.format(c2threshe))
        #parameters for bz component if ninp=3
        if nout==3:
            if c2threshe==0:
                s_lines.append('{0:d}'.format(0))
                s_lines.append('{0:.3f}'.format(c2threshe1))
            else:
                s_lines.append('{0:d}'.format(nz))
                s_lines.append('{0:.3f}'.format(c2threshe1))
        else:
            pass
        s_lines.append(ofil)
        s_lines.append('{0:d}'.format(nlev))
        
    elif ilev == 1:
        print 'Writing Advanced mode'
        s_lines.append('{0:d}'.format(ilev))
        s_lines.append('{0:d}'.format(nout))
        s_lines.append('{0:d}'.format(ninp))
        s_lines.append('{0:d}'.format(nref))
        if nref>3:
            nrrlist=np.array([len(rrlist) 
                            for rrlist in pdict['rrfn_list']])
            nr3=len(np.where(nrrlist==3)[0])
            nr2=len(np.where(nrrlist==2)[0])
            s_lines.append('{0:d},{1:d}'.format(nr3,nr2))
        s_lines.append('{0:d}'.format(nrr))
        #if remote referencing
        if int(nrr) == 0:
            s_lines.append('{0:.3f}'.format(tbw))
            s_lines.append('{0:.3f}'.format(deltat))
            s_lines.append('{0:d},{1:.2g},{2:d}'.format(nfft,nsctinc,nsctmax))
            s_lines.append('{0:d},{1:.2g},{2:d}'.format(nf1,nfinc,nfsect))
            s_lines.append('y')
            s_lines.append('{0:.2g}'.format(mfft))        
            s_lines.append('{0:.5g},{1:.5g},{2:.5g}'.format(uin,ainlin,ainuin))
            s_lines.append('{0:.3f}'.format(c2threshe))
            #parameters for bz component if ninp=3
            if nout == 3:
                if c2threshe != 0:
                    s_lines.append('{0:d}'.format(nz))
                    s_lines.append('{0:.3f}'.format(c2threshe1))
                else:
                    s_lines.append('{0:d}'.format(0))
                    s_lines.append('{0:.3f}'.format(c2threshe1))
                if c2threshe1 != 0.0 or c2threshe != 0.0:
                    s_lines.append('{0:.6g},{1:.6g}'.format(perlo,perhi))
            else:
                if c2threshe != 0.0:
                    s_lines.append('{0:.6g},{1:.6g}'.format(perlo,perhi))
            s_lines.append(ofil)
            s_lines.append('{0:d}'.format(nlev))
            s_lines.append('{0:d}'.format(nprej))
            if nprej!=0:
                if type(prej) is not list:
                    prej = [prej]
                s_lines.appendlines(['{0:.5g}'.format(nn) for nn in prej])
        #if 2 stage processing
        elif int(nrr) == 1:
            s_lines.append('{0:.5g}'.format(tbw))
            s_lines.append('{0:.5g}'.format(deltat))
            s_lines.append('{0:d},{1:.2g},{2:d}'.format(nfft,nsctinc,nsctmax))        
            s_lines.append('{0:d},{1:.2g},{2:d}'.format(nf1,nfinc,nfsect))
            s_lines.append('y')
            s_lines.append('{0:.2g}'.format(mfft))        
            s_lines.append('{0:.5g},{1:.5g},{2:.5g}'.format(uin,ainlin,ainuin))
            s_lines.append('{0:.3f}'.format(c2threshb))        
            s_lines.append('{0:.3f}'.format(c2threshe))
            if nout == 3:
                if c2threshb != 0 or c2threshe != 0:
                    s_lines.append('{0:d}'.format(nz))
                    s_lines.append('{0:.3f}'.format(c2threshe1))
                elif c2threshb == 0 and c2threshe == 0:
                    s_lines.append('{0:d}'.format(0))
                    s_lines.append('{0:.3f}'.format(0))
            if c2threshb != 0.0 or c2threshe != 0.0:
                s_lines.append('{0:.6g},{1:.6g}'.format(perlo,perhi))
            s_lines.append(ofil)
            s_lines.append('{0:d}'.format(nlev))
            s_lines.append('{0:d}'.format(nprej))
            if nprej!=0:
                if type(prej) is not list:
                    prej = [prej]
                s_lines.appendlines(['{0:.5g}'.format(nn) for nn in prej])
        
    s_lines.append('{0:d}'.format(npcs))    
    s_lines.append('{0:d}'.format(nar))    
    s_lines.append('{0:d}'.format(imode))    
    s_lines.append('{0:d}'.format(jmode)) 
    
    #!!!NEED TO SORT FILE NAMES SUCH THAT EX, EY, HZ, HX, HY or EX, EY, HX, HY
    
    #write in filenames 
    if npcs != 1:
        if jmode == 0:
            s_lines.append(str(nread[0]))
            #--> write filenames to process with other information for first
            #    time section
            for tt, tfile in enumerate(pdict['fn_list'][0]):
                #write in calibration files if given
                if tt == 2:
                    if hx_cal is not None:
                        s_lines.append('-2')
                        s_lines.append(hx_cal)
                    else:
                        s_lines.append(str(nfil))
                elif tt == 3:
                    if hy_cal is not None:
                        s_lines.append('-2')
                        s_lines.append(hy_cal)
                    else:
                        s_lines.append(str(nfil))
                elif tt == 4:
                    if hz_cal is not None:
                        s_lines.append('-2')
                        s_lines.append(hz_cal)
                    else:
                        s_lines.append(str(nfil))
                else:
                    s_lines.append(str(nfil))
                s_lines.append(tfile)
                s_lines.append(str(nskip[0]))
                
            #--> write remote reference time series
            for rr, rfile in enumerate(pdict['rrfn_list'][0]):
                if rr == 0:
                    if rrhx_cal is not None:
                        s_lines.append('-2')
                        s_lines.append(rrhx_cal)
                    else:
                        s_lines.append(str(nfil))
                elif rr == 1:
                    if rrhy_cal is not None:
                        s_lines.append('-2')
                        s_lines.append(rrhy_cal)
                    else:
                        s_lines.append(str(nfil))
                s_lines.append(rfile)
                s_lines.append(str(nskipr[0]))
            
            #--> write in other pieces if there are more, note calibrations 
            #    are only given for the first time block so it is assumed
            #    that the same remote referenc is used for all time blocks
            for nn in range(1,npcs):
                s_lines.append(str(nread[nn]))            
                #write filenames
                for tfile in pdict['fn_list'][nn]:
                    s_lines.append(tfile)
                    s_lines.append(str(nskip[nn]))
                for rfile in pdict['rrfn_list'][nn]:
                    s_lines.append(rfile)
                    s_lines.append(str(nskipr[nn]))
                    
        #--> if start and end time are give write in those
        elif jmode == 1:
            #write filenames
            for tt, tfile in enumerate(pdict['fn_list'][0]):
                #write in calibration files if given
                if tt == 2:
                    if hx_cal is not None:
                        s_lines.append('-2')
                        s_lines.append(hx_cal)
                    else:
                        s_lines.append(str(nfil))
                elif tt == 3:
                    if hy_cal is not None:
                        s_lines.append('-2')
                        s_lines.append(hy_cal)
                    else:
                        s_lines.append(str(nfil))
                elif tt == 4:
                    if hz_cal is not None:
                        s_lines.append('-2')
                        s_lines.append(hz_cal)
                    else:
                        s_lines.append(str(nfil))
                else:
                    s_lines.append(str(nfil))
                s_lines.append(tfile)
                s_lines.append(dstim)
                s_lines.append(wstim)
                s_lines.append(wetim)
                
            #--> write remote referenc information
            for rr, rfile in enumerate(pdict['rrfn_list'][0]):
                if rr == 0:
                    if rrhx_cal is not None:
                        s_lines.append('-2')
                        s_lines.append(rrhx_cal)
                    else:
                        s_lines.append(str(nfil))
                if rr == 1:
                    if rrhy_cal is not None:
                        s_lines.append('-2')
                        s_lines.append(rrhy_cal)
                    else:
                        s_lines.append(str(nfil))
                s_lines.append(rfile)
                s_lines.append(dstim)
                s_lines.append(wstim)
                s_lines.append(wetim)
                
            #--> write other time blocks
            for nn in range(1,npcs):
                s_lines.append(str(nread[nn]))            
                #write filenames
                for tfile in pdict['fn_list'][nn]:
                    s_lines.append(tfile)
                    s_lines.append(dstim)
                    s_lines.append(wstim)
                    s_lines.append(wetim)
                for rfile in pdict['rrfn_list'][nn]:
                    s_lines.append(rfile)
                    s_lines.append(dstim)
                    s_lines.append(wstim)
                    s_lines.append(wetim)
    else:
        if jmode == 0:
            if type(nread) is list:
                s_lines.append(str(nread[0]))
            else:
                s_lines.append(str(nread))
            #--> write filenames for first block
            if nds==0:
                for tt, tfile in enumerate(pdict['fn_list']):
                    if tt == 2:
                        if hx_cal is not None:
                            s_lines.append('-2')
                            s_lines.append(hx_cal)
                        else:
                            s_lines.append(str(nfil))
                    elif tt == 3:
                        if hy_cal is not None:
                            s_lines.append('-2')
                            s_lines.append(hy_cal)
                        else:
                            s_lines.append(str(nfil))
                    elif tt == 4:
                        if hz_cal is not None:
                            s_lines.append('-2')
                            s_lines.append(hz_cal)
                        else:
                            s_lines.append(str(nfil))
                    else:
                        s_lines.append(str(nfil))
                    s_lines.append(tfile)
                    if type(nskip) is list:
                        s_lines.append(str(nskip[0]))
                    else:
                        s_lines.append(str(nskip))
                for rr, rfile in enumerate(pdict['rrfn_list']):
                    if rr == 0:
                        if rrhx_cal is not None:
                            s_lines.append('-2')
                            s_lines.append(rrhx_cal)
                        else:
                            s_lines.append(str(nfil))
                    if rr == 1:
                        if rrhy_cal is not None:
                            s_lines.append('-2')
                            s_lines.append(rrhy_cal)
                        else:
                            s_lines.append(str(nfil))
                    s_lines.append(rfile)
                    if type(nskipr) is list:
                        s_lines.append(str(nskipr[0]))
                    else:
                        s_lines.append(str(nskipr))
            else:
                for tfile in pdict['fn_list'][0]:
                    s_lines.append(str(nfil))
                    s_lines.append(tfile)
                    if type(nskip) is list:
                        s_lines.append(str(nskip[0]))
                    else:
                        s_lines.append(str(nskip))
                for rfile in pdict['rrfn_list'][0]:
                    s_lines.append(str(nfil))
                    s_lines.append(rfile)
                    if type(nskipr) is list:
                        s_lines.append(str(nskipr[0]))
                    else:
                        s_lines.append(str(nskipr))
                        
        elif jmode == 1:
            #write filenames
            if nds==0:
                for tt, tfile in enumerate(pdict['fn_list']):
                    if tt == 2:
                        if hx_cal is not None:
                            s_lines.append('-2')
                            s_lines.append(hx_cal)
                        else:
                            s_lines.append(str(nfil))
                    elif tt == 3:
                        if hy_cal is not None:
                            s_lines.append('-2')
                            s_lines.append(hy_cal)
                        else:
                            s_lines.append(str(nfil))
                    elif tt == 4:
                        if hz_cal is not None:
                            s_lines.append('-2')
                            s_lines.append(hz_cal)
                        else:
                            s_lines.append(str(nfil))
                    else:
                        s_lines.append(str(nfil))
                    s_lines.append(tfile)
                    s_lines.append(dstim)
                    s_lines.append(wstim)
                    s_lines.append(wetim)
                for rr, rfile in enumerate(pdict['rrfn_list']):
                    if rr == 0:
                        if rrhx_cal is not None:
                            s_lines.append('-2')
                            s_lines.append(rrhx_cal)
                        else:
                            s_lines.append(str(nfil))
                    if rr == 1:
                        if rrhy_cal is not None:
                            s_lines.append('-2')
                            s_lines.append(rrhy_cal)
                        else:
                            s_lines.append(str(nfil))
                    s_lines.append(rfile)
                    s_lines.append(dstim)
                    s_lines.append(wstim)
                    s_lines.append(wetim)
            else:
                for tt, tfile in enumerate(pdict['fn_list'][0]):
                    if tt == 2:
                        if hx_cal is not None:
                            s_lines.append('-2')
                            s_lines.append(hx_cal)
                        else:
                            s_lines.append(str(nfil))
                    elif tt == 3:
                        if hy_cal is not None:
                            s_lines.append('-2')
                            s_lines.append(hy_cal)
                        else:
                            s_lines.append(str(nfil))
                    elif tt == 4:
                        if hz_cal is not None:
                            s_lines.append('-2')
                            s_lines.append(hz_cal)
                        else:
                            s_lines.append(str(nfil))
                    else:
                        s_lines.append(str(nfil))
                    s_lines.append(tfile)
                    s_lines.append(dstim)
                    s_lines.append(wstim)
                    s_lines.append(wetim)
                for rr, rfile in enumerate(pdict['rrfn_list'][0]):
                    if rr == 0:
                        if rrhx_cal is not None:
                            s_lines.append('-2')
                            s_lines.append(rrhx_cal)
                        else:
                            s_lines.append(str(nfil))
                    if rr == 1:
                        if rrhy_cal is not None:
                            s_lines.append('-2')
                            s_lines.append(rrhy_cal)
                        else:
                            s_lines.append(str(nfil))
                    s_lines.append(rfile)
                    s_lines.append(dstim)
                    s_lines.append(wstim)
                    s_lines.append(wetim)
                    
    #write rotation angles
    s_lines.append(thetae.replace(',',' '))
    s_lines.append(thetab.replace(',',' '))
    s_lines.append(thetaf.replace(',',' '))    
    
    with open(script_fn, 'w') as fid:
        fid.write('\n'.join(s_lines))
    
    birrp_dict = {}
    
    if ilev == 0:
        birrp_dict['ilev'] = ilev
        birrp_dict['nout'] = nout
        birrp_dict['ninp'] = ninp
        birrp_dict['tbw'] = tbw
        birrp_dict['sampling_rate'] = deltat
        birrp_dict['nfft'] = nfft
        birrp_dict['nsctmax'] = nsctmax
        birrp_dict['uin'] = uin
        birrp_dict['ainuin'] = ainuin
        birrp_dict['c2threshe'] = c2threshe
        birrp_dict['nz'] = nz
        birrp_dict['c2threshe1'] = c2threshe1
        birrp_dict['ofil'] = ofil
        birrp_dict['nlev'] = nlev
        birrp_dict['npcs'] = npcs
        birrp_dict['n_samples'] = nread
        birrp_dict['nar'] = nar
        birrp_dict['imode'] = imode
        birrp_dict['jmode'] = jmode
        birrp_dict['nfil'] = nfil
        birrp_dict['nskip'] = nskip
        birrp_dict['nskipr'] = nskipr
        birrp_dict['thetae'] = thetae
        birrp_dict['thetab'] = thetab
        birrp_dict['thetaf'] = thetaf
    elif ilev == 1:
        birrp_dict['ilev'] = ilev
        birrp_dict['nout'] = nout
        birrp_dict['ninp'] = ninp
        birrp_dict['nref'] = nref
        birrp_dict['nrr'] = nrr
        birrp_dict['tbw'] = tbw
        birrp_dict['sampling_rate'] = deltat
        birrp_dict['nfft'] = nfft
        birrp_dict['nsctinc'] = nsctinc
        birrp_dict['nsctmax'] = nsctmax
        birrp_dict['nf1'] = nf1
        birrp_dict['nfinc'] = nfinc
        birrp_dict['nfsect'] = nfsect
        birrp_dict['uin'] = uin
        birrp_dict['ainlin'] = ainlin
        birrp_dict['ainuin'] = ainuin
        if nrr == 1:
            birrp_dict['c2threshb'] = c2threshb
            birrp_dict['c2threshe'] = c2threshe
            if c2threshe == 0 and c2threshb == 0:
                birrp_dict['nz'] = 0
            else:
                birrp_dict['nz'] = 0
                birrp_dict['perlo'] = perlo
                birrp_dict['perhi'] = perhi
        elif nrr == 0:
            birrp_dict['c2threshb'] = 0
            birrp_dict['c2threshe'] = c2threshe
        birrp_dict['nprej'] = nprej
        birrp_dict['prej'] = prej
        birrp_dict['c2threshe1'] = c2threshe1
        birrp_dict['ofil'] = ofil
        birrp_dict['npcs'] = npcs
        birrp_dict['n_samples'] = nread
        birrp_dict['nlev'] = nlev
        birrp_dict['nar'] = nar
        birrp_dict['imode'] = imode
        birrp_dict['jmode'] = jmode
        birrp_dict['nfil'] = nfil
        if jmode == 0:
            birrp_dict['nskip'] = nskip
            birrp_dict['nskipr'] = nskipr
        elif jmode == 1:
            birrp_dict['dstim'] = dstim
            birrp_dict['wstim'] = wstim
            birrp_dict['wetim'] = wetim
        birrp_dict['thetae'] = thetae
        birrp_dict['thetab'] = thetab
        birrp_dict['thetaf'] = thetaf

    print 'Wrote BIRRP script file: {0}'.format(script_fn)
    
    return script_fn, birrp_dict
    
#==============================================================================
# run birrp
#==============================================================================
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

    # get an input string for communicating with the birrp executable
    with open(script_file, 'r') as sfid:
        input_string = ''.join(sfid.readlines())

    #correct inputstring for potential errorneous line endings due to strange
    #operating systems:
    temp_string = input_string.split()
    temp_string = [i.strip() for i in temp_string]
    input_string = '\n'.join(temp_string)
    input_string += '\n'

    #open a log file to catch process and errors of BIRRP executable
    log_file = open('birrp_logfile.log','w')

    print 'Starting Birrp processing at {0}...'.format(time.ctime())

    birrp_process = subprocess.Popen(birrp_exe, 
                                     stdin=subprocess.PIPE, 
                                     stdout=log_file,
                                     stderr=log_file)

    out, err = birrp_process.communicate(input_string)
    
    log_file.close()

    print 'Ending Birrp processing at   {0}...'.format(time.ctime())
    print 'Closed log file: {0}'.format(log_file.name)
 
    print 'Outputs: {0}'.format(out)
    print 'Errors: {0}'.format(err)
    
    #go back to initial directory
    os.chdir(current_dir)

    print '\n{0} DONE !!! {0}\n'.format('='*20)
#==============================================================================
# Class to read j_file
#==============================================================================
class JFile(object):
    """
    be able to read and write a j-file
    """
    
    def __init__(self, j_fn=None):
        self.j_fn = j_fn
        self.header_dict = None
        self.metadata_dict = None
        self._j_lines = None
        self.Z = None
        self.Tipper = None
        
    def _get_j_lines(self, j_fn=None):
        """
        read in the j_file as a list of lines, put the lines in attribute
        _j_lines
        """
        if j_fn is not None:
            self.j_fn = j_fn
            
        if os.path.isfile(os.path.abspath(self.j_fn)) is False:
            raise IOError('Could not find {0}, check path'.format(self.j_fn))
            
        with open(self.j_fn, 'r') as fid:
            self._j_lines = fid.readlines()
        
    def read_header(self, j_lines=None, j_fn=None):
        """
        Parsing the header lines of a j-file to extract processing information.
    
        Input:
        - j-file as list of lines (output of readlines())
    
        Output:
        - Dictionary with all parameters found

        """
        if j_lines is not None:
            self._j_lines = j_lines
            
        if j_fn is not None:
            self.j_fn = j_fn
            
        if self._j_lines is None:
            self._get_j_lines()
            
        header_lines = [j_line for j_line in self._j_lines if '#' in j_line]
        header_dict = {'title':header_lines[0][1:].strip()}
        
        fn_count = 0
        theta_count = 0
        # put the information into a dictionary 
        for h_line in header_lines[1:]:
            # replace '=' with a ' ' to be sure that when split is called there is a
            # split, especially with filenames
            h_list = h_line[1:].strip().replace('=', ' ').split()
            # skip if there is only one element in the list
            if len(h_list) == 1:
                continue
            # get the key and value for each parameter in the given line
            for h_index in range(0, len(h_list), 2):
                h_key = h_list[h_index]
                # if its the file name, make the dictionary value be a list so that 
                # we can append nread and nskip to it, and make the name unique by
                # adding a counter on the end
                if h_key == 'filnam':
                    h_key = '{0}_{1:02}'.format(h_key, fn_count)
                    fn_count += 1
                    h_value = [h_list[h_index+1]]
                    header_dict[h_key] = h_value
                    continue
                elif h_key == 'nskip' or h_key == 'nread':
                    h_key = 'filnam_{0:02}'.format(fn_count-1)
                    h_value = int(h_list[h_index+1])
                    header_dict[h_key].append(h_value)
                    
                # if its the line of angles, put them all in a list with a unique key
                elif h_key == 'theta1':
                    h_key = '{0}_{1:02}'.format(h_key, theta_count)
                    theta_count += 1
                    h_value = float(h_list[h_index+1])
                    header_dict[h_key] = [h_value]
                elif h_key == 'theta2' or h_key == 'phi':
                    h_key = '{0}_{1:02}'.format('theta1', theta_count-1)
                    h_value = float(h_list[h_index+1])
                    header_dict[h_key].append(h_value)
                    
                else:
                    try:
                        h_value = float(h_list[h_index+1])
                    except ValueError:
                        h_value = h_list[h_index+1]
                    
                    header_dict[h_key] = h_value
            
        self.header_dict = header_dict
        
    def read_metadata(self, j_lines=None, j_fn=None):
        """
        read in the metadata of the station, or information of station 
        logistics like: lat, lon, elevation
        
        Not really needed for a birrp output since all values are nan's
        """
        
        if j_lines is not None:
            self._j_lines = j_lines
            
        if j_fn is not None:
            self.j_fn = j_fn
            
        if self._j_lines is None:
            self._get_j_lines()
        
        metadata_lines = [j_line for j_line in self._j_lines if '>' in j_line]
    
        metadata_dict = {}
        for m_line in metadata_lines:
            m_list = m_line.strip().split('=')
            m_key = m_list[0][1:].strip().lower()
            try:
                m_value = float(m_list[0].strip())
            except ValueError:
                m_value = 0.0
                
            metadata_dict[m_key] = m_value
            
        self.metadata_dict = metadata_dict
        
    def read_j_file(self, j_lines=None, j_fn=None):
        """
        read_j_file will read in a *.j file output by BIRRP (better than reading lots of *.<k>r<l>.rf files)
    
        Input:
        j-filename
    
        Output: 4-tuple
        - periods : N-array
        - Z_array : 2-tuple - values and errors
        - tipper_array : 2-tuple - values and errors
        - processing_dict : parsed processing parameters from j-file header
    
        """   
    
        # read data
        z_index_dict = {'zxx':(0, 0),
                        'zxy':(0, 1),
                        'zyx':(1, 0),
                        'zyy':(1, 1)}
        t_index_dict = {'tzx':(0, 0),
                        'tzy':(0, 1)}
                        
        if j_lines is not None:
            self._j_lines = j_lines
            
        if j_fn is not None:
            self.j_fn = j_fn
            
        if self._j_lines is None:
            self._get_j_lines()
            
        self.header_dict = self.read_header()
        self.metadata_dict = self.read_metadata()
        
        data_lines = [j_line for j_line in self._j_lines 
                      if not '>' in j_line and not '#' in j_line][1:]
                
        # sometimes birrp outputs some missing periods, so the best way to deal with 
        # this that I could come up with was to get things into dictionaries with 
        # key words that are the period values, then fill in Z and T from there
        # leaving any missing values as 0
        
        # make empty dictionary that have keys as the component 
        z_dict = dict([(z_key, {}) for z_key in z_index_dict.keys()])
        t_dict = dict([(t_key, {}) for t_key in t_index_dict.keys()])
        for d_line in data_lines:
            # check to see if we are at the beginning of a component block, if so 
            # set the dictionary key to that value
            if 'z' in d_line.lower():
                d_key = d_line.strip().split()[0].lower()
            
            # if we are at the number of periods line, skip it
            elif len(d_line.strip().split()) == 1:
                continue
            # get the numbers into the correct dictionary with a key as period and
            # for now we will leave the numbers as a list, which we will parse later
            else:
                # split the line up into each number
                d_list = d_line.strip().split()
                
                # make a copy of the list to be sure we don't rewrite any values,
                # not sure if this is necessary at the moment
                d_value_list = list(d_list)
                for d_index, d_value in enumerate(d_list):
                    # check to see if the column number can be converted into a float
                    # if it can't, then it will be set to 0, which is assumed to be
                    # a masked number when writing to an .edi file
                    try:
                        d_value_list[d_index] = float(d_value)
                    except ValueError:
                        d_value_list[d_index] = 0.0
                
                # put the numbers in the correct dictionary as:
                # key = period, value = [real, imaginary, error]
                if d_key in z_index_dict.keys():
                    z_dict[d_key][d_value_list[0]] = d_value_list[1:4]
                elif d_key in t_index_dict.keys():
                    t_dict[d_key][d_value_list[0]] = d_value_list[1:4]
                    
        # now we need to get the set of periods for all components
        all_periods = sorted(list(set(np.array([z_dict[z_key].keys() for z_key in z_index_dict.keys()]+\
                                    [t_dict[t_key].keys() for t_key in t_index_dict.keys()]).flatten())))
        
        num_per = len(all_periods)
        
        # fill arrays using the period key from all_periods
        z_arr = np.zeros((num_per, 2, 2), dtype=np.complex)
        z_err_arr = np.zeros((num_per, 2, 2), dtype=np.float)
        
        t_arr = np.zeros((num_per, 1, 2), dtype=np.complex)
        t_err_arr = np.zeros((num_per, 1, 2), dtype=np.float)
        
        for p_index, per in enumerate(all_periods):
            for z_key in sorted(z_index_dict.keys()):
                kk = z_index_dict[z_key][0]
                ll = z_index_dict[z_key][1]
                z_value = z_dict[z_key][per][0]+1j*z_dict[z_key][per][1]
                z_arr[p_index, kk, ll] = z_value
                z_err_arr[p_index, kk, ll] = z_dict[z_key][per][2]
            for t_key in sorted(t_index_dict.keys()):
                kk = t_index_dict[t_key][0]
                ll = t_index_dict[t_key][1]
                t_value = t_dict[t_key][per][0]+1j*t_dict[t_key][per][1]
                t_arr[p_index, kk, ll] = t_value
                t_err_arr[p_index, kk, ll] = t_dict[t_key][per][2]
        
        # put the results into mtpy objects
        freq = 1./np.array(all_periods)    
        self.Z = mtz.Z(z_arr, z_err_arr, freq)
        self.Tipper = mtz.Tipper(t_arr, t_err_arr, freq)    

#==============================================================================
# Write edi file from birrp outputs
#==============================================================================
def write_edi_file_from_birrp(station, birrp_dir, survey_config_fn, 
                              birrp_config_fn=None, copy_path=None):
                                  
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
    if not os.path.isdir(birrp_dir):
        raise mtex.MTpyError_inputarguments('Could not find {0}, check path'.format(birrp_dir))
        
    if not os.path.isfile(survey_config_fn):
        raise mtex.MTpyError_inputarguments('Could not find {0}, check path'.format(survey_config_fn))
        
    if birrp_config_fn is not None:
        if not os.path.isfile(birrp_config_fn):
            raise mtex.MTpyError_inputarguments('Could not find {0}, check path'.format(birrp_config_fn))
            
    else:
        birrp_config_fn = os.path.join(birrp_dir, 
                                       '{0}_birrp_params.cfg'.format(station))
        if not os.path.isfile(birrp_config_fn):
            print 'Could not find birrp_config_fn'
            birrp_config_fn = None
            
    # read in survey information        
    survey_cfg_dict = mtcfg.read_survey_configfile(survey_config_fn)
    
    if not station in survey_cfg_dict:
        raise mtex.MTpyError_config_file('Could not find information for {0} in {1}'.format(station, survey_config_fn))
                                       
        