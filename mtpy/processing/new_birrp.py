# -*- coding: utf-8 -*-
"""
Created on Tue Sep 20 14:33:20 2016

@author: jrpeacock
"""

import numpy as np

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
                     'should be 0 or 1 not {0}'.format(self.nrr)
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
            
class BIRRP_Parameter_Error(Exception):
    pass