# -*- coding: utf-8 -*-
"""
Created on Mon May 03 13:44:51 2010

@author: Jared Peacock
"""

import numpy as np
import os

import matplotlib.pyplot as plt
import matplotlib.pylab as pylab
from matplotlib.ticker import MultipleLocator,FormatStrFormatter
import matplotlib.gridspec as gridspec
from matplotlib.patches import Ellipse
from matplotlib.colors import LinearSegmentedColormap,Normalize
from matplotlib.colorbar import *
import mtpy.utils.latlongutmconversion as utm2ll

#make a custom colormap to use for plotting
ptcmapdict={'red':((0.0,1.0,1.0),(1.0,1.0,1.0)),
            'green':((0.0,0.0,1.0),(1.0,0.0,1.0)),
            'blue':((0.0,0.0,0.0),(1.0,0.0,0.0))}
ptcmap=LinearSegmentedColormap('ptcmap',ptcmapdict,256)

#resistivity tensor map for calcluating apparent resistivity
rtcmapdict={'red':((0.0,1.0,1.0),(0.5,1.0,1.0),(1.0,0.0,0.0)),
            'green':((0.0,0.0,0.0),(0.5,1.0,1.0),(1.0,0.0,0.0)),
            'blue':((0.0,0.0,0.0),(0.5,1.0,1.0),(1.0,1.0,1.0))}
rtcmap=LinearSegmentedColormap('rtcmap',rtcmapdict,256)

#spacing for text files
tsp='   '

class Edi(object):
    """
    This is a data type for edi files to read, write and rewrite any edi file,
    well, at least that is the plan.
        
    """
    
    def __init__(self,edifilename,ncol=5):
        self.edifn=edifilename
        self.header={'text':[]}
        self.info={'text':[]}
        self.measurement={'text':[]}
        self.lat=0
        self.lon=0
        self.elevation=0
        self.frequency=0
        self.z=0
        self.zvar=0
        self.zrot=0
        self.tipper=0
        self.tippervar=0
        self.ncol=ncol

    def readEDI(self,verbose=None):
        """
        readEDI will read in the information of the edi file
        
        Arguments:
        ----------
            **Edi.edifn** : string
                            full path to .edi file to be read in.
            
            **verbose** : boolean (True,False)
                          * True to print out diagnositic messages
                          * False to not print out messages
                            
        Returns:
        --------
            data type Edi
            
        :Example: ::
            
            >>> import mtpy.core.z as Z
            >>> edi1 = Z.Edi(edifile)
            >>> edi1.readEDI(verbose=True)
        """
        #open up file
        edifid=file(self.edifn,'r')
        #read all the lines from that file
        edilines=edifid.readlines()
        edifid.close()
        
        #---Read in the information about the measurement---------
        gdict={'head':self.header,
               'info':self.info,
               'definemeas':self.measurement,
               'mtsect':{'text':[]},
               'emapsect':{'text':[]},
               'comp':{}}
        
        #get indexes of important information
        edict={}
        for ii,eline in enumerate(edilines):

            if eline.find('FREQ')>0:
                edict['freq']=ii
                continue
            if eline.find('SPECTRASECT')>0:
                edict['spectra']=ii
                continue
            if eline.find('>ZROT')>=0:
                if not edict.has_key('zrot'):
                    edict['zrot']=ii
                continue
            #elif eline.find('IMPEDANCE')>0:
            #    if eline.find('ROTATION')>0:
            #        edict['zrot']=ii
            #    else:
            #        edict['z']=ii
            if eline.find('>Z')>=0:
                if not edict.has_key('z'):
                    edict['z']=ii
                continue

            if eline.find('TIPPER')>0:
                if eline.find('ROTATION')>0:
                    edict['trot']=ii
                else:
                    edict['tipper']=ii

        #-------Read in header information------------------------
        ii=0
        while type(ii) is int:
        
            eline=edilines[ii]
            if eline.find('SPECTRA')>=0:      
                ii=None
            elif eline.find('IMPEDANCE')>=0:      
                ii=None
            elif eline.find('FREQUENCIES')>=0:     
                ii=None
            elif eline.find('>END')>=0:
                ii=None
                continue
            else:
                #get the block header
                if eline.find('>')==0 or eline.find('>')==1:
                    es=eline.find('>')
                    #For definemeas
                    if eline.find('=')>0 and eline.find('ACQCHAN')==-1:
                        if eline.find('INFO')>0:
                            gkey=eline[es+1:].strip().split()[0].lower()
                        else:
                            gkey=eline[es+2:].strip().lower()
                    #for a comment block
                    elif eline.find('!')>0:
                        pass
                    #for the measurement channels
                    elif eline.find('ACQCHAN')>0:
                        gkey='comp'
                        eline=eline.strip().split()
                        mid=eline[1].split('=')[1]
                        mtype=eline[2].split('=')[1]
                        gdict['comp'][mid]=mtype
                    #for the normal block header
                    else:
                        gkey=eline[eline.find('>')+1:].strip().lower()
                #for each block header put the information into a dictionary
                else:
                    eline=eline.strip()
                    #for a line with that can be separated by an =
                    if eline.find('=')>0 and eline.find('CHTYPE')==-1 and \
                        eline.find('mV')==-1:
                        eline=eline.split('=')
                        gdict[gkey][eline[0]]=eline[1].replace('"','')
                    #all other lines
                    elif eline.find('=')==-1 and eline.find(':')>0:
                        gdict[gkey]['text'].append(eline.strip())
                ii+=1
        
        #get latitude from the header
        try:
            latstr=gdict['head']['LAT']
            
        #if it isn't in the header it should be in the definemeas
        except KeyError:
            try:
                latstr=gdict['definemeas']['REFLAT']
            #if it isn't there then need to find it else where
            except KeyError:
                print 'Did not find Latitude'
                
        #change from hh:mm:ss to decimal deg
        if latstr.find(':')>=0:
            latlst=latstr.split(':')
            latstr=str(int(latlst[0]))+'.'+str(float(latlst[1])/60
                +float(latlst[2])/3600)[2:]
            self.lat=float(latstr)
        else:
            self.lat=float(latstr)
        
        #get longitude from header
        try:
            lonstr=gdict['head']['LONG']
            
        #if it isn't in the header it should be in the definemeas
        except KeyError:
            try:
                lonstr=gdict['definemeas']['REFLONG']
            except KeyError:
                print 'Did not find Longitude'
                
        #change from hh:mm:ss to decimal deg
        if lonstr.find(':')>=0:
            lonlst=lonstr.split(':')
            lonstr=str(int(lonlst[0]))+'.'+str(float(lonlst[1])/60
                +float(lonlst[2])/3600)[2:]
            self.lon=float(lonstr)
        else:
            self.lon=float(lonstr)
        
        #get elevation
        try:
            self.elevation=float(gdict['head']['ELEV'])
        except KeyError:
            if verbose:
                print 'Did not find elevation for ',self.edifn
        
        #======================================================================
        #         Get frequency, impedance and tipper
        #======================================================================
        #-------------------Get Frequencies------------------------------------        
        try:
            ii=edict['freq']
            #get number of frequencies
            nf=int(edilines[ii].strip().split('//')[1])
            #initialize some arrays    
            self.frequency=np.zeros(nf)
            
            #get frequencies
            kk=0
            ii+=1
            while type(kk) is int:
                if edilines[ii].find('>')>=0 or edilines[ii].find('*')>0 :
                    kk=None
                else:
                    eline=edilines[ii].strip().split()
                    if kk==0:
                        self.ncol=len(eline)
                    for nn,estr in enumerate(eline):
                        try:
                            self.frequency[kk*(self.ncol-1)+nn+kk]=float(estr)
                        except IndexError:
                            pass
                    kk+=1
                    ii+=1
        except KeyError:
            if verbose:
                print 'Did not find frequencies in ',self.edifn

            
        #-------------Get impedance from spectra---------------------------
        try:
            ii=edict['spectra']+1
            #make a dictionary with a list to put the component order into
            gdict['spectraset']={'comporder':{}}
            
            cc=0
            #get header information about how the spectra is stored
            while type(ii) is int:
                eline=edilines[ii]
                #stop at the first spectra block
                if eline.find('>SPECTRA')==0:
                    jj=ii
                    ii=None
                #get the spectraset information
                else:
                    gkey='spectraset'
                    #general information
                    if eline.find('=')>0:
                        eline=eline.strip().split('=')
                        gdict[gkey][eline[0]]=eline[1].replace('"','')
                    #get component order
                    else:
                        try:
                            estr=eline.strip()
                            try:
                                gdict[gkey]['comporder'][gdict['comp'][estr]]
                                gdict[gkey]['comporder'][gdict['comp'][estr]+'R']=cc
                            except KeyError:
                                gdict[gkey]['comporder'][gdict['comp'][estr]]=cc
                                
                            cc+=1
                        except KeyError:
                            pass
                    ii+=1
            #defind values of number of frequencies and number of channels
            nf=int(gdict['spectraset']['NFREQ'])
            nc=int(gdict['spectraset']['NCHAN'])
        
            #set some empty arrays to put stuff into    
            self.z=np.zeros((nf,2,2),dtype='complex')
            self.zvar=np.zeros((nf,2,2))
            self.tipper=np.zeros((nf,2),dtype='complex')
            self.tippervar=np.zeros((nf,2))
            self.frequency=np.zeros(nf)
            bw=np.zeros(nf)
            avgt=np.zeros(nf)
            self.zrot=np.zeros(nf)
            self.trot=np.zeros(nf)
            
            kk=0
            while type(kk) is int:
                eline=edilines[jj]
                #get information from the spectra line block
                if eline.find('>SPECTRA')==0:
                    eline=eline.strip().split()
                    for ee in eline:
                        estr=ee.split('=')
                        if estr[0]=='FREQ':
                            self.frequency[kk]=float(estr[1])
                        elif estr[0]=='ROTSPEC':
                            self.zrot[kk]=float(estr[1])
                        elif estr[0]=='BW':
                            bw[kk]=float(estr[1])
                        elif estr[0]=='AVGT':
                            avgt[kk]=float(estr[1])
                        else:
                            pass
                    #make and empty array to put spectra into for easy 
                    #calculation
                    spectra=np.zeros((nc,nc))
                    kk+=1
                    jj+=1
                #stop once all the spectra have been read
                elif eline.find('>')==0 and eline.find('SPECTRA')==-1:
                    kk=None
                #put the spectra values into the pre-defined array
                else:
                    for ll in range(nc):
                        eline=edilines[jj].strip().split()
                        jj+=1
                        for nn in range(nc):
                            spectra[ll,nn]=float(eline[nn])       
                    #get spectra in a useable order such that all values are
                    #complex in the upper right hand triangle
                    spect=np.zeros((nc,nc),dtype='complex')
                    for ll in range(nc):
                        for nn in range(ll,nc):
                            spect[ll,nn]=spectra[nn,ll]-1j*spectra[ll,nn]
                            
                    
                    #calculate the impedance from the spectra
                    cdict=gdict['spectraset']['comporder']
                    bx=cdict['HX']            
                    by=cdict['HY']
                    ex=cdict['EX']            
                    ey=cdict['EY']
                    try:
                        bz=cdict['HZ']
                    except KeyError:
                        print 'No HZ'
                    try: 
                        bxr=cdict['HXR']            
                        byr=cdict['HYR']
                    except KeyError:
                        print 'No remote reference'
                        
                    #get impedance from the spectra
                    zdet=(spect[bx,bxr]*spect[by,byr])-\
                            (spect[bx,byr]*spect[by,bxr])
                    self.z[kk-1,0,0]=((spect[ex,bxr]*spect[by,byr])-
                                (spect[ex,byr]*spect[by,bxr]))/zdet
                    self.z[kk-1,0,1]=((spect[ex,byr]*spect[bx,bxr])-
                                (spect[ex,bxr]*spect[bx,byr]))/zdet
                    self.z[kk-1,1,0]=((spect[ey,bxr]*spect[by,byr])-
                                (spect[ey,byr]*spect[by,bxr]))/zdet
                    self.z[kk-1,1,1]=((spect[ey,byr]*spect[bx,bxr])-
                                (spect[ey,bxr]*spect[bx,byr]))/zdet
                                
                    #get tipper from the spectra
                    self.tipper[kk-1,0]=((spect[bx,bz]*spect[by,by].real)-
                                    (spect[by,bz]*spect[bx,by]))/zdet
                    
                    #need to put in conjugate because of coordinate system
                    self.tipper[kk-1,1]=((spect[by,bz].conj()*spect[bx,bx].real)-
                                    (spect[bx,bz].conj()*spect[bx,by]))/zdet
        
        except KeyError:
            if verbose:
                print 'Did not find spectra information for ',self.edifn
        
        #--------------Get Impedance Rotation angles----------------------
        try:
            ii=edict['zrot']
            
            #get number of frequencies
            nf=int(edilines[ii].strip().split('//')[1])
            
            #initialize some arrays    
            self.zrot=np.zeros(nf)
            
            #get rotation angles
            kk=0
            ii+=1
            while type(kk) is int:
                if edilines[ii].find('>')>=0 or edilines[ii].find('*')>0:
                    kk=None
                else:     
                    eline=edilines[ii].strip().split()
                    if kk==0:
                        self.ncol=len(eline)
                    for nn,estr in enumerate(eline): 
                        try:
                            self.zrot[kk*(self.ncol-1)+nn+kk]=float(estr)
                        except IndexError:
                            pass
                    kk+=1
                    ii+=1
        
        except KeyError:
            if verbose:
                print 'Did not find impedance rotation block for ',self.edifn
        
        #--------------Get impedance--------------------------------------                           
        try:
            ii=edict['z']
            
            #define a dictionary of indecies to put information into z array
            zdict=dict([('ZXX',(0,0)),('ZXY',(0,1)),('ZYX',(1,0)),
                        ('ZYY',(1,1))])
            
            #initialize some arrays
            nf=int(edilines[ii].strip().split('//')[1])
            self.z=np.zeros((nf,2,2),dtype='complex')
            self.zvar=np.zeros((nf,2,2))
            self.tipper=np.zeros((nf,2),dtype='complex')
            self.tippervar=np.zeros((nf,2))
            self.trot=np.zeros(nf)
            try:
                edict['zrot']
            except KeyError:
                self.zrot=np.zeros(nf)
        
            kk=0
            while type(kk) is int:
                eline=edilines[ii]
                if eline.find('>')==0 or eline.find('>')==1:
                    es=eline.find('>')
                    if eline.find('Z')==-1 and eline.find('IMP')==-1:
                            kk=None
                            jj=ii
                    else:
                        eline=eline.strip().split()
                        zkey=eline[0][1:4]
                        zcomp=eline[0][4:]
                        if zkey!='ZRO' and zkey.find('!')==-1:
                            z0=zdict[zkey][0]
                            z1=zdict[zkey][1]
                            
                        kk=0
                        ii+=1
                else:
                    eline=eline.strip().split()
                    if zkey=='ZRO':
                        for nn,estr in enumerate(eline):
                            try:
                                self.zrot[kk*(self.ncol-1)+nn+kk]=float(estr)
                            except IndexError:
                                pass
                    elif zcomp=='R':
                        for nn,estr in enumerate(eline):
                            try:
                                self.z[kk*(self.ncol-1)+nn+kk,z0,z1]=float(estr)
                            except IndexError:
                                pass
                    elif zcomp=='I':
                        for nn,estr in enumerate(eline):
                            try:
                                self.z[kk*(self.ncol-1)+nn+kk,z0,z1]=\
                                        self.z[kk*(self.ncol-1)+nn+kk,z0,z1]+\
                                        1j*float(estr)
                            except IndexError:
                                pass
                    elif zcomp=='.VAR':
                        for nn,estr in enumerate(eline):
                            try:
                                self.zvar[kk*(self.ncol-1)+nn+kk,z0,z1]=\
                                                                    float(estr)
                            except IndexError:
                                pass
                    ii+=1
                    kk+=1
        except KeyError:
            if verbose:
                print 'Did not find impedance information for ',self.edifn
        
        #--------------Get Tipper Rotation angles----------------------
        try:
            ii=edict['trot']+1
            
            #get number of frequencies
            nf=int(edilines[ii].strip().split('//')[1])
            
            #initialize some arrays    
            self.trot=np.zeros(nf)
            
            #get rotation angles
            kk=0
            ii+=1
            while type(kk) is int:
                if edilines[ii].find('!')>0 or edilines[ii].find('*')>0:
                    kk=None
                else:
                    eline=edilines[ii].strip().split()
                    if kk==0:
                        self.ncol=len(eline)
                    for nn,estr in enumerate(eline):
                        try:
                            self.trot[kk*(self.ncol-1)+nn+kk]=float(estr)
                        except IndexError:
                            pass
                    kk+=1
                    ii+=1
        
        except KeyError:
            if verbose:
                print 'Did not find Tipper rotation block for ',self.edifn

        #-----Get Tipper Information-----------
        try:
            ii=edict['tipper']+1
            tdict=dict([('TX',0),('TY',1)])
            
            kk=0
            while type(kk) is int:
                eline=edilines[ii]
                if eline.find('>')==0 or eline.find('>')==1:
                    es=eline.find('>')
                    if eline.find('T',es,es+2)==-1:
                            kk=None
                            jj=ii
                    else:
                        eline=eline.strip().split()
                        tkey=eline[0][1:3]
                        tcomp=eline[0][3:]
                        if tkey!='TR' and tkey.find('!')==-1:
                            z0=tdict[tkey]
                            
                        kk=0
                        ii+=1
                else:
                    eline=eline.strip().split()
                    if tkey=='TR':
                        for nn,estr in enumerate(eline):
                            try:
                                self.trot[kk*(self.ncol-1)+nn+kk]=float(estr)
                            except IndexError:
                                pass
                    if tcomp=='R' or tcomp=='R.EXP':
                        for nn,estr in enumerate(eline):
                            try:
                                self.tipper[kk*(self.ncol-1)+nn+kk,z0]=float(estr)
                            except IndexError:
                                pass
                    elif tcomp=='I' or tcomp=='I.EXP':
                        for nn,estr in enumerate(eline):
                            try:
                                self.tipper[kk*(self.ncol-1)+nn+kk,z0]+=\
                                                                 1j*float(estr)
                            except IndexError:
                                pass
                    elif tcomp=='VAR.EXP' or '.VAR':
                        for nn,estr in enumerate(eline):
                            try:
                                self.tippervar[kk*(self.ncol-1)+nn+kk,z0]=\
                                                                float(estr)
                            except IndexError:
                                pass
                    ii+=1
                    kk+=1

        except KeyError:
            if verbose:
                print 'Did not find Tipper information for ',self.edifn 
        
        self.period=1./self.frequency        
        if self.period[0]>self.period[-1]:
            self.period=self.period[::-1]
            self.frequency=self.frequency[::-1]
            self.z=self.z[::-1,:,:]
            self.zvar=self.zvar[::-1,:,:]
            self.tipper=self.tipper[::-1,:]
            self.tipvar=self.tippervar[::-1,:]
                        
            print 'Flipped to descending frequency'
    
    def rewriteedi(self,znew=None,zvarnew=None,freqnew=None,newfile='y',
               tipnew=None,tipvarnew=None,thetar=0,ext='dr'):
        """
        rewriteedi(edifile) will rewrite an edifile say if it needs to be 
        rotated or distortion removed.
    
        Arguments:
        ----------
            **Edi.edifile** : string
                              full path to edifile to be rewritten
            
            **znew** : complex np.array with shape (nf,2,2) 
                        [[zxx,zxy],[zyx,zyy]]
                       impedance tensor if a new one has been created
            
            **zvarnew** : real np.array with shape (nf,2,2)
                          [[zxx_err,zxy_err],[zyx_err,zyy_err]]
                          errors in impedance tensor if a new one has been 
                          created
           
            **freqnew** : np.array with shape (nf)
                          new frequency list if one has been created
            
            **newfile** : string ('y','n')
                          * 'y' to create a new edi file
                          * 'n' to rewrite the existing edifile
            
            **tipnew** : complex np.array with shape (nf,2)
                         [[tx,ty]]
                         new tipper array
            
            **tipvarnew** : real np. array with shape (nf,2)
                            [[tx_err,ty_err]]
                            rnew tipper error array
            
            **thetar** : float (angle in degrees)
                         rotation angle clockwise (N=0, E=90)
            
            **ext** : string
                      extension on the end of the file name before .edi
        
        Returns:
        --------
            **Edi.nedifn** : string
                            full path to rewritten edi file
                            dirpath(edifile)+basename(edifile)+ext
                            
        :Example: ::
            
            >>> import mtpy.core.z as Z
            >>> edi1 = Z.Edi(edifile)
            >>> edi1.rewriteedi(thetar=65,ext='Rot65')
            Made directory:  c:\Rot65
            Made file:  c:\Rot65\pb01cRot65.edi
        
        """
    
        #get direcotry path make one if not there    
        dirpath=os.path.dirname(self.edifn)
        if newfile=='y':
            drdirpath=os.path.join(dirpath,ext.upper())
            if not os.path.exists(drdirpath):
                os.mkdir(drdirpath)
                print 'Made directory: ',drdirpath
            else:
                pass
            #copy files into a common directory
            drpath=os.path.dirname(dirpath)
            count=1
            while count<10:
                drpath=os.path.dirname(drpath)
                for folder in os.listdir(drpath):
                    if folder.find('EDIfiles')>=0:
                        edifolderyn='y'
                        drcopypath=os.path.join(drpath,'EDIfiles')
                        if os.path.exists(drcopypath):
                            if not os.path.exists(os.path.join(drcopypath,
                                                               ext.upper())):
                                os.mkdir(os.path.join(drcopypath,ext.upper()))
                                print 'Made directory ',os.path.join(drcopypath,
                                                                     ext.upper())
                    else:
                        edifolderyn='n'
                        drcopypath=os.path.dirname(drdirpath)
                count+=1
            
            #get new file name
            nedibn=os.path.basename(self.edifn)[:-4]+ext.lower()+'.edi'
            #open a new edifile
            newedifid=open(os.path.join(drdirpath,nedibn),'w')
            
        else:
            nedibn=os.path.basename(self.edifn)[:-4]+'c.edi'
            newedifid=open(os.path.join(dirpath,nedibn),'w')
            edifolderyn='n'
        
        #if there is no znew then rotate the data
        if znew==None:
            #read in the data
            self.readEDI()
            #rotate data if desired
            if thetar!=0:
                self.znew=np.zeros_like(self.z)
                self.zvarnew=np.zeros_like(self.zvar)
                self.tippernew=np.zeros_like(self.tipper)
                self.tippervarnew=np.zeros_like(self.tippervar)
                
                #convert thetar into radians
                if abs(thetar)>2*np.pi:
                    self.thetar=thetar*(np.pi/180)
                else:
                    self.thetar=thetar
                print self.thetar
                #make rotation matrix
                rotmatrix=np.array([[np.cos(self.thetar), np.sin(self.thetar)],
                             [-np.sin(self.thetar), np.cos(self.thetar)]])
                trotmatrix=np.array([np.cos(self.thetar),np.sin(self.thetar)])
                
                #fill in rotation for impedances (zrot)
                self.zrot=np.zeros_like(self.frequency)
                self.zrot[:]=thetar
                
                #rotate the data
                for rr in range(len(self.frequency)):
                    self.znew[rr]=np.dot(rotmatrix,np.dot(self.z[rr],
                                         rotmatrix.T))
                    self.zvarnew[rr]=np.dot(rotmatrix,np.dot(self.zvar[rr],
                                            rotmatrix.T))
                    self.tippernew[rr]=np.dot(trotmatrix,
                                            np.dot(self.tipper[rr],
                                                   trotmatrix.T))
                    self.tippervarnew[rr]=np.dot(trotmatrix,
                                                 np.dot(self.tippervar[rr],
                                                        trotmatrix.T))
            #otherwise read in the new Z
            else:
                self.znew=self.z
                self.zvarnew=self.zvar
                self.tippernew=self.tipper
                self.tippervarnew=self.tippervar
                #fill in rotation for impedances (zrot)
                self.zrot=np.zeros_like(self.frequency)
                self.zrot[:]=thetar
                self.thetar=0.0
        else:
            self.znew=znew
            self.zvarnew=zvarnew
            self.tippernew=self.tipper
            self.tippervarnew=self.tippervar
            #fill in rotation for impedances (zrot)
            self.zrot=np.zeros_like(self.frequency)
            self.zrot[:]=thetar
            self.thetar=0.0
        
        #open existing edi file
        edifid=open(self.edifn,'r')
        
        #read existing edi file and find where impedances start and number of 
        #frequencies
        edilines=edifid.readlines()
        spot=[]
        for ii,line in enumerate(edilines):
            if line.find('>FREQ')>=0:
                spot.append(ii)
            if line.find('IMPEDANCES')>=0:
                spot.append(ii)
            if line.find('SPECTRA')>=0:
                spot.append(ii)
        
        #write lines until the frequency list or spectra
        for ll,line in enumerate(edilines[0:spot[0]-1]):
            newedifid.write(line)
        
        newedifid.write('>!****FREQUENCIES****!'+'\n')
        #write in frequencies
        #if there is a new frequency list
        if freqnew!=None:
            #get number of frequencies    
            nfreq=len(freqnew)
            
            #get order of frequencies        
            if freqnew[0]<freqnew[-1]:
                order='INC'
            else:
                order='DEC'
            
            #write frequency header
            newedifid.write('>FREQ'+tsp+'NFREQ='+str(int(nfreq))+tsp+\
                            'ORDER='+order+tsp+'// '+str(int(nfreq))+'\n')
            for kk in range(int(nfreq)):
                newedifid.write(tsp+'{0:+.6E}'.format(freqnew[kk]))
                if np.remainder(float(kk)+1,5.)==0:
                    newedifid.write('\n')
            newedifid.write('\n')

        #from the start of file to where impedances start write from existing 
        #edifile
        else:
            #get order of frequencies        
            if self.frequency[0]<self.frequency[-1]:
                order='INC'
            else:
                order='DEC'
                
            nfreq=len(self.frequency)
            #write frequency header
            newedifid.write('>FREQ'+tsp+'NFREQ='+str(int(nfreq))+tsp+'ORDER='+\
                            order+tsp+'// '+str(int(nfreq))+'\n')
            for kk in range(int(nfreq)):
                newedifid.write(tsp+'{0:+.6E}'.format(self.frequency[kk]))
                if np.remainder(float(kk)+1,5.)==0:
                    newedifid.write('\n')
            newedifid.write('\n')
        
        #write in rotation block for impedance tensor
        newedifid.write('>!****IMPEDANCE ROTATION ANGLES****!\n')
        newedifid.write('>ZROT // {0}\n'.format(int(nfreq)))
        for ss in range(nfreq):
            newedifid.write(tsp+'{0:+.4f}'.format(self.thetar*180/np.pi))
            if np.remainder(ss+1.0,5.)==0:
                newedifid.write('\n')
        if np.remainder(ss+1.,5.)!=0:
            newedifid.write('\n')
        
        #write in impedance information        
        newedifid.write('>!****IMPEDANCES****!'+\
                       ' Rotated {0:.0f} clockwise\n'.format(self.thetar*\
                                                             180/np.pi))
        
        #create an implst to make code simpler
        implst=[['ZXXR',0,0],['ZXXI',0,0],['ZXX.VAR',0,0],['ZXYR',0,1],
                ['ZXYI',0,1],['ZXY.VAR',0,1],['ZYXR',1,0],['ZYXI',1,0], 
                ['ZYX.VAR',1,0],['ZYYR',1,1],['ZYYI',1,1],['ZYY.VAR',1,1]]
        #write new impedances and variances
        for jj,imp in enumerate(implst):
            mm=imp[1]
            nn=imp[2]
            newedifid.write('>'+imp[0]+' ROT=ZROT // '+str(int(nfreq))+'\n')
            if imp[0].find('.')==-1:
                if imp[0].find('R')>=0:
                    for kk in range(int(nfreq)):
                        newedifid.write(tsp+'{0:+.6E}'.format(self.znew.real[kk,mm,nn]))
                        if kk>0:
                            if np.remainder(float(kk)+1,5.)==0:
                                newedifid.write('\n')
                        else:
                            pass
                elif imp[0].find('I')>=0:
                    for kk in range(int(nfreq)):
                        newedifid.write(tsp+'{0:+.6E}'.format(self.znew.imag[kk,mm,nn]))
                        if kk>0:
                            if np.remainder(float(kk)+1,5.)==0:
                                newedifid.write('\n')
                        else:
                            pass
            else:
                for kk in range(int(nfreq)):
                    newedifid.write(tsp+'{0:+.6E}'.format(self.zvarnew[kk,mm,nn]))
                    if kk>0:
                        if np.remainder(float(kk)+1,5.)==0:
                            newedifid.write('\n')
                    else:
                        pass
            if np.remainder(kk+1.,5.)!=0:
                newedifid.write('\n')
        
        newedifid.write('>!****TIPPER****!'+'\n')
        tiplst=[['TXR',0],['TXI',0],['TX.VAR',0],['TYR',1],['TYI',1],
                ['TY.VAR',1]]

        for jj,tip in enumerate(tiplst):
            mm=tip[1]
            newedifid.write('>'+tip[0]+' // '+str(int(nfreq))+'\n')
            if tip[0].find('.')==-1:
                if tip[0].find('R')>=0:
                    for kk in range(int(nfreq)):
                        newedifid.write(tsp+'{0:+.6E}'.format(self.tippernew.real[kk,mm]))
                        if kk>0:
                            if np.remainder(float(kk)+1,5.)==0:
                                newedifid.write('\n')
                        else:
                            pass
                    if np.remainder(kk+1.,5.)!=0:
                        newedifid.write('\n')
                elif tip[0].find('I')>=0:
                    for kk in range(int(nfreq)):
                        newedifid.write(tsp+'{0:+.6E}'.format(self.tippernew.imag[kk,mm]))
                        if kk>0:
                            if np.remainder(float(kk)+1,5.)==0:
                                newedifid.write('\n')
                        else:
                            pass
                    if np.remainder(kk+1.,5.)!=0:
                        newedifid.write('\n')
            else:
                for kk in range(int(nfreq)):
                    newedifid.write(tsp+'{0:+.6E}'.format(self.tippervarnew.real[kk,mm]))
                    if kk>0:
                        if np.remainder(float(kk)+1,5.)==0:
                            newedifid.write('\n')
                    else:
                        pass
                if np.remainder(kk+1.,5.)!=0:
                        newedifid.write('\n')
        newedifid.write('\n')
        newedifid.write('>END')
        
        edifid.close()
        newedifid.close()
        #copy file to a common folder
        if edifolderyn=='y':
            if newfile=='y':
                shutil.copy(os.path.join(drdirpath,newedifile),
                            os.path.join(drcopypath,ext.upper(),newedifile))
                print 'Copied new edifile to: ',os.path.join(drcopypath,
                                                             ext.upper(),
                                                             newedifile)
            else:
                shutil.copy(os.path.join(dirpath,newedifile),
                        os.path.join(drcopypath,newedifile))
                print 'Copied new edifile to: ',os.path.join(drcopypath,
                                                             newedifile)
        else:
            pass
        if newfile=='y':
            self.nedifn=os.path.join(drdirpath,nedibn)
            
        else:
            self.nedifn=os.path.join(dirpath,nedibn)
            
        print 'Made file: ', self.nedifn

class Z(Edi):
    """
    Z is a data type to deal with edifiles and manipulate them to do 
    informative characterization.
    
    Inherits the class Edi so edi files can be read and rewritten.
            
    """
    def __init__(self,edifn,ncol=5):
        super(Edi,self).__init__()
        
        #define some parameters to be filled
        self.edifn=edifn
        self.header={'text':[]}
        self.info={'text':[]}
        self.measurement={'text':[]}
        self.lat=0
        self.lon=0
        self.elevation=0
        self.frequency=0
        self.z=0
        self.zvar=0
        self.zrot=0
        self.tipper=0
        self.tippervar=0
        self.ncol=ncol
        
        #get information from edifile as a datatype Edi
        Edi.readEDI(self)
        
        #define some parameters needed from older version
        self.period=1./self.frequency
        self.station=self.header['DATAID']
        
                    
    def getInvariants(self,thetar=0):
        """
        Calculate the invariants of the impedance tensor according to 
        Weaver et al. [2003].
        
        Arguments:
        ----------
            **thetar** : float (angle in degrees)
                         rotation angle clockwise positive
        
        Returns:
        --------
            **Zinvariants** : data type Zinvariants
            
        :Example: ::
            
            >>> z1 = Z.Z(edifile)
            >>> zinv = z1.getInvariants(thetar=0)
        """
        
        return Zinvariants(self.z,rotz=thetar)
    
    def getPhaseTensor(self,rotate=180,thetar=0):
        """
        Calculate phase tensor elements from the impedance tensor following 
        Caldwell et al. [2004].
        
        Arguments:
        ----------
            **rotate** : [ 90 | 180 | 270 ]  
                         Caldwell et al, [2004] assume the coordinate axis is 
                         Y North and X East.  If the data is in X North and Y 
                         East than rotation = 180. *Default* is 180
                         
            **thetar** : float (angle in degrees)
                         rotation angle clockwise positive assuming 0 is North.
            
        Returns:
        --------
            **PhaseTensor** : data type Phase Tensor

        :Example: ::
            
            >>> z1 = Z.Z(edifile)
            >>> pt = z1.getPhaseTensor()
        
        """
        
        pt=PhaseTensor(self.z,self.zvar,rotate=rotate,rotz=thetar)
        return pt

    def getTipper(self,thetar=0):
        """
        Get tipper information and return a type with attributes:
            
        Arguments:
        ----------
            **thetar** : float (angle in degrees)
                         rotation angle clockwise positive assuming 0 is North. 
        
        Returns:
        --------
            **Tipper** : data type Tipper
            
        :Example: ::
            
            >>> z1 = Z.Z(edifile)
            >>> tip = z1.getTipper(thetar=10)
            
        """
        
           
        return Tipper(self.tipper,rott=thetar)
   
    def removeDistortion(self,thetar=0):
        """
        removeDistortion(self) will remove the distortion from the impedance
        tensor as prescribed by Bibby et al. [2005] for the 1D case.  
        
        Argumens:
        ---------
            **thetar** : float (angle in degrees)
                         rotation angle clockwise positive assuming 0 is North. 
        
        Returns:
        --------
            **Z.D** : real np.array (2,2)
                    estimated distortion tensor
                    
            **Z.newedifn** : string
                           full path to new edifile the edifile as 
                           station+dr.edi.
                           
        :Example: ::
            
            >>> z1 = Z.Z(edifile)
            >>> z1.removeDistortion()
            >>> z1.D
            np.array([[0.01,0.9],[0.98,0.05]])
            
        
        """

        z=np.array(self.z)
        zvar=np.array(self.zvar)
        #calculate ellipticity to figure out the distortion from 1d structure
        pt=Z.getPhaseTensor(self,thetar=thetar)
        zd=[]
        #rotated identity matrix for static shift removal
        im=np.array([[0,-1],[1,0]])
        for ii,ellip in enumerate(pt.ellipticity):
            if ellip<.1:
                X=z[ii].real
                Y=z[ii].imag
                #calculate static shift factor
                gx=np.sqrt(abs(np.linalg.det(X)))
                gy=np.sqrt(abs(np.linalg.det(Y)))
                #remove static shift from real and imaginary parts
                #need to use the dot for matric multiplication
                Xs=np.dot(X,im)/gx
                Ys=np.dot(Y,im)/gy
                #append real and imaginary parts sequentially
                zd.append(Xs)
                zd.append(Ys)
            else:
                pass
        if len(zd)==0:
            print 'There is no 1D structure for this station'
            self.D=np.zeros((2,2))
            self.nedifn=self.edifn
        else:
            #estimate distortion matrix 
            zd=np.array(zd)
            D=np.array([[zd[:,0,0].mean(),zd[:,0,1].mean()],
                         [zd[:,1,0].mean(),zd[:,1,1].mean()]])
            Dvar=np.array([[zd[:,0,0].std(),zd[:,0,1].std()],
                            [zd[:,1,0].std(),zd[:,1,1].std()]])
            Dinv=np.linalg.inv(D)
            #remove distortion as (D^-1)Z[ii]
            zdr=np.array([np.dot(np.linalg.inv(D),z[ii]) for ii in range(len(z))])
            
            #estimate errors
            zvardr=np.zeros_like(z)
            for jj in range(len(z)):
                X=z[jj].real
                Xvar=zvar[jj].real
                zvardr[jj,0,0]=np.sqrt((Dinv[0,0]*X[0,0]*
                            np.sqrt((Dvar[1,1]/Dinv[0,0])**2+
                            (Xvar[0,0]/X[0,0])**2))**2+(Dinv[0,1]*X[1,0]
                            *np.sqrt((Dvar[0,1]/Dinv[0,1])**2+
                            (Xvar[0,1]/X[0,1])**2))**2)
                zvardr[jj,0,1]=np.sqrt((Dinv[0,0]*X[0,1]*
                            np.sqrt((Dvar[1,1]/Dinv[0,0])**2+
                            (Xvar[0,1]/X[0,1])**2))**2+(Dinv[0,1]*X[1,1]
                            *np.sqrt((Dvar[0,1]/Dinv[0,1])**2+
                            (Xvar[1,1]/X[1,1])**2))**2)
                zvardr[jj,1,0]=np.sqrt((Dinv[1,0]*X[0,0]*
                            np.sqrt((Dvar[1,0]/Dinv[1,0])**2+
                            (Xvar[0,0]/X[0,0])**2))**2+(Dinv[1,1]*X[1,0]
                            *np.sqrt((Dvar[0,0]/Dinv[1,1])**2+
                            (Xvar[1,0]/X[1,0])**2))**2)
                zvardr[jj,1,1]=np.sqrt((Dinv[1,0]*X[1,0]*
                            np.sqrt((Dvar[1,0]/Dinv[1,0])**2+
                            (Xvar[0,1]/X[0,1])**2))**2+(Dinv[1,1]*X[1,1]
                            *np.sqrt((Dvar[1,1]/Dinv[1,1])**2+
                            (Xvar[1,1]/X[1,1])**2))**2)
            #make new edi file. need to flip zdr and zvardr back to normal order 
            #for edi file.
            self.rewriteedi(znew=zdr,zvarnew=zvardr.real)
            self.D=D

    def removeStaticShift(self,stol=.2,dm=1000,fspot=20):
        """
        Will remove static shift by calculating the median of respones of near
        by stations, within dm.  If the ratio of the station response to the 
        median on either side of 1+-stol then the impedance tensor for that 
        electric component will be corrected for static shift.  
        
        .. note:: **Usually only works with closely spaced (<1km) stations.**
        
        .. note:: It is important to have all the edifiles in one folder.
        
        Arguments:
        ----------
            **stol** : float (ratio tolerance).  
                        If there is no static shift the ratio between the 
                        response and the median response should be 1, but noise
                        and other factors can be present so a tolerance 
                        around 1 is assumed. *Default* is 0.2
                   
            **dm** : float (nearby station radius in meters)  
                     If there is no station within that radius then no static 
                     shift will be corrected for. *Default* is 1000m
                 
            **fspot** : int (the last index of frequencies to look at) 
                        So if you want to look for static shift in the first 
                        20 frequencies fspot=20.  *Default* is 20
        
        Returns:
        --------
            **Z.nedifn** : string
                           full path to new edifile
        
        :Example: ::
            
            #to estimate the static shift for a radius of 1500m for the first 
            #50 frequencies enter:
            >>> z1 = Z.Z(edifile)
            >>> z1.removeStaticShift(stol=0.25,dm=1500,fspot=50)
        """
        
        znewss=np.copy(self.z)
        
        #get directory name where all edi files should be
        edipath=os.path.dirname(self.edifn)
        
        #make a list of filename from edipath
        edilst=[os.path.join(edipath,dedifile) 
                for dedifile in os.listdir(edipath) if dedifile.find('.edi')>0]    
        
        rp=ResPhase(self.z,self.period)
        #loop over files to find nearby stations
        resxlst=[]
        resylst=[]
        statlst=[]
        zone,northing,easting=utm2ll.LLtoUTM(23,self.lat,self.lon)   
        for kk,kedi in enumerate(edilst):
            zk=Edi(kedi)
            zk.readEDI()
            zone,dn,de=utm2ll.LLtoUTM(23,zk.lat,zk.lon)        
            deltad=np.sqrt((dn-northing)**2+(de-easting)**2)
            if deltad<=dm:
                zkrp=ResPhase(zk.z,1./zk.frequency)
                resxlst.append(zkrp.resxy[0:fspot])
                resylst.append(zkrp.resyx[0:fspot])
                statlst.append(kk)
        
        #make arrays for easy manipulation
        resxlst=np.array(resxlst)
        resylst=np.array(resylst)
        
        #calculate static shift of x-components
        staticx=np.mean(rp.resxy[0:fspot]/np.median(resxlst))
        
        #see if it is within the tolerance level    
        if staticx<1-stol or staticx>1+stol:
            znewss[:,0,:]=self.z[:,0,:]/np.sqrt(staticx)
            print 'X-Shifted '+self.station+' by '+str(1/np.sqrt(staticx))
            xyn='y'
        else:
            xyn='n'
        
        #calculate static shift of y-components
        staticy=np.mean(rp.resyx[0:fspot]/np.median(resylst))
        
        #see if it is within the tolerance level
        if staticy<1-stol or staticy>1+stol:
            znewss[:,1,:]=self.z[:,1,:]/np.sqrt(staticy)
            print 'Y-Shifted '+self.station+' by '+str(1/np.sqrt(staticy)) 
            yyn='y'
        else:yyn='n'
        
        #if there was a correction write a new edi file
        if xyn=='y' or yyn=='y':
            self.rewriteedi(znew=znewss,zvarnew=self.zvar,ext='SS')
            
        #if no correction was made return the same edifile
        else:
            print 'No Static Shift Correction for ',self.station
 
    def getResPhase(self,ffactor=1,thetar=0):
        """
        getResPhase will return a ResPhase class that has attributes of all 4
        components of resistivity and phase as well as the errors and 
        determinants.
        
        Arguments:
        ----------
            **ffactor** : float
                          a factor if the calibration or gains aren't quite 
                          correct
                          
            **thetar** : float (angle in degrees)
                         rotation angle clockwise positive assuming 0 is North.
        
        Returns:
        --------
            **ResPhase** : data type ResPhase
            
        :Example: ::
            
            >>> z1 = Z.Z(edifile)
            >>> rp = z1.getResPhase()
        
        """
    
        return ResPhase(self.z,self.period,zvar=self.zvar,rotz=thetar,
                        ffactor=ffactor)

    def getResTensor(self,thetar=0,rotate=180):
        """
        getResTensor will return a data type that describes the reistivity 
        tensor defined by O'reilly [1978] and Weckmann et al. [2003]
        
        Arguments:
        ----------
            **thetar** : float (angle in degrees)
                         rotation angle clockwise positive assuming 0 is North.
                         *Default* is 0
                         
            **rotate** : int [ 90 | 180 | 270 ] 
                        rotation of coordinate system, default is Y North and 
                        X East.  So if X is North rotate is 180. 
                        *Default* is 180
                        
        Returns:
            **ResistivityTensor** : data type Resistivity Tensor
            
        :Example: ::
            
            >>> z1 = Z.Z(edifile)
            >>> rt = z1.getResTensor()
            
        """
        
        return ResistivityTensor(self.z,self.frequency,rotz=thetar,
                                 rotate=rotate)
    
                
    def plotResPhase(self,thetar=0,fignum=1,plottype=1,title=None,ffactor=1,
                     savefigfilename=None,dpi=None,fmt=None,orientation=None,
                     phaselimits=(0,90),reslimits=None):
        """
        Will plot the apparent resistivity and phase for TE and TM modes or all
        modes.  If there is tipper data it will be plotted at the bottom as 
        real and imaginary parts.
    
        Arguments:
        ----------
            **fignum** : int (figure number). *Default* is 1
            
            **thetar** : float (angle in degrees)
                         rotation angle clockwise positive assuming 0 is North.
                         *Default* is 0
                         
            **ffactor** : float
                          a factor if the calibration or gains aren't quite 
                          correct.  *Default* is 1
                          
            **plottype** : int [ 1 | 2 | 3 ]
                            * 1 for just Ex/By and Ey/Bx
                            * 2 for all 4 components
                            * 3 for off diagonal plus the determinant
                            * Default is 1
                            
            **title** : string 
                        title of plot
                        
            **savefigfilename** : string
                                 supply full path filename to save figure 
                                 *Default* is None 
                                 
            **dpi** : int
                      Dots-per-inch resolution of figure.
                      *Default* is 100
                      
            **format** : string [ 'pdf' | 'eps' | 'svg' | 'png' | 'jpeg' ]
                         if savefigfilename!=None, save figure to format
            
            **orientation** : string [ 'landscape' | 'portrait' ]
                              orientation of figure on A4 paper
                              
            **phaselimits** : tuple (min,max)
                              min and max of phase limits in degrees. 
                              *Default* is (0,90)
            
            **reslimits** : tuple (min,max)
                            min and max of resistivity limits in log10
                            *Default* is None to automatically detect limits
            
        :Example: ::
            
            # to plot all 4 components
            >>> z1 = Z.Z(edifile)
            >>> z1.plotResPhase(plottype=2)
            
        """
        
        if dpi==None:
            dpi=100
        
        rp=ResPhase(self.z,self.period,zvar=self.zvar,rotz=thetar,
                       ffactor=ffactor)
        tp=Tipper(self.tipper,rott=thetar)
        if tp.magreal[0]==0.0 and tp.magimag[0]==0.0:
            tpyn='n'
        else:
            tpyn='y'               
        
        plt.rcParams['font.size']=10
        plt.rcParams['figure.subplot.left']=.13
        plt.rcParams['figure.subplot.right']=.98
        plt.rcParams['figure.subplot.bottom']=.1
        plt.rcParams['figure.subplot.top']=.95
        plt.rcParams['figure.subplot.wspace']=.25
        plt.rcParams['figure.subplot.hspace']=.05
        #plt.rcParams['font.family']='helvetica'
        fontdict={'size':14,'weight':'bold'}
        
        if tpyn=='y':
            gs=gridspec.GridSpec(3,2,height_ratios=[2,1.5,1],hspace=.05)
        if tpyn=='n':
            gs=gridspec.GridSpec(2,2,height_ratios=[2,1.5],hspace=.05)
        
            
        #make figure for xy,yx components
        if plottype==1 or plottype==3: 
            fig=plt.figure(fignum,[8,10],dpi=dpi)
            plt.clf()
            gs.update(hspace=.05,wspace=.15,left=.1)
        elif plottype==2:
            fig=plt.figure(fignum,[10,10],dpi=dpi)
            plt.clf()
            gs.update(hspace=.05,wspace=.15,left=.07)
        
        
        #---------plot the apparent resistivity-----------------------------------
        if plottype==1  or plottype==3:
            if tpyn=='n':
                ax=fig.add_subplot(gs[0,:])
                ax2=fig.add_subplot(gs[1,:],sharex=ax)
                ax.yaxis.set_label_coords(-.055, 0.5)
                ax2.yaxis.set_label_coords(-.055, 0.5)
            if tpyn=='y':
                ax=fig.add_subplot(gs[0,:])
                ax2=fig.add_subplot(gs[1,:],sharex=ax)
                ax5=fig.add_subplot(gs[2,:],sharex=ax)
                ax.yaxis.set_label_coords(-.055, 0.5)
                ax2.yaxis.set_label_coords(-.055, 0.5)
                ax5.yaxis.set_label_coords(-.055, 0.5)
        elif plottype==2:
            if tpyn=='n':
                ax=fig.add_subplot(gs[0,0])
                ax2=fig.add_subplot(gs[1,0],sharex=ax)
                ax.yaxis.set_label_coords(-.075, 0.5)
                ax2.yaxis.set_label_coords(-.075, 0.5)
            if tpyn=='y':
                ax=fig.add_subplot(gs[0,0])
                ax2=fig.add_subplot(gs[1,0],sharex=ax)
                ax5=fig.add_subplot(gs[2,:],sharex=ax)
                ax.yaxis.set_label_coords(-.075, 0.5)
                ax2.yaxis.set_label_coords(-.075, 0.5)
                ax5.yaxis.set_label_coords(-.075, 0.5)
        
        #--------Plot Resistivity-----------------------------------------------
        erxy=ax.errorbar(self.period,rp.resxy,marker='s',ms=4,mfc='None',
                         mec='b',mew=1,ls='None',yerr=rp.resxyerr,ecolor='b')
        eryx=ax.errorbar(self.period,rp.resyx,marker='o',ms=4,mfc='None',
                         mec='r',mew=1,ls='None',yerr=rp.resyxerr,ecolor='r')
        #ax.set_xlabel('Period (s)',fontdict=fontdict)
        pylab.setp( ax.get_xticklabels(), visible=False)
        ax.set_ylabel('app. res. ($\mathbf{\Omega \cdot m}$)',
                   fontdict=fontdict)
        ax.set_yscale('log')
        ax.set_xscale('log')
        ax.set_xlim(xmin=10**(np.floor(np.log10(self.period[0]))),
                    xmax=10**(np.ceil(np.log10(self.period[-1]))))
        if reslimits!=None:
            ax.set_ylim(ymin=10**reslimits[0],ymax=10**reslimits[1])
        ax.grid(True)
        ax.legend((erxy[0],eryx[0]),('$E_x/B_y$','$E_y/B_x$'),loc=3,
                    markerscale=1,borderaxespad=.01,labelspacing=.07,
                    handletextpad=.2,borderpad=.02)
        
        #-----Plot the phase----------------------------------------------------
        
        ax2.errorbar(self.period,rp.phasexy,marker='s',ms=4,mfc='None',mec='b',
                     mew=1,ls='None',yerr=rp.phasexyerr,ecolor='b')
        ax2.errorbar(self.period,np.array(rp.phaseyx)+180,marker='o',ms=4,
                     mfc='None',mec='r',mew=1,ls='None',yerr=rp.phaseyxerr,
                     ecolor='r')
        if tpyn=='y':
            pylab.setp(ax2.get_xticklabels(),visible=False)
        elif tpyn=='n':
            ax2.set_xlabel('period (s)',fontdict)
        ax2.set_ylabel('phase (deg)',fontdict)
        ax2.set_xscale('log')

        #check the phase to see if any point are outside of [0:90]    
        ax2.set_ylim(phaselimits)        
        ax2.yaxis.set_major_locator(MultipleLocator(30))
        ax2.yaxis.set_minor_locator(MultipleLocator(1))
        ax2.grid(True)
        
        #Add xx and yy components to the plot
        if plottype==2:
            #---------plot the apparent resistivity-----------------------------------
            ax3=fig.add_subplot(gs[0,1],sharex=ax)
            ax3.yaxis.set_label_coords(-.1, 0.5)
            erxx=ax3.errorbar(self.period,rp.resxx,marker='s',ms=4,mfc='None',
                              mec='b',mew=1,ls='None',yerr=rp.resxxerr,
                              ecolor='b')
            eryy=ax3.errorbar(self.period,rp.resyy,marker='o',ms=4,mfc='None',
                              mec='r',mew=1,ls='None',yerr=rp.resyyerr,
                              ecolor='r')
            ax3.set_yscale('log')
            ax3.set_xscale('log')
            pylab.setp( ax3.get_xticklabels(), visible=False)
            
            
            ax3.set_xlim(xmin=10**(np.floor(np.log10(self.period[0]))),
                     xmax=10**(np.ceil(np.log10(self.period[-1]))))
            ax3.grid(True)
            ax3.legend((erxx[0],eryy[0]),('$E_x/B_x$','$E_y/B_y$'),loc=3,
                        markerscale=1,borderaxespad=.01,labelspacing=.07,
                        handletextpad=.2,borderpad=.02)
            
            #-----Plot the phase----------------------------------------------------
            ax4=fig.add_subplot(gs[1,1],sharex=ax3)
            
            ax4.yaxis.set_label_coords(-.1, 0.5)
            ax4.errorbar(self.period,rp.phasexx,marker='s',ms=4,mfc='None',
                         mec='b',mew=1,ls='None',yerr=rp.phasexxerr,ecolor='b')
            ax4.errorbar(self.period,np.array(rp.phaseyy),marker='o',ms=4,
                         mfc='None',mec='r',mew=1,ls='None',yerr=rp.phaseyyerr,
                         ecolor='r')
            if tpyn=='y':
                pylab.setp(ax4.get_xticklabels(), visible=False)
            else:
                ax4.set_xlabel('period (s)',fontdict)
            ax4.set_xscale('log')
            ax4.set_ylim(ymin=-180,ymax=180)        
            ax4.yaxis.set_major_locator(MultipleLocator(30))
            ax4.yaxis.set_minor_locator(MultipleLocator(5))
            ax4.grid(True)
        
        #Add determinant to the plot
        if plottype==3:
                
            #-------Plot resistivity--------------------------------------------
            erdet=ax.errorbar(self.period,rp.resdet,marker='d',ms=4,mfc='None',
                              mec='g',mew=1,ls='None',yerr=rp.resdeterr,
                              ecolor='g')
        
            ax.legend((erxy[0],eryx[0],erdet[0]),('$E_x/B_y$','$E_y/B_x$',
                      '$\det(\mathbf{\hat{Z}})$'),loc=3,markerscale=1,
                        borderaxespad=.01,labelspacing=.07,handletextpad=.2,
                        borderpad=.02)
            
            #-----Plot the phase------------------------------------------------
            
            ax2.errorbar(self.period,rp.phasedet,marker='d',ms=4,mfc='None',
                         mec='g',mew=1,ls='None',yerr=rp.phasedeterr,ecolor='g')
        
        if tpyn=='y':
            txr=tp.magreal*np.cos(tp.anglereal*np.pi/180)
            tyr=tp.magreal*np.sin(tp.anglereal*np.pi/180)
    
            txi=tp.magimag*np.cos(tp.angleimag*np.pi/180)
            tyi=tp.magimag*np.sin(tp.angleimag*np.pi/180)
            
            nt=len(txr)
            
            for aa in range(nt):
                if np.log10(self.period[aa])<0:
                     hwidth=.1*self.period[aa]
                     hheight=.01*self.period[aa]
                else:
                    hwidth=.2/self.period[aa]
                    hheight=.01/self.period[aa]
                overhang=1
                ax5.arrow(self.period[aa],0,txr[aa],tyr[aa],lw=2,facecolor='k',
                          edgecolor='k',head_width=hwidth,head_length=hheight,
                          length_includes_head=False,overhang=overhang)
                ax5.arrow(self.period[aa],0,txi[aa],tyi[aa],lw=2,facecolor='b',
                          edgecolor='b',length_includes_head=False)
                          
            ax5.plot(self.period,[0]*nt,'k')
            ax5.set_xscale('log')              
            line1=ax5.plot(0,0,'k')
            line2=ax5.plot(0,0,'b')
#            ax5.set_xlim(xmin=10**np.floor(np.log10(period[0])),
#                         xmax=10**np.ceil(np.log10(period[-1])))
            ax5.yaxis.set_major_locator(MultipleLocator(.2))
                         
            ax5.legend([line1[0],line2[0]],['real','imag'],loc='upper left',markerscale=1,
                       borderaxespad=.01,labelspacing=.07,handletextpad=.2,borderpad=.02)
            ax5.set_xlabel('period (s)',fontdict={'size':12,'weight':'bold'})
            ax5.set_ylabel('tipper',fontdict={'size':12,'weight':'bold'})    
            
            ax5.set_ylim([-1,1])
            ax5.grid(True)
        
        #make title and show
        if title!=None:
            plt.suptitle(title,fontsize=14,fontweight='bold')
        else:
            plt.suptitle(self.station,fontsize=14,fontweight='bold')
        plt.show()
        if savefigfilename!=None:
                if dpi==None:
                    dpi=100
                if format==None:
                    format='pdf'
                if orientation==None:
                    orientation='landscape'
                fig.savefig(savefigfilename,dpi=dpi,format=format,
                            orientation=orientation)
                plt.clf()
                plt.close(fig)
                
    def plotPTAll(self,xspacing=6,esize=5,fignum=1,thetar=0, save='n',
                  savepath=None,fmt='pdf',coordrot=180,ptmm=None,rpmm=None,
                  dpi=300,restensor='n'):
        """
        Will plot phase tensor, strike angle, min and max phase angle, 
        azimuth, skew, and ellipticity as subplots on one plot.  It can plot
        the resistivity tensor along side the phase tensor for comparison.
        
        Arguments:
        ----------
            **xspacing** : float 
                            spacing of tensors along x direction.
                            *Default* is 6
                            
            **esize** : float
                        size of tensor ellipses.
                        *Default* is 5
                        
            **fignum** : int (figure number)
            
            **thetar** : float (angle in degrees)
                         rotation angle clockwise positive assuming 0 is North.
                         *Default* is 0
                         
            **save** : [ 'y' | 'n' ]
                       * 'y' to save figure
                       * 'n' to not save figure
                       * *Default* is 'n'
            
            **savepath** : string
                           path to save to, saves as savepath\statioAll.fmt
                           
            **fmt** : [ 'pdf' | 'eps' | 'svg' | 'png' | 'jpeg' ]
                      format of save figure pdf,svg,eps,ps,png
            
            **dpi** : int
                      Dots-per-inch resolution of figure.
                      *Default* is 300
            
            **coordrot** : [ 90 | 180 | 270 ]
                          rotation of coordinate directions assuming Y North
                          and X East.  If data measured in X North, Y East
                          coordrot=180.  *Default* is 180
                          
            **rpmm** : tuple (min,max)
                       min and max of resistivity tensor on log10 scale
            
            **ptmm** : tuple (min,max)
                       min an max of phase tensor
            
            **restensor** : [ 'y' | 'n' ]
                            * 'y' to plot the resistivity tensors along with 
                               the phase tensors
                            * 'n' to not plot resistivity tensors
                            * *Default* is 'n'
                            
        :Example: ::
            
            #To plot just the phase tensor components
            >>> z1 = Z.Z(edifile)
            >>> z1.plotPTAll()
            
        """
        
        #Set plot parameters
        plt.rcParams['font.size']=8
        plt.rcParams['figure.subplot.left']=.1
        plt.rcParams['figure.subplot.right']=.98
        plt.rcParams['figure.subplot.bottom']=.1
        plt.rcParams['figure.subplot.top']=.95
        plt.rcParams['figure.subplot.wspace']=.21
        plt.rcParams['figure.subplot.hspace']=.5
        
        fs=8
        tfs=10
        
        #begin plotting
        fig=plt.figure(fignum,[8,10],dpi=dpi)
        plt.clf()
        
        pt=PhaseTensor(self.z,self.zvar,rotate=coordrot,rotz=thetar)
        rp=ResistivityTensor(self.z,self.frequency,rotate=coordrot,rotz=thetar)
        zinv=Zinvariants(self.z,rotz=thetar)
        #tipper=Tipper(rott=thetar)
        period=self.period
        n=len(period)
        
        pperiod=np.logspace(np.floor(np.log10(period[0])),
                            np.ceil(np.log10(period[-1])),10*n)

        if rpmm==None:
            rpmin=min(np.log10(rp.rhodet))
            rpmax=max(np.log10(rp.rhodet))
        else:
            rpmin=float(rpmm[0])
            rpmax=float(rpmm[1])
        
        if ptmm==None:
            ptmin=min(pt.phimin)
            ptmax=max(pt.phimin)
        else:
            ptmin=float(ptmm[0])*np.pi/180
            ptmax=float(ptmm[1])*np.pi/180
            
        #-------------plotPhaseTensor-----------------------------------
        ax1=fig.add_subplot(3,1,1,aspect='equal')
        for ii in range(n):
            #make sure the ellipses will be visable
            scaling=esize/pt.phimax[ii]
            eheight=pt.phimin[ii]*scaling
            ewidth=pt.phimax[ii]*scaling
            
            rheight=rp.rhomin[ii]/max([rp.rhomin[ii],rp.rhomax[ii]])*esize
            rwidth=rp.rhomax[ii]/max([rp.rhomin[ii],rp.rhomax[ii]])*esize
                
            #create an ellipse scaled by phimin and phimax and oriented along
            #the azimuth    
            ellip=Ellipse((xspacing*ii,0),width=ewidth,
                          height=eheight,
                          angle=pt.azimuth[ii])
                          
            if pt.phimin[ii]<0 or pt.phimin[ii]=='nan':
                cvars=0
            else:
                cvars=(pt.phimin[ii]-ptmin)/(ptmax-ptmin)
                if cvars>1.0:
                    cvars=1
                
            ellip.set_facecolor((1,1-cvars,.1))                          
                          
            ax1.add_artist(ellip)
            
            #resistivity tensor
            if restensor=='y':
                rellip=Ellipse((xspacing*ii,3+esize),width=rwidth,
                              height=rheight,
                              angle=rp.rhoazimuth[ii])
                             
                cvar=(np.log10(rp.rhodet[ii])-rpmin)/(rpmax-rpmin)
                if cvar>.5:
                    if cvar>1:
                        rellip.set_facecolor((0,0,1))
                    else:
                        rellip.set_facecolor((1-abs(cvar),1-abs(cvar),1))
                else:
                    if cvar<-1:
                        rellip.set_facecolor((1,0,0))
                    else:
                        rellip.set_facecolor((1,1-abs(cvar),1-abs(cvar)))
                
                ax1.add_artist(rellip)
            else:
                pass
        
        #set tick labels and limits
        xticklabels=[]
        for xx in np.arange(start=0,stop=n,step=5):
            if period[xx]<100:
                xticklabels.append('{0:.2g}'.format(period[xx]))
            elif period[xx]>100 and period[xx]<1000:
                xticklabels.append('{0:.3g}'.format(period[xx]))
            elif period[xx]>1000 and period[xx]<10000:
                xticklabels.append('{0:.4g}'.format(period[xx]))
        plt.xlabel('Period (s)',fontsize=8,fontweight='bold')
        #plt.title('Phase Tensor Ellipses for '+stationstr,fontsize=14)
        plt.xticks(np.arange(start=0,stop=xspacing*n,step=5*xspacing),
                   xticklabels)
        if restensor=='y':
            ax1.set_ylim(-esize,2*esize+3)
        else:
            ax1.set_ylim(-esize,esize)
        ax1.set_xlim(-xspacing,n*xspacing+3)
        ax1.grid(alpha=.3)
        plt.setp(ax1.get_yticklabels(),visible=False)
        
        #add colorbar for PT
        cbpt=fig.add_subplot(3,32,1)
        cbpt.set_axis_off()
        if restensor=='y':
            ax1cbpt=make_axes(ax1,shrink=.7,orientation='vertical',pad=.005)
        else:
            ax1cbpt=make_axes(ax1,shrink=.7,orientation='vertical',pad=.01)
        cb1=ColorbarBase(ax1cbpt[0],cmap=ptcmap,
                        norm=Normalize(vmin=ptmin*180/np.pi,
                                       vmax=ptmax*180/np.pi),
                        orientation='vertical')
        cb1.set_ticks([ptmin*180/np.pi,ptmax*180/np.pi])
        cb1.set_ticklabels(['{0:.0f}'.format(ptmin*180/np.pi),
                            '{0:.0f}'.format(ptmax*180/np.pi)])
        
        #add color bar for RT
        if restensor=='y':
            cbrt=fig.add_subplot(3,32,1)
            cbrt.set_axis_off()
            ax1cbrt=make_axes(ax1,shrink=.7,orientation='vertical',pad=.01)
            cb2=ColorbarBase(ax1cbrt[0],cmap=rtcmap,
                            norm=Normalize(vmin=10**rpmin,
                                           vmax=10**rpmax),
                            orientation='vertical')
            cb2.draw_all()
            cb2.set_ticks([10**rpmin,10**rpmax])
            cb2.set_ticklabels(['{0:.2g}'.format(10**rpmin),
                                '{0:.5g}'.format(10**rpmax)]) 
        
        #---------------plotStrikeAngle-----------------------------------
        
        az=90-np.array(pt.azimuth)
        azvar=np.array(pt.azimuthvar)
        #put the strike into a coordinate system that goes from -90 to 90
        az[np.where(az>90)]=az[np.where(az>90)]-180
        az[np.where(az<-90)]=az[np.where(az<-90)]+180
        
        strike=zinv.strike
        #put the strike into a coordinate system that goes from -90 to 90
        strike[np.where(strike>90)]=strike[np.where(strike>90)]-180
        strike[np.where(strike<-90)]=strike[np.where(strike<-90)]+180
        
        ax2=plt.subplot(3,2,3)
        
        #plot invariant strike
        erxy=ax2.errorbar(period,strike,
                          marker='s',ms=4,mfc='None',
                          mec='c',mew=1,ls='None',
                          yerr=zinv.strikeerr,
                          ecolor='c')
        #plot phase tensor strike
        eraz=ax2.errorbar(period,az,marker='o',ms=4,
                          mfc='None',mec='purple',mew=1,
                          ls='None',yerr=azvar,ecolor='purple')
                          
        ax2.legend((erxy[0],eraz[0]),('Strike','Azimuth'),loc='lower left',
                   markerscale=.2,borderaxespad=.01,labelspacing=.1,
                   handletextpad=.2,ncol=2,borderpad=.1,columnspacing=.1)
        leg = plt.gca().get_legend()
        ltext  = leg.get_texts()  # all the text.Text instance in the legend
        plt.setp(ltext, fontsize=6)    # the legend text fontsize

        
        ax2.set_yscale('linear')
        ax2.set_xscale('log')
        ax2.set_xlim(xmax=10**(np.ceil(np.log10(period[-1]))),
                 xmin=10**(np.floor(np.log10(period[0]))))
        ax2.set_ylim(ymin=-95,ymax=95)
        ax2.yaxis.set_major_locator(MultipleLocator(20))
        ax2.yaxis.set_minor_locator(MultipleLocator(5))
        ax2.grid(True,alpha=.3)
        #plt.xlabel('Period (s)',fontsize=fs,fontweight='bold')
        ax2.set_ylabel('Angle (deg)',fontsize=fs,fontweight='bold')
        ax2.set_title('Strike Angle, Azimuth',fontsize=tfs,
                      fontweight='bold')
        
        #---------plot Min & Max Phase-----------------------------------------
        
        minphi=pt.phiminang
        minphivar=pt.phiminangvar
        maxphi=pt.phimaxang
        maxphivar=pt.phimaxangvar

        ax3=plt.subplot(3,2,4,sharex=ax2)
        ermin=plt.errorbar(period,minphi,marker='o',ms=4,
                           mfc='None',mec='r',mew=1,ls='None',
                           yerr=minphivar,ecolor='r')
        ermax=plt.errorbar(period,maxphi,marker='s',ms=4,
                           mfc='None',mec='b',mew=1,
                           ls='None',yerr=maxphivar,
                           ecolor='b')
        ax3.set_xscale('log')
        ax3.set_yscale('linear')
        plt.legend((ermin[0],ermax[0]),('$\phi_{min}$','$\phi_{max}$'),
                   loc='lower left',markerscale=.2,borderaxespad=.01,
                   labelspacing=.1,handletextpad=.2,ncol=2,borderpad=.01,
                   columnspacing=.01)
        leg = plt.gca().get_legend()
        ltext  = leg.get_texts()  # all the text.Text instance in the legend
        plt.setp(ltext, fontsize=6.5)    # the legend text fontsize
        plt.xlim(xmax=10**(np.ceil(np.log10(period[-1]))),
                 xmin=10**(np.floor(np.log10(period[0]))))
        plt.ylim(ymin=0,ymax=90)
        plt.grid(True,alpha=.3)
        #plt.xlabel('Period (s)',fontsize=fs,fontweight='bold')
        plt.ylabel('Phase (deg)',fontsize=fs,fontweight='bold')
        plt.title('$\mathbf{\phi_{min}}$ and $\mathbf{\phi_{max}}$',
                  fontsize=tfs,fontweight='bold')

        
        #-----------------------plotSkew---------------------------------------
        
        skew=pt.beta
        skewvar=pt.betavar

        ax5=plt.subplot(3,2,5,sharex=ax2)
        erskew=plt.errorbar(period,skew,marker='s',ms=4,
                            mfc='None',mec='g',mew=1,
                            ls='None',yerr=skewvar,
                            ecolor='g')
        ax5.plot(pperiod,[2]*(n*10),'--k',lw=1)
        ax5.plot(pperiod,[-2]*(n*10),'--k',lw=1)
        ax5.set_xscale('log')
        ax5.set_yscale('linear')
        ax5.yaxis.set_major_locator(MultipleLocator(2))
        plt.xlim(xmax=10**(np.ceil(np.log10(period[-1]))),xmin=10**(
                            np.floor(np.log10(period[0]))))
        plt.ylim(ymin=skew.min()-2,ymax=skew.max()+2)
        plt.grid(True,alpha=.3)
        plt.xlabel('Period (s)',fontsize=fs,fontweight='bold')
        plt.ylabel('Skew Angle (deg)',fontsize=fs,fontweight='bold')
        plt.title('Skew Angle',fontsize=tfs,fontweight='bold')
        
        #----------------------plotEllipticity--------------------------------
        
        ellipticity=pt.ellipticity
        ellipticityvar=pt.ellipticityvar

        ax6=plt.subplot(3,2,6,sharex=ax2)
        erskew=plt.errorbar(period,ellipticity,marker='s',
                            ms=4,mfc='None',mec='orange',mew=1,
                            ls='None',yerr=ellipticityvar,
                            ecolor='orange')
        ax6.set_xscale('log')
        ax6.set_yscale('linear')
        ax6.plot(pperiod,[.2]*(n*10),'--k',lw=1)
        ax6.yaxis.set_major_locator(MultipleLocator(.1))
        plt.xlim(xmax=10**(np.ceil(np.log10(period[-1]))),
                 xmin=10**(np.floor(np.log10(period[0]))))
        plt.ylim(ymin=0,ymax=1)
        #plt.yticks(range(10),np.arange(start=0,stop=1,step=.1))
        plt.grid(True,alpha=.3)
        plt.xlabel('Period (s)',fontsize=fs,fontweight='bold')
        plt.ylabel('$\mathbf{\phi_{max}-\phi_{min}/\phi_{max}+\phi_{min}}$',
                   fontsize=fs,fontweight='bold')
        plt.title('Ellipticity',fontsize=tfs,fontweight='bold')
        #plt.suptitle(self.z.station,fontsize=tfs,fontweight='bold')
        plt.suptitle('Phase Tensor Elements for: '+self.station,fontsize=12,
                     fontweight='bold')
            
        if save=='y':
            if savepath.find('.')==-1:
                if not os.path.exists(savepath):
                    os.mkdir(self.savepath)
                print 'Made Directory: '+ savepath
                fig.savefig(os.path.join(savepath,self.station+'All.'+fmt),fmt=fmt)
                print 'Saved figure to: '+os.path.join(self.savepath,
                                               self.z[0].station+'All.'+fmt)
                plt.close()
            else:
                fig.savefig(savepath,fmt=fmt)
        elif save=='n':
            pass
        
    def plotTipper(self,thetar=0,fignum=1,plotnum=1,dpi=100):
        """
        plotTipper will plot the resistivity, phase and tipper
        
        Arguments:
        ----------
            **fignum** : int (figure number). *Default* is 1
            
            **thetar** : float (angle in degrees)
                         rotation angle clockwise positive assuming 0 is North.
                         *Default* is 0
        
            **dpi** : int
                      Dots-per-inch resolution of figure.
                      *Default* is 100
                      
            **plotnum** : [ 1 | 2 | 3 ]
                          * 1 for just Ex/By and Ey/Bx
                          * 2 for all 4 components
                          * 3 for off diagonal plus the determinant
                          * Default is 1
        :Example: ::
            
            #To plot all 4 components and the tipper
            >>> z1 = Z.Z(edifile)
            >>> z1.plotTipper(plotnum=2)
        """

        rp=ResPhase(self.z,self.period,zvar=self.zvar,rotz=thetar,
                       ffactor=1)
        tip=Tipper(self.tipper,rott=thetar)
        
        period=self.period
        nt=len(tip.tipper)

        txr=tip.magreal*np.cos(tip.anglereal*np.pi/180)
        tyr=tip.magreal*np.sin(tip.anglereal*np.pi/180)

        txi=tip.magimag*np.cos(tip.angleimag*np.pi/180)
        tyi=tip.magimag*np.sin(tip.angleimag*np.pi/180)
        
        plt.rcParams['font.size']=10
        plt.rcParams['figure.subplot.left']=.13
        plt.rcParams['figure.subplot.right']=.98
        plt.rcParams['figure.subplot.bottom']=.1
        plt.rcParams['figure.subplot.top']=.95
        plt.rcParams['figure.subplot.wspace']=.25
        plt.rcParams['figure.subplot.hspace']=.05
        #plt.rcParams['font.family']='helvetica'
        fontdict={'size':14,'weight':'bold'}
        
        gs=gridspec.GridSpec(3,2,height_ratios=[2,1.5,1],hspace=.05)
        
            
        #make figure for xy,yx components
        if plotnum==1 or plotnum==3: 
            fig=plt.figure(fignum,[8,10],dpi=dpi)
            gs.update(hspace=.05,wspace=.15,left=.1)
        elif plotnum==2:
            fig=plt.figure(fignum,[10,10],dpi=dpi)
            gs.update(hspace=.05,wspace=.15,left=.07)
        
        #---------plot the apparent resistivity-----------------------------------
        if plotnum==1  or plotnum==3:
            ax=plt.subplot(gs[0,:])
            ax2=plt.subplot(gs[1,:])
            ax3=plt.subplot(gs[2,:])
            ax.yaxis.set_label_coords(-.055, 0.5)
            ax2.yaxis.set_label_coords(-.055, 0.5)
            ax3.yaxis.set_label_coords(-.055, 0.5)
        elif plotnum==2:
            ax=plt.subplot(gs[0,0])
            ax2=plt.subplot(gs[1,0])
            ax.yaxis.set_label_coords(-.075, 0.5)
            ax2.yaxis.set_label_coords(-.075, 0.5)
        
        
        erxy=ax.errorbar(period,rp.resxy,marker='s',ms=4,mfc='None',mec='b',
                          mew=1,ls='None',yerr=rp.resxyerr,ecolor='b')
        eryx=ax.errorbar(period,rp.resyx,marker='o',ms=4,mfc='None',mec='r',
                          mew=1,ls='None',yerr=rp.resyxerr,ecolor='r')
        #ax.set_xlabel('Period (s)',fontdict=fontdict)
        pylab.setp( ax.get_xticklabels(), visible=False)
        ax.set_ylabel('App. Res. ($\mathbf{\Omega \cdot m}$)',
                   fontdict=fontdict)
        ax.set_yscale('log')
        ax.set_xscale('log')
        ax.set_xlim(xmin=10**(np.floor(np.log10(period[0]))),
                 xmax=10**(np.ceil(np.log10(period[-1]))))
        ax.grid(True)
        ax.legend((erxy[0],eryx[0]),('$E_x/B_y$','$E_y/B_x$'),loc=3,
                    markerscale=1,borderaxespad=.01,labelspacing=.07,
                    handletextpad=.2,borderpad=.02)
        
        #-----Plot the phase----------------------------------------------------
        
        ax2.errorbar(period,rp.phasexy,marker='s',ms=4,mfc='None',mec='b',
                     mew=1,ls='None',yerr=rp.phasexyerr,ecolor='b')
        ax2.errorbar(period,np.array(rp.phaseyx)+180,marker='o',ms=4,
                     mfc='None',mec='r',mew=1,ls='None',yerr=rp.phaseyxerr,
                     ecolor='r')
        ax2.set_xlabel('Period (s)',fontdict)
        ax2.set_ylabel('Phase (deg)',fontdict)
        ax2.set_xscale('log')
        #ax2.set_xlim(xmin=10**(np.floor(np.log10(period[0]))),
        #         xmax=10**(np.ceil(np.log10(period[-1]))))
        #check the phase to see if any point are outside of [0:90]    
        if min(rp.phasexy)<0 or min(rp.phaseyx+180)<0:
            pymin=min([min(rp.phasexy),min(rp.phaseyx+180)])
            if pymin>0:
                pymin=0
        else:
            pymin=0
        
        if max(rp.phasexy)>90 or max(rp.phaseyx+180)>90:
            pymax=min([max(rp.phasexy),max(rp.phaseyx+180)])
            if pymax<91:
                pymax=90
        else:
            pymax=90
        
        pylab.setp( ax2.get_xticklabels(), visible=False)
        ax2.set_ylim(ymin=pymin,ymax=pymax)        
        ax2.yaxis.set_major_locator(MultipleLocator(30))
        ax2.yaxis.set_minor_locator(MultipleLocator(1))
        ax2.grid(True)
        
        #------Plot Tipper-----------------------------------------------------
        
        
        ml=2*(np.ceil(np.log10(period[-1]))-np.floor(np.log10(period[0])))/nt
        
        for aa in range(nt):
            if np.log10(period[aa])<0:
                 hwidth=.1*period[aa]
                 hheight=.01*period[aa]
            else:
                hwidth=.2/period[aa]
                hheight=.01/period[aa]
            overhang=1
            ax3.arrow(period[aa],0,txr[aa],tyr[aa],lw=2,facecolor='k',
                      edgecolor='k',head_width=hwidth,head_length=hheight,
                      length_includes_head=False,overhang=overhang)
            ax3.arrow(period[aa],0,txi[aa],tyi[aa],lw=2,facecolor='b',
                      edgecolor='b',length_includes_head=False)
        
        
        ax3.plot(period,[0]*nt,'k')
        ax3.set_xscale('log')              
        line1=ax3.plot(0,0,'k')
        line2=ax3.plot(0,0,'b')
        ax3.set_xlim(xmin=10**np.floor(np.log10(period[0])),
                     xmax=10**np.ceil(np.log10(period[-1])))
        ax3.yaxis.set_major_locator(MultipleLocator(.1))
                     
        ax3.legend([line1[0],line2[0]],['Real','Imag'],loc='upper left',markerscale=1,
                   borderaxespad=.01,labelspacing=.07,handletextpad=.2,borderpad=.02)
        ax3.set_xlabel('Period (s)',fontdict={'size':12,'weight':'bold'})
        ax3.set_ylabel('Tipper',fontdict={'size':12,'weight':'bold'})    
        
        ax3.set_ylim([tyi.min()-.1,tyr.max()+.1])
        ax3.grid(True)
        
        plt.show()        
            

                    
class PhaseTensor:
    """
    PhaseTensor calculates the components of the phase tensor following 
    Caldwell et al. [2004].
    
    Arguments:
    ----------
        **z** : complex np.array (nf,2,2)
                impedance tensor
        
        **zvar** : np.array (nf,2,2)
                   impedance tensor variances
                   
        **rotate** : [ 90 | 180 | 270 ] 
                     Caldwell et al, [2004] assume the coordinate axis is 
                     Y North and X East.  If the data is in X North and Y 
                     East than rotation = 180. *Default* is 180
        
        **rotz** : float (angle in degrees)
                   rotation angle clockwise positive assuming 0 is North.
        
    """

    def __init__(self,z,zvar=None,rotate=180,rotz=0):
        #make z an array
        try:
            nz=len(z)
        except TypeError:
            z=np.array([z])
            nz=len(z)
                
        #make a copy of z    
        self.z=z.copy()
        if zvar==None:
            self.zvar=np.zeros_like(z)
        else:
            try:
                self.zvar=zvar.copy()
            except AttributeError:
                zvar=np.array([zvar])                
                self.zvar=zvar.copy()
                
        self.rot=rotate
        self.thetar=rotz
        self.phi=np.zeros((nz,2,2))
        self.phivar=np.zeros_like(self.phi)
        self.phimin=np.zeros(nz)
        self.phiminvar=np.zeros(nz)
        self.phimax=np.zeros(nz)
        self.phimaxvar=np.zeros(nz)
        self.alpha=np.zeros(nz)
        self.alphavar=np.zeros(nz)
        self.beta=np.zeros(nz)
        self.betavar=np.zeros(nz)
        self.azimuth=np.zeros(nz)
        self.azimuthvar=np.zeros(nz)
        self.phiminang=np.zeros(nz)
        self.phiminangvar=np.zeros(nz)
        self.phimaxang=np.zeros(nz)
        self.phimaxangvar=np.zeros(nz)
        self.ellipticity=np.zeros(nz)
        self.ellipticityvar=np.zeros(nz)
        self.phidet=np.zeros(nz)
        self.phidetvar=np.zeros(nz)
        
        #rotate data if desired
        if self.thetar!=0:
            #convert thetar into radians
            if abs(self.thetar)>2*np.pi:
                self.thetar=self.thetar*(np.pi/180)
            #make rotation matrix
            rotmatrix=np.array([[np.cos(self.thetar), np.sin(self.thetar)],
                         [-np.sin(self.thetar), np.cos(self.thetar)]])
            #rotate the data
            for rr in range(nz):
                self.z[rr]=np.dot(rotmatrix,np.dot(self.z[rr],rotmatrix.T))
                self.zvar[rr]=np.dot(rotmatrix,np.dot(self.zvar[rr],
                                        rotmatrix.T))
        
            
        for ii in range(nz):
            X=self.z[ii].real
            Y=self.z[ii].imag
            #calculate phase tensor and rotate  to align x pointing 
            #east and y pointing north following convention of Caldwell (2004)
            if rotate==90:
                try:
                    phi=np.rot90(np.dot(np.linalg.inv(X),Y),1)
                except np.linalg.LinAlgError:
                    phi=np.zeros((2,2))
            elif rotate==180:
                try:
                    phi=np.rot90(np.dot(np.linalg.inv(X),Y),2)
                except np.linalg.LinAlgError:
                    phi=np.zeros((2,2))
            elif rotate==270:
                try:
                    phi=np.rot90(np.dot(np.linalg.inv(X),Y),3)
                except np.linalg.LinAlgError:
                    phi=np.zeros((2,2))
            else:
                phi=np.dot(np.linalg.inv(X),Y)
            #put variance into standard deviation
            zvar=self.zvar[ii]**2
            #create a matrix for errors to be calculated
            dphi=np.zeros(np.shape(z[ii]))
            #compute determinate of X
            detX=np.linalg.det(X)
            #calculate the deteriminate of the error matrix 
            ddet=np.sqrt((X[0,0]*X[1,1])**2*((zvar[0,0]/X[0,0])**2+
                (zvar[1,1]/X[1,1])**2)+(X[1,0]*X[0,1])**2*(
                (zvar[0,1]/X[0,1])**2+(zvar[1,0]/X[1,0])**2))
            #calculate errors for each component of the matrix ala Caldwell 
            #2004
            dphi[0,0]=np.sqrt((X[1,1]*Y[0,0])**2*((zvar[1,1]/X[1,1])**2+
                (zvar[0,0]/Y[0,0])**2)+(X[0,1]*Y[1,0])**2*(
                (zvar[0,1]/X[0,1])**2+(zvar[1,0]/X[1,0])**2))
            dphi[0,1]=np.sqrt((X[1,1]*Y[0,1])**2*((zvar[1,1]/X[1,1])**2+
                (zvar[0,1]/Y[0,1])**2)+(X[0,1]*Y[1,1])**2*(
                (zvar[0,1]/X[0,1])**2+(zvar[1,1]/X[1,1])**2))
            dphi[1,0]=np.sqrt((X[0,1]*Y[1,0])**2*((zvar[0,0]/X[0,0])**2+
                (zvar[1,0]/Y[1,0])**2)+(X[1,0]*Y[0,0])**2*(
                (zvar[1,0]/X[1,0])**2+(zvar[0,0]/X[0,0])**2))
            dphi[1,1]=np.sqrt((X[0,0]*Y[1,1])**2*((zvar[0,0]/X[0,0])**2+
                (zvar[1,1]/Y[1,1])**2)+(X[1,0]*Y[0,1])**2*(
                (zvar[1,0]/X[1,0])**2+(zvar[0,1]/X[0,1])**2))
            #rotate the error matrix
            dphi=np.rot90(dphi,2)
            #finish calculating the errors
            dphi[0,0]=(phi[0,0]/detX)**2*np.sqrt((dphi[0,0]/phi[0,0])**2+
                                                            (ddet/detX)**2)
            dphi[0,1]=(phi[0,1]/detX)**2*np.sqrt((dphi[0,1]/phi[0,1])**2+
                                                            (ddet/detX)**2)
            dphi[1,0]=(phi[1,0]/detX)**2*np.sqrt((dphi[1,0]/phi[1,0])**2+
                                                            (ddet/detX)**2)
            dphi[1,1]=(phi[1,1]/detX)**2*np.sqrt((dphi[1,1]/phi[1,1])**2+
                                                            (ddet/detX)**2)
            
            #Calculate Trace of Phi and error of trace of phi
            tr=phi[0,0]+phi[1,1]
            trvar=np.sqrt(dphi[0,0]**2+dphi[1,1]**2)
            
            #Calculate skew of phi and the cooresponding error
            skew=phi[0,1]-phi[1,0]
            skewvar=np.sqrt(dphi[0,1]**2+dphi[1,1]**2)
            
            #calculate the determinate and determinate error of phi
            phidet=abs(np.linalg.det(phi))
            phidetvar=np.sqrt((np.sqrt((dphi[0,0]/phi[0,0])**2+(
                dphi[1,1]/phi[1,1])**2)*phi[0,0]*phi[1,1])**2+(
                np.sqrt((dphi[0,1]/phi[0,1])**2+(
                dphi[1,0]/phi[1,0])**2)*phi[0,1]*phi[1,0])**2)
            
            #calculate reverse trace and error
            revtr=phi[0,0]-phi[1,1]
            revtrvar=np.sqrt(dphi[0,0]**2+dphi[1,1]**2)
            
            #calculate reverse skew and error
            revskew=phi[1,0]+phi[0,1]
            revskewvar=np.sqrt(phi[0,1]**2+dphi[1,0]**2)
            
            #calculate skew angle beta and error
            beta=.5*(np.arctan2(skew,tr)*180/np.pi)
            betavar=abs(np.arctan(skew*tr*np.sqrt((skewvar/skew)**2+(
                                                trvar/tr)**2))*180/np.pi)
            
            #calculate angle alpha corresponding to phase tensor's 
            #dependence on coordinate system
            alpha=.5*(np.arctan2(revskew,revtr)*180/np.pi)
            alphavar=abs(.5*np.arctan(revskew*revtr*np.sqrt(
                    (revskewvar/revskew)**2+(revtrvar/revtr)**2))*180/np.pi)
            
            #calculate azimuth as angle between major axis and x-axis
            azimuth=alpha-beta
            azimuthvar=np.sqrt(alphavar**2+betavar**2)
            
            #calulate maximum value for phi
            phimax=np.sqrt((.5*tr)**2+(.5*skew)**2)+\
                np.sqrt(abs((.5*tr)**2+(.5*skew)**2-np.sqrt(phidet)**2))
            phimaxvar=.5*np.sqrt(2*trvar**2+2*skewvar**2)+.5*np.sqrt(
                                        2*trvar**2+2*skewvar**2+phidetvar)
            
            #calculate minimum value for phi
            if np.linalg.det(phi)>=0:
                phimin=np.sqrt((.5*tr)**2+(.5*skew)**2)-\
                np.sqrt((.5*tr)**2+(.5*skew)**2-np.sqrt(phidet)**2)
                
            elif np.linalg.det(phi)<0:
                phimin=np.sqrt((.5*tr)**2+(.5*skew)**2)-np.sqrt(abs(
                                    (.5*tr)**2+(.5*skew)**2-np.sqrt(phidet)**2))
                            
            #set the variance in phimin to that of phimax
            phiminvar=phimaxvar
            
            #calculate ellipticity
            phiminang=(180/np.pi)*np.arctan(np.array(phimin))
            phiminangvar=(180/np.pi)*np.arctan(np.array(phiminvar))
            phimaxang=(180/np.pi)*np.arctan(np.array(phimax))
            phimaxangvar=(180/np.pi)*np.arctan(np.array(phimaxvar))
            
            ellipticity=(phimaxang-phiminang)/(phimaxang+phiminang)
            ellipticityvar=ellipticity*np.sqrt(phimaxangvar**2+
                    phiminangvar**2)*np.sqrt((1/(phimaxang-phiminang))**2+
                    (1/(phimaxang+phiminang))**2)
            
            phidet=phimin*phimax
            phidetvar=phiminvar*phimaxvar
            #append things to a dictionary of lists
            
            self.phi[ii]=phi
            self.phivar[ii]=np.array(dphi)
            self.phimin[ii]=float(phimin)
            self.phiminvar[ii]=float(phiminvar)
            self.phimax[ii]=float(phimax)
            self.phimaxvar[ii]=float(phimaxvar)
            self.alpha[ii]=float(alpha)
            self.alphavar[ii]=float(alphavar)
            self.beta[ii]=float(beta)
            self.betavar[ii]=float(betavar)
            self.azimuth[ii]=float(azimuth)
            self.azimuthvar[ii]=float(azimuthvar)
            self.phiminang[ii]=float(phiminang)
            self.phiminangvar[ii]=float(phiminangvar)
            self.phimaxang[ii]=float(phimaxang)
            self.phimaxangvar[ii]=float(phimaxangvar)
            self.ellipticity[ii]=float(ellipticity)
            self.ellipticityvar[ii]=float(ellipticityvar)
            self.phidet[ii]=float(phidet)
            self.phidetvar[ii]=float(phidetvar)

class ResPhase:
    """
    ResPhase is a resistivity and phase class
    
    Arguments:
    ----------
        **z** : complex np.array (nf,2,2)
                impedance tensor [[zxx,zxy],[zyx,zyy]]
        
        **period** : np.array (nf)
                    array of periods corresponding to the impedance tensor
        
        **zvar** : np.array (nf,2,2)
                   impedance tensor variances [[zxx,zxy],[zyx,zyy]]
        
        **rotz** : float (angle in degrees)
                   rotation angle clockwise positive assuming 0 is North.
        
        **ffactor** : float
                      factor to scale the apparent resistivities by.
                
    """
    
    def __init__(self,z,period,zvar=None,rotz=0,ffactor=1):
        try:
            nz=len(z)
        except TypeError:
            z=np.array([z])
            nz=len(z)
            period=np.array([period])
        
        #make a copy of z    
        self.z=z.copy()
        if zvar==None:
            self.zvar=np.zeros_like(z)
        else:
            try:
                self.zvar=zvar.copy()
            except AttributeError:
                zvar=np.array([zvar])                
                self.zvar=zvar.copy()
        self.z=z.copy()
        self.period=period.copy()
        self.thetar=rotz
        self.ffactor=ffactor
        
        nz=len(period)
        
        self.resxy=np.zeros(nz)
        self.resxyerr=np.zeros(nz)
        self.resyx=np.zeros(nz)
        self.resyxerr=np.zeros(nz)
        self.resxx=np.zeros(nz)
        self.resxxerr=np.zeros(nz)
        self.resyy=np.zeros(nz)
        self.resyyerr=np.zeros(nz)
        
        self.phasexy=np.zeros(nz)
        self.phasexyerr=np.zeros(nz)
        self.phaseyx=np.zeros(nz)
        self.phaseyxerr=np.zeros(nz)
        self.phasexx=np.zeros(nz)
        self.phasexxerr=np.zeros(nz)
        self.phaseyy=np.zeros(nz)
        self.phaseyyerr=np.zeros(nz)
        
        self.resdet=np.zeros(nz)
        self.resdeterr=np.zeros(nz)
        self.phasedet=np.zeros(nz)
        self.phasedeterr=np.zeros(nz)
        
        #rotate the data if desired
        if self.thetar!=0:
            #convert thetar into radians
            if abs(self.thetar)>2*np.pi:
                self.thetar=self.thetar*(np.pi/180)
            #make rotation matrix
            rotmatrix=np.array([[np.cos(self.thetar), np.sin(self.thetar)],
                         [-np.sin(self.thetar), np.cos(self.thetar)]])
            #rotate the data
            for rr in range(nz):
                self.z[rr]=np.dot(rotmatrix,np.dot(self.z[rr],rotmatrix.T))
                self.zvar[rr]=np.dot(rotmatrix,np.dot(self.zvar[rr],rotmatrix.T))
        
        #calculate resistivity and phase components
        for jj in range(nz):
            wt=(.2*self.period[jj])*ffactor
            self.resxx[jj]=wt*abs(self.z[jj,0,0])**2
            self.resxy[jj]=wt*abs(self.z[jj,0,1])**2
            self.resyx[jj]=wt*abs(self.z[jj,1,0])**2
            self.resyy[jj]=wt*abs(self.z[jj,1,1])**2
            
            self.resxxerr[jj]= 2 * wt * abs(self.z[jj,0,0]) * np.sqrt(self.zvar[jj,0,0])
            self.resxyerr[jj]= 2 * wt * abs(self.z[jj,0,1]) * np.sqrt(self.zvar[jj,0,1])
            self.resyxerr[jj]= 2 * wt * abs(self.z[jj,1,0]) * np.sqrt(self.zvar[jj,1,0])
            self.resyyerr[jj]= 2 * wt * abs(self.z[jj,1,1]) * np.sqrt(self.zvar[jj,1,1])
            
            self.phasexx[jj]=(np.arctan2(self.z[jj,0,0].imag,
                                        self.z[jj,0,0].real)*(180/np.pi))%360
            self.phasexy[jj]=(np.arctan2(self.z[jj,0,1].imag,
                                        self.z[jj,0,1].real)*(180/np.pi))%360
            # if not 0 <= self.phasexy[jj] <= 90:
            #     print self.phasexy[jj]
            #     self.phasexy[jj] = ((self.phasexy[jj]))%90
            self.phaseyx[jj]=(np.arctan2(self.z[jj,1,0].imag,
                                        self.z[jj,1,0].real)*(180/np.pi))%360
            # if not 180 <= self.phaseyx[jj] <= 270:
            #     self.phaseyx[jj] = ((self.phaseyx[jj])+180)%360

            self.phaseyy[jj]=(np.arctan2(self.z[jj,1,1].imag,
                                        self.z[jj,1,1].real)*(180/np.pi))%360
            
            self.phasexxerr[jj]=np.sqrt(self.zvar[jj,0,0]/((self.z[jj,0,0].imag)**2+(self.z[jj,0,0].real)**2))*(180/np.pi)
            self.phasexyerr[jj]=np.sqrt(self.zvar[jj,0,1]/((self.z[jj,0,1].imag)**2+(self.z[jj,0,1].real)**2))*(180/np.pi)
            self.phaseyxerr[jj]=np.sqrt(self.zvar[jj,1,0]/((self.z[jj,1,0].imag)**2+(self.z[jj,1,0].real)**2))*(180/np.pi)
            self.phaseyyerr[jj]=np.sqrt(self.zvar[jj,1,1]/((self.z[jj,1,1].imag)**2+(self.z[jj,1,1].real)**2))*(180/np.pi)


            #old
            # self.phasexxerr[jj]=np.arcsin(self.zvar[jj,0,0]/
            #                                 abs(self.z[jj,0,0]))*(180/np.pi)
            # self.phasexyerr[jj]=np.arcsin(self.zvar[jj,0,1]/
            #                                 abs(self.z[jj,0,1]))*(180/np.pi)
            # self.phaseyxerr[jj]=np.arcsin(self.zvar[jj,1,0]/
            #                                 abs(self.z[jj,1,0]))*(180/np.pi)
            # self.phaseyyerr[jj]=np.arcsin(self.zvar[jj,1,1]/
            #                                 abs(self.z[jj,1,1]))*(180/np.pi)
            
            #calculate determinant values
            #apparent resistivity
            zdet=np.sqrt(np.linalg.det(self.z[jj]))
            zdetvar=np.linalg.det(self.zvar[jj])
            self.resdet[jj]=wt*abs(zdet)**2
            self.resdeterr[jj]=wt*np.abs(zdet+zdetvar)**2-self.resdet[jj]
            
            #phase
            self.phasedet[jj]=np.arctan2(zdet.imag,zdet.real)*(180/np.pi)
            self.phasedeterr[jj]=np.arcsin(zdetvar/abs(zdet))*(180/np.pi)

class Tipper:
    """
    Tipper Class
    
    Arguments:
    ----------
        **tipper** : complex np.array (nf,2)
                     tipper array [tx,ty]
                     
        **roty** : float (angle in degrees)
                   rotation angle clockwise positive assuming 0 is North.        
        
    """

    def __init__(self,tipper,rott=0):
        
        try:
            nt=len(tipper)
        except TypeError:
            tipper=np.array([tipper])
            nt=len(tipper)
            
        self.tipper=tipper.copy()
        self.thetar=rott
        
         #rotate the data if desired
        if self.thetar!=0:
            #convert thetar into radians
            if abs(self.thetar)>2*np.pi:
                self.thetar=self.thetar*(np.pi/180)
            #make rotation matrix
            rotmatrix=np.array([[np.cos(self.thetar), np.sin(self.thetar)],
                         [-np.sin(self.thetar), np.cos(self.thetar)]])
            #rotate the data
            for rr in range(nt):
                self.tipper[rr]=np.dot(rotmatrix,np.dot(self.tipper[rr],
                                        rotmatrix.T))
        #get the magnitude
        self.magreal=np.sqrt(self.tipper[:,0].real**2+self.tipper[:,1].real**2)
        self.magimag=np.sqrt(self.tipper[:,0].imag**2+self.tipper[:,1].imag**2)
        
        #get the angle, need to make both parts negative to get it into the
        #parkinson convention where the arrows point towards the conductor
        self.anglereal=np.arctan2(-self.tipper[:,1].real,
                                        -self.tipper[:,0].real)*180/np.pi
        self.angleimag=np.arctan2(-self.tipper[:,1].imag,
                                            -self.tipper[:,0].imag)*180/np.pi

class Zinvariants:
    """
    calculates invariants from Weaver et al. [2003]
    
    Arguments:
    ----------
    
        **z** : complex np.array (nf,2,2)
                impedance tensor
    """
    
    def __init__(self,z,rotz=0):
        
        try:
            nz=len(z)
        except TypeError:
            z=np.array([z])
            nz=len(z)
            
        self.z=np.array(z)
        self.thetar=rotz
        self.inv1=np.zeros(nz)
        self.inv2=np.zeros(nz)
        self.inv3=np.zeros(nz)
        self.inv4=np.zeros(nz)
        self.inv5=np.zeros(nz)
        self.inv6=np.zeros(nz)
        self.inv7=np.zeros(nz)
        self.q=np.zeros(nz)
        self.strike=np.zeros(nz)
        self.strikeerr=np.zeros(nz)
        
        #rotate the data if desired
        if self.thetar!=0:
            #convert thetar into radians
            if abs(self.thetar)>2*np.pi:
                self.thetar=self.thetar*(np.pi/180)
            #make rotation matrix
            rotmatrix=np.array([[np.cos(self.thetar), np.sin(self.thetar)],
                         [-np.sin(self.thetar), np.cos(self.thetar)]])
            #rotate the data
            for rr in range(nz):
                self.z[rr]=np.dot(rotmatrix,np.dot(self.z[rr],rotmatrix.T))
        
        for ii in range(nz):
            #compute invariant parameters
            x1=.5*(self.z[ii,0,0].real+self.z[ii,1,1].real)
            x2=.5*(self.z[ii,0,1].real+self.z[ii,1,0].real)
            x3=.5*(self.z[ii,0,0].real-self.z[ii,1,1].real)
            x4=.5*(self.z[ii,0,1].real-self.z[ii,1,0].real)
            e1=.5*(self.z[ii,0,0].imag+self.z[ii,1,1].imag)
            e2=.5*(self.z[ii,0,1].imag+self.z[ii,1,0].imag)
            e3=.5*(self.z[ii,0,0].imag-self.z[ii,1,1].imag)
            e4=.5*(self.z[ii,0,1].imag-self.z[ii,1,0].imag)
            ex=x1*e1-x2*e2-x3*e3+x4*e4
            d12=(x1*e2-x2*e1)/ex
            d34=(x3*e4-x4*e3)/ex
            d13=(x1*e3-x3*e1)/ex
            d24=(x2*e4-x4*e2)/ex
            d41=(x4*e1-x1*e4)/ex
            d23=(x2*e3-x3*e2)/ex
            inv1=np.sqrt(x4**2+x1**2)
            inv2=np.sqrt(e4**2+e1**2)
            inv3=np.sqrt(x2**2+x3**2)/inv1
            inv4=np.sqrt(e2**2+e3**2)/inv2
            s12=(x1*e2+x2*e1)/ex
            s34=(x3*e4+x4*e3)/ex
            s13=(x1*e3+x3*e1)/ex
            s24=(x2*e4+x4*e2)/ex
            s41=(x4*e1+x1*e4)/ex
            s23=(x2*e3+x3*e2)/ex
            inv5=s41*ex/(inv1*inv2)
            inv6=d41*ex/(inv1*inv2)
            q=np.sqrt((d12-d34)**2+(d13+d24)**2)
            inv7=(d41-d23)/q
            strikeang=.5*np.arctan2(d12-d34,d13+d24)*(180/np.pi)
            strikeangerr=abs(.5*np.arcsin(inv7))*(180/np.pi)
            self.inv1[ii]=inv1
            self.inv2[ii]=inv2
            self.inv3[ii]=inv3
            self.inv4[ii]=inv4
            self.inv5[ii]=inv5
            self.inv6[ii]=inv6
            self.inv7[ii]=inv7
            self.q[ii]=q
            self.strike[ii]=strikeang
            self.strikeerr[ii]=strikeangerr
            
class ResistivityTensor:
    """
    gets components of the resistivity tensor defined by 
    O'reilly [1978] and Weckmann et al. [2003]
    
    Arguments:
    ----------
        **z** : complex np.array (nf,2,2)
                impedance tensor [[zxx,zxy],[zyx,zyy]]
        
        **frequency** : np.array (nf)
                    array of frequencies corresponding to the impedance tensor
        
        **rotate** : [ 90 | 180 | 270 ] 
                     Caldwell et al, [2004] assume the coordinate axis is 
                     Y North and X East.  If the data is in X North and Y 
                     East than rotation = 180. *Default* is 180
                     
        **rotz** : float (angle in degrees)
                   rotation angle clockwise positive assuming 0 is North.
    
    """
    
    def __init__(self,z,frequency,rotate=180,rotz=0):
        
        try:
            nz=len(z)
        except TypeError:
            z=np.array([z])
            frequency=np.array([frequency])
            nz=len(z)
        self.thetar=rotz
        self.frequency=frequency.copy()
        
        #rotate impedance tensor to align x pointing 
        #east and y pointing north following convention of Weckmann [2003]
        if rotate==90:
            self.z=np.array([np.rot90(zz,1) for zz in z])
        elif rotate==180:
            self.z=np.array([np.rot90(zz,2) for zz in z])
        elif rotate==270:
            self.z=np.array([np.rot90(zz,3) for zz in z])
        else:
            self.z=z.copy()    
        
        #rotate the data if desired
        if self.thetar!=0:
            #convert thetar into radians
            if abs(self.thetar)>2*np.pi:
                self.thetar=self.thetar*(np.pi/180)
            #make rotation matrix
            rotmatrix=np.array([[np.cos(self.thetar), np.sin(self.thetar)],
                         [-np.sin(self.thetar), np.cos(self.thetar)]])
            #rotate the data
            for rr in range(nz):
                self.z[rr]=np.dot(rotmatrix,np.dot(self.z[rr],rotmatrix.T))

        Y=np.zeros_like(self.z)
        for ii,zz in enumerate(self.z):
            try:            
                Y[ii]=np.linalg.inv(zz)
            except np.linalg.LinAlgError:
                pass
            
        
        self.gamma=np.zeros_like(self.z)
        self.gamma[:,0,0]=-frequency**2*(Y[:,1,1]*Y[:,0,0]-Y[:,1,0]*Y[:,1,0])
        self.gamma[:,0,1]=-frequency**2*(Y[:,1,1]*(Y[:,0,1]-Y[:,1,0]))
        self.gamma[:,1,0]=-frequency**2*(Y[:,0,0]*(Y[:,1,0]-Y[:,0,1]))
        self.gamma[:,1,1]=-frequency**2*(Y[:,1,1]*Y[:,0,0]-Y[:,0,1]*Y[:,0,1])
        
        self.rho=np.zeros_like(self.gamma.imag)
        for ii,gg in enumerate(self.gamma):
            try:                
                self.rho[ii]=.2*frequency[ii]*np.linalg.inv(gg.imag) 
            except np.linalg.LinAlgError:
                pass
                                        
        self.rhodet=np.sqrt(abs(np.array([np.linalg.det(rr) 
                            for rr in self.rho])))
        
        pi1=.5*np.sqrt((self.rho[:,0,0]-self.rho[:,1,1])**2+\
            (self.rho[:,0,1]+self.rho[:,1,0])**2)
        pi2=.5*np.sqrt((self.rho[:,0,0]+self.rho[:,1,1])**2+\
            (self.rho[:,0,1]-self.rho[:,1,0])**2)
        
        self.rhomax=pi1+pi2
        self.rhomin=pi2-pi1
		
        self.rhoalpha=.5*np.arctan((self.rho[:,0,1]+self.rho[:,1,0])/
                                 (self.rho[:,0,0]-self.rho[:,1,1]))*(180/np.pi)
        self.rhobeta=.5*np.arctan((self.rho[:,0,1]-self.rho[:,1,0])/
                                 (self.rho[:,0,0]-self.rho[:,1,1]))*(180/np.pi)
        self.rhoazimuth=self.rhoalpha-self.rhobeta
        
        self.eps=np.zeros_like(self.gamma.real)
        for ii,gg in enumerate(self.gamma):
            try:                
                self.eps[ii]=.2*frequency[ii]*np.linalg.inv(gg.real) 
            except np.linalg.LinAlgError:
                pass
            
        self.epsdet=np.sqrt(abs(np.array([np.linalg.det(rr) 
                            for rr in self.eps])))
        
        self.epsmax=.5*np.sqrt((self.eps[:,0,0]-self.eps[:,1,1])**2+
                                (self.eps[:,0,1]+self.eps[:,1,0])**2)
        self.epsmin=.5*np.sqrt((self.eps[:,0,0]+self.eps[:,1,1])**2+
                                (self.eps[:,0,1]-self.eps[:,1,0])**2)
        self.epsalpha=.5*np.arctan((self.eps[:,0,1]+self.eps[:,1,0])/
                                 (self.eps[:,0,0]-self.eps[:,1,1]))
        self.epsbeta=.5*np.arctan((self.eps[:,0,1]-self.eps[:,1,0])/
                                 (self.eps[:,0,0]-self.eps[:,1,1]))
        self.epsazimuth=self.epsalpha-self.epsbeta
        
class PhaseTensorResidual:
    """
    calculates the tensor residual
    as defined by 
    
    .. math:: 
        \Delta \Phi_{1,2} = \hat{I} - \Phi_1^{-1} \Phi_2
    """
    
    def __init__(self,z1,z2,rotz=0,rotate=180):
        z1=np.array(z1)        
        z2=np.array(z2)        
        if len(z1.shape)==2:        
            self.z1=np.array([z1])
            self.z2=np.array([z2])
        else:
            self.z1=z1
            self.z2=z2
            
        nf,xx,yy=z1.shape
        
        self.thetar=rotz
        self.rotatecoord=rotate
        self.phi=np.zeros((nf,xx,yy))
        self.phimin=np.zeros(nf)
        self.phimax=np.zeros(nf)
        self.azimuth=np.zeros(nf)
        self.beta=np.zeros(nf)
        self.ecolor=np.zeros(nf)
        
        pt1=PhaseTensor(z1,rotz=rotz,rotate=rotate)
        pt2=PhaseTensor(z1,rotz=rotz,rotate=rotate)
        
        #calculate the difference between the two phase tensor ellipses
        for ii in range(nf):
#            phi=np.eye(2)-\
#                (np.dot(np.linalg.inv(pt2.phi[ii]),pt1.phi[ii])+
#                np.dot(pt1.phi[ii],np.linalg.inv(pt2.phi[ii])))/2             
            phi=np.eye(2)-np.dot(np.linalg.inv(pt1.phi[ii]),pt2.phi[ii])             
            
            #compute the trace        
            tr=phi[0,0]+phi[1,1]
            #Calculate skew of phi and the cooresponding error
            skew=phi[0,1]-phi[1,0]
            #calculate the determinate and determinate error of phi
            phidet=abs(np.linalg.det(phi))
            
            #calculate reverse trace and error
            revtr=phi[0,0]-phi[1,1]
            
            #calculate reverse skew and error
            revskew=phi[1,0]+phi[0,1]
            
            beta=.5*np.arctan2(skew,tr)*(180/np.pi)
            alpha=.5*np.arctan2(revskew,revtr)*(180/np.pi)
            
            #need to figure out why azimuth is off by 90 deg
            azimuth=(alpha-beta)                   
           
            #calculate phimax
            phimax=np.sqrt(abs((.5*tr)**2+(.5*skew)**2))+\
                    np.sqrt(abs((.5*tr)**2+(.5*skew)**2-np.sqrt(phidet)**2))
                
            #calculate minimum value for phi
            if phidet>=0:
                phimin=np.sqrt(abs((.5*tr)**2+(.5*skew)**2))-\
                np.sqrt(abs((.5*tr)**2+(.5*skew)**2-np.sqrt(phidet)**2))
            elif phidet<0:
                phimin=-1*np.sqrt(abs((.5*tr)**2+(.5*skew)**2))-np.sqrt(abs(
                            (.5*tr)**2+(.5*skew)**2-(np.sqrt(phidet))**2))
                            
            ecolor=(abs(phi.min())+abs(phi.max()))/2
            
            #put things into arrays
            self.phimax[ii]=phimax
            self.phimin[ii]=phimin
            self.azimuth[ii]=azimuth
            self.beta[ii]=beta
            self.ecolor[ii]=ecolor
            self.phi[ii]=phi
        
class ResistivityTensorResidual:
    """
    will calculate the resistivity residual between two tensors
    
    .. math:: 
        \Delta \Phi_{1,2} = \hat{I} - \Phi_1^{-1} \Phi_2
    """
    
    def __init__(self,z1,z2,frequency,rotz=0,rotatecoord=180):
        z1=np.array(z1)        
        z2=np.array(z2)        
        if len(z1.shape)==2:        
            self.z1=np.array([z1])
            self.z2=np.array([z2])
        else:
            self.z1=z1
            self.z2=z2
            
        nf,xx,yy=z1.shape
        self.thetar=rotz
        self.rotatecoord=rotatecoord
        self.rho=np.zeros((nf,xx,yy))
        self.rhomin=np.zeros(nf)
        self.rhomax=np.zeros(nf)
        self.azimuth=np.zeros(nf)
        self.rhodet=np.zeros(nf)
        self.ecolor=np.zeros(nf)
        
        rt1=ResistivityTensor(z1,frequency,rotz=rotz,rotate=rotatecoord)
        rt2=ResistivityTensor(z2,frequency,rotz=rotz,rotate=rotatecoord)
        
        for ii in range(nf):
            rho=np.eye(2)-\
                    (np.dot(np.linalg.inv(rt2.rho[ii]),rt1.rho[ii])+
                    np.dot(rt1.rho[ii],np.linalg.inv(rt2.rho[ii])))/2
                    
            pi1=.5*np.sqrt((rho[0,0]-rho[1,1])**2+(rho[0,1]+rho[1,0])**2)
            pi2=.5*np.sqrt((rho[0,0]+rho[1,1])**2+(rho[0,1]-rho[1,0])**2)
            
            rhomax=pi1+pi2
            rhomin=pi2-pi1
            
            alpha=.5*np.arctan((rho[0,1]+rho[1,0])/(rho[0,0]-rho[1,1]))
            beta=.5*np.arctan((rho[0,1]-rho[1,0])/(rho[0,0]+rho[1,1]))
            
            azimuth=(alpha-beta)*180/np.pi
            
            ecolor=np.sign(rt1.rhomax[ii]-rt2.rhomin[ii])*\
                    (abs(rho.min())+abs(rho.max()))/2
                    
            rhodet=rt1.rhodet[ii]-rt2.rhodet[ii]
                    
            #put things into arrays
            self.rhomax[ii]=rhomax
            self.rhomin[ii]=rhomin
            self.azimuth[ii]=azimuth
            self.ecolor[ii]=ecolor
            self.rho[ii]=rho
            self.rhodet[ii]=rhodet
