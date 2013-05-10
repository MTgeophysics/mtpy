#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
This file will combine two edi files into one a number of files

If you want to see where the edi files match in frequency use:
mtplot.resPhasePlots([edi1,edi2],plottype=3)

where edi1 is the full path to the high frequency edi file and 
edi2 is the full path to the low frequency edi file

Then you can match up the frequencies, noting that the program counts
from the high frequencies down for the high frequency file and from low
frequencies to high for the low frequency edi file.

so if you want the first 12 frequencies from edi1 and the last 17 
frequencies from edi2, then n1=12 and n2=17
"""

import os
import os.path as op
import sys
import mtpy.imaging.mtplottools as mtplot
import mtpy.core.mttools as mt


def main():

    arglist = sys.argv[1:]

    if len(arglist) < 2:
        sys.exit('ERROR -- provide 2 arguments: <high frequency edi directory path> <low frequency edi directory path >')

    #enter the high frequency edi directory path 
    edipath1 = op.abspath(arglist[0])

    #enter the low frequency edi directory path 
    edipath2 = op.abspath(arglist[1])

    #station id is the index of the first and last character in the station
    #name which will be matched
    #for station 'pb01' the station id is 0 for 'p' and 3 for '1'
    #but python only indexes to that last number minus 1 ergo stationid
    #for 'pb01' would be (0,4)
    stationid=(0,4)
    s1,s2=stationid[0],stationid[1]
    #number of frequencies to read from high frequency data, starting from
    #high frequencies down
    n1=12

    #number of frequencies to read from low frequency data, starting from 
    #low frequencies and going up
    n2=12

    #make a new list to put the new edi file names into
    nedilst=[]
    for edihf in os.listdir(edipath1):
        cedilst=[]
        for edibb in os.listdir(edipath2):
            if edihf[s1:s2]==edibb[s1:s2]:
                cedilst.append(os.path.join(edipath1,edihf))
                cedilst.append(os.path.join(edipath2,edibb))
                break
        nedi=mt.combineEdifiles(cedilst[0],cedilst[1],n1,n2)
        nedilst.append(nedi)

    #if you want to plot the new edi file to see if the frequencies match
    #use (uncomment the following for a loop
    #for ii,nedi in enumerate(nedilst,1):
    #	mtplot.plotResPhase(nedi,fignum=ii)

    #or for just one find the index in edilst and use
    #mtplot.plotResPhase(nedilst[index])

if __name__ == '__main__':
    main()
