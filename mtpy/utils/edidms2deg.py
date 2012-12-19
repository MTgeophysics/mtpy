#!/usr/bin/env python

import sys
import os
import os.path as op
import glob


def dms2deg(dms_in):
    """
    ...

    """
    raw_in = dms_in.split(':')

    sec_fracs = float(raw_in[2])/60.
    min_fracs = (float(raw_in[1])+sec_fracs)/60.
    degs =  abs(float(raw_in[0])) +  min_fracs
    
    if float(raw_in[0]) < 0:
        degs = -degs
    
    return str(degs)


def deg2dms(deg_in):
    """
    ...

    """
 

    raw_in = float(deg_in)
    degs = int(raw_in)
    degfracs = abs(raw_in - degs )
    mins = int(degfracs*60.)
    minfracs = degfracs*60. - mins
    secs = minfracs *60

    dms_out = '%i:%i:%.2f'%(degs,mins,secs)
    
    return dms_out

def check_format(instring):
    """
    ...

    """
 
    try:
        a = float(instring)
        return 'deg'
    except: 
        stringformat = None

    try:
        b = instring.split(':')
        return 'dms'
    except:
        return None



def edidms2deg():
    """
    edidms2deg converts the LAT/LON information within given EDI files from/to the degree/deg:min:sec format.

    Call it with (a list of) files as arguments. Wildcards are
    allowed. The repective lines within the files are converted and
    the files written to the current working directory. The names are
    not changed, so existing files get overwritten!!
    
    Default: conversion from deg:min:sec into degrees

    To change the conversion direction, apend the letter 'i' at the end of the arguments

    Example:

    edidms2deg examplefile1.edi examplefile2.edi morefiles*.edi 

        or 
      
    edidms2deg otherwayround1.edi  moredegreefilesfiles*.edi  i  
    

    """


    invert = 0
    args = sys.argv
    usage = "\n usage: \n \t edidms2deg edi-file(s) [optional 'i' for inverted transformation]\n"
    if len(args) < 2 : 
        sys.exit(usage)

    lo_files_raw = args[1:]
    try:
        if lo_files_raw[-1].lower() == 'i':
            invert = 1
            lo_files_raw.pop()
    except:
        sys.exit(usage+'\n last argument must be filename or the letter "i"')

 
    lo_files =[]
    
    #walk through all arguments, using glob allows for wildcard usage
    
    for listarg in lo_files_raw:
        tmp_list = glob.glob(listarg)
        
        for arg in tmp_list:
            try:
                fn = op.abspath(arg)
                lo_files.append(fn)
            except:
                print '%s is no valid file'%arg

    for fn in lo_files:
        out_string = ''
        latlon = 0
        F = open(fn,'r')
        try:
            for curr_line_raw in  F.readlines():
                curr_line = curr_line_raw.strip()
                out_line  = curr_line
                raw_line1 = curr_line.split('=')


                if len(raw_line1)>1:
                    if  raw_line1[0].strip()[:3].lower() =='lat':
                        lat_in = raw_line1[1]
                        if check_format(lat_in) == 'dms' and invert == 0: 
                            lat_out = dms2deg(lat_in)
                            out_line = 'LAT=%s'%(lat_out)
                            latlon +=1
                        elif check_format(lat_in) == 'deg' and invert == 1:
                            lat_out = deg2dms(lat_in)
                            out_line = 'LAT=%s'%(lat_out)
                            latlon +=1

                        
                    elif raw_line1[0].strip()[:3].lower() =='lon': 
                        lon_in = raw_line1[1]
                        if check_format(lon_in) == 'dms' and invert == 0: 
                            lon_out = dms2deg(lon_in)
                            out_line = 'LON=%s'%(lon_out)
                            latlon +=1
                        elif check_format(lon_in) == 'deg' and invert == 1:
                            lon_out = deg2dms(lon_in)
                            out_line = 'LON=%s'%(lon_out)
                            latlon +=1
                out_string +=out_line
                out_string +='\n'

               
            F.close()

            if latlon == 2:

                F2 = open(op.basename(fn),'w')
                F2.write(out_string)
                F2.close()
            
            else:
                print "\t %s did not contain proper lat/lon information - maybe it's just not an EDI file... \n\t ...or it's in the correct form already...(file left unchanged)"%(fn)
                continue

        except:
            F.close()
            print "\t could not parse file %s ... maybe it's not an EDI file...\n"%(fn)
            continue
        
        

 

if __name__ == '__main__':

    edidms2deg()
