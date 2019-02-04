#!/usr/bin/env python

#-------------------------------------------------------------------------
"""
Loop for running birrp processing on one station over all days in folder:

The BIRRP string can be edited based on the BIRRP mode of interest (basic/advanced)


The input data directory must contain subdirectories, sorted and named by days

the name of the subdirectories must end with a date, separated from the rest
of the name by an underscore
date format: YYMMDD (e.g. 140320 for the 20th of March 2014)

"""
#-------------------------------------------------------------------------


import os
import sys
import os.path as op
import numpy as np
import subprocess
import pdb

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------

inputdatadir = 'indataL206'
outputdir = 'L206_birrpoutput'

station = 'L206'

#birrp_exe = 'birrp5_linux32_v2.exe'
birrp_exe = '/stash/Working_Scripts/BIRRP/birrp52'

notch_frequencies = [-46.875, -93.750, -156.25]

channels = ['n', 'e']

#------------------------------------------------------------------------------

# set up generic BIRRP input string - advanced mode - simple settings
birrp_string = """1
2
2
2
0
2
-500
65536,2,12
3,1,3
y
2
0,0.0001,0.9999
0
%s
0
%d
%s1
15
0
0
1000000
0
%s/%s.staA.%s
0
0
%s/%s.staA.%s
0
0
%s/%s.staB.north
0
0
%s/%s.staB.east
0
0
%s/%s.staC.north
0
0
%s/%s.staC.east
0
%d,%d,180
0,90,0
0,90,0
"""


#------------------------------------------------------------------------------
#  do not edit below this line-------------------------------------------------
#------------------------------------------------------------------------------


print()

basedir = op.abspath(os.curdir)
outdir = op.join(basedir, outputdir)

if not op.isdir(outdir):
    os.makedirs(outdir)

indir = op.join(basedir, inputdatadir)
if not op.isdir(indir):
    sys.exit('\n\tERROR - input data directory does not exist: %s' % (indir))

channel_dict = {
    'n': 0,
    'e': 1,
    's': 2,
    'w': 3,
    '0': 'n',
    '1': 'e',
    '2': 's',
    '3': 'w'}

longchannels = ['north', 'east', 'south', 'west']

xchannel_idx = channel_dict[channels[0]]
ychannel_idx = channel_dict[channels[1]]

xchannel = longchannels[xchannel_idx]
ychannel = longchannels[ychannel_idx]

xdegrees = int(90 * xchannel_idx)
ydegrees = int(90 * ychannel_idx)

notchstring = ''
for freq in notch_frequencies:
    notchstring += '%.3f\n' % (freq)


os.chdir(indir)

lo_subdirs = os.listdir('.')
lo_subdirs = [i for i in lo_subdirs if op.isdir(i)]

os.chdir(basedir)

for block in lo_subdirs:
    date = block.split('_')[-1]
    try:
        int(float(date))

    except:
        print('\t Warning - not a valid naming format for subdirectory %s \n' % (block))
        continue

    current_indir = op.join(indir, block)

    current_outdir = op.join(outdir, date)
    if not op.isdir(current_outdir):
        os.makedirs(current_outdir)

    os.chdir(current_outdir)

    relative_indir = '../../' + inputdatadir + '/' + block
    outdata_name = station + '_' + date

    # find filebase:
    current_indir_files = os.listdir(current_indir)
    files_of_interest = [
        i for i in current_indir_files if i.lower().endswith('timestamps')]
    filebasename = op.splitext(files_of_interest[0])[0]

    # abbreviation for easier handling:
    bn = filebasename

    print('...processing station %s, date: %s ' % (station, date))

    current_birrp_string = birrp_string % (outdata_name, len(notch_frequencies), notchstring,
                                           relative_indir, bn, xchannel,
                                           relative_indir, bn, ychannel,
                                           relative_indir, bn,
                                           relative_indir, bn,
                                           relative_indir, bn,
                                           relative_indir, bn,
                                           xdegrees, ydegrees)
    # print current_birrp_string
    if 1:
        P = subprocess.Popen(
            [birrp_exe],
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE)
        o, e = P.communicate(current_birrp_string)
        print('\t...Done!\n')
    # except:
    #     print '\t...ERROR - processing failed!\n'
    # pdb.set_trace()

    os.chdir(outdir)
print('Processing outputs in directory %s' % (outdir))
print()

os.chdir(basedir)
