#!/usr/bin/env python

"""

    qel_convert_allstations.py

Convert the outputs of BIRRP. Loop over all stations for a single day.
The directory structure must be:
- input dir (to be set in the header of this file)
 - stationname (at least starting with it and then separated by underscore)
     -date (at least ending with it and  separated by underscore)
       - outputfiles (only one of each file type *.2c2 and *.j will be processed)

Output:
- new directory (to be set in the header of this file)
    - date (to be set in the header of this file)
        - subdirectories 'coh columns edi plots'


Requires:
survey_configfile,
calibrationfile,
plot component dictionary (can be empty)

 all to be set in the header of this file


The EDI file prefix is set in  "qel_monitoring_j2edi.py"

"""


import os
import sys
import shutil
import os.path as op

import mtpy.core.edi as MTedi
import mtpy.utils.convert_birrp_output as MTbp
import numpy as np
import mtpy.utils.exceptions as MTex
import mtpy.uofa.qel_monitoring_j2edi as qel2edi
import mtpy.utils.edi2columnsonly as edi2col
import mtpy.uofa.simpleplotEDI as smplplt
import mtpy.uofa.simpleplotCOH as smplpltCOH

#import pdb

#=========================================================================

#indir = 'BIRRP_Outtape'
#indir ='L2_All_stations_March_Basic_18Mar'
#indir = 'birrp_output'
#indir ='testin'
indir = '.'

#outdir = 'qel_collected_L2_All_stations_March_Basic_18Mar'
outdir = 'testout'

date = '140318'

# 20
plot_component_dict = {}  # 'L209':'e','L213':'e','L224':'e','L218':'n'}
# 21
# plot_component_dict={'L209':'e','L213':'e','L224':'e','L202':'e','L204':'e','L220':'n','L223':'n'}
# 18
# plot_component_dict={'L209':'e','L213':'e','L224':'e','L208':'n','L200':'e','L210':'n','L214':'n'}


survey_configfile = op.abspath('/data/temp/nigel/romasurvey.cfg')
instr_resp = op.abspath(
    '/data/mtpy/mtpy/uofa/lemi_coils_instrument_response_freq_real_imag_normalised.txt')
#instr_resp = op.abspath('/data/mtpy/mtpy/uofa/lemi_coils_instrument_response_freq_real_imag_microvolts.txt')

outdir_prefix = ''

string2strip = ['_RR', '_B125']

#=========================================================================

outdir = op.abspath(outdir)

indir = op.abspath(indir)

dirs = os.listdir(indir)

dirs = sorted([i for i in dirs if op.isdir(op.join(indir, i))])

print(dirs)
basedir = op.abspath(os.curdir)

for stationdir in dirs:

    station = stationdir.split('_')[0]

    stationbase = op.abspath(op.join(indir, stationdir))

    os.chdir(stationbase)

    daydirs = os.listdir('.')
    daydirs = sorted([i for i in daydirs if op.isdir(i)])

    for daydir in daydirs:
        fullday = daydir.split('_')[-1]

        day = int(float(fullday[-2:]))
        month = int(float(fullday[-4:-2]))
        year = int(float(fullday[-6:-4]))

        os.chdir(daydir)
        print(op.abspath(os.curdir))

        donefiles = os.listdir('.')
        donefiles = [i for i in donefiles if i.lower().endswith('.j.done')]
        for df in donefiles:
            os.rename(df, df.replace('.done', ''))

        lo_old_coh_files = os.listdir('.')
        lo_old_coh_files = [
            i for i in lo_old_coh_files if i.lower().endswith('.coh')]
        for i in lo_old_coh_files:
            os.remove(i)
        try:
            outfn, outfn_coh = qel2edi.convert2edi(
                station, '.', survey_configfile, instr_resp, string2strip=string2strip, datestring=fullday)

        except:
            print('no information found in folder {0}'.format(op.abspath(os.curdir)))
            continue
        try:
            colfile = edi2col.convert2columns(op.basename(outfn))
        except:
            pass

        outdir_edi = op.join(
            basedir, outdir, '{0}{1:02d}{1:02d}{2:02d}'.format(
                outdir_prefix, year, month, day), 'edi')

        print(outfn, outfn_coh, colfile)

        if not op.isdir(outdir_edi):
            os.makedirs(outdir_edi)
        try:
            shutil.copy(op.basename(outfn), outdir_edi)
            print('copied EDI file to %s' % (outdir_edi))
        except:
            pass

        outdir_coh = op.join(
            basedir, outdir, '{0}{1:02d}{1:02d}{2:02d}'.format(
                outdir_prefix, year, month, day), 'coh')
        if not op.isdir(outdir_coh):
            os.makedirs(outdir_coh)

        try:
            shutil.copy(op.basename(outfn_coh), outdir_coh)
        except:
            pass

        outdir_cols = op.join(
            basedir, outdir, '{0}{1:02d}{1:02d}{2:02d}'.format(
                outdir_prefix, year, month, day), 'columns')
        if not op.isdir(outdir_cols):
            os.makedirs(outdir_cols)

        try:
            shutil.copy(op.basename(colfile), outdir_cols)
        except:
            pass

        outdir_plots = op.join(
            basedir, outdir, '{0}{1:02d}{1:02d}{2:02d}'.format(
                outdir_prefix, year, month, day), 'plots')
        if not op.isdir(outdir_plots):
            os.makedirs(outdir_plots)

        try:
            plot_component = 'ne'
            if station.upper() in plot_component_dict:
                plot_component = plot_component_dict[station.upper()]

            plotfn = smplplt.plotedi(
                outfn, saveplot=True, component=plot_component)
            shutil.copy(op.basename(plotfn), outdir_plots)
            print('copied res/phase plot %s' % (plotfn))

        except:
            pass

        try:
            plotfncoh = smplpltCOH.plotcoh(outfn_coh, saveplot=True)
            shutil.copy(op.basename(plotfncoh), outdir_plots)
            print('copied coherence plot %s' % (plotfncoh))

        except:
            pass

        donefiles = os.listdir('.')
        donefiles = [i for i in donefiles if i.lower().endswith('.j.done')]
        for df in donefiles:
            os.rename(df, df.replace('.done', ''))

        # pdb.set_trace()
        os.chdir(stationbase)

    os.chdir(indir)

os.chdir(basedir)
print()
