#!/usr/bin/env python

import os,sys,shutil
import os.path as op

import mtpy.core.edi as MTedi
import mtpy.utils.convert_birrp_output as MTbp
import pdb
import numpy as np
import mtpy.utils.exceptions as MTex
import mtpy.uofa.qel_monitoring_j2edi as qel2edi
import mtpy.utils.edi2columnsonly as edi2col
import mtpy.uofa.simpleplotEDI as smplplt


#indir = 'L09_before_23Feb_birrpoutput'
#indir = 'test'
indir = 'L102_after_23Feb_birrpoutput'

outdir = 'qel_collected_L102_after'

station = 'L102'

survey_configfile= op.abspath('/data/temp/nigel/romasurvey.cfg')
instr_resp = op.abspath('/data/mtpy/mtpy/uofa/lemi_coils_instrument_response_freq_real_imag_microvolts.txt')

outdir = op.join(op.abspath(outdir),station)

indir = op.abspath(indir)

dirs = os.listdir(indir)

dirs = sorted([i for i in dirs if op.isdir(op.join(indir,i))])

basedir = op.abspath(os.curdir)

for date in dirs:

    daybase = op.abspath(op.join(indir,date))
    os.chdir(daybase)

    if 1:
        fullday = date.split('-')[0]
        day_try = int(float(fullday))
        day = int(float(fullday[-2:]))
        month_num = int(float(fullday[-4:-2]))
        #month_num = {'jan':1,'feb':2,'mar':3,'apr':4,'may':5,'jun':6,
        #                'jul':7,'aug':8,'sep':9,'oct':10,'nov':11,'dec':12,}[month]
    # except:
    #     continue
        
    if 1:
        outfn,outfn_coh = qel2edi.convert2edi(station,'.',survey_configfile,instr_resp,string2strip=['_before','_23Feb'], datestring=fullday)

    # except:
    #     print 'no information found in folder {0}'.format(op.abspath(os.curdir))
    #     pass
    try:
        colfile = edi2col.convert2columns(op.basename(outfn))
    except:
        pass


    outdir_edi = op.join(basedir,outdir,'edi')

    print outfn,outfn_coh,colfile

    if not op.isdir(outdir_edi):
        os.makedirs(outdir_edi)
    try:
        shutil.copy(op.basename(outfn),outdir_edi)
    except:
        pass

    outdir_coh = op.join(basedir,outdir,'coh')
    if not op.isdir(outdir_coh):
        os.makedirs(outdir_coh)

    try:
        shutil.copy(op.basename(outfn_coh),outdir_coh)
    except:
        pass

    outdir_cols = op.join(basedir,outdir,'columns')
    if not op.isdir(outdir_cols):
        os.makedirs(outdir_cols)

    try:
        shutil.copy(op.basename(colfile),outdir_cols)
    except:
        pass

    outdir_plots = op.join(basedir,outdir,'plots')
    if not op.isdir(outdir_plots):
        os.makedirs(outdir_plots)

    try:
        plotfn = smplplt.plotedi(outfn,saveplot=True)
        shutil.copy(op.basename(plotfn),outdir_plots)
    except:
        pass

    donefiles = os.listdir('.')
    donefiles = [i for i in donefiles if i.lower().endswith('.j.done')]
    for df in donefiles:
        os.rename(df,df.replace('.done',''))
    
    #pdb.set_trace()

    os.chdir(indir)

os.chdir(basedir)







