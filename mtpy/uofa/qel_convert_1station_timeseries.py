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
import mtpy.uofa.simpleplotCOH as smplpltCOH


#indir = 'L09_before_23Feb_birrpoutput'
indir = 'test'
indir = 'L224_all_days_birrpoutput'

outdir = 'qel_collected_L224_all_days_birrpoutput'
#outdir = 'testout'

station = 'L224'

plot_component_dict={'0111':'ne'}

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

    donefiles = os.listdir('.')
    donefiles = [i for i in donefiles if i.lower().endswith('.j.done')]
    for df in donefiles:
        os.rename(df,df.replace('.done',''))

    lo_old_coh_files = os.listdir('.')
    lo_old_coh_files = [i for i in lo_old_coh_files if i.lower().endswith('.coh')]
    for i in lo_old_coh_files:
        os.remove(i)

    if 1:
        fullday = date.split('-')[0]
        day_try = int(float(fullday))
        day = int(float(fullday[-2:]))
        month_num = int(float(fullday[-4:-2]))
        #month_num = {'jan':1,'feb':2,'mar':3,'apr':4,'may':5,'jun':6,
        #                'jul':7,'aug':8,'sep':9,'oct':10,'nov':11,'dec':12,}[month]
    # except:
    #     continue
        
    try:
        outfn,outfn_coh = qel2edi.convert2edi(station,'.',survey_configfile,instr_resp,string2strip=['_before','_23Feb'], datestring=fullday)

    except:
        print 'no information found in folder {0}'.format(op.abspath(os.curdir))
        pass
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
        plot_component = 'ne'
        if fullday[-4:] in plot_component_dict:
            plot_component = plot_component_dict[fullday[-4:]]

        plotfn = smplplt.plotedi(outfn,saveplot=True,component=plot_component)
        shutil.copy(op.basename(plotfn),outdir_plots)
        print 'copied res/phase plot %s'%(plotfn)
    except:
        pass

    try:
        plotfncoh = smplpltCOH.plotcoh(outfn_coh,saveplot=True)
        shutil.copy(op.basename(plotfncoh),outdir_plots)
        print 'copied coherence plot %s'%(plotfncoh)
    except:
        pass


    donefiles = os.listdir('.')
    donefiles = [i for i in donefiles if i.lower().endswith('.j.done')]
    for df in donefiles:
        os.rename(df,df.replace('.done',''))
    
    #pdb.set_trace()

    os.chdir(indir)

os.chdir(basedir)







