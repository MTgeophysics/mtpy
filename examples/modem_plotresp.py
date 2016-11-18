# -*- coding: utf-8 -*-
"""
Created on Mon Oct 17 10:31:59 2016

@author: u64125
"""

import mtpy.modeling.modem_new as mtmn
import os.path as op
import os
import matplotlib.pyplot as plt
import numpy as np

workdir = r'C:\Git\mtpydev\examples\data'
modeldir = op.join(workdir,'ModEM_files')


read_data = False
respfn = op.join(modeldir,max([ff for ff in os.listdir(modeldir) if (ff.endswith('.dat') and (ff.startswith('Modular')))]))
datafn = op.join(modeldir,'ModEM_Data.dat')
if read_data:
    doo = mtmn.Data()
    doo.read_data_file(op.join(modeldir,'ModEM_Data.dat'))
    roo = mtmn.Data()
    roo.read_data_file(op.join(modeldir,respfn))

#mtmn.PlotResponse(data_fn = datafn,
#                  resp_fn = respfn, 
#                  plot_type=['VIC028'],
##                  res_limits=(1,100),
##                  phase_limits=(0,90),
#                  plot_z=False)
def resistivity(z,period,zerr=None):
    res = 0.2*np.abs((z.transpose(1,2,0)**2.*period).transpose(2,0,1))
    if zerr is None:
        reserr = np.zeros_like(res)
    else:
        reserr = res*2.*zerr/np.abs(z)
    return res,reserr

def phase(z,zerr=None):
    phs = np.degrees(np.arctan(np.imag(z)/np.real(z)))
    if zerr is None:
        phserr = np.zeros_like(phs)
    else:
        phserr = np.rad2deg(np.arctan(np.abs(zerr)/np.abs(z)))
    return phs,phserr
    
    
def plotres(period, resistivity_data=None, resistivity_error=None, resistivity_mod=None, c='b', label=''):
    if resistivity_mod is not None:
        plt.plot(period,resistivity_mod,c+'-',label=label)
    if resistivity_data is not None:
        if resistivity_error is None:
            resistivity_error = np.zeros_like(resistivity_data)
        plt.errorbar(period,resistivity_data,yerr=resistivity_error,color=c,label=label,fmt='.')
    plt.xscale('log')
    plt.yscale('log')    
    
def plotphs(period, phase_data=None, phase_error = None, phase_mod=None, c='b', label=''):
    if phase_mod is not None:
        plt.plot(period,phase_mod,c+'-',label=label)
    if phase_data is not None:
        if phase_error is None:
            phase_error=np.zeros_like(phase_data)
        plt.errorbar(period,phase_data,yerr=phase_error,color=c,label=label,fmt='.')
    
    plt.xscale('log')

def plotresponse(station,data_object=None,response_object=None,reslim_diags=None,
                 phaselim_diags=None,reslim_offdiags=None,phaselim_offdiags=None):
    
    if data_object is not None:
        data_array = data_object.data_array
        period = data_object.period_list
    if response_object is not None:
        response_array = response_object.data_array
        period = data_object.period_list        
    
    if ((response_object is None) and (data_object is None)):
        print "No data to plot"
        return
    else:
        plt.figure(figsize=(8,4))
        ss = np.where(data_array['station'] == station)[0][0]
        resdata,reserr = resistivity(data_array['z'][ss],period,zerr=data_array['z_err'][ss])
        phsdata,phserr = phase(data_array['z'][ss],zerr=data_array['z_err'][ss])
        resmod = resistivity(response_array['z'][ss],period)[0]
        phsmod = phase(response_array['z'][ss])[0]
        sp = 1
        for i in range(2):
            for j in range(2):
                c='br'[i]
                plt.subplot(2,4,sp)
                plotres(period,resistivity_data=resdata[:,i,j],
                        resistivity_error=reserr[:,i,j],
                        resistivity_mod=resmod[:,i,j],
                        c=c)
                if i == j:
                    if reslim_diags is not None:
                        plt.ylim(*reslim_diags)
                else:
                    if reslim_offdiags is not None:
                        plt.ylim(*reslim_offdiags)                    
                plt.gca().xaxis.set_visible(False)
                plt.subplot(2,4,sp+4)
                plotphs(period,phase_data=phsdata[:,i,j],
                        phase_error=phserr[:,i,j],
                        phase_mod=phsmod[:,i,j],
                        c=c)
                if i == j:
                    if phaselim_diags is not None:
                        plt.ylim(*phaselim_diags)
                else:
                    if phaselim_offdiags is not None:
                        plt.ylim(*phaselim_offdiags)  
                sp += 1
        plt.subplots_adjust(hspace=0.05)

plotresponse('GB01',doo,roo,reslim_offdiags=(100,10000),reslim_diags=(1,1000),
             phaselim_offdiags=(0,90),phaselim_diags=(-180,180))