# -*- coding: utf-8 -*-
"""
Created on Fri Aug 08 11:16:46 2014

@author: a1655681
"""

import matplotlib.pyplot as plt
import numpy as np

class PlotResponses():
    """
    plot responses and data from any 2d model, requires two dictionaries
    containing period, resistivity, phase arrays with the following dimensions:
    period: (nperiod)
    resistivity: (num_stations, num_periods, 2, 2)
    phase: (num_stations, num_periods, 2, 2)
    
    """
    def __init__(self, measured, modelled, **input_parameters):
                
        
        self.period_mod = modelled['period']
        self.period_meas = measured['period']
        self.resistivity_mod = modelled['resistivity']
        self.resistivity_meas = measured['resistivity']
        self.phase_mod = modelled['phase']
        self.phase_meas = measured['phase']
        self.subplot_rows = 1 # number of rows to plot
        self.plot_tf = True # whether or not to plot on creation of an object
        self.symbol_dict = {'measured':'--','modelled':'-','colors':['b','r']}
        self.xlim = [1e-3,1e3]
        self.ylim = {'resistivity':[[0.01,1],[1,100]],'phase':[[0,90],[0,90]]}
    
        for key in input_parameters.keys():
            setattr(self,key,input_parameters[key])   

        if self.plot_tf:
            self.plot()

    def __set_axis_params(self,parameter,component=1,labels='xy'):
        """
        component = diagonals or offdiagonals, 1 for diags, 0 for offdiags
        """
        ax = plt.gca()
        plt.xscale('log')
        
        if parameter == 'resistivity':
            plt.yscale('log')
        
        plt.xlim(self.xlim)
        plt.ylim(self.ylim[parameter][component])
        
        if 'x' not in labels:
            ax.xaxis.set_visible(False)
        if 'y' not in labels:
            ax.yaxis.set_visible(False)
    
#    def prepare_dictionaries(self):


    def plot(self):
        """
        """

        nsp = len(self.resistivity_mod)
        
        # rows, columns in subplots
        # multiply rows by 2 as we are plotting res and phase
        spr = self.subplot_rows*2
        # number of columns
        spc = 2*int(np.ceil(float(nsp)/float(self.subplot_rows)))
        # starting subplot number
        row = -1        

        for i in range(nsp):
            if (2*i)%spc == 0:
                row += 1
            spn = 2*i + 1 + spc*row
            for parameter,meas,mod in [['resistivity',self.resistivity_meas[i],self.resistivity_mod[i]],
                                       ['phase',self.phase_meas[i]%180,self.phase_mod[i]%180]]:
                print i
                for ii in range(2):
                    for jj in range(2):
                        labels = ''
                        if row == int(spr/2.) - 1:
                            if parameter == 'phase':
                                labels += 'x' 
                        if (2*i)%spc == 0:
                            if (ii != jj):
                                labels += 'y'
                        if ii == jj:
                            cpt = 0
                            plt.subplot(spr,spc,spn+1)
                        else:
                            cpt = 1
                            plt.subplot(spr,spc,spn)
                        if ii == 0:
                            c = self.symbol_dict['colors'][0]
                        else:
                            c = self.symbol_dict['colors'][1]
                        plt.plot(self.period_meas[i],meas[:,ii,jj],
                                 self.symbol_dict['measured'],
                                 color = c)
                        plt.plot(self.period_mod[i],mod[:,ii,jj],
                                 self.symbol_dict['modelled'],
                                 color = c)
                        self.__set_axis_params(parameter,
                                               component=cpt,
                                               labels=labels)

                spn += spc
        plt.subplots_adjust(hspace=0.02,wspace=0.02)
            

        