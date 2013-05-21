# -*- coding: utf-8 -*-
"""
Created on Mon May 20 16:00:13 2013

@author: jpeacock-pr
"""

import mtpy.modeling.ws3dtools as ws
import matplotlib.pyplot as plt
import matplotlib.widgets as widgets
import numpy as np
from matplotlib.ticker import MultipleLocator
import matplotlib.colorbar as mcb
import matplotlib.colors as colors


class WS3DModelManipulator(object):
    """
    will plot a model from wsinv3d or init file so the user can manipulate the 
    resistivity values relatively easily.  At the moment only plotted
    in map view.
    
    
    """

    def __init__(self, model_fn=None, init_fn=None, data_fn=None,
                 reslst=None, mapscale='km', plot_yn='y', xlimits=None, 
                 ylimits=None, cbdict={}):
        
        self.model_fn = model_fn
        self.init_fn = init_fn
        self.data_fn = data_fn
        
        #make a default resistivity list to change values
        if reslst is None:
            self.reslst = np.array([.3, 1, 10, 50, 100, 500, 1000, 5000],
                                   dtype=np.float)
        
        else:
            self.reslst = reslst
            
        #make a dictionary of values to write to file.
        self.resdict = dict([(res, ii) for ii,res in enumerate(self.reslst,1)])
        
        #--> set map scale
        self.mapscale = mapscale
        
        #--> set map limits
        self.xlimits = xlimits
        self.ylimits = ylimits
        
        self.cb_dict = cbdict

        self.font_size = 7
        self.dpi = 300
        self.fignum = 1
        self.figsize = [6,6]
        self.cmap = 'jet_r'
        self.depth_index = 0
        
        #station locations in relative coordinates read from data file
        self.station_x = None
        self.station_y = None
        
        #plot on initialization
        self.plot_yn = plot_yn
        if self.plot_yn=='y':
            self.plot()
    
    #---read files-------------------------------------------------------------    
    def read_file(self):
        """
        reads in initial file or model file and set attributes:
            -resmodel
            -xgrid
            -ygrid
            -zgrid
            -reslst if initial file
            
        """

        if self.model_fn is not None and self.init_fn is None:
            mtuple = ws.readModelFile(self.model_fn)
            self.xg = mtuple[0]
            self.yg = mtuple[1]
            self.zg = mtuple[2]
            self.res = mtuple[3]
            
        elif self.init_fn is not None and self.model_fn is None:
            mtuple = ws.readInit3D(self.init_fn)
            self.xg = mtuple[0]
            self.yg = mtuple[1]
            self.zg = mtuple[2]
            self.res = mtuple[5]
            self.reslst = mtuple[3]
            
            #need to convert index values to resistivity values
            rdict = dict([(ii,res) for ii,res in enumerate(self.reslst,1)])
            
            for ii in range(len(self.reslst)):
                self.res[np.where(self.res==ii+1)] = rdict[ii+1]
                
        elif self.init_fn is None and self.model_fn is None:
            print 'Need to input either an initial file or model file to plot'
        else:
            print 'Input just initial file or model file not both.'
            
        if self.data_fn is not None:
            dtuple = ws.readDataFile(self.data_fn)
            self.station_x = dtuple[3]
            self.station_y = dtuple[4]
            
            
            
    #---plot model-------------------------------------------------------------    
    def plot(self):
        """
        plots the model with:
            -a radio dial for depth slice 
            -radio dial for resistivity value
            
        """
        
        self.read_file()
        
        cmin = np.floor(np.log10(min(self.reslst)))
        cmax = np.ceil(np.log10(max(self.reslst)))
        print cmin, cmax
        
        #-->Plot properties
        plt.rcParams['font.size'] = self.font_size
        
        fdict = {'size':self.font_size+2, 'weight':'bold'}
    
        cblabeldict = {-5:'$10^{-5}$',
                       -4:'$10^{-4}$',
                       -3:'$10^{-3}$',
                       -2:'$10^{-2}$',
                       -1:'$10^{-1}$',
                        0:'$10^{0}$',
                        1:'$10^{1}$',
                        2:'$10^{2}$',
                        3:'$10^{3}$',
                        4:'$10^{4}$',
                        5:'$10^{5}$',
                        6:'$10^{6}$',
                        7:'$10^{7}$',
                        8:'$10^{8}$'}
        
        
        #--> scale the map coordinates
        if self.mapscale=='km':
            dscale = 1000.
        if self.mapscale=='m':
            dscale = 1.
            
        self.dscale = dscale

        
        #make a mesh grid for plotting
        xgrid, ygrid = np.meshgrid(self.xg/dscale, self.yg/dscale)
        
        self.fig = plt.figure(self.fignum, figsize=self.figsize, dpi=self.dpi)
        self.ax1 = self.fig.add_subplot(1, 1, 1, aspect='equal')
        
        self.ax1.pcolormesh(ygrid, xgrid, 
                            np.log10(np.rot90(self.res[:,:,self.depth_index],3)),
                            cmap=self.cmap,
                            vmin=cmin,
                            vmax=cmax)
                       
        #plot the stations
        if self.station_x is not None:
            for ee,nn in zip(self.station_x, self.station_y):
                self.ax1.text(ee/dscale, nn/dscale,
                              '*',
                              verticalalignment='center',
                              horizontalalignment='center',
                              fontdict={'size':self.font_size-2,
                                        'weight':'bold'})

        #set axis properties
        if self.xlimits is not None:
            self.ax1.set_xlim(self.xlimits)
        else:
            self.ax1.set_xlim(xmin=self.xg.min()/dscale, 
                              xmax=self.xg.max()/dscale)
        
        if self.ylimits is not None:
            self.ax1.set_ylim(self.ylimits)
        else:
            self.ax1.set_ylim(ymin=self.yg.min()/dscale,
                              ymax=self.yg.max()/dscale)
            
        self.ax1.xaxis.set_minor_locator(MultipleLocator(100*1./dscale))
        self.ax1.yaxis.set_minor_locator(MultipleLocator(100*1./dscale))
        
        self.ax1.set_ylabel('Northing ('+self.mapscale+')',fontdict=fdict)
        self.ax1.set_xlabel('Easting ('+self.mapscale+')',fontdict=fdict)
        
        self.ax1.set_title('Depth = {:.3f} '.format(
                                            self.zg[self.depth_index]/dscale)+\
                                            '('+self.mapscale+')',
                           fontdict=fdict)
        
        #plot the grid if desired              
        for xx in self.yg:
            self.ax1.plot([self.xg.min()/dscale, self.xg.max()/dscale],
                           [xx/dscale, xx/dscale],
                           lw=.25,
                           color='k')
        for yy in self.xg:
            self.ax1.plot([yy/dscale, yy/dscale], 
                          [self.yg.min()/dscale, self.yg.max()/dscale],
                           lw=.25,
                           color='k')
        
        #plot the colorbar
        self.ax2 = mcb.make_axes(self.ax1, orientation='vertical', shrink=.5)
        
        self.cb = mcb.ColorbarBase(self.ax2[0],cmap=self.cmap,
                                   norm=colors.Normalize(vmin=cmin,
                                                         vmax=cmax))
                            
        self.cb.set_label('Resistivity ($\Omega \cdot$m)',
                     fontdict={'size':self.font_size})
        self.cb.set_ticks(np.arange(cmin, cmax+1))
        self.cb.set_ticklabels([cblabeldict[cc] 
                            for cc in np.arange(cmin,cmax+1)])
                            
        #make a resistivity radio button
        resrb = self.fig.add_axes([.85,.1,.1,.15])
        reslabels = ['{0:.4g}'.format(res) for res in self.reslst]
        self.radio_res = widgets.RadioButtons(resrb, reslabels,active=0)
        


        
        plt.show()
        
        #needs to go after show()
        self.radio_res.on_clicked(self.set_res_value)
        
        #on plus or minus change depth slice
        self.cid_depth = \
                    self.ax1.figure.canvas.mpl_connect('button_press_event',
                                                        self._change_depth)
        


    def redraw_plot(self):
        """
            redraws the plot
        """
        
        self.fig.clear()
        self.plot()
        
    def set_res_value(self, label):
        print 'set resistivity to ',label
        
        
    def _change_depth(self,event):
        self.event_change_depth = event
        if self.event_change_depth.key=='=':
            self.depth_index += 1
            if self.depth_index<0:
                self.depth_index = 0
                print 'Already at shallowest depth try "=" button'
            print 'Plotting Depth {0:.1f}'.format(self.zg[self.depth_index/\
                    self.dscale])+'('+self.mapscale+')'
            self.redraw_plot()
        
            
        elif self.event_change_depth.key=='-':
            self.depth_index -= 1
            if self.depth_index>len(self.zg)-1:
                self.depth_index = len(self.zg)-1
                print 'already at deepest depth'
            print 'Plotting Depth {0:.1f} '.format(self.zg[self.depth_index/\
                    self.dscale])+'('+self.mapscale+')'
            self.redraw_plot()
            
        elif self.event_change_depth.key=='q':
            self.event_change_depth.canvas.mpl_disconnect(self.cid_depth)
            plt.close(self.event_change_depth.canvas.figure)
        
        
 
        
                        

        
        
            
        
        