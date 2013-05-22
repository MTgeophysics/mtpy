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
import matplotlib.cm as cm

class WS3DModelManipulator(object):
    """
    will plot a model from wsinv3d or init file so the user can manipulate the 
    resistivity values relatively easily.  At the moment only plotted
    in map view.
    
    
    """

    def __init__(self, model_fn=None, init_fn=None, data_fn=None,
                 res_lst=None, mapscale='km', plot_yn='y', xlimits=None, 
                 ylimits=None, cbdict={}):
        
        self.model_fn = model_fn
        self.init_fn = init_fn
        self.data_fn = data_fn
        
        #station locations in relative coordinates read from data file
        self.station_x = None
        self.station_y = None
        
        #make a default resistivity list to change values
        if res_lst is None:
            self.res_lst = np.array([.3, 1, 10, 50, 100, 500, 1000, 5000],
                                   dtype=np.float)
        
        else:
            self.res_lst = res_lst        
        
        self.read_file()
   
        #make a dictionary of values to write to file.
        self.res_dict = dict([(res, ii) 
                              for ii,res in enumerate(self.res_lst,1)])
        
        #--> set map scale
        self.mapscale = mapscale
        self.res_value = self.res_lst[0]
        
        #--> set map limits
        self.xlimits = xlimits
        self.ylimits = ylimits
        
        self.cb_dict = cbdict

        self.font_size = 7
        self.dpi = 300
        self.fignum = 1
        self.figsize = [6,6]
        self.cmap = cm.jet_r
        self.depth_index = 0
        
        self.fdict = {'size':self.font_size+2, 'weight':'bold'}
    
        self.cblabeldict = {-5:'$10^{-5}$',
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
        

        
        #plot on initialization
        self.plot_yn = plot_yn
        if self.plot_yn=='y':
            self.plot()
    
    #---read files-------------------------------------------------------------    
    def read_file(self):
        """
        reads in initial file or model file and set attributes:
            -resmodel
            -northrid
            -eastrid
            -zgrid
            -res_lst if initial file
            
        """

        if self.model_fn is not None and self.init_fn is None:
            mtuple = ws.readModelFile(self.model_fn)
            self.north = mtuple[0]
            self.east = mtuple[1]
            self.zg = mtuple[2]
            self.res = mtuple[3]
            self.north_nodes = mtuple[5]
            self.east_nodes = mtuple[6]
            self.z_nodes = mtuple[7]
            
            self.convert_res_to_model()
            
        elif self.init_fn is not None and self.model_fn is None:
            mtuple = ws.readInit3D(self.init_fn)
            self.north = mtuple[0]
            self.east = mtuple[1]
            self.zg = mtuple[2]
            self.res = mtuple[5]
            self.res_lst = mtuple[3]
            self.north_nodes = mtuple[6]
            self.east_nodes = mtuple[7]
            self.z_nodes = mtuple[8]
            
            #need to convert index values to resistivity values
            rdict = dict([(ii,res) for ii,res in enumerate(self.res_lst,1)])
            
            for ii in range(len(self.res_lst)):
                self.res[np.where(self.res==ii+1)] = rdict[ii+1]
                
        elif self.init_fn is None and self.model_fn is None:
            print 'Need to input either an initial file or model file to plot'
        else:
            print 'Input just initial file or model file not both.'
            
        if self.data_fn is not None:
            dtuple = ws.readDataFile(self.data_fn)
            self.station_x = dtuple[3]
            self.station_y = dtuple[4]
            
        #make a copy of original in case there are unwanted changes
        self.res_copy = self.res.copy()
            
            
            
    #---plot model-------------------------------------------------------------    
    def plot(self):
        """
        plots the model with:
            -a radio dial for depth slice 
            -radio dial for resistivity value
            
        """
        
        self.cmin = np.floor(np.log10(min(self.res_lst)))
        self.cmax = np.ceil(np.log10(max(self.res_lst)))
        
        #-->Plot properties
        plt.rcParams['font.size'] = self.font_size
    

        #--> scale the map coordinates
        if self.mapscale=='km':
            dscale = 1000.
        if self.mapscale=='m':
            dscale = 1.
            
        self.dscale = dscale

        
        #make a mesh grid for plotting
        self.northgrid, self.eastgrid = np.meshgrid(self.north/dscale, 
                                                    self.east/dscale)
        
        self.fig = plt.figure(self.fignum, figsize=self.figsize, dpi=self.dpi)
        self.ax1 = self.fig.add_subplot(1, 1, 1, aspect='equal')
        
        plot_res = np.log10(np.rot90(self.res[:,:,self.depth_index],3))
        
        self.mesh_plot = self.ax1.pcolormesh(self.eastgrid, self.northgrid, 
                                             plot_res,
                                             cmap=self.cmap,
                                             vmin=self.cmin,
                                             vmax=self.cmax)
                                             
        #on plus or minus change depth slice
        self.cid_depth = \
                    self.mesh_plot.figure.canvas.mpl_connect('key_press_event',
                                                        self._on_key_callback)
                                    
                       
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
            self.ax1.set_xlim(xmin=self.north.min()/dscale, 
                              xmax=self.north.max()/dscale)
        
        if self.ylimits is not None:
            self.ax1.set_ylim(self.ylimits)
        else:
            self.ax1.set_ylim(ymin=self.east.min()/dscale,
                              ymax=self.east.max()/dscale)
            
        self.ax1.xaxis.set_minor_locator(MultipleLocator(100*1./dscale))
        self.ax1.yaxis.set_minor_locator(MultipleLocator(100*1./dscale))
        
        self.ax1.set_ylabel('Northing ('+self.mapscale+')',
                            fontdict=self.fdict)
        self.ax1.set_xlabel('Easting ('+self.mapscale+')',
                            fontdict=self.fdict)
        
        depth_title = self.zg[self.depth_index]/self.dscale
                                                        
        self.ax1.set_title('Depth = {:.3f} '.format(depth_title)+\
                           '('+self.mapscale+')',
                           fontdict=self.fdict)
        
        #plot the grid if desired              
        for xx in self.east:
            self.ax1.plot([self.north.min()/dscale, self.north.max()/dscale],
                           [xx/dscale, xx/dscale],
                           lw=.25,
                           color='k')

        for yy in self.north:
            self.ax1.plot([yy/dscale, yy/dscale], 
                          [self.east.min()/dscale, self.east.max()/dscale],
                           lw=.25,
                           color='k')
        
        #plot the colorbar
        self.ax2 = mcb.make_axes(self.ax1, orientation='vertical', shrink=.5)
        seg_cmap = cmap_discretize(self.cmap, len(self.res_lst))
        self.cb = mcb.ColorbarBase(self.ax2[0],cmap=seg_cmap,
                                   norm=colors.Normalize(vmin=self.cmin,
                                                         vmax=self.cmax))
                                                         
        
#        #make bounds so that the middle is white
#        bounds = np.arange(self.cmin, self.cmax, len(self.res_lst))
#        
#        #normalize the colors
#        norms = colors.BoundaryNorm(bounds, seg_cmap.N)
#        
#        #make the colorbar
#        self.cb = mcb.ColorbarBase(self.ax2,
#                                   cmap=seg_cmap,
#                                   norm=norms,
#                                   orientation='vertical',
#                                   ticks=bounds[1:-1])
                            
        self.cb.set_label('Resistivity ($\Omega \cdot$m)',
                     fontdict={'size':self.font_size})
        self.cb.set_ticks(np.arange(self.cmin, self.cmax+1))
        self.cb.set_ticklabels([self.cblabeldict[cc] 
                            for cc in np.arange(self.cmin, self.cmax+1)])
                            
        #make a resistivity radio button
        resrb = self.fig.add_axes([.85,.1,.1,.15])
        reslabels = ['{0:.4g}'.format(res) for res in self.res_lst]
        self.radio_res = widgets.RadioButtons(resrb, reslabels,active=0)
        
        #make a rectangular selector
        self.rect_selector = widgets.RectangleSelector(self.ax1, 
                                                       self.rect_onselect,
                                                       drawtype='box',
                                                       useblit=True)

        
        plt.show()
        
        #needs to go after show()
        self.radio_res.on_clicked(self.set_res_value)


    def redraw_plot(self):
        """
            redraws the plot
        """
        
        self.ax1.cla()
        
        plot_res = np.log10(np.rot90(self.res[:,:,self.depth_index],3))
        
        self.mesh_plot = self.ax1.pcolormesh(self.eastgrid, self.northgrid, 
                                             plot_res,
                                             cmap=self.cmap,
                                             vmin=self.cmin,
                                             vmax=self.cmax)
                                             
         #plot the stations
        if self.station_x is not None:
            for ee,nn in zip(self.station_x, self.station_y):
                self.ax1.text(ee/self.dscale, nn/self.dscale,
                              '*',
                              verticalalignment='center',
                              horizontalalignment='center',
                              fontdict={'size':self.font_size-2,
                                        'weight':'bold'})

        #set axis properties
        if self.xlimits is not None:
            self.ax1.set_xlim(self.xlimits)
        else:
            self.ax1.set_xlim(xmin=self.north.min()/self.dscale, 
                              xmax=self.north.max()/self.dscale)
        
        if self.ylimits is not None:
            self.ax1.set_ylim(self.ylimits)
        else:
            self.ax1.set_ylim(ymin=self.east.min()/self.dscale,
                              ymax=self.east.max()/self.dscale)
            
        self.ax1.xaxis.set_minor_locator(MultipleLocator(100*1./self.dscale))
        self.ax1.yaxis.set_minor_locator(MultipleLocator(100*1./self.dscale))
        
        self.ax1.set_ylabel('Northing ('+self.mapscale+')',
                            fontdict=self.fdict)
        self.ax1.set_xlabel('Easting ('+self.mapscale+')',
                            fontdict=self.fdict)
        
        depth_title = self.zg[self.depth_index]/self.dscale
                                                        
        self.ax1.set_title('Depth = {:.3f} '.format(depth_title)+\
                           '('+self.mapscale+')',
                           fontdict=self.fdict)
                     
        #plot the grid if desired              
        for xx in self.east:
            self.ax1.plot([self.north.min()/self.dscale, 
                           self.north.max()/self.dscale],
                           [xx/self.dscale, xx/self.dscale],
                           lw=.25,
                           color='k')

        for yy in self.north:
            self.ax1.plot([yy/self.dscale, yy/self.dscale], 
                          [self.east.min()/self.dscale, 
                           self.east.max()/self.dscale],
                           lw=.25,
                           color='k')
        
        #be sure to redraw the canvas                  
        self.fig.canvas.draw()
        
    def set_res_value(self, label):
        self.res_value = float(label)
        print 'set resistivity to ', label
        print self.res_value
        
        
    def _on_key_callback(self,event):
        """
        on pressing a key do something
        
        """
        
        self.event_change_depth = event

        if self.event_change_depth.key=='=':
            self.depth_index += 1
            
            if self.depth_index>len(self.zg)-1:
                self.depth_index = len(self.zg)-1
                print 'already at deepest depth'
                
            print 'Plotting Depth {0:.3f}'.format(self.zg[self.depth_index]/\
                    self.dscale)+'('+self.mapscale+')'
            
            self.redraw_plot()

        elif self.event_change_depth.key=='-':
            self.depth_index -= 1
            
            if self.depth_index<0:
                self.depth_index = 0
                
            print 'Plotting Depth {0:.3f} '.format(self.zg[self.depth_index]/\
                    self.dscale)+'('+self.mapscale+')'
            
            self.redraw_plot()

        elif self.event_change_depth.key == 'q':
            self.event_change_depth.canvas.mpl_disconnect(self.cid_depth)
            plt.close(self.event_change_depth.canvas.figure)
            
        #copy the layer above
        elif self.event_change_depth.key == 'a':
            try:
                if self.depth_index == 0:
                    print 'No layers above'
                else:
                    self.res[:,:,self.depth_index] = \
                                              self.res[:,:,self.depth_index-1]
            except IndexError:
                print 'No layers above'
                
            self.redraw_plot()
        
        #copy the layer below
        elif self.event_change_depth.key == 'b':
            try:
                self.res[:,:,self.depth_index] = \
                                              self.res[:,:,self.depth_index+1]
            except IndexError:
                print 'No more layers below'
                
            self.redraw_plot() 
            
        #undo
        elif self.event_change_depth.key == 'u':
            self.res[self.ni0:self.ni1, self.ei0:self.ei1, self.depth_index] =\
            self.res_copy[self.ni0:self.ni1,self.ei0:self.ei1,self.depth_index]
            
            self.redraw_plot()
            
           
    def rect_onselect(self, eclick, erelease):
        """
        on selecting a rectangle change the colors to the resistivity values
        """
        x1, y1 = eclick.xdata, eclick.ydata
        x2, y2 = erelease.xdata, erelease.ydata
        
        
        ei0, ei1 = self._get_east_index(x1, x2)
        ni0, ni1 = self._get_north_index(y1, y2)
        
        #reset values of resistivity
        
        self.res[ni0:ni1, ei0:ei1, self.depth_index] = self.res_value
        
        self.redraw_plot()
        
        self.ei0 = ei0
        self.ei1 = ei1
        self.ni0 = ni0
        self.ni1 = ni1
            
        
        
    def _get_east_index(self, x1, x2):
        """
        get the index value of the points to be changed
        
        """
        
        x1 *= 1.
        x2 *= 1.
        
        if x1 < x2:
            p1 = 0
        elif x2 < x1:
            p1 = 0
        elif x1 == x2:
            p1 = 0
        
        nx = len(self.east)
        xreturn=[]
        for xx in [x1, x2]:
            for ii in range(nx-1):
                if xx == self.east[ii]/self.dscale:
                    xreturn.append(ii)
                    break
                
                elif xx >= self.east[ii]/self.dscale and \
                        xx <= self.east[ii+1]/self.dscale:
                    xreturn.append(ii+p1)
                    break
                
                elif xx > self.east[-1]/self.dscale:
                    xreturn.append(nx-1)
                    break
                
                elif xx<self.east[0]/self.dscale:
                    xreturn.append(0)
                    break
                
        if xreturn[0] > xreturn[1]:
            return xreturn[1], xreturn[0]
            
        elif xreturn[0] < xreturn[1]:
            return xreturn[0], xreturn[1]
            
        elif xreturn[0]==xreturn[1]:
            return xreturn[0], xreturn[0]+1
                
    def _get_north_index(self, y1, y2):
        """
        get the index value of the points to be changed in north direction
        
        need to flip the index because the plot is flipped
        
        """
        
        y1 *= -1
        y2 *= -1
        
        if y1 < y2:
            p1 = 1
        elif y2 < y1:
            p1 = 1
        elif y1 == y2:
            p1 = 1
            
        ny = len(self.north)
        yreturn=[]
        for yy in [y1, y2]:
            for ii in range(ny-1):
                if yy==self.north[ii]/self.dscale:
                    yreturn.append(ii)
                    break
                
                elif yy>=self.north[ii]/self.dscale and \
                        yy<=self.north[ii+1]/self.dscale:
                    yreturn.append(ii+p1)
                    break
                
                elif yy>self.north[-1]/self.dscale:
                    yreturn.append(ny-1)
                    break
                
                elif yy<self.north[0]/self.dscale:
                    yreturn.append(0)
                    break
                    
        if yreturn[0] > yreturn[1]:
            return yreturn[1], yreturn[0]
        
        elif yreturn[0] < yreturn[1]:
            return yreturn[0], yreturn[1]
            
        elif yreturn[0]==yreturn[1]:
            return yreturn[0], yreturn[0]+1
            
    def convert_model_to_int(self):
        """
        convert the resistivity model that is in ohm-m to integer values
        corresponding to res_lst
        
        """
 
        self.res_model = np.ones_like(self.res)
        
        for key in self.res_dict.keys():
            self.res_model[np.where(self.res==key)] = self.res_dict[key]
            
    def convert_res_to_model(self):
        """
        converts an output model into an array of segmented valued according
        to res_lst.        
        
        """
        
        #make values in model resistivity array a value in res_lst
        resm = np.ones_like(self.res)
        resm[np.where(self.res<self.res_lst[0])] = \
                                            self.res_dict[self.res_lst[0]]
        resm[np.where(self.res)>self.res_lst[-1]] = \
                                            self.res_dict[self.res_lst[-1]]
        
        for zz in range(self.res.shape[2]):
            for yy in range(self.res.shape[1]):
                for xx in range(self.res.shape[0]):
                    for rr in range(len(self.res_lst)-1):
                        if self.res[xx,yy,zz]>=self.res_lst[rr] and \
                            self.res[xx,yy,zz]<=self.res_lst[rr+1]:
                            resm[xx,yy,zz] = self.res_dict[self.res_lst[rr]]
                            break
                        elif self.res[xx,yy,zz]<=self.res_lst[0]:
                            resm[xx,yy,zz] = self.res_dict[self.res_lst[0]]
                            break
                        elif self.res[xx,yy,zz]>=self.res_lst[-1]:
                            resm[xx,yy,zz] = self.res_dict[self.res_lst[-1]]
                            break
    
        self.res = resm
            
        
    def write_init_file(self, savepath, north_nodes=None, east_nodes=None,
                        z_nodes=None, title='Initial Model for wsinv3d'):
        """
        write an initial file for wsinv3d from the model created.
        """
        
        self.convert_model_to_int()
        
        try:
            init_new = ws.writeInit3DFile(self.north_nodes, 
                                          self.east_nodes,
                                          self.z_nodes, 
                                          savepath, 
                                          reslst=self.res_lst,
                                          title=title,
                                          resmodel=self.res_model)
            return init_new
            
        except AttributeError:
            if north_nodes is not None:
                init_new = ws.writeInit3DFile(north_nodes, 
                                              east_nodes,
                                              z_nodes, 
                                              savepath, 
                                              reslst=self.res_lst,
                                              title=title,
                                              resmodel=self.res_model)
                return init_new
            else:
                print 'Need to input the starting grid'
                                              
                

def cmap_discretize(cmap, N):
    """Return a discrete colormap from the continuous colormap cmap.
      
         cmap: colormap instance, eg. cm.jet. 
         N: number of colors.
     
     Example
         x = resize(arange(100), (5,100))
         djet = cmap_discretize(cm.jet, 5)
         imshow(x, cmap=djet)
    """

    colors_i = np.concatenate((np.linspace(0, 1., N), (0.,0.,0.,0.)))
    colors_rgba = cmap(colors_i)
    indices = np.linspace(0, 1., N+1)
    cdict = {}
    for ki,key in enumerate(('red','green','blue')):
        cdict[key] = [ (indices[i], colors_rgba[i-1,ki], colors_rgba[i,ki]) for i in xrange(N+1) ]
    # Return colormap object.
    return colors.LinearSegmentedColormap(cmap.name + "_%d"%N, cdict, 1024)

        
        
            
        
        