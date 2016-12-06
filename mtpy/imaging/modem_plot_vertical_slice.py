# -*- coding: utf-8 -*-
"""
Visualize Vertical Slices of the ModEM Model

Created on Tue 2016-12-05
by
fei.zhang@ga.gov.au
"""

import os.path as op

import matplotlib.pyplot as plt
import numpy as np

from mtpy.modeling.modem_new import Model
from mtpy.modeling.modem_new import Data

class ModemPlotVerticalSlice():
    
    def __init__(self,filedat, filerho, plot_orient='ew', **kwargs):
        """Constructor
        :param filedat: path2file.dat
        :param filerho: path2file.rho
        :param  plot_orient: plot orientation ['ew','ns', 'z']
        """
        self.datfile= filedat
        self.rhofile= filerho
        
        # plot orientation ('ns' (north-south),'ew' (east-west) or 'z' (horizontal slice))
        self.plot_orientation=plot_orient

        #plotdir = 'z' #'ns' #'ew'
        # slice location, in local grid coordinates (if it is a z slice, this is slice depth)
        self.slice_location = 10000
        # maximum distance in metres from vertical slice location and station
        self.stationdist = 50000
        # z limits (positive down so order is reversed)
        self.zlim = (1e5, -5e3)
        # colour limits
        self.clim = [0.3, 3.7]

        self.fig_size = kwargs.pop('fig_size', [16, 12])
        self.font_size = kwargs.pop('font_size', 10)
        self.map_scale='km'
        
        # read in the model data
        self._read_data()
        
        return

    def _read_data(self):

        self.datObj = Data()
        self.datObj.read_data_file(data_fn=self.datfile)
        
        self.modObj = Model(model_fn= self.rhofile)
        self.modObj.read_model_file()
        
        return

    def set_plot_orientation(self, orient):
        """set a new plot orientation for plotting
        :param orient: z, ew, ns
        :return:
        """
        self.plot_orientation=orient


    def make_plot(self):
        """ create a plot based on the input data and parameters       
        :return: 
        """

        fdict = {'size': self.font_size + 2, 'weight': 'bold'}

        # get grid centres
        gcz = np.mean([self.modObj.grid_z[:-1], self.modObj.grid_z[1:]], axis=0)
        gceast, gcnorth = [np.mean([arr[:-1], arr[1:]], axis=0) for arr in [self.modObj.grid_east, self.modObj.grid_north]]
        
        # distance from slice to grid centre locations
        if self.plot_orientation == 'ew':
            sdist = np.abs(gcnorth - self.slice_location)
        elif self.plot_orientation == 'ns':
            sdist = np.abs(gceast - self.slice_location)
        elif self.plot_orientation == 'z':
            sdist = np.abs(gcz - self.slice_location)
        
        # find closest slice index to specified location
        sno = np.where(sdist == np.amin(sdist))[0][0]
        
        # get data for plotting
        if self.plot_orientation == 'ew':
            X, Y, res = self.modObj.grid_east, self.modObj.grid_z, np.log10(self.modObj.res_model[sno, :, :].T)
            ss = np.where(np.abs(self.datObj.station_locations['rel_north'] - np.median(gcnorth)) < self.stationdist)[0]
        
            sX, sY = self.datObj.station_locations['rel_east'][ss], self.datObj.station_locations['elev'][ss]
            xlim = (self.modObj.grid_east[self.modObj.pad_east], self.modObj.grid_east[-self.modObj.pad_east - 1])
            ylim = self.zlim
            title = 'East-west slice at {}km north'.format(gcnorth[sno])
        elif self.plot_orientation == 'ns':
            X, Y, res = self.modObj.grid_north, self.modObj.grid_z, np.log10(self.modObj.res_model[:, sno, :].T)
            # indices for selecting stations close to profile
            ss = np.where(np.abs(self.datObj.station_locations['rel_east'] - np.median(gceast)) < self.stationdist)[0]
        
            sX, sY = self.datObj.station_locations['rel_north'][ss], self.datObj.station_locations['elev'][ss]
            xlim = (self.modObj.grid_north[self.modObj.pad_north], self.modObj.grid_north[-self.modObj.pad_north - 1])
            ylim = self.zlim
            title = 'North-south slice at {}km east'.format(gceast[sno])
        elif self.plot_orientation == 'z':
            X, Y, res = self.modObj.grid_east, self.modObj.grid_north, np.log10(self.modObj.res_model[:, :, sno])
            sX, sY = self.datObj.station_locations['rel_east'], self.datObj.station_locations['rel_north']
            xlim = (self.modObj.grid_east[self.modObj.pad_east], self.modObj.grid_east[-self.modObj.pad_east - 1])
            ylim = (self.modObj.grid_north[self.modObj.pad_north], self.modObj.grid_north[-self.modObj.pad_north - 1])
            title = 'Horizontal Slice at Depth {}km'.format(gcz[sno])
        
        # make the plot
        plt.figure(figsize=self.fig_size)
        mesh_plot=plt.pcolormesh(X, Y, res, cmap='bwr_r')
        plt.xlim(*xlim)
        plt.ylim(*ylim)
        
        # plot station locations
        plt.plot(sX, sY, 'o')  #station marker:'kv'
        
        # set title
        plt.title(title)
        
        if self.plot_orientation == 'z':
            plt.gca().set_aspect('equal')

        plt.clim(*self.clim)
        #plt.colorbar()

        # FZ: fix miss-placed colorbar
        from mpl_toolkits.axes_grid1 import make_axes_locatable
        ax = plt.gca()

        # create an axes on the right side of ax. The width of cax will be 5%
        # of ax and the padding between cax and ax will be fixed at 0.05 inch.
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.05)

        mycb = plt.colorbar(mesh_plot, cax=cax, label='Resistivity ($\Omega \cdot$m)')

        ax.set_ylabel('Northing (' + self.map_scale + ')', fontdict=fdict)
        ax.set_xlabel('Easting (' + self.map_scale + ')', fontdict=fdict)

        plt.show()

#########################################################################
# Example code how to use the class and in-situ testing of the class
# Usage:
#   export PYTHONPATH="E:/Githubz/mtpy2"
#   python  mtpy/imaging/modem_plot_vertical_slice.py
#
#   (Future: python this_script.py  path2datfile path2rhofile plot_orient)
#-----------------------------------------------------------------------
if __name__=="__main__":
    # INPUTS #
    # define a workdir for your environ
    workdir = r'V:\Geology\conductivity_modelling'
    workdir = r'E:\Githubz\mtpy2\examples\data\ModEM_files'
    #workdir = r'/Softlab/Githubz/mtpy2/examples/data/ModEM_files'
    # workdir = r'/g/data/ha3/fxz547/Githubz/mtpy2/examples/data/ModEM_files'

    modeldir = op.join(workdir, 'VicSynthetic07')

    datafn = 'ModEM_Data_noise10inv.dat'
    iterfn = 'Modular_MPI_NLCG_019.rho'

    datf=op.join(modeldir, datafn)
    rhof=op.join(modeldir,iterfn)

    myObj=ModemPlotVerticalSlice(datf, rhof, plot_orient='ew')
    myObj.make_plot()

    myObj.set_plot_orientation('ns')
    myObj.make_plot()

    myObj.set_plot_orientation('z')
    myObj.make_plot()
