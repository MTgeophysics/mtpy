"""
Description:
    Extract data from modem output files (.dat and .rho) and produce csv files for further visualization
    The output look like: StationName, Lat, Long, X, Y, Z, Log(Resistivity)
    where (X,Y,Z) are relative distances in meters from the mesh's origin.
    Projection/Coordinate system must be known in order to associate (Lat, Long) to (X, Y)

References:
    https://gajira.atlassian.net/browse/ALAMP-30
    https://gajira.atlassian.net/browse/ALAMP-31

CreationDate:   8/09/2017
CreatedBy:      SysUser='u25656'

Developer:      fei.zhang@ga.gov.au
LastUpdate:     8/09/2017
"""

import os
import sys
import glob
import logging
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import MultipleLocator
from mpl_toolkits.axes_grid1 import make_axes_locatable

from mtpy.modeling.modem import Data
from mtpy.modeling.modem import Model

from mtpy.utils.mtpylog import MtPyLog

logger = MtPyLog().get_mtpy_logger(__name__)
logger.setLevel(logging.DEBUG)

class ModemSlices():

    def __init__(self, filedat, filerho, plot_orient='ew', **kwargs):
        """Constructor
        :param filedat: path2file.dat
        :param filerho: path2file.rho
        :param  plot_orient: plot orientation ['ew','ns', 'z']
        """
        self.datfile = filedat
        self.rhofile = filerho

        # plot orientation 'ns' (north-south),'ew' (east-west) or
        # 'z' (horizontal slice))
        self.plot_orientation = plot_orient

        # slice location, in local grid coordinates (if it is a z slice, this
        # is slice depth)
        # self.slice_location = kwargs.pop('slice_location', 1000)
        # maximum distance in metres from vertical slice location and station
        self.station_dist = kwargs.pop('station_dist', 50000)
        # z limits (positive down so order is reversed)
        self.zlim = kwargs.pop('zlim', (200000, -2000))
        # colour limits
        self.clim = kwargs.pop('clim', [0.3, 3.7])
        self.fig_size = kwargs.pop('fig_size', [12, 10])
        self.font_size = kwargs.pop('font_size', 16)
        self.border_linewidth = 2

        self.map_scale = kwargs.pop('map_scale', 'm')
        # make map scale
        if self.map_scale == 'km':
            self.dscale = 1000.
        elif self.map_scale == 'm':
            self.dscale = 1.
        else:
            print("Unknown map scale:", self.map_scale)

        self.xminorticks = kwargs.pop('xminorticks', 10000)
        self.yminorticks = kwargs.pop('yminorticks', 10000)

        # read in the model data-file and rho-file
        self._read_model_data()

        return

    def _read_model_data(self):

        self.datObj = Data()
        self.datObj.read_data_file(data_fn=self.datfile)

        self.modObj = Model(model_fn=self.rhofile)
        self.modObj.read_model_file()

        return

    def set_plot_orientation(self, orient):
        """set a new plot orientation for plotting
        :param orient: z, ew, ns
        :return:
        """
        if orient in ['z','ew','ns']:
            self.plot_orientation = orient
        else:
            raise Exception("Error: unknown orientation value= %s"%orient)


    def get_slice_data(self, slice_location):
        """
        get the resistivity slices at the specified location
        :param slice_location:
        :return: slice data
        """

        # get grid centres (finite element cells centres)
        gcz = np.mean([self.modObj.grid_z[:-1],
                       self.modObj.grid_z[1:]], axis=0)
        gceast, gcnorth = [np.mean([arr[:-1], arr[1:]], axis=0) for arr in
                           [self.modObj.grid_east, self.modObj.grid_north]]

          # distance from slice to grid centre locations
        if self.plot_orientation == 'ew':
            sdist = np.abs(gcnorth - slice_location)
        elif self.plot_orientation == 'ns':
            sdist = np.abs(gceast - slice_location)
        elif self.plot_orientation == 'z':
            sdist = np.abs(gcz - slice_location)
        # find the closest slice index to specified location
        snos = np.where(sdist == np.amin(sdist))
        print (type(snos), len(snos))  # ((index1), (index2), (index3))

        # unpack the index tupple, and get the integer value as index number
        sno=snos[0][0]
        print("the slice index number=", sno)
        # get data for plotting
        if self.plot_orientation == 'ew':
            X, Y, res = self.modObj.grid_east, self.modObj.grid_z, np.log10(
                self.modObj.res_model[sno, :, :].T)
            ss = np.where(np.abs(self.datObj.station_locations['rel_north'] - np.median(gcnorth)) < self.station_dist)[
                0]

            sX, sY = self.datObj.station_locations['rel_east'][
                         ss], self.datObj.station_locations['elev'][ss]
            xlim = (self.modObj.grid_east[
                        self.modObj.pad_east], self.modObj.grid_east[-self.modObj.pad_east - 1])
            ylim = self.zlim
            title = 'East-west slice at {} meters north'.format(gcnorth[sno])
        elif self.plot_orientation == 'ns':
            X, Y, res = self.modObj.grid_north, self.modObj.grid_z, np.log10(
                self.modObj.res_model[:, sno, :].T)
            # indices for selecting stations close to profile
            ss = np.where(
                np.abs(
                    self.datObj.station_locations['rel_east'] -
                    np.median(gceast)) < self.station_dist)[0]

            sX, sY = self.datObj.station_locations['rel_north'][
                         ss], self.datObj.station_locations['elev'][ss]
            xlim = (self.modObj.grid_north[
                        self.modObj.pad_north], self.modObj.grid_north[-self.modObj.pad_north - 1])
            ylim = self.zlim
            title = 'North-south slice at {} meters east'.format(gceast[sno])
        elif self.plot_orientation == 'z':
            X, Y, res = self.modObj.grid_east, self.modObj.grid_north, np.log10(
                self.modObj.res_model[:, :, sno])
            sX, sY = self.datObj.station_locations[
                         'rel_east'], self.datObj.station_locations['rel_north']
            xlim = (self.modObj.grid_east[
                        self.modObj.pad_east], self.modObj.grid_east[-self.modObj.pad_east - 1])
            ylim = (self.modObj.grid_north[
                        self.modObj.pad_north], self.modObj.grid_north[-self.modObj.pad_north - 1])
            title = 'Horizontal Slice at Depth {} meters'.format(gcz[sno])


        return (X,Y,res,sX,sY,xlim,ylim,title)


    def create_csv(self):
        """
        dump into CSV file with the output columns:
            StationName, Lat, Long, X, Y, Z, Log(Resistivity)
        where (X,Y,Z) are relative distances in meters from the mesh's origin.
        Projection/Coordinate system must be known in order to associate (Lat, Long) to (X, Y)
        :return:
        """
        self.set_plot_orientation('z')

        (X, Y, res, sX, sY, xlim, ylim, title) = self.get_slice_data(1000)

        #print (X,Y,res)
        print(sX,sY)

    def make_plot(self,slice_location=1000):
        """ create a plot based on the input data and parameters
        :return:
        """

        (X, Y, res, sX, sY, xlim, ylim, title)= self.get_slice_data(slice_location)

        # make the plot

        fdict = {'size': self.font_size, 'weight': 'bold'}
        plt.figure(figsize=self.fig_size)
        plt.rcParams['font.size'] = self.font_size

        # plot station locations
        # print("station locations sX:", sX)
        # print("station locations sY:", sY)

        plt.plot(sX, sY, 'kv')  # station marker:'kv'

        mesh_plot = plt.pcolormesh(X, Y, res, cmap='bwr_r')

        xlim2 = (xlim[0] / self.dscale, xlim[1] / self.dscale)
        ylim2 = (ylim[0] / self.dscale, ylim[1] / self.dscale)

        plt.xlim(*xlim)
        plt.ylim(*ylim)

        # set title
        plt.title(title, fontdict=fdict)

        #if self.plot_orientation == 'z':
        plt.gca().set_aspect('equal', 'datalim')

        plt.clim(*self.clim)
        # plt.colorbar()

        # FZ: fix miss-placed colorbar
        ax = plt.gca()
        ax.xaxis.set_minor_locator(
            MultipleLocator(
                self.xminorticks))  # /self.dscale
        ax.yaxis.set_minor_locator(
            MultipleLocator(
                self.yminorticks))  # /self.dscale
        ax.tick_params(axis='both', which='minor', width=2, length=5)
        ax.tick_params(
            axis='both',
            which='major',
            width=3,
            length=15,
            labelsize=20)
        for axis in ['top', 'bottom', 'left', 'right']:
            ax.spines[axis].set_linewidth(self.border_linewidth)
        # ax.tick_params(axis='both', which='major', labelsize=20)
        # ax.tick_params(axis='both', which='minor', labelsize=20)

        # http://stackoverflow.com/questions/10171618/changing-plot-scale-by-a-factor-in-matplotlib
        xticks = ax.get_xticks() / self.dscale
        ax.set_xticklabels(xticks)
        yticks = ax.get_yticks() / self.dscale
        ax.set_yticklabels(yticks)

        # create an axes on the right side of ax. The width of cax will be 5%
        # of ax and the padding between cax and ax will be fixed at 0.05 inch.
        divider = make_axes_locatable(ax)
        # pad = separation from figure to colorbar
        cax = divider.append_axes("right", size="5%", pad=0.2)

        mycb = plt.colorbar( mesh_plot,cax=cax)
        mycb.outline.set_linewidth(self.border_linewidth)
        mycb.set_label('Resistivity ($\Omega \cdot$m)', fontdict=fdict)

        if self.plot_orientation == 'z':
            ax.set_ylabel('Northing (' + self.map_scale + ')', fontdict=fdict)
            ax.set_xlabel('Easting (' + self.map_scale + ')', fontdict=fdict)
        if self.plot_orientation == 'ew':
            ax.set_ylabel('Depth (' + self.map_scale + ')', fontdict=fdict)
            ax.set_xlabel('Easting (' + self.map_scale + ')', fontdict=fdict)
        if self.plot_orientation == 'ns':
            ax.set_ylabel('Depth (' + self.map_scale + ')', fontdict=fdict)
            ax.set_xlabel('Northing (' + self.map_scale + ')', fontdict=fdict)

        plt.show()

        return


#########################################################################
if __name__ == "__main__":
    """ Usage:
    python mtpy/modeling/modem_outfiles_to_csv.py
    /e/Data/Modeling/Isa/100hs_flat_BB/Isa_run3_NLCG_048.dat /e/Data/Modeling/Isa/100hs_flat_BB/Isa_run3_NLCG_048.rho 20

    python mtpy/modeling/modem_outfiles_to_csv.py
    /e/tmp/GA_UA_edited_10s-10000s_16/ModEM_Data.dat  /e/tmp/GA_UA_edited_10s-10000s_16/ModEM_Model.ws 5000
    """

    # Take commandline input
    if (len(sys.argv) == 2):  # A model dir provided
        modeldir = sys.argv[1]
        datf = os.path.join(modeldir, 'ModEM_Data.dat')
        rhofiles = glob.glob(os.path.join(modeldir, '*.rho'))

        print(rhofiles)

        if len(rhofiles) < 1:
            print ("No rho files found in the dir %s", modeldir)
            sys.exit(1)
        else:
            # the file with highest numbers in the last 3 numbers before *.rho
            rhof = sorted(rhofiles)[-1]

        print("Effective Files Used in Plot: ", datf, rhof)

    # dat and rho file both provided
    if len(sys.argv) >= 3:
        datf = sys.argv[1]
        rhof = sys.argv[2]

    if len(sys.argv) >= 4:
        slice_locs = sys.argv[3:]
    else:
        # a list of depth where h-slice to be visualized
        slice_locs=[-2000, -1900, -1700, -1500, -1200, -1000, -800, -600, -400, -200,
                    0, 20, 50, 80, 100,150, 200, 400, 600,800,1000,
                    2000,3000,4000,5000,6000,7000,8000,9000,10000]

    # construct plot object
    myObj = ModemSlices(datf, rhof)  # ,map_scale='m')

    myObj.create_csv()

    # myObj.set_plot_orientation('ew')
    # myObj.set_plot_orientation('ns')
    # horizontal at a given depth
    myObj.set_plot_orientation('z')

    for dist in slice_locs:
        sdist=int(dist)

        print("**** plotting slice at location: ****", sdist)
        print("**** actual location will be at nearest cell centre ****")

        # plot resistivity image at slices in three orientations at a given slice_location=sdist

        myObj.make_plot(slice_location=sdist)  # actual location will be nearest cell centre

        # myObj.set_plot_orientation('ew')
        # myObj.make_plot(slice_location=sdist)
        #
        # myObj.set_plot_orientation('ns')
        # myObj.make_plot(slice_location=sdist)

        plt.show()
