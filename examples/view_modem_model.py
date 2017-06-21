"""
Description:
    what does this script module do? How to do it.

Author: fei.zhang@ga.gov.au

Date:
"""
import sys
import glob

import matplotlib.pyplot as plt


def plot_model_mesh(modelfile, east_limits=None, north_limits=None, z_limits=None, **kwargs):
    """ Plot the model grid/mesh

    :param modelfile:  path2/ModEM_Model.ws

    :param east_limits: tuple (xmin,xmax)
                         plot min and max distances in meters for the
                         E-W direction.  If None, the east_limits
                         will be set to furthest stations east and west.
                         *default* is None
    :param north_limits: tuple (ymin,ymax)
                         plot min and max distances in meters for the
                         N-S direction.  If None, the north_limits
                         will be set to furthest stations north and south.
                         *default* is None
    :param z_limits:tuple (zmin,zmax)
                        plot min and max distances in meters for the
                        vertical direction.  If None, the z_limits is
                        set to the number of layers.  Z is positive down
                        *default* is None
    :param kwargs:
    :return:
    """

    fig_size = kwargs.pop('fig_size', [6, 4])
    fig_dpi = kwargs.pop('fig_dpi', 300)
    fig_num = kwargs.pop('fig_num', 1)

    station_marker = kwargs.pop('station_marker', 'v')
    marker_color = kwargs.pop('station_color', 'b')
    marker_size = kwargs.pop('marker_size', 2)

    line_color = kwargs.pop('line_color', 'k')
    line_width = kwargs.pop('line_width', .5)

    plt.rcParams['figure.subplot.hspace'] = .3
    plt.rcParams['figure.subplot.wspace'] = .3
    plt.rcParams['figure.subplot.left'] = .12
    plt.rcParams['font.size'] = 7

    fig = plt.figure(fig_num, figsize=fig_size, dpi=fig_dpi)
    plt.clf()

    # make a rotation matrix to rotate data
    # cos_ang = np.cos(np.deg2rad(self.mesh_rotation_angle))
    # sin_ang = np.sin(np.deg2rad(self.mesh_rotation_angle))

    # turns out ModEM has not accomodated rotation of the grid, so for
    # now we will not rotate anything.
    cos_ang = 1
    sin_ang = 0

    # --->plot map view
    ax1 = fig.add_subplot(1, 2, 1, aspect='equal')

    # 1. plot station locations
    # plot_east = station_locations['rel_east']
    # plot_north = station_locations['rel_north']
    #
    # ax1.scatter(plot_east,
    #             plot_north,
    #             marker=station_marker,
    #             c=marker_color,
    #             s=marker_size)

    # 2. plot east lines and north lines
    east_line_xlist = [-10,-6, -4, -2,0,2,4,6,10, None, -10,-6, -4, -2,0,2,4,6,10 ]
    east_line_ylist = [1,1,1,1,1,1,1,1,1, None, 20,20,20,20,20,20,20,20,20]


    ax1.scatter(east_line_xlist,
             east_line_ylist,
             lw=line_width,
             color=line_color)

    #
    # if east_limits == None:
    #     ax2.set_xlim(plot_east.min() - 50 * self.cell_size_east,
    #                  plot_east.max() + 50 * self.cell_size_east)
    # else:
    #     ax2.set_xlim(east_limits)

    # ax2.set_ylabel('Depth (m)', fontdict={'size': 9, 'weight': 'bold'})
    # ax2.set_xlabel('Easting (m)', fontdict={'size': 9, 'weight': 'bold'})
    plt.show()


    if modelfile is not None:
        with open(modelfile) as modfn:
            lines= modfn.readlines()

        print("How many lines in %s ?== %s "%(modelfile, len(lines)))

        for aline in lines[:5]:
            if aline.startswith('#'):
                print ("Header line skipped: ", aline)
            elif aline.endswith('LOGE\n'):
                (ew, ns, depth)= aline.split()[:3]
                print(ew, ns, depth)
            else:  # cell sizes
                values = aline.split()
                print(len(values))
                print(values)

                plt.plot(values, 'o',  lw=line_width, color=line_color)
                plt.show()


if __name__ == "__main__":

    if len(sys.argv)>1:
        modelfile= sys.argv[1]
    else:
        modelfile = None

    plot_model_mesh(modelfile)

