"""
Description
    view a model mesh (x y z) cell sizes.
    Input files format .ws, .rho, .prm

How to run
    python examples/view_modem_model.py /e/Data/Modeling/Isa/100hs_flat_BB/Isa_run3_NLCG_048.rho

Author: fei.zhang@ga.gov.au
Date: 2017-06-21
"""
import sys
import glob
import numpy as np
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

    line_color = kwargs.pop('line_color', 'k')
    line_width = kwargs.pop('line_width', .5)


    # make a rotation matrix to rotate data
    # cos_ang = np.cos(np.deg2rad(self.mesh_rotation_angle))
    # sin_ang = np.sin(np.deg2rad(self.mesh_rotation_angle))

    # turns out ModEM has not accomodated rotation of the grid, so for
    # now we will not rotate anything.
    cos_ang = 1
    sin_ang = 0


    # east_line_xlist = [-10,-6, -4, -2,0,2,4,6,10, None, -10,-6, -4, -2,0,2,4,6,10 ]
    # east_line_ylist = [1,1,1,1,1,1,1,1,1, None, 20,20,20,20,20,20,20,20,20]
    # ax1.scatter(east_line_xlist, east_line_ylist, lw=line_width, color=line_color)
    # plt.show()


    if modelfile is not None:
        with open(modelfile) as modfn:
            lines= modfn.readlines()

        print("How many lines in %s ?== %s "%(modelfile, len(lines)))

        for aline in lines[:5]:
            aline = aline.strip() # remove leading and trailing white spaces and invisible chars \n
            if aline.startswith('#'):
                print ("Header line skipped: ", aline)
            elif aline.endswith('LOGE'):
                (ew, ns, depth)= aline.split()[:3]
                print("cells = ", ew, ns, depth)
            else:  # cell sizes
                values = aline.split()
                print("cells: ",len(values))
                print('Cell sizes: ', values)

                nvalues= np.array(values).astype('float')
                plt.plot(nvalues, 'o',  lw=line_width, color=line_color)
                plt.show()

                # Check the changes in the cell sizes?
                diffval = np.diff(nvalues)
                print (diffval)
                plt.plot(diffval, '*',  lw=line_width, color=line_color)
                plt.show()

# ----------------------------------------------------------------------------
# view a model mesh (x y z) cell sizes
# python examples/view_modem_model.py /e/Data/Modeling/Isa/100hs_flat_BB/Isa_run3_NLCG_048.rho
if __name__ == "__main__":

    if len(sys.argv)>1:
        modelfile= sys.argv[1]
    else:
        modelfile = None
        print ("USAGE: python %s path2_model_file.ws|.mod|.rho" % sys.argv[0])
        sys.exit(1)

    plot_model_mesh(modelfile)

