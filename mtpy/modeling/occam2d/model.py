# -*- coding: utf-8 -*-
"""
Created on Tue Mar  7 19:11:57 2023

@author: jpeacock
"""
# =============================================================================
# Imports
# =============================================================================
import numpy as np

from mtpy.modeling.occam2d import Startup, Regularization

# =============================================================================
class Model(Startup):
    """
    Read .iter file output by Occam2d.  Builds the resistivity model from
    mesh and regularization files found from the .iter file.  The resistivity
    model is an array(x_nodes, z_nodes) set on a regular grid, and the values
    of the model response are filled in according to the regularization grid.
    This allows for faster plotting.

    Inherets Startup because they are basically the same object.

    Argument:
    ----------
        **iter_fn** : string
                      full path to .iter file to read. *default* is None.

        **model_fn** : string
                       full path to regularization file. *default* is None
                       and found directly from the .iter file.  Only input
                       if the regularization is different from the file that
                       is in the .iter file.

        **mesh_fn** : string
                      full path to mesh file. *default* is None
                      Found directly from the model_fn file.  Only input
                      if the mesh is different from the file that
                      is in the model file.

    ===================== =====================================================
    Key Words/Attributes  Description
    ===================== =====================================================
    data_fn               full path to data file
    iter_fn               full path to .iter file
    mesh_fn               full path to mesh file
    mesh_x                np.ndarray(x_nodes, z_nodes) mesh grid for plotting
    mesh_z                np.ndarray(x_nodes, z_nodes) mesh grid for plotting
    model_values          model values from startup file
    plot_x                nodes of mesh in horizontal direction
    plot_z                nodes of mesh in vertical direction
    res_model             np.ndarray(x_nodes, z_nodes) resistivity model
                          values in linear scale
    ===================== =====================================================


    ===================== =====================================================
    Methods               Description
    ===================== =====================================================
    build_model           get the resistivity model from the .iter file
                          in a regular grid according to the mesh file
                          with resistivity values according to the model file
    read_iter_file        read .iter file and fill appropriate attributes
    write_iter_file       write an .iter file incase you want to set it as the
                          starting model or a priori model
    ===================== =====================================================

    :Example: ::
        >>> model = occam2D.Model(r"/home/occam/line1/inv1/test_01.iter")
        >>> model.build_model()

    """

    def __init__(self, iter_fn=None, model_fn=None, mesh_fn=None, **kwargs):
        super().__init__(**kwargs)
        self.iter_fn = iter_fn
        self.model_fn = model_fn
        self.mesh_fn = mesh_fn
        self.data_fn = None
        self.model_values = None
        self.res_model = None
        self.plot_x = None
        self.plot_z = None
        self.mesh_x = None
        self.mesh_z = None

        for key, value in kwargs.items():
            setattr(self, key, value)

    def read_iter_file(self, iter_fn=None):
        """
        Read an iteration file.

        Arguments:
        ----------
            **iter_fn** : string
                        full path to iteration file if iterpath=None.  If
                        iterpath is input then iterfn is just the name
                        of the file without the full path.

        Returns:
        --------

        :Example: ::

            >>> import mtpy.modeling.occam2d as occam2d
            >>> itfn = r"/home/Occam2D/Line1/Inv1/Test_15.iter"
            >>> ocm = occam2d.Model(itfn)
            >>> ocm.read_iter_file()

        """

        if iter_fn is not None:
            self.iter_fn = iter_fn

        if self.iter_fn is None:
            raise ValueError("iter_fn is None, input iteration file")

        # check to see if the file exists
        if not self.iter_fn.exists():
            raise ValueError(f"Can not find {self.iter_fn}")

        self.save_path = self.iter_fn.parent

        # open file, read lines, close file

        with open(self.iter_fn, "r") as ifid:
            ilines = ifid.readlines()

        ii = 0
        # put header info into dictionary with similar keys
        while ilines[ii].lower().find("param") != 0:
            iline = ilines[ii].strip().split(":")
            key = iline[0].strip().lower()
            if key.find("!") != 0:
                key = (
                    key.replace(" ", "_")
                    .replace("file", "fn")
                    .replace("/", "_")
                )
                value = iline[1].strip()
                try:
                    setattr(self, key, float(value))
                except ValueError:
                    setattr(self, key, value)
            ii += 1

        # get number of parameters
        iline = ilines[ii].strip().split(":")
        key = iline[0].strip().lower().replace(" ", "_")
        value = int(iline[1].strip())
        setattr(self, key, value)

        self.model_values = np.zeros(self.param_count)
        kk = int(ii + 1)

        jj = 0
        mv_index = 0
        while jj < len(ilines) - kk:
            iline = np.array(ilines[jj + kk].strip().split(), dtype="float")
            self.model_values[mv_index : mv_index + iline.shape[0]] = iline
            jj += 1
            mv_index += iline.shape[0]

        # make sure data file is full path
        if not self.data_fn.is_file():
            self.data_fn = self.save_path.joinpath(self.data_fn)

        # make sure model file is full path
        if not self.model_fn.is_file():
            self.model_fn = self.save_path.joinpath(self.model_fn)

    def write_iter_file(self, iter_fn=None):
        """
        write an iteration file if you need to for some reason, same as
        startup file
        """
        if iter_fn is not None:
            self.iter_fn = iter_fn

        self.write_startup_file(iter_fn)

    def build_model(self):
        """
        build the model from the mesh, regularization grid and model file

        """
        # first read in the iteration file
        self.read_iter_file()

        # read in the regulariztion file
        r1 = Regularization()
        r1.read_regularization_file(self.model_fn)
        r1.model_rows = np.array(r1.model_rows)
        self.model_rows = r1.model_rows
        self.model_columns = r1.model_columns

        # read in mesh file
        r1.read_mesh_file(r1.mesh_fn)

        # get the binding offset which is the right side of the furthest left
        # block, this helps locate the model in relative space
        bndgoff = r1.binding_offset

        # make sure that the number of rows and number of columns are the same
        assert len(r1.model_rows) == len(r1.model_columns)

        # initiate the resistivity model to the shape of the FE mesh
        self.res_model = np.zeros((r1.z_nodes.shape[0], r1.x_nodes.shape[0]))

        # read in the model and set the regularization block values to map onto
        # the FE mesh so that the model can be plotted as an image or regular
        # mesh.
        mm = 0
        for ii in range(len(r1.model_rows)):
            # get the number of layers to combine
            # this index will be the first index in the vertical direction
            ny1 = r1.model_rows[:ii, 0].sum()
            # the second index  in the vertical direction
            ny2 = ny1 + r1.model_rows[ii][0]
            # make the list of amalgamated columns an array for ease
            lc = np.array(r1.model_columns[ii])
            # loop over the number of amalgamated blocks
            for jj in range(len(r1.model_columns[ii])):
                # get first in index in the horizontal direction
                nx1 = lc[:jj].sum()
                # get second index in horizontal direction
                nx2 = nx1 + lc[jj]
                # put the apporpriate resistivity value into all the amalgamated
                # model blocks of the regularization grid into the forward model
                # grid
                self.res_model[ny1:ny2, nx1:nx2] = self.model_values[mm]
                mm += 1

        # make some arrays for plotting the model
        self.plot_x = np.array(
            [r1.x_nodes[: ii + 1].sum() for ii in range(len(r1.x_nodes))]
        )
        self.plot_z = np.array(
            [r1.z_nodes[: ii + 1].sum() for ii in range(len(r1.z_nodes))]
        )

        # center the grid onto the station coordinates
        x0 = bndgoff - self.plot_x[r1.model_columns[0][0]]
        self.plot_x += x0

        # flip the arrays around for plotting purposes
        # plotx = plotx[::-1] and make the first layer start at zero
        self.plot_z = self.plot_z[::-1] - self.plot_z[0]

        # make a mesh grid to plot in the model coordinates
        self.mesh_x, self.mesh_z = np.meshgrid(self.plot_x, self.plot_z)

        # flip the resmodel upside down so that the top is the stations
        self.res_model = np.flipud(self.res_model)
