# -*- coding: utf-8 -*-
"""
Created on Tue Mar  7 18:13:52 2023

@author: jpeacock
"""
# =============================================================================
# Imports
# =============================================================================
from pathlib import Path

import numpy as np

from mtpy.modeling.occam2d import Mesh

# =============================================================================
class Regularization(Mesh):
    """
    Creates a regularization grid based on Mesh.  Note that Mesh is inherited
    by Regularization, therefore the intended use is to build a mesh with
    the Regularization class.

    The regularization grid is what Occam calculates the inverse model on.
    Setup is tricky and can be painful, as you can see it is not quite fully
    functional yet, as it cannot incorporate topography yet.  It seems like
    you'd like to have the regularization setup so that your target depth is
    covered well, in that the regularization blocks to this depth are
    sufficiently small to resolve resistivity structure at that depth.
    Finally, you want the regularization to go to a half space at the bottom,
    basically one giant block.

    Arguments:
    -----------
        **station_locations** : np.ndarray(n_stations)
                                array of station locations along a profile
                                line in meters.

    ======================= ===================================================
    Key Words/Attributes    Description
    ======================= ===================================================
    air_key                 letter associated with the value of air
                            *default* is 0
    air_value               value given to an air cell, *default* is 1E13
    binding_offset          offset from the right side of the furthest left
                            hand model block in meters.  The regularization
                            grid is setup such that this should be 0.
    cell_width              width of cells with in station area in meters
                            *default* is 100
    description             description of the model for the model file.
                            *default* is 'simple inversion'
    elevation_profile       elevation profile along the profile line.
                            given as np.ndarray(nx, 2), where the elements
                            are x_location, elevation.  If elevation profile
                            is given add_elevation is called automatically.
                            *default* is None
    mesh_fn                 full path to mesh file.
    mesh_values             letter values of each triangular mesh element
                            if the cell is free value is ?
    model_columns
    model_name
    model_rows

    min_block_width         [ float ] minimum model block width in meters,
                            *default* is 2*cell_width
    n_layers                number of vertical layers in mesh
                            *default* is 90
    num_free_param          [ int ] number of free parameters in the model.
                            this is a tricky number to estimate apparently.
    num_layers              [ int ] number of regularization layers.
    num_x_pad_cells         number of horizontal padding cells outside the
                            the station area that will increase in size
                            by x_pad_multiplier. *default* is 7
    num_x_pad_small_cells   number of horizonal padding cells just outside
                            the station area with width cell_width.  This is
                            to extend the station area if needed.
                            *default* is 2
    num_z_pad_cells         number of vertical padding cells below
                            z_target_depth down to z_bottom. *default* is 5
    prejudice_fn            full path to prejudice file
                            *default* is 'none'
    reg_basename            basename of regularization file (model file)
                            *default* is 'Occam2DModel'
    reg_fn                  full path to regularization file (model file)
                            *default* is save_path/reg_basename
    rel_station_locations   relative station locations within the mesh.  The
                            locations are relative to the center of the station
                            area.  *default* is None, filled later
    save_path               full path to save mesh and model file to.
                            *default* is current working directory.
    statics_fn              full path to static shift file
                            Static shifts in occam may not work.
                            *default* is 'none'
    station_locations       location of stations in meters, can be on a
                            relative grid or in UTM.
    trigger                 [ float ] multiplier to merge model blocks at
                            depth.  A higher number increases the number of
                            model blocks at depth.  *default* is .65
    x_grid                  location of horizontal grid nodes in meters
    x_nodes                 relative spacing between grid nodes
    x_pad_multiplier        horizontal padding cells will increase by this
                            multiple out to the edge of the grid.
                            *default* is 1.5
    z1_layer                thickness of the first layer in the model.
                            Should be at least 1/4 of the first skin depth
                            *default* is 10
    z_bottom                bottom depth of the model (m).  Needs to be large
                            enough to be 1D at the edge.
                            *default* is 200000.0
    z_grid                  location of vertical nodes in meters
    z_nodes                 relative distance between vertical nodes in meters
    z_target_depth          depth to deepest target of interest.  Below this
                            depth cells will be padded to z_bottom
    ======================= ===================================================

    .. note:: regularization does not work with topography yet.  Having
              problems calculating the number of free parameters.

    ========================= =================================================
    Methods                   Description
    ========================= =================================================
    add_elevation             adds elevation to the mesh given elevation
                              profile.
    build_mesh                builds the mesh given the attributes of Mesh.  If
                              elevation_profile is not None, add_elevation is
                              called inside build_mesh
    build_regularization      builds the regularization grid from the build mesh
                              be sure to plot the grids before starting the
                              inversion to make sure coverage is appropriate.
    get_num_free_param        estimate the number of free parameters.
                              **This is a work in progress**
    plot_mesh                 plots the built mesh with station location.
    read_mesh_file            reads in an existing mesh file and populates the
                              appropriate attributes.
    read_regularization_file  read in existing regularization file, populates
                              apporopriate attributes
    write_mesh_file           writes a mesh file to save_path
    write_regularization_file writes a regularization file
    ======================= ===================================================

    :Example: ::

        >>> edipath = r"/home/mt/edi_files"
        >>> profile = occam2d.Profile(edi_path=edi_path)
        >>> profile.generate_profile()
        >>> reg = occam2d.Regularization(profile.station_locations)
        >>> reg.build_mesh()
        >>> reg.build_regularization()
        >>> reg.save_path = r"/home/occam2d/Line1/Inv1"
        >>> reg.write_regularization_file()

    """

    def __init__(self, station_locations=None, **kwargs):
        # Be sure to initialize Mesh
        super().__init__(station_locations, **kwargs)

        self.min_block_width = 2 * np.median(self.cell_width)
        self.trigger = 0.75
        self.model_columns = None
        self.model_rows = None
        self.binding_offset = None
        self.reg_fn = None
        self.reg_basename = "Occam2DModel"
        self.model_name = "model made by mtpy.modeling.occam2d"
        self.description = "simple Inversion"
        self.num_param = None
        self.num_free_param = None
        self.statics_fn = "none"
        self.prejudice_fn = "none"
        self.num_layers = None

        # --> build mesh
        if self.station_locations is not None:
            self.build_mesh()
            self.build_regularization()

    def __str__(self):
        lines = []
        lines.append("=" * 55)
        lines.append(f"{'Regularization Parameters':^55}")
        lines.append("=" * 55)
        lines.append(f"\tbinding offset        = {self.binding_offset:.1f}")
        lines.append(f"\tnumber layers        = {len(self.model_columns)}")
        lines.append(f"\tnumber of parameters = {self.num_param}")
        lines.append(f"\tnumber of free param = {self.num_free_param}")
        lines.append("=" * 55)

        return "\n".join(lines)

    def __repr__(self):
        return self.__str__()

    def build_regularization(self):
        """
        Builds larger boxes around existing mesh blocks for the regularization.
        As the model deepens the regularization boxes get larger.

        The regularization boxes are merged mesh cells as prescribed by the
        Occam method.

        """
        # list of the mesh columns to combine
        self.model_columns = []
        # list of mesh rows to combine
        self.model_rows = []

        # At the top of the mesh model blocks will be 2 combined mesh blocks
        # Note that the padding cells are combined into one model block
        station_col = [2] * int(
            (self.x_nodes.shape[0] - 2 * self.num_x_pad_cells + 1) / 2
        )
        model_cols = (
            [self.num_x_pad_cells] + station_col + [self.num_x_pad_cells]
        )
        station_widths = [
            self.x_nodes[ii] + self.x_nodes[ii + 1]
            for ii in range(
                self.num_x_pad_cells,
                self.x_nodes.shape[0] - self.num_x_pad_cells,
                2,
            )
        ]

        pad_width = self.x_nodes[0 : self.num_x_pad_cells].sum()
        model_widths = [pad_width] + station_widths + [pad_width]
        num_cols = len(model_cols)

        model_thickness = np.hstack(
            [
                self.z_nodes[:2].sum(),
                self.z_nodes[2 : self.z_nodes.shape[0] - self.num_z_pad_cells],
                self.z_nodes[-self.num_z_pad_cells :].sum(),
            ]
        )

        self.num_param = 0
        # --> now need to calulate model blocks to the bottom of the model
        columns = list(model_cols)
        widths = list(model_widths)
        for zz, thickness in enumerate(model_thickness):
            # index for model column blocks from first_row, start at 1 because
            # 0 is for padding cells
            block_index = 1
            num_rows = 1
            if zz == 0:
                num_rows += 1
            if zz == len(model_thickness) - 1:
                num_rows = self.num_z_pad_cells
            while block_index + 1 < num_cols - 1:
                # check to see if horizontally merged mesh cells are not larger
                # than the thickness times trigger
                if thickness < self.trigger * (
                    widths[block_index] + widths[block_index + 1]
                ):
                    block_index += 1
                    continue
                # merge 2 neighboring cells to avoid vertical exaggerations
                else:
                    widths[block_index] += widths[block_index + 1]
                    columns[block_index] += columns[block_index + 1]
                    # remove one of the merged cells
                    widths.pop(block_index + 1)
                    columns.pop(block_index + 1)

                    num_cols -= 1
            self.num_param += num_cols

            self.model_columns.append(list(columns))
            self.model_rows.append([num_rows, num_cols])

        # calculate the distance from the right side of the furthest left
        # model block to the furthest left station which is half the distance
        # from the center of the mesh grid.
        self.binding_offset = (
            self.x_grid[self.num_x_pad_cells + 1]
            + self.station_locations.mean()
        )

        self.get_num_free_params()

    def get_num_free_params(self):
        """
        estimate the number of free parameters in model mesh.

        I'm assuming that if there are any fixed parameters in the block, then
        that model block is assumed to be fixed. Not sure if this is right
        cause there is no documentation.

        **DOES NOT WORK YET**
        """

        self.num_free_param = 0

        row_count = 0
        # loop over columns and rows of regularization grid
        for col, row in zip(self.model_columns, self.model_rows):
            rr = row[0]
            col_count = 0
            for ii, cc in enumerate(col):
                # make a model block from the index values of the regularization
                # grid
                model_block = self.mesh_values[
                    row_count : row_count + rr, col_count : col_count + cc, :
                ]

                # find all the free triangular blocks within that model block
                find_free = np.where(model_block == "?")
                try:
                    # test to see if the number of free parameters is equal
                    # to the number of triangular elements with in the model
                    # block, if there is the model block is assumed to be free.
                    if find_free[0].size == model_block.size:
                        self.num_free_param += 1
                except IndexError:
                    pass
                col_count += cc
            row_count += rr

    def write_regularization_file(
        self,
        reg_fn=None,
        reg_basename=None,
        statics_fn="none",
        prejudice_fn="none",
        save_path=None,
    ):
        """
        Write a regularization file for input into occam.

        Calls build_regularization if build_regularization has not already
        been called.

        if reg_fn is None, then file is written to save_path/reg_basename

        Arguments:
        ----------
            **reg_fn** : string
                         full path to regularization file. *default* is None
                         and file will be written to save_path/reg_basename

            **reg_basename** : string
                               basename of regularization file

            **statics_fn** : string
                             full path to static shift file
                             .. note:: static shift does not always work in
                                       occam2d.exe
            **prejudice_fn** : string
                               full path to prejudice file

            **save_path** : string
                            path to save regularization file.
                            *default* is current working directory

        """
        if save_path is not None:
            self.save_path = Path(save_path)
        if reg_basename is not None:
            self.reg_basename = reg_basename
        if reg_fn is None:
            if self.save_path is None:
                self.save_path = Path()
            self.reg_fn = self.save_path.joinpath(self.reg_basename)

        self.statics_fn = statics_fn
        self.prejudice_fn = prejudice_fn

        if self.model_columns is None:
            if self.binding_offset is None:
                self.build_mesh()
            self.build_regularization()

        reg_lines = []

        # --> write out header information
        reg_lines.append(f"{'Format:':<18}{'occam2mtmod_1.0'.upper()}")
        reg_lines.append(f"{'Model Name:':<18}{self.model_name.upper()}")
        reg_lines.append(f"{'Description':<18}{self.description.upper()}")

        if self.mesh_fn.parent == self.save_path:
            reg_lines.append(f"{'Mesh File:':<18}{self.mesh_fn.name}")
        else:
            reg_lines.append(f"{'Mesh File:':<18}{self.mesh_fn}")
        reg_lines.append(f"{'Mesh Type':<18}{'pw2d'.upper()}")
        if self.statics_fn.parent == self.save_path:
            reg_lines.append(f"{'Statics File':<18}{self.statics_fn.name}")

        else:
            reg_lines.append(f"{'Statics File':<18}{self.statics_fn}")
        if self.prejudice_fn.parent == self.save_path:
            reg_lines.append(f"{'Prejudice File':<18}{self.prejudice_fn.name}")
        else:
            reg_lines.append(f"{'Prejudice File':<18}{self.prejudice_fn}")
        reg_lines.append(f"{'Binding Offset':<20}{self.binding_offset:.1f}")
        reg_lines.append(f"{'Num Layers':<20}{len(self.model_columns)}")

        # --> write rows and columns of regularization grid
        for row, col in zip(self.model_rows, self.model_columns):
            reg_lines.append("".join([f" {rr:>5}" for rr in row]) + "\n")
            reg_lines.append("".join([f"{cc:>5}" for cc in col]) + "\n")

        reg_lines.append(f"{'NO. EXCEPTIONS:':<18}0")

        with open(self.reg_fn, "w") as rfid:
            rfid.write("\n".join(reg_lines))

        print("Wrote Regularization file to {self.reg_fn}")

    def read_regularization_file(self, reg_fn):
        """
        Read in a regularization file and populate attributes:
            * binding_offset
            * mesh_fn
            * model_columns
            * model_rows
            * prejudice_fn
            * statics_fn

        """
        self.reg_fn = Path(reg_fn)
        self.save_path = self.reg_fn.parent

        with open(self.reg_fn, "r") as rfid:
            rlines = rfid.readlines()

        self.model_rows = []
        self.model_columns = []
        ncols = []

        for ii, iline in enumerate(rlines):
            # read header information
            if iline.find(":") > 0:
                iline = iline.strip().split(":")
                key = iline[0].strip().lower()
                key = key.replace(" ", "_").replace("file", "fn")
                value = iline[1].strip()
                try:
                    setattr(self, key, float(value))
                except ValueError:
                    setattr(self, key, value)

                # append the last line
                if key.find("exception") > 0:
                    self.model_columns.append(ncols)

            # get mesh values
            else:
                iline = iline.strip().split()
                iline = [int(jj) for jj in iline]
                if len(iline) == 2:
                    if len(ncols) > 0:
                        self.model_columns.append(ncols)
                    self.model_rows.append(iline)
                    ncols = []
                elif len(iline) > 2:
                    ncols = ncols + iline

        # set mesh file name
        if not self.mesh_fn.is_file():
            self.mesh_fn = self.save_path.joinpath(self.mesh_fn)
            self.statics_fn = self.save_path.joinpath(self.statics_fn)
            self.prejudice_fn = self.save_path.joinpath(self.prejudice_fn)
