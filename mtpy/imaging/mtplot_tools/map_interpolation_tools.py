# -*- coding: utf-8 -*-
"""
Created on Sun Sep 25 15:06:58 2022

@author: jpeacock
"""

# =============================================================================
# Imports
# =============================================================================
import numpy as np

import matplotlib.tri as tri

from scipy.spatial import Delaunay, cKDTree
from scipy import interpolate

# =============================================================================
def in_hull(p, hull):
    """
    Test if points in p are within the convex hull
    """

    try:
        if not isinstance(hull, Delaunay):
            hull = Delaunay(hull)
        return hull.find_simplex(p) >= 0
    except:
        from scipy.optimize import linprog

        # Delaunay triangulation will fail if there are collinear points;
        # in those instances use linear programming (much slower) to define
        # a convex hull.
        def in_hull_lp(points, x):
            """
            :param points:
            :param x:
            :return:
            """
            n_points = len(points)
            c = np.zeros(n_points)
            A = np.r_[points.T, np.ones((1, n_points))]
            b = np.r_[x, np.ones(1)]
            lp = linprog(c, A_eq=A, b_eq=b)
            return not lp.success

        result = []
        for cp in p:
            result.append(in_hull_lp(hull, cp))

        return np.array(result)


def get_plot_xy(plot_array, cell_size, n_padding_cells):
    """
    Get plot x and plot y from a plot array to interpolate on to

    :param plot_array: DESCRIPTION
    :type plot_array: TYPE
    :param cell_size: DESCRIPTION
    :type cell_size: TYPE
    :param n_padding_cells: DESCRIPTION
    :type n_padding_cells: TYPE
    :return: DESCRIPTION
    :rtype: TYPE

    """

    # create uniform x, y to plot on.
    ds = cell_size * n_padding_cells
    n_points = int(
        abs(
            plot_array["longitude"].max()
            - plot_array["longitude"].min()
            + 2 * ds
        )
        / cell_size
    )
    plot_x = np.linspace(
        plot_array["longitude"].min() - ds,
        plot_array["longitude"].max() + ds,
        n_points,
    )

    n_points = int(
        abs(
            plot_array["latitude"].max() - plot_array["latitude"].min() + 2 * ds
        )
        / cell_size
    )
    plot_y = np.linspace(
        plot_array["latitude"].min() - ds,
        plot_array["latitude"].max() + ds,
        n_points,
    )

    return plot_x, plot_y


def griddata_interpolate(
    x, y, values, new_x, new_y, interpolation_method="cubic"
):
    """

    :param x: DESCRIPTION
    :type x: TYPE
    :param y: DESCRIPTION
    :type y: TYPE
    :param value: DESCRIPTION
    :type value: TYPE
    :param new_x: DESCRIPTION
    :type new_x: TYPE
    :param new_y: DESCRIPTION
    :type new_y: TYPE
    :param interpolation_method: DESCRIPTION, defaults to "cubic"
    :type interpolation_method: TYPE, optional
    :return: DESCRIPTION
    :rtype: TYPE

    """
    points = np.array([x, y]).T
    grid_x, grid_y = np.meshgrid(new_x, new_y)

    return (
        grid_x,
        grid_y,
        interpolate.griddata(
            points,
            values,
            (grid_x, grid_y),
            method=interpolation_method,
        ),
    )


def interpolate_to_map_griddata(
    plot_array,
    component,
    cell_size=0.002,
    n_padding_cells=10,
    interpolation_method="cubic",
):
    """
    Interpolate using scipy.interpolate.griddata

    :param plot_array: DESCRIPTION
    :type plot_array: TYPE
    :param component: DESCRIPTION
    :type component: TYPE
    :param cell_size: DESCRIPTION, defaults to .002
    :type cell_size: TYPE, optional
    :param n_padding_cells: DESCRIPTION, defaults to 10
    :type n_padding_cells: TYPE, optional
    :return: DESCRIPTION
    :rtype: TYPE

    """

    plot_x, plot_y = get_plot_xy(plot_array, cell_size, n_padding_cells)

    grid_x, grid_y, image = griddata_interpolate(
        plot_array["longitude"],
        plot_array["latitude"],
        plot_array[component],
        plot_x,
        plot_y,
    )

    if "res" in component:
        image = np.log10(image)

    return grid_x, grid_y, image


def triangulate_interpolation(
    x,
    y,
    values,
    padded_x,
    padded_y,
    new_x,
    new_y,
    nearest_neighbors=7,
    interp_pow=4,
):
    """

    :param x: DESCRIPTION
    :type x: TYPE
    :param y: DESCRIPTION
    :type y: TYPE
    :param values: DESCRIPTION
    :type values: TYPE
    :param new_x: DESCRIPTION
    :type new_x: TYPE
    :param new_y: DESCRIPTION
    :type new_y: TYPE
    :param : DESCRIPTION
    :type : TYPE
    :return: DESCRIPTION
    :rtype: TYPE

    """

    grid_x, grid_y = np.meshgrid(new_x, new_y)
    grid_x = grid_x.flatten()
    grid_y = grid_y.flatten()

    triangulation = tri.Triangulation(grid_x, grid_y)

    mean_x = grid_x[triangulation.triangles].mean(axis=1)
    mean_y = grid_y[triangulation.triangles].mean(axis=1)

    points_mean = np.array([mean_x, mean_y]).T
    padded_points = np.array([padded_x, padded_y]).T

    inside_indices = in_hull(points_mean, padded_points)
    inside_indices = np.bool_(inside_indices)
    triangulation.set_mask(~inside_indices)

    tree = cKDTree(np.array([x, y]).T)

    xy = np.array([grid_x, grid_y]).T
    d, l = tree.query(xy, k=nearest_neighbors)

    image = None

    if nearest_neighbors == 1:
        # extract nearest neighbour values
        image = values[l]
    else:
        image = np.zeros((xy.shape[0]))

        # field values are directly assigned for coincident locations
        coincident_indices = d[:, 0] == 0
        image[coincident_indices] = values[l[coincident_indices, 0]]

        # perform idw interpolation for non-coincident locations
        idw_indices = d[:, 0] != 0
        w = np.zeros(d.shape)
        w[idw_indices, :] = 1.0 / np.power(d[idw_indices, :], interp_pow)

        image[idw_indices] = np.sum(
            w[idw_indices, :] * values[l[idw_indices, :]], axis=1
        ) / np.sum(w[idw_indices, :], axis=1)

    return triangulation, image, inside_indices


def interpolate_to_map_triangulate(
    plot_array,
    component,
    cell_size=0.002,
    n_padding_cells=10,
    nearest_neighbors=7,
    interp_pow=4,
):
    """

    `plot_array` must have key words:

        - **latitude**: latitude in decimal degrees of measured points
        - **longitude**: longitude in decimal degrees of measured points
        - values should have the proper name of the input component. For
          example if you are plotting the resistivity of the xy component then
          the keyword should be 'res_xy'

    :param plot_x: desired regular x locations to interpolate to
    :type plot_x: np.ndarray
    :param plot_y: desired regular x locations to interpolate to
    :type plot_y: np.ndarray
    :param plot_array: structured array, see above
    :type plot_array: np.ndarray
    :param component: component or keyword of the plot_array to plot
    :type component: string
    :param cell_size: size of cells , defaults to 0.002
    :type cell_size: TYPE, optional
    :param n_padding_cells: DESCRIPTION, defaults to 10
    :type n_padding_cells: TYPE, optional
    :param nearest_neighbors: DESCRIPTION, defaults to 7
    :type nearest_neighbors: TYPE, optional
    :param interp_pow: DESCRIPTION, defaults to 4
    :type interp_pow: TYPE, optional
    :return: DESCRIPTION
    :rtype: TYPE

    """

    # add padding to the locations
    ds = cell_size * n_padding_cells

    padded_x = plot_array["longitude"].copy()
    padded_y = plot_array["latitude"].copy()

    padded_x[np.argmin(padded_x)] -= ds
    padded_x[np.argmax(padded_x)] += ds
    padded_y[np.argmin(padded_y)] -= ds
    padded_y[np.argmax(padded_y)] += ds

    plot_x, plot_y = get_plot_xy(plot_array, cell_size, n_padding_cells)

    triangulation, image, inside_indices = triangulate_interpolation(
        plot_array["longitude"],
        plot_array["latitude"],
        plot_array[component],
        padded_x,
        padded_y,
        plot_x,
        plot_y,
        nearest_neighbors=nearest_neighbors,
        interp_pow=interp_pow,
    )

    if "res" in component:
        image = np.log10(image)

    return triangulation, image, inside_indices


def interpolate_to_map(
    plot_array,
    component,
    cell_size=0.002,
    n_padding_cells=10,
    interpolation_method="delaunay",
    interpolation_power=5,
    nearest_neighbors=7,
):
    """

    :param plot_array: DESCRIPTION
    :type plot_array: TYPE
    :param component: DESCRIPTION
    :type component: TYPE
    :param cell_size: DESCRIPTION, defaults to .002
    :type cell_size: TYPE, optional
    :param n_padding_cells: DESCRIPTION, defaults to 10
    :type n_padding_cells: TYPE, optional
    :param interpolation_method: DESCRIPTION, defaults to "delaunay"
    :type interpolation_method: TYPE, optional
    :return: DESCRIPTION
    :rtype: TYPE

    """

    if interpolation_method in ["nearest", "linear", "cubic"]:
        return interpolate_to_map_griddata(
            plot_array,
            component,
            cell_size=cell_size,
            n_padding_cells=n_padding_cells,
            interpolation_method=interpolation_method,
        )

    elif interpolation_method in [
        "fancy",
        "delaunay",
        "triangulate",
    ]:
        return interpolate_to_map_triangulate(
            plot_array,
            component,
            cell_size=cell_size,
            n_padding_cells=n_padding_cells,
            nearest_neighbors=nearest_neighbors,
            interp_pow=interpolation_power,
        )
