import numpy as np
from scipy.spatial import cKDTree as KDTree
from scipy.interpolate import griddata

def IDW(source_points, source_values, query_points, k=6, p=5):
    tree = KDTree(source_points.points, k)
    distances, indices = tree.query(query_points.points, k=k)
    inv_dist = 1. / np.power(distances, p)
    weights = inv_dist / inv_dist.sum(axis=1)[:, np.newaxis]
    return (weights * source_values[indices]).sum(axis=1)


def nearest(source_points, source_values, query_points):
    tree = KDTree(source_points.points, 1)
    _, indices = tree.query(query_points.points, k=1, p=2)
    return source_values[indices]


def linear(source_points, source_values, query_grid):
    assert query_grid.shape[-1] == 3
    grids = (query_grid[..., 0], query_grid[..., 1], query_grid[..., 2])
    print(grids[0].shape)
    return griddata(source_points.points, source_values, grids)
