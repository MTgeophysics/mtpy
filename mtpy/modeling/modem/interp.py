import numpy as np
from scipy.spatial import cKDTree as KDTree

def IDW(source_points, source_values, query_points, k=6, p=5):
    tree = KDTree(source_points.points, k)
    distances, indices = tree.query(query_points.points, k=k, p=p)
    inv_dist = 1. / distances
    weights = inv_dist / inv_dist.sum(axis=1)[:, np.newaxis]
    return (weights * source_values[indices]).sum(axis=1)


def nearest(source_points, source_values, query_points):
    tree = KDTree(source_points.points, 1)
    _, indices = tree.query(query_points.points, k=1, p=2)
    return source_values[indices]
