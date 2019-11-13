#!/usr/bin/env python

"""
Convert raster/gridded population data into discrete
polygon population clusters.

Must be run with Python 3.6+ with the following packages available:
    numpy, pandas, geopandas, sklearn, rasterio, shapely
"""

from statistics import median
from argparse import ArgumentParser

import numpy as np
import pandas as pd
import geopandas as gpd
from sklearn.neighbors import NearestNeighbors
import rasterio
from rasterio.transform import xy
import shapely.wkt
from shapely.geometry import MultiPoint


def read_raster(rast_in):
    """Read file and return numpy array, affine and crs."""

    rast_rd = rasterio.open(rast_in)
    arr = rast_rd.read(1)
    affine = rast_rd.transform
    crs = rast_rd.crs
    nodata = rast_rd.nodata
    arr[arr == nodata] = 0

    return arr, affine, crs


def extract_points(arr, min_val=1):
    """Extract points with value==true_val into a list of coordinates."""

    rr = np.where(arr >= min_val)
    X = [(x, y) for x, y in zip(rr[0], rr[1])]
    print("Num points:", len(X))
    print("Points sample:", X[:4])

    return X


def neighbors(
    X,
    method="radius",
    radius=3,
    n_neighbors=10,
    algorithm="auto",
    min_neighbors=5,
    max_dist=5,
):
    """
    Run sklearn nearest neighbor algorithm.

    Parameters
    ----------
    X : array-like, shape = [n_samples, n_features]
        List of points.
    method : string, default "radius"
        Can be either "radius" or "k" for the radius_neighbors or kneighbors
        algorithms from sklearn.neighbors.NearestNeighbors.
    radius : int, default 3
        Radius to use in radius_neighbors method.
    n_neighbors : int, default 10
        n_neighbors to use in kneighbors method.
    algorithm : string, default "auto"
        Algorithm to use, can be "ball_tree", "kd_tree", "brute" or "auto".
    min_neighbors : int, default 5
        Minimum number of neighbors for a grouping to be kept.
    max_dist : int, default 5
        Maximum distance (in cells) for a certain neighbor to be kept.

    Returns
    -------
    groups : list of sets
        List of sets of nearest neighbor groups.
    """

    nbrs = NearestNeighbors(algorithm=algorithm).fit(X)

    if method == "radius":
        distances, indices = nbrs.radius_neighbors(X, radius=radius)
    elif method == "k":
        distances, indices = nbrs.kneighbors(X, n_neighbors=n_neighbors)
    else:
        raise ValueError("method must be either radius or k")

    # Extract all values into a list of sets where the two conditions are met
    groups = []
    for ids, dists in zip(indices, distances):
        keepers = set([i for i, d in zip(ids, dists) if d < max_dist])
        if len(keepers) >= min_neighbors:
            groups.append(keepers)
    print("Num groups:", len(groups))

    return groups


def make_clusters(groups):
    """Transform groups into clusters containing all shared points."""

    # Combine all values that share nearest neighbors into clusters
    clu = []
    for gr in groups:
        found = False
        for i, c in enumerate(clu):
            if len(gr & c) > 0:  # If there are shared elements
                clu[i] |= gr  # Add the elements of k to that cluster
                found = True
                break
        if not found:
            clu.append(gr)

    # There will be adjacent clusters that aren't merged,
    # but the geometry merge later will fix that and seems more robust.

    print("Num clusters:", len(clu))
    print("Median size:", median(len(c) for c in clu))

    return clu


def raster_out(file_out, clu, X, arr, affine, crs):
    """Save  clustered cells as raster."""

    # Assign cluster_num to each cluster and burn into raster
    rast_out = np.zeros_like(arr, dtype="int32")
    cluster_num = 1  # avoid 0 because nodatafilename
    for c in clu:
        for loc in c:
            rast_out[X[loc]] = cluster_num

    # Export values
    filtered_out = rasterio.open(
        file_out,
        "w",
        driver="GTiff",
        height=rast_out.shape[0],
        width=rast_out.shape[1],
        count=1,
        dtype=rast_out.dtype,
        crs=crs,
        transform=affine,
        nodata=0,
    )
    filtered_out.write(rast_out, 1)
    filtered_out.close()


def merge_overlap(gdf):
    """Merge overlapping geometries."""

    gdf["same"] = 1
    gdf = gdf.dissolve(by="same")
    gdf = gdf.explode()
    gdf = gdf.reset_index()
    gdf = gdf.drop(columns=["same", "level_1"])

    return gdf


def make_geometry(clu, X, affine, crs, buffer_amount=100):
    """
    Convert clusters cells into contigious convex hulls.
    These are then buffered and overlapping polygons are merged.

    Parameters
    ----------
    clu : list of sets
        The array of clusters with indices.
    X : array-like, shape = [n_samples, n_features]
        The list of points with coordinates.
    affine : affine.Affine
        Raster affine transformation.
    crs : CRS
        Coordinate reference system.
    buffer_amount : int, default 100
        Amount in metres by which to buffer polygons before merging.

    Returns
    -------
    clusters : GeoDataFrame
        The geometry-fied clusters.
    """

    # This is the Africa Albers Equal Area Conic EPSG: 102022
    EPSG102022 = """+proj=aea
                    +lat_1=20 +lat_2=-23 +lat_0=0 +lon_0=25
                    +x_0=0 +y_0=0
                    +ellps=WGS84 +datum=WGS84 +units=m +no_defs"""
    EPSG3857 = {"init": "epsg:3857"}

    clusters = []
    for c in clu:
        coords = [X[loc] for loc in c]
        coords_real = [xy(affine, loc[0], loc[1]) for loc in coords]
        m = MultiPoint(coords_real)
        clusters.append(m.wkt)

    gdf = pd.DataFrame(clusters)
    geometry = gdf[0].map(shapely.wkt.loads)
    gdf = gdf.drop(0, axis=1)
    gdf = gpd.GeoDataFrame(gdf, crs=crs, geometry=geometry)

    gdf["geometry"] = gdf.geometry.convex_hull
    # gdf = gdf.to_crs(EPSG3857)
    buffer_amount /= 1e5
    gdf["geometry"] = gdf.geometry.buffer(buffer_amount)
    gdf = gdf.to_crs(crs)

    gdf = merge_overlap(gdf)
    gdf["geometry"] = gdf.geometry.convex_hull
    gdf = merge_overlap(gdf)

    print("Number of clusters:", len(gdf))

    return gdf
