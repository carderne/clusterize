#!/usr/bin/env python

"""
Attempt at clustering population raster using KMeans.
Doesn't work -- I don't think KMeans is designed for thousands of clusters.
"""

import argparse

import numpy as np
import rasterio
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans


def clusterize(raster_in, raster_out, n_clusters):
    rast_rd = rasterio.open(raster_in)
    rast = rast_rd.read(1)

    rast = rast.astype("int32")

    rr = np.where(rast == 1)
    points = [(x, y) for x, y in zip(rr[0], rr[1])]

    y_pred = KMeans(n_clusters=n_clusters).fit_predict(points)

    for loc, y in zip(points, y_pred):
        rast[loc] = y

    rast_write = rasterio.open(
        raster_out,
        "w",
        driver="GTiff",
        height=rast.shape[0],
        width=rast.shape[1],
        count=1,
        dtype=rast.dtype,
        crs=rast_rd.crs,
        transform=rast_rd.transform,
    )
    rast_write.write(rast, 1)
    rast_write.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("raster_in")
    parser.add_argument("raster_out")
    parser.add_argument("n_clusters")
    args = parser.parse_args()

    clusterize(args.raster_in, args.raster_out, int(args.n_clusters))

