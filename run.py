#!/usr/bin/env python

"""
Convert raster/gridded population data into discrete
polygon population clusters.

Must be run with Python 3.6+ with the following packages available:
    numpy, pandas, geopandas, matplotlib, sklearn, rasterio, shapely
"""

import os
import sys
from argparse import ArgumentParser
from pathlib import Path

import yaml
import geopandas as gpd

from clusterize import extract_points, neighbors, make_clusters, make_geometry
from parser_feat import add_raster_layer, add_vector_layer, fix_column, read_raster

script_dir = Path(os.path.dirname(__file__))
cfg_default = script_dir / "features.yml"

EPSG102022 = "+proj=aea +lat_1=20 +lat_2=-23 +lat_0=0 +lon_0=25 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"


def prep(args):
    pass

# Downsample and filter a provided population raster.
# Arg 1: Input raster name
# Arg 2: New resolution in degrees
# Arg 3: Minimum value after resampling
# Arg 4: Output file name

# gdal_translate -ot Byte -a_nodata none $1 /vsistdout/ | gdalwarp -tr $2 $2 -r average /vsistdin/ /vsistdout/ | gdal_calc.py -A /vsistdin/ --outfile=$4 --calc="A>$3" --NoDataValue=0


def make(args):
    arr, affine, crs = read_raster(args.raster_in)
    X = extract_points(arr)
    groups = neighbors(
        X,
        method=args.method,
        radius=args.radius,
        n_neighbors=args.n_neighbors,
        min_neighbors=args.min_neighbors,
        max_dist=args.max_dist,
    )
    clu = make_clusters(groups)
    gdf = make_geometry(clu, X, affine, crs, buffer_amount=args.buffer)
    gdf.to_file(args.gpkg_out, driver="GPKG")


def feat(args):
    cfg_in = Path(args.config)
    with open(cfg_in, "r") as ymlfile:
        cfg = yaml.safe_load(ymlfile)

    clusters = gpd.read_file(args.clusters_in)

    for f in cfg:
        print(f["name"])
        if f["type"] == "raster":
            clusters = add_raster_layer(
                clusters=clusters,
                raster=Path(f["file"]).expanduser(),
                operation=f["operation"],
                col_name=f["name"],
                crs=f["crs"],
                decimals=f["decimals"],
            )

        elif f["type"] == "vector":
            clusters = add_vector_layer(
                clusters=clusters,
                vector=Path(f["file"]).expanduser(),
                operation=f["operation"],
                col_name=f["name"],
                shape="shape",
                affine="affine",
                raster_crs="raster_crs",
                decimals=f["decimals"],
            )

        else:
            raise ValueError("Only 'raster' or 'vector' supported for 'type'.")

        if "fix" in f:
            clusters = fix_column(
                clusters=clusters,
                col_name=f["name"],
                factor=f["factor"],
                minimum=f["minimum"],
                maximum=f["maximum"],
                no_value=f["no_value"],
                per_capita=f["per_capita"],
            )

    clusters = clusters.to_crs(EPSG102022)
    clusters["area"] = clusters.geometry.area
    clusters["x"] = clusters.geometry.centroid.x
    clusters["y"] = clusters.geometry.centroid.y
    clusters = clusters.to_crs(epsg=4326)

    clusters = clusters.sort_values(by="area", ascending=False)
    clusters["fid"] = clusters.index
    clusters = clusters.dropna(axis=0, subset=["geometry"])
    clusters = clusters.fillna(0)  # There were NaN in NTL (at least)
    clusters = clusters.loc[clusters["area"] > 0]

    clusters.geometry = clusters.simplify(tolerance=0.001, preserve_topology=False)

    clusters.to_file(args.clusters_out, driver="GeoJSON")


# aoi = gpd.read_file("/home/chris/Documents/GIS/gadm.gpkg")
# aoi = aoi.loc[aoi["NAME_0"] == "Tanzania"]
# clipped, affine, crs = clip_raster(gis/"GHS_POP_250.tif", aoi)


if __name__ == "__main__":
    parser = ArgumentParser()
    subparsers = parser.add_subparsers(dest="tool", title="Subcommands")
    parser_prep = subparsers.add_parser("prep", help="Prep raster for clustering")
    parser_prep.set_defaults(func=prep)

    parser_make = subparsers.add_parser("make", help="Create clusters")
    parser_make.add_argument("raster_in", help="Input raster file to be processed")
    parser_make.add_argument("gpkg_out", help="Filename for output GPKG")
    parser_make.add_argument(
        "--min_val",
        type=int,
        default=1,
        help="Minimum raster value to be considered a population site",
    )
    parser_make.add_argument(
        "-m",
        "--method",
        default="radius",
        help="Options are 'radius' or 'k' for the two different sklean methods",
    )
    parser_make.add_argument(
        "-r", "--radius", type=int, default=3, help="Radius in cells for radius method"
    )
    parser_make.add_argument(
        "--n_neighbors",
        type=int,
        default=10,
        help="Number of neighbors for kneighbors method",
    )
    parser_make.add_argument(
        "--min_neighbors",
        type=int,
        default=5,
        help="Discard grouping below this cutoff",
    )
    parser_make.add_argument(
        "--max_dist", type=int, default=5, help="Discard neighbors above this cutoff"
    )
    parser_make.add_argument(
        "-b",
        "--buffer",
        type=int,
        default=100,
        help="Amount in metres by which to buffer clusters before merging",
    )
    parser_make.set_defaults(func=make)

    parser_feat = subparsers.add_parser("feat", help="Add features")
    parser_feat.add_argument(
        "-c", "--config", default=cfg_default, help="Path to config file"
    )
    parser_feat.add_argument("clusters_in")
    parser_feat.add_argument("clusters_out")
    parser_feat.set_defaults(func=feat)
    args = parser.parse_args()
    args.func(args)
