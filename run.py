#!/usr/bin/env python

"""
Convert raster/gridded population data into discrete
polygon population clusters.

Must be run with Python 3.6+ with the following packages available:
    numpy, pandas, geopandas, sklearn, rasterio, shapely
"""

import os
from pathlib import Path

import click
import yaml
import geopandas as gpd

from clusterize import (
    extract_points,
    neighbors,
    make_clusters,
    make_geometry,
    read_raster,
)
from clusterize.features import add_raster_layer, add_vector_layer, fix_column

script_dir = Path(os.path.dirname(__file__))
cfg_default = script_dir / "features.yml"


@click.group(help="Utility for creating population clusters.")
def cli():
    pass


@cli.command()
@click.argument("raster_in")
@click.argument("raster_out")
@click.option(
    "--res", default=0.1, type=click.FLOAT, help="New resolution in CRS units"
)
@click.option(
    "--min_val", default=5, type=click.INT, help="Minimum value (after resampling)"
)
def prep(raster_in, raster_out, res, min_val):
    """Resample and threshold a population raster."""

    command = (
        f"gdal_translate -ot Byte -a_nodata none {raster_in} /vsistdout/ | "
        f"gdalwarp -tr {res} {res} -r average /vsistdin/ /vsistdout/ | "
        f"gdal_calc.py -A /vsistdin/ --outfile={raster_out} --calc='A>{min_val}' --NoDataValue=0"
    )
    os.system(command)


@cli.command()
@click.argument("raster_in")
@click.argument("gpkg_out")
@click.option(
    "--min-val",
    default=1,
    type=click.INT,
    help="Minimum raster value to be considered a population site",
)
@click.option(
    "--method",
    default="radius",
    type=click.STRING,
    help="Options are 'radius' or 'k' for the two different sklean methods",
)
@click.option(
    "--radius", default=3, type=click.INT, help="Radius in cells for radius method"
)
@click.option(
    "--n-neighbors",
    default=10,
    type=click.INT,
    help="Number of neighbors for kneighbors method",
)
@click.option(
    "--min-neighbors",
    default=5,
    type=click.INT,
    help="Discard grouping below this cutoff",
)
@click.option(
    "--max-dist", default=5, type=click.INT, help="Discard neighbors above this cutoff"
)
@click.option(
    "--buffer",
    default=100,
    type=click.INT,
    help="Amount in metres by which to buffer clusters before merging",
)
def make(
    raster_in,
    gpkg_out,
    min_val,
    method,
    radius,
    n_neighbors,
    min_neighbors,
    max_dist,
    buffer,
):
    """Create clusters."""

    arr, affine, crs = read_raster(raster_in)
    X = extract_points(arr, min_val)
    groups = neighbors(
        X,
        method=method,
        radius=radius,
        n_neighbors=n_neighbors,
        min_neighbors=min_neighbors,
        max_dist=max_dist,
    )
    clu = make_clusters(groups)
    gdf = make_geometry(clu, X, affine, crs, buffer_amount=buffer)
    gdf.to_file(gpkg_out, driver="GPKG")


@cli.command()
@click.argument("clusters_in")
@click.argument("clusters_out")
@click.option(
    "--config", default=cfg_default, type=click.STRING, help="Path to config file"
)
def feat(clusters_in, clusters_out, config):
    """Add features."""

    cfg_in = Path(config)
    with open(cfg_in, "r") as ymlfile:
        cfg = yaml.safe_load(ymlfile)

    clusters = gpd.read_file(clusters_in)

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
    clusters.to_file(clusters_out, driver="GeoJSON")


if __name__ == "__main__":
    cli()
