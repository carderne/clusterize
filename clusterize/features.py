"""
clusters module for openelec

Provides functions to read in a raster population dataset and convert to discrete
vector polgons, each with a set population value. Additionally calculate each
polygon's distance from a provided grid infrastructure vector.

Functions:

- clip_rasters
- add_raster_layer
- add_vector_layer
- fix_column
- save_clusters
"""

import json
from pathlib import Path

import numpy as np
from scipy import ndimage
import geopandas as gpd

import rasterio
from rasterio.mask import mask
from rasterio.features import shapes, rasterize
from rasterstats import zonal_stats


def clip_raster(raster, boundary, boundary_layer=None):
    """
    Clip the raster to the given administrative boundary.

    Parameters
    ----------
    raster: string, pathlib.Path or rasterio.io.DataSetReader
        Location of or already opened raster.
    boundary: string, pathlib.Path or geopandas.GeoDataFrame
        The poylgon by which to clip the raster.
    boundary_layer: string, optional
        For multi-layer files (like GeoPackage), specify the layer to be used.

    Returns
    -------
    tuple
        Three elements:
            clipped: numpy.ndarray
                Contents of clipped raster.
            affine: affine.Affine()
                Information for mapping pixel coordinates
                to a coordinate system.
            crs: dict
                Dict of the form {'init': 'epsg:4326'} defining the coordinate
                reference system of the raster.
    """

    if isinstance(raster, Path):
        raster = str(raster)
    if isinstance(raster, str):
        raster = rasterio.open(raster)

    crs = raster.crs

    if isinstance(boundary, Path):
        boundary = str(boundary)
    if isinstance(boundary, str):
        if ".gpkg" in boundary:
            driver = "GPKG"
        else:
            driver = None  # default to shapefile
            boundary_layer = ""  # because shapefiles have no layers

        boundary = gpd.read_file(boundary, layer=boundary_layer, driver=driver)

    boundary = boundary.to_crs(crs=raster.crs)
    coords = [json.loads(boundary.to_json())["features"][0]["geometry"]]

    # mask/clip the raster using rasterio.mask
    clipped, affine = mask(dataset=raster, shapes=coords, crop=True)

    return clipped, affine, crs


def add_raster_layer(
    clusters, raster, operation, col_name, affine=None, crs=None, decimals=2
):
    """
    The filter_merge_clusters() process loses the underlying raster values.
    So we need to use rasterstats.zonal_stats() to get it back.

    Parameters
    ----------
    clusters: geopandas.GeoDataFrame
        The processed clusters.
    raster: str, pathlib.Path or numpy.ndarray
        Either a path to the raster, or numpy.ndarray with the data.
    operation: str
        The operation to perform when extracting the raster data.
        Either 'sum', 'max', or 'mean'
    col_name: str
        Name of the column to add.
    affine: affine.Affine(), optional
        If a numpy ndarray is passed above, the affine is also needed.
    crs: proj.crs, optional
        Override raster's reported crs

    Returns
    -------
    clusters: geopandas.GeoDataFrame
        The processed clusters with new column.
    """

    if isinstance(raster, Path):
        raster = str(raster)
    if isinstance(raster, str):
        # rasterstats doesn't check for same CRS
        # Throws memory error if don't ensure they are same
        if not crs:
            crs = rasterio.open(raster).crs
        clusters_proj = clusters.to_crs(crs)
        stats = zonal_stats(clusters_proj, raster, stats=operation)

        clusters_proj[col_name] = [x[operation] for x in stats]

        clusters = clusters_proj.to_crs(clusters.crs)
        clusters[col_name] = clusters[col_name].round(decimals)

        return clusters

    else:
        raise NotImplementedError("Only implemented for path input.")


def add_vector_layer(
    clusters, vector, operation, col_name, shape, affine, raster_crs, decimals=2
):
    """
    Use a vector containing grid infrastructure to determine
    each cluster's distance from the grid.

    Parameters
    ----------
    clusters: geopandas.GeoDataFrame
        The processed clusters.
    vector: str, pathlib.Path or geopandas.GeoDataFrame
        Path to or already imported grid dataframe.
    operation: str
        Operation to perform in extracting vector data.
        Currently only 'distance' supported.
    shape: tuple
        Tuple of two integers representing the shape of the data
        for rasterizing grid. Sould match the clipped raster.
    affine: affine.Affine()
        As above, should match the clipped raster.

    Returns
    -------
    clusters: geopandas.GeoDataFrame
        The processed clusters with new column.
    """

    if isinstance(vector, Path):
        vector = str(vector)
    if isinstance(vector, str):
        vector = gpd.read_file(vector)

    vector = vector.to_crs(crs=raster_crs)
    clusters = clusters.to_crs(crs=raster_crs)

    if operation == "distance":
        vector = vector.loc[vector["geometry"].length > 0]

        grid_raster = rasterize(
            vector.geometry,
            out_shape=shape,
            fill=1,
            default_value=0,
            all_touched=True,
            transform=affine,
        )
        dist_raster = ndimage.distance_transform_edt(grid_raster) * affine[0]

        dists = zonal_stats(
            vectors=clusters, raster=dist_raster, affine=affine, stats="min"
        )
        clusters[col_name] = [x["min"] for x in dists]
        clusters[col_name] = clusters[col_name].round(decimals)
        clusters = clusters.to_crs(epsg=4326)

        return clusters

    else:
        raise NotImplementedError('Currently only "distance" is supported.')


def fix_column(
    clusters,
    col_name,
    factor=1,
    minimum=0,
    maximum=None,
    no_value=None,
    per_capita=False,
):
    """
    A number of operations to apply to a columns values to get desired output.

    Parameters
    ----------
    clusters : GeoDataFrame
        The clusters object.
    col_name : str
        The column to apply the operation to.
    factor : float, optional (default 1.)
        Factor by which to multiply the column vales.
    minimum : float, optional (default 0.)
        Apply a minimum threshold to the values.
    maximum : str, optional
        Currently only supported for 'largest'.
        Limits the values to double the value of the cluster with the highest
        population.
    no_value : str, optional
        Currently only supported for 'median'.
        Replaces NaN instances with the median value.
    per_capita : boolean, optional (default False.)
        Divide values by cluster population.

    Returns
    -------
    clusters : GeoDataFrame
        The 'fixed' clusters.
    """

    # multiply the column by a fixed factor
    if factor is not None and factor != 1:
        clusters[col_name] = clusters[col_name] * factor

    # remove negative values
    if minimum is not None:
        clusters.loc[clusters[col_name] < minimum, col_name] = minimum

    if per_capita:
        clusters[col_name] = clusters[col_name] / clusters["pop"]

    # apply a cutoff maximum value
    if maximum is not None:
        if maximum == "largest":
            limit = 2 * float(
                clusters.loc[
                    clusters["pop"] == clusters["pop"].max(), col_name
                ].tolist()[0]
            )
            clusters.loc[clusters[col_name] > limit, col_name] = limit

        else:
            raise NotImplementedError("maximum only implemented for largest.")

    # replace nan values
    if no_value is not None:
        if no_value == "median":
            replace = {col_name: clusters[col_name].median()}
            clusters = clusters.fillna(value=replace)

        else:
            raise NotImplementedError("no_value only implemented for median.")

    return clusters


def save_clusters(clusters, out_path):
    """
    Convert to EPSG:4326 and save to the specified file.
    clusters: geopandas.GeoDataFrame
        The processed clusters.
    out_path: str or pathlib.Path
        Where to save the clusters file.
    """

    if isinstance(out_path, Path):
        out_path = str(out_path)
    if ".gpkg" in out_path:
        driver = "GPKG"
    elif ".geojson" in out_path or ".json" in out_path:
        driver = "GeoJSON"
    else:
        driver = None

    clusters = clusters.to_crs(epsg=4326)
    clusters.to_file(out_path, driver=driver)
