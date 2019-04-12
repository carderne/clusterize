# clusterize

Tool for converting raster population data into polygon settlement clusters, primarily for use with [openelec](https://github.com/carderne/openelec).

## Usage
Use `./run.py` with one of the three sub-commads:
1. `prep` prepares a raster for clusterizing
2. `make` convert a raster into polygon clusters
3. `feat` adds features from other data sources

Example:

    ```
    # Downsample to 0.001 degrees and remove average values below 1
    ./run.py prep pop_in.tif pop_out.tif -s 0.001 -f 1

    # Create clusters using radiusneighbors method, radius 3
    # and buffering by 100 metres before merging
    ./run.py make -m radius -r 3 -b 100 pop_prepped.tif clusters.gpkg

    # Add features as using a config file
    # Example in features.yml
    ./run.py feat -c features.yml clusters.gpkg clusters_out.gpkg

    ```

## Requirements
- `pyyaml`
- `numpy`
- `pandas`
- `scipy`
- `scikit-learn`
- `rasterio`
- `rasterstats`
- `shapely`
- `geopandas`

## Installation
Clone from GitHub, install requirements into a virtual environment, and use the `./run.py` script as described above.

