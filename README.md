# clusterize

Tool for converting raster population data into polygon settlement clusters, primarily for use with [openelec](https://github.com/carderne/openelec).

Steps to use:
1. Use `prep_data.sh` with a population raster to downsample and filter low population pixels.
Example usage to downsample to 0.001 degrees and rsubsequently remove pixels with an average value below 1:
    ```
    ./prep_data.sh pop_in.tif 0.001 1 pop_prepped.tif
    ```

2. Use `clusterize.py` to create polygon clusters from this file.
Example usage to use the radiusneighbors method with a radius of 3, subsequently buffering by 100 metres and merging:
    ```
    ./clusterize.py -m radius -r 3 -b 100 pop_prepped.tif clusters.gpkg
    ```

3. Use `add_features.py` to add population, grid distance and whatever other raster and vector features to the clusters.
(Simple example still to come.)

