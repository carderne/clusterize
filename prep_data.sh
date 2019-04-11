#!/bin/bash

# Downsample and filter a provided population raster.
# Arg 1: Input raster name
# Arg 2: New resolution in degrees
# Arg 3: Minimum value after resampling
# Arg 4: Output file name

gdal_translate -ot Byte -a_nodata none $1 /vsistdout/ | gdalwarp -tr $2 $2 -r average /vsistdin/ /vsistdout/ | gdal_calc.py -A /vsistdin/ --outfile=$4 --calc="A>$3" --NoDataValue=0

