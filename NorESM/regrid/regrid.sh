#!/bin/bash

# Directories
INDIR="in"
OUTDIR="out"
TEMP="/tmp"

# The remapping weights file
MAP_FILE="map_tnx1v1_to_woa09_aave_20120501.nc"

for fullpath in $INDIR/*.hd.*
do
	file=$(basename $fullpath)

	# The model output duplicates the last latitude (index 385),
	# so we need to remove it
	sed -e 's!%%INFILE%%!'"$INDIR/$file"'!' 384.jnl > $TEMP/384.jnl
	pyferret -script /tmp/384.jnl

	# Convert to a 1x1° grid
	ncremap -m $MAP_FILE -R '--rgr lat_nm=Y1_384 --rgr lon_nm=X' -i /tmp/384.nc -o /tmp/1deg.nc

	# The output of remap has values and fractions in separate variables
	# Here we combine them
	sed -e 's!%%OUTFILE%%!'"$OUTDIR/$file"'!' frac_b.jnl > $TEMP/frac_b.jnl
	pyferret -script /tmp/frac_b.jnl

	# Make the variables in the output file more sensible
	ncrename -O -d Y1_384,latitude -v Y1_384,latitude -d X,longitude -v X,longitude -v SST_FIXED,sst /tmp/sst_fixed.nc $OUTDIR/$file

	# Cleanup
	rm $TEMP/384.jnl
	rm $TEMP/384.nc
	rm $TEMP/1deg.nc
	rm $TEMP/frac_b.jnl
	rm $TEMP/sst_fixed.nc
done