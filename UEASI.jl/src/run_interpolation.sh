#!/bin/bash

echo "Deleting..."
rm -r interpolation_data
echo "Copying..."
cp -r interpolation_data.new interpolation_data

echo "Running..."
#julia -p 8 interpolation.jl
julia --project=.. interpolation.jl
