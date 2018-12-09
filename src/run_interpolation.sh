#!/bin/bash

read -p "Copy? " copy

if [ $copy == "y" ]
then 
    echo "Deleting..."
    rm -r interpolation_data
    echo "Copying..."
    cp -r interpolation_data.new interpolation_data
fi

echo "Running..."
julia -p 8 interpolation.jl
