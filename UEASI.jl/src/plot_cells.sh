#!/bin/bash

rm -r plots
mkdir plots

export JULIA_PROJECT=".."
#julia -p 8 plot_cells.jl
julia plot_cells.jl
