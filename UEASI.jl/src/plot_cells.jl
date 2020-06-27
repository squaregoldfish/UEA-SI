using Distributed

@everywhere include("InterpolationData.jl")
@everywhere using .InterpolationData

using NCDatasets
using ProgressMeter

const SEA_FILE = "sea.nc"

function run()

	#####################################
	## LOAD DATA
	local lons::Vector{Float32} = Dataset(SEA_FILE)["LON"][:]
	local lats::Vector{Float32} = Dataset(SEA_FILE)["LAT"][:]

	cells::Vector{Cell} = makecells(length(lons), length(lats))

	@showprogress 1 pmap(x -> plotcell(x, "plots"), cells)

  #cell::Cell = cells[25]
  #plotcell(cell, "plots")
end

run()
