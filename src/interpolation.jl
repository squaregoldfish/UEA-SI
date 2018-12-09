using Distributed
using NCDatasets
using ProgressMeter
using Serialization

@everywhere include("InterpolationData.jl")
@everywhere using .InterpolationData

const FCO2_FILE = "daily.nc"
const SPATIAL_VARIATION_FILE = "fco2_spatial_variation.jldata"
const SEA_FILE = "sea.nc"

function run()
	
	#####################################
	## LOAD DATA
	loadprogress::Progress = Progress(6, 1, "Loading data")

	# fCO2 values
	local lons::Array{Union{Missing, Float32},1} = Dataset(FCO2_FILE)["longitude"][:]
	next!(loadprogress)
	local lats::Array{Union{Missing, Float32},1} = Dataset(FCO2_FILE)["latitude"][:]
	next!(loadprogress)
	local times::Array{Union{Missing, Float32},1} = Dataset(FCO2_FILE)["time"][:]
	next!(loadprogress)
	local fco2::Array{Union{Missing, Float64}, 3} = Dataset(FCO2_FILE)["fCO2"][:,:,:]
	next!(loadprogress)
	local uncertainty::Array{Union{Missing, Float64}, 3} = Dataset(FCO2_FILE)["uncertainty"][:,:,:]
	next!(loadprogress)

	# Spatial variation
	local inchan::IOStream = open(SPATIAL_VARIATION_FILE, "r")
	local spatialvariation::Array{Float64, 4} = deserialize(inchan)
	close(inchan)
	next!(loadprogress)

	# Sea mask
	local seamask::Array{Int8, 2} = convert.(Int8, Dataset(SEA_FILE)["SEA"][:,:])
	finish!(loadprogress)

	######################################
	## SET UP DATA STRUCTURES
	local cells::Array{Cell, 1} = makecells(length(lons), length(lats), length(times), seamask, fco2, uncertainty)

	print("\n")
end

@time run()
