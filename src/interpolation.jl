using NCDatasets
using ProgressMeter
using Serialization
include("InterpolationData.jl")
using .InterpolationData

const PCO2_FILE = "daily.nc"
const SPATIAL_VARIATION_FILE = "pco2_spatial_variation.jldata"
const SEA_FILE = "sea.nc"

function run()
	
	#####################################
	## LOAD DATA
	loadprogress::Progress = Progress(6, 1, "Loading data")

	# pCO2 values
	local lons::Array{Union{Missing, Float32},1} = Dataset(PCO2_FILE)["longitude"][:]
	next!(loadprogress)
	local lats::Array{Union{Missing, Float32},1} = Dataset(PCO2_FILE)["latitude"][:]
	next!(loadprogress)
	local times::Array{Union{Missing, Float32},1} = Dataset(PCO2_FILE)["time"][:]
	next!(loadprogress)
	local pco2::Array{Union{Missing, Float64}, 3} = Dataset(PCO2_FILE)["pCO2"][:,:,:]
	next!(loadprogress)
	local uncertainty::Array{Union{Missing, Float64}, 3} = Dataset(PCO2_FILE)["uncertainty"][:,:,:]
	next!(loadprogress)

	# Spatial variation
	local inchan::IOStream = open(SPATIAL_VARIATION_FILE, "r")
	local spatialvariation::Array{Float64, 4} = deserialize(inchan)
	close(inchan)
	next!(loadprogress)

	# Sea mask
	local seamask::Array{Int64, 2} = convert.(Int64, Dataset(SEA_FILE)["SEA"][:,:])
	finish!(loadprogress)

	######################################
	## SET UP DATA STRUCTURES
	local interpolationbase = makeinterpolationbase(length(lons), length(lats), length(times), seamask)


	print("\n")
end

@time run()
