using NCDatasets
using ProgressMeter
using Serialization

const PCO2_FILE = "daily.nc"
const SPATIAL_VARIATION_FILE = "pco2_spatial_variation.jldata"
const SEA_FILE = "sea.nc"

function run()
	loadprogress::Progress = Progress(6, 1, "Loading data")

	# pCO2 values
	local lons::Array{Union{Missing, Float32},1} = Dataset(PCO2_FILE)["longitude"][:]
	next!(loadprogress)
	local lats::Array{Union{Missing, Float32},1} = Dataset(PCO2_FILE)["latitude"][:]
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
	next!(loadprogress)

	print(seamask)

	print("\n")
end

@time run()
