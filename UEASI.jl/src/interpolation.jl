# Indicates whether the processing data should be initialised from scratch
# For use during development. In production, this should be true
const __INIT_DATA__ = false

# Indicates whether the program should stop after initialising the data
# In production this should be false.
# If __INIT_DATA__ = false and this = true, nothing will happen
const __INIT_ONLY__ = false

using Distributed

@everywhere include("InterpolationData.jl")
@everywhere using .InterpolationData

using NCDatasets
using DataStructures
using ProgressMeter
using Serialization
using DelimitedFiles
using SharedArrays
using Missings

const FCO2_FILE = "daily.nc"
const SPATIAL_VARIATION_FILE = "fco2_spatial_variation.jldata"
const SEA_FILE = "sea.nc"
const SPATIAL_ACFS_FILE = "mean_directional_acfs.nc"
const TEMPORAL_ACFS_FILE = "mean_temporal_acf.csv"

function run()

	#####################################
	## LOAD DATA
	local loadcount::UInt8 = 0
	if __INIT_DATA__
		loadcount = 9
	else
		loadcount = 7
	end
	local loadprogress::Progress = Progress(loadcount, 0.001, "Loading data...")

	# fCO2 values
	local lons::Array{Union{Missing, Float32},1} = Dataset(FCO2_FILE)["longitude"][:]
	next!(loadprogress)
	local lats::Array{Union{Missing, Float32},1} = Dataset(FCO2_FILE)["latitude"][:]
	next!(loadprogress)

	if __INIT_DATA__
		local times::Array{Union{Missing, Float32},1} = Dataset(FCO2_FILE)["time"][:]
		next!(loadprogress)
		local fco2::Array{Union{Missing, Float64}, 3} = Dataset(FCO2_FILE)["fCO2"][:,:,:]
		next!(loadprogress)
		local uncertainty::Array{Union{Missing, Float64}, 3} = Dataset(FCO2_FILE)["uncertainty"][:,:,:]
		next!(loadprogress)
	end

	# Sea mask
	local seamask::SharedArray{UInt8, 2} = SharedArray(convert.(UInt8, Dataset(SEA_FILE)["SEA"][:,:]))
	next!(loadprogress)

	# Spatial variation
	local inchan::IOStream = open(SPATIAL_VARIATION_FILE, "r")
	local spatialvariation::SharedArray{Float64, 4} = SharedArray(deserialize(inchan))
	close(inchan)
	next!(loadprogress)

	# Autocorrelation data
	local spatialacfs::SharedArray{Float64, 4} = SharedArray(Float64.(collect(Missings.replace(Dataset(SPATIAL_ACFS_FILE)["mean_directional_acfs"][:,:,:,:], NaN))))
	next!(loadprogress)

	local temporalacf::SharedArray{Float64, 1} = SharedArray(readdlm(TEMPORAL_ACFS_FILE, ',', Float64, '\n')[:,2])

	finish!(loadprogress)

	######################################
	## SET UP DATA STRUCTURES

	local cells::Vector{Cell} = []
	if __INIT_DATA__
		cells = makecells(length(lons), length(lats), length(times), seamask, fco2, uncertainty)
	else
		cells = makecells(length(lons), length(lats))
	end

	######################################
	## REMOVE UNNEEDED DATA
	if __INIT_DATA__
		fco2 = zeros(1, 1, 1)
		uncertainty = zeros(1, 1, 1)
	end

	if __INIT_ONLY__
		exit()
	end

	######################################
	# NOW THE PROCESSING
	#
	# Repeatedly process all cells until the number of finished cells stablises
	local finishedarray::Vector{Int8} = zeros(length(cells))
	local lastfinishedcount::Int64 = 0
	local interpolationstep::UInt8 = 0

  while lastfinishedcount == 0 || lastfinishedcount != sum(finishedarray .≥ 1)
		lastfinishedcount = sum(finishedarray .≥ 1)
		interpolationstep += 1
  	finishedarray = @showprogress 1 "Interpolating cells..." pmap((x, y) ->
  		interpolatecell(x, y, interpolationstep, temporalacf, spatialacfs, spatialvariation, seamask), cells, finishedarray)

		println("Finished cells: $(sum(finishedarray .≥ 1))")
	end

	println("Final finished count: $lastfinishedcount")

	writefinished(finishedarray)
end

function writefinished(finished::Vector{Int8})

	lons = Dataset(SEA_FILE)["LON"][:]
	lats = Dataset(SEA_FILE)["LAT"][:]

	nc = Dataset("finished_steps.nc", "c")
	defDim(nc, "lon", length(lons))
	defDim(nc, "lat", length(lats))

	lonVar = defVar(nc, "lon", Float32, ("lon",), attrib = OrderedDict(
    "units" => "degrees_east"))
	lonVar[:] = lons

	latVar = defVar(nc, "lat", Float32, ("lat",), attrib = OrderedDict(
    "units" => "degrees_north"))
	latVar[:] = lats

	var = defVar(nc, "step", Int8, ("lat", "lon"), fillvalue = -1)
	var[:,:] = reshape(finished, (length(lats), length(lons)))

	close(nc)
end

run()
