# Indicates whether the processing data should be initialised from scratch
# For use during development. In production, this should be true
const __INIT_DATA__ = false

# Indicates whether the program should stop after initialising the data
# In production this should be false.
# If __INIT_DATA__ = false and this = true, nothing will happen
const __INIT_ONLY__ = false

using Distributed
using NCDatasets
using ProgressMeter
using Serialization
using DelimitedFiles

@everywhere include("InterpolationData.jl")
@everywhere using .InterpolationData

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
	local loadprogress::Progress = Progress(loadcount, 1, "Loading data...")

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
	local seamask::Array{UInt8, 2} = convert.(UInt8, Dataset(SEA_FILE)["SEA"][:,:])
	next!(loadprogress)

	# Spatial variation
	local inchan::IOStream = open(SPATIAL_VARIATION_FILE, "r")
	local spatialvariation::Array{Union{Missing, Float64}, 4} = deserialize(inchan)
	close(inchan)
	replace!(spatialvariation, -1e35=>missing)
	next!(loadprogress)

	# Autocorrelation data
	local spatialacfs::Array{Union{Missing, Float64}, 4} = Dataset(SPATIAL_ACFS_FILE)["mean_directional_acfs"][:,:,:,:]
	next!(loadprogress)

	local temporalacf::Array{Float64, 1} = readdlm(TEMPORAL_ACFS_FILE, ',', Float64, '\n')[:,2]

	# Initialise InterpolationData module, incl. loading bathymetry
	InterpolationData.init(length(lons), length(lats))

	finish!(loadprogress)

	######################################
	## SET UP DATA STRUCTURES

	local cells::Array{Cell, 1} = []
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
#	local finishedarray::Array{Bool, 1} = falses(length(cells))
#	local lastfinishedcount::Int64 = 0
#	local interpolationstep::Int8 = 0

	#while lastfinishedcount == 0 || lastfinishedcount != sum(finishedarray .== true)
	#	lastfinishedcount = sum(finishedarray .== true)
	#	interpolationstep += 1
	#	finishedarray = @showprogress 1 "Interpolating cells..." pmap(x -> interpolatecell(x, interpolationstep), cells)
	#	println("Finished cells: $(sum(finishedarray .== true))")
	#end

	testcell::Cell = (lon=69, lat=25)
	interpolatecell(testcell, convert(UInt8, 1), temporalacf, spatialacfs, spatialvariation, seamask)

#	println("Final finished count: $lastfinishedcount")

	print("\n")
end

@time run()