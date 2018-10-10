using NCDatasets
using Serialization

const INFILE = "daily.nc"
const OUTFILE = "pco2_spatial_variation.jldata"

function run()
	print("Loading input data...")
	lons = Dataset(INFILE)["longitude"][:]
	lats = Dataset(INFILE)["latitude"][:]
	pco2 = Dataset(INFILE)["pCO2"][:,:,:]

	lonsize::Int64 = length(lons)
	latsize::Int64 = length(lats)
	timesize::Int64 = size(pco2)[3]

	print("\033[1K\rInitialising data structures...")
	variationtotal::Array{Float64, 4} = zeros(lonsize, latsize, lonsize, latsize)
	variationcount::Array{Int64, 4} = zeros(lonsize, latsize, lonsize, latsize)

	for timeloop::Int64 in 1:1000
		print("\033[1K\r$timeloop / $timesize")
		for lonloop::Int64 in 1:lonsize
			for latloop::Int64 in 1:latsize
				currentvalue::Union{Missing, Float64} = pco2[lonloop, latloop, timeloop]
				if !ismissing(currentvalue)

	                # Search for other cells with values at the same time.
	                # Note that we include values from a 7 timesteps either side too,
	                # just to increase the amount of data we can use. It should also
	                # temper the most optimistic change values we see.
	                compare::Array{Union{Missing, Float64}, 3} = pco2[:, :, (timeloop - 7):(timeloop + 7)]
	                diffs::Array{Union{Missing, Float64}, 3} = similar(compare)
	                @. diffs = abs(compare - currentvalue)

	                for difflon in 1:lonsize
	                	for difflat in 1:latsize
	                		celldiffs = diffs[difflon, difflat, :]
	                		variationtotal[lonloop, latloop, difflon, difflat] = 
	                			variationtotal[lonloop, latloop, difflon, difflat] + sum(skipmissing(celldiffs))
	                		variationcount[lonloop, latloop, difflon, difflat] =
	                			variationcount[lonloop, latloop, difflon, difflat] + sum(@. count(!ismissing(celldiffs)))
	                	end
	                end
				end
			end
		end
	end

	print("\033[1K\rCalculating means...")
	meanvar::Array{Float64, 4} = similar(variationtotal)
	meanvar = variationtotal ./ variationcount
	meanvar = replace(meanvar, NaN=>-1e35)

	out::IOStream = open(OUTFILE, "w")
	serialize(out, meanvar)
	close(out)

	print("\n")
end

@time run()
