using NCDatasets
using Serialization

const INFILE = "daily.nc"
const OUTFILE = "pco2_spatial_variation.jldata"

function run()
	print("Loading input data...")
	local lons::Array{Union{Missing, Float32},1} = Dataset(INFILE)["longitude"][:]
	local lats::Array{Union{Missing, Float32},1} = Dataset(INFILE)["latitude"][:]
	local pco2::Array{Union{Missing, Float64}, 3} = Dataset(INFILE)["pCO2"][:,:,:]

	local lonsize::Int64 = length(lons)
	local latsize::Int64 = length(lats)
	local timesize::Int64 = size(pco2)[3]

	print("\033[1K\rInitialising data structures...")
	local variationtotal::Array{Float64, 4} = zeros(lonsize, latsize, lonsize, latsize)
	local variationcount::Array{Int64, 4} = zeros(lonsize, latsize, lonsize, latsize)

	local processedcount::Int64 = 0

	Threads.@threads for timeloop::Int64 in 1:timesize
		local timemin::Int64 = timeloop - 7
		if timemin < 1
			timemin = 1
		end

		local timemax::Int64 = timeloop + 7
		if timemax > timesize
			timemax = timesize
		end

		for lonloop::Int64 in 1:lonsize
			for latloop::Int64 in 1:latsize

				@inbounds local currentvalue::Union{Missing, Float64} = pco2[lonloop, latloop, timeloop]
				if !ismissing(currentvalue)

	                # Search for other cells with values at the same time.
	                # Note that we include values from a 7 timesteps either side too,
	                # just to increase the amount of data we can use. It should also
	                # temper the most optimistic change values we see.
	                local compare::Array{Union{Missing, Float64}, 3} = pco2[:, :, timemin:timemax]
	                local diffs::Array{Union{Missing, Float64}, 3} = similar(compare)
	                @inbounds @. diffs = abs(compare - currentvalue)
	                for difflon in 1:lonsize
	                	for difflat in 1:latsize
	                		local celldiffs = diffs[difflon, difflat, :]
                			@sync begin
	                			@inbounds variationtotal[lonloop, latloop, difflon, difflat] =
	                				variationtotal[lonloop, latloop, difflon, difflat] + sum(skipmissing(celldiffs))
	                			@inbounds variationcount[lonloop, latloop, difflon, difflat] =
	                				variationcount[lonloop, latloop, difflon, difflat] + sum(@. count(!ismissing(celldiffs)))
	                		end
	                	end
	                end
				end
			end
		end

		@sync begin
		    processedcount = processedcount + 1
			Core.print("\033[1K\rCalculating spatial variation $processedcount / $timesize")
		end
	end

	print("\033[1K\rCalculating means...")
	local meanvar::Array{Float64, 4} = similar(variationtotal)
	@inbounds meanvar = variationtotal ./ variationcount
	@inbounds meanvar = replace(meanvar, NaN=>-1e35)

	local out::IOStream = open(OUTFILE, "w")
	serialize(out, meanvar)
	close(out)

	print("\n$(length(meanvar)) $(count(meanvar .< 0))")
	print("\n")
end

@time run()
