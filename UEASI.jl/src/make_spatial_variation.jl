using NCDatasets
using Serialization
using SharedArrays
using ProgressMeter

# I had a go at parallelising this, but @threads doesn't work
# because array updates aren't thread safe.
# We'd have to make a proper @distributed/pmap
# to make separate small arrays and combine them.
# Since it's debatable how much time it saves anyway
# in this program I've not bothered.

const INFILE = "daily.nc"
const OUTFILE = "fco2_spatial_variation.jldata"

function run()
	print("Loading input data...")

	local lons::Array{Union{Missing, Float32},1} = Dataset(INFILE)["longitude"][:]
	local lats::Array{Union{Missing, Float32},1} = Dataset(INFILE)["latitude"][:]
	local fco2::Array{Union{Missing, Float64}, 3} = Dataset(INFILE)["fCO2"][:,:,:]

	local lonsize::Int64 = length(lons)
	local latsize::Int64 = length(lats)
	local timesize::Int64 = size(fco2)[3]

	print("\033[1K\rInitialising data structures...")
	local variationtotal::Array{Float64, 4} = zeros(lonsize, latsize, lonsize, latsize)
	local variationcount::Array{Int64, 4} = zeros(lonsize, latsize, lonsize, latsize)

	local p::Progress = Progress(timesize, 1, "Calculating spatial variation")

	for timeloop::Int64 in 1:timesize
		local timemin::Int64 = timeloop - 7
		if timemin < 1
			timemin = 1
		end

		local timemax::Int64 = timeloop + 7
		if timemax > timesize
			timemax = timesize
		end

		for lonloop::Int64 in 1:lonsize
			@inbounds for latloop::Int64 in 1:latsize

				local currentvalue::Union{Missing, Float64} = fco2[lonloop, latloop, timeloop]
				if !ismissing(currentvalue)

					for timeoffset in -7:7
						innertimeindex = timeloop + timeoffset
						if innertimeindex > 0 && innertimeindex ≤ timesize
							for innerlonloop::Int64 in 1:lonsize
								for innerlatloop::Int64 in 1:latsize

									if innerlonloop != lonloop || innerlatloop != latloop
										if !ismissing(fco2[innerlonloop, innerlatloop, innertimeindex])
											diff::Float64 = abs(currentvalue - fco2[innerlonloop, innerlatloop, innertimeindex])
											variationtotal[lonloop, latloop, innerlonloop, innerlatloop] =
											  variationtotal[lonloop, latloop, innerlonloop, innerlatloop] + diff

											variationcount[lonloop, latloop, innerlonloop, innerlatloop] =
											  variationcount[lonloop, latloop, innerlonloop, innerlatloop] + 1
										end
									end
								end
							end
						end
				  end
				end
			end
		end
		next!(p)
	end

	finish!(p)

	print("\033[1K\rCalculating means...")
	local meanvar::Array{Float64, 4} = similar(variationtotal)
	@inbounds meanvar = variationtotal ./ variationcount

	local out::IOStream = open(OUTFILE, "w")
	serialize(out, meanvar)
	close(out)

	print("\n")
end

@time run()