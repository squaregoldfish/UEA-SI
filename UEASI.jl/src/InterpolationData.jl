module InterpolationData
using Distributed
using ProgressMeter
using Serialization
using Logging
using Statistics
using LsqFit
using NCDatasets

export Cell
export makecells
export interpolatecell
export InterpolationCellData

const NO_CELL = 65535
const ETOPO_FILE = "etopo60.cdf"
const ETOPO_URL = "https://github.com/NOAA-PMEL/FerretDatasets/blob/master/data/etopo60.cdf"
const EARTH_RADIUS = 6367.5
const SPATIAL_ACF_LAG_STEP = 25

# Interpolation Limits
const MAX_STDEV = 75.0
const MIN_TIME_SPAN = 1825 # 5 years
const MONTH_END_DAYS = [31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334, 365]
const MIN_POPULATED_MONTHS = 8
const MAX_LINEAR_TREND = 4.5
const MIN_LINEAR_TREND = -2.5
const MAX_PEAKS_PER_YEAR = 2
const MAX_PEAK_RATIO = 0.33
const MIN_CURVE_RATIO = 0.5
const MAX_CURVE_RATIO = 1.5
const MAX_LIMIT_DIFFERENCE = 75
const TEMPORAL_INTERPOLATION_LIMIT = 7
const MAX_SPATIAL_INTERPOLATION_STEPS = 10
const SPATIAL_INTERPOLATION_WEIGHT_LIMIT = 0.367879


######################################################
#
# Global stuff
#
const Cell = NamedTuple{(:lon, :lat), Tuple{UInt16, UInt16}}
const InterpolationCell = NamedTuple{(:lon, :lat, :weight, :uncertainty),
    Tuple{UInt16, UInt16, Union{Missing, Float64}, Union{Missing, Float64}}}
const INTERPOLATION_DATA_DIR = "interpolation_data"

# These aren't fixed, but they must be set based on the data
# grid being processed
LON_SIZE = 0
LAT_SIZE = 0
GRID_SIZE = 0

ETOPO = nothing
ETOPO_LONS = nothing
ETOPO_LATS = nothing

# Harmonic fit functions
@. harmonicfit1(x, p) = p[1] + p[2]x +
  p[3]*sin(2*π*(x/365)) + p[4]*cos(2*π*(x/365))

@. harmonicfit2(x, p) = p[1] + p[2]x +
  p[3]*sin(2*π*(x/365)) + p[4]*cos(2*π*(x/365)) +
  p[5]*sin(2*π*2*(x/365)) + p[6]*cos(2*π*2*(x/365))

@. harmonicfit3(x, p) = p[1] + p[2]x +
  p[3]*sin(2*π*(x/365)) + p[4]*cos(2*π*(x/365)) +
  p[5]*sin(2*π*2*(x/365)) + p[6]*cos(2*π*2*(x/365)) +
  p[7]*sin(2*π*3*(x/365)) + p[8]*cos(2*π*3*(x/365))

@. harmonicfit4(x, p) = p[1] + p[2]x +
  p[3]*sin(2*π*(x/365)) + p[4]*cos(2*π*(x/365)) +
  p[5]*sin(2*π*2*(x/365)) + p[6]*cos(2*π*2*(x/365)) +
  p[7]*sin(2*π*3*(x/365)) + p[8]*cos(2*π*3*(x/365)) +
  p[9]*sin(2*π*4*(x/365)) + p[10]*cos(2*π*4*(x/365))


function init(lonsize::Int64, latsize::Int64)
    global LON_SIZE, LAT_SIZE, GRID_SIZE
    LON_SIZE = lonsize
    LAT_SIZE = latsize
    GRID_SIZE = 360.0 / LON_SIZE

    if !isfile(ETOPO_FILE)
        println("\nMissing bathymetry file $ETOPO_FILE")
        println("Available from $ETOPO_URL")
        exit()
    else
        global ETOPO, ETOPO_LONS, ETOPO_LATS
        ETOPO = convert.(Float64, Dataset(ETOPO_FILE)["ROSE"][:,:])
        ETOPO_LONS = convert.(Float64, Dataset(ETOPO_FILE)["ETOPO60X"][:,:])
        ETOPO_LATS = convert.(Float64, Dataset(ETOPO_FILE)["ETOPO60Y"][:,:])
    end
end

######################################################
#
# Data

mutable struct InterpolationCellData
    cell::Cell
    land::Bool
    finished::Bool
    fitparams::Vector{Float64}
    originalinputseries::Vector{Union{Missing, Float64}}
    originalinputuncertainties::Vector{Union{Missing, Float64}}
    originalinputweights::Vector{Union{Missing, Float64}}
    paraminputseries::Vector{Union{Missing, Float64}}
    paraminputuncertainties::Vector{Union{Missing, Float64}}
    paraminputweights::Vector{Union{Missing, Float64}}

    function InterpolationCellData(cell::Cell, timesize::Int64, island::Bool, fco2::Vector{Union{Missing, Float64}},
        uncertainty::Vector{Union{Missing, Float64}})

        newobj::InterpolationCellData = new(cell)

        newobj.land = island
        newobj.finished = island # If the cell is on land, we've already finished
        newobj.fitparams = Vector{Float64}(undef, 5)

        newobj.originalinputseries = fco2
        newobj.originalinputuncertainties = uncertainty
        newobj.originalinputweights = newobj.originalinputuncertainties
        replace!(newobj.originalinputweights, !ismissing=>1.0)

        newobj.paraminputseries = copy(newobj.originalinputseries)
        newobj.paraminputuncertainties = copy(newobj.originalinputuncertainties)
        newobj.paraminputweights = copy(newobj.originalinputweights)

        newobj
    end
end #InterpolationCellData

# Create the cell objects without initialising the corresponding data structures
# Useful during development
function makecells(lonsize::Int64, latsize::Int64)::Vector{Cell}
    local cells::Vector{Cell} = Vector{Cell}(undef, lonsize * latsize)
    for i in 1:(lonsize * latsize)
        cells[i] = _makecell(i, lonsize, latsize)
    end

    return cells
end

# Create the array of base data structures for the interpolations
function makecells(lonsize::Int64, latsize::Int64, timesize::Int64, seamask::Array{UInt8, 2},
    fco2::Array{Union{Missing, Float64}, 3}, uncertainty::Array{Union{Missing, Float64}, 3})::Vector{Cell}

    if isdir(INTERPOLATION_DATA_DIR)
        rm(INTERPOLATION_DATA_DIR, recursive=true)
    end

    mkdir(INTERPOLATION_DATA_DIR)

    local cells::Vector{Cell} = Vector{Cell}(undef, lonsize * latsize)
    @showprogress 1 "Initialising working data..." for i in 1:(lonsize * latsize)
        cells[i] = makeinterpolationdata(i, lonsize, latsize, seamask, timesize, fco2, uncertainty)
    end

    return cells
end #_makecells

# Generate the basic interpolation data
function makeinterpolationdata(cellindex::Int64, lonsize::Int64, latsize::Int64, seamask::Array{UInt8, 2}, timesize::Int64,
    fco2::Array{Union{Missing, Float64}, 3}, uncertainty::Array{Union{Missing, Float64}, 3})::Cell

    local cell::Cell = _makecell(cellindex, lonsize, latsize)
    local land::Bool = seamask[cell.lon, cell.lat] == 0
    local interpolationdata::InterpolationCellData =
        InterpolationCellData(cell, timesize, land, fco2[cell.lon, cell.lat, :], uncertainty[cell.lon, cell.lat, :])
    _saveinterpolationdata(interpolationdata)

    return cell
end

# Perform the main interpolation for a cell
function interpolatecell(cell::Cell, interpolationstep::UInt8, temporalacf::Vector{Float64},
    spatialacfs::Array{Union{Missing, Float64}, 4}, spatialvariation::Array{Union{Missing, Float64}, 4}, seamask::Array{UInt8, 2})

    data::InterpolationCellData = _loadinterpolationdata(cell)

    if !data.finished
        local logger::Tuple{SimpleLogger, IOStream} = _makelogger(cell, interpolationstep)
        interpolate!(data, interpolationstep, temporalacf, spatialacfs, spatialvariation, seamask, logger[1])
        _closelogger(logger[2])
    end

    return data.finished
end


###############################################################
#
# Interpolation stuff

# Structure to hold a series and its related data
mutable struct SeriesData
    measurements::Vector{Union{Missing, Float64}}
    uncertainties::Vector{Union{Missing, Float64}}
    weights::Vector{Union{Missing, Float64}}
    curve::Vector{Float64}
    curveparams::Vector{Float64}

    function SeriesData()
        newobj::SeriesData = new()
        newobj.measurements = Vector{Union{Missing, Float64}}()
        newobj.uncertainties = Vector{Union{Missing, Float64}}()
        newobj.weights = Vector{Union{Missing, Float64}}()
        newobj.curve = Vector{Float64}()
        newobj.curveparams = Vector{Float64}()
        newobj
    end

    function SeriesData(measurements::Vector{Union{Missing, Float64}}, uncertainties::Vector{Union{Missing, Float64}}, weights::Vector{Union{Missing, Float64}})
        newobj::SeriesData = new(deepcopy(measurements), deepcopy(uncertainties), deepcopy(weights))
        newobj.curve = Vector{Float64}()
        newobj.curveparams = Vector{Float64}()
        newobj
    end

    function SeriesData(series::SeriesData)
        newobj::SeriesData = new()
        newobj.measurements = deepcopy(series.measurements)
        newobj.uncertainties = deepcopy(series.uncertainties)
        newobj.weights = deepcopy(series.weights)
        newobj.curve = deepcopy(series.curve)
        newobj.curveparams = deepcopy(series.curveparams)
        newobj
    end
end

# Perform interpolation
function interpolate!(data::InterpolationCellData, step::UInt8, temporalacf::Vector{Float64},
    spatialacfs::Array{Union{Missing, Float64}, 4}, spatialvariation::Array{Union{Missing, Float64}, 4},
    seamask::Array{UInt8, 2}, logger::SimpleLogger)

    with_logger(logger) do

        local fitfound::Bool = false
        local continuefit::Bool = true

        local spatialinterpolationcount::UInt8 = 0
        local spatialinterpolationcells::Vector{InterpolationCell} = Vector{InterpolationCell}()

        # Variables for processing (empty placeholders)
        local currentseries::SeriesData = SeriesData()

        # Temporary storage used while trying extra interpolations
        local storedseries::SeriesData = SeriesData()

        while continuefit
            local attemptcurvefit::Bool = true

            # Put together the series for the curve fit
            if spatialinterpolationcount == 0
                # Just use the original series as is
                currentseries.measurements = deepcopy(data.paraminputseries)
                currentseries.uncertainties = deepcopy(data.paraminputuncertainties)
                currentseries.weights = deepcopy(data.paraminputweights)
            else
                println("Interpolation $spatialinterpolationcount of $(length(spatialinterpolationcells))")
                # Need to do spatial interpolation
                dospatialinterpolation!(data, spatialinterpolationcells, spatialinterpolationcount,
                    spatialacfs, spatialvariation, seamask)
            end

            if continuefit
                if attemptcurvefit
                    fitfound = attemptfit!(currentseries, temporalacf)

                    if fitfound
                        if length(storedseries.curve) == 0
                            # No previous successful fit - store this fit
                            # Store the current fit for later comparison
                            storedseries = SeriesData(currentseries)

                            @debug "Trying extra interpolation"

                            # Now we attempt a spatial interpolation to see if this makes a difference.
                            # I think this won't happen if we get a valid fit at the first attempt, and it should (Issue #5)
                            if length(spatialinterpolationcells) == 0
                                continuefit = false
                            elseif spatialinterpolationcount > 0 && spatialinterpolationcount == length(spatialinterpolationcells)
                                continuefit = false
                            end
                        else
                            # There is a previous fit
                            if significantdifference(storedseries, currentseries)
                                # The interpolated fit is different to the previous fit.
                                # Continue the interpolation to refine the fit further
                                storedseries = SeriesData(currentseries)

                                if spatialinterpolationcount > 0
                                    if length(spatialinterpolationcells) == 0
                                        continuefit = false
                                    elseif spatialinterpolationcount > 0 && spatialinterpolationcount == length(spatialinterpolationcells)
                                        continuefit = false
                                    end
                                end
                            else
                                # The spatial interpolation made no difference. We can stop
                                # with the previous fit
                                @debug "Curve difference not significant - reverting"
                                continuefit = false
                            end
                        end
                    end
                end

                if !attemptcurvefit || length(currentseries.curve) == 0
                    # We didn't get a successful fit. Do a spatial interpolation to get more
                    # measurements in play. Unless, of course, there are no more candidate cells
                    # for interpolation
                    if spatialinterpolationcount > 0
                        if length(spatialinterpolationcells) == 0
                            continuefit = false
                        elseif spatialinterpolationcount > 0 && spatialinterpolationcount == length(spatialinterpolationcells)
                            continuefit = false
                        end
                    end
                end

                if continuefit
                    spatialinterpolationcount += 1
                    if spatialinterpolationcount > length(spatialinterpolationcells)
                        continuefit = false
                    end
                end
            end

            # Store data from this attempt, successful or otherwise
            data.paraminputseries = currentseries.measurements
            data.paraminputweights = currentseries.weights
            data.paraminputuncertainties = currentseries.uncertainties
            data.fitparams = currentseries.curveparams

            if fitfound
                @info "FIT SUCCESS"
                data.finished = true
            else
                @info "FIT FAILED"
            end

            _saveinterpolationdata(data)
        end
    end
end

# Try to fit a curve to the supplied time series
function attemptfit!(series::SeriesData, temporalacf::Vector{Float64})::Bool

    local fitsuccess::Bool = false

    # Attempt a curve fit
    fitsuccess = fitcurve!(series)

    if length(series.curveparams) == 0
        dotemporalinterpolation!(series, temporalacf)
        fitsuccess = fitcurve!(series)
    end

    fitsuccess
end

# Fit a harmonic curve to a time series
function fitcurve!(series::SeriesData)::Bool

    local fitsuccess::Bool = false

    local fitseries::Vector{Union{Missing, Float64}} = removeoutliers(series.measurements)

    if !doprefitcheck(fitseries)
        @info "PREFIT CHECKS FAILED"
    else
        local harmoniccount::UInt8 = 4

        while !fitsuccess && harmoniccount > 0
            @debug "Fitting $harmoniccount harmonic(s)"

            # TODO Handle fit failure
            local days::Vector{UInt16} = findall((!ismissing).(fitseries))
            local p0::Vector{Float64} = zeros(harmoniccount * 2 + 2)
            local functionname::String = "harmonicfit$harmoniccount"

            local fit::LsqFit.LsqFitResult = curve_fit(getfield(InterpolationData, Symbol(functionname)), days, collect(skipmissing(fitseries)), p0)
            local fitparams::Vector{Float64} = fit.param
            local fittedcurve::Vector{Float64} = makecurve(fitparams, length(fitseries))
            @debug "Fitted params: $fitparams"

            # TODO Handle fit failure
            if true
                # Check that the fit is OK

                # Slope
                local slope::Float64 = fitparams[2]
                @debug "Fitted slope is $slope"
                if slope > MAX_LINEAR_TREND || slope < MIN_LINEAR_TREND
                    @info "Curve slope is outside $MIN_LINEAR_TREND - $MAX_LINEAR_TREND range"
                    fitsuccess = false
                else
                    fitsuccess = true
                end

                if fitsuccess && harmoniccount > 1
                    fitsuccess = checkcurvepeaks(fitparams)
                end

                if fitsuccess
                    fitsuccess = checkcurvefit(fitseries, fittedcurve)
                end
            end

            if fitsuccess
                series.curve = fittedcurve
                series.curveparams = fitparams
            else
                harmoniccount -= 1
            end
        end
    end

    fitsuccess
end

# Remove outliers from a time series. The outlier limit
# is defined in MAX_STDEV
function removeoutliers(series::Vector{Union{Missing, Float64}})::Vector{Union{Missing, Float64}}
    local filteredseries::Vector{Union{Missing, Float64}} = Vector{Union{Missing, Float64}}(missing, length(series))

    local meanvalue::Float64 = mean(skipmissing(series))
    local stdev::Float64 = std(skipmissing(series), mean=meanvalue)

    for i in 1:length(series)
        if !ismissing(series[i])
            local stdevsfrommean::Float64 = abs(series[i] - meanvalue) / stdev
            if stdevsfrommean < 3.0
                filteredseries[i] = series[i]
            end
        end
    end

    @debug "Removed $(_countvalues(series) - _countvalues(filteredseries)) outliers"
    filteredseries
end

# Check that a time series is suitable to have a curve fit attempted on it
function doprefitcheck(series::Vector{Union{Missing, Float64}})::Bool
    local ok::Bool = true

    # Standard deviation
    local stdev = std(skipmissing(series))
    @debug "Series standard deviation = $stdev"
    if stdev > MAX_STDEV
        @info "Standard too large, should be ≤ $MAX_STDEV"
        ok = false
    end

    # Minimum time coverage
    local measurementdays::Vector{Int64} = findall(!ismissing, series)
    if length(measurementdays) == 0
        @info "No measurements"
        ok = false
    end

    if ok
        local dayspan::Int64 = measurementdays[end] - measurementdays[1]
        @debug "Measurements span $dayspan days"
        if dayspan < MIN_TIME_SPAN
            @info "Measurements must span at least $MIN_TIME_SPAN days"
            ok = false
        end
    end

    # Populated months
    if ok
        local populatedmonths::Vector{Bool} = falses(12)
        for i in 1:length(series)
            if !ismissing(series[i])
                local day::Int16 = mod(i, 365)
                if day == 0
                    day = 365
                end

                populatedmonths[findall(day .≤ MONTH_END_DAYS)[1]] = true
            end
        end

        local populatedmonthcount::UInt8 = count(populatedmonths)
        @debug "$populatedmonthcount populated months"
        if populatedmonthcount < MIN_POPULATED_MONTHS
            @info "Should have at least $MIN_POPULATED_MONTHS populated months"
            ok = false
        end
    end

    ok
end

# Perform temporal interpolation
function dotemporalinterpolation!(series::SeriesData, temporalacf::Vector{Float64})

    local interpolatedmeasurements::Vector{Union{Missing, Float64}} = Array{Union{Missing, Float64}}(missing, length(series.measurements))
    local interpolatedweights::Vector{Union{Missing, Float64}} = Array{Union{Missing, Float64}}(missing, length(series.measurements))
    local interpolateduncertainties::Vector{Union{Missing, Float64}} = Array{Union{Missing, Float64}}(missing, length(series.measurements))

    for i in 1:length(series.measurements)

        if ismissing(series.measurements[i])
            local interpstart::Int64 = i - TEMPORAL_INTERPOLATION_LIMIT
            if interpstart < 1
                interpstart = 1
            end

            local interpend::Int64 = i + TEMPORAL_INTERPOLATION_LIMIT
            if interpend > length(series.measurements)
                interpend = length(series.measurements)
            end

            local interptotal::Float64 = 0.0
            local interpweightsum::Float64 = 0.0
            local interpcount::Int16 = 0

            for j in interpstart:interpend
                if !ismissing(series.measurements[j])
                    local lag::Int16 = abs(i - j)
                    interptotal += series.measurements[j] * temporalacf[lag]
                    interpweightsum += temporalacf[lag]
                    interpcount += 1
                end
            end

            if interpcount > 0
                interpolatedmeasurements[i] = interptotal / interpweightsum
                interpolatedweights[i] = interpweightsum / interpcount

                # The uncertainty doesn't matter, since the temporally interpolated
                # values are only temporary anyway - they don't make it to the final
                # interpolated product
                interpolateduncertainties[i] = 0
            end
        end
    end

    mergevalues!(series.measurements, interpolatedmeasurements)
    mergevalues!(series.weights, interpolatedweights)
    mergevalues!(series.uncertainties, interpolateduncertainties)
end

# Merge two sets of values. Where two values overlap,
# use the original
function mergevalues!(original::Vector{Union{Missing, Float64}}, incoming::Vector{Union{Missing, Float64}})
    for i in 1:length(original)
        if ismissing(original[i]) && !ismissing(incoming[i])
            original[i] = incoming[i]
        end
    end
end

# Perform one step of spatial interpolation
function dospatialinterpolation!(data::InterpolationCellData, interpolationcells::Vector{InterpolationCell},
    interpolationcount::UInt8, spatialacfs::Array{Union{Missing, Float64}, 4},
    spatialvariation::Array{Union{Missing, Float64}, 4}, seamask::Array{UInt8, 2})

    # Build the list of interpolation cells if it hasn't been done yet
    if length(interpolationcells) == 0
        interpolationcandidates::Vector{Cell} = calculateinterpolationcells(data.cell, MAX_SPATIAL_INTERPOLATION_STEPS)
        append!(interpolationcells, filterandsortinterpolationcells(data.cell,
            interpolationcandidates, data.paraminputweights, data.paraminputuncertainties,
            spatialacfs, spatialvariation, seamask))
    end

    # Only interpolate if there are interpolation cells
    if length(interpolationcells) > 0 && interpolationcount < length(interpolationcells)
        @info "Performing spatial interpolation for cell $(data.cell.lon), $(data.cell.lat), cell $interpolationcount of $(length(interpolationcells))"

        local interpolatedmeasurements::Array{Union{Missing, Float64}, 2} = Array{Union{Missing, Float64}}(missing, interpolationcount, length(data.originalinputseries))
        local interpolatedweights::Array{Union{Missing, Float64}, 2} = Array{Union{Missing, Float64}}(missing, interpolationcount, length(data.originalinputseries))
        local interpolateduncertainties::Array{Union{Missing, Float64}, 2} = Array{Union{Missing, Float64}}(missing, interpolationcount, length(data.originalinputseries))

        # Interpolate the specified number of cells
        for i in 1:interpolationcount
            local interpcell::InterpolationCell = interpolationcells[i]
            @info "Interpolating to $(data.cell.lon), $(data.cell.lat) from $(interpcell.lon), $(interpcell.lat)"
            interpdata::InterpolationCellData = _loadinterpolationdata(interpcell)
            interpolatedmeasurements[i,:] = interpdata.originalinputseries
            interpolatedweights[i,:] = interpdata.originalinputweights
            interpolateduncertainties[i,:] = interpdata.originalinputuncertainties
        end

        # Loop through the time steps combining values from the original input
        # and the interpolated cells
        for i in 1:length(data.originalinputseries)

            # Count of values for the time step
            local valuecount::Int16 = 0

            # Measurements
            local weightedmeasurementtotal::Float64 = 0
            local weighttotal::Float64 = 0
            local uncertaintysquaretotal::Float64 = 0

            if !ismissing(data.originalinputseries[i])
                valuecount += 1
                weightedmeasurementtotal += data.originalinputseries[i] * data.originalinputweights[i]
                weighttotal += data.originalinputweights[i]
                uncertaintysquaretotal += data.originalinputuncertainties[i]^2
            end

            for j in 1:interpolationcount
                if !ismissing(interpolatedmeasurements[j, i])
                    valuecount += 1
                    weightedmeasurementtotal += interpolatedmeasurements[j, i] * interpolatedweights[j, i]
                    weighttotal += interpolatedweights[j, i]
                    uncertaintysquaretotal += interpolateduncertainties[j, i]^2
                end
            end

            if valuecount == 0
                data.paraminputseries[i] = missing
                data.paraminputweights[i] = missing
                data.paraminputuncertainties[i] = missing
            else
                data.paraminputseries[i] = weightedmeasurementtotal / weighttotal
                data.paraminputweights[i] = weighttotal / valuecount
                data.paraminputuncertainties[i] = √(uncertaintysquaretotal / valuecount)
            end
        end
    end
end

# Determine candidate cells for spatial interpolation of a given cell
function calculateinterpolationcells(cell::Cell, maxstep::Int64)::Vector{Cell}

    # Work out how many interpolation cells there will be.
    # Each step means (8 * step) cells are added to the list.
    # The number of cells is therefore the triangular number of steps * 8
    # I have discovered a truly marvellous proof of why this is the case,
    # but this comment is too brief to contain it.
    local interpolationcells::Vector{Cell} = Vector{Cell}()

    for step_loop in 1:maxstep
        for y in (step_loop * -1):step_loop
            local celly::UInt16 = calcspatialinterpycell(cell.lat, y)
            if celly != NO_CELL

                # For the top and bottom rows of the step grid, add all horizontal cells
                if abs(y) == step_loop
                    for x in (step_loop * -1):step_loop
                        local cellx::UInt16 = calcspatialinterpxcell(cell.lon, x)
                        push!(interpolationcells, (lon = cellx, lat = celly))
                    end
                else
                    # For all other rows, just add the left and right edges
                    local cellx::UInt16 = calcspatialinterpxcell(cell.lon, (step_loop * -1))
                    push!(interpolationcells, (lon = cellx, lat = celly))

                    cellx = calcspatialinterpxcell(cell.lon, step_loop)
                    push!(interpolationcells, (lon = cellx, lat = celly))
                end
            end
        end
    end

    interpolationcells
end

# Filter a set of interpolation cell candidates to remove unusable cells and sort them by
# lowest uncertainty/highest weight
function filterandsortinterpolationcells(cell::Cell, interpolationcandidates::Vector{Cell},
    weights::Vector{Union{Missing, Float64}}, uncertainties::Vector{Union{Missing, Float64}},
    spatialacfs::Array{Union{Missing, Float64}, 4}, spatialvariation::Array{Union{Missing, Float64}, 4},
    seamask::Array{UInt8, 2})::Vector{InterpolationCell}

    # The output of this is a table of cell, weight, uncertainty.
    # Sorted by uncertainty, distance between cells, weight, cell.lon, cell.lat
    local interpolationcells::Vector{InterpolationCell} = Vector{InterpolationCell}()

    local cellcount::Int16 = 0

    for candidate in interpolationcandidates
        # Only process sea cells
        if seamask[candidate.lon, candidate.lat] == 1
            if !islandbetween(cell, candidate)
                local weight::Union{Float64, Missing} = getspatialinterpolationweight(cell, candidate, spatialacfs)

                if !ismissing(weight) && weight ≥ SPATIAL_INTERPOLATION_WEIGHT_LIMIT
                    cellcount += 1

                    # Now get the uncertainty. The list of candidate cells will be sorted by this
                    local uncertainty::Union{Float64, Missing} = getspatialinterpolationuncertainty(cell, candidate, spatialvariation)
                    push!(interpolationcells, _makeinterpolationcell(candidate, weight, uncertainty))
                end
            end
        end
    end

    sort!(interpolationcells, by = x -> (x.uncertainty, x.weight))
    interpolationcells
end

# Retrieve the uncertainty for the interpolation between two cells
function getspatialinterpolationuncertainty(source::Cell, dest::Cell,
    spatialvariation::Array{Union{Missing, Float64}, 4})::Union{Float64, Missing}

    local uncertainty::Union{Float64, Missing} = missing

    # Here we use the pco2_spatial_variation object loaded at the start of the program.
    # This is a 4D array of <target_lon, target_lat, interp_lon, interp_lat>
    # The values in this object are the mean difference in pCO2 between two cells
    # within a small time period.

    # We can use the uncertainty from both the interp and target cells, since they are
    # are in the same direction and therefore equivalent
    local forwarduncertainty::Union{Missing, Float64} = spatialvariation[source.lon, source.lat, dest.lon, dest.lat]
    local reverseuncertainty::Union{Missing, Float64} = spatialvariation[dest.lon, dest.lat, source.lon, source.lat]

    if length(collect(skipmissing([forwarduncertainty, reverseuncertainty]))) > 0
        uncertainty = mean(skipmissing([forwarduncertainty, reverseuncertainty]))
    else
        # If we can't find an uncertainty for the specific cell combination, use the
        # mean uncertainty from the cells surrounding the target cell
        local uncertaintycells::Vector{Cell} = calculateinterpolationcells(source, 1)
        local uncertainties::Vector{Union{Missing, Float64}} = []

        for uncertaintycell in uncertaintycells
            push!(uncertainties, spatialvariation[source.lon, source.lat, uncertaintycell.lon, uncertaintycell.lat])
        end

        if length(collect(skipmissing(uncertainties))) > 0
            uncertainty = mean(skipmissing(uncertainties))
        end
    end

    uncertainty
end

# Calculate the weighting to be used for spatial interpolation between two cells
function getspatialinterpolationweight(source::Cell, dest::Cell, spatialacfs::Array{Union{Missing, Float64}, 4})::Union{Float64, Missing}

    local weight::Union{Float64, Missing} = missing

    # Here we use the mean_directional_acfs object loaded at the start of the program.
    # This is a 4D array of lat, lon, direction, distance
    # Distances are in bins of ACF_LAG_STEP km. We calculate the weighting from the
    # centre of the two cells.
    local distance::Float64 = calccelldistance(source, dest)
    local distancebin::Int16 = floor(distance / SPATIAL_ACF_LAG_STEP) + 1
    local bearing::Float64 = calccellbearing(source, dest)
    local acfdirection::UInt8 = getacfdirection(bearing)

    local sourcecentrelon::Float64 = (source.lon * GRID_SIZE) - (GRID_SIZE / 2)
    local sourcecentrelat::Float64 = (source.lat * GRID_SIZE) - (90 + (GRID_SIZE / 2))
    local destcentrelon::Float64 = (dest.lon * GRID_SIZE) - (GRID_SIZE / 2)
    local destcentrelat::Float64 = (dest.lat * GRID_SIZE) - (90 + (GRID_SIZE / 2))


    # We can use the ACF from both the source and dest cells, since they are
    # are in the same direction and therefore equivalent
    local acfvalue::Union{Missing, Float64} = spatialacfs[source.lon, source.lat, acfdirection, distancebin]
    local reverseacfvalue::Union{Missing, Float64} = spatialacfs[dest.lon, dest.lat, acfdirection, distancebin]

    if length(collect(skipmissing([acfvalue, reverseacfvalue]))) > 0
        weight = mean(skipmissing([acfvalue, reverseacfvalue]))
    else
        # If there are no directional ACF values, use the omnidirectional value
        # for the target cell only. If this doesn't exist, we have to return a zero
        # weight because we can't guess at whether or not the two cells' values are related
        local omnidirectionalacfvalue::Union{Missing, Float64} = spatialacfs[source.lon, source.lat, 5, distancebin]
        if !ismissing(omnidirectionalacfvalue)
            weight = omnidirectionalacfvalue
        end
    end

    # Anything below the weighting limit is given a weighting of zero
    if !ismissing(weight) && weight < SPATIAL_INTERPOLATION_WEIGHT_LIMIT
        weight = missing
    end

    weight
end

# Get the ACF direction of a bearing
#
# 1 = NORTH/SOUTH
# 2 = EAST/WEST
# 3 = NE/SW
# 4 = SE/NW
function getacfdirection(bearing::Float64)::UInt8

    direction::UInt8 = 0

    if bearing ≤ 22.5 ||
       bearing ≥ 337.5 ||
       (bearing ≥ 157.5 && bearing ≤ 202.5)

       direction = 1

    elseif (bearing ≥ 67.5 && bearing ≤ 112.5) ||
           (bearing ≥ 247.5 && bearing ≤ 292.5)

        direction = 2

    elseif (bearing ≥ 22.5 && bearing ≤ 67.5) ||
           (bearing ≥ 202.5 && bearing ≤ 247.5)

        direction = 3

    elseif (bearing ≥ 112.5 && bearing ≤ 157.5) ||
            (bearing ≥ 292.5 && bearing ≤ 337.5)

        direction = 4

    end

    direction
end

# Calculate the distance between two cells in km
#
# θ = lat, λ = lon (radians)
function calccelldistance(cell1::Cell, cell2::Cell)::Float64

    # Calculate lats and lons in radians
    local λ1::Float64, θ1::Float64 = _getcellradians(cell1)
    local λ2::Float64, θ2::Float64 = _getcellradians(cell2)

    local ∆θ::Float64 = θ2 - θ1
    local ∆λ::Float64 = λ2 - λ1

    local a::Float64 = sin(∆θ/2)^2 + cos(θ1) * cos(θ2) * sin(∆λ/2)^2
    local c::Float64 = 2 * atan(√a, √(1-a))

    c * EARTH_RADIUS
end

# Calculate the bearing between two cells, on rhumb lines
# θ = lat, λ = lon (radians)
# Bearing is returned in degrees
function calccellbearing(source::Cell, dest::Cell)::Float64

    # Calculate lats and lons in radians
    local λ1::Float64, θ1::Float64 = _getcellradians(source)
    local λ2::Float64, θ2::Float64 = _getcellradians(dest)

    local ∆λ::Float64 = λ2 - λ1
    local ∆ψ::Float64 = log(tan(π/4 + θ2/2) / tan(π/4 + θ1/2))

    if abs(∆λ) > π
        ∆λ = ∆λ > 0 ? -2π - ∆λ : 2π + ∆λ
    end

    local bearing::Float64 = atan(∆λ, ∆ψ) * 180/π
    if bearing < 0
        bearing = 360.0 - abs(bearing)
    end

    bearing
end

# Determine whether there is any land between two cells
function islandbetween(source::Cell, target::Cell)::Bool

    local sourcex::Int16 = getetopolonindex(source.lon)
    local sourcey::Int16 = getetopolatindex(source.lat)
    local targetx::Int16 = getetopolonindex(target.lon)
    local targety::Int16 = getetopolatindex(target.lat)

    # Step through all cells in etopo between the two cells passed in
    # This uses an extension of the Bresenham algorithm by Eugen Dedu
    # http://lifc.univ-fcomte.fr/~dedu/projects/bresenham/index.html
    local foundland::Bool = false

    local i::Int16 = 0
    local ystep::Int16 = 0
    local xstep::Int16 = 0
    local currenterror::Int16 = 0
    local preverror::Int16 = 0
    local y::Int16 = sourcey
    local x::Int16 = sourcex
    local ∂∂y::Int16 = 0
    local ∂∂x::Int16 = 0
    local ∂y::Int16 = targety - sourcey
    local ∂x::Int16 = targetx - sourcex
    local currentx::Int16 = sourcex
    local currenty::Int16 = sourcey

    if ∂y < 0
        ystep = -1
        ∂y = -∂y
    else
        ystep = 1
    end

    if ∂x < 0
        xstep = -1
        ∂x = -∂x
    else
        xstep = 1
    end

    ∂∂y = 2∂y
    ∂∂x = 2∂x

    if ∂∂x ≥ ∂∂y
        currenterror = ∂x
        preverror = ∂x

        for i in 0:(∂x - 1)
            if foundland
                break
            end

            x += xstep
            currenterror += ∂∂y

            if currenterror > ∂∂x
                y += ystep
                currenterror -= ∂∂x

                if currenterror + preverror < ∂∂x
                    currenty = y - ystep
                    currentx = x
                    foundland = isetopoland(currentx, currenty)
                elseif currenterror + preverror > ∂∂x
                    currenty = y
                    currentx = x - xstep
                    foundland = isetopoland(currentx, currenty)
                else
                    currenty = y - ystep
                    currentx = x
                    foundland = isetopoland(currentx, currenty)

                    currenty = y
                    currentx = x - xstep
                    foundland = isetopoland(currentx, currenty)
                end
            end

            currenty = y
            currentx = x
            foundland = isetopoland(currentx, currenty)
            preverror = currenterror
        end
    else

        currenterror = ∂y
        preverror = ∂y

        for i in 0:(∂y - 1)
            if foundland
                break
            end

            y += ystep
            currenterror += ∂∂x

            if currenterror > ∂∂y
                x += xstep
                currenterror -= ∂∂y

                if currenterror + preverror < ∂∂y
                    currenty = y
                    currentx -= xstep
                    foundland = isetopoland(currentx, currenty)
                elseif currenterror + preverror > ∂∂y
                    currenty -= ystep
                    currentx = x
                    foundland = isetopoland(currentx, currenty)
                else
                    currenty = y
                    currentx -= xstep
                    foundland = isetopoland(currentx, currenty)

                    currenty -= ystep
                    currentx = x
                    foundland = isetopoland(currentx, currenty)
                end
            end

            currenty = y
            currentx = x
            foundland = isetopoland(currentx, currenty)

            preverror = currenterror
        end

    end

    foundland
end

function isetopoland(x::Int16, y::Int16)::Bool
    ETOPO[x, y] ≥ 0
end

# Convert an interpolation grin lon index to an ETOPO lon index
function getetopolonindex(lonindex::UInt16)::Int16
    local lon_degrees::Float64 = (lonindex * GRID_SIZE) - (GRID_SIZE / 2)

    # ETOPO longitudes are 20 - 380 for some reason
    if lon_degrees ≤ 20
        lon_degrees += 360
    end

    convert(Int16, findfirst(ETOPO_LONS .≥ floor(lon_degrees)))
end

# Convert an interpolation grin lat index to an ETOPO lat index
function getetopolatindex(latindex::UInt16)::Int16
    local lat_degrees::Float64 = (latindex * GRID_SIZE) - (90 + GRID_SIZE / 2)
    convert(Int16, findfirst(ETOPO_LATS .≥ floor(lat_degrees)))
end

# Get the longitude cell index that is x cells from the specified lon cell index
function calcspatialinterpxcell(lon::UInt16, x::Int64)::UInt16
    local newx::Int64 = lon + x
    if newx > LON_SIZE
        newx = newx - LON_SIZE
    elseif newx < 1
        newx = LON_SIZE - abs(newx)
    end

    convert(UInt16, newx)
end

# Get the latitude cell index that is y cells from the specified lat cell index
function calcspatialinterpycell(lat::UInt16, y::Int64)::UInt16
    local newy::UInt16 = convert(UInt16, lat + y)
    if newy > LAT_SIZE || newy < 1
        newy = NO_CELL
    end

    newy
end

# Make a curve using the supplied parameters
function makecurve(params::Vector{Float64}, curvelength::Int64)::Vector{Float64}

    local curve::Vector{Float64} = zeros(curvelength)

    for i in 1:length(curve)

        local value::Float64 = params[1] + params[2] * (i - 1)

        local term::UInt8 = 2
        for trigloop in 1:((length(params) - 2) / 2)
            term += 1
            value = value + params[term] * sin(2 * π * trigloop * ((i - 1) / 365))
            term += 1
            value = value + params[term] * cos(2 * π * trigloop * ((i - 1) / 365))
        end

        curve[i] = value
    end

    curve
end

# Make a curve using the supplied parameters
function makeseasonalcycle(curveparams::Vector{Float64})::Vector{Float64}

    local seasonalcycle::Vector{Float64} = zeros(365)

    for i in 1:length(seasonalcycle)

        local value::Float64 = 0.0

        local term::UInt8 = 2
        for trigloop in 1:((length(curveparams) - 2) / 2)
            term += 1
            value = value + curveparams[term] * sin(2 * π * trigloop * ((i - 1) / 365))
            term += 1
            value = value + curveparams[term] * cos(2 * π * trigloop * ((i - 1) / 365))
        end

        seasonalcycle[i] = value
    end

    seasonalcycle
end

# Check multiple harmonics in a seasonal cycle
function checkcurvepeaks(curveparams::Vector{Float64})::Bool

    local peaksok::Bool = true

    local seasonalcycle::Vector{Float64} = makeseasonalcycle(curveparams)

    local maxima::Vector{Float64} = Vector{Float64}()
    local maximapos::Vector{Int16} = Vector{Int16}()
    local minima::Vector{Float64} = Vector{Float64}()
    local minimapos::Vector{Int16} = Vector{Int16}()

    local lastvalue::Float64 = seasonalcycle[1]
    local slopedirection::Int8 = 0 # Starting value; 1 = up, -1 = down
    local startdirection::Int8 = 0
    local enddirection::Int8 = 0
    local curveamplitude::Float64 = maximum(seasonalcycle) - minimum(seasonalcycle)

    for i in 2:365
        if seasonalcycle[i] ≥ lastvalue
            if slopedirection == -1
                # Found a trough!
                push!(minima, lastvalue)
                push!(minimapos, i - 1)

            end

            slopedirection = 1

            # Record the start direction
            if startdirection == 0
                startdirection = slopedirection
            end
        else
            if slopedirection > 0
                # We found a peak!
                push!(maxima, lastvalue)
                push!(maximapos, i - 1)
            end

            slopedirection = -1

            # Record the start direction
            if startdirection == 0
                startdirection = slopedirection
            end
        end

        lastvalue = seasonalcycle[i]
    end

    # Record the end direction
    enddirection = slopedirection

    # If the start and end directions are different, there's a peak/trough at zero.
    if startdirection == 1 && enddirection == -1
        insert!(minima, 1, seasonalcycle[1])
        insert!(minimapos, 1, 1)
    elseif startdirection == -1 && enddirection == 1
        push!(maxima, seasonalcycle[365])
        push!(maximapos, 365)
    end

    # Check the number of peaks
    @debug "$(length(maxima)) peaks in the seasonal cycle"
    if length(maxima) > MAX_PEAKS_PER_YEAR
        @info "Peak check FAILED - too many peaks in seasonal cycle"
        peaksok = false
    else
        # Calculate the peak sizes
        local peaksizes::Vector{Float64} = zeros(length(maxima))

        for i in 1:length(maxima)
            local maxvalue::Float64 = maxima[i]
            local maxpos::UInt16 = maximapos[i]

            # For each peak, we find the position of the preceding minimum
            # We loop round to the end of the year if necessary
            local minentry::UInt16 = 0
            local minposcandidates::Vector{UInt16} = findall(x -> (x < maxpos), minimapos)
            if length(minposcandidates) == 0
                minentry = length(minima)
            else
                minentry = minposcandidates[end]
            end

            local minvalue::Float64 = minima[minentry]
            peaksizes[i] = abs(maxvalue - minvalue)
        end

        local peaksizelimit::Float64 = curveamplitude * MAX_PEAK_RATIO
        local largepeakcount::UInt8 = length(findall(x -> (x > peaksizelimit), peaksizes))
        @debug "Peak sizes = $(peaksizes ./ curveamplitude)"
        if largepeakcount > 1 # The 'main' cycle is detected as a peak, so ignore it
            @info "Peak check FAILED - secondary peak(s) too large"
            peaksok = false
        end
    end

    peaksok
end

function checkcurvefit(series::Vector{Union{Missing, Float64}}, fittedcurve::Vector{Float64})::Bool

    local curveok::Bool = true

    # Compare the curve range with the series range
    local seriesmax::Float64 = maximum(skipmissing(series))
    local seriesmin::Float64 = minimum(skipmissing(series))
    local curvemax::Float64 = maximum(fittedcurve)
    local curvemin::Float64 = minimum(fittedcurve)
    local seriesrange::Float64 = seriesmax - seriesmin
    local curverange::Float64 = curvemax - curvemin
    local rangeratio::Float64 = curverange / seriesrange

    @debug "Range ratio = $rangeratio"
    if rangeratio < MIN_CURVE_RATIO || rangeratio > MAX_CURVE_RATIO
        curveok = false
        @info "Ratio FAILED should be $MIN_CURVE_RATIO  to $MAX_CURVE_RATIO"
    end

    local maxdiff::Float64 = abs(seriesmax - curvemax)
    local mindiff::Float64 = abs(seriesmin - curvemin)

    @debug "Limit differences = $maxdiff, $mindiff"
    if maxdiff > MAX_LIMIT_DIFFERENCE || mindiff > MAX_LIMIT_DIFFERENCE
        curveok = false
        @info "Limit difference FAILED should be $MAX_LIMIT_DIFFERENCE or less"
    end

    curveok
end

function significantdifference(series1::SeriesData, series2::SeriesData)::Bool
    local corr::Float64 = Statistics.cor(series1.curve, series2.curve)^2
    @debug "Curve difference level = $corr"
    corr < 0.99
end

###############################################################
#
# Simple utility methods

# Initialise a logger for a cell and interpolation run
function _makelogger(cell::Cell, interpolationstep::UInt8)::Tuple{SimpleLogger, IOStream}
    logfile::String = "$(INTERPOLATION_DATA_DIR)/$(cell.lon)_$(cell.lat)_$(interpolationstep).log"
    io::IOStream = open(logfile, "w")
    logger::SimpleLogger = SimpleLogger(io, Logging.Debug, Dict{Any,Int64}())
    return (logger, io)
end

# Close a logger's IO stream
function _closelogger(io::IOStream)
    close(io)
end

# Generate the filename for an InterpolationCellData object
function _getdatafilename(data::InterpolationCellData)::String
    return _getdatafilename(data.cell)
end

# Generate the filename for an InterpolationCellData object using its cell
function _getdatafilename(cell::Cell)::String
    return "$(INTERPOLATION_DATA_DIR)/$(cell.lon)_$(cell.lat).jldata"
end

# Save an InterpolationCellData object to disk
function _saveinterpolationdata(data::InterpolationCellData)
    local out::IOStream = open(_getdatafilename(data), "w")
    serialize(out, data)
    close(out)
end

# Load an InterpolationCellData object
function _loadinterpolationdata(cell::InterpolationCell)::InterpolationCellData
    _loadinterpolationdata((lon=cell.lon, lat=cell.lat))
end

# Load an InterpolationCellData object
function _loadinterpolationdata(cell::Cell)::InterpolationCellData
    local inchan::IOStream = open(_getdatafilename(cell), "r")
    local data::InterpolationCellData = deserialize(inchan)
    close(inchan)
    return data
end

# Count the number of non-missing values in a series
function _countvalues(series::Vector{Union{Missing, Float64}})::Int64
    length(collect(skipmissing(series)))
end

# Get a cell's lon/lat in radians
function _getcellradians(cell::Cell)::Tuple{Float64, Float64}
    local lon::Float64 = (cell.lon * GRID_SIZE - GRID_SIZE / 2) * π/180
    local lat::Float64 = (cell.lat * GRID_SIZE - 90 - GRID_SIZE / 2) * π/180
    (lon, lat)
end

# Create a cell without generating the interpolation data
function _makecell(cellindex::Int64, lonsize::Int64, latsize::Int64)::Cell
    local lonindex::Int64 = floor(cellindex / latsize) + 1
    local latindex::Int64 = mod(cellindex, latsize)
    if latindex == 0
        latindex = 72
        lonindex -= 1
    end

    return (lon=lonindex, lat=latindex)
end

# Create an InterpolationCell
function _makeinterpolationcell(cell::Cell, weight::Union{Missing, Float64}, uncertainty::Union{Missing, Float64})::InterpolationCell
    (lon=cell.lon, lat=cell.lat, weight=weight, uncertainty=uncertainty)
end

end #module
