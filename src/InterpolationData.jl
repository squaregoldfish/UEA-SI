module InterpolationData
using Distributed
using ProgressMeter
using Serialization
using Logging
using Statistics
using LsqFit

export Cell
export makecells
export interpolatecell
export InterpolationCellData


const Cell = NamedTuple{(:lon, :lat), Tuple{UInt16, UInt16}}
const INTERPOLATION_DATA_DIR = "interpolation_data"

mutable struct InterpolationCellData
    cell::Cell
    land::Bool
    finished::Bool
    fitparams::Array{Float64, 1}
    originalinputseries::Array{Union{Missing, Float64}, 1}
    originalinputuncertainties::Array{Union{Missing, Float64}, 1}
    originalinputweights::Array{Union{Missing, Float64}, 1}
    paraminputseries::Array{Union{Missing, Float64}, 1}
    paraminputuncertainties::Array{Union{Missing, Float64}, 1}
    paraminputweights::Array{Union{Missing, Float64}, 1}

    function InterpolationCellData(cell::Cell, timesize::Int64, island::Bool, fco2::Array{Union{Missing, Float64}, 1},
        uncertainty::Array{Union{Missing, Float64}, 1})

        newobj::InterpolationCellData = new(cell)

        newobj.land = island
        newobj.finished = island # If the cell is on land, we've already finished
        newobj.fitparams = Array{Float64, 1}(undef, 5)

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
function makecells(lonsize::Int64, latsize::Int64)::Array{Cell, 1}
    local cells::Array{Cell, 1} = Array{Cell, 1}(undef, lonsize * latsize)
    for i in 1:(lonsize * latsize)
        cells[i] = _makecell(i, lonsize, latsize)
    end

    return cells
end

# Create the array of base data structures for the interpolations
function makecells(lonsize::Int64, latsize::Int64, timesize::Int64, seamask::Array{UInt8, 2},
    fco2::Array{Union{Missing, Float64}, 3}, uncertainty::Array{Union{Missing, Float64}, 3})::Array{Cell, 1}

    if isdir(INTERPOLATION_DATA_DIR)
        rm(INTERPOLATION_DATA_DIR, recursive=true)
    end

    mkdir(INTERPOLATION_DATA_DIR)

    local cells::Array{Cell, 1} = Array{Cell, 1}(undef, lonsize * latsize)
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
    local interpolationdata::InterpolationCellData = InterpolationCellData(cell, timesize, land, fco2[cell.lon, cell.lat, :], uncertainty[cell.lon, cell.lat, :])
    _saveinterpolationdata(interpolationdata)

    return cell
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

# Perform the main interpolation for a cell
function interpolatecell(cell::Cell, interpolationstep::UInt8)
    data::InterpolationCellData = _loadinterpolationdata(cell)

    if !data.finished
        local logger::Tuple{SimpleLogger, IOStream} = _makelogger(cell, interpolationstep)
        interpolate!(data, interpolationstep, logger[1])
        _closelogger(logger[2])
    end

    return data.finished
end


###############################################################
#
# Interpolation stuff

# Limits
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


# Structure to hold a series and its related data
mutable struct SeriesData
    measurements::Array{Union{Missing, Float64}, 1}
    uncertainties::Array{Union{Missing, Float64}, 1}
    weights::Array{Union{Missing, Float64}, 1}
    curve::Array{Float64, 1}
    curveparams::Array{Float64, 1}

    function SeriesData()
        newobj::SeriesData = new()
        newobj.measurements = Array{Union{Missing, Float64}, 1}()
        newobj.uncertainties = Array{Union{Missing, Float64}, 1}()
        newobj.weights = Array{Union{Missing, Float64}, 1}()
        newobj.curve = Array{Float64, 1}()
        newobj.curveparams = Array{Float64, 1}()
        newobj
    end

    function SeriesData(measurements::Array{Union{Missing, Float64}, 1}, uncertainties::Array{Union{Missing, Float64}, 1}, weights::Array{Union{Missing, Float64}, 1})
        newobj::SeriesData = new(deepcopy(measurements), deepcopy(uncertainties), deepcopy(weights))
        newobj.curve = Array{Float64, 1}()
        newobj.curveparams = Array{Float64, 1}()
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
function interpolate!(data::InterpolationCellData, step::UInt8, logger::SimpleLogger)
    with_logger(logger) do

        local fitfound::Bool = false
        local continuefit::Bool = true

        local spatialinterpolationcount::UInt8 = 0
        local spatialinterpolationcells::Array{Cell, 1} = Array{Cell, 1}()

        # Variables for processing (empty placeholders)
        local currentseries::SeriesData = SeriesData()

        # Temporary storage used while trying extra interpolations
        local storedseries::SeriesData = SeriesData()

        while continuefit
            local attemptcurvefit::Bool = true

            # Put together the series for the curve fit
            if spatialinterpolationcount == 0
                # Just use the original series as is
                currentseries.measurements = deepcopy(data.originalinputseries)
                currentseries.uncertainties = deepcopy(data.originalinputuncertainties)
                currentseries.weights = deepcopy(data.originalinputweights)
            else
                # Need to do spatial interpolation
                println("Must do spatial interpolation.")
                exit()
            end

            if continuefit
                if attemptcurvefit
                    attemptfit!(currentseries)

                    if length(currentseries.curve) > 0
                        fitfound = true

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
                            else
                                spatialinterpolationcount += 1
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
                                    else
                                        spatialinterpolationcount += 1
                                    end
                                else
                                    spatialinterpolationcount += 1
                                end
                            else
                                # The spatial interpolation made no difference. We can stop
                                # with the previous fit
                                continuefit = false
                            end
                        end
                    end
                end
            end

            # Safety cutoff
            exit()
        end
    end
end

# Try to fit a curve to the supplied time series
function attemptfit!(series::SeriesData)

    # Attempt a curve fit
    fitcurve!(series)

    if length(series.curveparams) == 0
        println("Interpolate")
    end
end

# Fit a harmonic curve to a time series
function fitcurve!(series::SeriesData)
    local fitseries::Array{Union{Missing, Float64}, 1} = removeoutliers(series.measurements)
    if !doprefitcheck(fitseries)
        @info "PREFIT CHECKS FAILED"
    else
        local fitsuccess::Bool = false
        local harmoniccount::UInt8 = 4

        while !fitsuccess && harmoniccount > 0
            @debug "Fitting $harmoniccount harmonic(s)"

            local termcount::UInt8 = 0

            local formula::String = "p[1] + p[2]x"
            termcount = 2

            for i in 1:harmoniccount
                formula = "$formula + p[$(termcount + 1)]*sin(2*π*$i*(x/365)) + p[$(termcount + 2)]*cos(2*π*$i*(x/365))"
                termcount += 2
            end

            fitfunction = eval(Meta.parse("@. (x, p) -> " * formula))
            model(x, p) = Base.invokelatest(fitfunction, x, p)

            local p0::Array{Float64, 1} = zeros(termcount)
            local days::Array{UInt16, 1} = findall((!ismissing).(fitseries))

            # TODO Handle fit failure
            local fit::LsqFit.LsqFitResult = curve_fit(model, days, collect(skipmissing(fitseries)), p0)
            local fitparams::Array{Float64, 1} = fit.param
            local fittedcurve::Array{Float64, 1} = makecurve(fitparams, length(fitseries))
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
end

# Remove outliers from a time series. The outlier limit
# is defined in MAX_STDEV
function removeoutliers(series::Array{Union{Missing, Float64}, 1})::Array{Union{Missing, Float64}, 1}
    local filteredseries::Array{Union{Missing, Float64}, 1} = Array{Union{Missing, Float64}, 1}(missing, length(series))

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
function doprefitcheck(series::Array{Union{Missing, Float64}, 1})::Bool
    local ok::Bool = true

    # Standard deviation
    local stdev = std(skipmissing(series))
    @debug "Series standard deviation = $stdev"
    if stdev > MAX_STDEV
        @info "Standard too large, should be ≤ $MAX_STDEV"
        ok = false
    end

    # Minimum time coverage
    local missingdays::Array{Int64, 1} = findall(!ismissing, series)
    local dayspan::Int64 = missingdays[end] - missingdays[1]
    @debug "Measurements span $dayspan days"
    if dayspan < MIN_TIME_SPAN
        @info "Measurements must span at least $MIN_TIME_SPAN days"
        ok = false
    end

    # Populated months
    if ok
        local populatedmonths::Array{Bool, 1} = falses(12)
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
        end
    end

    ok
end

# Make a curve using the supplied parameters
function makecurve(params::Array{Float64, 1}, curvelength::Int64)::Array{Float64, 1}

    local curve::Array{Float64, 1} = zeros(curvelength)

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
function makeseasonalcycle(curveparams::Array{Float64, 1})::Array{Float64, 1}

    local seasonalcycle::Array{Float64, 1} = zeros(365)

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
function checkcurvepeaks(curveparams::Array{Float64, 1})::Bool

    local peaksok::Bool = true

    local seasonalcycle::Array{Float64, 1} = makeseasonalcycle(curveparams)

    local maxima::Array{Float64, 1} = Array{Float64, 1}()
    local maximapos::Array{Int16, 1} = Array{Int16, 1}()
    local minima::Array{Float64, 1} = Array{Float64, 1}()
    local minimapos::Array{Int16, 1} = Array{Int16, 1}()

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
                startdirection = slope_direction
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
        local peaksizes::Array{Float64, 1} = zeros(length(maxima))

        for i in 1:length(maxima)
            local maxvalue::Float64 = maxima[i]
            local maxpos::UInt16 = maximapos[i]

            # For each peak, we find the position of the preceding minimum
            # We loop round to the end of the year if necessary
            local minentry::UInt16 = 0
            local minposcandidates::Array{UInt16, 1} = findall(x -> (x < maxpos), minimapos)
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

function checkcurvefit(series::Array{Union{Missing, Float64}, 1}, fittedcurve::Array{Float64, 1})::Bool

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
        curveok = true
        @info "Ratio FAILED should be $MIN_CURVE_RATIO  to $MAX_CURVE_RATIO"
    end

    local maxdiff::Float64 = abs(seriesmax - curvemax)
    local mindiff::Float64 = abs(seriesmin - curvemin)

    @debug "Limit differences = $maxdiff, $mindiff"
    if maxdiff > MAX_LIMIT_DIFFERENCE || mindiff > MAX_LIMIT_DIFFERENCE
        curveok = false
        @info "Limit difference FAILED should be in be $MAX_LIMIT_DIFFERENCE or less"
    end

    curveok
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
function _loadinterpolationdata(cell::Cell)::InterpolationCellData
    local inchan::IOStream = open(_getdatafilename(cell), "r")
    local data::InterpolationCellData = deserialize(inchan)
    close(inchan)
    return data
end

# Count the number of non-missing values in a series
function _countvalues(series::Array{Union{Missing, Float64}, 1})::Int64
    length(collect(skipmissing(series)))
end

end #module
