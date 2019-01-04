module InterpolationData
using Distributed
using ProgressMeter
using Serialization
using Logging
using Statistics

export Cell
export makecells
export interpolatecell

const Cell = NamedTuple{(:lon, :lat), Tuple{Int64, Int64}}
const INTERPOLATION_DATA_DIR = "interpolation_data"

# Limits
const MAX_STDEV = 75.0
const MIN_TIME_SPAN = 1825 # 5 years
const MONTH_END_DAYS = [31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334, 365]
const MIN_POPULATED_MONTHS = 8

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

        return newobj
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
function makecells(lonsize::Int64, latsize::Int64, timesize::Int64, seamask::Array{Int8, 2},
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
function makeinterpolationdata(cellindex::Int64, lonsize::Int64, latsize::Int64, seamask::Array{Int8, 2}, timesize::Int64,
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
function interpolatecell(cell::Cell, interpolationstep::Int8)
    data::InterpolationCellData = _loadinterpolationdata(cell)
    
    if !data.finished
        local logger::Tuple{SimpleLogger, IOStream} = _makelogger(cell, interpolationstep)
        interpolate!(data, interpolationstep, logger[1])
        _closelogger(logger[2])
    end

    return data.finished
end

# Perform interpolation
function interpolate!(data::InterpolationCellData, step::Int8, logger::SimpleLogger)
    with_logger(logger) do
        
        local fitfound::Bool = false
        local continuefit::Bool = true
        local spatialinterpolationcount::Int8 = 0

        # Variables for processing (empty placeholders)
        local currentseries::Array{Union{Missing, Float64}, 1} = Array{Union{Missing, Float64}, 1}()
        local currentuncertainties::Array{Union{Missing, Float64}, 1} = currentseries
        local currentweights::Array{Union{Missing, Float64}, 1} = currentseries
        local currentcurve::Array{Union{Missing, Float64}, 1} = currentseries

        while continuefit
            local attemptcurvefit::Bool = true

            # Put together the series for the curve fit
            if spatialinterpolationcount <= 0
                # Just use the original series as is
                currentseries = deepcopy(data.originalinputseries)
                currentuncertainties = deepcopy(data.originalinputuncertainties)
                currentweights = deepcopy(data.originalinputweights)
            else
                # Need to do spatial interpolation
                println("Must do spatial interpolation.")
                exit()
            end

            if continuefit
                if attemptcurvefit
                    currentcurve = attemptfit(currentseries, currentweights, currentuncertainties)
                    if length(currentcurve) > 0
                        println("Fit achieved")
                    else
                        println("Fit failed")
                    end
                end
            end
        end
    end
end

function attemptfit(series::Array{Union{Missing, Float64}, 1}, weights::Array{Union{Missing, Float64}, 1}, uncertainties::Array{Union{Missing, Float64}, 1})::Array{Union{Missing, Float64}, 1}
    local successfulparams::Array{Float64, 1} = Array{Float64, 1}()

    local fittedparams::Array{Float64, 1} = fitcurve(series, weights)


end

function fitcurve(series::Array{Union{Missing, Float64}, 1}, weights::Array{Union{Missing, Float64}, 1})::Array{Float64, 1}
    local fitseries::Array{Union{Missing, Float64}, 1} = removeoutliers(series)
    if !doprefitcheck(fitseries)
        @info("PREFIT CHECKS FAILED")
    else
        println("I can continue")
    end

    exit()

end

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

    @debug("Removed $(_countvalues(series) - _countvalues(filteredseries)) outliers")
    filteredseries
end

function doprefitcheck(series::Array{Union{Missing, Float64}, 1})::Bool
    local ok::Bool = true

    # Standard deviation
    local stdev = std(skipmissing(series))
    if stdev > MAX_STDEV
        @info("Standard deviation is $stdev, should be <= $MAX_STDEV")
        ok = false
    end

    # Minimum time coverage
    local missingdays::Array{Int64, 1} = findall(!ismissing, series)
    local dayspan::Int64 = missingdays[end] - missingdays[1]
    if dayspan < MIN_TIME_SPAN
        @info("Measurements must span at least $MIN_TIME_SPAN days (only has $dayspan)")
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

        local populatedmonthcount::Int8 = sum(populatedmonths)
        if populatedmonthcount < MIN_POPULATED_MONTHS
            @info("Should have at least $MIN_POPULATED_MONTHS populated months (only has $populatedmonthcount)")
        end
    end

    ok
end

###############################################################
#
# Simple utility methods

# Initialise a logger for a cell and interpolation run
function _makelogger(cell::Cell, interpolationstep::Int8)::Tuple{SimpleLogger, IOStream}
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
