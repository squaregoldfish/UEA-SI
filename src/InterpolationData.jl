module InterpolationData
using Distributed
using ProgressMeter
using Serialization

export Cell
export makecells

const Cell = NamedTuple{(:lon, :lat), Tuple{Int64, Int64}}
const INTERPOLATION_DATA_DIR = "interpolation_data"

mutable struct InterpolationCellData
    cell::Cell
    land::Bool
    finished::Bool
    fitparams::Array{Float64, 1}
    originalinputseries::Array{Union{Missing, Float64}, 1}
    originalinputuncertainties::Array{Union{Missing, Float64}, 1}
    paraminputseries::Array{Union{Missing, Float64}, 1}
    paraminputuncertainties::Array{Union{Missing, Float64}, 1}

    function InterpolationCellData(cell::Cell, timesize::Int64, island::Bool, fco2::Array{Union{Missing, Float64}, 1},
        uncertainty::Array{Union{Missing, Float64}, 1})
        
        newobj::InterpolationCellData = new(cell)
        
        newobj.land = island
        newobj.finished = island # If the cell is on land, we've already finished
        newobj.fitparams = Array{Float64, 1}(undef, 5)
        newobj.originalinputseries = fco2
        newobj.originalinputuncertainties = uncertainty
        newobj.paraminputseries = newobj.originalinputseries
        newobj.paraminputuncertainties = newobj.originalinputuncertainties

        return newobj
    end
end #InterpolationCellData

# Create the array of base data structures for the interpolations
function makecells(lonsize::Int64, latsize::Int64, timesize::Int64, seamask::Array{Int8, 2},
    fco2::Array{Union{Missing, Float64}, 3}, uncertainty::Array{Union{Missing, Float64}, 3})::Array{Cell, 1}

    #if isdir(INTERPOLATION_DATA_DIR)
    #    rm(INTERPOLATION_DATA_DIR, recursive=true)
    #end

    #mkdir(INTERPOLATION_DATA_DIR)

    local cells::Array{Cell, 1} = Array{Cell, 1}(undef, lonsize * latsize)
    @showprogress 1 "Initialising working data" for i in 1:(lonsize * latsize)
        cells[i] = makeinterpolationdata(i, lonsize, latsize, seamask, timesize, fco2, uncertainty)
    end

    return cells
end #makecells

function makeinterpolationdata(cellindex::Int64, lonsize::Int64, latsize::Int64, seamask::Array{Int8, 2}, timesize::Int64,
    fco2::Array{Union{Missing, Float64}, 3}, uncertainty::Array{Union{Missing, Float64}, 3})

    local lonindex::Int64 = floor(cellindex / latsize) + 1
    local latindex::Int64 = mod(cellindex, latsize)
    if latindex == 0
        latindex = 72
        lonindex -= 1
    end

    local cell::Cell = makecell(lonindex, latindex)
    local land::Bool = seamask[cell.lon, cell.lat] == 0
    local interpolationdata::InterpolationCellData = InterpolationCellData(cell, timesize, land, fco2[cell.lon, cell.lat, :], uncertainty[cell.lon, cell.lat, :])
    #saveinterpolationdata(interpolationdata)

    return cell
end

# Generate the filename for and InterpolationCellData object
function getdatafilename(data::InterpolationCellData)::String
    return "$(INTERPOLATION_DATA_DIR)/$(data.cell.lon)_$(data.cell.lat).jldata"
end

# Save an InterpolationCellData object to disk
function saveinterpolationdata(data::InterpolationCellData)
    local out::IOStream = open(getdatafilename(data), "w")
    serialize(out, data)
    close(out)
end

# Load an InterpolationCellData object
function loadinterpolationdata(cell::Cell)::InterpolationCellData
    local inchan::IOStream = open(getdatafilename(base), "r")
    local data::InterpolationCellData = deserialize(inchan)
    close(inchan)
    return data
end

function makecell(lonindex::Int64, latindex::Int64)::Cell
    (lon=lonindex, lat=latindex)
end
end #module