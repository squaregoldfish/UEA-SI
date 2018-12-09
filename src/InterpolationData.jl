module InterpolationData
using Serialization
using ProgressMeter
using Distributed

export Cell
export makecells

const Cell = NamedTuple{(:lon, :lat), Tuple{Int64, Int64}}
const INTERPOLATION_DATA_DIR = "interpolation_data"

mutable struct InterpolationCellData
    cell::Cell
    finished::Bool
    fitparams::Array{Float64, 1}
    paraminputseries::Array{Float64, 1}
    paraminputuncertainties::Array{Float64, 1}

    function InterpolationCellData(cell::Cell, timesize::Int64)
        newobj::InterpolationCellData = new(cell)
        
        newobj.finished = false
        newobj.fitparams = Array{Float64, 1}(undef, 5)
        newobj.paraminputseries = Array{Float64, 1}(undef, timesize)
        newobj.paraminputuncertainties = Array{Float64, 1}(undef, timesize)

        return newobj
    end
end #InterpolationCellData

# Create the array of base data structures for the interpolations
function makecells(lonsize::Int64, latsize::Int64, timesize::Int64, seamask::Array{Int64, 2})::Array{Cell, 1}

    if isdir(INTERPOLATION_DATA_DIR)
        rm(INTERPOLATION_DATA_DIR, recursive=true)
    end

    mkdir(INTERPOLATION_DATA_DIR)

    println("Initialising working data...")
    local cells::Array{Cell, 1} = @showprogress pmap((x) -> makeinterpolationdata(x, lonsize, latsize, timesize), (range(1, length=lonsize*latsize)))

    return cells
end #makecells

function makeinterpolationdata(cellindex::Int64, lonsize::Int64, latsize::Int64, timesize::Int64)
    local lonindex::Int64 = floor(cellindex / latsize) + 1
    local latindex::Int64 = mod(cellindex, latsize)
    if latindex == 0
        latindex = 72
        lonindex -= 1
    end

    local cell::Cell = makecell(lonindex, latindex)
    local interpolationdata::InterpolationCellData = InterpolationCellData(cell, timesize)
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