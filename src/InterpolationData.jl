module InterpolationData
using Serialization
using ProgressMeter

export Cell
export makecells

const Cell = Pair{Int64, Int64}
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

    function lonindex()
        cell[0]
    end

    function latindex()
        cell[1]
    end
end #InterpolationCellData

# Create the array of base data structures for the interpolations
function makecells(lonsize::Int64, latsize::Int64, timesize::Int64, seamask::Array{Int64, 2})::Array{Pair, 1}

    if isdir(INTERPOLATION_DATA_DIR)
        rm(INTERPOLATION_DATA_DIR, recursive=true)
    end

    mkdir(INTERPOLATION_DATA_DIR)

    local cells::Array{Cell, 1} = Array{Cell, 1}(undef, lonsize * latsize)

    local counter::Int64 = 0
    local prog::Progress = Progress(lonsize, 1, "Preparing data structures")
    for i in 1:lonsize
        for j in 1:latsize
            counter = counter + 1

            cells[counter] = Cell(i, j)
            interpolationdata::InterpolationCellData = InterpolationCellData(cells[counter], timesize)
            saveinterpolationdata(interpolationdata)
        end
        next!(prog)
    end
    finish!(prog)

    return cells
end #makeinterpolationbase

# Generate the filename for and InterpolationCellData object
function getdatafilename(data::InterpolationCellData)::String
    return "$(INTERPOLATION_DATA_DIR)/$(lonindex(data))_$(latindex(data)).jldata"
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

function lonindex(data::InterpolationCellData)::Int64
    data.cell.first
end

function latindex(data::InterpolationCellData)::Int64
    data.cell.second
end

end #module