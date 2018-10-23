module InterpolationData
using Serialization
using ProgressMeter

export makeinterpolationbase

const INTERPOLATION_DATA_DIR = "interpolation_data"

struct InterpolationCellBase
    lonindex::Int64
    latindex::Int64
    land::Bool
end

mutable struct InterpolationCellData

    lonindex::Int64
    latindex::Int64
    finished::Bool
    fitparams::Array{Float64, 1}
    paraminputseries::Array{Float64, 1}
    paraminputuncertainties::Array{Float64, 1}

    function InterpolationCellData(lonindex::Int64, latindex::Int64, timesize::Int64)
        newobj::InterpolationCellData = new(lonindex, latindex)
        
        newobj.finished = false
        newobj.fitparams = Array{Float64, 1}(undef, 5)
        newobj.paraminputseries = Array{Float64, 1}(undef, timesize)
        newobj.paraminputuncertainties = Array{Float64, 1}(undef, timesize)

        return newobj
    end

end

# Create the array of base data structures for the interpolations
function makeinterpolationbase(lonsize::Int64, latsize::Int64, timesize::Int64, seamask::Array{Int64, 2})

    if isdir(INTERPOLATION_DATA_DIR)
        rm(INTERPOLATION_DATA_DIR, recursive=true)
    end

    mkdir(INTERPOLATION_DATA_DIR)

    local interpolationbase::Array{InterpolationCellBase, 1} = Array{InterpolationCellBase, 1}(undef, lonsize * latsize)

    local counter::Int64 = 0
    local prog::Progress = Progress(lonsize, 1, "Preparing data structures")
    for i in 1:lonsize
        for j in 1:latsize
            counter = counter + 1

            interpolationbase[counter] = InterpolationCellBase(i, j, seamask[i, j] == 1)
            interpolationdata::InterpolationCellData = InterpolationCellData(i, j, timesize)
            saveinterpolationdata(interpolationdata)
        end
        next!(prog)
    end
    finish!(prog)

end #makeinterpolationbase

# Generate the filename for and InterpolationCellData object
function getdatafilename(data::InterpolationCellData)
    return "$(INTERPOLATION_DATA_DIR)/$(data.lonindex)_$(data.latindex).jldata"
end

# Generate the filename for and InterpolationCellData object
function getdatafilename(base::InterpolationCellBase)
    return "$(INTERPOLATION_DATA_DIR)/$(base.lonindex)_$(base.latindex).jldata"
end

# Save an InterpolationCellData object to disk
function saveinterpolationdata(data::InterpolationCellData)
    local out::IOStream = open(getdatafilename(data), "w")
    serialize(out, data)
    close(out)
end

# Load an InterpolationCellData object
function loadinterpolationdata(base::InterpolationCellBase)
    local inchan::IOStream = open(getdatafilename(base), "r")
    local data::InterpolationCellData = deserialize(inchan)
    close(inchan)
    return data
end


end #module