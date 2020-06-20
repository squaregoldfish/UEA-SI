using NCDatasets
using Plots
include("InterpolationData.jl")
using .InterpolationData

function missingtozero(series)
  series[ismissing.(series)] .= 0
end

FCO2_FILE = "daily.nc"

lons = Dataset(FCO2_FILE)["longitude"][:]
lats = Dataset(FCO2_FILE)["latitude"][:]
cells = makecells(length(lons), length(lats))

cell = cells[25]
data = _loadinterpolationdata(cell)

# Replace missings in uncertaintes with zeros
missingtozero(data.originalinputuncertainties)

timesteps = size(data.originalinputseries)[1]

interpolatedpoints = data.paraminputseries
interpolateduncertainties = data.paraminputuncertainties
originalindices = findall((!ismissing).(data.originalinputseries))

for i in originalindices
    interpolatedpoints[i] = missing
    interpolateduncertainties[i] = missing
end

missingtozero(interpolateduncertainties)

scatter(1:timesteps, data.originalinputseries, yerr=data.originalinputuncertainties, label="Original", size=(800,800), mc=:red, ms=8, msc=:red)
#scatter!(1:timesteps, interpolatedpoints, yerr=interpolateduncertainties, label="Interpolated", mc=:blue, ms=8)



#scatter(1:timesteps, [interpolatedpoints data.originalinputseries], size=(1800,1200),
#  markerstrokewidth=0, markersize=8, markercolor=[:blue :red],
#  label=["Interpolated" "Original"])

#scatter(1:timesteps, data.originalinputseries, size=(800,800), markerstrokewidth=0.1, markersize=8, markercolor=[:red], label="Original")