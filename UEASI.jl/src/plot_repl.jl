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

cell = cells[5747]
data = _loadinterpolationdata(cell)

# Replace missings in uncertaintes with zeros
missingtozero(data.originalinputuncertainties)

timesteps = size(data.originalinputseries)[1]

# Make two copies of the param input series.
# The copy for temporal interpolations does not need uncertainties
temporalinterpolatedpoints = copy(data.paraminputseries)
spatialinterpolatedpoints = copy(data.paraminputseries)
spatialinterpolateduncertainties = copy(data.paraminputuncertainties)

originalindices = findall((!ismissing).(data.originalinputseries))

for i in 1:timesteps
  if !ismissing(data.originalinputseries[i])
    temporalinterpolatedpoints[i] = missing
    spatialinterpolatedpoints[i] = missing
    spatialinterpolateduncertainties[i] = missing
  elseif !ismissing(spatialinterpolatedpoints[i])

    # If the interpolated uncertainty is 0.001, then it's from a temporal
    # interpolation. Otherwise it's a spatial one.
    if spatialinterpolateduncertainties[i] == 0.001
      spatialinterpolatedpoints[i] = missing
      spatialinterpolateduncertainties[i] = missing
    else
      temporalinterpolatedpoints[i] = missing
    end
  end
end


missingtozero(spatialinterpolateduncertainties)

if data.finished

  plot()

  if sum((!ismissing).(temporalinterpolatedpoints)) > 0
    scatter!(1:timesteps, temporalinterpolatedpoints, label="Temporal Interpolated", seriescolor=RGB(.75, .75, 1), ms=8, msw=0.1, size=(2400,800))
  end

  if sum((!ismissing).(spatialinterpolatedpoints)) > 0
    scatter!(1:timesteps, spatialinterpolatedpoints, yerr=spatialinterpolateduncertainties, label="Spatial Interpolated", seriescolor=RGB(0, 0, 1), ms=8, msw=0.1, size=(2400,800))
  end

  if sum((!ismissing).(data.originalinputseries)) > 0
      scatter!(1:timesteps, data.originalinputseries, yerr=data.originalinputuncertainties, label="Original", c=RGB(1, 0, 0), ms=8, msw=0.1, size=(2400,800))
  end

  fittedcurve = makecurve(data.fitparams, timesteps)
  plot!(fittedcurve, label="Fitted curve", seriescolor=RGB(.5, .5, .5))
end