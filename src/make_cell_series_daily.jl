using NCDatasets

const INFILE = "SOCATv5.tsv"
const OUTFILE = "daily.nc"
const CELLSIZE = 2.5
const STARTYEAR = 1985
const ENDYEAR = 2016

# Zero-based day of year of each month
const MONTHSTARTS = [0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334]
const LEAPMONTHSTARTS = [0, 31, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335]

# Days of the year
DAYSTARTS = range(1, step=1, length=365)

# Calculate the conversion of leap year days to 365 'normal' year days
# Each leap year day lasts 1 1/365 normal days 
LEAPDAYSTARTS = Array{Float64}(undef, 365)

for i in 1:365
    LEAPDAYSTARTS[i] = i + (((366 / 365) / 365) * (i - 1))
end

# Determine whether or not a given year is a leap year
function isleapyear(year::Int64)::Bool
    local leapyear::Bool = false
    if rem(year, 4) == 0
        leapyear = true
        if rem(year, 100) == 0
            if rem(year, 400) != 0
                leapyear = false
            end
        end
    end

    leapyear
end

# Calculate the Julian day from a year/month/day set
function calcjdate(year::Int64, month::Int64, day::Int64)::Int64
    local result::Int64 = 0

    if isleapyear(year)
        result = LEAPMONTHSTARTS[month]
    else
        result = MONTHSTARTS[month]
    end

    result = result + day
end

# Calculate the day of the year for a given date
function getdayindex(date, days)::Float64
    return findall(days .<= date)[end]
end

# Calcaulte the nth day of the complete data set
# Days are only calculated between the start and end years
function getdateindex(year::Int64, month::Int64, day::Int64)::Int64
    
    local jdate::Float64 = calcjdate(year, month, day)    
    local index::Float64 = -1
    local dayindex::Float64 = 0

    if year >= STARTYEAR && year <= ENDYEAR
        if isleapyear(year)
            dayindex = getdayindex(jdate, DAYSTARTS)
        else
            dayindex = getdayindex(jdate, LEAPDAYSTARTS)
        end
    
        index = ((year - STARTYEAR) * 365) + dayindex
    end

    return convert(Int64, floor(index))
end

# Calculate the cell indices from a longitude and latitude
function getcellindex(longitude::Float64, latitude::Float64)::Tuple{Int64, Int64}
    local loncell::Int64 = floor(longitude / CELLSIZE) + 1
    local latcell::Int64 = floor((latitude + 90) / CELLSIZE) + 1

    if loncell == 145
        loncell = 1
    end

    return loncell, latcell
end

function applydataset(datasettotals::Array{Float64, 3}, datasetcounts::Array{Int64, 3}, overalltotals::Array{Float64, 3}, overallcounts::Array{Int64, 3})
    datasetmean::Array{Float64, 3} = datasettotals ./ datasetcounts
    datasetmean[isnan.(datasetmean)] .= 0

    overalltotals .= overalltotals .+ datasetmean
    datasetmean[datasetmean .> 0] .= 1
    overallcounts .= overallcounts .+ datasetmean
end

function run()

    totaldays::Int64 = (ENDYEAR - STARTYEAR + 1) * 365

    # Output data set
    overallcelltotals::Array{Float64, 3} = zeros(convert(Int64, 360 / CELLSIZE), convert(Int64, 180 / CELLSIZE), totaldays)
    overallcellcounts::Array{Int64, 3} = zeros(convert(Int64, 360 / CELLSIZE), convert(Int64, 180 / CELLSIZE), totaldays)


    # Open input file
    inchan::IOStream = open(INFILE)

    currentdataset::String = ""
    datasetcelltotals::Array{Float64, 3} = zeros(convert(Int64, 360 / CELLSIZE), convert(Int64, 180 / CELLSIZE), totaldays)
    datasetcellcounts::Array{Int64, 3} = zeros(convert(Int64, 360 / CELLSIZE), convert(Int64, 180 / CELLSIZE), totaldays)

    currentline::String = readline(inchan)
    linecount::Int64 = 1
    while length(currentline) > 0

        fields::Array{String, 1} = split(currentline, "\t")

        dataset::String = fields[1]

        if dataset != currentdataset
            if length(currentdataset) > 0
                applydataset(datasetcelltotals, datasetcellcounts, overallcelltotals, overallcellcounts)

                print("\033[1K\r$currentdataset ($linecount)")

                datasetcelltotals .= 0
                datasetcellcounts .= 0
            end

            currentdataset = dataset
        end

        year::Int64 = parse(Int64, fields[5])
        month::Int64 = parse(Int64, fields[6])
        day::Int64 = parse(Int64, fields[7])
        longitude::Float64 = parse(Float64, fields[11])
        latitude::Float64 = parse(Float64, fields[12])
        fco2::Float64 = parse(Float64, fields[24])

        cellindex::Tuple{Int64, Int64} = getcellindex(longitude, latitude)

        dateindex::Int64 = getdateindex(year, month, day)

        if dateindex != -1
            datasetcelltotals[cellindex[1], cellindex[2], dateindex] = datasetcelltotals[cellindex[1], cellindex[2], dateindex] + fco2
            datasetcellcounts[cellindex[1], cellindex[2], dateindex] = datasetcellcounts[cellindex[1], cellindex[2], dateindex] + 1
        end

        currentline = readline(inchan)
        linecount = linecount + 1
    end

    close(inchan)

    # The last dataset
    applydataset(datasetcelltotals, datasetcellcounts, overallcelltotals, overallcellcounts)

    # Overall cell means
    local meanfco2::Array{Float64, 3} = overallcelltotals ./ overallcellcounts
    meanfco2 = replace(meanfco2, NaN=>-1e35)

    # Write NetCDF
    nc = Dataset(OUTFILE, "c")
    defDim(nc, "longitude", trunc(Int, (360 / CELLSIZE)))
    defDim(nc, "latitude", trunc(Int, (180 / CELLSIZE)))
    defDim(nc, "time", totaldays)

    nclon = defVar(nc, "longitude", Float32, ("longitude",))
    nclon.attrib["units"] = "degrees_east"

    nclat = defVar(nc, "latitude", Float32, ("latitude",))
    nclat.attrib["units"] = "degrees_north"

    nctime = defVar(nc, "time", Float32, ("time",))
    nctime.attrib["calendar"] = "noleap"

    ncfco2 = defVar(nc, "fCO2", Float64, ("longitude", "latitude", "time"))
    ncfco2.attrib["_FillValue"] = -1e35

    nclon[:] = collect(range(CELLSIZE / 2, step=CELLSIZE, stop=360))
    nclat[:] = collect(range(-90 + CELLSIZE / 2, step=CELLSIZE, stop=90))
    nctime[:] = collect(range(STARTYEAR, step=(1/365), stop=(ENDYEAR + 1) - (1/365)))
    ncfco2[:,:,:] = meanfco2

    close(nc)
end

@time run()
