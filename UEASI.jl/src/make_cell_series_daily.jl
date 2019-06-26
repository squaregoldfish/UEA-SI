using NCDatasets
using ProgressMeter

const INFILES = ["SOCATv6.tsv" "SOCATv6_FlagE.tsv"]
const UNCERTAINTIES = [2.5 10]
const OUTFILE = "daily.nc"
const CELLSIZE = 2.5
const STARTYEAR = 1985
const ENDYEAR = 2017

# Zero-based day of year of each month
const MONTHSTARTS = [0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334]
const LEAPMONTHSTARTS = [0, 31, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335]

# Days of the year
const DAYSTARTS = range(1, step=1, length=365)

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
function getdayindex(date, days)::Int64
    return findall(days .<= date)[end]
end

# Calcaulte the nth day of the complete data set
# Days are only calculated between the start and end years
function getdateindex(year::Int64, month::Int64, day::Int64)::Int64

    local jdate::Int64 = calcjdate(year, month, day)
    local index::Int64 = -1
    local dayindex::Int64 = 0

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

function run()

    local totaldays::Int64 = (ENDYEAR - STARTYEAR + 1) * 365

    # Output data set
    print("Initialising...")
    local overallcelltotals::Array{Float64, 3} = zeros(convert(Int64, 360 / CELLSIZE), convert(Int64, 180 / CELLSIZE), totaldays)
    local overalluncertaintytotals::Array{Float64, 3} = zeros(convert(Int64, 360 / CELLSIZE), convert(Int64, 180 / CELLSIZE), totaldays)
    local overallcellcounts::Array{Int64, 3} = zeros(convert(Int64, 360 / CELLSIZE), convert(Int64, 180 / CELLSIZE), totaldays)

    # Loop over all files
    for fileindex = 1:length(INFILES)
        local file::String = INFILES[fileindex]
        local uncertainty::Float64 = UNCERTAINTIES[fileindex]

        # Open input file
        local inchan::IOStream = open(file)

        seekend(inchan)
        local chansize::Int64 = position(inchan)
        seekstart(inchan)

        p::Progress = Progress(chansize, 1, "Reading $file")

        # Skip the header
        local inheader::Bool = true
        while inheader
            if findfirst(r"Expocode.*SOCAT_DOI", readline(inchan)) !== nothing
                inheader = false
            end
        end

        update!(p, position(inchan))

        local currentdataset::String = ""
        local currentcell::Tuple{Int64, Int64} = -1, -1
        local currentdate::Int = -1

        local currenttotal::Float64 = 0
        local currentuncertaintytotal::Float64 = 0
        local currentcount::Int64 = 0

        local currentline::String = readline(inchan)
        update!(p, position(inchan))

        while length(currentline) > 0
            local fields::Array{String, 1} = split(currentline, "\t")

            local dataset::String = fields[1]
            local year::Int64 = parse(Int64, fields[5])
            local month::Int64 = parse(Int64, fields[6])
            local day::Int64 = parse(Int64, fields[7])
            local dateindex::Int64 = getdateindex(year, month, day)

            local longitude::Float64 = parse(Float64, fields[11])
            local latitude::Float64 = parse(Float64, fields[12])
            local cellindex::Tuple{Int64, Int64} = getcellindex(longitude, latitude)

            local fco2::Float64 = parse(Float64, fields[24])

            if dataset != currentdataset ||
                cellindex != currentcell ||
                dateindex != currentdate

                if currentdate != -1
                    @inbounds overallcelltotals[currentcell[1], currentcell[2], currentdate] =
                        overallcelltotals[currentcell[1], currentcell[2], currentdate] + currenttotal / currentcount

                    @inbounds overalluncertaintytotals[currentcell[1], currentcell[2], currentdate] =
                        overalluncertaintytotals[currentcell[1], currentcell[2], currentdate] + currentuncertaintytotal / currentcount

                    @inbounds overallcellcounts[currentcell[1], currentcell[2], currentdate] =
                        overallcellcounts[currentcell[1], currentcell[2], currentdate] + 1
                end

                currentdataset = dataset
                currentcell = cellindex
                currentdate = dateindex

                currenttotal = 0
                currentuncertaintytotal = 0
                currentcount = 0
            end


            if dateindex != -1
                currenttotal = currenttotal + fco2
                currentuncertaintytotal = currentuncertaintytotal + uncertainty
                currentcount = currentcount + 1
            end

            currentline = readline(inchan)
            update!(p, position(inchan))
        end

        close(inchan)

        # The last dataset
        update!(p, position(inchan))
        if currentdate != -1
            @inbounds overallcelltotals[currentcell[1], currentcell[2], currentdate] =
                overallcelltotals[currentcell[1], currentcell[2], currentdate] + currenttotal / currentcount

            @inbounds overalluncertaintytotals[currentcell[1], currentcell[2], currentdate] =
                overalluncertaintytotals[currentcell[1], currentcell[2], currentdate] + currentuncertaintytotal / currentcount

            @inbounds overallcellcounts[currentcell[1], currentcell[2], currentdate] =
                overallcellcounts[currentcell[1], currentcell[2], currentdate] + 1
        end

        finish!(p)
    end

    # Overall cell means
    print("\033[1K\rCalculating means...")
    local meanfco2::Array{Float64, 3} = overallcelltotals ./ overallcellcounts
    meanfco2 = replace(meanfco2, NaN=>-1e35)

    local meanuncertainty::Array{Float64, 3} = overalluncertaintytotals ./ overallcellcounts
    meanuncertainty = replace(meanuncertainty, NaN=>-1e35)

    # Write NetCDF
    print("\033[1K\rWriting output...")
    local nc = Dataset(OUTFILE, "c")
    defDim(nc, "longitude", trunc(Int, (360 / CELLSIZE)))
    defDim(nc, "latitude", trunc(Int, (180 / CELLSIZE)))
    defDim(nc, "time", totaldays)

    local nclon = defVar(nc, "longitude", Float32, ("longitude",))
    nclon.attrib["units"] = "degrees_east"

    local nclat = defVar(nc, "latitude", Float32, ("latitude",))
    nclat.attrib["units"] = "degrees_north"

    local nctime = defVar(nc, "time", Float32, ("time",))
    nctime.attrib["units"] = "year"

    local ncfco2 = defVar(nc, "fCO2", Float64, ("longitude", "latitude", "time"))
    ncfco2.attrib["_FillValue"] = -1e35

    local ncuncertainty = defVar(nc, "uncertainty", Float64, ("longitude", "latitude", "time"))
    ncuncertainty.attrib["_FillValue"] = -1e35

    nclon[:] = collect(range(CELLSIZE / 2, step=CELLSIZE, stop=360))
    nclat[:] = collect(range(-90 + CELLSIZE / 2, step=CELLSIZE, stop=90))
    nctime[:] = collect(range(STARTYEAR, step=(1/365), stop=(ENDYEAR + 1) - (1/365)))
    ncfco2[:,:,:] = meanfco2
    ncuncertainty[:,:,:] = meanuncertainty

    close(nc)
    print("\n")
end

@time run()
