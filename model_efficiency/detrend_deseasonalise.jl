using NetCDF

function run()
    nanvalue = parse(Float64, ARGS[6])

    nc = NetCDF.open(ARGS[1], readdimvar=true)

    lons = nc[ARGS[2]]
    lats = nc[ARGS[3]]
    times = nc[ARGS[4]]
    fco2 = nc[ARGS[5]]

    localfco2 = copy(fco2)

    for i in eachindex(localfco2)
        if localfco2[i] == nanvalue
            localfco2[i] = NaN
        end
    end

    NetCDF.ncclose(ARGS[1])
end

@time run()