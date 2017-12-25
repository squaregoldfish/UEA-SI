# Take a full pCO2 field and sample it to match SOCAT
using NetCDF

if size(ARGS) != (13,)
    println("Usage: julia sample_full_field.jl [start_year] [full_field] [var] [start] [end] [leap_years?] [sample_map] [var] [start] [end] [leap_years?] [lon_var] [lat_var]")
    quit()
end

start_year = parse(Int64, ARGS[1])
full_file = ARGS[2]
full_var_name = ARGS[3]
full_start = parse(Int64, ARGS[4])
full_end = parse(Int64, ARGS[5])
full_leap_years = parse(Bool, ARGS[6])
sample_file = ARGS[7]
sample_var_name = ARGS[8]
sample_start = parse(Int64, ARGS[9])
sample_end = parse(Int64, ARGS[10])
sample_leap_years = parse(Bool, ARGS[11])
lon_var = ARGS[12]
lat_var = ARGS[13]

if !isfile(full_file)
    println("Cannot find file $full_file")
    quit()
end

if !isfile(sample_file)
    println("Cannot find file $sample_file")
    quit()
end

function run()
    
    # Open the two files
    full_nc::NcFile = NetCDF.open(full_file)
    sample_nc::NcFile = NetCDF.open(sample_file)

    # Get the file dimensions
    full_size = size(full_nc[full_var_name])
    sample_size = size(sample_nc[sample_var_name])

    # Check that the grid sizes match
    if full_size[1] != sample_size[1] || full_size[2] != sample_size[2]
        error("Grid sizes do not match")
    end

    # Calculate the final number of time steps
    time_steps::Integer = calc_time_steps(sample_size[3], start_year, sample_leap_years)


    # Create the output file
    lons = NetCDF.readvar(sample_nc, lon_var)
    lats = NetCDF.readvar(sample_nc, lat_var)
    times = build_times(start_year, time_steps)

    lon_dim::NcDim = NetCDF.NcDim(lon_var, length(lons), values=lons)
    lat_dim::NcDim = NetCDF.NcDim(lat_var, length(lats), values=lats)
    time_dim::NcDim = NetCDF.NcDim("time", length(times), values=times)

    output_var::NcVar = NetCDF.NcVar(full_var_name, [lon_dim, lat_dim, time_dim], t=Float64)


    output_filename::String = full_file[1:length(full_file) - 3] * "_sampled.nc"
    output_nc::NcFile = NetCDF.create(output_filename, output_var)

    # Step through both files
    full_var = full_nc[full_var_name]
    sample_var = sample_nc[sample_var_name]

    current_year::Integer = start_year
    full_step::Integer = full_start
    sample_step::Integer = sample_start
    day_of_year::Integer = 1

    # Loop through the time steps
    while full_step < full_end

        # Get the data for the time step
        full_data = NetCDF.readvar(full_var, start=[1, 1, full_step], count=[size(lons)[1], size(lats)[1], 1])
        sample_data = NetCDF.readvar(sample_var, start=[1, 1, sample_step], count=[size(lons)[1], size(lats)[1], 1])

        println(sample_data)

        quit()

    end


    # Close the output file
    NetCDF.close(output_nc)

    # Close the input files
    NetCDF.close(full_nc)
    NetCDF.close(sample_nc)    
end

#=
    Calculate the true number of time steps.
    If there are no leap years, the is the same
    as the supplied number of time steps.
    Otherwise the leap year days are removed
=#
function calc_time_steps(all_time_steps::Integer, start_year::Integer, has_leap_years::Bool)

    time_steps::Integer = all_time_steps

    if has_leap_years

        time_steps = 0
        current_step = 0
        current_year = start_year
        
        while current_step < all_time_steps

            time_steps += 365

            if Dates.isleapyear(Date(current_year))
                current_step += 366
            else
                current_step += 365
            end

            current_year += 1
        end
    end

    return time_steps
end

function build_times(start_year::Integer, days_count::Integer)

    times = zeros(Float64, days_count)

    current_year::Float64 = start_year
    day_of_year::Integer = 0
    
    for day_index in 1:days_count

        day_of_year += 1
        if day_of_year > 365
            current_year += 1
            day_of_year = 1
        end

        times[day_index] = current_year + (day_of_year - 1) / 365
    end

    return times
end

run()