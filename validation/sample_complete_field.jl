# Take a complete pCO2 field and sample it to match SOCAT

if size(ARGS) != (2,)
    println("Usage: julia sample_complete_field.jl [compelte_field] [sample_map]")
    quit()
end

