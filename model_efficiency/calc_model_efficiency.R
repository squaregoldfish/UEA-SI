# Calculate the model efficiency on a 1° version
# using Bastani's code (converted from Matlab)
library(ncdf4)

# Pre-made 1° biome 9 mask
nc <- nc_open("Time_Varying_Biomes.nc")
biome10 <- ncvar_get(nc, "MeanBiomes")
nc_close(nc)

# Rotate the biomes map so it's lon*lat
biome10 <- aperm(biome10, c(2,1))


biome10[biome10 != 9] <- NA
biome10[is.nan(biome10)] <- NA

nc <- nc_open("socat_monthly.nc")
socat <- ncvar_get(nc, "fco2")
socat_times <- ncvar_get(nc, "time")
nc_close(nc)

for (i in 1:length(socat_times)) {
	step_data <- socat[,,i]
	step_data[is.na(biome10)] <- NA
	socat[,,i] <- step_data
}
