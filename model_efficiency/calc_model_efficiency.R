# Calculate the model efficiency on a 1° version
# using Bastani's code (converted from Matlab)

# Currently restricted to the North Atlantic for comparison

library(ncdf4)

# Pre-made 1° biome 9 mask
cat("\rLoading biomes      ")
nc <- nc_open("mean_biomes_fix.nc")
biome10 <- ncvar_get(nc, "MEANBIOMES")
nc_close(nc)

cat("\rPreprocessing biomes      ")

# Clear out everything except Biome 10
biome10[biome10 != 9] <- NA
biome10[is.nan(biome10)] <- NA

# Load SOCAT data for biome 10
cat("\rLoading SOCAT            ")
nc <- nc_open("socat_monthly.nc")
socat <- ncvar_get(nc, "fco2")
socat_times <- ncvar_get(nc, "time")
nc_close(nc)

for (i in 1:length(socat_times)) {
	step_data <- socat[,,i]
	step_data[is.na(biome10)] <- NA
	socat[,,i] <- step_data
}

# Mean of SOCAT observations
cat("\rPreprocessing SOCAT      ")
socat_mean <- apply(socat, c(1,2), mean, na.rm=TRUE)

# Load the interpolation output
cat("\rLoading interpolated fCO₂      ")
nc <- nc_open("fco2_1x1.nc")
interp <- ncvar_get(nc, "fco2")
nc_close(nc)

# Calculate the efficiency
efficiency <- vector(mode="numeric", length=360*180*372)
dim(efficiency) <- c(360, 180, 372)
efficiency[efficiency == 0] <- NA

for (i in 1:360) {
	cat("\rCalculating efficiency", i, "     ")
	for (j in 1:180) {

		socat_interp_diff_sum <- 0
		socat_mean_diff_sum <- 0

		for (k in 1:length(socat_times)) {
			if (!is.na(socat[i, j, k]) && !is.na(interp[i, j, k])) {
				socat_interp_diff <- (socat[i, j, k] - interp[i, j, k])^2
				socat_mean_diff <- (socat[i, j, k] - socat_mean[i, j])^2

				if (socat_mean_diff != 0) {
					efficiency[i, j, k] <- socat_interp_diff / socat_mean_diff
				}
			}
		}
	}
}

cat("\rWriting output            ")
lon_dim <- ncdim_def("lon", "degrees_east", seq(0.5, 359.5, 1))
lat_dim <- ncdim_def("lat", "degrees_north", seq(-89.5, 89.5, 1))
time_dim <- ncdim_def("time", "year", socat_times)

efficiency_var <- ncvar_def("efficiency", "", list(lon_dim, lat_dim, time_dim), -1e35, prec="double")

nc <- nc_create("efficiency.nc", list(efficiency_var))
ncvar_put(nc, efficiency_var, efficiency)
nc_close(nc)
cat("\n")