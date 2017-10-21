# Calculate the model efficiency on a 1° version
# using Bastani's code (converted from Matlab)

library(ncdf4)

# Load SOCAT data for biome 10
cat("\rLoading SOCAT            ")
nc <- nc_open("socat_monthly.nc")
socat <- ncvar_get(nc, "fco2")
socat_times <- ncvar_get(nc, "time")
nc_close(nc)

# Mean of SOCAT observations
cat("\rPreprocessing SOCAT      ")
socat_mean <- apply(socat, c(1,2), mean, na.rm=TRUE)

# Load the interpolation output
cat("\rLoading interpolated fCO₂      ")
nc <- nc_open("fco2.nc")
interp <- ncvar_get(nc, "fco2")
nc_close(nc)

# Calculate the efficiency
efficiency <- vector(mode="numeric", length=144*72*372)
dim(efficiency) <- c(144, 72, 372)
efficiency[efficiency == 0] <- NA

for (i in 1:144) {
	cat("\rCalculating efficiency", i, "     ")
	for (j in 1:72) {

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
lon_dim <- ncdim_def("lon", "degrees_east", seq(1.25, 358.75, 2.5))
lat_dim <- ncdim_def("lat", "degrees_north", seq(-88.75, 88.75, 2.5))
time_dim <- ncdim_def("time", "year", socat_times)

efficiency_var <- ncvar_def("efficiency", "", list(lon_dim, lat_dim, time_dim), -1e35, prec="double")

nc <- nc_create("efficiency.nc", list(efficiency_var))
ncvar_put(nc, efficiency_var, efficiency)
nc_close(nc)
cat("\n")