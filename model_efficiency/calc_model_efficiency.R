library(ncdf4)

nc <- nc_open("socat_monthly.nc")
lons <- ncvar_get(nc, "lon")
lats <- ncvar_get(nc, "lat")
socat <- ncvar_get(nc, "fco2")
nc_close(nc)

nc <- nc_open("fco2.nc")
interp <- ncvar_get(nc, "fco2")
nc_close(nc)

efficiency <- vector(mode="numeric", length=(length(lons) * length(lats)))
dim(efficiency) <- c(length(lons), length(lats))
efficiency[efficiency == 0] <- NA

for (lon_loop in 1:length(lons)) {
	for (lat_loop in 1:length(lats)) {

		socat_series <- socat[lon_loop, lat_loop, ]
		if (sum(!is.na(socat_series)) > 0) {

			interp_series <- interp[lon_loop, lat_loop, ]

			mean_obs <- mean(socat_series, na.rm=TRUE)
			interp_diff_sq_sum <- 0
			obs_diff_sq_sum <- 0

			for (time_loop in 1:length(socat_series)) {
				if (!is.na(socat_series[time_loop])) {

					interp_diff_sq_sum <- interp_diff_sq_sum + (socat_series[time_loop] - interp_series[time_loop]) ** 2
					obs_diff_sq_sum <- obs_diff_sq_sum + (socat_series[time_loop] - mean_obs) ** 2

				}
			}

			if (obs_diff_sq_sum == 0) {
				efficiency[lon_loop, lat_loop] <- 1
			} else {
				efficiency[lon_loop, lat_loop] <- interp_diff_sq_sum / obs_diff_sq_sum
			}
		}
	}
}

lon_dim <- ncdim_def("lon", "degrees_east", lons)
lat_dim <- ncdim_def("lat", "degrees_north", lats)

efficiency_var <- ncvar_def("efficiency", "", list(lon_dim, lat_dim), -1e35, prec="double")

nc <- nc_create("efficiency.nc", list(efficiency_var))
ncvar_put(nc, efficiency_var, efficiency)
nc_close(nc)
