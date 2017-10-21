# Calculate the R-squared error for each grid cell

library(ncdf4)

# Load SOCAT data
cat("\rLoading SOCAT            ")
nc <- nc_open("socat_monthly.nc")
socat <- ncvar_get(nc, "fco2")
nc_close(nc)

# Load the interpolation output
cat("\rLoading interpolated fCO₂      ")
nc <- nc_open("fco2.nc")
lons <- ncvar_get(nc, "lon")
lats <- ncvar_get(nc, "lat")
interp <- ncvar_get(nc, "fco2")
nc_close(nc)

rsquared <- vector(mode="numeric", length=(length(lons) * length(lats)))
dim(rsquared) <- c(length(lons), length(lats))
rsquared[rsquared == 0] <- NA
socat_count <- rsquared

for (lon in 1:length(lons)) {
	cat("\rCalculating r²", lon, "      ")

	for (lat in 1:length(lats)) {

		socat_series <- socat[lon, lat, ]
		interp_series <- interp[lon, lat, ]

		if (sum(!is.na(socat_series)) > 2) {
			correlation <- cor(socat_series, interp_series, use="pairwise.complete.obs")
			
			if (!is.na(correlation)) {
				sign <- 1
				if (correlation < 0) {
					sign <- -1
				}

				rsquared[lon, lat] <- (correlation ^ 2) * sign
				socat_count[lon, lat] <- sum(!is.na(socat_series))
			}
		}
	}
}

cat("\rWriting output            ")
lon_dim <- ncdim_def("lon", "degrees_east", lons)
lat_dim <- ncdim_def("lat", "degrees_north", lats)

r2_var <- ncvar_def("r2", "", list(lon_dim, lat_dim), -1e35, prec="double")
count_var <- ncvar_def("socat_count", "", list(lon_dim, lat_dim), -9999, prec="integer")

nc <- nc_create("rsquared.nc", list(r2_var, count_var))
ncvar_put(nc, r2_var, rsquared)
ncvar_put(nc, count_var, socat_count)
nc_close(nc)
cat("\n")