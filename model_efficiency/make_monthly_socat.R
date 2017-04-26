library(ncdf4)

MONTH_END_DAYS <- c(31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334, 365)

nc <- nc_open("daily.nc")
lons <- ncvar_get(nc, "LON")
lats <- ncvar_get(nc, "LAT")
times <- ncvar_get(nc, "TIME")
daily <- ncvar_get(nc, "pco2")
nc_close(nc)

month_count <- length(times) / 365 * 12

monthly <- vector(mode="numeric", length=(length(lons) * length(lats) * month_count))
dim(monthly) <- c(length(lons), length(lats), month_count)
monthly[monthly == 0] <- NA

for (lon_loop in 1:length(lons)) {
    cat("\r",lon_loop, "            ")
    for (lat_loop in 1:length(lats)) {

        current_month <- 1
        month_of_year <- 1
        day_of_year <- 0
        month_total <- 0
        days_used <- 0

        for (time_loop in 1:length(times)) {
            day_of_year <- day_of_year + 1

            if (!is.na(daily[lon_loop, lat_loop, time_loop])) {
                month_total <- month_total + daily[lon_loop, lat_loop, time_loop]
                days_used <- days_used + 1
            }

            if (day_of_year == MONTH_END_DAYS[month_of_year]) {
                if (days_used > 0) {
                    cat("\r", lon_loop, lat_loop, time_loop, "                    ")
                    monthly[lon_loop, lat_loop, current_month] <- month_total / days_used
                }

                current_month <- current_month + 1
                month_of_year <- month_of_year + 1
                if (month_of_year == 13) {
                    month_of_year <- 1
                    day_of_year <- 0
                }
                month_total <- 0
                days_used <- 0
            }
        }
    }
}

cat("\rWriting...                    ")

lon_dim <- ncdim_def("lon", "degrees_east", lons)
lat_dim <- ncdim_def("lat", "degrees_north", lats)
time_dim <- ncdim_def("time", "year", seq(1, month_count), unlim=TRUE)

fco2_var <- ncvar_def("fco2", "uatm", list(lon_dim, lat_dim, time_dim), -1e35, prec="double")

nc <- nc_create("socat_monthly.nc", list(fco2_var))
ncvar_put(nc, fco2_var, monthly)
nc_close(nc)

cat("\n")
