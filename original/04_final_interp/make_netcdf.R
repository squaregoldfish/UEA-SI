# R!
library(ncdf4)

MONTH_STARTS <- c(0.000, 0.085, 0.162, 0.247, 0.329, 0.414, 0.496, 0.581, 0.666, 0.748, 0.833, 0.915)
START_YEAR <- 1985

LON_SIZE <- 144
LAT_SIZE <- 72
VALID_LAT_SIZE <- 72
TIME_SIZE <- 324


lons <- vector(mode="numeric", length=LON_SIZE)
for (i in 1:length(lons)) {
    lons[i] <- (i * 2.5) - 1.25
}

lats <- vector(mode="numeric", length=LAT_SIZE)
for (i in 1:length(lats)) {
    lats[i] <- (i * 2.5) - 91.25
}

times <- vector(mode="numeric", length=TIME_SIZE)
current_month <- 0
current_year <- START_YEAR
for (i in 1:length(times)) {
    current_month <- current_month + 1
    if (current_month > 12) {
        current_month <- 1
        current_year <- current_year + 1
    }

    times[i] <- current_year + MONTH_STARTS[current_month]
}

pco2 <- vector(mode="numeric", length=LON_SIZE * LAT_SIZE * TIME_SIZE)
pco2[pco2 == 0] <- NA
dim(pco2) <- c(LON_SIZE, LAT_SIZE, TIME_SIZE)

uncertainties <- pco2


for (lon_loop in 1:LON_SIZE) {
    cat("\r",lon_loop,"   ")
    for (lat_loop in 1:VALID_LAT_SIZE) {

        pco2_file <- paste("output/spline_",lon_loop,"_",lat_loop,".csv",sep="")
        if (file.exists(pco2_file)) {
            cell_pco2 <- read.csv(pco2_file, header=F)[[2]]

            for (time_loop in 1:TIME_SIZE) {
                pco2[lon_loop, lat_loop, time_loop] <- cell_pco2[time_loop]
            }
        }

        uncertainty_file <- paste("output/uncertainty_",lon_loop,"_",lat_loop,".csv",sep="")
        if (file.exists(uncertainty_file)) {
            cell_uncertainty <- read.csv(uncertainty_file, header=F)[[2]]

            for (time_loop in 1:TIME_SIZE) {
                uncertainties[lon_loop, lat_loop, time_loop] <- cell_uncertainty[time_loop]
            }
        }


    }
}

cat("\n")

lon_dim <- ncdim_def("lon", "degrees_east", lons)
lat_dim <- ncdim_def("lat", "degrees_north", lats)
time_dim <- ncdim_def("time", "year", times, unlim=TRUE)

pco2_var <- ncvar_def("fco2", "uatm", list(lon_dim, lat_dim, time_dim), -1e35, prec="double")
uncertainty_var <- ncvar_def("uncertainty", "uatm", list(lon_dim, lat_dim, time_dim), -1e35, prec="double")

nc <- nc_create("fco2.nc", list(pco2_var, uncertainty_var))
ncvar_put(nc, pco2_var, pco2)
ncvar_put(nc, uncertainty_var, uncertainties)
#ncatt_put(nc, "time", "calendar", "noleap")
ncatt_put(nc, 0, "Title", "Statistically gap-filled fCO2 data south of 70N, 1985-2011")
ncatt_put(nc, 0, "Version", "1.0 - 2015-09-10")
ncatt_put(nc, 0, "Info", "Derived from SOCAT v2")
ncatt_put(nc, 0, "Reference", "A statistical gap-filling method to interpolate global monthly surface ocean carbon dioxide data, doi to come")
ncatt_put(nc, 0, "Contact", "Steve Jones <s.d.jones@exeter.ac.uk>")
nc_close(nc)



