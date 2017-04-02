# Libraries
library(ncdf4)

# PARAMETERS
IN_DIR <- "cell_series_daily"
LON_CELL_SIZE <- 2.5
LAT_CELL_SIZE <- 2.5
DATES_FILE <- "days.csv"
OUT_FILE <- "daily.nc"

# Build lons
lons <- vector(mode="numeric", length=(360 / LON_CELL_SIZE))
for (i in 1:length(lons)) {
    lons[i] <- i * LON_CELL_SIZE - (LON_CELL_SIZE / 2)
}

lats <- vector(mode="numeric", length=(180 / LAT_CELL_SIZE))
for (i in 1:length(lats)) {
    lats[i] <- i * LAT_CELL_SIZE - (90 + LON_CELL_SIZE / 2)
}

dates <- read.csv(DATES_FILE, header=FALSE)[[1]]

pco2 <- vector(mode="numeric", length=(length(lons) * length(lats) * length(dates)))
dim(pco2) <- c(length(lons), length(lats), length(dates))

for (lon in 1:(360 / LON_CELL_SIZE)) {
    cat("\rReading",lon," ")
    for (lat in 1:(180 / LAT_CELL_SIZE)) {

        in_file <- paste(IN_DIR, "/", "cell_series_", lon, "_", lat, ".csv", sep="")
        series_data <- read.csv(in_file, header=F)[[2]]

        pco2[lon, lat, ] <- series_data

    }
}

cat("\n")

cat("Writing...\n")
lon_dim <- ncdim_def("LON", "degrees_east", lons)
lat_dim <- ncdim_def("LAT", "degrees_north", lats)
time_dim <- ncdim_def("TIME", "year", dates)

pco2_var <- ncvar_def("pco2", "uatm", list(lon_dim, lat_dim, time_dim), -999, prec="double")

nc <- nc_create(OUT_FILE, list(pco2_var))

ncvar_put(nc, "pco2", pco2)

nc_close(nc)
