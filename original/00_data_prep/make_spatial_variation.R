# Libraries
library(ncdf4)

# PARAMETERS
PCO2_FILE <- "daily.nc"
PCO2_VAR <- "pco2"
TIME_LENGTH <- 12045

LON_CELL_SIZE <- 2.5
LAT_CELL_SIZE <- 2.5
LON_LENGTH <- 360 / LON_CELL_SIZE
LAT_LENGTH <- 180 / LAT_CELL_SIZE

##############################################################################

getChangeValues <- function(lon, lat, time, source_value) {
    for (lon_loop in 1:LON_LENGTH) {
        for (lat_loop in 1:LAT_LENGTH) {

            if (lon_loop == lon && lat_loop == lat) {
                # NOOP
            } else {
                if (!is.na(pco2[lon_loop, lat_loop, time])) {
                    diff <- abs(source_value - pco2[lon_loop, lat_loop, time])
                    spatial_variations_total[lon, lat, lon_loop, lat_loop] <<- spatial_variations_total[lon, lat, lon_loop, lat_loop] + diff
                    spatial_variations_count[lon, lat, lon_loop, lat_loop] <<- spatial_variations_count[lon, lat, lon_loop, lat_loop] + 1

                }
            }
        }
    }
}

############################################################################


# Construct the output data structure
# This is a 4-D structure. For each cell (lat/lon), there's the mean
# difference between that cell and all the others (also lat/lon).
cat("Initialising data structures\n")
spatial_variations_total = array(0, c(LON_LENGTH, LAT_LENGTH, LON_LENGTH, LAT_LENGTH))
spatial_variations_count = array(0, c(LON_LENGTH, LAT_LENGTH, LON_LENGTH, LAT_LENGTH))

# Open the source pCO2 file#
cat("Reading pco2 data\n")
nc <- nc_open(PCO2_FILE)
pco2 <- ncvar_get(nc, PCO2_VAR)
nc_close(nc)


cat("\n")

for (time_loop in 1:TIME_LENGTH) {
#for (time_loop in 83:83) {

    for (lon_loop in 1:LON_LENGTH) {
    #for (lon_loop in 18:18) {
        for (lat_loop in 1:LAT_LENGTH) {
        #for (lat_loop in 21:21) {
            current_value <- pco2[lon_loop, lat_loop, time_loop]
            if (!is.na(current_value)) {
                cat("\r",time_loop, lon_loop, lat_loop,"       ")

                # Search for other cells with values at the same time.
                # Note that we include values from a 7 timesteps either side too,
                # just to increase the amount of data we can use. It should also
                # temper the most optimistic change values we see.
                for (i in -7:7) {
                    var_index <- time_loop + i
                    if (var_index > 0 && var_index <= TIME_LENGTH) {
                        getChangeValues(lon_loop, lat_loop, var_index, current_value)
                    }
                }
            }
        }
    }
}

cat("\n")

cat("Calculating means\n")
pco2_spatial_variation <- array(NA, c(LON_LENGTH, LAT_LENGTH, LON_LENGTH, LAT_LENGTH))
for (lon_loop1 in 1:LON_LENGTH) {
#for (lon_loop1 in 18:18) {
    cat("Calculating mean for", lon_loop1, "\n")
    for (lat_loop1 in 1:LAT_LENGTH) {
        for (lon_loop2 in 1:LON_LENGTH) {
            for (lat_loop2 in 1:LAT_LENGTH) {
                if (spatial_variations_count[lon_loop1, lat_loop1, lon_loop2, lat_loop2] > 0) {
                    pco2_spatial_variation[lon_loop1, lat_loop1, lon_loop2, lat_loop2] <- spatial_variations_total[lon_loop1, lat_loop1, lon_loop2, lat_loop2] / spatial_variations_count[lon_loop1, lat_loop1, lon_loop2, lat_loop2]
cat(lon_loop1, lat_loop1, lon_loop2, lat_loop2,pco2_spatial_variation[lon_loop1, lat_loop1, lon_loop2, lat_loop2],"\n")
                }
            }
        }
    }
}

cat("Writing output data\n")
save(pco2_spatial_variation, file="pco2_spatial_variation.R")

