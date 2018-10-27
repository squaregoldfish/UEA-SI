# Add uncertainty time series files alongside the actual time series
SERIES_DIR <- "cell_series_daily"
SERIES_LENGTH <- 12045

for (lon in 1:144) {
#for (lon in 72:72) {
    for (lat in 1:72) {
    #for (lat in 36:36) {
        cat("\r",lon,lat,"   ")

        series_file <- paste(SERIES_DIR, "/cell_series_", lon, "_", lat, ".csv", sep="")
        series <- as.double(unlist(read.csv(series_file,header=F)[2]))

        uncertainties <- array(NA, c(SERIES_LENGTH))
        uncertainties[!is.na(series)] <- 2.5
   
        weights <- array(NA, c(SERIES_LENGTH))
        weights[!is.na(series)] <- 1.0

        uncertainty_file <- paste(SERIES_DIR, "/cell_uncertainties_", lon, "_", lat, ".csv", sep="")
        sink(uncertainty_file)
        for (i in 1:SERIES_LENGTH) {
            cat(i, ",", uncertainties[i], "\n", sep="")
        }
        sink()

        weights_file <- paste(SERIES_DIR, "/cell_weights_", lon, "_", lat, ".csv", sep="")
        sink(weights_file)
        for (i in 1:SERIES_LENGTH) {
            cat(i, ",", weights[i], "\n", sep="")
        }
        sink()
    }
}

cat("\n")

