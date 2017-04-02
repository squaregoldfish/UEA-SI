# R!

getPrevValue <- function(series, start) {

    prev_value <- NA
    pos <- start

    while (is.na(prev_value)) {
        pos <- pos - 1
        if (pos == 0) {
            pos <- length(series)
        }

        if (!is.na(series[pos])) {
            prev_value <- series[pos]
        }
    }

    return (prev_value)
}

getNextValue <- function(series, start) {
    next_value <- NA
    pos <- start

    while (is.na(next_value)) {
        pos <- pos + 1
        if (pos > length(series)) {
            pos <- 1
        }

        if (!is.na(series[pos])) {
            next_value <- series[pos]
        }
    }

    return (next_value)
}


for (lon in 1:144) {
#for (lon in 22:22) {
    cat("\r",lon,"    ")
    for (lat in 1:72) {
    #for (lat in 28:28) {


        # Load data
        measurements_file <- paste("spline_output/measurements_",lon,"_",lat,".csv",sep="")
        uncertainties_file <- paste("spline_output/uncertainty_",lon,"_",lat,".csv",sep="")
        curve_file <- paste("spline_output/curve_",lon,"_",lat,".csv",sep="")
        spline_file <- paste("spline_output/spline_",lon,"_",lat,".csv",sep="")

        if (file.exists(curve_file)) {
            measurements <- read.csv(measurements_file,header=F)[[2]]
            uncertainties <- read.csv(uncertainties_file,header=F)[[2]]
            curve<- read.csv(curve_file,header=F)[[2]]
            spline<- read.csv(spline_file,header=F)[[2]]

            # Calculate the curve fit uncertainties.
            month_uncertainties <- vector(mode="numeric", length=12)
            month_uncertainties[month_uncertainties == 0] <- NA

            # Collect the measurements for each calendar month in turn.
            for (month in 1:12) {
                sequence <- seq(month,length(measurements),by=12)

                collected_uncertainties <- vector(mode="numeric", length=0)

                for (i in sequence) {
                    if (!is.na(measurements[i])) {
                        # The curve uncertainty is the measurement uncertainty plus the measurement's distance from the curve
                        measurement_uncertainty <- uncertainties[i]
                        curve_anomaly <- abs(measurements[i] - curve[i])
                        collected_uncertainties[length(collected_uncertainties) + 1] <- curve_anomaly
                    }
                }

                # The uncertainty for the curve is the mean of all the uncertainties calculated above
                if (length(collected_uncertainties) > 0) {
                    month_uncertainties[month] <- sqrt(sum(collected_uncertainties ** 2) / length(collected_uncertainties))
                }
            }


            # For some months, there are no measurements from which we can caluculate the curve uncertainty.
            # Find the uncertainties either side of the gap and add them together.
            completed_month_uncertainties <- month_uncertainties
            for (month in 1:12) {
                if (is.na(month_uncertainties[month])) {
                    prev_value = getPrevValue(month_uncertainties, month)
                    next_value = getNextValue(month_uncertainties, month)

                    completed_month_uncertainties[month] <- sqrt(prev_value ** 2 + next_value ** 2)
                }
            }

            # Now we have uncertainties for the full seasonal cycle.
            # We can apply all known uncertainties to the spline fit
            # to create a time series of uncertainties.
            all_uncertainties <- vector(mode="numeric", length=length(spline))
            all_uncertainties[all_uncertainties == 0] <- NA

            # Build a time series of the curve's uncertainty
            for (i in 1:length(spline)) {
                month <- i %% 12
                if (month == 0) {
                    month <- 12
                }

                all_uncertainties[i] <- completed_month_uncertainties[month]
            }

            # Now insert uncertainties from any measurements
            for (i in 1:length(spline)) {
                month <- i %% 12
                if (month == 0) {
                    month <- 12
                }

                if (!is.na(measurements[i])) {

                    measurement_value <- measurements[i]
                    measurement_uncertainty <- uncertainties[i]
                    spline_value <- spline[i]
                    spline_uncertainty <- sqrt(measurement_uncertainty ** 2 + abs(spline_value - measurement_value) ** 2)

                    all_uncertainties[i] <- spline_uncertainty
                }
            }

            # Copy the spline file to the output
            file.copy(spline_file, paste("output/spline_",lon,"_",lat,".csv",sep=""))
            file.copy(measurements_file, paste("output/measurements_",lon,"_",lat,".csv",sep=""))
            
            # Write the uncertainties
            sink(paste("output/uncertainty_",lon,"_",lat,".csv",sep=""))
            for (i in 1:length(all_uncertainties)) {
                cat(i,",",all_uncertainties[i],"\n",sep="")
            }
            sink()
        }
    }
}

cat("\n")
