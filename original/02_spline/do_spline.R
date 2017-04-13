# R
MONTH_START_DAYS <- c(0, 32, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335, 366)

makeCurve <- function(params, series_length, per_year) {

    curve <- vector(mode="numeric", length=series_length)

    for (curve_step in 1:series_length) {
        curve_value <- params[1] + params[2] * (curve_step - 1) * (365 / 12)
        term <- 2
        for (trig_loop in 1:((length(params) - 2) / 2)) {
            term <- term + 1
            curve_value <- curve_value + params[term] * sin(2*pi*trig_loop*((curve_step - 1)/per_year))
            term <- term + 1
            curve_value <- curve_value + params[term] * cos(2*pi*trig_loop*((curve_step - 1)/per_year))
        }

        curve[curve_step] <- curve_value
    }

    return(curve)
}

makeMonthly <- function(series) {

    series_month <- vector(mode="numeric", length=0)
    monthly <- vector(mode="numeric", length=(length(series) / 365 * 12))

    index <- 0
    month_index <- 0
    current_month <- 1

    while (index < length(series)) {
        month_finished <- FALSE
        year_finished <- FALSE

        index <- index + 1

        day_of_year <- index %% 365
        if (day_of_year == 0) {
            day_of_year <- 365
        }

        if (day_of_year == 1 && current_month != 1) {
            month_finished <- TRUE
            year_finished <- TRUE
        } else if (day_of_year >= MONTH_START_DAYS[current_month + 1]) {
            month_finished <- TRUE
        }

        if (month_finished) {
            month_index <- month_index + 1

            series_mean <- NA
            if (length(series_month) > 0) {
                series_mean <- mean(series_month)
            }
            monthly[month_index] <- series_mean

            current_month <- current_month + 1
            if (year_finished) {
                current_month <- 1
                year_finished <- FALSE
            }

            month_finished <- FALSE

            series_month <- vector(mode="numeric", length=0)
        }

        if (!is.na(series[index])) {
            series_month[length(series_month) + 1] <- series[index]
        }

    }

    month_index <- month_index + 1
    series_mean <- NA
    if (length(series_month) > 0) {
        series_mean <- mean(series_month)
    }
    monthly[month_index] <- series_mean

    return (monthly)
}

makeMonthlyUncertainty <- function(series) {

    series_month <- vector(mode="numeric", length=0)
    monthly <- vector(mode="numeric", length=(length(series) / 365 * 12))

    index <- 0
    month_index <- 0
    current_month <- 1

    while (index < length(series)) {
        month_finished <- FALSE
        year_finished <- FALSE

        index <- index + 1

        day_of_year <- index %% 365
        if (day_of_year == 0) {
            day_of_year <- 365
        }

        if (day_of_year == 1 && current_month != 1) {
            month_finished <- TRUE
            year_finished <- TRUE
        } else if (day_of_year >= MONTH_START_DAYS[current_month + 1]) {
            month_finished <- TRUE
        }

        if (month_finished) {
            month_index <- month_index + 1
            series_mean <- NA
            if (length(series_month) > 0) {
                series_mean <- sqrt(sum(series_month ** 2) / length(series_month))
            }
            monthly[month_index] <- series_mean

            current_month <- current_month + 1
            if (year_finished) {
                current_month <- 1
                year_finished <- FALSE
            }

            month_finished <- FALSE

            series_month <- vector(mode="numeric", length=0)
        }

        if (!is.na(series[index])) {
            series_month[length(series_month) + 1] <- series[index]
        }

    }

    month_index <- month_index + 1
    series_mean <- NA
    if (length(series_month) > 0) {
        series_mean <- sqrt(sum(series_month ** 2) / length(series_month))
    }
    monthly[month_index] <- series_mean


    return (monthly)
}

##################################################################

for (lon in 1:144) {
    for (lat in 1:72) {

        cat("\r",lon,lat,"     ")

        spline_file <- paste("output/spline_",lon,"_",lat,".csv",sep="")
        month_curve_file <- paste("output/curve_",lon,"_",lat,".csv",sep="")
        month_meas_file <- paste("output/measurements_",lon,"_",lat,".csv",sep="")
        uncertainty_file <- paste("output/uncertainty_",lon,"_",lat,".csv",sep="")

        file.remove(spline_file)
        file.remove(month_curve_file)
        file.remove(month_meas_file)
        file.remove(uncertainty_file)

        params_file <- paste("curve_fit_results/cell_params_",lon,"_",lat,".csv",sep="")
        if (file.exists(params_file)) {

            measurements <- read.csv(paste("curve_fit_results/cell_series_",lon,"_",lat,".csv",sep=""),header=F)[[2]]
            uncertainties <- read.csv(paste("curve_fit_results/cell_uncertainties_",lon,"_",lat,".csv",sep=""),header=F)[[2]]
            curve_params <- read.csv(paste("curve_fit_results/cell_params_",lon,"_",lat,".csv",sep=""),header=F)[[1]]

            curve <- makeCurve(curve_params, length(measurements) / 365 * 12, 12)
            monthly_measurements <- makeMonthly(measurements)
            monthly_uncertainties <- makeMonthlyUncertainty(uncertainties)

            combined <- curve
            for (i in 1:length(combined)) {
                if (!is.na(monthly_measurements[i])) {
                    combined[i] <- monthly_measurements[i]
                }
            }

            spline <- smooth.spline(combined,all.knots=TRUE,spar=0.3)$y

            sink(month_meas_file, append=T)
            for (i in 1:length(spline)) {
                cat(i,",",monthly_measurements[i],"\n",sep="")
            }
            sink()

            sink(uncertainty_file, append=T)
            for (i in 1:length(spline)) {
                cat(i,",",monthly_uncertainties[i],"\n",sep="")
            }
            sink()

            sink(month_curve_file, append=T)
            for (i in 1:length(spline)) {
                cat(i,",",curve[i],"\n",sep="")
            }
            sink()

            sink(spline_file, append=T)
            for (i in 1:length(spline)) {
                cat(i,",",spline[i],"\n",sep="")
            }
            sink()
        }
    }
}

cat("\n")

