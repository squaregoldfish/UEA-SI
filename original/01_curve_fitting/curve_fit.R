# Curve fitting functions for do_interpolation.R

# Attempt to perform a curve fit. Includes temporal interpolation of data.
# Temporal interpolation only adds extra values to give more data to the curve
# fitting algorithm; the interpolated values are never kept.
attemptFit <- function(measurements, weights, uncertainties) {
    fitted_params <- NULL
    used_measurements <- measurements

    series <<- measurements # Need to do this for the fitCurve function call - I don't know why!
    fit_params <- fitCurve(measurements, weights)

    # If the fit was successful, return the result
    if (!is.null(fit_params)) {
        fitted_params <- fit_params
    } else {

        # Apply temporal interpolation
        cat("Performing temporal interpolation\n")
        interpolated_values <- calcTemporalInterpolationValues(measurements, weights, uncertainties)
        interpolated_measurements <- mergeValues(current_measurements, interpolated_values[, 1])
        interpolated_weights <- mergeValues(current_weights, interpolated_values[, 2])
        interpolated_uncertainties <- mergeValues(current_uncertainties, interpolated_values[, 3])

        series <<- interpolated_measurements # Need to do this for the fitCurve function call - I don't know why!
        fit_params <- fitCurve(interpolated_measurements, interpolated_weights)

        if (!is.null(fit_params)) {
            used_measurements <- interpolated_measurements
        }
    }

    if (!is.null(fit_params)) {
        # We had a successful fit. Store the measurents (excluding
        # temporally interpolated ones) for later.
        stored_measurements <<- measurements
        stored_weights <<- weights
        stored_uncertainties <<- uncertainties
        stored_fit_params <<- fit_params
    }

    fitted_curve <- NULL
    if (!is.null(fitted_params)) {
        fitted_curve <- makeCurve(fitted_params, length(measurements))
    }
    
    return (fitted_curve)
}

removeOutliers <- function(series) {

    filtered_series <- vector(mode="numeric",length=length(series))
    filtered_series[filtered_series == 0] <- NA

    if (sum(!is.na(series)) > 0) {
        mean_value <- mean(series,na.rm=T)
        stdev <- sd(series,na.rm=T)

        for (i in 1:length(series)) {
            if (!is.na(series[i])) {
                stdevs_from_mean <- floor((abs(series[i] - mean_value)) / stdev)
                if (!is.na(stdevs_from_mean) && stdevs_from_mean < 3) {
                    filtered_series[i] <- series[i]
                }

            }
        }
    }

    return (filtered_series)
}


# Fit a curve to a time series
fitCurve <- function(series, weights) {
    fitted_params <- NULL
    result <- NULL

    # Remove outliers
    fit_series <<- removeOutliers(series)

    prefit_ok <- doPrefitCheck(fit_series)
    if (!prefit_ok) {
        cat("PREFIT CHECKS FAILED\n")
    } else {
        # Set up the initial curve parameters
        col_names <- c("a", "b", "c", "d", "e", "f", "g", "h", "i", "j")
        fit_params <- c(a=340, b=30, c=-1, d=39, e=4, f=17, g=-8.3, h=-8.1, i=0.3, j=-1)

        fit_success <- FALSE
        harmonic_count <- 4

        while (!fit_success && harmonic_count > 0) {
            cat("Fitting",harmonic_count," harmonics\n")

            # Construct the formula
            x <<- 1:length(fit_series) # For some reason, x has to live in the global space so nls can find it.
            formula <- "a + b*(x)"
            term_count <- 2

            for (trig_loop in 1:harmonic_count) {
                formula <- sprintf("%s + %s*sin(2*%s*%d*(x/365)) + %s*cos(2*%s*%d*(x/365))",formula,col_names[term_count + 1],pi,trig_loop,col_names[term_count + 2],pi,trig_loop)
                term_count <- term_count + 2
            }

            # Perform the curve fit
            start_values <- .subset(fit_params, 1:term_count)

            model <- NULL
            tryCatch({
                model <- nls(paste("fit_series ~ ",formula), weights=weights, start=start_values, control=nls.control(minFactor=1/(8192*8192)))
            }, error=function(err){print(err)})

            if (!is.null(model)) {
                # The fit was successful. Now check it to see if it's a good fit

                # Check that the linear slope is OK
                slope <- model$m$getPar()[2]
                if (slope > MAX_LINEAR_TREND || slope < MIN_LINEAR_TREND) {
                    cat("Curve slope is too steep (",slope,"): should be ",MAX_LINEAR_TREND," or less\n",sep="")
                    fit_success <- FALSE
                } else {
                    fit_success <- TRUE
                }

                # Now we need the complete curve
                fitted_params <- model$m$getPar()
                fit_curve <- makeCurve(fitted_params, length(fit_series))

                if (fit_success && harmonic_count != 1) {
                    fit_success <- checkCurvePeaks(fitted_params)
                }

                if (fit_success) {
                    fit_success <- checkCurveFit(fit_series, fit_curve, term_count)
                }
            }

            if (fit_success) {
                result <- fitted_params
            } else {
                harmonic_count <- harmonic_count - 1
            }
        }
    }

    return (result)
}

# Perform pre-curve-fitting checks on the time series
doPrefitCheck <- function(series) {

    ok <- TRUE

    # Check the standard deviation of the measurements
    stdev <- sd(series, na.rm=T)
    if (!is.na(stdev) && stdev > MAX_STDEV) {
        cat("Standard deviation is ", stdev, ", should be <= ", MAX_STDEV, "\n", sep="")
        ok <- FALSE
    }

    # Check the measurement coverage
    measurements <- which(!is.na(series))
    if (length(measurements) < 2) {
        cat("Measurements must span at least ", MIN_TIME_SPAN, " days\n", sep="")
        ok <- FALSE
    } else {
        measurement_span <- measurements[length(measurements)] - measurements[1]
        if (measurement_span < MIN_TIME_SPAN) {
            cat("Measurements only span ", measurement_span, " days, should be ", MIN_TIME_SPAN, "\n", sep="")
            ok <- FALSE
        }

        populated_months <- array(FALSE, c(12))
        for (i in 1:length(series)) {
            if (!is.na(series[i])) {
                day <- i %% 365
                if (day == 0) {
                    day <- 365
                }
                month <- tail(which(MONTH_STARTS <= (day / 365)), 1)
                populated_months[month] <- TRUE
            }
        }

        if (sum(populated_months == TRUE) < MIN_POPULATED_MONTHS) {
            cat("No. of populated months = ", sum(populated_months == TRUE), ", should be ", MIN_POPULATED_MONTHS, "\n", sep="")
            ok <- FALSE
        }
    }

    return (ok)
}

checkCurveFit <- function(series, fit_curve, degrees) {
    curve_ok <- TRUE

    # Compare the curve range with the series range
    series_max <- max(series, na.rm=TRUE)
    series_min <- min(series, na.rm=TRUE)
    curve_max <- max(fit_curve, na.rm=TRUE)
    curve_min <- min(fit_curve, na.rm=TRUE)
    series_range <- series_max - series_min
    curve_range <- curve_max - curve_min
    range_ratio <- curve_range / series_range

    cat("Curve fit to series ratio = ",range_ratio)

    if (range_ratio < MIN_CURVE_RATIO || range_ratio > MAX_CURVE_RATIO) {
        curve_ok <- FALSE
        cat(", Ratio FAILED should be",MIN_CURVE_RATIO,"to",MAX_CURVE_RATIO,"\n")
    } else {
        cat("\n")
    }

    max_diff <- abs(series_max - curve_max)
    min_diff <- abs(series_min - curve_min)

    cat("Limit differences =",max_diff,",",min_diff)
    if (max_diff > MAX_LIMIT_DIFFERENCE || min_diff > MAX_LIMIT_DIFFERENCE) {
        curve_ok <- FALSE
        cat(", Limit difference FAILED should be",MAX_LIMIT_DIFFERENCE,"or less\n")
    } else {
        cat("\n")
    }

    return(curve_ok)
}

checkCurvePeaks <- function(curve_params) {

    peaks_ok <- TRUE

    curve <- makeSeasonalCurve(curve_params)


    # Loop through the curve, finding maxima and minima
    maxima <- vector(mode="numeric", length=0)
    maxima_pos <- vector(mode="numeric", length=0)
    minima <- vector(mode="numeric", length=0)
    minima_pos <- vector(mode="numeric", length=0)

    last_value <- curve[1]
    slope_direction <- 0 # Starting value; 1 = up, -1 = down
    start_direction <- 0
    end_direction <- 0
    curve_amplitude <- max(curve) - min(curve)


    for (i in 2:365) {
        if (curve[i] >= last_value) {
            if (slope_direction == -1) {
                # Found a trough!
                minima[length(minima) + 1] <- last_value
                minima_pos[length(minima_pos) + 1] <- i - 1

            }
            slope_direction <- 1
            
            # Record the start direction
            if (start_direction == 0) {
                start_direction <- slope_direction
            }

        } else {
            if (slope_direction > 0) {
                # We found a peak!
                maxima[length(maxima) + 1] <- last_value
                maxima_pos[length(maxima_pos) + 1] <- i - 1
            }
            slope_direction <- -1

            # Record the start direction
            if (start_direction == 0) {
                start_direction <- slope_direction
            }
        }

        last_value <- curve[i]
    }

    # Record the end direction
    end_direction <- slope_direction

    # If the start and end directions are different, there's a peak/trough at zero.
    if (start_direction == 1 && end_direction == -1) {
        minima <- c(curve[1],minima)
        minima_pos <- c(1, minima_pos)
    }

    if (start_direction == -1 && end_direction == 1) {
        maxima[length(maxima) + 1] <- curve[365]
        maxima_pos[length(maxima_pos) + 1] <- 365
    }

    # Check the number of peaks
    if (length(maxima) > MAX_PEAKS_PER_YEAR) {
        cat("Peak check FAILED - too many peaks in seasonal cycle\n")
        peaks_ok <- FALSE
    } else {
        # Calculate the peak sizes
        peak_sizes <- vector(mode="numeric", length=length(maxima))

        for (i in 1:length(maxima)) {
            max_value <- maxima[i]
            max_pos <- maxima_pos[i]

            # For each peak, we find the position of the preceding minimum
            # We loop round to the end of the year if necessary
            min_entry <- 0
            min_pos_candidates <- which(minima_pos < max_pos)
            if (length(min_pos_candidates) == 0) {
                min_entry <- length(minima)
            } else {
                min_entry <- min_pos_candidates[length(min_pos_candidates)]
            }

            min_value <- minima[min_entry]
            peak_sizes[i] <- abs(max_value - min_value)
        }

        peak_size_limit <- curve_amplitude * MAX_PEAK_RATIO
        large_peak_count <- length(which(peak_sizes > peak_size_limit))
        if (large_peak_count > 1) {
            cat("Peak check FAILED - secondary peak(s) are too large\n")
            peaks_ok <- FALSE
        }
    }

    
    return(peaks_ok)
}

# Make a curve from the supplied parameters
makeCurve <- function(params, series_length) {
    curve <- vector(mode="numeric", length=series_length)

    for (curve_step in 1:series_length) {
        curve_value <- params[1] + params[2] * (curve_step - 1)
        term <- 2
        for (trig_loop in 1:((length(params) - 2) / 2)) {
            term <- term + 1
            curve_value <- curve_value + params[term] * sin(2*pi*trig_loop*((curve_step - 1)/365))
            term <- term + 1
            curve_value <- curve_value + params[term] * cos(2*pi*trig_loop*((curve_step - 1)/365))
        }

        curve[curve_step] <- curve_value
    }

    return(curve)
}

# Make a seasonal cycle curve from the supplied parameters
makeSeasonalCurve <- function(params) {
    curve <- vector(mode="numeric", length=365)

    for (curve_step in 1:365) {
        curve_value <- 0
        term <- 2
        for (trig_loop in 1:((length(params) - 2) / 2)) {
            term <- term + 1
            curve_value <- curve_value + params[term] * sin(2*pi*trig_loop*((curve_step - 1)/365))
            term <- term + 1
            curve_value <- curve_value + params[term] * cos(2*pi*trig_loop*((curve_step - 1)/365))
        }

        curve[curve_step] <- curve_value
    }

    return(curve)
}

