# Interpolation
library(ncdf4)

# External source files
source("curve_fit.R")
source("spatial_interpolation.R")

# PARAMETERS
SERIES_FILE_START <- "cell_series_"
UNCERTAINTIES_FILE_START <- "cell_uncertainties_"
WEIGHTS_FILE_START <- "cell_weights_"
TEMPORAL_INTERPOLATION_LIMIT <- 7
SPATIAL_INTERPOLATION_WEIGHT_LIMIT <- 0.367879
START_DATE <- 1985

MAKE_PLOTS <- FALSE
DIAGNOSTIC <- TRUE

# Calendar points
MONTH_STARTS <- c(0.0000, 0.0849, 0.1616, 0.2466, 0.3288, 0.4134, 0.4959, 0.5808, 0.6678, 0.7479, 0.8329, 0.9151)

# Curve fitting restrictions
MIN_TIME_SPAN <- 1825
MAX_STDEV <- 75
MIN_POPULATED_MONTHS <- 8
MAX_LINEAR_TREND <- 4.5
MIN_LINEAR_TREND <- -2.5
MAX_PEAKS_PER_YEAR <- 2
MAX_PEAK_RATIO <- 0.33
MIN_CURVE_RATIO <- 0.5
MAX_CURVE_RATIO <- 1.5
MAX_LIMIT_DIFFERENCE <- 75

# Other bits
LON_CELL_SIZE <- 2.5
LAT_CELL_SIZE <- 2.5
LON_SIZE <- 360 / LON_CELL_SIZE
LAT_SIZE <- 180 / LAT_CELL_SIZE
SPATIAL_ACF_LAG_STEP <- 25
EARTH_RADIUS <- 6367.5
SPATIAL_INTERP_WEIGHT_LIMIT <- 0.367879

##########################################################

# Calculate the values for temporal interpolation.
# The returned value is a list of <time step, interpolated value, weight, uncertainty>
# to be added to an existing series.
#
# We only ever interpolate one time step
calcTemporalInterpolationValues <- function(series, series_weights, series_uncertainties) {
    interpolated_points <- array(NA, c(length(series), 3))

    for (i in 1:length(series)) {
    #for (i in 6451:6451) {
        if (is.na(series[i])) {
            interp_start <- i - TEMPORAL_INTERPOLATION_LIMIT
            if (interp_start < 1) {
                interp_start <- 1
            }

            interp_end <- i + TEMPORAL_INTERPOLATION_LIMIT
            if (interp_end > length(series)) {
                interp_end <- length(series)
            }

            interp_total <- 0
            interp_weight_sum <- 0
            interp_count <- 0

            for (interp_loop in interp_start:interp_end) {
                if (!is.na(series[interp_loop])) {
                    lag <- abs(i - interp_loop)

                    interp_total <- interp_total + series[interp_loop] * temporal_acf[lag]
                    interp_weight_sum <- interp_weight_sum + temporal_acf[lag]
                    interp_count <- interp_count + 1
                }
            }

            if (interp_count > 0) {
                new_value <- interp_total / interp_weight_sum
                new_weight <- interp_weight_sum / interp_count

                # There's no need for uncertainties, as these values get discarded later.
                new_uncertainty <- 2.5

                interpolated_points[i,] <- c(new_value, new_weight, new_uncertainty)
            }
        }
    }

    return (interpolated_points)
}

# Merge two series of values together. Values from the original always override new values
mergeValues <- function(original, new) {
    result <- vector(mode="numeric", length=length(original))
    result[result == 0] <- NA

    for (i in 1:length(original)) {
        if (!is.na(original[i])) {
            result[i] <- original[i]
        } else if (!is.na(new[i])) {
            result[i] <- new[i]
        }
    }

    return (result)
}

loadCellMeasurements <- function(indir, lon, lat) {
    series_file <- paste(indir, "/", SERIES_FILE_START, lon, "_", lat, ".csv", sep="")
    loaded_series <- as.double(unlist(read.csv(series_file, header=FALSE)[2]))
}

loadCellWeights <- function(indir, lon, lat) {
    series_file <- paste(indir, "/", WEIGHTS_FILE_START, lon, "_", lat, ".csv", sep="")
    loaded_series <- as.double(unlist(read.csv(series_file, header=FALSE)[2]))
}

loadCellUncertainties <- function(indir, lon, lat) {
    series_file <- paste(indir, "/", UNCERTAINTIES_FILE_START, lon, "_", lat, ".csv", sep="")
    loaded_series <- as.double(unlist(read.csv(series_file, header=FALSE)[2]))
}

writeSeries <- function(outdir, lon, lat, series, name) {
    filename <- paste(outdir, "/cell_", name, "_", lon, "_", lat, ".csv", sep="")
    sink(filename)
    for (i in 1:length(series)) {
        cat(i, ",", series[i], "\n", sep="")
    }
    sink()
}

writeParams <- function(outdir, lon, lat, params) {
    filename <- paste(outdir, "/cell_params_", lon, "_", lat, ".csv", sep="")
    sink(filename)
    for (i in 1:length(params)) {
        cat(params[i],"\n")
    }
    sink()
}

# Function to see if two fitted curves are significantly different.
# Significance is assessed as the correlation between the two curves over one year
# If r^2 < 0.99, the curves are deemed to be different. This limit has been established
# empirically
significantDifference <- function(curve1, curve2) {
    cat("Curve difference level =",cor(curve1, curve2)**2,"\n")
    return (cor(curve1, curve2)**2 < 0.99)
}

etopoLonIndex <- function(lon) {
    if (lon <= 20) {
        lon <- lon + 360
    }

    return (which(etopo_lons >= lon)[1])
}

etopoLatIndex <- function(lat) {
    return (which(etopo_lats >= lat)[1])
}


##########################################################

# Command line parameters
for (arg in commandArgs()) {
    argset <- strsplit(arg, "=", fixed=TRUE)
    if (!is.na(argset[[1]][2])) {
        if (argset[[1]][1] == "lon") {
            assign("lon",argset[[1]][2])
        } else if (argset[[1]][1] == "lat") {
            assign("lat",argset[[1]][2])
        } else if (argset[[1]][1] == "indir") {
            assign("indir",argset[[1]][2])
        } else if (argset[[1]][1] == "outdir") {
            assign("outdir",argset[[1]][2])
        }
    }
}

lon <- as.integer(lon)
lat <- as.integer(lat)

sink(paste(outdir,"/",lon,"_",lat,".log",sep=""))

cat("CELL",lon,lat,"\n")

nc <- nc_open("sea.nc")
sea <- ncvar_get(nc, "SEA")
nc_close(nc)

# If there's a final curve file, this means that the cell has been
# successfully processed. Copy the files to the new output directory.
final_curve_file <- paste(indir, "/cell_curve_", lon, "_", lat, ".csv", sep="")
if (file.exists(final_curve_file)) {
    cat("Cell already has SUCCESSFUL fit - copying files\n")

    file.copy(final_curve_file, outdir)
    file.copy(paste(indir, "/cell_series_", lon, "_", lat, ".csv", sep=""), outdir)
    file.copy(paste(indir, "/cell_weights_", lon, "_", lat, ".csv", sep=""), outdir)
    file.copy(paste(indir, "/cell_uncertainties_", lon, "_", lat, ".csv", sep=""), outdir)
    file.copy(paste(indir, "/cell_params_", lon, "_", lat, ".csv", sep=""), outdir)
} else if (sea[lon, lat] == 0) {
    cat("Land cell - skipping\n")
} else {
    # Load background data
    load("mean_directional_acfs.R")
    load("pco2_spatial_variation.R")
    temporal_acf <- read.csv("mean_temporal_acf.csv")[[2]]

    nc <- nc_open("etopo60.cdf")
    etopo <- ncvar_get(nc, "ROSE")
    etopo_lons <- ncvar_get(nc, "ETOPO60X")
    etopo_lats <- ncvar_get(nc, "ETOPO60Y")
    nc_close(nc)

    # Load the data
    in_series <- loadCellMeasurements(indir, lon, lat)
    in_weights <- loadCellWeights(indir, lon, lat)
    in_uncertainties <- loadCellUncertainties(indir, lon, lat)

    # Interpolation controls
    fit_found <- FALSE
    continue_fit <- TRUE
    spatial_interpolation_count <- 0

    # Stores a previously successful curve fit, and the measurements on which it was based
    # (which may include interpolated values)
    stored_measurements <- NULL
    stored_weights <- NULL
    stored_uncertainties <- NULL
    stored_fit_params <- NULL
    stored_curve <- NULL

    # Variables for processing
    current_measurements <- NULL
    current_weights <- NULL
    current_uncertainties <- NULL
    interpolated_measurements <- NULL
    interpolated_weights <- NULL
    interpolated_uncertainties <- NULL
    current_curve <- NULL

    spatial_interpolation_cells <- NULL

    while (continue_fit) {

        attempt_curve_fit <- TRUE

        # If this is the first run, simply use the input measurements
        if (spatial_interpolation_count <= 0) {
            current_measurements <- in_series
            current_weights <- in_weights
            current_uncertainties <- in_uncertainties
        } else {
            # Otherwise perform the spatial interpolation. The spatial interpolation is always based off the original measurements
            # The spatial interpolation method directly sets the interpolated value variables.
            doSpatialInterpolation(indir, lon, lat, in_series, in_weights, in_uncertainties, spatial_interpolation_count)

            # Add the interpolated measurements to the original measurements
            interpolated_measurements <- mergeValues(in_series, interpolated_measurements)
            interpolated_weights <- mergeValues(in_weights, interpolated_weights)
            interpolated_uncertainties <- mergeValues(in_uncertainties, interpolated_uncertainties)

            # If the interpolation hasn't added anything, don't bother attempting the fit. Simply loop round.
            if (sum(!is.na(current_measurements)) == sum(!is.na(interpolated_measurements))) {
                attempt_curve_fit <- FALSE
            } else {
                current_measurements <- mergeValues(current_measurements, interpolated_measurements)
                current_weights <- mergeValues(current_weights, interpolated_weights)
                current_uncertainties <- mergeValues(current_uncertainties, interpolated_uncertainties)

                # Remove any outliers from the interpolated time series
                if (sum(!is.na(current_measurements)) > 2) {

                    outliers_removed <- TRUE
                    while (outliers_removed == TRUE) {
                        outliers_removed <- FALSE
                        old_count <- sum(!is.na(current_measurements))

                        series_mean <- mean(current_measurements, na.rm=T)
                        stdev <- sd(current_measurements, na.rm=T)

                        current_measurements[current_measurements > (series_mean + (stdev * 3))] <- NA
                        current_measurements[current_measurements < (series_mean - (stdev * 3))] <- NA
                        current_weights[is.na(current_measurements)] <- NA
                        current_weights[is.na(current_measurements)] <- NA
                        current_uncertainties[is.na(current_measurements)] <- NA
                        current_uncertainties[is.na(current_measurements)] <- NA

                        new_count <- sum(!is.na(current_measurements))

                        if (new_count != old_count) {
                            outliers_removed <- TRUE
                        }
                    }
                }
            }
        }

        if (continue_fit) {

            if (attempt_curve_fit) {

                # Attempt to fit a curve to the data. If it succeeds, the stored values will have been
                # updated automatically, so we don't need to do it here.
                current_curve <- attemptFit(current_measurements, current_weights, current_uncertainties)

                if (!is.null(current_curve)) {
                    fit_found <- TRUE

                    if (is.null(stored_curve)) {

                        # No previous successful fit - store this fit
                        # Note that the measurements, weights and uncertainties
                        # have already been stored by the attemptFit function above.
                        stored_curve <- current_curve

                        # Now we attempt a spatial interpolation to see if this makes a difference.
                        if (is.null(spatial_interpolation_cells) || length(spatial_interpolation_cells) == 0) {
                            continue_fit <- FALSE
                        } else if (spatial_interpolation_count > 0 && spatial_interpolation_count == nrow(spatial_interpolation_cells)) {
                            continue_fit <- FALSE
                        } else {
                            spatial_interpolation_count <- spatial_interpolation_count + 1
                        }

                    } else {
                        if (significantDifference(stored_curve, current_curve)) {
                            # The interpolated fit is different to the previous fit.
                            # Continue the interpolation to refine the fit further
                            stored_measurements <- current_measurements
                            stored_weights <- current_weights
                            stored_uncertainties <- current_uncertainties
                            stored_curve <- current_curve

                            if (spatial_interpolation_count > 0) {
                                if (is.null(spatial_interpolation_cells) || length(spatial_interpolation_cells) == 0) {
                                    continue_fit <- FALSE
                                } else if (spatial_interpolation_count > 0 && spatial_interpolation_count == nrow(spatial_interpolation_cells)) {
                                    continue_fit <- FALSE
                                } else {
                                    spatial_interpolation_count <- spatial_interpolation_count + 1
                                }
                            } else {
                                spatial_interpolation_count <- spatial_interpolation_count + 1
                            }
                        } else {
                            # The spatial interpolation made no difference. We can stop
                            # with the previous fit
                            continue_fit <- FALSE
                        }
                    }
                }
            }

            if (!attempt_curve_fit || is.null(current_curve)) {
                # We didn't get a successful fit. Do a spatial interpolation to get more
                # measurements in play. Unless, of course, there are no more candidate cells
                # for interpolation
                if (spatial_interpolation_count > 0) {
                    if (is.null(spatial_interpolation_cells) || length(spatial_interpolation_cells) == 0) {
                        continue_fit <- FALSE
                    } else if (spatial_interpolation_count > 0 && spatial_interpolation_count == nrow(spatial_interpolation_cells)) {
                        continue_fit <- FALSE
                    } else {
                        spatial_interpolation_count <- spatial_interpolation_count + 1
                    }
                } else {
                    spatial_interpolation_count <- spatial_interpolation_count + 1
                }
            }
        }

        if (fit_found) {
            cat("SUCCESS\n")

            # Write the cell details to disk
            writeSeries(outdir, lon, lat, stored_measurements, "series")
            writeSeries(outdir, lon, lat, stored_weights, "weights")
            writeSeries(outdir, lon, lat, stored_uncertainties, "uncertainties")
            writeParams(outdir, lon, lat, stored_fit_params)
            writeSeries(outdir, lon, lat, stored_curve, "curve")
        } else {
            if (!continue_fit) {
                cat("FITTING FAILED. COPYING FILES FOR NEXT ITERATION\n")
                file.copy(paste(indir, "/cell_series_", lon, "_", lat, ".csv", sep=""), outdir)
                file.copy(paste(indir, "/cell_weights_", lon, "_", lat, ".csv", sep=""), outdir)
                file.copy(paste(indir, "/cell_uncertainties_", lon, "_", lat, ".csv", sep=""), outdir)
            }
        }
    }
}

sink()
