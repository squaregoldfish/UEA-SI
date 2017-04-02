# R!
library(ncdf4)

LON_SIZE <- 144
LAT_SIZE <- 72
LON_CELL_SIZE <- 2.5
LAT_CELL_SIZE <- 2.5
SERIES_LENGTH <- 324
EARTH_RADIUS <- 6367.5
SPATIAL_ACF_LAG_STEP <- 25
SPATIAL_INTERPOLATION_WEIGHT_LIMIT <- 0.367879

# Calculate the bearing between two cells, on rhumb lines
# Answer returned in degrees
calcCellBearing <- function(from_lon, from_lat, to_lon, to_lat) {

    # First, work out the cell sizes
    cell_lat_size <- pi / LAT_SIZE
    cell_lon_size <- (2 * pi) / LON_SIZE

    # Convert the cell indeces to radians
    lat1_rad <- (from_lat * cell_lat_size) - (pi / 2) - (cell_lat_size / 2)
    lon1_rad <- (from_lon * cell_lon_size) - (cell_lon_size / 2)
    lat2_rad <- (to_lat * cell_lat_size) - (pi / 2) - (cell_lat_size / 2)
    lon2_rad <- (to_lon * cell_lon_size) - (cell_lon_size / 2)

    # And calculate the distance between them
    d_lat <- lat2_rad - lat1_rad
    d_lon <- lon2_rad - lon1_rad

    d_phi <- log(tan(lat2_rad / 2 + pi / 4) / tan(lat1_rad / 2 + pi / 4))
    q <- 0
    if (abs(d_lat) > 1e-10) {
        q <- d_lat / d_phi
    } else {
        q <- cos(lat1_rad)
    }

    if (abs(d_lon) > pi) {
        if (d_lon > 0) {
            d_lon <- (2 * pi - d_lon) * -1
        } else {
            d_lon <- 2 * pi * d_lon
        }
    }

    bearing <- as.integer(((atan2(d_lon, d_phi)) / (pi / 180)) + 360) %% 360

    return (bearing)
}

# Get the spatial ACF direction for a bearing
getACFDirection <- function(bearing) {
    direction <- 0

    if (bearing <= 22.5 || bearing >= 337.5 || (bearing >= 157.5 && bearing <= 202.5)) {
        direction <- 1
    } else if ((bearing >= 67.5 && bearing <= 112.5) || (bearing >= 247.5 && bearing <= 292.5)) {
        direction <- 2
    } else if ((bearing >= 22.5 && bearing <= 67.5) || (bearing >= 202.5 && bearing <= 247.5)) {
        direction <- 3
    } else if ((bearing >= 112.5 && bearing <= 157.5) || (bearing >= 292.5 && bearing <= 337.5)) {
        direction <- 4
    }

    return (direction)
}

                        
# Calculate the distance between two cells in km
# Everything here is done in radians
calcCellDistance <- function(lon1, lat1, lon2, lat2) {

    # First, work out the cell sizes
    cell_lat_size <- pi / LAT_SIZE
    cell_lon_size <- (2 * pi) / LON_SIZE

    # Now convert the cell indeces to positions
    lat1_rad <- (lat1 * cell_lat_size) - (pi / 2) - (cell_lat_size / 2)
    lon1_rad <- (lon1 * cell_lon_size) - (cell_lon_size / 2)
    lat2_rad <- (lat2 * cell_lat_size) - (pi / 2) - (cell_lat_size / 2)
    lon2_rad <- (lon2 * cell_lon_size) - (cell_lon_size / 2)

    # And calculate the distance between them
    d_lat <- lat2_rad - lat1_rad
    d_lon <- lon2_rad - lon1_rad

    a <- sin(d_lat / 2) * sin(d_lat / 2) + cos(lat1_rad) * cos(lat2_rad) * sin(d_lon / 2) * sin(d_lon / 2)
    c <- 2 * atan2(sqrt(a), sqrt(1 - a))
    distance <- c * EARTH_RADIUS

    return (distance)
}


makeMeanValue <- function(values, weights) {

    weighted_total <- 0
    total_weight <- 0

    for (i in 1:length(values)) {
        if (!is.na(values[i])) {
            weighted_total <- weighted_total + (values[i] * weights[i])
            total_weight <- total_weight + weights[i]
        }
    }

    return (weighted_total / total_weight)

}

makeMeanUncertainty <- function(uncertainties) {

    total <- 0
    count <- 0

    for (i in 1:length(uncertainties)) {
        if (!is.na(uncertainties[i])) {
            total <- total + uncertainties[i] ** 2
            count <- count + 1
        }
    }

    return (sqrt(total / count))

}

# Retrieve the uncertainty for the interpolation between two cells
getSpatialInterpolationUncertainty <- function(target_lon, target_lat, interp_lon, interp_lat) {

    uncertainty <- NA

    # Here we use the pco2_spatial_variation object loaded at the start of the program.
    # This is a 4D array of <target_lon, target_lat, interp_lon, interp_lat>
    # The values in this object are the mean difference in pCO2 between two cells
    # within a small time period.

    # We can use the uncertainty from both the interp and target cells, since they are
    # are in the same direction and therefore equivalent
    forward_uncertainty <- pco2_spatial_variation[target_lon, target_lat, interp_lon, interp_lat]
    reverse_uncertainty <- pco2_spatial_variation[interp_lon, interp_lat, target_lon, target_lat]
    if (!is.na(forward_uncertainty) && !is.na(reverse_uncertainty)) {
        uncertainty <- (forward_uncertainty + reverse_uncertainty) /2
    } else if (!is.na(forward_uncertainty)) {
        uncertainty <- forward_uncertainty
    } else if (!is.na(reverse_uncertainty)) {
        uncertainty <- reverse_uncertainty
    } else {

        # If we can't find an uncertainty for the specific cell combination, use the
        # mean uncertainty from the cells surrounding the interpolation cell.
        # This keeps expanding until we get something we can use. It may result
        # in higher uncertainties than are strictly necessary, but there's not much
        # to be done about that.
        #
        interpolation_step <- 0
        uncertainties <- NA

        while (sum(!is.na(uncertainties)) == 0) {
            interpolation_step <- interpolation_step + 1


            uncertainty_cells <- calculateInterpolationCells(target_lon, target_lat, interpolation_step)

            uncertainties <- vector(mode="numeric", length=nrow(uncertainty_cells))

            for (cell_loop in 1:nrow(uncertainty_cells)) {
                if (!is.na(uncertainty_cells[cell_loop,1])) {
                    cell_lon <- uncertainty_cells[cell_loop,1]
                    cell_lat <- uncertainty_cells[cell_loop,2]

                    forward_uncertainty <- pco2_spatial_variation[target_lon, target_lat, interp_lon, interp_lat]
                    reverse_uncertainty <- pco2_spatial_variation[interp_lon, interp_lat, target_lon, target_lat]
                    if (!is.na(forward_uncertainty) && !is.na(reverse_uncertainty)) {
                        uncertainties[cell_loop] <- (forward_uncertainty + reverse_uncertainty) /2
                    } else if (!is.na(forward_uncertainty)) {
                        uncertainties[cell_loop] <- forward_uncertainty
                    } else if (!is.na(reverse_uncertainty)) {
                        uncertainties[cell_loop] <- reverse_uncertainty
                    }
                }
            }
        }

        uncertainty <- mean(uncertainties,na.rm=T)
    }
    return (uncertainty)
}

# Calculate the weighting to be used for spatial interpolation between two cells
getSpatialInterpolationWeight <- function(target_lon, target_lat, interp_lon, interp_lat) {

    weight <- 0

    # Here we use the mean_directional_acfs object loaded at the start of the program.
    # This is a 4D array of lat, lon, direction, distance
    # Distances are in bins of ACF_LAG_STEP km. We calculate the weighting from the
    # centre of the two cells.
    target_centre_lon <- (target_lon * LON_CELL_SIZE) - (LON_CELL_SIZE / 2)
    target_centre_lat <- (target_lat * LAT_CELL_SIZE) - (90 + (LAT_CELL_SIZE / 2))
    interp_centre_lon <- (interp_lon * LON_CELL_SIZE) - (LON_CELL_SIZE / 2)
    interp_centre_lat <- (interp_lat * LAT_CELL_SIZE) - (90 + (LAT_CELL_SIZE / 2))

    # We need to know the distance and bearing between the two grid cells
    distance <- calcCellDistance(target_lon, target_lat, interp_lon, interp_lat)

    # There's a bug near the poles which means we sometimes disappear off the globe.
    # This only happens if the ACF weight is also off the scale, so we can just ignore
    # this eventuality.
    if (!is.na(distance)) {
        distance_bin <- trunc(distance / SPATIAL_ACF_LAG_STEP) + 1
        bearing <- calcCellBearing(target_lon, target_lat, interp_lon, interp_lat)
        acf_direction <- getACFDirection(bearing)

        # We can use the ACF from both the interp and target cells, since they are
        # are in the same direction and therefore equivalent
        acf_value <- mean_directional_acfs[target_lat, target_lon, acf_direction, distance_bin]
        reverse_acf_value <- mean_directional_acfs[interp_lat, interp_lon, acf_direction, distance_bin]

        if (!is.na(acf_value) && !is.na(reverse_acf_value)) {
            weight <- (acf_value + reverse_acf_value) / 2
        } else if (!is.na(acf_value)) {
            weight <- acf_value
        } else if (!is.na(reverse_acf_value)) {
            weight <- reverse_acf_value
        } else {
            # If there are no directional ACF values, use the omnidirectional value
            # for the target cell only. If this doesn't exist, we have to return a zero
            # weight because we can't guess at whether or not the two cells' values are related
            weight <- mean_directional_acfs[target_lat, target_lon, 5, distance_bin]
            if (is.na(weight)) {
                weight <- 0
            }
        }
    }

    # Anything below the weighting limit is given a weighting of zero
    if (weight < SPATIAL_INTERPOLATION_WEIGHT_LIMIT) {
        weight <- 0
    }

    return (weight)
}


calcSpatialInterpXCell <- function(lon, x) {
    new_x <- lon + x
    if (new_x > LON_SIZE) {
        new_x <- new_x - LON_SIZE
    } else if (new_x < 1) {
        new_x <- LON_SIZE - abs(new_x)
    }

    return (new_x)
}


calcSpatialInterpYCell <- function(lat, y) {
    new_y <- lat + y
    if (new_y > LAT_SIZE || new_y < 1) {
        new_y <- NA
    }
    return (new_y)
}


calculateInterpolationCells <- function(lon, lat, step) {

    # Work out how many interpolation cells there will be.
    # Each step means (8 * step) cells are added to the list.
    # The number of cells is therefore the triangular number of steps * 8
    # I have discovered a truly marvellous proof of why this is the case,
    # but this comment is too brief to contain it.
    cell_count <- ((step * (step + 1)) / 2) * 8

    cells <- vector(mode="numeric", length=cell_count * 2)
    cells[cells == 0] <- NA
    dim(cells) <- c(cell_count,2)

    cell_count <- 0
    for (step_loop in 1:step) {
        for (y in (step_loop * -1):step_loop) {

            cell_y <- calcSpatialInterpYCell(lat, y)
            if (!is.na(cell_y)) {

                # For the top and bottom rows of the step grid, add all horizontal cells
                if (abs(y) == step_loop) {
                    for (x in (step_loop * -1):step_loop) {
                        cell_count <- cell_count + 1
                        cells[cell_count,1] <- calcSpatialInterpXCell(lon, x)
                        cells[cell_count,2] <- cell_y
                    }
                } else {
                    # For all other rows, just add the left and right edges
                    cell_count <- cell_count + 1
                    cells[cell_count,1] <- calcSpatialInterpXCell(lon, (step_loop * -1))
                    cells[cell_count,2] <- cell_y

                    cell_count <- cell_count + 1
                    cells[cell_count,1] <- calcSpatialInterpXCell(lon, step_loop)
                    cells[cell_count,2] <- cell_y
                }
            }
        }
    }

    return (cells)
}





# Load data
cat("Loading data...")
nc <- nc_open("sea.nc")
sea <- ncvar_get(nc, "SEA")
nc_close(nc)


load("mean_directional_acfs.R")
load("pco2_spatial_variation.R")

cat("done\n")


failed_count <- 1000000
last_failed_count <- 0

while (failed_count > 0) {

    cat("ROUND\n")

    # Reset the failure counter
    last_failed_count <- failed_count
    failed_count <- 0
    dir.create("output")

    for (lon in 1:144) {
    #for (lon in 2:2) {
        for (lat in 1:72) {
        #for (lat in 36:36) {

            # Set up filenames
            in_spline_file <- paste("input/spline_",lon,"_",lat,".csv",sep="")
            in_uncertainties_file <- paste("input/uncertainty_",lon,"_",lat,".csv",sep="")

            out_series_file <- paste("output/spline_",lon,"_",lat,".csv",sep="")
            out_uncertainties_file <- paste("output/uncertainty_",lon,"_",lat,".csv",sep="")

            cat("\r",lon,lat," ")


            # If the output file already exists, we don't need to process this cell
            if (!file.exists(out_series_file)) {

                if (file.exists(in_spline_file)) {
                    # We already have a complete time series for this cell - just copy the files
                    file.copy(in_spline_file, out_series_file)
                    file.copy(in_uncertainties_file, out_uncertainties_file)
                    cat("already done    ")
                } else {

                    # If this isn't a sea cell, skip it
                    if (sea[lon, lat] == 1) {


                        # Work out which interpolation cells we need to look at
                        neighbour_cells <- calculateInterpolationCells(lon, lat, 1)

                        # Set up variables for interpolated data values
                        interpolation_spline <- vector(mode="numeric", length=(nrow(neighbour_cells) * SERIES_LENGTH))
                        interpolation_spline[interpolation_spline == 0] <- NA
                        dim(interpolation_spline) <- c(nrow(neighbour_cells), SERIES_LENGTH)

                        interpolation_weight <- interpolation_spline
                        interpolation_uncertainty <- interpolation_spline


                        # Get the data from the neighbour cells
                        for (cell_loop in 1:nrow(neighbour_cells)) {
                            neighbour_lon <- neighbour_cells[cell_loop, 1]
                            neighbour_lat <- neighbour_cells[cell_loop, 2]

                            neighbour_spline_file <- paste("input/spline_",neighbour_lon,"_",neighbour_lat,".csv",sep="")
                            neighbour_uncertainty_file <- paste("input/uncertainty_",neighbour_lon,"_",neighbour_lat,".csv",sep="")

                            if (file.exists(neighbour_spline_file)) {
                                # Get the interpolation weight and uncertainty for this cell
                                #
                                interp_weight <- getSpatialInterpolationWeight(lon, lat, neighbour_lon, neighbour_lat)
                                if (interp_weight == 0) {
                                    interp_weight <- 0.01
                                }
                                interp_uncertainty <- getSpatialInterpolationUncertainty(lon, lat, neighbour_lon, neighbour_lat)

                                neighbour_data <- read.csv(neighbour_spline_file,header=F)[[2]]
                                neighbour_uncertainty <- read.csv(neighbour_uncertainty_file,header=F)[[2]]

                                for (i in 1:SERIES_LENGTH) {
                                    interpolation_spline[cell_loop, i] <- neighbour_data[i]
                                    interpolation_uncertainty[cell_loop, i] <- neighbour_uncertainty[i] + interp_uncertainty
                                    interpolation_weight[cell_loop, i] <- interp_weight
                                }
                            }
                        }

                        if (sum(!is.na(interpolation_spline)) == 0) {
                            # There is no data to interpolate. Log a failure.
                            failed_count <- failed_count + 1
                            cat(" failed")
                        } else {
                            cat(" interpolated")
                            sink(out_series_file)
                            for (i in 1:SERIES_LENGTH) {
                                cat(i,",",makeMeanValue(interpolation_spline[, i], interpolation_weight[, i]),"\n",sep="")
                            }
                            sink()

                            sink(out_uncertainties_file)
                            for (i in 1:SERIES_LENGTH) {
                                cat(i,",",makeMeanUncertainty(interpolation_uncertainty[, i]),"\n",sep="")
                            }
                            sink()
                        }
                    }
                }
            }
        }
    }

    cat("\nFailed cells =",failed_count,"\n")

    if (failed_count == last_failed_count) {
        cat("Cannot fix any more cells. Stopping.\n")
        break
    } else {

        # Copy the output files back into the input for the next round
        for (i in list.files("output")) {
            file.rename(paste("output/",i,sep=""), paste("input/",i,sep=""))
        }
        unlink("output", recursive=TRUE)
    }
}



# Now we've done the best we can,
# fill in the remaining cells as follows:
#
# Copy the values from the neighbour cells as a weighted mean.
# Set the uncertainty to equal that of the surrounding cells.
# This is because we can't tell what the uncertainty should actually be,
# and doing anything else just gives silly results
#
# In the process, this will put all output files in the output directory.
cat("Fixing up the cells that failed\n")
for (lon in 1:144) {
    for (lat in 1:72) {
        cat("\r",lon,lat," ")

        # Set up filenames
        in_spline_file <- paste("input/spline_",lon,"_",lat,".csv",sep="")
        in_uncertainties_file <- paste("input/uncertainty_",lon,"_",lat,".csv",sep="")

        out_series_file <- paste("output/spline_",lon,"_",lat,".csv",sep="")
        out_uncertainties_file <- paste("output/uncertainty_",lon,"_",lat,".csv",sep="")

        if (file.exists(in_spline_file)) {
            #We already have a complete time series for this cell - just copy the files
            file.copy(in_spline_file, out_series_file)
            file.copy(in_uncertainties_file, out_uncertainties_file)
            cat("already done")
        } else {

            if (sea[lon, lat] == 1) {

                cat("fixing")

                failed <- TRUE
                step_count <- 0

                while (failed) {

                    step_count <- step_count + 1

                    neighbour_cells <- calculateInterpolationCells(lon, lat, step_count)

                    # Set up variables for interpolated data values
                    interpolation_spline <- vector(mode="numeric", length=(nrow(neighbour_cells) * SERIES_LENGTH))
                    interpolation_spline[interpolation_spline == 0] <- NA
                    dim(interpolation_spline) <- c(nrow(neighbour_cells), SERIES_LENGTH)

                    interpolation_weight <- interpolation_spline
                    interpolation_uncertainty <- interpolation_spline

                    # Get the data from the neighbour cells
                    for (cell_loop in 1:nrow(neighbour_cells)) {
                        neighbour_lon <- neighbour_cells[cell_loop, 1]
                        neighbour_lat <- neighbour_cells[cell_loop, 2]

                        neighbour_spline_file <- paste("input/spline_",neighbour_lon,"_",neighbour_lat,".csv",sep="")
                        neighbour_uncertainty_file <- paste("input/uncertainty_",neighbour_lon,"_",neighbour_lat,".csv",sep="")

                        if (file.exists(neighbour_spline_file)) {

                            # Get the interpolation weight and uncertainty for this cell
                            interp_weight <- getSpatialInterpolationWeight(lon, lat, neighbour_lon, neighbour_lat)
                            if (interp_weight == 0) {
                                interp_weight <- 0.01
                            }


                            neighbour_data <- read.csv(neighbour_spline_file,header=F)[[2]]
                            neighbour_uncertainty <- read.csv(neighbour_uncertainty_file,header=F)[[2]]

                            for (i in 1:SERIES_LENGTH) {
                                interpolation_spline[cell_loop, i] <- neighbour_data[i]
                                interpolation_uncertainty[cell_loop, i] <- neighbour_uncertainty[i]
                                interpolation_weight[cell_loop, i] <- interp_weight
                            }
                        }
                    }

                    if (sum(!is.na(interpolation_spline)) == 0) {
                        # There is no data to interpolate. Log a failure.
                        failed <- TRUE
                    } else {
                        sink(out_series_file)
                        for (i in 1:SERIES_LENGTH) {
                            cat(i,",",makeMeanValue(interpolation_spline[, i], interpolation_weight[, i]),"\n",sep="")
                        }
                        sink()

                        sink(out_uncertainties_file)
                        for (i in 1:SERIES_LENGTH) {
                            cat(i,",",makeMeanUncertainty(interpolation_uncertainty[, i]),"\n",sep="")
                        }
                        sink()
                        failed <- FALSE
                        cat(step_count)
                    }
                }
            }
        }
    }
}

cat("\n")

