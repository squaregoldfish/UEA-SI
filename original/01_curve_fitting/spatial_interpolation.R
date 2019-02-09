# Spatial interpolation methods for do_interpolation.R

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

filterAndSortInterpolationCells <- function(lon, lat, candidates, weights, uncertainties) {

    # The output of this is a table of lon, lat, weight, uncertainty.
    # Sorted by uncertainty, distance between cells, weight, lon, lat
    lons <- vector(mode="numeric", length=nrow(candidates))
    lons[lons == 0] <- NA

    lats <- lons
    weights <- lons
    uncertainties <- lons

    cell_count <- 0

    for (i in 1:nrow(candidates)) {

        cell_lon <- candidates[i, 1]
        cell_lat <- candidates[i, 2]

        if (!is.na(cell_lon) && !is.na(cell_lat)) {

            # If the cell is land, skip it
            if (sea[cell_lon, cell_lat] == 1) {
                # If there's land between the two cells, we can't go any further
                if (!landBetween(lon, lat, cell_lon, cell_lat)) {
                    # Get the spatial interpolation weight. If this is below the interpolation threshold,
                    # we discard the cell.
                    weight <- getSpatialInterpolationWeight(lon, lat, cell_lon, cell_lat)
                    if (weight >= SPATIAL_INTERPOLATION_WEIGHT_LIMIT) {
                    #if (weight >= 0) {

                        # Now get the uncertainty. The list of candidate cells will be sorted by this
                        uncertainty <- getSpatialInterpolationUncertainty(lon, lat, cell_lon, cell_lat)
                        cell_count <- cell_count + 1

                        lons[cell_count] <- cell_lon
                        lats[cell_count] <- cell_lat
                        weights[cell_count] <- weight
                        uncertainties[cell_count] <- uncertainty
                    }
                }
            }
        }
    }

    if (sum(!is.na(uncertainties)) == 0) {
        result <- vector(mode="numeric", length=0)
    } else {

        sorted_indices <- order(uncertainties, -weights, lons, lats)

        result <- vector(mode="numeric", length=(sum(!is.na(lons)) * 4))
        dim(result) <- c(sum(!is.na(lons)), 4)

        result[,1] <- lons[sorted_indices][!is.na(lons)]
        result[,2] <- lats[sorted_indices][!is.na(lons)]
        result[,3] <- weights[sorted_indices][!is.na(lons)]
        result[,4] <- uncertainties[sorted_indices][!is.na(lons)]
    }

    return (result)
}

landBetween <- function(lon, lat, interp_lon, interp_lat) {

    x1 <- getEtopoLon(lon)
    y1 <- getEtopoLat(lat)
    x2 <- getEtopoLon(interp_lon)
    y2 <- getEtopoLat(interp_lat)

    # Step through all cells in etopo between the two cells passed in
    # This uses an extension of the Bresenham algorithm by Eugen Dedu
    # http://lifc.univ-fcomte.fr/~dedu/projects/bresenham/index.html

    found_land <- FALSE

    i <- NULL
    ystep <- NULL
    xstep <- NULL
    error <- NULL
    errorprev <- NULL
    y <- y1
    x <- x1
    ddy <- NULL
    ddx <- NULL
    dx <- x2 - x1
    dy <- y2 - y1
    current_x <- x1
    current_y <- y1

    if (dy < 0) {
        ystep <- -1
        dy <- -dy
    } else {
        ystep <- 1
    }

    if (dx < 0) {
        xstep <- -1
        dx <- -dx
    } else {
        xstep <- 1
    }

    ddy <- 2 * dy
    ddx <- 2 * dx

    if (ddx >= ddy) {
        error <- dx
        errorprev <- dx

        for (i in 0:(dx - 1)) {

            # If we've already found a land cell, stop
            if (found_land) {
                break
            }

            x <- x + xstep
            error <- error + ddy

            if (error > ddx) {
                y <- y + ystep
                error <- error - ddx

                if (error + errorprev < ddx) {
                    current_y <- y - ystep
                    current_x <- x
                    found_land <- isEtopoLand(current_x, current_y)
                } else if (error + errorprev > ddx) {
                    current_y <- y
                    current_x <- x - xstep
                    found_land <- isEtopoLand(current_x, current_y)
                } else {
                    current_y <- y - ystep
                    current_x <- x
                    found_land <- isEtopoLand(current_x, current_y)

                    current_y <- y
                    current_x <- x - xstep
                    found_land <- isEtopoLand(current_x, current_y)
                }
            }

            current_y <- y
            current_x <- x
            found_land <- isEtopoLand(current_x, current_y)

            errorprev <- error
        }
    } else {

        error <- dy
        errorprev <- dy

        for (i in 0:(dy - 1)) {

            # If we've already found a land cell, stop
            if (found_land) {
                break
            }

            y <- y + ystep
            error <- error + ddx

            if (error > ddy) {
                x <- x + xstep
                error <- error - ddy

                if (error + errorprev < ddy) {
                    current_y <- y
                    current_x <- x - xstep
                    found_land <- isEtopoLand(current_x, current_y)
                } else if (error + errorprev > ddy) {
                    current_y <- y - ystep
                    current_x <- x
                    found_land <- isEtopoLand(current_x, current_y)
                } else {
                    current_y <- y
                    current_x <- x - xstep
                    found_land <- isEtopoLand(current_x, current_y)

                    current_y <- y - ystep
                    current_x <- x
                    found_land <- isEtopoLand(current_x, current_y)
                }
            }

            current_y <- y
            current_x <- x
            found_land <- isEtopoLand(current_x, current_y)

            errorprev <- error
        }
    }

    return (found_land)
}

isEtopoLand <- function(lon, lat) {
    return (etopo[lon, lat] >= 0)
}

getEtopoLon <- function(lon) {
    lon_degrees <- (lon * 2.5) - 1.25
    return (etopoLonIndex(floor(lon_degrees)))
}

getEtopoLat <- function(lat) {
    lat_degrees <- (lat * 2.5) - 91.25
    return (etopoLatIndex(floor(lat_degrees)))
}


doSpatialInterpolation <- function(indir, lon, lat, measurements, weights, uncertainties, cell_count) {

    cat("Performing spatial interpolation for cell",lon,lat,", cell",cell_count, "/", nrow(spatial_interpolation_cells), "\n")

    if (is.null(spatial_interpolation_cells)) {
        cat("Getting spatial interpolation cells\n")

        interpolation_candidates <- calculateInterpolationCells(lon, lat, 10)
        spatial_interpolation_cells <<- filterAndSortInterpolationCells(lon, lat, interpolation_candidates, weights, uncertainties)
    }

    if (length(spatial_interpolation_cells) == 0) {

        # There are no interpolation candidates. Return the input data
        interpolated_measurements <<- measurements
        interpolated_weights <<- weights
        interpolated_uncertainties <<- uncertainties
    } else {

        # Prepare data structures to hold the interpolated data. They will all be combined to
        # work out the final interpolated values
        interpolation_measurements <- vector(mode="numeric", length=(nrow(spatial_interpolation_cells) * length(measurements)))
        interpolation_measurements[interpolation_measurements == 0] <- NA
        dim(interpolation_measurements) <- c(nrow(spatial_interpolation_cells), length(measurements))

        interpolation_weights <- interpolation_measurements
        interpolation_uncertainties <- interpolation_measurements


        # Retrieve the data for each interpolation cell
        current_cell <- 0
        for (i in 1:cell_count) {
            interp_lon <- spatial_interpolation_cells[i, 1]
            interp_lat <- spatial_interpolation_cells[i, 2]
            interp_weight <- spatial_interpolation_cells[i, 3]
            interp_uncertainty <- spatial_interpolation_cells[i, 4]

            cat("INTERPOLATING TO",lon,lat,"from",interp_lon,interp_lat,"\n")

            # The measurements are easiest - they simply get loaded
            cell_measurements <- loadCellMeasurements(indir, interp_lon, interp_lat)

            # All measurements have a weight assigned to them already. They must be further weighted by the
            # distance over which they are being interpolated according to the spatial autocorrelation.
            cell_weights <- loadCellWeights(indir, interp_lon, interp_lat)

            current_cell <- current_cell + 1
            cell_weights <- cell_weights * interp_weight

            # All measurements have an uncertainty associated with them already. The uncertainty
            # must be increased as the values are interpolated.
            cell_uncertainties <- loadCellUncertainties(indir, interp_lon, interp_lat)

            cell_uncertainties <- sqrt(cell_uncertainties^2 + interp_uncertainty^2)

            # Add the interpolated cell's details to the data structure
            for (j in 1:length(measurements)) {
                interpolation_measurements[current_cell,j] <- cell_measurements[j]
                interpolation_weights[current_cell,j] <- cell_weights[j]
                interpolation_uncertainties[current_cell,j] <- cell_uncertainties[j]
            }
        }

        # Now combine all the interpolated cells with the original measurements to make a single time series
        output_measurements <- vector(mode="numeric", length=length(measurements))
        output_measurements[output_measurements == 0] <- NA

        output_weights <- output_measurements
        output_uncertainties <- output_measurements

        for (i in 1:length(measurements)) {

            # Calculate the weighted mean measurement and mean uncertainty from all interpolated cells
            weighted_mean <- sum(interpolation_measurements[,i] * interpolation_weights[,i],na.rm=T) / sum(interpolation_weights[,i],na.rm=T)
            mean_weight <- mean(interpolation_weights[,i],na.rm=T)
            mean_uncertainty <- sqrt(sqrt(sum(interpolation_uncertainties[,i]**2,na.rm=T)) / sum(!is.na(interpolation_uncertainties[,i]))^2 + 2.5^2)
            output_measurements[i] <- weighted_mean
            output_weights[i] <- mean_weight
            output_uncertainties[i] <- mean_uncertainty
        }

        # Push the interpolated values up to the calling function
        interpolated_measurements <<- output_measurements
        interpolated_weights <<- output_weights
        interpolated_uncertainties <<- output_uncertainties
    }
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


# Retrieve the uncertainty for the interpolation between two cells
getSpatialInterpolationUncertainty <- function(target_lon, target_lat, interp_lon, interp_lat) {

    uncertainty <- NA

    # Here we use the pco2_spatial_variation obect loaded at the start of the program.
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
        # mean uncertainty from the cells surrounding the target cell
        uncertainty_cells <- calculateInterpolationCells(target_lon, target_lat, 1)
        uncertainties <- vector(mode="numeric", length=nrow(uncertainty_cells))

        for (cell_loop in 1:nrow(uncertainty_cells)) {
            if (!is.na(uncertainty_cells[cell_loop,1])) {
                cell_lon <- uncertainty_cells[cell_loop,1]
                cell_lat <- uncertainty_cells[cell_loop,2]

                uncertainties[cell_loop] <- pco2_spatial_variation[target_lon, target_lat, cell_lon, cell_lat]
            }
        }
        uncertainty <- mean(uncertainties,na.rm=T)
    }
    return (uncertainty)
}


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

