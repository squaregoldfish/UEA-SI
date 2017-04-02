# NB REQUIRES THAT ALL HEADER LINES ARE REMOVED.

MONTH_STARTS <- c(0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334)
LEAP_MONTH_STARTS <- c(0, 31, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335)


# Fractions of a year for the start date of each week
DAY_STARTS <- seq(1, 365)

LEAP_DAY_STARTS <- vector(mode="numeric", length=365)

for (i in 1:365) {
    LEAP_DAY_STARTS[i] <- i + (((366 / 365) / 365) * (i - 1))
}

calcJDate <- function(year, month, day) {
    result <- 0

    if (isLeapYear(year)) {
        result <- LEAP_MONTH_STARTS[month]
    } else {
        result <- MONTH_STARTS[month]
    }

    result <- result + day

    return (result)
}


getDateIndex <- function(year, month, day) {

    jdate <- calcJDate(year, month, day)

    index <- NULL
    day_index = 0

    if (year >= 1985 && year <= 2015) {
        if (isLeapYear(year)) {
            day_index <- getDayIndex(jdate, DAY_STARTS)
        } else {
            day_index <- getDayIndex(jdate, LEAP_DAY_STARTS)
        }
    
        index <- ((year - 1985) * 365) + day_index
    }

    return (index)
}

getDayIndex <- function(date, days) {
    return  (tail(which(days < date), 1))
}

isLeapYear <- function(year) {
    leap_year <- FALSE
    if (year %% 4 == 0) {
        leap_year <- TRUE
        if (year %% 100 == 0) {
            if (year %% 400 != 0) {
                leap_year <- FALSE
            }
        }
    }

    return (leap_year)
}

getLatCell <- function(lat) {
    boundary <- trunc(lat / 2.5)
    cell <- boundary + 36
    if (lat > 0) {
        cell <- cell + 1
    }
    return (cell)
}

getLonCell <- function(lon) {
    result <- (trunc(lon / 2.5) + 1)
    if (result == 145) {
        result <- 144
    }
    return (result)
}

###############################################

totals <- array(0, c(144,72,11315))
counts <- array(0, c(144,72,11315))

#conn <- file("SOCATv4.tsv",open="r")
conn <- file("short.tsv",open="r")

line <- readLines(conn,n=1)
line_count <- 1

while (length(line) > 0) {

    if (line_count %% 1000 == 0) {
        cat("\r",line_count)
    }

    fields <- unlist(strsplit(line, "\t"))
    
    year <- as.numeric(fields[5])
    month <- as.numeric(fields[6])
    day <- as.numeric(fields[7])

    lat <- as.numeric(fields[12])
    lon <- as.numeric(fields[11])
    fco2 <- as.numeric(fields[24])

    if (!is.null(fco2) && !is.na(fco2)) {

        if (!is.null(year) && !is.null(month) && !is.null(day) && !is.null(lat) && !is.null(lon) && !is.null(fco2)) {
            if (!is.na(year) && !is.null(month) && !is.na(day) && !is.na(lat) && !is.na(lon) && !is.na(fco2)) {

                if (lon < 0) {
                    lon <- 360 - abs(lon)
                }

                date_index <- getDateIndex(year, month, day)

                if (!is.null(date_index)) {
                    lat_cell <- getLatCell(lat)
                    lon_cell <- getLonCell(lon)

                    totals[lon_cell, lat_cell, date_index] <- totals[lon_cell, lat_cell, date_index] + fco2
                    counts[lon_cell, lat_cell, date_index] <- counts[lon_cell, lat_cell, date_index] + 1
                }
            }
        }
    }

    line <- readLines(conn,n=1)
    line_count <- line_count + 1
}

cat("\n")

for (lon in 1:144) {
#for (lon in 80:80) {
    for (lat in 1:72) {
    #for (lat in 37:37) {
        cat("\r","Output file",lon,lat)


        output <- vector(mode="numeric", length=11315)
        output[output == 0] <- NA

        for (i in 1:11315) {
            if (counts[lon, lat, i] > 0) {
                output[i] <- totals[lon, lat, i] / counts[lon, lat, i]
            }
        }

        if (sum(!is.na(output)) > 2) {

            outliers_removed <- TRUE
            while (outliers_removed == TRUE) {
                outliers_removed <- FALSE
                old_count <- sum(!is.na(output))

                series_mean <- mean(output, na.rm=T)
                stdev <- sd(output, na.rm=T)

                output[output > (series_mean + (stdev * 3))] <- NA
                output[output < (series_mean - (stdev * 3))] <- NA

                new_count <- sum(!is.na(output))

                if (new_count != old_count) {
                    cat("  Outliers removed: ",lon,lat,old_count - new_count,"\n")
                    outliers_removed <- TRUE
                }
            }
        }

        
        out_file <- paste("cell_series_daily/cell_series_",lon,"_",lat,".csv",sep="")
        sink(out_file)
        for (i in 1:11315) {
            cat(i,",",output[i],"\n",sep="")
        }

        sink()
    }
}
cat("\n")
