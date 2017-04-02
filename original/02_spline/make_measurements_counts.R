# R!
library(ncdf)

REGION_NAMES <- c("so", "sa", "ea", "na", "io", "sp", "ep", "np")

nc <- open.ncdf("mask.nc")
mask <- get.var.ncdf(nc, "mask")
close.ncdf(nc)


for (r in 1:8) {

    counts <- vector(mode="numeric", length=324)

    for (lon in 1:144) {
        cat("\r",r,lon,"   ")
        for (lat in 1:72) {
            if (mask[lon, lat] == r) {
            measurements_file <- paste("output/measurements_",lon,"_",lat,".csv",sep="")
                if (file.exists(measurements_file)) {
                    measurements <- read.csv(measurements_file,header=F)[[2]]
                    for (t in 1:length(measurements)) {
                        if (!is.na(measurements[t])) {
                            counts[t] <- counts[t] + 1
                        }
                    }
                }
            }
        }
    }

    sink(paste(REGION_NAMES[r],"_counts.csv",sep=""))
    for (j in 1:length(counts)) {
        cat(j - 0.5,",-2.375,",counts[j],"\n",sep="")  
    }
    sink()

    cat("\n",max(counts),"\n")
}

cat("\n")
