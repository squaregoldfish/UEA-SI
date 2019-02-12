#!/usr/bin/python3

import shutil
import math
import numpy as np
from netCDF4 import Dataset


def calc_gas_xfer(wind_speed, sst):
    schmidt = 2073.1 - 125.62 * sst + 3.6276 * sst**2 - 0.043219 * sst**3
    return 0.266 * wind_speed**2 * (schmidt / 660)**-0.5

def calc_solubility(sst, salinity):
    kelvin_offset = 273.15
    a1 = -60.2409
    a2 = 93.4517
    a3 = 23.3585
    b1 = 0.023517
    b2 = -0.023656
    b3 = 0.0047036

    kelvin = sst + kelvin_offset

    calc_value = a1 + a2 * (100 / kelvin) + a3 * math.log(kelvin / 100) + \
        salinity * (b1 + b2 * (kelvin / 100) + b3 * (kelvin / 100)**2)

    return math.exp(calc_value)


# Copy the fco2 file ready to add new variables
shutil.copyfile("fco2.nc", "fco2_with_flux.nc")

# Load all input data
fco2 = Dataset("fco2.nc").variables["fco2"][:, :, :]
atmos_co2 = Dataset("atmos_co2.nc").variables["ATMOS_CO2"][:, :, :]
sst = Dataset("sst.nc").variables["SST"][:, :, :]
salinity = Dataset("salinity.nc").variables["salinity"][:, :, :]
wind = Dataset("wind.nc").variables["si10"][:, :, :]

flux = np.empty(shape=(396, 72, 144))
flux.fill(-9999.9)

area = np.empty(shape=(72, 144))


for lon in range(0, 144):
    lon1 = lon * 2.5
    lon2 = lon * 2.5 + 2.5
    for lat in range(0, 72):
        lat1 = math.radians(lat * 2.5 - 90)
        lat2 = math.radians(lat * 2.5 - 90 + 2.5)

        area[lat, lon] = (math.pi / 180) * 6371**2 * abs(math.sin(lat1) - math.sin(lat2)) * \
                         abs(lon1 - lon2)

        print("\r" + str(lon) + " " + str(lat) + "    ", end="")
        for time in range(0, 396):

            if sst[time, lat, lon] > 0:
                delta_pco2 = fco2[time, lat, lon] - atmos_co2[time, lat, lon]
                gas_xfer = calc_gas_xfer(wind[time, lat, lon], sst[time, lat, lon])
                solubility = calc_solubility(sst[time, lat, lon], salinity[time, lat, lon])
                gas_xchange = gas_xfer * solubility * 1.027 * 24 * 365 * 1e-5
                flux[time, lat, lon] = gas_xchange * delta_pco2
                if np.isnan(flux[time, lat, lon]):
                    flux[time, lat, lon] = -9999.9


flux_nc = Dataset("fco2_with_flux.nc", mode="a")

fluxvar = flux_nc.createVariable("flux", "d", ("time","lat", "lon"), fill_value = -9999.9)
fluxvar[:, :, :] = flux

areavar = flux_nc.createVariable("area", "d", ("lat", "lon"), fill_value = -9999.9)
areavar[:, :] = area

flux_nc.close()