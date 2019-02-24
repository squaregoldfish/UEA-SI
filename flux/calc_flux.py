#!/usr/bin/python3

import shutil
import math
import numpy as np
from netCDF4 import Dataset


def calc_schmidt(sst):
    return 2073.1 - 125.62 * sst + 3.6276 * sst**2 - 0.043219 * sst**3

def calc_solubility(sst, salinity):
    kelvin_offset = 273.15
    a1 = -58.0931
    a2 = 90.5069
    a3 = 22.2940
    b1 = 0.027766
    b2 = -0.025888
    b3 = 0.0050578

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

delta_pco2_nc = np.empty(shape=(396, 72, 144))
delta_pco2_nc.fill(-9999.9)

wind_sq_nc = np.empty(shape=(396, 72, 144))
wind_sq_nc.fill(-9999.9)

solubility_nc = np.empty(shape=(396, 72, 144))
solubility_nc.fill(-9999.9)

schmidt_nc = np.empty(shape=(396, 72, 144))
schmidt_nc.fill(-9999.9)

gas_xfer_nc = np.empty(shape=(396, 72, 144))
gas_xfer_nc.fill(-9999.9)

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
               delta_pco2_nc[time, lat, lon] = delta_pco2

               wind_sq = wind[time, lat, lon]**2
               wind_sq_nc[time, lat, lon] = wind_sq

               solubility = calc_solubility(sst[time, lat, lon], salinity[time, lat, lon])
               solubility_nc[time, lat, lon] = solubility

               schmidt = calc_schmidt(sst[time, lat, lon])
               schmidt_nc[time, lat, lon] = schmidt

               gas_xfer = 0.585 * solubility * (schmidt)**0.5 * wind_sq
               gas_xfer_nc[time, lat, lon] = gas_xfer

               flux[time, lat, lon] = gas_xfer * delta_pco2
               if np.isnan(flux[time, lat, lon]):
                   flux[time, lat, lon] = -9999.9


flux_nc = Dataset("fco2_with_flux.nc", mode="a")

dpco2var = flux_nc.createVariable("dpco2", "d", ("time","lat", "lon"), fill_value = -9999.9)
dpco2var[:, :, :] = delta_pco2_nc

windsqvar = flux_nc.createVariable("wind_sq", "d", ("time","lat", "lon"), fill_value = -9999.9)
windsqvar[:, :, :] = wind_sq_nc

solubilityvar = flux_nc.createVariable("solubility", "d", ("time","lat", "lon"), fill_value = -9999.9)
solubilityvar[:, :, :] = solubility_nc

schmidtvar = flux_nc.createVariable("schmidt", "d", ("time","lat", "lon"), fill_value = -9999.9)
schmidtvar[:, :, :] = schmidt_nc

gasxfervar = flux_nc.createVariable("gas_xfer", "d", ("time","lat", "lon"), fill_value = -9999.9)
gasxfervar[:, :, :] = gas_xfer_nc

fluxvar = flux_nc.createVariable("flux", "d", ("time","lat", "lon"), fill_value = -9999.9)
fluxvar[:, :, :] = flux

areavar = flux_nc.createVariable("area", "d", ("lat", "lon"), fill_value = -9999.9)
areavar[:, :] = area

flux_nc.close()