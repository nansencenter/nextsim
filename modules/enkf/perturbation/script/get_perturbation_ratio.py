
#  this script is used to calculate the ratio of averaged original velocity/perturbed velocity
import netCDF4
import datetime as dt
import os
import sys
import numpy as np
sys.path.append(os.path.join(os.getcwd(), 'src'))
import mod_perturb

def Get_wind(filename):
    """Read wind from the ECMWF files
    """
    f = netCDF4.Dataset(filename, 'r')
    print(f['x_wind_10m'].shape)
    u10 = np.array(f['x_wind_10m'])
    v10 = np.array(f['y_wind_10m'])
    f.close()
    print(u10.shape)
    return u10, v10

def Get_perturb(nx, ny, synforc, randfld, previous_perturbation_exist):
    """Read wind from the ECMWF files
    """
    mod_perturb.generate_perturbation(nx, 
                                    ny, 
                                       synforc,
                                       randfld, 
                                       previous_perturbation_exist)

    du10 = synforc[:, 0]
    dv10 = synforc[:, 1]
    return du10, dv10

def Get_ratio(u10, v10, du10, dv10):
    u_norm = u10*u10 + v10*v10
    u10_purterb = u10 + du10
    v10_purterb = v10 + dv10
    u_purterb_norm = u10_purterb*u10_purterb + v10_purterb*v10_purterb
    return np.sum(np.sqrt(u_norm))/np.sum(np.sqrt(u_purterb_norm))

def main(dir_wind, ndays, N):

    # get dimension size from files
    starting_date = dt.date(2019, 9, 14)
    t = starting_date.strftime("%Y%m%d")
    filename_wind = f"generic_atm_{t}.nc"
    filename_wind = os.path.join(dir_wind, filename_wind)
    f = netCDF4.Dataset(filename_wind, 'r')
    nt = f.dimensions['time'].size
    nx = f.dimensions['x'].size
    ny = f.dimensions['y'].size
    f.close()

    #define array size
    du10 = np.zeros((N, ny*nx))
    dv10 = np.zeros((N, ny*nx))
    ratio = np.zeros(nt*ndays)

    synforc = np.zeros((ny*nx, 4), order='F')
    randfld = np.zeros((ny*nx, 4), order='F')

    for iday in range(ndays):
        # get wind
        t = starting_date + dt.timedelta(days=iday)
        t = t.strftime("%Y%m%d")
        filename_wind = f"generic_atm_{t}.nc"
        filename_wind = os.path.join(dir_wind, filename_wind)
        u10, v10 = Get_wind(filename_wind)
        
        # get perturbation
        for it in range(ndays):
            for ie in range(N):
                du10[ie], dv10[ie] = Get_perturb(nx, ny, synforc, randfld, previous_perturbation_exist=4*iday + it)
        
            ratio[4*iday + it] = Get_ratio(u10[it].flatten(order='F'), 
                                               v10[it].flatten(order='F'), 
                                               du10, dv10)

            print('the ratio is:', ratio[4*iday + it])

main(dir_wind='/cluster/projects/nn2993k/sim/sukun_test/nextsim_data_dir/GENERIC_ATM',  
    ndays=2, N=40)
