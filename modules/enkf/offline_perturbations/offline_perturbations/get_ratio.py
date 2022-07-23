import netCDF4
import datetime as dt
import os
import numpy as np
import pandas as pd

def Get_wind(filename):
    """Read wind from the ECMWF files
    """
    f = netCDF4.Dataset(filename, 'r')
    u10 = np.array(f['x_wind_10m'])
    v10 = np.array(f['y_wind_10m'])
    f.close()
    return u10, v10

def Get_perturb(filename):
    """Read wind from the perturbation files
    """
    perturb = pd.read_csv(filename, header=None, delim_whitespace=True).to_numpy()#, sep=r'\s{2,}', engine='python')
    du10 = perturb[:, 0]
    dv10 = perturb[:, 1]

    return du10, dv10

def Get_ratio(u10, v10, du10, dv10, r):
    u_norm = u10*u10 + v10*v10
    du10_view = du10.view()
    u10_purterb = u10 + (du10.view().reshape(len(r), -1)*r[:, np.newaxis]).view().reshape(du10.shape)
    v10_purterb = v10 + (dv10.view().reshape(len(r), -1)*r[:, np.newaxis]).view().reshape(dv10.shape)
    u_purterb_norm = u10_purterb*u10_purterb + v10_purterb*v10_purterb
    return np.mean((np.sqrt(u_norm)/np.mean(np.sqrt(u_purterb_norm), 1)).view().reshape(len(r), -1), -1)

def main(dir_wind, dir_perturb, ndays, N):

    # get dimension size from files
    starting_date = dt.date(2019, 9, 4)
    t = starting_date.strftime("%Y%m%d")
    filename_wind = f"generic_atm_{t}.nc"
    filename_wind = os.path.join(dir_wind, filename_wind)
    f = netCDF4.Dataset(filename_wind, 'r')
    nt = f.dimensions['time'].size
    nx = f.dimensions['x'].size
    ny = f.dimensions['y'].size
    f.close()

    # r = np.arange(0.5, 5.5, 0.5)
    r =np.array( [1.] )
    nr = len(r)
    #define array size
    du10 = np.zeros((nr, N, ny*nx))
    dv10 = np.zeros((nr, N, ny*nx))
    ratio = np.zeros((nr, nt*ndays))
    for iday in range(ndays):
        # get wind
        t = starting_date + dt.timedelta(days=iday)
        t = t.strftime("%Y%m%d")
        filename_wind = f"generic_atm_{t}.nc"
        filename_wind = os.path.join(dir_wind, filename_wind)
        u10, v10 = Get_wind(filename_wind)
        
        # get perturbation
        for it in range(nt):
            for ie in range(N):
                filename_perturb = os.path.join(dir_perturb, f"mem{ie+1}")
                # filename_perturb = dir_perturb
                filename_perturb = os.path.join(filename_perturb, f"synforc_randfld{4*iday + it + 1}.dat")
                du10[..., ie, :], dv10[..., ie, :] = Get_perturb(filename_perturb)
        
            ratio[:, 4*iday + it] = Get_ratio(u10[it].view().ravel(order='F'), 
                                               v10[it].view().ravel(order='F'), 
                                               du10, dv10, r)

            print(ratio[:, 4*iday + it])
    print(np.std(ratio))

main(dir_wind='/cluster/projects/nn2993k/sim/sukun_test/nextsim_data_dir/GENERIC_ATM', 
    dir_perturb='/cluster/projects/nn2993k/sim/sukun_test/nextsim_data_dir/perturbation/result', 
    ndays=245, N=40)
