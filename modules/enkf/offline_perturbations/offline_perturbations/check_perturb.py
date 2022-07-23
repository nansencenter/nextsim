import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def read_and_plot(filename, varnames, figname, matrix_size):
    df = pd.read_csv(filename, header=None, sep='\s+')
    fig = plt.figure(1)
    for i, varname in enumerate(varnames):
        ax = fig.add_subplot(2, 3, i+1)
        data = np.reshape(df.iloc[:,i].to_numpy(), matrix_size);
        pc = ax.pcolor(data, shading='flat')
        fig.colorbar(pc, ax = ax)
        ax.set_title(f"noise on {varname}")
    fig.tight_layout() 
    fig.savefig(figname)

if __name__ == "__main__":
    filename = "synforc_2.dat"
    varnames=['u', 'v', 'snowfall', 'longwave radiation', 'SSS', 'SST']
    figname = 'display_result.png'
    matrix_size = (1024, 1024)
    read_and_plot(filename, varnames, figname, matrix_size)