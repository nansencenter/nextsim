# Instruction _7-Feb-2020_
This is a stand-alone FORTRAN version that can create sequential perturbations offline. The code is a subset of the original perturbation code obtained by Sukun. 
- A more complete version is in nextsim/modules/enkf/perturbation.
The online perturbation in neXtSIM modules/enkf/perturbation. The latest maintainence is in nextsim's git repository: IOperturbation/bf1686de6b8757da947448aaf9a0c4b5e22cd51c

- In ./src/main_pseudo2D.F90,     
    - Set the length of a sequential perturbations for one member to variable i_step.
    - Set domain size to xdim = 1024, ydim = 1024,xy_full = xdim*ydim

- THEN, compile the code by makefile in ./src.

- Configuration of perturbations are set i pseudo2D.nml

- In ./result folder, pertubation series are saved in subfolders distincted by ensemble id. For examples,
    -mem1 containts perturbations in netcdf as synforc_i.nc
    ncdump -h synforc_1.nc shows
        netcdf synforc_1 {
        dimensions:
            xy = 1048576 ;       (! xy =xdim*ydim)
        variables:               (! the following variables related variances are defined in pesudo.nml. uwind,vwind are horizontal wind speed in u,v directions. The other variables are independent (not correlated). One can add/reduce variables on the specific needs.)
            float uwind(xy) ;    
            float vwind(xy) ;
            float snowfall(xy) ;
            float Qlw_in(xy) ;
            float sss(xy) ;
            float sst(xy) ;
        }
- In ./report folder, it saves a document records previous studies. The estimation of the amplification is in ./report/get_ratio.m or .py

Use run_script.sh for a fresh compilation and generating perturbations, where ensemble size is given. The code can run in sequential or parallel on HPC
