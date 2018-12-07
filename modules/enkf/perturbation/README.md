# neXtSIM_FORCING

**Synoptic forcing routines from TOPAZ adapted to neXtSIM to generate ensemble**

How to:

_Generates a 360x360 2-dimensional random perturbation_

- A random seed is generated based on the machine clock to avoid generating the same field each
  time. 

- __Compiling__: the platform independent commands are in the makefile and expect machine-dependent
  options from make.inc 

- __Execution__: ./bin/p_pseudo2D

- __Postproc__: Ferret scripts are written to plot the outputs.

ToDo:

- utilize a namelist file for parameter input (grid size [idm,jdm], variance, etc...) 

- Write outputs in netcdf format 

- Non-gaussian perturbation

- Sensitivity analysis using only wind fields

__28.11.18__

- an main cpp  gen_pseudo_fields.cpp is added and complied with the fortran object

- it will probably be called in checkReloadDatasets in finiteelement.cpp or loadDataset in
  externaldata.cpp

- in checkReloadDatasets: load pseudo windfield and add it to the &M_wind


__03.12.18__

-  A branch in neXtSIM repo namely **ensemble** is created.

-  A new module in nextsim/modules/enkf for neXtSIM-EnKF

-  Scripts running an ensemble is in scripts/ensemble


__04.12.18__

-  Forcing and IC perturbations are in modules/enkf/perturbation

-  A cpp class called ensemble is being written. Whole EnKF related objects and functions can be
   constructed here.

-  The code now works only for ASR forcing. Write also for ECMWF.

-  Initiate an issue on github


__05.12.18__

- Two methods to compute_minmax and compute_vmean are added

- Struct data type is used for v.data, v.name and v.id 

- infile2.in should be read in limits_ranf routine in mod_random_forcing.F90


__06.12.18__

- neXtSIM compiled with libpseudo2D.so.1.0. Need to add its path to
  LD_LIBRARY_PATH in nextsim.src and link in the Makefile in the model.

__07.12.18__

- Compile with ensemble.cpp
