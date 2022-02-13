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

- an main cpp  gen_pseudo_fields.cpp is added and compiled with the fortran object

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

- Two methods to computeMinMax and computeVecMean are added

- Struct data type is used for v.data, v.name and v.id

- infile2.in should be read in limits_ranf routine in mod_random_forcing.F90


__06.12.18__

- neXtSIM compiled with libpseudo2D.so.1.0. Need to add its path to
  LD_LIBRARY_PATH in nextsim.src and link in the Makefile in the model.

__07.12.18__

- Compile with ensemble.cpp


__10.12.18__

- To compile nextsim with any module including enkf, module should be copied to
  the image in Dockerfile i.e. COPY module $NEXTSIMDIR/module should be added.

- Added FFTW into Dockerfile

__12.12.18__

- Ensemble is run successfully for the first time. Time-correlation should be
  added to two perturbation fields. (done)

__13.12.18__

- Time correlated perturbations are generated in ran_update_ran1 in
  mod_random_forcing.F90. Write the new perturbations in there and save two
outputs to feed t-1 and t+1 of forcing. (done)

__14.12.18__

- grep 'CHeCK' to see where to improve in the code for generating a coupled time-correlated perturbation


__14.01.19__

- Doesn't write ranfld_prev and ranfld_next (done)

- Distribute read data to tiles (done)

__18.01.20__

- Read ranfld_prev.dat when restart do not generate from the scratch

- Check if the order of the vectors is the same in fortran output and c++ input


