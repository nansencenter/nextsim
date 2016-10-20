# TODO #

## develop branch: serial code ##

  * DONE! > support hexagon
  * benchmarking: computational time in function of number of threads used (multi-threads) and the mesh resolution --- in progress
  * add documentayion: using for example doxygen --- priority (1 month)
  * forcast platform: work with Sylvain --- priority (1 month)
  * support cmake/automake: needed to make code compilation more generic --- middle term goal
  * support ParaView/VTK: visullization tools used in the most of complex physical models developed in C/C++ --- long term goal

## parallel branch: parallel code ##

  * DONE! > comparison between serial and parallel codes
  * DONE! > support hexagon
  * DONE! > support gmsh mesh partitioning in memory (avoid system call). It works well on my mac but crases on hexagon (to be investigated)
  * DONE! > merge develop into parallel (get recent changes from develop branch, mainly the physics)
  * benchmarking: peformance of the parallel code by increasing the number of processor cores (strong and weak scalability analyses) --- in progress
  * generalize the parallel code to be able to use it in serial (mpirun -np 1 bin/nextsim.exec --config-file=nextsim.cfg) --- in progress

## time estimate ##

Task list for Abdoulaye:

neXtWIM:
- get it to agree with matlab coupled code (2 weeks) --> check with Tim
- DONE! > reduce diffusion: at each wim time step interp 3 wave fields onto mesh, then loop over elements to do breaking on them using the wave fields (2 weeks)
- make sure it is running efficiently for large scale simulations (MPI needed? --> important to answer this question)
- DONE! > external forcing/realistic domains (1 month) --> Sylvain should lead, Abdoulaye should support

Publication:
write a paper related to nextsim C++ (2 months, in parallel with another works) --> start to write the outline of the paper
benchmarking: peformance of the parallel code by increasing the number of processor cores (strong and weak scalability analyses) --- in progress

Kara:
forcast platform: work with Sylvain --- priority (1 month) --> to be clarified during the Sylvain's visit

neXtSIM:
- generalisation of the code to be able to use it in serial (mpirun -np 1 bin/nextsim.exec --config-file=nextsim.cfg) --- in progress
- DONE! > support gmsh mesh partitioning in memory (avoid system call). It works well on my mac but crases on hexagon (to be investigated)
- add documentayion: using for example doxygen --- priority (1 month)
- support cmake/automake: needed to make code compilation more generic --- middle term goal
- support ParaView/VTK: visullization tools used in the most of complex physical models developed in C/C++ --- long term goal
