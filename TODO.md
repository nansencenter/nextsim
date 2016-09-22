# TODO #

## develop branch: serial code ##

  * support hexagon --- done
  * benchmarking: computational time in function of number of threads used (multi-threads) and the mesh resolution --- in progress
  * add documentayion: using for example doxygen --- priority
  * forcast platform: work with Sylvain --- priority
  * support cmake/automake: needed to make code compilation more generic --- middle term goal
  * support ParaView/VTK: visullization tools used in the most of complex physical models developed in C/C++ --- long term goal

## parallel branch: parallel code ##

  * comparison between serial and parallel codes --- done
  * support hexagon --- done
  * support gmsh mesh partitioning in memory (avoid system call). It works well on my mac but crases on hexagon (to be investigated) --- priority
  * merge develop into parallel (get recent changes from develop branch, mainly the physics) --- done
  * benchmarking: peformance of the parallel code by increasing the number of processor cores (strong and weak scalability analyses) --- priority
  * generalize the parallel code to be able to use it in serial (mpirun -np 1 bin/nextsim.exec --config-file=nextsim.cfg) --- priority

## coupling code nextwim ##

  * get it to agree with matlab coupled code
  * reduce diffusion: at each wim time step interp 3 wave fields onto mesh, then loop over elements to do breaking on them using the wave fields
  * external forcing/realistic domains --- in progress
  * make sure it is running efficiently for large scale simulations (MPI needed?)
