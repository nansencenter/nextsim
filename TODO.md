# TODO #

## develop branch: serial code ##

  * support hexagon --- done
  * benchmarking: computational time in function of number of threads used (multi-threads) and the mesh resolution --- in progress
  * add documentayion: using for example doxygen --- middle term goal
  * support cmake/automake: needed to make code compilation more generic --- middle term goal
  * support ParaView/VTK: visullization tools used in the most of complex physical models developed in C/C++ --- long term goal

## parallel branch: parallel code ##

  * comparison between serial and parallel codes --- done
  * support hexagon --- done
  * support gmsh mesh partitioning in memory (avoid system call). It works well on my mac but crases on hexagon (to be investigated) --- priority
  * merge develop into parallel (get recent changes from develop branch, mainly the physics) --- priority
  * benchmarking: peformance of the parallel code by increasing the number of processor cores (strong and weel scalability analyses) --- priority
  * generalize the parallel code to be able to use it in serial (mpirun -np 1 bin/nextsim.exec --config-file=nextsim.cfg) --- priority
