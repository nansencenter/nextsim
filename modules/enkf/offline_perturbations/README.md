# Instruction _7-Feb-2020_
This is a stand-alone FORTRAN version that can create sequential perturbations offline. The code is a subset of the original perturbation code obtained by Sukun. 
- A more complete version is in nextsim/modules/enkf/perturbation.
The online perturbation in neXtSIM modules/enkf/perturbation. The latest maintainence is in nextsim's git repository: IOperturbation/bf1686de6b8757da947448aaf9a0c4b5e22cd51c

- In ./report folder, it saves previous relevant studies. The estimation of the amplification is in ./report/calculate_ratio.m 

- In ./src folder,     
    - Set the length of a sequential perturbations for one member to variable i_step in ./src/main_pseudo2D.F90.
    - Set domain size in main_pseudo2D.F90.
    - Compile the code by makefile.

- In ./result folder, pertubation series are saved in seperated subfolders distincted by ensemble id

- Configuration of perturbations are set i pseudo2D.nml

Use run_script.sh for a fresh compilation the fortran codes and generating perturbations in parallel by sbatch, where ensemble size is set.

slurm_script.sh is missing