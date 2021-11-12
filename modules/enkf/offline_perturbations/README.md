# Instruction
This is an offline FORTRAN version that can create sequential perturbations. It is derived from the online perturbation in neXtSIM modules/enkf/perturbation

__The size of the perturbation series is set in main_pseudo2D.F90, i_step.__

__Use run_script.sh to compile the fortran codes and execute to create perturbation__

To the FORTRAN code, the origianl purpose is to estimate the amplification of wind velocity by perturbing wind fields in ensemble simulations. See details in  run_script.sh.

The estimation of the amlification is in calculate_ratio.m 
Dignostic output of an perturbation is given in display_perturbation.m
