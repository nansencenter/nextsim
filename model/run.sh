#!/bin/bash

#mpirun -np 4 bin/nextsim.exec --config-file=nextsim.cfg -log_summary -sub_pc_type cholesky -sub_pc_factor_mat_solver_package cholmod -sub_ksp_type preonly #-pc_factor_shift_type POSITIVE_DEFINITE -pc_factor_shift_type NONZERO

#mpirun -np 4 bin/nextsim.exec --config-file=nextsim.cfg #-log_summary #-sub_pc_type cholesky -sub_pc_factor_mat_solver_package cholmod -sub_ksp_type preonly

mpirun -np 4 bin/nextsim.exec --config-file=nextsim.cfg -sub_pc_type cholesky -sub_pc_factor_mat_solver_package cholmod -sub_ksp_type preonly
