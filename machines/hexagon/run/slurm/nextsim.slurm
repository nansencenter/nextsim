#!/bin/bash
#SBATCH --job-name="mpijob"                   # Give a name (optional) to the job
#SBATCH --account=nn2993k                     # Specify the project the job should be accounted on
#SBATCH --ntasks=64                           # Number of cores needed for the job
#SBATCH --ntasks-per-node=16                  # Number of tasks to be launched per node
##SBATCH --mem=16000                          # Specify the real memory required per node (in MB)
##SBATCH --mem-per-cpu=2000                   # Specify mimimum memory required per allocated CPU (in MB)
##SBATCH --nodes=32                           # Request a minimum and maximum nodes be allocated
#SBATCH --time=60:00:00                       # Total run time (in format HH:MM:SS) for the job
#SBATCH --begin=now+0hour                     # Defer the job allocation until the specified time
#SBATCH --mail-type=ALL                       # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=abdoulaye.samake@nersc.no # email to the user
#SBATCH --output=nextsim_%j.out               # Standard output
#SBATCH --error=nextsim_%j.err                # Standard error
source $HOME/.bash_profile
export CONFIG_DIR=`pwd`
cd $NEXTSIMDIR/model/bin
mpirun --bind-to-core -x LD_LIBRARY_PATH nextsim.exec \
	   --config-file=$CONFIG_DIR/nextsim.cfg \
	   -sub_pc_type cholesky \
	   -sub_pc_factor_mat_solver_package cholmod \
	   -sub_ksp_type preonly \
	   -mat_mumps_icntl_14 20 \
	   2>&1 > $WORK/Logs
