# Instruction _31-March-2022_
This is a stand-alone FORTRAN version that can create sequential perturbations offline. 

- This offline version only keeps the necessary code for generating perturbations. The online perturbation in neXtSIM modules/enkf/perturbation. The latest maintainence is in nextsim's git repository: IOperturbation/bf1686de6b8757da947448aaf9a0c4b5e22cd51c. Both code are subsets of the perturbation code, originally from HYCOM https://svn.nersc.no/hycom/browser/HYCOM_2.2.37/CodeOnly/src_2.2.37/nersc/mod_random_forcing.F. About details about the developmet in the original code could contact Laurent Bertino laurent.bertino@nersc.no

The core algorithm is in mod_pseudo.F90, which doesn't need to modify.

For operational purpose,
- In ./src/main_pseudo2D.F90,     
    - Set the length of a sequential perturbations for one member to variable i_step.
    - Set domain size to xdim = 1024, ydim = 1024,xy_full = xdim*ydim.  It uses FFT, which generates faster when using power of 2.

- THEN, compile the code by makefile in ./src.

- Configuration of perturbations are set in pseudo2D.nml
- set mod_random_forcing.F90/rdtime as time step of forcing update. Also check the consistency with tcorr in pseudo2D.nml.

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

For a fresh compilation and generating perturbations, where ensemble size is given. The code can run in sequential (slow) or parallel on HPC.
    sbatch run_script.sh


# Q&A
**(a) Regarding the compilation of the fortran libraries (in the makefile)**:
Your compilation line (with mpiifort) includes:

-L/cluster/software/netCDF-Fortran/4.5.2-iimpi-2019b/lib -lfftw3 -lnetcdf -lnetcdff -lhdf5_hl -lhdf5 -lsz -lz -ldl -lm -liomp5 -lpthread -lcurl
It works for me (compilation and run) only if I remove the -lsz option.
—> Do you have any idea if/what this option -lsz is needed for and if it matters if I remove it? (Again, by removing it I am able to compile and run the code apparently successfully).
A: I think you can remove it as you can compile it successfully. I am not good at compilation on unix. It could depends the specific machine.

**(b) input grid size:**
I need to understand better what is the input grid size (xdim=1024, ydim=1024 in fortran routine main_pseudo2D.F90).:
Are we working on the space grid of the input forcing files or in the space grid of the model?
And why 1024x1024? Is that the exact grid size or an approximation to use a power of 2?
A: We want to add perturbation to the forcing grid. Note the grid has be a polar stereograph grid otherwise it needs to interpolate the raw forcing data to a polar stereograph grid. The perturbation grid size is power of 2, said xdim=2^N, because the perturbation code uses FFT algorithm, which runs fastest if the domain size is power of 2. N should be the integer so that 2^(N-1)< size of forcing grid<= 2^N. 

**(c) output grid size:**
Why have the output perturbations become 1-dimension in the output files? (I mean dimension xy in the syncforc.nc files) And how is this 1-dimension ‘xy’ then handled by the model?
A: Using 1-dimension data is easily applied for nextsim data structure and easily saved in this fortran code. The data can export the data in 2-dimension, refering to function save_randfld_synforc() in mod_random_forcing.f90


**(d) usage of the perturbation files:**
More generally, how are these perturbation files then used in a nextsim run? Where could I look at an example code where you have applied those perturbations to the atmospheric forcing? Which GitHub branch should I look at? What part of the code?
A: You can search for code segment within #ifdef ENSEMBLE label or function perturbLoadedData() directly in nextsim/model/externaldata.cpp. For example, https://github.com/nansencenter/nextsim/blob/IOPerturbation-fram-compile/model/externaldata.cpp#L233


**1/** I’d be interested in seeing the kind of script you use to run an ensemble simulation with nextsim (+ an example of associated namelist). In particular, i’m curious to see how you run N instances of the executable nextsim, how you specify the ensemble members, and how you handle the outputs and output names for each member if they are all run in parallele.

The scripts Sukun used was from a complex DA script (his archive), I suggest you start from this simplified script https://github.com/myying/nextsim-dart that I prepared (still work in progress, I haven't added perturbation code to this yet). For each member, we prepare a separate working directory, in which we link the necessary input/boundary files, and then call the nextsim in each directory. Later when the runs finish, we can rename/move the output to a common directory and tagged each with a member id (001, 002...)

job_submit.sh is the one that calls mpirun on HPC

make_nextsim_cfg.sh is the script for namelist template

**2/** As far as i understand, each member seems to be referenced to in the model code by the variable M_ensemble_member, defined from the namelist parameter ensemble_member. Am i right?

In the nextsim option, ensemble_member is just used to locate the correct input file (previously all members are run in the same directory).

With the new script, nextsim just runs a single member in separately, so this option can be removed.


**3/** Also, Sukun, there’s one thing i am not sure to understand at this point: If I use your “off-line-perturbation” code to generate perturbation files offline (which i was able to run already), then how should I specify to the model where are those files, when running the experiment. I understand that it is all in externaldata.cpp and then perturbLoadedData(), and then it seems to me that it is hardcoded in the routine, am I right? Like around line 346 :

std::string filename = Environment::vm()["restart.input_path"].as<std::string>() + "/Perturbations/AtmospherePerturbations_mem" + std::to_string(M_ensemble_member) +"_series" + std::to_string(M_dataset->N_wind+it) + ".nc";
Do you confirm that this is the path where the perturbation files are looked for by the model?

the offline perturbations are stored as synforc nc files

nextsim reads in the perturbation through addPerturbation routine in externaldata.cpp

the file name is hard coded in nextsim and there is no namelist option, you link synforc_$j.nc files to Perturbations/AtmospherePerturbations_mem$n_series$j.nc


**4/** Can you also confirm that if i use the model code in branch IOPerturbation-fram-compile and if i switch-on the compilation key ENSEMBLE, i will activate the perturbation of the forcing without data assimilation ? (i’m a bit confused because of the enkf module…)

The IOPerturbation branch has both perturbation and data assimilation interfaces, but you don't need to worry about the enkf module, for now you only need the offline perturbation module (where nextsim reads in the synforc files),

I'm in the process of removing ENSEMBLE from nextsim code, since the addPerturbation can be done outside, nextsim just reads in nc file through gridoutput.cpp anyway. If you would like, we can work together to clean up the code a bit.