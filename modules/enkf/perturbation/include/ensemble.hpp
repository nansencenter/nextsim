/*
Header file for ensemble.cpp
@author= -aLi- <ali.aydogdu@nersc.no>
- This class calls fortran objects in order to create synoptic perturbations to generate an ensemble using neXtSIM
- It reads the output file from ran_update and ports it to related forcing data in neXtSIM 
*/

#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <numeric>                                                              
#include <vector>
#include <fstream>
#include <string>
#include <boost/program_options.hpp>

extern "C" {
	double p_pseudo2D_fld_sub();
/* Other subroutines in libpseudo2D.dylib fortran library
	void init_fvars_();
        void m_set_random_seed_(); 
	void memory_stack_(); 
	void mod_forcing_nersc_(); 
	void mod_pseudo_(); 
	void mod_random_forcing_(); 
	void assign_force_(); 
	void assign_vars_(); 
	void var_sqrt_(); 
*/
}

using namespace std;

namespace Nextsim 
{
class ensemble
{
public:

	struct synoptic_forcing { 

		struct zonal_wind {
		std::vector<double> data; 
		const int id  = 2; 
		static constexpr const char* name  = "uwind"; 
		} uwind;


		struct merid_wind {
		std::vector<double> data; 
		const int id  = 3; 
		static constexpr const char* name  = "vwind"; 
		} vwind;

		struct temperature {
		std::vector<double> data; 
		const int id  = 4; 
		static constexpr const char* name  = "t2air"; 
		} t2air;

		struct pressure {
		std::vector<double> data; 
		const int id  = 5; 
		static constexpr const char* name  = "slp"; 
		} slp;

		struct precipitation {
		std::vector<double> data; 
		const int id  = 6; 
		static constexpr const char* name  = "precip"; 
		} precip;

		struct humidity {
		std::vector<double> data; 
		const int id  = 7; 
		static constexpr const char* name  = "relhum"; 
		} relhum;

	} synoptic;

private: 
	

	std::vector<std::string> ranfile = { "/docker_io/synforc.00", "/docker_io/synforc.01" };

public:

	void synopticPerturbation();

	void synopticPerturbation(int);

	void addPerturbation();

	void addPerturbation(int);

	void addPerturbation(std::vector<double>&, std::vector<double>&, int, int);
	
	void computeMinMax(const std::vector<double> &, const char*); 

	void computeVecMean(const std::vector<double> &, const char*);

};
};
