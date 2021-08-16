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
//	void p_pseudo2D_fld_sub(int const *xdm, int const *ydm, int const *previous_perturbation_exist );   
  void p_pseudo2D_fld_sub(int const *xdm, int const *ydm, double *synforc, double *randfld,int const *previous_perturbation_exist);
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

class ensemble
{
    // public:

    //     struct synoptic_forcing {

    //         struct zonal_wind {
    //             std::vector<double> data;
    //             const int id  = 2;
    //             static constexpr const char* name  = "uwind";
    //         } uwind;


    //         struct merid_wind {
    //             std::vector<double> data;
    //             const int id  = 3;
    //             static constexpr const char* name  = "vwind";
    //         } vwind;

    //         struct temperature {
    //             std::vector<double> data;
    //             const int id  = 4;
    //             static constexpr const char* name  = "t2air";
    //         } t2air;

    //         struct pressure {
    //             std::vector<double> data;
    //             const int id  = 5;
    //             static constexpr const char* name  = "slp";
    //         } slp;

    //         struct precipitation {
    //             std::vector<double> data;
    //             const int id  = 6;
    //             static constexpr const char* name  = "precip";
    //         } precip;

    //         struct humidity {
    //             std::vector<double> data;
    //             const int id  = 7;
    //             static constexpr const char* name  = "relhum";
    //         } relhum;

    //     } synoptic;

    // private:

    //     std::vector<std::string> M_ranfile;
    //     std::string M_ranpath;

    public:

        void synopticPerturbation(std::vector<double> &synforc,std::vector<double> &randfld,int const& ydim, int const& xdim, int const& previous_perturbation_exist); 

        void addPerturbation(std::vector<double>& perturbed_field, std::vector<double>& synforc, int M_full, int N_full, int x_start, int y_start, int x_count, int y_count, int opr);
        // void addPerturbation(std::vector<double>& perturbed_field, std::vector<double>& synforc, int M_full, int N_full, int x_start, int y_start, int x_count, int y_count,int index, int opr);


        // void synopticPerturbation();

        // void synopticPerturbation(int);

        // void synopticPerturbation(int const& ydim, int const& xdim,std::vector<std::vector<double>>& synforc,std::vector<std::vector<double>>& randfld, int const& previous_perturbation_exist);

        // void synopticPerturbation(int const& ydim, int const& xdim, std::vector<std::vector<double> > &synforc,std::vector<std::vector<double> > &randfld, int const& previous_perturbation_exist, double *synforc_p);


        // void addPerturbation();

        // void addPerturbation(int);

        // void addPerturbation(std::vector<double>&, std::vector<double>&, int, int);

        // void addPerturbation(std::vector<double>&, std::vector<double>&, std::vector<double>& synforc, int M_full, int N_full, int, int, int, int);

        //void addPerturbation(Dataset *dataset, std::vector<std::vector<double> >&, std::vector<std::vector<double> >&);

        // void computeMinMax(const std::vector<double> &, const char*);

        // void computeVecMean(const std::vector<double> &, const char*);

        // void getpath(std::string);
//
       // void loadPerturbation(double *synforc,int rdm, int ranid);
       // void loadPerturbation(std::vector<std::vector<float> >& ,int, int);
       // std::vector<std::vector<double> > loadPerturbation(int, int);
};
