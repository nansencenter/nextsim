/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim: set fenc=utf-8 ft=cpp et sw=4 ts=4 sts=4: */

/**
 * @file   externaldata.hpp
 * @author Sylvain Bouillon <sylvain.bouillon@nersc.no>
 * @date   Tue May  3 09:22:00 2016
 */

#ifndef __ExternalData_H
#define __ExternalData_H 1

#include <environment.hpp>
#include <assert.hpp>
#include <boost/ptr_container/ptr_vector.hpp>

#include <gmshmesh.hpp>
#include <Bamgx.h>
#include <InterpFromGridToMeshx.h>
#include <InterpFromMeshToMesh2dx.h>
#include <ConservativeRemapping.hpp>
#include <netcdf>
#include <dataset.hpp>
#include "debug.hpp"

#if defined ENSEMBLE 
#include "ensemble.hpp"
#endif

#ifdef OASIS
#include <InterpFromMeshToMesh2dx_apply.hpp>
#endif

/**
 * @class ExternalData
 * @brief Manage the external data
 *
 * @see
 *
 */

namespace Nextsim
{

class ExternalData: public DataSet
{

public:

    typedef std::size_t size_type;
    typedef double value_type;
    typedef DataSet Dataset;
    typedef void (* vFunctionCall)(value_type args);

    //constructor for empty object
    ExternalData();

    //normal constructors
    ExternalData(Dataset *dataset, GmshMesh const& mesh, int VariableId, bool is_vector,
            double StartingTime ); // 5 args
    ExternalData(Dataset *dataset, GmshMesh const& mesh, int VariableId, bool is_vector,
            double StartingTime, double SpinUpDuration ); //6 args
    ExternalData(Dataset *dataset, GmshMesh const& mesh, int VariableId,
            bool is_vector, double StartingTime, double SpinUpDuration,
            double bias_correction );// 7 args
    ExternalData(Dataset *dataset, GmshMesh const& mesh, int VariableId,
            bool is_vector, double StartingTime, double SpinUpDuration,
            double bias_correction, int const& ensemble_member );// 8 args

    //constructors for constant forcing
    ExternalData(double ConstantValue );
    ExternalData(double ConstantValue, double ConstantValuebis );
    ExternalData(double ConstantValue,
        double StartingTime, double SpinUpDuration );
    ExternalData(double ConstantValue, double ConstantValuebis,
        double StartingTime, double SpinUpDuration );

    ~ExternalData();

    void check_and_reload(std::vector<double> const& RX,
            std::vector<double> const& RY, const double current_time, const int ensemble_member );

    value_type operator[] (const size_type i);
    value_type get(const size_type i);
    size_type size();

    std::vector<double> getVector();
    void clear();
    void loadDataset(Dataset *dataset, std::vector<double> const& RX,
            std::vector<double> const& RY );

#if defined OASIS
    void check_and_reload(std::vector<double> const& RX_in,
            std::vector<double> const& RY_in, const double current_time,
            Communicator comm, const int cpl_time, const int cpl_dt);
#endif

    void transformData(Dataset *dataset);
    void interpolateDataset(Dataset *dataset, std::vector<double> const& RX,
            std::vector<double> const& RY );

#if defined OASIS
    void receiveCouplingData(Dataset *dataset, int cpl_time, Communicator comm);
#endif

    bool isInitialized() const { return M_initialized; }
    bool isVector() const { return M_is_vector; }

    std::string getDatasetName() const {
        if(M_is_constant)
            return "constant";
        return M_dataset->name;
    }
    std::string getVariableName() const {
        if(M_is_constant)
            return "constant";
        return M_dataset->variables[M_VariableId].name;
    }
#ifdef ENSEMBLE    
    Dataset* get_M_dataset(){
        return M_dataset;
    }
#endif
 
private:
    double fdt;
    std::vector<double> fcoeff;
    Dataset *M_dataset;
    bool M_is_constant;
    bool M_initialized;
    bool M_is_vector;
    int M_target_size;

    int M_VariableId;
    double M_current_time;
    double M_constant_value;
    double M_constant_valuebis;
    double M_bias_correction;
    int M_VariableIdbis;
    double M_StartingTime;
    double M_SpinUpDuration;
    double M_factor;
    int M_ensemble_member;

    LogLevel M_log_level;
    bool M_log_all;
    Communicator M_comm;
};

} // Nextsim
#endif // __ExternalData_H
