//
// Created by Einar Ólason on 15/09/2022.
//

#ifndef NEXTSIM_SCATTERGATHER_HPP
#define NEXTSIM_SCATTERGATHER_HPP

#include <vector>
#include <boost/mpi.hpp>

#include "variablemanager.hpp"
#include "datasetmanager.hpp"
#include "meshpartition.hpp"
#include "debug.hpp"

namespace Nextsim {

    class ScatterGather : public VariableManager, public DataSetManager
    {
    public:
        ScatterGather(Timer::timer *timer)
                : DataSetManager(timer), vm(Environment::vm())
        {}

    protected:
        void gatherFieldsElement(std::vector<double>& interp_in_elements);
        void scatterFieldsElement(double* interp_elt_out);

        //void gatherUM(std::vector<double>& um);
        void gatherNodalField(std::vector<double> const& field_local, std::vector<double>& field_root);
        void scatterNodalField(std::vector<double> const& field_root, std::vector<double>& field_local);

        // other interfaces
        void gatherNodalField(std::vector<double> const& field1_local, std::vector<double> const& field2_local,
                              std::vector<double>& field1_root, std::vector<double>& field2_root);

        void gatherElementField(std::vector<double> const& field_local, std::vector<double>& field_root, int nb_fields = 1);
        void scatterElementField(std::vector<double> const& field_root, std::vector<double>& field_local, int nb_fields = 1);

        void gatherFieldsNode(std::vector<double>& interp_in_elements, std::vector<int> const& rmap_nodes, std::vector<int> sizes_nodes);
        void scatterFieldsNode(double* interp_nd_out);

        void partitionMeshRestart();
        void collectNodesRestart(std::vector<double>& interp_nd_out);
        void collectElementsRestart(std::vector<double>& interp_elt_out,
                                    std::vector<std::vector<double>*> &data_elements_root);

        void collectVariables(std::vector<double>& interp_elt_in_local, bool ghosts);
        void redistributeVariables(std::vector<double> const& out_elt_values, bool const& apply_maxima);

        // IO
        void redistributeVariablesIO(std::vector<double> const& out_elt_values,
                                     std::vector<ModelVariable*> &vars_elements);
        void scatterFieldsElementIO(std::vector<double> const& interp_elt_out,
                                    std::vector<ModelVariable*> &vars_elements);

        void collectVariablesIO(std::vector<double>& elt_values_local,
                                std::vector<ModelVariable*> const& vars_elements,
                                std::vector<ExternalData*> const& ext_data_elements,
                                bool const& ghosts);

        void gatherFieldsElementIO(std::vector<double>& elt_values_root,
                                   std::vector<ModelVariable*> const& vars_elements,
                                   std::vector<ExternalData*> const& ext_data_elements);
        void gatherFieldsElementIO(std::vector<double>& elt_values_root,
                                   std::vector<ModelVariable*> const& vars_elements)
        {
            std::vector<ExternalData*> ext_data_elements = {};// add a place-holder
            this->gatherFieldsElementIO(elt_values_root, vars_elements, ext_data_elements);
        }

        // Variables
        int M_nb_var_node;

        int M_prv_local_ndof;
        int M_prv_num_nodes;
        int M_prv_num_elements;
        int M_prv_global_num_nodes;
        int M_prv_global_num_elements;

    private:
        po::variables_map vm;
        boost::mpi::timer chrono, chrono_tot;

    };

}

#endif //NEXTSIM_SCATTERGATHER_HPP
