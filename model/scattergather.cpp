//
// Created by Einar Ólason on 15/09/2022.
//

#include "scattergather.hpp"

namespace Nextsim
{

//------------------------------------------------------------------------------------------------------
//! Gathers field values (velocities, displacements) at the nodes into a single structure, interp_node_in_local.
//! Called by the interpFields() function.
    void
    ScatterGather::gatherFieldsNode(std::vector<double>& interp_in_nodes, std::vector<int> const& rmap_nodes, std::vector<int> sizes_nodes)
    {
        chrono.restart();

        LOG(DEBUG) <<"----------GATHER NODE starts\n";

        M_nb_var_node = 6;
        std::vector<double> interp_node_in_local(M_nb_var_node*M_prv_local_ndof,0.);

        chrono.restart();
        //std::cout<<"Nodal Interp starts\n";
        //std::cout<<"NODAL: Interp starts\n";

        for (int i=0; i<M_prv_local_ndof; ++i)
        {
            int tmp_nb_var = 0;

            // VT
            interp_node_in_local[M_nb_var_node*i] = M_VT[i];
            tmp_nb_var++;
            interp_node_in_local[M_nb_var_node*i+tmp_nb_var] = M_VT[i+M_prv_num_nodes];
            tmp_nb_var++;

            // UM
            interp_node_in_local[M_nb_var_node*i+tmp_nb_var] = M_UM[i];
            tmp_nb_var++;
            interp_node_in_local[M_nb_var_node*i+tmp_nb_var] = M_UM[i+M_prv_num_nodes];
            tmp_nb_var++;

            // UT
            interp_node_in_local[M_nb_var_node*i+tmp_nb_var] = M_UT[i];
            tmp_nb_var++;
            interp_node_in_local[M_nb_var_node*i+tmp_nb_var] = M_UT[i+M_prv_num_nodes];
            tmp_nb_var++;

            if ( tmp_nb_var != M_nb_var_node )
                throw std::logic_error("tmp_nb_var not equal to nb_var");
        }

        std::for_each(sizes_nodes.begin(), sizes_nodes.end(), [&](int& f){ f = M_nb_var_node*f; });

        if (M_comm.rank() == 0)
        {
            interp_in_nodes.resize(M_nb_var_node*M_prv_global_num_nodes);
            boost::mpi::gatherv(M_comm, interp_node_in_local, &interp_in_nodes[0], sizes_nodes, 0);
        }
        else
        {
            boost::mpi::gatherv(M_comm, interp_node_in_local, 0);
        }

        if (M_comm.rank() == 0)
        {
            auto interp_in_nodes_nrd = interp_in_nodes;

            for (int i=0; i<M_prv_global_num_nodes; ++i)
            {
                int ri =  rmap_nodes[i];

                for (int j=0; j<M_nb_var_node; ++j)
                {
                    interp_in_nodes[M_nb_var_node*i+j] = interp_in_nodes_nrd[M_nb_var_node*ri+j];
                }
            }
        }

        LOG(DEBUG) <<"----------GATHER NODE done in "<< chrono.elapsed() <<"s\n";
    }//gatherFieldsNode


//------------------------------------------------------------------------------------------------------
//! Scatters field values (velocities, displacements) at the field nodes from the root
//! Called by the interpFields() and readRestart() functions.
    void
    ScatterGather::scatterFieldsNode(double* interp_nd_out)
    {
        chrono.restart();

        LOG(DEBUG) <<"----------SCATTER NODE starts\n";

        std::vector<double> in_nd_values;

        if (M_comm.rank() == 0)
        {
            in_nd_values.resize(M_nb_var_node*M_id_nodes.size());

            for (int i=0; i<M_id_nodes.size(); ++i)
            {
                //int ri = rmap_nodes.right.find(id_nodes[i])->second-1;
                int ri = M_id_nodes[i]-1;

                for (int j=0; j<M_nb_var_node; ++j)
                {
                    in_nd_values[M_nb_var_node*i+j] = interp_nd_out[M_nb_var_node*ri+j];
                }
            }
        }

        std::vector<double> out_nd_values(M_nb_var_node*M_num_nodes);
        std::vector<int> sizes_nodes = M_sizes_nodes_with_ghost;

        if (M_comm.rank() == 0)
        {
            std::for_each(sizes_nodes.begin(), sizes_nodes.end(), [&](int& f){ f = M_nb_var_node*f; });
            boost::mpi::scatterv(M_comm, in_nd_values, sizes_nodes, &out_nd_values[0], 0);
        }
        else
        {
            boost::mpi::scatterv(M_comm, &out_nd_values[0], M_nb_var_node*M_num_nodes, 0);
        }


        M_VT.resize(2*M_num_nodes);
        M_UM.resize(2*M_num_nodes);
        M_UT.resize(2*M_num_nodes);

        D_tau_w.assign(2*M_num_nodes,0.);
        D_tau_a.assign(2*M_num_nodes,0.);

        for (int i=0; i<M_num_nodes; ++i)
        {
            // VT
            M_VT[i] = out_nd_values[M_nb_var_node*i];
            M_VT[i+M_num_nodes] = out_nd_values[M_nb_var_node*i+1];

            // UM
            M_UM[i] = out_nd_values[M_nb_var_node*i+2];
            M_UM[i+M_num_nodes] = out_nd_values[M_nb_var_node*i+3];

            // UT
            M_UT[i] = out_nd_values[M_nb_var_node*i+4];
            M_UT[i+M_num_nodes] = out_nd_values[M_nb_var_node*i+5];
        }


        LOG(DEBUG) <<"----------SCATTER NODE done in "<< chrono.elapsed() <<"s\n";
    }//scatterFieldsNode


//------------------------------------------------------------------------------------------------------
//! Sends displacement vector to the root process.
//! Called by the regrid(), checkUpdateDrifters(), exportResults() functions.
    void
    ScatterGather::gatherNodalField(std::vector<double> const& field_local, std::vector<double>& field_root)
    {
        int const num_components = field_local.size() / M_num_nodes;
        std::vector<double> um_local(num_components * M_local_ndof, 0.);
        for (int i=0; i<M_local_ndof; ++i)
            for (int k=0; k<num_components; ++k)
                um_local[num_components * i + k] = field_local[i + k*M_num_nodes];

        std::vector<int> sizes_nodes = M_sizes_nodes;
        std::for_each(sizes_nodes.begin(), sizes_nodes.end(),
                      [&](int& f){ f = num_components * f; });

        // send displacement vector to the root process (rank 0)
        if (M_comm.rank() == 0)
        {
            field_root.resize(num_components * M_ndof);
            boost::mpi::gatherv(M_comm, um_local, &field_root[0], sizes_nodes, 0);
        }
        else
        {
            boost::mpi::gatherv(M_comm, um_local, 0);
        }

        if (M_comm.rank() == 0)
        {
            int global_num_nodes = M_mesh.numGlobalNodes();

            auto field_root_nrd = field_root;

            for (int i=0; i<global_num_nodes; ++i)
            {
                int ri =  M_rmap_nodes[i];
                for(int k=0; k<num_components; k++)
                    field_root[i + k * global_num_nodes]
                            = field_root_nrd[num_components * ri + k];
            }
        }
    }//gatherNodalField

//------------------------------------------------------------------------------------------------------
//! Gathers nodal fields.
//! Only called by the advecRoot() function, which is not used anymore.
#if 0
    void
ScatterGather::gatherNodalField(std::vector<double> const& field1_local, std::vector<double> const& field2_local,
                                std::vector<double>& field1_root, std::vector<double>& field2_root)
{
    std::vector<double> um_local(4*M_local_ndof,0.);
    for (int i=0; i<M_local_ndof; ++i)
    {
        // field1
        um_local[4*i] = field1_local[i];
        um_local[4*i+1] = field1_local[i+M_num_nodes];

        // field2
        um_local[4*i+2] = field2_local[i];
        um_local[4*i+3] = field2_local[i+M_num_nodes];
    }

    std::vector<int> sizes_nodes = M_sizes_nodes;
    std::for_each(sizes_nodes.begin(), sizes_nodes.end(), [&](int& f){ f = 4*f; });

    // send displacement vector to the root process (rank 0)

    std::vector<double> um_root;

    if (M_comm.rank() == 0)
    {
        um_root.resize(4*M_ndof);
        boost::mpi::gatherv(M_comm, um_local, &um_root[0], sizes_nodes, 0);
    }
    else
    {
        boost::mpi::gatherv(M_comm, um_local, 0);
    }

    if (M_comm.rank() == 0)
    {
        int global_num_nodes = M_mesh.numGlobalNodes();

        field1_root.resize(2*M_ndof);
        field2_root.resize(2*M_ndof);

        auto um_root_nrd = um_root;

        for (int i=0; i<global_num_nodes; ++i)
        {
            int ri =  M_rmap_nodes[i];

            field1_root[i] = um_root_nrd[4*ri];
            field1_root[i+global_num_nodes] = um_root_nrd[4*ri+1];

            field2_root[i] = um_root_nrd[4*ri+2];
            field2_root[i+global_num_nodes] = um_root_nrd[4*ri+3];
        }
    }
}//gatherNodalField
#endif


//------------------------------------------------------------------------------------------------------
//! Scatter nodal fields.
//! Called by the advect() function.
    void
    ScatterGather::scatterNodalField(std::vector<double> const& field_root, std::vector<double>& field_local)
    {
        std::vector<double> in_nd_values;

        if (M_comm.rank() == 0)
        {
            int global_num_nodes = M_mesh.numGlobalNodes();

            in_nd_values.resize(2*M_id_nodes.size());

            for (int i=0; i<M_id_nodes.size(); ++i)
            {
                int ri = M_id_nodes[i]-1;

                in_nd_values[2*i]   = field_root[ri];
                in_nd_values[2*i+1] = field_root[ri+global_num_nodes];

                // for (int j=0; j<2; ++j)
                // {
                //     in_nd_values[2*i+j] = field_root[2*ri+j];
                // }
            }
        }

        field_local.resize(2*M_num_nodes);
        std::vector<int> sizes_nodes = M_sizes_nodes_with_ghost;

        if (M_comm.rank() == 0)
        {
            std::for_each(sizes_nodes.begin(), sizes_nodes.end(), [&](int& f){ f = 2*f; });
            boost::mpi::scatterv(M_comm, in_nd_values, sizes_nodes, &field_local[0], 0);
        }
        else
        {
            boost::mpi::scatterv(M_comm, &field_local[0], 2*M_num_nodes, 0);
        }

        std::vector<double> field_local_copy = field_local;

        for (int i=0; i<M_num_nodes; ++i)
        {
            // U component
            field_local[i] = field_local_copy[2*i];
            // V component
            field_local[i+M_num_nodes] = field_local_copy[2*i+1];
        }
    }//scatterNodalField


//------------------------------------------------------------------------------------------------------
//! Gather field values over elements.
//! Called by the advect(), diffuse(), checkUpdateDrifters() functions.
    void
    ScatterGather::gatherElementField(std::vector<double> const& field_local, std::vector<double>& field_root, int nb_fields)
    {
        std::vector<double> field_local_copy(nb_fields*M_local_nelements);

        for (int i=0; i<M_local_nelements; ++i)
        {
            // copy values without ghosts
            for (int j=0; j<nb_fields; ++j)
            {
                field_local_copy[nb_fields*i+j] = field_local[nb_fields*i+j];
            }
        }

        std::vector<int> sizes_elements = M_sizes_elements;

        if (nb_fields != 1)
        {
            std::for_each(sizes_elements.begin(), sizes_elements.end(), [&](int& f){ f = nb_fields*f; });
        }

        if (M_comm.rank() == 0)
        {
            field_root.resize(nb_fields*M_mesh_root.numTriangles());
            boost::mpi::gatherv(M_comm, field_local_copy, &field_root[0], sizes_elements, 0);
        }
        else
        {
            boost::mpi::gatherv(M_comm, field_local_copy, 0);
        }

        if (M_comm.rank() == 0)
        {
            auto field_root_nrd = field_root;

            for (int i=0; i<M_mesh_root.numTriangles(); ++i)
            {
                int ri = M_rmap_elements[i];

                for (int j=0; j<nb_fields; ++j)
                {
                    field_root[nb_fields*i+j] = field_root_nrd[nb_fields*ri+j];
                }
            }
        }

    }//gatherElementField


//------------------------------------------------------------------------------------------------------
//! Scatters back vector of field values at the elements from root to all processes.
//! Called by the advectRoot() function.
    void
    ScatterGather::scatterElementField(std::vector<double> const& field_root, std::vector<double>& field_local, int nb_fields)
    {
        std::vector<double> field_root_extended;

        if (M_comm.rank() == 0)
        {
            field_root_extended.resize(nb_fields*M_id_elements.size());

            for (int i=0; i<M_id_elements.size(); ++i)
            {
                int ri = M_id_elements[i]-1;

                for (int j=0; j<nb_fields; ++j)
                {
                    field_root_extended[nb_fields*i+j] = field_root[nb_fields*ri+j];
                }
            }
        }

        field_local.resize(nb_fields*M_num_elements);

        std::vector<int> sizes_elements = M_sizes_elements_with_ghost;
        if (nb_fields != 1)
        {
            std::for_each(sizes_elements.begin(), sizes_elements.end(), [&](int& f){ f = nb_fields*f; });
        }

        if (M_comm.rank() == 0)
        {
            boost::mpi::scatterv(M_comm, field_root_extended, sizes_elements, &field_local[0], 0);
        }
        else
        {
            boost::mpi::scatterv(M_comm, &field_local[0], nb_fields*M_num_elements, 0);
        }
    }//scatterElementField

//! collect the restart elemental variables (already on the root)
//! and put them into 1 long vector to be scattered
//! called by readRestart()
    void
    ScatterGather::collectElementsRestart(std::vector<double>& interp_elt_out,
                                          std::vector<std::vector<double>*> &data_elements_root)
{

    // get the variables (only on the root processor so far)
    // from data and put it in interp_elt_out
    // TODO do something similar for the nodes
    std::vector<double> out_elt_values;

    int num_elements_root = M_id_elements.size();
    int const nb_var_element = data_elements_root.size();
    interp_elt_out.resize(nb_var_element*num_elements_root);

    for (int i=0; i<num_elements_root; ++i)
{
    int ri = M_id_elements[i]-1;
    for(int j=0; j<data_elements_root.size(); j++)
{
    auto ptr = data_elements_root[j];
    interp_elt_out[nb_var_element*i+j] = (*ptr)[ri];
}
}

}//collectElementsRestart


//------------------------------------------------------------------------------------------------------
//! Gets the variables (only on the root processor so far) from data and store it in a structure (interp_elt_out)
//! Called by the readRestart() function.
void
ScatterGather::collectNodesRestart(std::vector<double>& interp_nd_out)
{
    // * output: interp_nd_out is vector containing all the variables
    //   on the nodes to be scattered from root during readRestart

    M_nb_var_node = 6;
    if (M_comm.rank() == 0)
    {
        int num_nodes_root = M_mesh_root.numNodes();
        interp_nd_out.resize(M_nb_var_node*num_nodes_root,0.);

        int tmp_nb_var = 0;

        for (int i=0; i<num_nodes_root; ++i)
        {
            tmp_nb_var = 0;

            // VT
            interp_nd_out[M_nb_var_node*i+tmp_nb_var] = M_VT[i];
            tmp_nb_var++;

            interp_nd_out[M_nb_var_node*i+tmp_nb_var] = M_VT[i+num_nodes_root];
            tmp_nb_var++;

            // UM
            interp_nd_out[M_nb_var_node*i+tmp_nb_var] = M_UM[i];
            tmp_nb_var++;

            interp_nd_out[M_nb_var_node*i+tmp_nb_var] = M_UM[i+num_nodes_root];
            tmp_nb_var++;

            // UT
            interp_nd_out[M_nb_var_node*i+tmp_nb_var] = M_UT[i];
            tmp_nb_var++;

            interp_nd_out[M_nb_var_node*i+tmp_nb_var] = M_UT[i+num_nodes_root];
            tmp_nb_var++;

            if(tmp_nb_var>M_nb_var_node)
            {
                throw std::logic_error("tmp_nb_var not equal to nb_var");
            }
        }
    }
}//collectRootRestart

//------------------------------------------------------------------------------------------------------
//! Gathers information about the fields for interpolation onto the mesh grid.
//! Called by interpFields() function.
    void
    ScatterGather::gatherFieldsElement(std::vector<double>& interp_in_elements)
    {
        int nb_var_element = M_prognostic_variables_elt.size();

        chrono.restart();
        LOG(DEBUG) <<"----------GATHER ELEMENT starts\n";

        std::vector<int> sizes_elements = M_sizes_elements;
        std::for_each(sizes_elements.begin(), sizes_elements.end(), [&](int& f){ f = nb_var_element*f; });

        std::vector<double> interp_elt_in_local;
        bool ghosts = false;
        this->collectVariables(interp_elt_in_local, ghosts);

        if (M_comm.rank() == 0)
        {
            interp_in_elements.resize(nb_var_element*M_mesh_previous_root.numTriangles());
            boost::mpi::gatherv(M_comm, interp_elt_in_local, &interp_in_elements[0], sizes_elements, 0);
        }
        else
            boost::mpi::gatherv(M_comm, interp_elt_in_local, 0);

        if (M_comm.rank() == 0)
        {
            auto interp_in_elements_nrd = interp_in_elements;
            for (int i=0; i<M_mesh_previous_root.numTriangles(); ++i)
            {
                int ri = M_rmap_elements[i];
                for (int j=0; j<nb_var_element; ++j)
                    interp_in_elements[nb_var_element*i+j] = interp_in_elements_nrd[nb_var_element*ri+j];
            }
        }

        LOG(DEBUG) <<"----------GATHER ELEMENT done in "<< chrono.elapsed() <<"s\n";
    }//gatherFieldsElement


//------------------------------------------------------------------------------------------------------
//! Gathers information about the fields for outputting.
//! Called by the writeRestart() and exportResults() function.
    void
    ScatterGather::gatherFieldsElementIO( std::vector<double>& elt_values_root,
                                          std::vector<ModelVariable*> const& vars_elements,
                                          std::vector<ExternalData*> const& ext_data_elements)
    {

        chrono.restart();
        LOG(DEBUG) <<"----------IO: GATHER ELEMENT starts\n";

        int const nb_var_element = vars_elements.size() + ext_data_elements.size();
        std::vector<double> elt_values_local;
        bool const ghosts = false;
        this->collectVariablesIO(elt_values_local, vars_elements,
                                 ext_data_elements, ghosts);

        std::vector<int> sizes_elements = M_sizes_elements;
        std::for_each(sizes_elements.begin(), sizes_elements.end(), [&](int& f){ f = nb_var_element*f; });

        if (M_comm.rank() == 0)
        {
            elt_values_root.resize(nb_var_element*M_mesh_root.numTriangles());
            boost::mpi::gatherv(M_comm, elt_values_local, &elt_values_root[0], sizes_elements, 0);
        }
        else
        {
            boost::mpi::gatherv(M_comm, elt_values_local, 0);
        }

        LOG(DEBUG) <<"----------IO: GATHER ELEMENT done in "<< chrono.elapsed() <<"s\n";
    }//gatherFieldsElementsIO


//------------------------------------------------------------------------------------------------------
//! Scatters (redistributes) P0 (elemental) field values to the subdomains (parallel computing).
//! Called by the interpFields() function.
    void
    ScatterGather::scatterFieldsElement(double* interp_elt_out)
    {
        chrono.restart();
        LOG(DEBUG) <<"----------SCATTER ELEMENT starts\n";

        int nb_var_element = M_prognostic_variables_elt.size();
        std::vector<int> sizes_elements = M_sizes_elements_with_ghost;
        std::vector<double> in_elt_values;

        if (M_comm.rank() == 0)
        {
            in_elt_values.resize(nb_var_element*M_id_elements.size());
            for (int i=0; i<M_id_elements.size(); ++i)
            {
                int ri = M_id_elements[i]-1;
                for (int j=0; j<nb_var_element; ++j)
                    in_elt_values[nb_var_element*i+j] = interp_elt_out[nb_var_element*ri+j];
            }
        }

        std::vector<double> out_elt_values(nb_var_element*M_num_elements);
        if (M_comm.rank() == 0)
        {
            std::for_each(sizes_elements.begin(), sizes_elements.end(), [&](int& f){ f = nb_var_element*f; });
            boost::mpi::scatterv(M_comm, in_elt_values, sizes_elements, &out_elt_values[0], 0);
        }
        else
            boost::mpi::scatterv(M_comm, &out_elt_values[0], nb_var_element*M_num_elements, 0);

        for(auto ptr: M_variables_elt)
        {
            if(ptr->isPrognostic())
                // resize prognostic variables
                // - they are set in redistributeVariables
                ptr->resize(M_num_elements);
            else
                // assign diagnostic variables
                // - they are not interpolated
                ptr->assign(M_num_elements, 0.);
        }
        this->redistributeVariables(out_elt_values, true);//apply maxima during interpolation

        LOG(DEBUG) <<"----------SCATTER ELEMENT done in "<< chrono.elapsed() <<"s\n";
    }//scatterFieldsElement


//------------------------------------------------------------------------------------------------------
//! scatter from root to local
//! * both input and output are 1 long vector containing all the variables
//! Called by the restartScatterElementVariables() function.
    void
    ScatterGather::scatterFieldsElementIO(std::vector<double> const& elt_values_root,
                                          std::vector<ModelVariable*> &vars_elements)
    {
        //! * elt_values_root is a vector containing all the variables to be
        //!   redistributed (eg after scattering from root) into the
        //!   individual variables (eg M_conc, M_thick,...)
        //!   - rearranged using M_id_elements and passed to
        //!     boost::mpi::scatterv
        //! * data is a vector of pointers to the variables to be assigned
        //!   values from elt_values_local
        //!   - passed to redistributeVariablesIO
        chrono.restart();

        LOG(DEBUG) <<"----------SCATTER ELEMENT starts\n";

        std::vector<int> sizes_elements = M_sizes_elements_with_ghost;
        int const nb_var_element = vars_elements.size();

        std::vector<double> elt_values_local(nb_var_element*M_num_elements);
        if (M_comm.rank() == 0)
        {
            std::for_each(sizes_elements.begin(), sizes_elements.end(),
                          [&](int& f){ f = nb_var_element*f; });
            boost::mpi::scatterv(M_comm, elt_values_root, sizes_elements,
                                 &elt_values_local[0], 0);
        }
        else
        {
            boost::mpi::scatterv(M_comm, &elt_values_local[0],
                                 nb_var_element*M_num_elements, 0);
        }

        // transfer data from elt_values_local to vars_elements
        this->redistributeVariablesIO(elt_values_local, vars_elements);

        LOG(DEBUG) <<"----------SCATTER ELEMENT done in "<< chrono.elapsed() <<"s\n";
    }//scatterFieldsElementIO

//------------------------------------------------------------------------------------------------------
//! Collects model variables and stores them into a single vector, interp_elt_in_local: called by the update() function,
//! before updating all variables after solving.
//! Called by the gatherFieldsElement() function.
    void
    ScatterGather::collectVariables(std::vector<double>& interp_elt_in_local, bool ghosts)
    {
        // ELEMENT INTERPOLATION With Cavities
        int nb_var_element = M_prognostic_variables_elt.size();
        int num_elements = M_local_nelements;
        if (ghosts)
            num_elements = M_num_elements;

        //! loop over elements and grab the different variables to be interpolated
        interp_elt_in_local.resize(nb_var_element*num_elements);
        for (int i=0; i<num_elements; ++i)
        {
            for (int j = 0; j<M_prognostic_variables_elt.size(); j++)
            {
                auto vptr = M_prognostic_variables_elt[j];
                double val = (*vptr)[i];
                switch (vptr->getInterpTransformation())
                {
                    case ModelVariable::interpTransformation::conc:
                        val *= M_conc[i];
                        break;
                    case ModelVariable::interpTransformation::thick:
                        val *= M_thick[i];
                        break;
                    case ModelVariable::interpTransformation::enthalpy:
                        val = ( val - physical::mu*physical::si*physical::Lf/(physical::C*val) ) * M_thick[i]; // (Winton, 2000, eq 39) times volume with f1=1
                        break;
                }
                interp_elt_in_local[nb_var_element*i+j] = val;
            }
        }
    }//collectVariables


//------------------------------------------------------------------------------------------------------
//! Collects model variables and stores them into a single vector, interp_elt_in_local, for outputting.
//! Called by the gatherFieldsElementIO() function.
    void
    ScatterGather::collectVariablesIO(std::vector<double>& elt_values_local,
                                      std::vector<ModelVariable*> const& vars_elements,
                                      std::vector<ExternalData*> const& ext_data_elements,
                                      bool const& ghosts)
    {

        int const nb_data = vars_elements.size();
        int const nb_ext_data = ext_data_elements.size();
        int const nb_var_element = nb_data + nb_ext_data;

        int num_elements = M_local_nelements;
        if (ghosts)
            num_elements = M_num_elements;
        elt_values_local.resize(nb_var_element*num_elements);

        for (int i=0; i<num_elements; ++i)
        {
            int k = 0;
            for(int j=0; j<nb_data; j++, k++)
            {
                auto ptr = vars_elements[j];
                elt_values_local[nb_var_element*i+k] = (*ptr)[i];
            }
            for(int j=0; j<nb_ext_data; j++, k++)
            {
                auto ptr = ext_data_elements[j];
                elt_values_local[nb_var_element*i+k] = (*ptr)[i];
            }
        }
    }//collectVariablesIO


//------------------------------------------------------------------------------------------------------
//! Interpolates variables (other than velocities and displacements) onto the mesh grid once updated.
//! \note apply_maxima can be true for interpolation, but should be false if calling for advection
//! - then we redistribute extra ice concentration with a ridging scheme.
//! Called by the interpFields() function, after the advect() function.
    void
    ScatterGather::redistributeVariables(std::vector<double> const& out_elt_values, bool const& apply_maxima)
    {
        //! -1) check if we are using the maximum values for the variables
        //!     - don't do this if calling during advection
        int nb_var_element = M_prognostic_variables_elt.size();
        std::vector<bool> has_max(nb_var_element, false);
        if(apply_maxima)
            for(int j=0; j<nb_var_element; j++)
                has_max[j] = M_prognostic_variables_elt[j]->hasMaxVal();

        //! -2) loop over elements and assign the values to the different variables
        for (int i=0; i<M_num_elements; ++i)
        {
            for (int j = 0; j<M_prognostic_variables_elt.size(); j++)
            {
                auto vptr = M_prognostic_variables_elt[j];
                if(M_comm.rank() + i==0)
                    LOG(DEBUG)<<"redistribute (none): variable "<<j << " = "<<vptr->name()<<"\n";

                double val = out_elt_values[nb_var_element*i+j];
                if(vptr->getInterpTransformation() == ModelVariable::interpTransformation::conc)
                {
                    if(M_conc[i]>0)
                        val /= M_conc[i];
                    else if (vptr->hasValueNoThickIce())
                        val = vptr->valueNoThickIce();
                }
                else if (vptr->getInterpTransformation() == ModelVariable::interpTransformation::thick)
                {
                    if(M_thick[i]>0)
                        val /= M_thick[i];
                    else if (vptr->hasValueNoThickIce())
                        val = vptr->valueNoThickIce();
                }
                else if (vptr->getInterpTransformation() == ModelVariable::interpTransformation::enthalpy)
                {
                    if(M_thick[i]>0)
                    {
                        double enth = out_elt_values[nb_var_element*i+j]/M_thick[i];//divide by volume to get enthalpy back
                        val = 0.5*(
                                enth - std::sqrt(enth*enth + 4*physical::mu*physical::si*physical::Lf/physical::C) ); // (Winton, 2000, eq 38)
                    }
                    else if (vptr->hasValueNoThickIce())
                        val = vptr->valueNoThickIce();
                }

                val = vptr->hasMinVal() ? std::max(vptr->minVal(), val ) : val ;
                if(apply_maxima)
                    val = has_max[j] ? std::min(vptr->maxVal(), val ) : val ;
                (*vptr)[i] = val;
            }

            if(apply_maxima && M_ice_cat_type==setup::IceCategoryType::YOUNG_ICE)
            {
                if ((M_conc[i] + M_conc_young[i]) > 1.) M_conc_young[i] = 1. - M_conc[i];
            }
        }//loop over i
    }//redistributeVariables


//------------------------------------------------------------------------------------------------------
//! Redistributes variables (parallel computing).
//! Called by function scatterFieldsElementIO().
//! * elt_values_local is vector containing all the variables to be redistributed (eg after scattering from root)
//    into the individual variables (eg M_conc, M_thick,...)
//! * data is a vector of pointers to the variables to be assigned values from elt_values_local
    void
    ScatterGather::redistributeVariablesIO(std::vector<double> const& elt_values_local,
                                           std::vector<ModelVariable*> &data)
    {
        int nb_var_element = data.size();
        for(int j=0; j<data.size(); j++)
        {
            for (int i=0; i<M_num_elements; ++i)
            {
                auto ptr = data[j];
                (*ptr)[i] = elt_values_local[nb_var_element*i+j];
            }
        }
    }//redistributeVariablesIO



}