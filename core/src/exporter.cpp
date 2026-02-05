/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim: set fenc=utf-8 ft=cpp et sw=4 ts=4 sts=4: */

/**
 * @file   exporter.cpp
 * @author Abdoulaye Samake <abdoulaye.samake@nersc.no>
 * @author Sylvain Bouillon <sylvain.bouillon@nersc.no>
 * @author Einar Olason <einar.olason@nersc.no>
 * @date   Mon Dec 28 15:59:43 2015
 */

#include <exporter.hpp>

namespace Nextsim
{
Exporter::Exporter()
	:
	M_mrecord(),
    M_frecord(),
    M_type_record(),
    M_name_record()
{}

// writeContainer in sequential
template<typename Type>
void
Exporter::writeContainer(std::fstream& out, std::vector<Type> const& container, std::string const precision)
{
    if (out.is_open())
	{
        int fsize = container.size();
        out.write((char*)&fsize, sizeof(fsize)); // write first the record length

        int typesize = sizeof(Type);
        if ( precision == "float" )
        {
            for (int i=0; i<container.size(); ++i)
            {
                // convert to float before writing into binary file
                float fvalue = (float)container[i];
                out.write((char*)&fvalue, sizeof(float));
            }
        }
        else
        {
            // write as original format
            out.write((char*)&container[0], container.size() * sizeof(Type));
        }
	}
	else
	{
		std::cerr << "writeContainer error: opening file for output failed!" <<"\n";
		std::abort();
	}
}

// writeContainer with parallel writing
template<typename Type>
void
Exporter::writeContainer(MPI_File& out, std::vector<Type> const& container, int is_vector, 
                         Communicator M_comm, MPI_Offset& base_offset, std::string const precision)
{
    int fsize = container.size();
    int global_fsize;
    boost::mpi::all_reduce(M_comm, fsize, global_fsize, std::plus<int>());

    // Write the size only on rank 0
    if (M_comm.rank() == 0) MPI_File_write_at(out, base_offset, &global_fsize, 1 ,MPI_INT, MPI_STATUS_IGNORE);

    // Offset is the location to write in the output file
    MPI_Offset offset = static_cast<MPI_Offset>(sizeof(global_fsize)) + base_offset;

    MPI_Datatype mpi_type;
    size_t elem_size = 0;
    if (precision == "float")      { mpi_type = MPI_FLOAT;  elem_size = sizeof(float); }
    else if (precision == "double"){ mpi_type = MPI_DOUBLE; elem_size = sizeof(double); }
    else if (precision == "int")   { mpi_type = MPI_INT;    elem_size = sizeof(int); }

    if (!is_vector)
    {
        // Find the offset for each processor
        int offset_count = 0;
        MPI_Exscan(&fsize, &offset_count, 1, MPI_INT, MPI_SUM, MPI_Comm(M_comm));
        if (M_comm.rank() == 0) offset_count = 0;

        offset += static_cast<MPI_Offset>(offset_count) * elem_size;

        std::vector<float> fvalues(fsize);
        if (precision == "float") {
            for (int i=0; i<fsize; ++i) fvalues[i] = (float)container[i];
            MPI_File_write_at(out, offset, fvalues.data(), fsize, mpi_type, MPI_STATUS_IGNORE);
        }
        else MPI_File_write_at(out, offset, container.data(), fsize, mpi_type, MPI_STATUS_IGNORE);

        base_offset += sizeof(global_fsize) + global_fsize * elem_size;

        return;
    }

    // else, the container is a vector with two components
    // it must be diveded in two parts and write separatly
    int fsize2 = fsize/2;
    int offset_count1 = 0;
    MPI_Exscan(&fsize2, &offset_count1, 1, MPI_INT, MPI_SUM, MPI_Comm(M_comm));
    if (M_comm.rank() == 0) offset_count1 = 0;

    offset += static_cast<MPI_Offset>(offset_count1) * elem_size;

    if (precision == "float")
    {
        std::vector<float> fvalues(fsize2);
        for (int i=0; i<fsize2; ++i) fvalues[i] = (float)container[2*i];
        MPI_File_write_at(out, offset, fvalues.data(), fsize2, mpi_type, MPI_STATUS_IGNORE);
    }
    else
    {
        std::vector<Type> fvalues(fsize2);
        for (int i=0; i<fsize2; ++i) fvalues[i] = container[2*i];
        MPI_File_write_at(out, offset, fvalues.data(), fsize2, mpi_type, MPI_STATUS_IGNORE);
    }

    MPI_Offset offset2 = static_cast<MPI_Offset>(sizeof(global_fsize)) + base_offset +
                         (global_fsize/2) * elem_size +
                         static_cast<MPI_Offset>(offset_count1) * elem_size;

    if (precision == "float")
    {
        std::vector<float> fvalues(fsize2);
        for (int i=0; i<fsize2; ++i) fvalues[i] = (float)container[2*i+1];
        MPI_File_write_at(out, offset2, fvalues.data(), fsize2, mpi_type, MPI_STATUS_IGNORE);
    }
    else
    {
        std::vector<Type> fvalues(fsize2);
        for (int i=0; i<fsize2; ++i) fvalues[i] = container[2*i+1];
        MPI_File_write_at(out, offset2, fvalues.data(), fsize2, mpi_type, MPI_STATUS_IGNORE);
    }

    // Update the base offset for the next fields to write in the same file
    base_offset += sizeof(global_fsize) + global_fsize * elem_size;

} //writeContainer

template<typename Type>
void
Exporter::orderField(std::vector<Type> const& field_local, std::vector<Type>& ordered_field, std::vector<int>& rmap, 
                     Communicator M_comm, int size_without_ghost, int num_nodes)
{
    int is_vector = 0;
    if (field_local.size() == 2*num_nodes) is_vector = 1;

    std::vector<int> rcounts(M_comm.size());
    boost::mpi::all_gather(M_comm, size_without_ghost, rcounts);
    std::vector<int> displs(M_comm.size());
    displs[0] = 0;
    for ( int k = 1; k < M_comm.size(); k++ ) {
        displs[k] = displs[k-1] + rcounts[k-1];
    }

    std::vector<std::vector<int>> list_indices_recv(M_comm.size());
    for (int i = 0; i < size_without_ghost; i++)
    {
        int ri = rmap[displs[M_comm.rank()] + i];
        int j = 0;
        while (j < M_comm.size()-1 && ri >= displs[j+1]) j++;

        if (j == M_comm.rank())
        {
            if (is_vector)
            {
                ordered_field[2*i] = field_local[ri - displs[M_comm.rank()]];
                ordered_field[2*i+1] = field_local[ri - displs[M_comm.rank()] + num_nodes];
            }
            else ordered_field[i] = field_local[ri - displs[M_comm.rank()]];
            continue;
        }

        list_indices_recv[j].push_back(displs[M_comm.rank()] + i);
    }

    std::vector<boost::mpi::request> requests;
    std::vector<std::vector<int>> list_indices_send(M_comm.size());
    for (int i = 0; i < M_comm.size(); i++)
    {
        requests.push_back(M_comm.isend(i, M_comm.rank(), list_indices_recv[i]));
        requests.push_back(M_comm.irecv(i, i, list_indices_send[i]));
    }
    boost::mpi::wait_all(requests.begin(), requests.end());
    requests.clear();

    std::vector<std::vector<Type>> field_send(M_comm.size());
    for (int i = 0; i < M_comm.size(); i++)
    {
        for (int j = 0; j < list_indices_send[i].size(); j++)
        {
            field_send[i].push_back(field_local[rmap[list_indices_send[i][j]] - displs[M_comm.rank()]]);
            if (is_vector) field_send[i].push_back(field_local[rmap[list_indices_send[i][j]] - displs[M_comm.rank()] + num_nodes]);
        }
    }

    std::vector<std::vector<Type>> field_recv(M_comm.size());
    for (int i = 0; i < M_comm.size(); i++)
    {
        requests.push_back(M_comm.isend(i, M_comm.rank(), field_send[i]));
        requests.push_back(M_comm.irecv(i, i, field_recv[i]));
    }
    boost::mpi::wait_all(requests.begin(), requests.end());

    for (int i = 0; i < M_comm.size(); i++)
    {
        if (M_comm.rank() == i) continue;

        for (int j = 0; j < list_indices_recv[i].size(); j++)
        {
            if (is_vector)
            {
                ordered_field[2*(list_indices_recv[i][j]-displs[M_comm.rank()])] = field_recv[i][2*j];
                ordered_field[2*(list_indices_recv[i][j]-displs[M_comm.rank()])+1] = field_recv[i][2*j+1];
            }
            else ordered_field[list_indices_recv[i][j]-displs[M_comm.rank()]] = field_recv[i][j];
        }
    }

}//orderField

// writeMesh function in parallel
void
Exporter::writeMesh(MPI_File& out, GmshMesh const& Mesh, std::vector<int>& rmap_nodes, std::vector<int>& rmap_elements)
{
    Communicator M_comm = Mesh.comm();

    // Local 
    std::vector<double> coordX = Mesh.coordXPartition();
    std::vector<double> coordY = Mesh.coordYPartition();
    std::vector<int> indexTr = Mesh.indexTrPartition();
    std::vector<int> id = Mesh.id();
    int num_nodes = coordX.size();
    id.resize(num_nodes); // Remove ghosts

    // Global coordinates
    std::vector<double> coordX_out(num_nodes);
    std::vector<double> coordY_out(num_nodes);
    std::vector<int> id_out(num_nodes);
    orderField(coordX, coordX_out, rmap_nodes, M_comm, num_nodes, 0);
    orderField(coordY, coordY_out, rmap_nodes, M_comm, num_nodes, 0);
    orderField(id, id_out, rmap_nodes, M_comm, num_nodes, 0);

    // Global indices of triangles
    std::vector<int> size_nodes;
    boost::mpi::all_gather(M_comm, Mesh.numLocalNodesWithoutGhost(), size_nodes);

    std::vector<int> displs(M_comm.size());
    displs[0] = 0;
    for ( int k = 1; k < M_comm.size(); k++ ) {
        displs[k] = displs[k-1] + size_nodes[k-1];
    }

    std::vector<int> id_root(Mesh.numGlobalNodes());
    std::vector<int> id_ghost(Mesh.numGlobalNodes()+1);
    int ier = MPI_Allgatherv(&id[0], Mesh.numLocalNodesWithoutGhost(), MPI_INT, &id_root[0], 
                             &size_nodes[0], &displs[0], MPI_INT, MPI_Comm(M_comm));

    for (int k = 0; k < Mesh.numGlobalNodes(); k++) id_ghost[id_root[k]] = k;

    // Correct the local indices in the triangles to have global indices
    std::map<int, Nextsim::entities::GMSHPoint > nodes = Mesh.nodes();
    auto it = nodes.begin();
    std::advance(it, num_nodes);

    std::vector<int> id_local_ghost(Mesh.numNodes() - num_nodes);
    int cpt = 0;
    for (; it != nodes.end(); ++it)
    {
        id_local_ghost[cpt] = it->second.id;
        ++cpt;
    }

    for (int k = 0; k < 3*Mesh.numTrianglesWithoutGhost(); k++)
    {
        if (indexTr[k] <= Mesh.numLocalNodesWithoutGhost())
        {
            indexTr[k] += displs[M_comm.rank()];
        }
        else // it is a ghost that belongs to another partition
        {
            indexTr[k] = id_ghost[id_local_ghost[indexTr[k]-num_nodes-1]]+1;
        }
    }

    // Correct the global indices with the root node ordering
    // And split the indexTr array to use the orderField function
    std::vector<int> node_correspondence(rmap_nodes.size());
    for (int i = 0; i < rmap_nodes.size(); ++i) node_correspondence[rmap_nodes[i]] = i;

    int n_triangles = indexTr.size()/3;
    std::vector<int> indexTr1(n_triangles);
    std::vector<int> indexTr2(n_triangles);
    std::vector<int> indexTr3(n_triangles);
    for (int i = 0; i < n_triangles; i++)
    {
        indexTr1[i] = node_correspondence[indexTr[3*i]-1]+1;
        indexTr2[i] = node_correspondence[indexTr[3*i+1]-1]+1;
        indexTr3[i] = node_correspondence[indexTr[3*i+2]-1]+1;
    }

    std::vector<int> indexTr1_out(n_triangles);
    std::vector<int> indexTr2_out(n_triangles);
    std::vector<int> indexTr3_out(n_triangles);
    orderField(indexTr1, indexTr1_out, rmap_elements, Mesh.comm(), n_triangles, 0);
    orderField(indexTr2, indexTr2_out, rmap_elements, Mesh.comm(), n_triangles, 0);
    orderField(indexTr3, indexTr3_out, rmap_elements, Mesh.comm(), n_triangles, 0);

    std::vector<int> indexTr_out(3*n_triangles);
    for (int i = 0; i < n_triangles; i++)
    {
        indexTr_out[3*i] = indexTr1_out[i];
        indexTr_out[3*i+1] = indexTr2_out[i];
        indexTr_out[3*i+2] = indexTr3_out[i];
    }

    // Write the mesh
    this->writeMesh(out, coordX_out, coordY_out, id_out, indexTr_out, M_comm);

}//writeMesh

template<typename Type>
void
Exporter::writeMesh(MPI_File& out, std::vector<Type> const &xnod, std::vector<Type> const &ynod,
                    std::vector<int> const &idnod, std::vector<int> const &elements, Communicator M_comm)
{
    std::string description;
    MPI_Offset base_offset = 0;

    std::vector<int> field_int = elements;
    writeContainer(out, field_int, 0, M_comm, base_offset, "int");

    int min_element = *std::min_element(field_int.begin(), field_int.end());
    int max_element = *std::max_element(field_int.begin(), field_int.end());
    int min_element_global, max_element_global, global_element_size;
    boost::mpi::reduce(M_comm, min_element, min_element_global, boost::mpi::minimum<int>(), 0);
    boost::mpi::reduce(M_comm, max_element, max_element_global, boost::mpi::maximum<int>(), 0);
    boost::mpi::reduce(M_comm, int(field_int.size()), global_element_size, std::plus<int>(), 0);

    if (M_comm.rank() == 0)
    {
        description = (boost::format( "%1% %2% %3$g %4$g %5$g" )
                       % "Elements"
                       % "int"
                       % global_element_size
                       % min_element_global
                       % max_element_global).str();
        M_mrecord.push_back(description);
    }

    field_int = idnod;
    writeContainer(out, field_int, 0, M_comm, base_offset, "int");

    int min_id = *std::min_element(field_int.begin(), field_int.end());
    int max_id = *std::max_element(field_int.begin(), field_int.end());
    int min_id_global, max_id_global, global_id_size;
    boost::mpi::reduce(M_comm, min_id, min_id_global, boost::mpi::minimum<int>(), 0);
    boost::mpi::reduce(M_comm, max_id, max_id_global, boost::mpi::maximum<int>(), 0);
    boost::mpi::reduce(M_comm, int(field_int.size()), global_id_size, std::plus<int>(), 0);

    if (M_comm.rank() == 0)
    {
        description = (boost::format( "%1% %2% %3$g %4$g %5$g" )
                       % "id"
                       % "int"
                       % global_id_size
                       % min_id_global
                       % max_id_global).str();
        M_mrecord.push_back(description);
    }

    std::string prec = "double";
    if (sizeof(Type)==sizeof(float))
        prec = "float";

    auto field = xnod;
    writeContainer(out, field, 0, M_comm, base_offset, prec);

    int min_x = *std::min_element(field.begin(), field.end());
    int max_x = *std::max_element(field.begin(), field.end());
    int min_x_global, max_x_global, global_x_size;
    boost::mpi::reduce(M_comm, min_x, min_x_global, boost::mpi::minimum<double>(), 0);
    boost::mpi::reduce(M_comm, max_x, max_x_global, boost::mpi::maximum<double>(), 0);
    boost::mpi::reduce(M_comm, int(field.size()), global_x_size, std::plus<int>(), 0);

    if (M_comm.rank() == 0)
    {
        description = (boost::format( "%1% %2% %3$g %4$g %5$g" )
                       % "Nodes_x"
                       % prec
                       % global_x_size
                       %  min_x_global
                       % max_x_global ).str();
        M_mrecord.push_back(description);
    }

    field = ynod;
    writeContainer(out, field, 0, M_comm, base_offset, prec);

    int min_y = *std::min_element(field.begin(), field.end());
    int max_y = *std::max_element(field.begin(), field.end());
    int min_y_global, max_y_global, global_y_size;
    boost::mpi::reduce(M_comm, min_y, min_y_global, boost::mpi::minimum<double>(), 0);
    boost::mpi::reduce(M_comm, max_y, max_y_global, boost::mpi::maximum<double>(), 0);
    boost::mpi::reduce(M_comm, int(field.size()), global_y_size, std::plus<int>(), 0);

    if (M_comm.rank() == 0)
    {
        description = (boost::format( "%1% %2% %3$g %4$g %5$g" )
                       % "Nodes_y"
                       % prec
                       % global_y_size
                       % min_y_global
                       % max_y_global ).str();
        M_mrecord.push_back(description);
    }
}

// writeMesh sequentially
void
Exporter::writeMesh(std::fstream& out, GmshMeshSeq const& Mesh)
{
    auto xnod = Mesh.coordX();
    auto ynod = Mesh.coordY();
    auto elements = Mesh.indexTr();
    auto idnod = Mesh.id();
    this->writeMesh(out, xnod, ynod, idnod, elements);
}

template<typename Type>
void
Exporter::writeMesh(std::fstream& out, std::vector<Type> const &xnod, std::vector<Type> const &ynod,
        std::vector<int> const &idnod, std::vector<int> const &elements)
{
    std::string description;

    std::vector<int> field_int = elements;
    writeContainer(out, field_int, "int");
    description = (boost::format( "%1% %2% %3$g %4$g %5$g" )
                   % "Elements"
                   % "int"
                   % field_int.size()
                   % *std::min_element(field_int.begin(), field_int.end())
                   % *std::max_element(field_int.begin(), field_int.end()) ).str();
    M_mrecord.push_back(description);

    field_int = idnod;
    writeContainer(out, field_int, "int");
    description = (boost::format( "%1% %2% %3$g %4$g %5$g" )
                   % "id"
                   % "int"
                   % field_int.size()
                   % *std::min_element(field_int.begin(), field_int.end())
                   % *std::max_element(field_int.begin(), field_int.end()) ).str();
    M_mrecord.push_back(description);

    std::string prec = "double";
    if (sizeof(Type)==sizeof(float))
        prec = "float";

    auto field = xnod;
    writeContainer(out, field, prec);
    description = (boost::format( "%1% %2% %3$g %4$g %5$g" )
                   % "Nodes_x"
                   % prec
                   % field.size()
                   % *std::min_element(field.begin(), field.end())
                   % *std::max_element(field.begin(), field.end()) ).str();
    M_mrecord.push_back(description);

    field = ynod;
    writeContainer(out, field, prec);
    description = (boost::format( "%1% %2% %3$g %4$g %5$g" )
                   % "Nodes_y"
                   % prec
                   % field.size()
                   % *std::min_element(field.begin(), field.end())
                   % *std::max_element(field.begin(), field.end()) ).str();

    M_mrecord.push_back(description);
}

//writeField in parallel
template<typename Type>
void
Exporter::writeField(MPI_File& out, std::vector<Type> const& field, std::string const& name, std::string const& precision, Communicator M_comm, 
                     MPI_Offset& base_offset, std::vector<int>& rmap, int size_without_ghost, int num_nodes, int root)
{
    std::string description;

    MPI_Datatype mpi_type;
    size_t elem_size = 0;
    if (precision == "float")      { mpi_type = MPI_FLOAT;  elem_size = sizeof(float); }
    else if (precision == "double"){ mpi_type = MPI_DOUBLE; elem_size = sizeof(double); }
    else if (precision == "int")   { mpi_type = MPI_INT;    elem_size = sizeof(int); }
    else {
        std::string const error_str = "Exporter: Unknown precision: " + precision;
        throw std::logic_error(error_str);
    }

    // For the fields written only by the root process
    if (root)
    {
        if (M_comm.rank() == 0) 
        {
            int size = field.size();
            MPI_File_write_at(out, base_offset, &size, 1, MPI_INT, MPI_STATUS_IGNORE);
            base_offset += sizeof(size);
            MPI_File_write_at(out, base_offset, field.data(), size, mpi_type, MPI_STATUS_IGNORE);
            base_offset += size * elem_size;

            description = (boost::format( "%1% %2% %3$g %4$g %5$g" )
                       % name
                       % precision
                       % field.size()
                       % ( (field.size() > 0) ? *std::min_element(field.begin(), field.end()) : 0 )
                       % ( (field.size() > 0) ? *std::max_element(field.begin(), field.end()) : 0 ) ).str();

            M_frecord.push_back(description);
        }

        boost::mpi::broadcast(M_comm, &base_offset, 1, 0);
        return;
    }

    std::vector<Type> order_field;
    int is_vector = 0;
    if (field.size() == 2*num_nodes) // it is a vector node field
    {
        order_field.resize(2*size_without_ghost);
        is_vector = 1;
    }
    else if (field.size() == num_nodes) // it is a scalar node field
    {
        order_field.resize(size_without_ghost);
    }
    else // it is an element field
    {
        order_field.resize(field.size());
    }

    orderField(field, order_field, rmap, M_comm, size_without_ghost, num_nodes);
    writeContainer(out, order_field, is_vector, M_comm, base_offset, precision);

    int min_element = *std::min_element(field.begin(), field.end());
    int max_element = *std::max_element(field.begin(), field.end());
    int min_element_global, max_element_global, global_element_size;
    boost::mpi::reduce(M_comm, min_element, min_element_global, boost::mpi::minimum<int>(), 0);
    boost::mpi::reduce(M_comm, max_element, max_element_global, boost::mpi::maximum<int>(), 0);
    boost::mpi::reduce(M_comm, int(order_field.size()), global_element_size, std::plus<int>(), 0);

    if (M_comm.rank() == 0)
    {
        description = (boost::format( "%1% %2% %3$g %4$g %5$g" )
                       % name
                       % precision
                       % global_element_size
                       % min_element_global
                       % max_element_global).str();

        M_frecord.push_back(description);
    }

}//writeField


void
Exporter::writeRecord(std::fstream& out, std::string const& rtype)
{
	if (out.is_open())
	{
        if (rtype == "field")
        {
            for (auto const& str : M_frecord)
            {
                //out << std::setw(15) << std::left << R.name;
                out << str <<"\n";
            }
        }
        else
        {
            for (auto const& str : M_mrecord)
            {
                //out << std::setw(15) << std::left << R.name;
                out << str <<"\n";
            }
        }
	}
	else
	{
		std::cerr << "writeRecord error: opening file for output failed!" <<"\n";
		std::abort();
	}
}

void
Exporter::readRecord(std::ifstream &in)
{
    std::string name, type, size, min_element, max_element;

    while ( in >> name >> type >> size >> min_element >> max_element)
    {
        M_name_record.push_back(name);
        M_type_record.push_back(type);
    }
}

void
Exporter::loadFile(std::fstream &in, boost::unordered_map<std::string, std::vector<int>> &field_map_int, boost::unordered_map<std::string, std::vector<double>> &field_map_dbl)
{
    int reclen;
    for ( int i=0; i<M_type_record.size(); ++i )
    {
        in.read((char*) &reclen, sizeof(reclen));

        if ( M_type_record[i] == "double" )
        {
            std::vector<double> dvec(reclen);
            in.read((char*) &dvec[0], reclen*sizeof(double));
            field_map_dbl.emplace(M_name_record[i], dvec);
        }
        else if ( M_type_record[i] == "int" )
        {
            std::vector<int> ivec(reclen);
            in.read((char*) &ivec[0], reclen*sizeof(int));
            field_map_int.emplace(M_name_record[i], ivec);
        }
        else
        {
            throw std::logic_error("unknown type in file");
        }
    }
}


template void Exporter::writeField(MPI_File& out, std::vector<double> const& field, std::string const& name, std::string const& precision, Communicator M_comm,
                                   MPI_Offset& base_offset, std::vector<int>& rmap, int size_without_ghost, int num_nodes, int root);
template void Exporter::writeField(MPI_File& out, std::vector<int> const& field, std::string const& name, std::string const& precision, Communicator M_comm,
                                   MPI_Offset& base_offset, std::vector<int>& rmap, int size_without_ghost, int num_nodes, int root);

template void Exporter::writeMesh<double>(std::fstream&, std::vector<double> const&, std::vector<double> const&,
        std::vector<int> const&,std::vector<int> const&);
template void Exporter::writeMesh<float>(std::fstream&, std::vector<float> const&, std::vector<float> const&,
        std::vector<int> const&,std::vector<int> const&);

} // Nextsim
