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
Exporter::Exporter(std::string const& precision)
	:
	M_mrecord(),
    M_frecord(),
    M_type_record(),
    M_name_record(),
    M_precision(precision)
{
    if ( precision != "double" & precision != "float" & precision != "int" )
    {
        std::string error_str = "Exporter: Unknown precision: " + precision;
        throw std::logic_error(error_str);
    }
}

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

void
Exporter::writeMesh(std::fstream& out, GmshMesh const& Mesh)
{
    Communicator M_comm = Mesh.comm();

    std::vector<double> coordX = Mesh.coordXPartition();
    std::vector<double> coordY = Mesh.coordYPartition();
    std::vector<int> indexTr = Mesh.indexTrPartition();
    std::vector<int> id = Mesh.id();

    std::vector<double> coordX_root, coordY_root;
    std::vector<int> indexTr_root, id_root;
    std::vector<int> id_ghost(Mesh.numGlobalNodes()+1);

    std::vector<int> size_elements;
    boost::mpi::all_gather(M_comm, 3*Mesh.numTrianglesWithoutGhost(), size_elements);

    std::vector<int> size_nodes;
    boost::mpi::all_gather(M_comm, Mesh.numLocalNodesWithoutGhost(), size_nodes);

    std::vector<int> displs(M_comm.size());
    displs[0] = 0;
    for ( int k = 1; k < M_comm.size(); k++ ) {
        displs[k] = displs[k-1] + size_nodes[k-1];
    }

    if (M_comm.rank() == 0)
    {
        coordX_root.resize(Mesh.numGlobalNodes());
        coordY_root.resize(Mesh.numGlobalNodes());
        indexTr_root.resize(3*Mesh.numGlobalElements());
        id_root.resize(Mesh.numGlobalNodes());
    }

    // Gather coordinates and id
    int ier = MPI_Gatherv(&coordX[0], Mesh.numLocalNodesWithoutGhost(), MPI_DOUBLE, &coordX_root[0], &size_nodes[0], &displs[0], MPI_DOUBLE, 0, MPI_Comm(M_comm));
    ier = MPI_Gatherv(&coordY[0], Mesh.numLocalNodesWithoutGhost(), MPI_DOUBLE, &coordY_root[0], &size_nodes[0], &displs[0], MPI_DOUBLE, 0, MPI_Comm(M_comm));
    ier = MPI_Gatherv(&id[0], Mesh.numLocalNodesWithoutGhost(), MPI_INT, &id_root[0], &size_nodes[0], &displs[0], MPI_INT, 0, MPI_Comm(M_comm));

    if (M_comm.rank() == 0)
    {
        for (int k = 0; k < Mesh.numGlobalNodes(); k++) id_ghost[id_root[k]] = k;
    }

    boost::mpi::broadcast(M_comm, &id_ghost[0], Mesh.numGlobalNodes()+1, 0);

    // Correct the indices in the triangles
    std::map<int, Nextsim::entities::GMSHPoint > nodes = Mesh.nodes();
    for (int k = 0; k < 3*Mesh.numTrianglesWithoutGhost(); k++)
    {
        if (indexTr[k] <= Mesh.numLocalNodesWithoutGhost())
        {
            indexTr[k] += displs[M_comm.rank()];
        }
        else // it is a ghost that belongs to another partition
        {
            indexTr[k] = id_ghost[nodes.find(indexTr[k])->second.id]+1;
        }
    }

    // Gather triangle indices
    displs[0] = 0;
    for (int k = 1; k < M_comm.size(); k++) displs[k] = displs[k-1] + size_elements[k-1];

    ier = MPI_Gatherv(&indexTr[0], 3*Mesh.numTrianglesWithoutGhost(), MPI_INT, &indexTr_root[0], &size_elements[0], &displs[0], MPI_INT, 0, MPI_Comm(M_comm));

    if (M_comm.rank() == 0)
    {
        std::vector<double> coordX_out(Mesh.numGlobalNodes());
        std::vector<double> coordY_out(Mesh.numGlobalNodes());
        std::vector<int> indexTr_out(indexTr_root.size());
        std::vector<int> id_out(Mesh.numGlobalNodes());
        std::vector<int> node_correspondence(Mesh.numGlobalNodes());

        std::vector<int> M_rmap_nodes = Mesh.mapNodes();
        std::vector<int> M_rmap_elements = Mesh.mapElements();

        // Reorder nodes to correspond to the node fields
        for (int i = 0; i < Mesh.numGlobalNodes(); ++i)
        {
            int ri =  M_rmap_nodes[i];
            node_correspondence[ri] = i;

            coordX_out[i] = coordX_root[ri];
            coordY_out[i] = coordY_root[ri];
            id_out[i] = id_root[ri];
        }

        for (int i = 0; i < indexTr_root.size(); i++) indexTr_out[i] = node_correspondence[indexTr_root[i]-1]+1;

        // Reorder the elements to correspond to the element fields
        auto indexTr_out_cpy = indexTr_out;

        for (int i = 0; i < indexTr_root.size(); i++)
        {
            int ri = M_rmap_elements[i/3];
            indexTr_out[i] = indexTr_out_cpy[3*ri+i%3];
        }

        // Write the mesh
        this->writeMesh(out, coordX_out, coordY_out, id_out, indexTr_out);
    }

}

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

template<typename Type>
void
Exporter::writeField(std::fstream& out, std::vector<Type> const& field, std::string const& name)
{
    std::string description;

    std::string precision;
    // We must choose either integer or floating point
    if (((boost::any)field[0]).type() == typeid(int))
        precision = "int";
    else
        precision = M_precision;

    // Time should always be in double
    if ( name == "Time" )
        precision = "double";

    writeContainer(out, field, precision);
    description = (boost::format( "%1% %2% %3$g %4$g %5$g" )
                   % name
                   % precision
                   % field.size()
                   % ( (field.size() > 0) ? *std::min_element(field.begin(), field.end()) : 0 )
                   % ( (field.size() > 0) ? *std::max_element(field.begin(), field.end()) : 0 ) ).str();

	M_frecord.push_back(description);
}

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

template void Exporter::writeField<double>(std::fstream&, std::vector<double> const&, std::string const&);
template void Exporter::writeField<int>(std::fstream&, std::vector<int> const&, std::string const&);

template void Exporter::writeMesh<double>(std::fstream&, std::vector<double> const&, std::vector<double> const&,
        std::vector<int> const&,std::vector<int> const&);
template void Exporter::writeMesh<float>(std::fstream&, std::vector<float> const&, std::vector<float> const&,
        std::vector<int> const&,std::vector<int> const&);

} // Nextsim
