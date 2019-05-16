/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim: set fenc=utf-8 ft=cpp et sw=4 ts=4 sts=4: */

/**
 * @file   gmshmodel.hpp
 * @author Abdoulaye Samake <abdoulaye.samake@nersc.no>
 * @date   Thu May  2 12:43:36 2019
 */

#ifndef __GmshModel_HPP
#define __GmshModel_HPP 1

class GModel;
#include <GModel.h>
#include <meshPartition.h>

namespace Nextsim
{

class GmshModel : public GModel
{
    typedef GModel super;

public:

GmshModel()
    :
    super()
    {}

GModel* createGModel(std::map<int, MVertex*>& vertexMap,
                     std::vector<int> const& elementNum,
                     std::vector<std::vector<int> >& vertexIndices,
                     std::vector<int> const& elementType,
                     std::vector<int> const& physical,
                     std::vector<int> const& elementary,
                     std::vector<int> const& partition)
{
    int numVertices = (int)vertexMap.size();
    int numElement = (int)elementNum.size();

    if(numElement != (int)vertexIndices.size()){
        Msg::Error("Dimension in vertices numbers");
        return 0;
    }
    if(numElement != (int)elementType.size()){
        Msg::Error("Dimension in elementType numbers");
        return 0;
    }
    if(numElement != (int)physical.size()){
        Msg::Error("Dimension in physical numbers");
        return 0;
    }
    if(numElement != (int)elementary.size()){
        Msg::Error("Dimension in elementary numbers");
        return 0;
    }
    if(numElement != (int)partition.size()){
        Msg::Error("Dimension in partition numbers");
        return 0;
    }

    // GModel *gm = new GModel();

    GModel *gm = this;

    std::map<int, std::vector<MElement*> > elements[11];
    std::map<int, std::map<int, std::string> > physicals[4];
    std::vector<MVertex*> vertexVector;

    std::map<int, MVertex*>::const_iterator it = vertexMap.begin();
    std::map<int, MVertex*>::const_iterator end = vertexMap.end();

    int maxVertex = std::numeric_limits<int>::min();
    int minVertex = std::numeric_limits<int>::max();
    int num;

    for(; it != end; ++it){
        num = it->first;
        minVertex = std::min(minVertex, num);
        maxVertex = std::max(maxVertex, num);
    }
    if(minVertex == std::numeric_limits<int>::max())
        Msg::Error("Could not determine the min index of vertices");

    // if the vertex numbering is dense, transfer the map into a vector to speed
    // up element creation
    if((minVertex == 1 && maxVertex == numVertices) ||
       (minVertex == 0 && maxVertex == numVertices - 1)){
        Msg::Info("Vertex numbering is dense");
        vertexVector.resize(vertexMap.size() + 1);
        if(minVertex == 1)
            vertexVector[0] = 0;
        else
            vertexVector[numVertices] = 0;
        std::map<int, MVertex*>::const_iterator it = vertexMap.begin();
        for(; it != vertexMap.end(); ++it)
            vertexVector[it->first] = it->second;
        vertexMap.clear();
    }
    int *indices;
    int nbVertices;
    for(int i = 0; i < numElement; ++i){
        num = elementNum[i];
        std::vector<MVertex*> vertices;
        nbVertices = (int)vertexIndices[i].size();
        indices = &vertexIndices[i][0];

        MElementFactory f;
        MElement *e = f.create(elementType[i], vertices, num, partition[i]);
        if(!e){
            Msg::Error("Unknown type of element %d", elementType[i]);
            return 0;
        }
        switch(e->getType()){
        case TYPE_PNT : elements[0][elementary[i]].push_back(e); break;
        case TYPE_LIN : elements[1][elementary[i]].push_back(e); break;
        case TYPE_TRI : elements[2][elementary[i]].push_back(e); break;
        case TYPE_QUA : elements[3][elementary[i]].push_back(e); break;
        case TYPE_TET : elements[4][elementary[i]].push_back(e); break;
        case TYPE_HEX : elements[5][elementary[i]].push_back(e); break;
        case TYPE_PRI : elements[6][elementary[i]].push_back(e); break;
        case TYPE_PYR : elements[7][elementary[i]].push_back(e); break;
        case TYPE_TRIH : elements[8][elementary[i]].push_back(e); break;
        case TYPE_POLYG: elements[9][elementary[i]].push_back(e); break;
        case TYPE_POLYH: elements[10][elementary[i]].push_back(e); break;
        default : Msg::Error("Wrong type of element"); return 0;
        }

        int dim = e->getDim();
        if(physical[i] && (!physicals[dim].count(elementary[i]) ||
                           !physicals[dim][elementary[i]].count(physical[i])))
            physicals[dim][elementary[i]][physical[i]] = "unnamed";
    }

    // store the elements in their associated elementary entity. If the
    // entity does not exist, create a new (discrete) one.
    for(int i = 0; i < (int)(sizeof(elements) / sizeof(elements[0])); i++)
        this->_storeElementsInEntities(elements[i]);

    // associate the correct geometrical entity with each mesh vertex
    this->_associateEntityWithMeshVertices();

    // store the vertices in their associated geometrical entity
    if(vertexVector.size())
        this->_storeVerticesInEntities(vertexVector);
    else
        this->_storeVerticesInEntities(vertexMap);

    // store the physical tags
    for(int i = 0; i < 4; i++)
        this->_storePhysicalTagsInEntities(i, physicals[i]);

    return gm;
}

};
} //Nextsim

#endif
