#include "./include.h"
#include "./shared.h"

/*Constructors/Destructors*/
BamgMesh::BamgMesh(){/*{{{*/

	this->VerticesSize[0]=0,                  this->VerticesSize[1]=0;                 this->Vertices=NULL;          this->PreviousNumbering = NULL;
	this->EdgesSize[0]=0,                     this->EdgesSize[1]=0;                    this->Edges=NULL;
	this->TrianglesSize[0]=0,                 this->TrianglesSize[1]=0;                this->Triangles=NULL;

	this->SubDomainsSize[0]=0,                this->SubDomainsSize[1]=0;               this->SubDomains=NULL;
	this->SubDomainsFromGeomSize[0]=0,        this->SubDomainsFromGeomSize[1]=0;       this->SubDomainsFromGeom=NULL;

	this->VerticesOnGeomVertexSize[0]=0,      this->VerticesOnGeomVertexSize[1]=0;     this->VerticesOnGeomVertex=NULL;
	this->VerticesOnGeomEdgeSize[0]=0,        this->VerticesOnGeomEdgeSize[1]=0;       this->VerticesOnGeomEdge=NULL;
	this->EdgesOnGeomEdgeSize[0]=0,           this->EdgesOnGeomEdgeSize[1]=0;          this->EdgesOnGeomEdge=NULL;

	this->IssmEdgesSize[0]=0,                 this->IssmEdgesSize[1]=0;                this->IssmEdges=NULL;
	this->IssmSegmentsSize[0]=0,              this->IssmSegmentsSize[1]=0;             this->IssmSegments=NULL;

	this->ElementConnectivitySize[0]=0,       this->ElementConnectivitySize[1]=0;      this->ElementConnectivity=NULL;
	this->NodalConnectivitySize[0]=0,         this->NodalConnectivitySize[1]=0;        this->NodalConnectivity=NULL;
	this->NodalElementConnectivitySize[0]=0,  this->NodalElementConnectivitySize[1]=0; this->NodalElementConnectivity=NULL;
}
/*}}}*/
BamgMesh::~BamgMesh(){/*{{{*/

	delete [] this->Vertices;
	delete [] this->PreviousNumbering;
	delete [] this->Edges;
	delete [] this->Triangles;

	delete [] this->SubDomains;
	delete [] this->SubDomainsFromGeom;

	delete [] this->VerticesOnGeomVertex;
	delete [] this->VerticesOnGeomEdge;
	delete [] this->EdgesOnGeomEdge;

	delete [] this->IssmEdges;
	delete [] this->IssmSegments;

	delete [] this->ElementConnectivity;
	delete [] this->NodalConnectivity;
	delete [] this->NodalElementConnectivity;
}
/*}}}*/
