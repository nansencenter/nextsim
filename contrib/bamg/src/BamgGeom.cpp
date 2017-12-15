#include "./include.h"
#include "./shared.h"

/*Constructors/Destructors*/
BamgGeom::BamgGeom(){/*{{{*/

	this->VerticesSize[0]=0,          this->VerticesSize[1]=0;          this->Vertices=NULL;
	this->EdgesSize[0]=0,             this->EdgesSize[1]=0;             this->Edges=NULL;
	this->TangentAtEdgesSize[0]=0,    this->TangentAtEdgesSize[1]=0;    this->TangentAtEdges=NULL;
	this->CornersSize[0]=0,           this->CornersSize[1]=0;           this->Corners=NULL;
	this->RequiredVerticesSize[0]=0,  this->RequiredVerticesSize[1]=0;  this->RequiredVertices=NULL;
	this->RequiredEdgesSize[0]=0,     this->RequiredEdgesSize[1]=0;     this->RequiredEdges=NULL;
	this->SubDomainsSize[0]=0,        this->SubDomainsSize[1]=0;        this->SubDomains=NULL;

}
/*}}}*/
BamgGeom::~BamgGeom(){/*{{{*/

	delete [] this->Vertices;
	delete [] this->Edges;
	delete [] this->TangentAtEdges;
	delete [] this->Corners;
	delete [] this->RequiredVertices;
	delete [] this->RequiredEdges;
	delete [] this->SubDomains;
}
/*}}}*/

void
BamgGeom::reset()
{
	this->VerticesSize[0]=0,          this->VerticesSize[1]=0;
	this->EdgesSize[0]=0,             this->EdgesSize[1]=0;
	this->TangentAtEdgesSize[0]=0,    this->TangentAtEdgesSize[1]=0;
	this->CornersSize[0]=0,           this->CornersSize[1]=0;
	this->RequiredVerticesSize[0]=0,  this->RequiredVerticesSize[1]=0;
	this->RequiredEdgesSize[0]=0,     this->RequiredEdgesSize[1]=0;
	this->SubDomainsSize[0]=0,        this->SubDomainsSize[1]=0;

	delete [] this->Vertices;
	delete [] this->Edges;
	delete [] this->TangentAtEdges;
	delete [] this->Corners;
	delete [] this->RequiredVertices;
	delete [] this->RequiredEdges;
	delete [] this->SubDomains;
}
