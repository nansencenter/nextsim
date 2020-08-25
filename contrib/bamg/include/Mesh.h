#ifndef _MESH_H_
#define _MESH_H_

#include "./BamgOpts.h"
#include "./BamgMesh.h"
#include "./BamgGeom.h"
#include "./Triangle.h"
#include "./VertexOnGeom.h"
#include "./VertexOnVertex.h"
#include "./VertexOnEdge.h"
#include "./ListofIntersectionTriangles.h"
#include "./typedefs.h"
#include "./macros.h"

namespace bamg {

	class Geometry;
	class BamgQuadtree;
	class SubDomain;

	class Mesh {

		public:

			Geometry                    & Gh;                    // Geometry
			Mesh                        & BTh;                   // Background Mesh Bth== *this =>no background
			BamgVertex                   *vertices;
			Triangle                     *triangles;
			Edge                         *edges;
			BamgQuadtree                 *quadtree;
			BamgVertex                  **orderedvertices;
			SubDomain                    *subdomains;
			long                          NbRef;                 // counter of ref on the this class if 0 we can delete
			long                          maxnbv,maxnbt;         // nombre max de sommets , de triangles
			long                          nbv,nbt,nbe;       // nb of vertices, of triangles, of edges
			long                          nbsubdomains;
			long                          nbtout;                // Nb of oudeside triangle

			R2                            pmin,pmax;             // extrema
			double                        coefIcoor;             // coef to integer Icoor1;
			ListofIntersectionTriangles   lIntTria;

			long                          NbVerticesOnGeomVertex;
			VertexOnGeom                 *VerticesOnGeomVertex;
			long                          NbVerticesOnGeomEdge;
			VertexOnGeom                 *VerticesOnGeomEdge;
			long                          NbVertexOnBThVertex;
			VertexOnVertex               *VertexOnBThVertex;
			long                          NbVertexOnBThEdge;
			VertexOnEdge                 *VertexOnBThEdge;

			//Constructors/Destructors
			Mesh(BamgGeom* bamggeom,BamgMesh* bamgmesh,BamgOpts* bamgopts);
			Mesh(int* index,double* x,double* y,int nods,int nels);/*MeshConvert*/
			Mesh(double* x,double* y,int nods); /*BamgTriangulate*/
			Mesh(Mesh &,Geometry * pGh=0,Mesh* pBTh=0,long maxnbv_in=0 ); //copy operator
			Mesh(const Mesh &,const int *flag,const int *bb,BamgOpts* bamgopts); // truncature
			Mesh(long maxnbv,Mesh & BT,BamgOpts* bamgopts,int keepBackVertices=1);
			Mesh(long maxnbv,Geometry & G,BamgOpts* bamgopts);
			~Mesh();

			//Operators
			const BamgVertex &operator[](unsigned long i) const { return vertices[i];  };
			BamgVertex       &operator[](unsigned long i) { return vertices[i];        };
			const Triangle   &operator()(unsigned long i) const { return triangles[i]; };
			Triangle         &operator()(unsigned long i) { return triangles[i];       };

			//Methods
			void SetIntCoor(const char * from =0);
			double MinimalHmin();
			double MaximalHmax();
			I2 R2ToI2(const R2 & P) const;
			R2 I2ToR2(const I2 & P) const;
			void AddVertex(BamgVertex & s,Triangle * t,Icoor2 *  =0) ;
			void Insert(bool random);
			void Echo(void);
			void ForceBoundary();
			void FindSubDomain(int OutSide=0);
			long TriangleReferenceList(long*) const;
			void TriangleIntNumbering(long* renumbering);
			void SmoothMetric(double raisonmax) ;
			void BoundAnisotropy(double anisomax,double hminaniso= 1e-100) ;
			void MaxSubDivision(double maxsubdiv);
			Edge** MakeGeomEdgeToEdge();
			long SplitInternalEdgeWithBorderVertices();
			void MakeBamgQuadtree();
			void NewPoints(Mesh &,BamgOpts* bamgopts,int KeepVertices=1);
			long InsertNewPoints(long nbvold,long & NbTSwap,bool random);
			void TrianglesRenumberBySubDomain(bool justcompress=false);
			Metric MetricAt (const R2 &) const;
			void SmoothingVertex(int =3,double=0.3);
			GeomEdge* ProjectOnCurve( Edge & AB, BamgVertex &  A, BamgVertex & B,double theta, BamgVertex & R,VertexOnEdge & BR,VertexOnGeom & GR);
			long GetId(const Triangle & t) const;
			long GetId(const Triangle * t) const;
			long GetId(const BamgVertex & t) const;
			long GetId(const BamgVertex * t) const;
			long GetId(const Edge & t) const;
			long GetId(const Edge * t) const;
			BamgVertex* NearestVertex(Icoor1 i,Icoor1 j) ;
			Triangle* TriangleFindFromCoord(const I2 & ,Icoor2 [3],Triangle *tstart=0) const;
			void ReadMesh(int* index,double* x,double* y,int nods,int nels);
			void ReadMesh(BamgMesh* bamgmesh, BamgOpts* bamgopts);
			void WriteMesh(BamgMesh* bamgmesh,BamgOpts* bamgopts);
			void ReadMetric(const BamgOpts* bamgopts);
			void WriteMetric(BamgOpts* bamgopts);
			void WriteIndex(int** pindex,int* pnels);
			void AddMetric(BamgOpts* bamgopts);
			void BuildMetric0(BamgOpts* bamgopts);
			void BuildMetric1(BamgOpts* bamgopts);
			void AddGeometryMetric(BamgOpts* bamgopts);
			void BuildGeometryFromMesh(BamgOpts* bamgopts=NULL);
			void ReconstructExistingMesh();

			void CreateSingleVertexToTriangleConnectivity();
			void UnMarkUnSwapTriangle();
			void SetVertexFieldOn();
			void SetVertexFieldOnBTh();

		private:
			void TriangulateFromGeom1(BamgOpts* bamgopts,int KeepVertices=1);// the real constructor mesh adaption
			void TriangulateFromGeom0(BamgOpts* bamgopts);// the real constructor mesh generator
			void Triangulate(double* x,double* y,int nods);
			void Init(long);
	};

	/*Intermediary*/
	AdjacentTriangle CloseBoundaryEdge(I2 ,Triangle *, double &,double &) ;
	void  swap(Triangle *t1,short a1,
				Triangle *t2,short a2,
				BamgVertex *s1,BamgVertex *s2,Icoor2 det1,Icoor2 det2);
	int SwapForForcingEdge(BamgVertex   *  & pva ,BamgVertex  * &   pvb ,
				AdjacentTriangle & tt1,Icoor2 & dets1,
				Icoor2 & detsa,Icoor2 & detsb, int & nbswap);
	int ForceEdge(BamgVertex &a, BamgVertex & b,AdjacentTriangle & taret) ;
	inline AdjacentTriangle Previous(const AdjacentTriangle & ta){
		return AdjacentTriangle(ta.t,PreviousEdge[ta.a]);
	}
	inline AdjacentTriangle Next(const AdjacentTriangle & ta){
		return AdjacentTriangle(ta.t,NextEdge[ta.a]);
	}
	inline  AdjacentTriangle Adj(const AdjacentTriangle & a){
		return  a.Adj();
	}
	inline void Adj(GeomEdge * & on,int &i){
		int j=i;i=on->AdjVertexIndex[i];on=on->Adj[j];
	}
	inline double qualite(const BamgVertex &va,const BamgVertex &vb,const BamgVertex &vc){
		double ret;
		I2 ia=va,ib=vb,ic=vc;
		I2 ab=ib-ia,bc=ic-ib,ac=ic-ia;
		Icoor2 deta=Det(ab,ac);
		if (deta <=0) ret = -1;
		else {
			double a = sqrt((double) (ac,ac)),
					 b = sqrt((double) (bc,bc)),
					 c = sqrt((double) (ab,ab)),
					 p = a+b+c;
			double h= Max(Max(a,b),c),ro=deta/p;
			ret = ro/h;
		}
		return ret;
	}

}
#endif
