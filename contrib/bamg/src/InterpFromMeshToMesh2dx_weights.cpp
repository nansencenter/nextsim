/*!\file InterpFromMeshToMesh2dx_weights
 */

#include "./InterpFromMeshToMesh2dx.h"

#include "./Mesh.h"
#include "./shared.h"
#include "./typedefs.h"
//#include "../../toolkits/toolkits.h"
//#include "../../classes/classes.h"
//#include "../../bamg/bamgobjects.h"

using namespace bamg;
using namespace std;

int InterpFromMeshToMesh2dx_weights(
      std::vector<std::vector<double>>& areacoord, std::vector<std::vector<int>>& vertex,
      std::vector<int> it,
      int* index_data,double* x_data,double* y_data,
      int nods_data,int nels_data,
      int M_data,
      const double* x_interp, const double* y_interp,int N_interp)
{
    /*Output*/
    if (M_data==nods_data)
    {
        areacoord.resize(N_interp);
        vertex.resize(N_interp);
    } else {
        it.resize(N_interp);
    }

    /*Intermediary*/
    R2     r;
    I2     I;
    int    i;
    Icoor2 dete[3];

    /*Checks*/
   if (M_data!=nods_data && M_data!=nels_data)
   {
        _error_("data provided should have either " << nods_data << " or " << nels_data << " lines (not " << M_data << ")");
    }

    /*read background mesh*/
    Mesh* Th=new Mesh(index_data,x_data,y_data,nods_data,nels_data);

    /*Get reference number (for subdomains)*/
    long* reft = xNew<long>(Th->nbt);
    Th->TriangleReferenceList(reft);
    Th->CreateSingleVertexToTriangleConnectivity();

    /*Loop over output points*/
#ifndef BAMG_NO_OMP
    int thread_id;
    int max_threads = omp_get_max_threads(); /*8 by default on MACOSX (2,5 GHz Intel Core i7)*/
#pragma omp parallel for num_threads(max_threads) private(thread_id)
#endif
    for (int i=0; i < N_interp; ++i)
    {
        /*Get current point coordinates*/
        r.x=x_interp[i]; r.y=y_interp[i];
        I2 I=Th->R2ToI2(r);

        /*Find triangle holding r/I*/
        Triangle &tb=*Th->TriangleFindFromCoord(I,dete);

        /*point inside convex*/
        if (tb.det>0)
        {
            int it_tmp=Th->GetId(tb);
            if (reft[it_tmp]<0)
                throw std::runtime_error("InterpFrmoMeshToMesh2dx_weights: Point outside of mesh (but inside convex)\n");

            if (M_data==nods_data)
            {
                /*Area coordinates*/
                areacoord[i].resize(3);
                areacoord[i][0]= (double) dete[0]/tb.det;
                areacoord[i][1]= (double) dete[1]/tb.det;
                areacoord[i][2]= (double) dete[2]/tb.det;
                /*3 vertices of the triangle*/
                vertex[i].resize(3);
                vertex[i][0]=Th->GetId(tb[0]);
                vertex[i][1]=Th->GetId(tb[1]);
                vertex[i][2]=Th->GetId(tb[2]);
            }
            else
            {
                /*triangle number*/
                it[i]=it_tmp;
            }
        }
        //external point
        else{
            //Get closest adjacent triangle (inside the mesh)
            double aa,bb;
            AdjacentTriangle ta=CloseBoundaryEdge(I,&tb,aa,bb).Adj();
            int k=ta;
            Triangle &tc=*(Triangle*)ta;
            if (M_data==nods_data)
            {
                //Area coordinate
                areacoord[i].resize(3);
                areacoord[i][VerticesOfTriangularEdge[k][1]] = aa;
                areacoord[i][VerticesOfTriangularEdge[k][0]] = bb;
                areacoord[i][OppositeVertex[k]] = 1 - aa -bb;
                //3 vertices of the triangle
                vertex[i].resize(3);
                vertex[i][0]=Th->GetId(tc[0]);
                vertex[i][1]=Th->GetId(tc[1]);
                vertex[i][2]=Th->GetId(tc[2]);
            }
            else
            {
                //triangle number
                it[i]=Th->GetId(tc);
            }
        }
    }

    /*clean-up and return*/
    delete Th;
    xDelete<long>(reft);
    return 1;
}

