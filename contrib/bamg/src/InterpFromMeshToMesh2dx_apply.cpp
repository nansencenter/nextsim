/*!\file InterpFromMeshToMesh2dx
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

int InterpFromMeshToMesh2dx_apply(double** pdata_interp,
      std::vector<std::vector<double>>& areacoord, std::vector<std::vector<int>>& vertex,
      std::vector<int> it,
      int nods_data,int nels_data,
      double* data,int M_data,int N_data,
      double* x_interp,double* y_interp)
{
    /* Set N_interp */
    int N_interp;
    if (M_data==nods_data)
        N_interp = areacoord.size();
    else
        N_interp = it.size();

    /*Output*/
    double* data_interp=NULL;

    /*Intermediary*/
    int    i,j;

    /*Checks*/
   if (M_data!=nods_data && M_data!=nels_data)
   {
        _error_("data provided should have either " << nods_data << " or " << nels_data << " lines (not " << M_data << ")");
    }

    /*Initialize output*/
    data_interp=xNew<double>(N_interp*N_data);

    /*Loop over output nodes*/
#ifndef BAMG_NO_OMP
    int thread_id;
    int max_threads = omp_get_max_threads(); /*8 by default on MACOSX (2,5 GHz Intel Core i7)*/
#pragma omp parallel for num_threads(max_threads) private(thread_id)
#endif
    for (int i=0; i < N_interp; ++i)
    {
        if (M_data==nods_data)
        {
            for (j=0;j<N_data;j++)
            {
                data_interp[i*N_data+j]=areacoord[i][0]*data[N_data*vertex[i][0]+j]
                    +areacoord[i][1]*data[N_data*vertex[i][1]+j]
                    +areacoord[i][2]*data[N_data*vertex[i][2]+j];
            }
        }
        else{
            for (j=0;j<N_data;j++)
            {
                if (it[i]<0 || it[i]>=nels_data)
                {
                    _error_("Triangle number " << it[i] << " not in [0 " << nels_data << "], report bug to developers");
                }
                data_interp[i*N_data+j]=data[N_data*it[i]+j];
            }
        }
    }

    /*clean-up and return*/
    *pdata_interp=data_interp;
    return 1;
}
