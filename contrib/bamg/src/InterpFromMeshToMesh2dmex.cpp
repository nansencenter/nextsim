/*!\file InterpFromMeshToMesh2dx
 */

#include "./InterpFromMeshToMesh2dmex.h"

#include "./Mesh.h"
#include "./shared.h"
#include "./typedefs.h"
//#include "../../toolkits/toolkits.h"
//#include "../../classes/classes.h"
//#include "../../bamg/bamgobjects.h"

using namespace bamg;
using namespace std;

int InterpFromMeshToMesh2dmex(double** pdata_interp,int* index_data,double* x_data,double* y_data,int nods_data,int nels_data,
			double* data,int M_data,int N_data,double* x_interp,double* y_interp,int N_interp,bool isdefault, double defaultvalue /*Options* options*/){

	/*Output*/
	double* data_interp=NULL;

	/*Intermediary*/
	double xmin,xmax,ymin,ymax;
	//bool   isdefault;
	//double defaultvalue;
	R2     r;
	I2     I;
	int    i,j;
	int    it;
	int    i0,i1,i2;
	double areacoord[3];
	double aa,bb;
	Icoor2 dete[3];

	/*Checks*/
	if (M_data!=nods_data && M_data!=nels_data){
		_error_("data provided should have either " << nods_data << " or " << nels_data << " lines (not " << M_data << ")");
	}

#if 0
	/*Get default*/
	if(options->GetOption("default")){
		isdefault=true;
		options->Get(&defaultvalue,"default");
	}
	else{
		isdefault=false;
	}
#endif

	/*Initialize output*/
	data_interp=xNew<double>(N_interp*N_data);

	/*read background mesh*/
	Mesh* Th=new Mesh(index_data,x_data,y_data,nods_data,nels_data);

	/*Get reference number (for subdomains)*/
	long* reft = xNew<long>(Th->nbt);
	Th->TriangleReferenceList(reft);
	Th->CreateSingleVertexToTriangleConnectivity();

	/*Get domain boundaries*/
	xmin=x_data[0]; ymin=y_data[0];
	xmax=x_data[0]; ymax=y_data[0];
	for(i=1;i<nods_data;i++){
		if(x_data[i]<xmin) xmin=x_data[i];
		if(x_data[i]>xmax) xmax=x_data[i];
		if(y_data[i]<ymin) ymin=y_data[i];
		if(y_data[i]>ymax) ymax=y_data[i];
	}

	/*Loop over output nodes*/
	//for(i=0;i<N_interp;i++){
    int thread_id;

    for (int i=0; i < N_interp; ++i)
    {
    
    	//if(i%100==0) _printf_("\r      interpolation progress: "<<setw(6)<<setprecision(2)<<double(i)/double(N_interp)*100.<<"%   ");

		if(isdefault){
			if(x_interp[i]<xmin || x_interp[i]>xmax || y_interp[i]<ymin || y_interp[i]>ymax){
				for(j=0;j<N_data;j++) data_interp[i*N_data+j]=defaultvalue;
				continue;
			}
		}

		/*Get current point coordinates*/
		r.x=x_interp[i]; r.y=y_interp[i];
		I2 I=Th->R2ToI2(r);

		/*Find triangle holding r/I*/
		Triangle &tb=*Th->TriangleFindFromCoord(I,dete);

		/*point inside convex*/
		if (tb.det>0){

			/*Area coordinates*/
			areacoord[0]= (double) dete[0]/tb.det;
			areacoord[1]= (double) dete[1]/tb.det;
			areacoord[2]= (double) dete[2]/tb.det;
			/*3 vertices of the triangle*/
			i0=Th->GetId(tb[0]);
			i1=Th->GetId(tb[1]);
			i2=Th->GetId(tb[2]);
			/*triangle number*/
			it=Th->GetId(tb);

			/*Inside convex but outside mesh*/
			if (reft[it]<0 && isdefault){
				for(j=0;j<N_data;j++) data_interp[i*N_data+j]=defaultvalue;
				continue;
			}
		}
		//external point
		else{
			if(isdefault){
				for(j=0;j<N_data;j++) data_interp[i*N_data+j]=defaultvalue;
				continue;
			}
			else{
				//Get closest adjacent triangle (inside the mesh)
				AdjacentTriangle ta=CloseBoundaryEdge(I,&tb,aa,bb).Adj();
				int k=ta;
				Triangle &tc=*(Triangle*)ta;
				//Area coordinate
				areacoord[VerticesOfTriangularEdge[k][1]] = aa;
				areacoord[VerticesOfTriangularEdge[k][0]] = bb;
				areacoord[OppositeVertex[k]] = 1 - aa -bb;
				//3 vertices of the triangle
				i0=Th->GetId(tc[0]);
				i1=Th->GetId(tc[1]);
				i2=Th->GetId(tc[2]);
				//triangle number
				it=Th->GetId(tc);
			}
		}

		if (M_data==nods_data){
			for (j=0;j<N_data;j++){
				data_interp[i*N_data+j]=areacoord[0]*data[N_data*i0+j]+areacoord[1]*data[N_data*i1+j]+areacoord[2]*data[N_data*i2+j];

			}
		}
		else{
			for (j=0;j<N_data;j++){
				if (it<0 || it>=nels_data){
					_error_("Triangle number " << it << " not in [0 " << nels_data << "], report bug to developers");
				}
				data_interp[i*N_data+j]=data[N_data*it+j];
			}
		}
	}
	//if(N_interp>=100) _printf_("\r      interpolation progress: "<<fixed<<setw(6)<<setprecision(2)<<100.<<"%  \n");

	/*clean-up and return*/
	delete Th;
	xDelete<long>(reft);
	*pdata_interp=data_interp;
	return 1;
}
