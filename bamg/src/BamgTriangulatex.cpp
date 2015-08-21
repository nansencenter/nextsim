/*!\file BamgTriangulatex
 */

#include "./BamgTriangulatex.h"
#include "./Mesh.h"

using namespace bamg;
using namespace std;

int BamgTriangulatex(int** pindex,int* pnels,double* x,double* y,int nods){

	Mesh Th(x,y,nods);
	Th.WriteIndex(pindex,pnels);
	//delete &Th;
	return 0;

}
