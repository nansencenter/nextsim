/*!\file OptionUtilities.cpp
 * \brief: implementation of the options utilities
 */

/*Headers:*/
/*{{{*/
//#include "../classes.h"
//#include "include.h"
#include "./OptionUtilities.h"
#include "./shared.h"
#include "recast.h"
#include "./exceptions.h"
/*}}}*/

int ColumnWiseDimsFromIndex(int* dims,int index,int* size,int ndims){/*{{{*/

	int   i;
	int   aprod=1;

	/*check for index too large  */
	for (i=0;i<ndims;i++) aprod*=size[i];
	if (index >= aprod) _error_("Index " << index << " exceeds number of elements " << aprod << ".");

	/*calculate the dimensions (being careful of integer division)  */
	for (i=ndims-1; i>=0; i--) {
		aprod=reCast<int>(((double)aprod+0.5)/(double)size[i]);
		dims[i]=(int)floor(((double)index+0.5)/(double)aprod);
		index-=dims[i]*aprod;
	}

	return(0);
}/*}}}*/
int IndexFromColumnWiseDims(int* dims, int* size, int ndims) {/*{{{*/

	int   i;
	int   index=0;

	/*check for any dimension too large  */
	for (i=0;i<ndims;i++){
		if (dims[i] >= size[i]) _error_("Dimension " << i << " of " << dims[i] << " exceeds size of " << size[i] << ".");
	}

	/*calculate the index  */
	for (i=ndims-1; i>=0; i--){
		index*=size[i];
		index+=dims[i];
	}

	return(index);
}/*}}}*/
int RowWiseDimsFromIndex(int* dims, int index, int* size, int ndims) {/*{{{*/

	int   i;
	int   aprod=1;

	/*check for index too large  */
	for (i=0; i<ndims; i++) aprod*=size[i];
	if (index >= aprod) _error_("Index " << index << " exceeds number of elements " << aprod << ".");

	/*calculate the dimensions (being careful of integer division)  */
	for (i=0; i<ndims; i++) {
		aprod=(int)(((double)aprod+0.5)/(double)size[i]);
		dims[i]=(int)floor(((double)index+0.5)/(double)aprod);
		index-=dims[i]*aprod;
	}

	return(0);
}/*}}}*/
int IndexFromRowWiseDims(int* dims, int* size, int ndims) {/*{{{*/

	int   i;
	int   index=0;

	/*check for any dimension too large  */
	for (i=0; i<ndims; i++){
		if (dims[i] >= size[i]) _error_("Dimension " << i << " of " << dims[i] << " exceeds size of " << size[i] << ".");
	}

	/*calculate the index  */
	for (i=0; i<ndims; i++) {
		index*=size[i];
		index+=dims[i];
	}

	return(index);
}/*}}}*/
int StringFromDims(char* cstr, int* dims, int ndims) {/*{{{*/

	sprintf(&cstr[0],"[");
	for(int i=0; i<ndims-1; i++) sprintf(&cstr[strlen(cstr)],"%d,",dims[i]);
	sprintf(&cstr[strlen(cstr)],"%d]",dims[ndims-1]);

	return(0);
}/*}}}*/
int StringFromSize(char* cstr, int* size, int ndims) {/*{{{*/

	sprintf(&cstr[0],"[");
	for(int i=0; i<ndims-1; i++) sprintf(&cstr[strlen(cstr)],"%dx",size[i]);
	sprintf(&cstr[strlen(cstr)],"%d]",size[ndims-1]);

	return(0);
}/*}}}*/
