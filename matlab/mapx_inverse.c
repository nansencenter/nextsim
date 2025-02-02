#include "mex.h"
#include "matrix.h"
#include "mapx.h"
#include "unistd.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
        /*---------- Input  ----------*/

        /* Check for proper number of input and output arguments */
        if (nrhs < 3)
                mexErrMsgTxt("Insufficient inputs");
        if (nrhs > 3)
                mexErrMsgTxt("Too many input arguments");

        /* Check for proper input type */
        if (!mxIsChar(prhs[0]) || (mxGetM(prhs[0]) != 1 ) )
                mexErrMsgTxt("First input argument must be a string.");
        if (!mxIsDouble(prhs[1]) || (mxGetM(prhs[1]) != 1 ) )
                mexErrMsgTxt("Second input argument must be a double and a column vector.");
        if (!mxIsDouble(prhs[2]) || (mxGetM(prhs[2]) != 1 ) )
                mexErrMsgTxt("Third input argument must be a double and a column vector.");

        size_t buflen  = mxGetN(prhs[0])*sizeof(mxChar)+1;
        char *filename = mxMalloc(buflen);

        /* Copy the string data into filename. */ 
        int status = mxGetString(prhs[0], filename, (mwSize)buflen);  

        double *x = (double *) mxGetPr(prhs[1]);
        double *y = (double *) mxGetPr(prhs[2]);

        mwSize const Np = mxGetNumberOfElements(prhs[1]);

        /*---------- Output ----------*/

        mxArray *plhs0  = mxCreateDoubleMatrix(Np, 1, mxREAL);
        double  *lon    = mxGetPr(plhs0);
        mxArray *plhs1  = mxCreateDoubleMatrix(Np, 1, mxREAL);
        double  *lat    = mxGetPr(plhs1);

        /*---------- Code ----------*/

        mapx_class *map;

        /* Check if file exists */
        if ( access(filename, R_OK) == -1 )
                mexErrMsgTxt("MPP file not readable");

        map = init_mapx(filename);

	int i;
        for ( i=0; i<Np; i++ )
                inverse_mapx(map, x[i], y[i], &lat[i], &lon[i]);

        close_mapx(map);

        /* Export to matlab */
        plhs[0]=plhs0;
        plhs[1]=plhs1;
}

