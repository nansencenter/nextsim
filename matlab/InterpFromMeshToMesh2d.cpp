#include "mex.h"
#include "matrix.h"
#include "InterpFromMeshToMesh2dmex.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
        /*---------- Input  ----------*/

        /* Check for proper number of input and output arguments */
        if (nrhs < 8) mexErrMsgTxt("Insufficient inputs");
        if (nrhs > 8) mexErrMsgTxt("Too many input arguments");

        int    *index_data    = (int *)    mxGetPr(prhs[0]);
        double *x_data        = (double *) mxGetPr(prhs[1]);
        double *y_data        = (double *) mxGetPr(prhs[2]);
        double *data          = (double *) mxGetPr(prhs[3]);
        double *x_interp      = (double *) mxGetPr(prhs[4]);
        double *y_interp      = (double *) mxGetPr(prhs[5]);
        bool   isdefault      = (bool)     mxGetScalar(prhs[6]);
        double defaultvalue   = (double)   mxGetScalar(prhs[7]);

        int nels_data = mxGetN(prhs[0]);
        int nods_data = mxGetM(prhs[1]);
        int M_data    = mxGetN(prhs[3]);
        int N_data    = mxGetM(prhs[3]);
        int N_interp  = mxGetN(prhs[4]);

        /*---------- Output ----------*/

        plhs[0]      = mxCreateDoubleMatrix(N_data, N_interp, mxREAL);
        double  *interp_out;
        double  *output = mxGetPr(plhs[0]);

        /*---------- Code ----------*/

        // Error checks
        if ( ! mxIsInt32(prhs[0]) ) mexErrMsgTxt("index must be int32(index)");
        if ( mxGetM(prhs[0]) != 3 ) mexErrMsgTxt("index should have 3 columns");
        if ( nods_data < 3 ) mexErrMsgTxt("there should be at least three points");
        if ( mxGetM(prhs[2]) != nods_data ) mexErrMsgTxt("vectors x and y should have the same length");
        if ( M_data*N_data < 1 ) mexErrMsgTxt("data is empty");
        if ( M_data != nods_data & M_data != nels_data ) mexErrMsgTxt("data has wrong size");
        if ( N_interp < 1 ) mexErrMsgTxt("no interpolation requested");
        if ( mxGetN(prhs[5]) != N_interp ) mexErrMsgTxt("vectors x_interp and y_interp should have the same length");

        // Function call
        InterpFromMeshToMesh2dmex(&interp_out,
                            index_data, x_data, y_data,
                            nods_data, nels_data,
                            data,
                            M_data, N_data,
                            x_interp, y_interp, N_interp,
                            isdefault, defaultvalue);

        // Copy things to MATLAB
        // Couldn't figure out how to do it without copying ... but it should be possible!
        for (int i=0; i<N_interp*N_data; ++i)
            output[i] = interp_out[i];
}
