// nsra_mex_omp.cpp
#include "mex.h"
#include "nsra.cpp"  // include the optimized core

/* How to compile:
    >> cd cpp
    >> mex nsra_mex.cpp
*/

void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[]) {

    if (nrhs < 2)
        mexErrMsgIdAndTxt("NSRA:InvalidInput","At least meas and pred required");

    size_t n = mxGetNumberOfElements(prhs[0]);
    double* meas = mxGetPr(prhs[0]);
    double* pred = mxGetPr(prhs[1]);

    if (mxGetNumberOfElements(prhs[1]) != n)
        mexErrMsgIdAndTxt("NSRA:InvalidInput","meas and pred must have same length");

    double epsilon = 0.0;
    double tie_score = 0.5;
    if (nrhs >= 3) epsilon = mxGetScalar(prhs[2]);
    if (nrhs >= 4) tie_score = mxGetScalar(prhs[3]);

    double score = nsra_core_omp(meas, pred, n, epsilon, tie_score);

    plhs[0] = mxCreateDoubleScalar(score);
}
