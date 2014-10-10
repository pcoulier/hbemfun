/*BEMMATCONV   Time domain boundary element matrix convolution.
 *
 *   B = BEMMATCONV(A,X) computes the matrix convolution of the
 *   matrices A and X, defined as:
 *
 *          N
 *   B  =  sum A(:,:,n) * X(:,N-n)
 *         n=1
 *
 *   A   Influence matrices (nRow * nDof * nTime).
 *   X   Stress or displacement history (nDof * nTime).
 *   B   Convolution (nRow * 1).
 */
 
/* $Make: mex -O -output bemmatconv bemmatconv_mex.cpp$*/ 

#include "mex.h"
#include "math.h"
#include <string.h>
#include "checklicense.h"

void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
  try
  {    
    checklicense();
    
    if (nrhs!=2) mexErrMsgTxt("Two input arguments required.");
    if (nlhs>1)  mexErrMsgTxt("Only one output argument allowed.");
    
    
    if (!mxIsNumeric(prhs[0])) throw("Input argument 'A' must be numeric.");
    if (mxIsSparse(prhs[0])) throw("Input argument 'A' must not be sparse.");
    if (mxIsComplex(prhs[0])) throw("Input argument 'A' must be real.");
    const size_t nDimA=mxGetNumberOfDimensions(prhs[0]);
    const size_t* const dimA=mxGetDimensions(prhs[0]);
    if (nDimA>3) throw("Input argument 'A' must have 3 dimensions at most.");
    int nRow=dimA[0];
    int nDof=dimA[1];
    int nTime;
    if (nDimA == 3) nTime=dimA[2]; else nTime=1;
    const double* const A = mxGetPr(prhs[0]);
    
    if (!mxIsNumeric(prhs[1])) throw("Input argument 'X' must be numeric.");
    if (mxIsSparse(prhs[1])) throw("Input argument 'X' must not be sparse.");
    if (mxIsComplex(prhs[1])) throw("Input argument 'X' must be real.");
    const int nDimX=mxGetNumberOfDimensions(prhs[1]);
    const size_t* const dimX=mxGetDimensions(prhs[1]);
    if (nDimX>2) throw("Input argument 'X' must have 2 dimensions at most.");
    if (dimX[0]!=(unsigned int)nDof) throw("Number of elements in input arguments 'A' and 'X' does not correspond.");
    if (dimX[1]!=(unsigned int)nTime) throw("Number of elements in input arguments 'A' and 'X' does not correspond.");
    const double* const X = mxGetPr(prhs[1]);
    
    // OUTPUT ARGUMENT AND POINTER TO OUTPUT ARGUMENT
    plhs[0] = mxCreateDoubleMatrix(nRow,1,mxREAL);
    double* const B = mxGetPr(plhs[0]);
    for (int iElt=0; iElt<nRow; iElt++) B[iElt]=0.0;
    
    for (int iTime=0; iTime<nTime; iTime++)
    {
      for(int iRow=0; iRow<nRow; iRow++)
      {
        for(int jDof=0; jDof<nDof; jDof++)
        {
          B[iRow]=B[iRow]+A[iRow+nRow*jDof+nRow*nDof*iTime]*X[jDof+nDof*(nTime-iTime-1)];
        }
      }
    }
  }
  catch (const char* exception)
  {
    mexErrMsgTxt(exception);
  }
}
