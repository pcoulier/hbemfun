/*GAUSSPW1D Gaussian points and weights
 *  [xi,H]=gausspw1d(nEltDiv,nGauss) computes the integration point coordinates
 *  and weigths for an integration scheme with nGauss integration points.
 *
 *  nEltDiv Number of Element divisions
 *  nGauss  Number of Gaussian points.
 *  xi      Natural coordinates of the Gaussian points (1 * nGauss).
 *  H       Gaussian weights (1 * nGauss).
 */
 
/* $Make: mex -O -output gausspw1d gausspw1d_mex.cpp gausspw.cpp $*/


#include "mex.h"
#include <string.h>
#include "gausspw.h"
#include "checklicense.h"
#include <new>

using namespace std;

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{
  try
  {
    checklicense();
    
    if (nrhs<2) throw("Not enough input arguments.");
    if (nrhs>2) throw("Too many input arguments.");
    if (nlhs>2) throw("Too many output arguments.");
    
    if (!mxIsNumeric(prhs[0])) throw("Input argument 'nEltDiv' must be numeric.");
    if (mxIsSparse(prhs[0])) throw("Input argument 'nEltDiv' must not be sparse.");
    if (mxIsComplex(prhs[0])) throw("Input argument 'nEltDiv' must be real.");
    if (!(mxGetNumberOfElements(prhs[0])==1)) throw("Input argument 'nEltDiv' should be a scalar.");
    const double* const Inarg=mxGetPr(prhs[0]);
    int nEltDiv = int(Inarg[0]);
                
    if (!mxIsNumeric(prhs[1])) throw("Input argument 'nGauss' must be numeric.");
    if (mxIsSparse(prhs[1])) throw("Input argument 'nGauss' must not be sparse.");
    if (mxIsComplex(prhs[1])) throw("Input argument 'nGauss' must be real.");
    if (!(mxGetNumberOfElements(prhs[1])==1)) throw("Input argument 'nGauss' should be a scalar.");
    const double* const Inarg2=mxGetPr(prhs[1]);
    int nGauss = int(Inarg2[0]);

    double* const Xi=new(nothrow) double[nGauss*nEltDiv];
    double* const H =new(nothrow) double[nGauss*nEltDiv];
    
    gausspw1D(nEltDiv,nGauss,Xi,H);

    // Output arguments
    plhs[0]=mxCreateDoubleMatrix(1,nGauss*nEltDiv,mxREAL);
    double* const xiOut = mxGetPr(plhs[0]);
    for (int iXi=0; iXi<nGauss*nEltDiv; iXi++) xiOut[iXi]=Xi[iXi];
    if (nlhs>1)
    {
       plhs[1]=mxCreateDoubleMatrix(1,nGauss*nEltDiv,mxREAL);
       double* const HOut = mxGetPr(plhs[1]);
       for (int iXi=0; iXi<nGauss*nEltDiv; iXi++)  HOut[iXi]=H[iXi];
    }
    delete [] Xi;
    delete [] H;
  }
  catch (const char* exception)
  {
    mexErrMsgTxt(exception);
  }
}
