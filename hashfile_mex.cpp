/*HASHFILE   Compute RMD-128 hash code for a file.
 *
 *   h = HASHFILE(file) computes the RMD-128 hash code for the specified file.
 *
 *   file  Name of the file.
 *   h     Hash code.

 * Mattias Schevenels
 * June 2008
 */

/* $Make: mex -largeArrayDims -O -output hashfile hashfile_mex.cpp ripemd128.cpp $ */

#include "mex.h"
#include "ripemd128.h"

using namespace std;

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{
  try
  {
    if (nrhs<1) throw("Not enough input arguments.");
    if (nrhs>1) throw("Too many input arguments.");
    if (nlhs>1) throw("Too many output arguments.");

    if (!mxIsChar(prhs[0])) throw("Input argument must be a char array.");

    const int n = mxGetNumberOfElements(prhs[0])+1;
    char* file = new(nothrow) char[n];
    if (file==0) throw("Out of memory.");
    mxGetString(prhs[0],file,n);

    string hash;

    rmdfile(file,hash);

    plhs[0] = mxCreateString(hash.c_str());

    delete [] file;
  }
  catch (const char* exception)
  {
    mexErrMsgTxt(exception);
  }
}
