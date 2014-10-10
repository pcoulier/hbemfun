/*BEMFUNLICENSE   Verify BEMFUN license.
 *
 *   BEMFUNLICENSE verifies the BEMFUN license.
 *
 *   BEMFUNLICENSE VERIFYONCE does the same, but only if the license has not been
 *   verified during the last 2 hours.  In order to determine whether the
 *   license has already been verified within this period, the verification
 *   state is stored and the function BEMFUNLICENSE is locked in memory.
 *
 *   BEMFUNLICENSE RESET unlocks the function BEMFUNLICENSE and clears the
 *   verification state.

 * Mattias Schevenels
 * July 2009
 */

/* $Make: mex GCCLIBS="$GCCLIBS -liphlpapi" -largeArrayDims -O -output bemfunlicense bemfunlicense_mex.cpp bemfunlicense.cpp checklicense.cpp bigint.cpp ripemd128.cpp rsa.cpp getmac.cpp $ */

#include "mex.h"
#include "bemfunlicense.h"
#include <string>

using namespace std;

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{
  try
  {
    if (nrhs>1) throw("Too many input arguments.");
    if (nlhs>1) throw("Too many output arguments.");

    char action[32];
    if (nrhs>=1) mxGetString(prhs[0],action,32); else action[0] = '\0';
    action[31] = '\0';

    int result = bemfunlicense(action);

    if (nlhs>0) plhs[0] = mxCreateDoubleScalar(result);

  }
  catch (const char* exception)
  {
    mexErrMsgTxt(exception);
  }
}
