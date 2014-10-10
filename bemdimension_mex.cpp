/*BEMDIMENSION   Boundary element problem dimension.
 *
 *   dim=BEMDIMENSION(elt,typ) returns the boundary element problem dimension
 *   for the boundary element mesh consisting of elements Typ. Returns 2 for
 *   2D/2.5D, and 3 for 3D/axisymmetric.
 *
 *   elt  Elements.
 *   typ  Element types.
 *   dim  Element dimension.
 */

/* $Make: mex -O -output bemdimension bemdimension_mex.cpp bemdimension.cpp eltdef.cpp$*/


#include "mex.h"
#include <string.h>
#include "bemdimension.h" 
#include <new>
#include "checklicense.h"

using namespace std;

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{
  try
  {
    checklicense();
   
    if (nrhs<2) throw("Not enough input arguments.");
    if (nrhs>2) throw("Too many input arguments.");
    if (nlhs>2) throw("Too many output arguments.");

    if (!mxIsNumeric(prhs[0])) throw("Input argument 'elt' must be numeric.");
    if (mxIsSparse(prhs[0])) throw("Input argument 'elt' must not be sparse.");
    if (mxIsComplex(prhs[0])) throw("Input argument 'elt' must be real.");
    if (mxGetN(prhs[0])<=2) throw("Input argument 'elt' should have at least 2 columns.");
    const double* const Elt=mxGetPr(prhs[0]);
    const unsigned int nElt=mxGetM(prhs[0]);

    bool keyOpts=true;
    if (mxGetN(prhs[1])==3) keyOpts=true;
    else if  (mxGetN(prhs[1])==2) keyOpts=false;
    else throw("Input argument 'typ' should have 2 or 3 columns.");
    if (!(mxIsCell(prhs[1]))) throw("Input argument 'typ' should be a cell array.");
    const unsigned int nEltType=mxGetM(prhs[1]);
    const unsigned int maxKeyOpts = 50;  // Maximum number of keyoptions per element type
    unsigned int* const TypeID=new(nothrow) unsigned int[nEltType];
      if (TypeID==0) throw("Out of memory.");
    unsigned int* const nKeyOpt=new(nothrow) unsigned int[nEltType];
      if (nKeyOpt==0) throw("Out of memory.");
    char** const TypeName=new(nothrow) char*[nEltType];
      if (TypeName==0) throw("Out of memory.");
    char** const TypeKeyOpts=new(nothrow) char*[nEltType*maxKeyOpts];
      if (TypeKeyOpts==0) throw("Out of memory.");
    for (unsigned int iTyp=0; iTyp<nEltType; iTyp++)
    { 
      // TypeID
      const mxArray* TypPtr0=mxGetCell(prhs[1],iTyp+nEltType*0);
      if (!mxIsNumeric(TypPtr0)) throw("Type ID should be numeric.");
      if (mxIsSparse(TypPtr0)) throw("Type ID should not be sparse.");
      if (mxIsComplex(TypPtr0)) throw("Type ID should not be complex.");
      if (!(mxGetNumberOfElements(TypPtr0)==1)) throw("Type ID should be a scalar.");
      TypeID[iTyp]= (unsigned int)(mxGetScalar(TypPtr0));
      
      // TypeName
      const mxArray* TypPtr1=mxGetCell(prhs[1],iTyp+nEltType*1);
      if (!mxIsChar(TypPtr1)) throw("Element types should be input as stings.");
      TypeName[iTyp] =  mxArrayToString(TypPtr1);
      
      // TypeKeyOpts
      if (keyOpts)
      {
        const mxArray* TypPtr2=mxGetCell(prhs[1],iTyp+nEltType*2); // Keyoptions cell array
        if (!mxIsCell(TypPtr2)) throw("Keyopts should be input as a cell array of stings.");
        nKeyOpt[iTyp]= mxGetNumberOfElements(TypPtr2);
        if (nKeyOpt[iTyp] > maxKeyOpts) throw("Number of keyoptions is too large.");
        for (unsigned int iKeyOpt=0; iKeyOpt<nKeyOpt[iTyp]; iKeyOpt++)
        {
          const mxArray* keyOptPtr=mxGetCell(TypPtr2,iKeyOpt);
          if (!mxIsChar(keyOptPtr)) throw("Keyopts should be input as a cell array of stings.");
          TypeKeyOpts[iTyp+nEltType*iKeyOpt] = mxArrayToString(keyOptPtr);
        }
      }
      else nKeyOpt[iTyp]=0;
    }
    
    const unsigned int probDim = bemDimension(Elt,nElt,TypeID,TypeName,TypeKeyOpts,nKeyOpt,nEltType);
    
    plhs[0]=mxCreateDoubleScalar(double(probDim));

    // DEALLOCATE MEMORY ALLOCATED BY "mxArrayToString" IN TYPE DEFINITIONS
    for (unsigned int iTyp=0; iTyp<nEltType; iTyp++)
    {
      mxFree(TypeName[iTyp]); 
      for (unsigned int iKeyOpt=0; iKeyOpt<nKeyOpt[iTyp]; iKeyOpt++) mxFree(TypeKeyOpts[iTyp+nEltType*iKeyOpt]);
    }
    delete [] TypeID;
    delete [] nKeyOpt;
    delete [] TypeName;
    delete [] TypeKeyOpts;
    
  }
  catch (const char* exception)
  {
    const char* error = "Error using ==> bemdimension\n";
    char msg[512];
    strcpy(msg,error);
    strcat(msg,exception);
    mexErrMsgTxt(msg);
  }
}
