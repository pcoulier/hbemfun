/*BEMELTDEF   Boundary element properties.
 *
 *   [map,nNod,nCol,nID,mID,nodDef,eltDim,axiSym,nGauss,nEltDiv]
 *                                                         =BEMELTDEF(typID,typ)
 *   returns various element properties.
 *
 *   typID   Element type ID.
 *   typ     Element types.
 *   map     Parent element mapping (for 2D elements only). 1 for a triangular 
 *           mapping, 2 for a rectangular mapping. 0 for 1D elements
 *   nNod    Number of element nodes.
 *   nCol    Number of element collocation points.
 *   nID     Shape function ID for element Geometry (see BEMSHAPE).
 *   mID     Boundary element interpolation function ID (see BEMSHAPE).
 *   nodDef  Natural coordinates of the nodes (Nnod * 2).
 *   eltDim  Element dimension. 1 for 1D element, 2 for 2D element.
 *   axiSym  1 for an axisymmetric geometry, 0 otherwise.
 *   nGauss  Number of Gaussian points.
 *   nEltDiv Number of element divisions for integration.
 */

/* $Make: mex -O -output bemeltdef bemeltdef_mex.cpp eltdef.cpp$ */

#include "mex.h"
#include <string.h>
#include "eltdef.h"
#include <new>
using namespace std;


void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{
  try
  {
    if (nrhs<2) throw("Not enough input arguments.");
    if (nrhs>2) throw("Too many input arguments.");
    if (nlhs>10) throw("Too many output arguments.");
    
    if (!mxIsNumeric(prhs[0])) throw("Input argument 'eltTyp' must be numeric.");
    if (mxIsSparse(prhs[0])) throw("Input argument 'eltTyp' must not be sparse.");
    if (mxIsComplex(prhs[0])) throw("Input argument 'eltTyp' must be real.");
    if (!(mxGetNumberOfElements(prhs[0])==1)) throw("Input argument 'eltTyp' should be a scalar.");
    const double* const Inarg=mxGetPr(prhs[0]);
    const unsigned int eltType = (unsigned int)(Inarg[0]);
    
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
    
    unsigned int EltParent;
    unsigned int nEltNod;
    unsigned int nEltCol;
    unsigned int EltShapeN;
    unsigned int EltShapeM;
    unsigned int EltDim;
    unsigned int AxiSym;
    unsigned int Periodic;
    unsigned int nGauss;
    unsigned int nEltDiv;
    unsigned int nGaussSing;
    unsigned int nEltDivSing;
    eltdef(eltType,TypeID,TypeName,TypeKeyOpts,nKeyOpt,nEltType,EltParent,
           nEltNod,nEltCol,EltShapeN,EltShapeM,EltDim,AxiSym,Periodic,nGauss,
           nEltDiv,nGaussSing,nEltDivSing);
           
    plhs[0]=mxCreateDoubleScalar(double(EltParent));
    if (nlhs>1) plhs[1]=mxCreateDoubleScalar(double(nEltNod));
    if (nlhs>2) plhs[2]=mxCreateDoubleScalar(double(nEltCol));
    if (nlhs>3) plhs[3]=mxCreateDoubleScalar(double(EltShapeN));
    if (nlhs>4) plhs[4]=mxCreateDoubleScalar(double(EltShapeM));
    if (nlhs>5)
    {
      plhs[5]=mxCreateDoubleMatrix(nEltNod,2,mxREAL);
      double* const eltNodXi=mxGetPr(plhs[5]);
      eltnoddef(eltType,TypeID,TypeName,nEltType,eltNodXi);
    }
    if (nlhs>6) plhs[6]=mxCreateDoubleScalar(double(EltDim));
    if (nlhs>7) plhs[7]=mxCreateDoubleScalar(double(AxiSym));
    if (nlhs>8) plhs[8]=mxCreateDoubleScalar(double(nGauss));
    if (nlhs>9) plhs[9]=mxCreateDoubleScalar(double(nEltDiv));
    
    
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
    mexErrMsgTxt(exception);
  }
}
