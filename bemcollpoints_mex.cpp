/*BEMCOLLPOINTS   Boundary element collocation points.
 *
 *   [col,colTyp,ID] = BEMCOLLPOINTS(nod,elt,typ) returns the boundary
 *   element collocation points.
 *
 *   nod    Nodes. 
 *   elt    Elements.
 *   typ    Element types.
 *   col    Collocation point coordinates (nCol * 3).
 *   colTyp Collocation type (nCol * 1). 1 for a centroid collocation, 2 for a
 *          nodal collocation point.
 *   ID     The corresponding element number (if colTyp=1) or node number 
 *          (if colTyp=2).
 */

/* $Make: mex -O -output bemcollpoints bemcollpoints_mex.cpp eltdef.cpp
                         bemcollpoints.cpp shapefun.cpp$*/
                         
#include "mex.h"
#include "string"
#include "shapefun.h"
#include "bemcollpoints.h"
#include <new>
using namespace std;

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{
   try
   {
    if (nrhs<3) throw("Not enough input arguments.");
    if (nrhs>3) throw("Too many input arguments.");
    if (nlhs>3) throw("Too many output arguments.");

    // PROCESS ARGUMENT "NOD"
    if (!mxIsNumeric(prhs[0])) throw("Input argument 'nod' must be numeric.");
    if (mxIsSparse(prhs[0])) throw("Input argument 'nod' must not be sparse.");
    if (mxIsComplex(prhs[0])) throw("Input argument 'nod' must be real.");
    if (!(mxGetN(prhs[0])==4)) throw("Input argument 'nod' should have 4 columns.");
    const double* const Nod=mxGetPr(prhs[0]);
    const unsigned int nNod=mxGetM(prhs[0]);

    if (!mxIsNumeric(prhs[1])) throw("Input argument 'elt' must be numeric.");
    if (mxIsSparse(prhs[1])) throw("Input argument 'elt' must not be sparse.");
    if (mxIsComplex(prhs[1])) throw("Input argument 'elt' must be real.");
    if (mxGetN(prhs[1])<=2) throw("Input argument 'elt' should have at least 2 columns.");
    const double* const Elt=mxGetPr(prhs[1]);
    const unsigned int nElt=mxGetM(prhs[1]);
    const unsigned int maxEltCol=mxGetN(prhs[1]);

    bool keyOpts=true;
    if (mxGetN(prhs[2])==3) keyOpts=true;
    else if  (mxGetN(prhs[2])==2) keyOpts=false;
    else throw("Input argument 'typ' should have 2 or 3 columns.");
    if (!(mxIsCell(prhs[2]))) throw("Input argument 'typ' should be a cell array.");
    const unsigned int nEltType=mxGetM(prhs[2]);
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
      const mxArray* TypPtr0=mxGetCell(prhs[2],iTyp+nEltType*0);
      if (!mxIsNumeric(TypPtr0)) throw("Type ID should be numeric.");
      if (mxIsSparse(TypPtr0)) throw("Type ID should not be sparse.");
      if (mxIsComplex(TypPtr0)) throw("Type ID should not be complex.");
      if (!(mxGetNumberOfElements(TypPtr0)==1)) throw("Type ID should be a scalar.");
      TypeID[iTyp]= (unsigned int)(mxGetScalar(TypPtr0));
      
      // TypeName
      const mxArray* TypPtr1=mxGetCell(prhs[2],iTyp+nEltType*1);
      if (!mxIsChar(TypPtr1)) throw("Element types should be input as stings.");
      TypeName[iTyp] =  mxArrayToString(TypPtr1);
      
      // TypeKeyOpts
      if (keyOpts)
      {
        const mxArray* TypPtr2=mxGetCell(prhs[2],iTyp+nEltType*2); // Keyoptions cell array
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
    
    unsigned int* const NodalColl=new(nothrow) unsigned int[nNod];
    if (NodalColl==0) throw("Out of memory.");
    unsigned int* const CentroidColl=new(nothrow) unsigned int[nElt];
    if (CentroidColl==0) throw("Out of memory.");
    unsigned int nCentroidColl;
    unsigned int nNodalColl;
    
    BemCollPoints(Elt,Nod,TypeID,nKeyOpt,TypeName,TypeKeyOpts,nEltType,nElt,maxEltCol,nNod,NodalColl,CentroidColl,nNodalColl,nCentroidColl);
    unsigned int nTotalColl = nNodalColl + nCentroidColl;
    
    double* const CollPoints = new(nothrow) double[nTotalColl*5];
    if (CollPoints==0) throw("Out of memory.");
    
    BemCollCoords(Elt,Nod,TypeID,nKeyOpt,TypeName,TypeKeyOpts,nEltType,CentroidColl,NodalColl,CollPoints,nTotalColl,nElt,nNod);
    
    // OUTPUT ARGUMENTS
    plhs[0]=mxCreateDoubleMatrix(nTotalColl,3,mxREAL);
    double* const CollCoords = mxGetPr(plhs[0]);
    for (unsigned int iColl=0; iColl<nTotalColl; iColl++)
    {
      CollCoords[nTotalColl*0+iColl]=CollPoints[nTotalColl*2+iColl];
      CollCoords[nTotalColl*1+iColl]=CollPoints[nTotalColl*3+iColl];
      CollCoords[nTotalColl*2+iColl]=CollPoints[nTotalColl*4+iColl];
    }
    if (nlhs>1)
    {
      plhs[1]=mxCreateDoubleMatrix(nTotalColl,1,mxREAL);
      double* const CollType = mxGetPr(plhs[1]);
      for (unsigned int iColl=0; iColl<nTotalColl; iColl++) CollType[iColl]=CollPoints[iColl];
    }
    if (nlhs>2)
    {
      plhs[2]=mxCreateDoubleMatrix(nTotalColl,1,mxREAL);
      double* const CollID = mxGetPr(plhs[2]);
      for (unsigned int iColl=0; iColl<nTotalColl; iColl++) CollID[iColl]=CollPoints[nTotalColl+iColl];
    }
    
    
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
    delete [] CollPoints;
    delete [] NodalColl;
    delete [] CentroidColl;
  }
  catch (const char* exception)
  {
    mexErrMsgTxt(exception);
  }
}
