/*BEMNORMAL   Element normals.
 *
 *   normal = BEMNORMAL(nod,elt,typ,xi) computes the normal in the specified
 *   points for all elements of the boundary element mesh.
 *
 *   nod    Nodes. 
 *   elt    Elements.
 *   typ    Element types.
 *   Xi     Natural coordinates of the sampling points (nXi * 2). 
 *          For 1D elements, the second column is not used.
 *   normal Element normals (3 * nXi * nElt).
 */

/* $Make: mex -O -output bemnormal bemnormal_mex.cpp eltdef.cpp shapefun.cpp 
                         bemnormal.cpp bemcollpoints.cpp $ */                                                                  

#include "mex.h"
#include <string.h>
#include "eltdef.h"
#include "shapefun.h"
#include "bemcollpoints.h"
#include "bemnormal.h"
#include "bemdimension.h"
#include "checklicense.h"
#include <new>

using namespace std;

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{
  try
  {
    checklicense();
    
    if (!nrhs==4) throw("4 input arguments required.");
    if (nlhs>1) throw("Too many output arguments.");

    if (!mxIsNumeric(prhs[0])) throw("Input argument 'Nod' must be numeric.");
    if (mxIsSparse(prhs[0])) throw("Input argument 'Nod' must not be sparse.");
    if (mxIsComplex(prhs[0])) throw("Input argument 'Nod' must be real.");
    if (!(mxGetN(prhs[0])==4)) throw("Input argument 'Nod' should have 4 columns.");
    const unsigned int nNod=mxGetM(prhs[0]);
    const double* const Nod=mxGetPr(prhs[0]);

    if (!mxIsNumeric(prhs[1])) throw("Input argument 'Elt' must be numeric.");
    if (mxIsSparse(prhs[1])) throw("Input argument 'Elt' must not be sparse.");
    if (mxIsComplex(prhs[1])) throw("Input argument 'Elt' must be real.");
    if (mxGetN(prhs[1])<=3) throw("Input argument 'Elt' should have at least 3 columns.");
    const double* const Elt=mxGetPr(prhs[1]);
    const unsigned int nElt=mxGetM(prhs[1]);
 
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
    //const unsigned int probDim=bemDimension(Elt,nElt,TypeID,TypeName,TypeKeyOpts,nKeyOpt,nEltType);/
    
    if (!mxIsNumeric(prhs[3])) throw("Input argument 'xi' must be numeric.");
    if (mxIsSparse(prhs[3])) throw("Input argument 'xi' must not be sparse.");
    if (mxIsComplex(prhs[3])) throw("Input argument 'xi' must be real.");
    if (!(mxGetN(prhs[3])==2)) throw("Input argument 'xi' should have 2 columns.");
    const unsigned int nXi=mxGetM(prhs[3]);
    const double* const xi=mxGetPr(prhs[3]);

    // OUTPUT ARGUMENT: NORMAL
    size_t* const outDims = new(nothrow) size_t[3];
    if (outDims==0) throw("Out of memory.");
    outDims[0]=3;
    outDims[1]=nXi;
    outDims[2]=nElt;
    plhs[0]=mxCreateNumericArray(3,outDims,mxDOUBLE_CLASS,mxREAL);
    double* const normal=mxGetPr(plhs[0]);

    for (unsigned int iElt=0; iElt<nElt; iElt++)
    {
      unsigned int eltType = (unsigned int)(Elt[nElt+iElt]);
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
      eltdef(eltType,TypeID,TypeName,TypeKeyOpts,nKeyOpt,nEltType,EltParent,nEltNod,
             nEltCol,EltShapeN,EltShapeM,EltDim,AxiSym,Periodic,nGauss,nEltDiv,nGaussSing,nEltDivSing);

      int NodIndex;
      unsigned int NodID;
      double* const EltNod =new(nothrow) double[3*nEltNod];
      if (EltNod==0) throw("Out of memory.");

      for (unsigned int iEltNod=0; iEltNod<nEltNod; iEltNod++)
      {
        NodID=(unsigned int)(Elt[(2+iEltNod)*nElt+iElt]);
        BemNodeIndex(Nod,nNod,NodID,NodIndex);
        EltNod[0*nEltNod+iEltNod]=Nod[1*nNod+NodIndex];
        EltNod[1*nEltNod+iEltNod]=Nod[2*nNod+NodIndex];
        EltNod[2*nEltNod+iEltNod]=Nod[3*nNod+NodIndex];
      }

      double* const dN=new(nothrow) double[2*nXi*nEltNod];
      if (dN==0) throw("Out of memory.");
      double* const nat=new(nothrow) double[6*nXi];
      if (nat==0) throw("Out of memory.");
      double* const eltNormal=new(nothrow) double[3*nXi];
      if (eltNormal==0) throw("Out of memory.");

      shapederiv(EltShapeN,nXi,xi,dN);
      shapenatcoord(dN,nEltNod,nXi,EltNod,nat,EltDim);

      bemnormal(nat,nXi,EltDim,eltNormal);

      for (unsigned int iXi=0; iXi<nXi; iXi++)
      {
        normal[3*nXi*iElt+3*iXi+0]=eltNormal[3*iXi+0];
        normal[3*nXi*iElt+3*iXi+1]=eltNormal[3*iXi+1];
        normal[3*nXi*iElt+3*iXi+2]=eltNormal[3*iXi+2];
      }

      delete [] eltNormal;
      delete [] EltNod;
      delete [] dN;
      delete [] nat;
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
    delete [] outDims;
  }
  catch (const char* exception)
  {
    mexErrMsgTxt(exception);
  }
}
