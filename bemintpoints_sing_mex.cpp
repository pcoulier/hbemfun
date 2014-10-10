/*BEMINTPOINTS  Boundary element integration points.
 *
 *   X = BEMINTPOINTSSING(nod,elt,typ) returns the coordinates of the boundary
 *   element integration points.
 *
 *   nod Nodes.
 *   elt Elements.
 *   typ Element types.
 *   X   Boundary element integration point coordinates (nPts * 3).
 *
 */

/* $Make: mex -O -output bemintpoints bemintpoints_mex.cpp eltdef.cpp gausspw.cpp 
                                      bemcollpoints.cpp shapefun.cpp bemdimension.cpp$*/
#include "mex.h"
#include "eltdef.h"
#include "gausspw.h"
#include "bemcollpoints.h"
#include "shapefun.h"
#include "bemdimension.h"
#include "checklicense.h"
#include <complex>

using namespace std;

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{
  try
  {
    checklicense();
    
    if (nrhs<3) throw("Not enough input arguments.");
    if (nrhs>3) throw("Too many input arguments.");
    if (nlhs>1) throw("Too many output arguments.");

    if (!mxIsNumeric(prhs[0])) throw("Input argument 'Nod' must be numeric.");
    if (mxIsSparse(prhs[0])) throw("Input argument 'Nod' must not be sparse.");
    if (mxIsComplex(prhs[0])) throw("Input argument 'Nod' must be real.");
    if (!(mxGetN(prhs[0])==4)) throw("Input argument 'Nod' should have 4 columns.");
    const int nNod=mxGetM(prhs[0]);
    const double* const Nod=mxGetPr(prhs[0]);

    if (!mxIsNumeric(prhs[1])) throw("Input argument 'Elt' must be numeric.");
    if (mxIsSparse(prhs[1])) throw("Input argument 'Elt' must not be sparse.");
    if (mxIsComplex(prhs[1])) throw("Input argument 'Elt' must be real.");
    if (mxGetN(prhs[1])<=3) throw("Input argument 'Elt' should have at least 3 columns.");
    const double* const Elt=mxGetPr(prhs[1]);
    const int nElt=mxGetM(prhs[1]);

    bool keyOpts=true;
    if (mxGetN(prhs[2])==3) keyOpts=true;
    else if  (mxGetN(prhs[2])==2) keyOpts=false;
    else throw("Input argument 'typ' should have 2 or 3 columns.");
    if (!(mxIsCell(prhs[2]))) throw("Input argument 'typ' should be a cell array.");
    const int nEltType=mxGetM(prhs[2]);
    const int maxKeyOpts = 50;  // Maximum number of keyoptions per element type
    int* const TypeID=new(nothrow) int[nEltType];
      if (TypeID==0) throw("Out of memory.");
    int* const nKeyOpt=new(nothrow) int[nEltType];
      if (nKeyOpt==0) throw("Out of memory.");
    char** const TypeName=new(nothrow) char*[nEltType];
      if (TypeName==0) throw("Out of memory.");
    char** const TypeKeyOpts=new(nothrow) char*[nEltType*maxKeyOpts];
      if (TypeKeyOpts==0) throw("Out of memory.");
    for (int iTyp=0; iTyp<nEltType; iTyp++)
    { 
      // TypeID
      const mxArray* TypPtr0=mxGetCell(prhs[2],iTyp+nEltType*0);
      if (!mxIsNumeric(TypPtr0)) throw("Type ID should be numeric.");
      if (mxIsSparse(TypPtr0)) throw("Type ID should not be sparse.");
      if (mxIsComplex(TypPtr0)) throw("Type ID should not be complex.");
      if (!(mxGetNumberOfElements(TypPtr0)==1)) throw("Type ID should be a scalar.");
      TypeID[iTyp]= int(mxGetScalar(TypPtr0));
      
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
        for (int iKeyOpt=0; iKeyOpt<nKeyOpt[iTyp]; iKeyOpt++)
        {
          const mxArray* keyOptPtr=mxGetCell(TypPtr2,iKeyOpt);
          if (!mxIsChar(keyOptPtr)) throw("Keyopts should be input as a cell array of stings.");
          TypeKeyOpts[iTyp+nEltType*iKeyOpt] = mxArrayToString(keyOptPtr);
        }
      }
      else nKeyOpt[iTyp]=0;
    }
    
    const int probDim=bemDimension(Elt,nElt,TypeID,TypeName,TypeKeyOpts,nKeyOpt,nEltType);
    
    // Determine total number of integration points
    int nIntPoint=0;
    for (int iElt=0; iElt<nElt; iElt++)
    {
      int EltType = int(Elt[nElt+iElt]);
      int EltParent;
      int nEltNod;
      int nEltColl;
      int EltShapeN;
      int EltShapeM;
      int EltDim;
      int AxiSym;
      int Periodic;
      int nGauss;
      int nEltDiv; 
      int nGaussSing;
      int nEltDivSing; 
      eltdef(EltType,TypeID,TypeName,TypeKeyOpts,nKeyOpt,nEltType,EltParent,
             nEltNod,nEltColl,EltShapeN,EltShapeM,EltDim,AxiSym,Periodic,nGauss,nEltDiv,
             nGaussSing,nEltDivSing);
             
      if (probDim==2) nIntPoint+=nGauss*nEltDiv;
      else
      {
        if (EltParent==1) nIntPoint+=nGauss*nEltDiv;
        else if (EltParent==2) nIntPoint+=nEltDiv*nEltDiv*nGauss*nGauss;
      }
    }
    
    plhs[0]=mxCreateDoubleMatrix(nIntPoint,3,mxREAL);
    double* const intPoints=mxGetPr(plhs[0]);
        
    // Determine integration point coordinates
    int iIntPoint=0;
    for (int iElt=0; iElt<nElt; iElt++)
    {
      int EltType = int(Elt[nElt+iElt]);
      int EltParent;
      int nEltNod;
      int nEltColl;
      int EltShapeN;
      int EltShapeM;
      int EltDim;
      int AxiSym;
      int Periodic;
      int nGauss;
      int nEltDiv; 
      int nGaussSing;
      int nEltDivSing;
      eltdef(EltType,TypeID,TypeName,TypeKeyOpts,nKeyOpt,nEltType,EltParent,
             nEltNod,nEltColl,EltShapeN,EltShapeM,EltDim,AxiSym,Periodic,nGauss,nEltDiv,
             nGaussSing,nEltDivSing);
             
      int nXi=0;
      if (probDim==2) nXi=nGauss*nEltDiv;
      else
      {
        if (EltParent==1) nXi=nGauss*nEltDiv;
        else if (EltParent==2) nXi=nEltDiv*nEltDiv*nGauss*nGauss;
      }
      
      int NodIndex;
      int NodID;
      double* const EltNod =new(nothrow) double[3*nEltNod];
      if (EltNod==0) throw("Out of memory.");
      
      // Determine coordinates of element nodes (element iElt).
      for (int iEltNod=0; iEltNod<nEltNod; iEltNod++)
      {
        NodID=int(Elt[(2+iEltNod)*nElt+iElt]);
        BemNodeIndex(Nod,nNod,NodID,NodIndex);
        EltNod[0*nEltNod+iEltNod]=Nod[1*nNod+NodIndex];
        EltNod[1*nEltNod+iEltNod]=Nod[2*nNod+NodIndex];
        EltNod[2*nEltNod+iEltNod]=Nod[3*nNod+NodIndex];
      }
      
      if (probDim==2)
      {
        double* const xi=new(nothrow) double[nXi];
          if (xi==0) throw("Out of memory.");
        double* const H=new(nothrow) double[nXi];
          if (H==0) throw("Out of memory.");
        double* const N=new(nothrow) double[nXi*nEltNod];
          if (N==0) throw("Out of memory.");  
        double* const xiCart=new(nothrow) double[2*nXi];
          if (xiCart==0) throw("Out of memory.");
        gausspw1D(nEltDiv,nGauss,xi,H);
        shapefun(EltShapeN,nXi,xi,N); 
        
        
        for (int icomp=0; icomp<2*nXi; icomp++) xiCart[icomp]=0.0;
        for (int iXi=0; iXi<nXi; iXi++)
        {
          for (int iEltNod=0; iEltNod<nEltNod; iEltNod++)
          {
            xiCart[2*iXi+0]+=N[nEltNod*iXi+iEltNod]*EltNod[0*nEltNod+iEltNod];
            xiCart[2*iXi+1]+=N[nEltNod*iXi+iEltNod]*EltNod[2*nEltNod+iEltNod];
          }
          intPoints[0*nIntPoint+iIntPoint]=xiCart[2*iXi+0];
          intPoints[2*nIntPoint+iIntPoint]=xiCart[2*iXi+1];
          iIntPoint+=1;
        }
        delete [] xi;
        delete [] H;
        delete [] N;
        delete [] xiCart;
      }
      else if (probDim==3)
      {
        double* const xi=new(nothrow) double[2*nXi];
          if (xi==0) throw("Out of memory.");
        double* const H=new(nothrow) double[nXi];
          if (H==0) throw("Out of memory.");
        double* const N=new(nothrow) double[nXi*nEltNod];
          if (N==0) throw("Out of memory.");  
        double* const xiCart=new(nothrow) double[3*nXi];
          if (xiCart==0) throw("Out of memory.");
          
        if (EltParent == 1) gausspwtri(nGauss,xi,H);
        else gausspw2D(nEltDiv,nGauss,xi,H);
        shapefun(EltShapeN,nXi,xi,N); 
        
        for (int icomp=0; icomp<3*nXi; icomp++) xiCart[icomp]=0.0;
        for (int iXi=0; iXi<nXi; iXi++)
        {
          for (int iEltNod=0; iEltNod<nEltNod; iEltNod++)
          {
            xiCart[3*iXi+0]+=N[nEltNod*iXi+iEltNod]*EltNod[0*nEltNod+iEltNod];
            xiCart[3*iXi+1]+=N[nEltNod*iXi+iEltNod]*EltNod[1*nEltNod+iEltNod];
            xiCart[3*iXi+2]+=N[nEltNod*iXi+iEltNod]*EltNod[2*nEltNod+iEltNod];
          }
          intPoints[0*nIntPoint+iIntPoint]=xiCart[3*iXi+0];
          intPoints[1*nIntPoint+iIntPoint]=xiCart[3*iXi+1];
          intPoints[2*nIntPoint+iIntPoint]=xiCart[3*iXi+2];
          iIntPoint+=1;
        }
        delete [] xi;
        delete [] H;
        delete [] N;
        delete [] xiCart;
      }
      delete [] EltNod;
    }

    // DEALLOCATE MEMORY ALLOCATED BY "mxArrayToString" IN TYPE DEFINITIONS
    for (int iTyp=0; iTyp<nEltType; iTyp++)
    {
      mxFree(TypeName[iTyp]);
      for (int iKeyOpt=0; iKeyOpt<nKeyOpt[iTyp]; iKeyOpt++) mxFree(TypeKeyOpts[iTyp+nEltType*iKeyOpt]);
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
