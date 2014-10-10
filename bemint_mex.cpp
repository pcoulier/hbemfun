/*BEMINT  Integration of displacements and tractions over the boundary element mesh.
 *
 *   I = BEMINT(nod,elt,typ,t,u) computes the integral of the product of
 *   stresses and displacements over the boundary, defined as:
 *
 *                    /
 *            I_ij =  |   t_i *u_j     dS
 *                    /
 *                  Sigma
 *
 *   where both the displacements and tractions are approximated by means
 *   of the boundary element shape functions. The vectors t and u may be used 
 *   to evaluate any regular integral of vectors over the boundary Sigma.
 *
 *   nod   Nodes (nNod * 4). Each row has the layout [nodID x y z] where
 *         nodID is the node number and x, y, and z are the nodal
 *         coordinates.
 *   elt   Elements (nElt * nColumn). Each row has the layout
 *         [eltID typID n1 n2 n3 ... ] where eltID is the element number,
 *         typID is the element type number and n1, n2, n3, ... are the node 
 *         numbers representing the nodal connectivity of the element.
 *   typ   Element type definitions. Cell array with the layout
 *         {{typID type keyOpts} ... } where typID is the element type number,
 *         type is the element type (string) and keyOpts is a cell array of
 *         strings with key options.
 *   t     Boundary element tractions (nDof * ntMode * ...) for ntMode modes.
 *         The number of degrees of freedom is equal to:
 *            nDof=1*nBemCol      2D out-of-plane
 *            nDof=2*nBemCol      2D in-plane or axisymmetric
 *            nDof=3*nBemCol      2.5D or 3D
 *         where nBemColl is the total number of collocation points.
 *   u     Boundary element displacements (nDof * nuMode * ...) for nuMode modes.
 *   I     The resulting boundary integral (ntMode * nuMode * ...).
 *
 */

/* $Make: mex -O -output bemint bemint_mex.cpp eltdef.cpp bemcollpoints.cpp
                      shapefun.cpp gausspw.cpp bemint.cpp bemisaxisym.cpp$*/

#include "mex.h"
#include "string"
#include "eltdef.h"
#include "bemcollpoints.h"
#include "bemint.h"
#include "bemisaxisym.h"
#include "checklicense.h"
#include <complex>
#include <new>
using namespace std;

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{
  try
  {
    checklicense();
    
    if (nrhs<5) throw("Not enough input arguments.");
    if (nrhs>5) throw("Too many input arguments.");
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
    
    const bool probAxi = isAxisym(Elt,nElt,TypeID,TypeName,TypeKeyOpts,nKeyOpt,nEltType);

    // DETERMINE COLLOCATION POINTS: NODAL AND CENTROID
    unsigned int* const NodalColl=new(nothrow) unsigned int[nNod];
    if (NodalColl==0) throw("Out of memory.");
    unsigned int* const CentroidColl=new(nothrow) unsigned int[nElt];
    if (CentroidColl==0) throw("Out of memory.");
    unsigned int nCentroidColl=0;
    unsigned int nNodalColl=0;
    BemCollPoints(Elt,Nod,TypeID,nKeyOpt,TypeName,TypeKeyOpts,nEltType, 
                   nElt,maxEltCol,nNod,NodalColl,CentroidColl,nNodalColl,nCentroidColl);


    unsigned int nTotalColl = nNodalColl + nCentroidColl;
    double* const CollPoints=new(nothrow) double[5*nTotalColl];
    if (CollPoints==0) throw("Out of memory.");
    BemCollCoords(Elt,Nod,TypeID,nKeyOpt,TypeName,TypeKeyOpts,nEltType,
                  CentroidColl,NodalColl,CollPoints,nTotalColl,nElt,nNod);
    

    if (!mxIsNumeric(prhs[3])) throw("Input argument 't' must be numeric.");
    if (mxIsSparse(prhs[3])) throw("Input argument 't' must not be sparse.");
    if (mxIsEmpty(prhs[3])) throw("Input argument 't' must not be empty.");
    const unsigned int ntdim=mxGetNumberOfDimensions(prhs[3]);
    const size_t* const tdim=mxGetDimensions(prhs[3]);
    const unsigned int nDof=tdim[0];
    unsigned int ntMode =((ntdim>1)? tdim[1]:1);
    unsigned int ntSet  =((ntdim>2)? tdim[2]:1);

    unsigned int nColDof=0;
    if (nDof==1*nTotalColl) nColDof=1;
    else if (nDof==2*nTotalColl) nColDof=2;
    else if (nDof==3*nTotalColl) nColDof=3;
    else throw("The first dimension of 't' does not correspond with the number of DOFs in the boundary element mesh.");

    const double* const tReal = mxGetPr(prhs[3]);
    const double* const tImag = mxGetPi(prhs[3]);
    complex<double>* const t = new(nothrow) complex<double>[mxGetNumberOfElements(prhs[3])];
    if (t==0) throw("Out of memory.");
    for (unsigned int i=0; (unsigned)i<mxGetNumberOfElements(prhs[3]); i++)
      t[i]=complex<double>(tReal[i],(mxIsComplex(prhs[3]) ? tImag[i] : 0.0));

    if (!mxIsNumeric(prhs[4])) throw("Input argument 'u' must be numeric.");
    if (mxIsSparse(prhs[4])) throw("Input argument 'u' must not be sparse.");
    if (mxIsEmpty(prhs[4])) throw("Input argument 'u' must not be empty.");
    const unsigned int nudim=mxGetNumberOfDimensions(prhs[4]);
    const size_t* const udim=mxGetDimensions(prhs[4]);
    if (!(udim[0]==(unsigned)nDof)) throw("The first dimension of 'u' does not correspond with the number of DOFs.");
    unsigned int nuMode =((nudim>1)? udim[1]:1);

    const double* const uReal = mxGetPr(prhs[4]);
    const double* const uImag = mxGetPi(prhs[4]);
    complex<double>* const u = new(nothrow) complex<double>[mxGetNumberOfElements(prhs[4])];
    if (u==0) throw("Out of memory.");
    for (unsigned int i=0; (unsigned)i<mxGetNumberOfElements(prhs[4]); i++)
      u[i]=complex<double>(uReal[i],(mxIsComplex(prhs[4]) ? uImag[i] : 0.0));

    const bool cmplx = (mxIsComplex(prhs[3]) || mxIsComplex(prhs[4]));

    complex<double>* const OutMat = new(nothrow) complex<double>[ntMode*nuMode*ntSet];
    if (OutMat==0) throw("Out of memory.");
    for (unsigned int i=0; i<ntMode*nuMode*ntSet; i++) OutMat[i]=0.0;

    for (unsigned int iElt=0; iElt<nElt; iElt++)
    {
      unsigned int EltType = (unsigned int)(Elt[nElt+iElt]);
      unsigned int EltParent;
      unsigned int nEltNod;
      unsigned int nEltColl;
      unsigned int EltShapeN;
      unsigned int EltShapeM;
      unsigned int EltDim;
      unsigned int AxiSym;
      unsigned int Periodic;
      unsigned int nGauss;
      unsigned int nEltDiv;
      unsigned int nGaussSing;
      unsigned int nEltDivSing;
      eltdef(EltType,TypeID,TypeName,TypeKeyOpts,nKeyOpt,nEltType,EltParent,nEltNod,
             nEltColl,EltShapeN,EltShapeM,EltDim,AxiSym,Periodic,nGauss,nEltDiv,
             nGaussSing,nEltDivSing);

      // element collocation points
      unsigned int* const eltCollIndex=new(nothrow) unsigned int[nEltColl];
      if (eltCollIndex==0) throw("Out of memory.");
      BemEltCollIndex(Elt,iElt,nElt,CollPoints,nCentroidColl,nTotalColl,
                                                 nEltColl,nEltNod,eltCollIndex);
      
      // boundary element integration
        bemint(Nod,nNod,Elt,iElt,nElt,TypeID,nKeyOpt,TypeName,TypeKeyOpts,nEltType,
               eltCollIndex,OutMat,nDof,t,u,ntMode,nuMode,ntSet,nColDof,probAxi);
            
        delete [] eltCollIndex;
    }

    // OUTPUT ARGUMENT POINTER
    size_t* const OutDim = new(nothrow) size_t[3];
    if (OutDim==0) throw("Out of memory.");
    OutDim[0]=nuMode;
    OutDim[1]=ntMode;
    OutDim[2]=ntSet;
    plhs[0] =mxCreateNumericArray(3,OutDim,mxDOUBLE_CLASS,(cmplx ? mxCOMPLEX : mxREAL));
    double* const OutPr=mxGetPr(plhs[0]);
    double* const OutPi=mxGetPi(plhs[0]);
    for (unsigned int i=0; i<nuMode*ntMode*ntSet; i++)
    {
      OutPr[i] = real(OutMat[i]);
      if (cmplx) OutPi[i] = imag(OutMat[i]);
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
    delete [] NodalColl;
    delete [] CentroidColl;
    delete [] CollPoints;
    delete [] t;
    delete [] u;
    delete [] OutDim;
    delete [] OutMat;
  }
  catch (const char* exception)
  {
    mexErrMsgTxt(exception);
  }
}
