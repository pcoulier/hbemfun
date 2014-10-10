/*BEMXFER   Boundary element transfer matrices.
 *
 *   [Up,Tp] = BEMXFER(nod,elt,typ,rec,green,...) returns boundary element
 *   transfer matrices. The Green's functions are specified as a full space
 *   solution (green='fs***') or a user specified solution (green='user') and
 *   are integrated over the boundary element mesh defined by nodes, elements
 *   and element types. The wave field scattered by the boundary element mesh
 *   in a number of receiver points is computed as:
 *
 *     u_rec = Up*t - Tp*u
 *
 *   where t (nDof * 1) and u (nDof * 1) are the tractions and displacements on
 *   the boundary element mesh.
 *
 *   The routine detemines whether the receivers are located on the boundary 
 *   element mesh or not. For receivers located on the interface, the matrices 
 *   Up and Tp are derived from the boundary element shape functions. For 
 *   receivers not located on the boundary element mesh, the boundary integral 
 *   theorem is used.
 *
 *   Depending on the Green's function, the following syntax is used:
 *
 *   [Up,Tp] = BEMXFER(nod,elt,typ,rec,green,...)
 *   [Up,Tp] = BEMXFER(nod,elt,typ,rec,'fsgreen2d_outofplane0',mu)
 *   [Up,Tp] = BEMXFER(nod,elt,typ,rec,'fsgreen2d_inplane0',E,nu)
 *   [Up,Tp] = BEMXFER(nod,elt,typ,rec,'fsgreen3d0',E,nu)
 *   [Up,Tp] = BEMXFER(nod,elt,typ,rec,'fsgreen2d_outofplane',Cs,Ds,rho,omega)
 *   [Up,Tp] = BEMXFER(nod,elt,typ,rec,'fsgreen2d_inplane',Cs,Cp,Ds,Dp,rho,omega)
 *   [Up,Tp] = BEMXFER(nod,elt,typ,rec,'fsgreen3d',Cs,Cp,Ds,Dp,rho,omega)
 *   [Up,Tp] = BEMXFER(nod,elt,typ,rec,'fsgreenf',Cs,Cp,Ds,Dp,rho,py,omega)
 *   [Up,Tp] = BEMXFER(nod,elt,typ,rec,'user',zs,r,z,ug,sg)
 *
 *
 *   nod      Nodes (nNod * 4). Each row has the layout [nodID x y z] where
 *            nodID is the node number and x, y, and z are the nodal
 *            coordinates.
 *   elt      Elements (nElt * nColumn). Each row has the layout
 *            [eltID typID n1 n2 n3 ... ] where eltID is the element number,
 *            typID is the element type number and n1, n2, n3, ... are the node
 *            numbers representing the nodal connectivity of the element.
 *   typ      Element type definitions. Cell array with the layout
 *            {{typID type keyOpts} ... } where typID is the element type number,
 *            type is the element type (string) and keyOpts is a cell array of
 *            strings with key options.
 *   rec      Boundary element receiver points (|nRec| * 3). Each row represents
 *            a receiver and contains the receiver coordinates x, y, and z.
 *   green    Green's function (string). 'fs***' for a full-space solution or
 *            'user' for a user specified Green's function.
 *   mu       Shear modulus (1 * 1).
 *   E        Young's modulus (1 * 1).
 *   nu       Poisson coefficient (1 * 1).
 *   Cs       Shear wave velocity (1 * 1).
 *   Cp       Dilatational wave velocity (1 * 1).
 *   Ds       Shear damping ratio (1 * 1).
 *   Dp       Dilatational damping ratio (1 * 1).
 *   rho      Density (1 * 1).
 *   omega    Circular frequency (nFreq * 1).
 *   py       Slowness in y-direction, logarithmically sampled (nyWave * 1).
 *            If omega ~= 0, the wavenumber sampling is given by ky = omega * py.
 *            If omega == 0, the wavenumber sampling is given by ky = py.
 *   zs       Source locations (vertical coordinate) (nzSrc * 1).
 *   r        Receiver locations (x-coordinate) (nxRec * 1).
 *   z        Receiver locations (z-coordinate) (nzRec * 1).
 *   ug       Green's displacements.
 *   sg       Green's stresses.
 *   Up       Boundary element displacement system matrix (nRecDof * nDof * nSet).
 *   Tp       Boundary element traction system matrix (nRecDof * nDof * nSet).
 */

/* $Make: mex -O -output bemxfer bemxfer_mex.cpp eltdef.cpp bemcollpoints.cpp 
                                 shapefun.cpp bemnormal.cpp gausspw.cpp search1.cpp 
                                 bemxfer3d.cpp bemxfer2d.cpp bemxferaxi.cpp 
                                 bemdimension.cpp bemisaxisym.cpp greeneval2d.cpp 
                                 fsgreenf.cpp fsgreen3d.cpp fsgreen3dt.cpp 
                                 fsgreen2d_inplane.cpp fsgreen2d_outofplane.cpp 
                                 besselh.cpp greeneval3d.cpp greenrotate2d.cpp 
                                 boundaryrec2d.cpp boundaryrec3d.cpp fminstep.cpp 
                                 greenrotate3d.cpp checklicense.cpp ripemd128.cpp$*/


/*
mexFunction----> IntegrateGreenUser --->  bemIntegrate
            |
            |--> IntegrateFsGreenf  --->  bemIntegrate
            |
            |--> IntegrateFsGreen3d  --->  bemIntegrate
            |
            |--> IntegrateFsGreen2d_inplane --->  bemIntegrate
            |
            |--> IntegrateFsGreenf2d_outofplane --->  bemIntegrate
            |
            |--> ...
*/

#include "mex.h"
#include <string.h>
#include "eltdef.h"
#include "bemcollpoints.h"
#include "bemdimension.h"
#include "bemisaxisym.h"
#include "bemisperiodic.h"
#include "bemxfer2d.h"
#include "bemxfer3d.h"
#include "bemxfer3dperiodic.h"
#include "bemxferaxi.h"
#include "boundaryrec2d.h"
#include "boundaryrec3d.h"     
#include "checklicense.h"
#include <math.h>
#include <new>


#ifndef __GNUC__
#define strcasecmp _strcmpi
#endif

using namespace std;

//==============================================================================
void bemIntegrate(mxArray* plhs[], const bool& probAxi, const bool& probPeriodic, const unsigned int& probDim,
                  const unsigned int& nColDof, const bool& TmatOut,
                  const double* const Nod, const unsigned int& nNod,
                  const double* const Elt, const unsigned int& nElt,
                  const double* const Rec, const unsigned int& nRec,
                  const unsigned int* const TypeID,
                  const char* const TypeName[], const char* const TypeKeyOpts[],
                  const unsigned int* const nKeyOpt,
                  const unsigned int& nEltType, const double* const CollPoints,
                  const unsigned int& nTotalColl, const unsigned int& nCentroidColl,
                  const void* const* const greenPtr, const unsigned int& nGrSet,
                  const unsigned int& nugComp, const bool& ugCmplx,
                  const bool& tgCmplx, const unsigned int* const greenDim,
                  const unsigned int nGreenDim,
                  const double L, const double* const ky, const unsigned int nWave, 
                  const unsigned int nmax)
/* BemIntegrate performs the actual boundary element integration.
 * This function is called from each separate Green's function integration
 * separately.
 */
//==============================================================================
{
  // DEGREES OF FREEDOM
  unsigned int nDof=nColDof*nTotalColl;
  unsigned int nRecDof=nColDof*nRec;

  // OUTPUT ARGUMENT POINTERS
  const unsigned int nMatDim=(probPeriodic ? 3+nGreenDim : 2+nGreenDim);
  
  size_t* const MatDim = new(nothrow) size_t[nMatDim];
  if (MatDim==0) throw("Out of memory.");
  MatDim[0]=nRecDof;
  MatDim[1]=nDof;
  for (unsigned int iDim=0; iDim<nGreenDim; iDim++)  MatDim[2+iDim]=greenDim[iDim];
  if (probPeriodic) MatDim[nMatDim-1]=nWave;
  

  plhs[0]=mxCreateNumericArray(nMatDim,MatDim,mxDOUBLE_CLASS,mxCOMPLEX);
  double* const URe=mxGetPr(plhs[0]);
  double* const UIm=mxGetPi(plhs[0]);
  
  double* TRe = 0;
  double* TIm = 0;
  if (TmatOut)
  {
    plhs[1]=mxCreateNumericArray(nMatDim,MatDim,mxDOUBLE_CLASS,mxCOMPLEX);
    TRe=mxGetPr(plhs[1]);
    TIm=mxGetPi(plhs[1]);
  }

  bool* const boundaryRec=new(nothrow) bool[nRec];
        if (boundaryRec==0) throw("Out of memory.");
  for (unsigned int iRec=0; iRec<nRec; iRec++) boundaryRec[iRec]=false;
  
  for (unsigned int iElt=0; iElt<nElt; iElt++)
  {
    if (probDim==3)
    {
     boundaryRec3d(Nod,nNod,Elt,nElt,iElt,TypeID,TypeName,TypeKeyOpts,
                   nKeyOpt,nEltType,CollPoints,nTotalColl,nCentroidColl,
                   Rec,nRec,nRecDof,boundaryRec,TRe,TmatOut,nDof,nGrSet);
    }
    else
    {
      boundaryRec2d(Nod,nNod,Elt,nElt,iElt,TypeID,TypeName,TypeKeyOpts,
                    nKeyOpt,nEltType,CollPoints,nTotalColl,nCentroidColl,
                    Rec,nRec,nRecDof,boundaryRec,TRe,TmatOut,nDof,nGrSet,
                    nugComp,nColDof);
    }
  }

// Apply Boundary integral theorem for receivers not on the interface.
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

    unsigned int* const eltCollIndex=new(nothrow) unsigned int[nEltColl];
    if (eltCollIndex==0) throw("Out of memory.");
    BemEltCollIndex(Elt,iElt,nElt,CollPoints,nCentroidColl,nTotalColl,
                    nEltColl,nEltNod,eltCollIndex);

    if (probDim==3)
    {
      if (probPeriodic)
        bemxfer3dperiodic(Nod,nNod,Elt,iElt,nElt,eltCollIndex,Rec,nRec,boundaryRec,
                          URe,UIm,TRe,TIm,1,TmatOut,nDof,nRecDof,TypeID,nKeyOpt,TypeName,
                          TypeKeyOpts,nEltType,greenPtr,nGrSet,ugCmplx,tgCmplx,L,ky,nWave,nmax);
      else
        bemxfer3d(Nod,nNod,Elt,iElt,nElt,eltCollIndex,Rec,nRec,boundaryRec,
                  URe,UIm,TRe,TIm,1,TmatOut,nDof,nRecDof,TypeID,nKeyOpt,TypeName,
                  TypeKeyOpts,nEltType,greenPtr,nGrSet,ugCmplx,tgCmplx);
    }
    else if ((probDim==2)&& probAxi)
    {
       bemxferaxi(Nod,nNod,Elt,iElt,nElt,eltCollIndex,Rec,nRec,boundaryRec,
                  URe,UIm,TRe,TIm,1,TmatOut,nDof,nRecDof,TypeID,nKeyOpt,TypeName,
                  TypeKeyOpts,nEltType,greenPtr,nGrSet,ugCmplx,tgCmplx);
    }
    else
    {
       bemxfer2d(Nod,nNod,Elt,iElt,nElt,eltCollIndex,Rec,nRec,boundaryRec,
                 URe,UIm,TRe,TIm,TmatOut,nDof,nRecDof,TypeID,nKeyOpt,TypeName,
                 TypeKeyOpts,nEltType,greenPtr,nGrSet,nugComp,ugCmplx,tgCmplx);
    }
    delete [] eltCollIndex;
  }
  delete [] boundaryRec;
  delete [] MatDim;
}


//==============================================================================
void IntegrateGreenUser(mxArray* plhs[], int nrhs,
                        const mxArray* prhs[], const bool probAxi, const bool& probPeriodic,
                        const unsigned int& probDim, const double* const Nod,
                        const unsigned int& nNod, const double* const Elt,
                        const unsigned int& nElt, const double* const Rec, const unsigned int& nRec,
                        const unsigned int* const TypeID,
                        const char* const TypeName[], const char* const TypeKeyOpts[],
                        const unsigned int* const nKeyOpt, const unsigned int& nEltType,
                        const double* const CollPoints, const unsigned int& nTotalColl,
                        const unsigned int& nCentroidColl,
                        const bool& TmatOut)
//==============================================================================
{
  // INPUT ARGUMENT PROCESSING
  if (probPeriodic)
  {
    if (!(nrhs==13)) throw("Wrong number of input arguments.");
  }
  else
  {
    if (nrhs>10) throw("Too many input arguments.");
    if (nrhs<9) throw("Not enough input arguments.");
  }

  if (!mxIsNumeric(prhs[5])) throw("Input argument 'zs' must be numeric.");
  if (mxIsSparse(prhs[5])) throw("Input argument 'zs' must not be sparse.");
  if (mxIsEmpty(prhs[5])) throw("Input argument 'zs' must not be empty.");
  const unsigned int nzs=mxGetNumberOfElements(prhs[5]);
  const double* const zs=mxGetPr(prhs[5]);
  for (unsigned int izs=1; izs<nzs; izs++) if (!(zs[izs-1]<zs[izs])) throw("Input argument 'zs' must be monotonically increasing.");

  if (!mxIsNumeric(prhs[6])) throw("Input argument 'r' must be numeric.");
  if (mxIsSparse(prhs[6])) throw("Input argument 'r' must not be sparse.");
  if (mxIsEmpty(prhs[6])) throw("Input argument 'r' must not be empty.");
  const unsigned int nr=mxGetNumberOfElements(prhs[6]);
  const double* const r=mxGetPr(prhs[6]);
  for (unsigned int ir=1; ir<nr; ir++) if (!(r[ir-1]<r[ir])) throw("Input argument 'r' must be monotonically increasing.");

  if (!mxIsNumeric(prhs[7])) throw("Input argument 'z' must be numeric.");
  if (mxIsSparse(prhs[7])) throw("Input argument 'z' must not be sparse.");
  if (mxIsEmpty(prhs[7])) throw("Input argument 'z' must not be empty.");
  const unsigned int nz=mxGetNumberOfElements(prhs[7]);
  const double* const z=mxGetPr(prhs[7]);
  for (unsigned int iz=1; iz<nz; iz++) if (!(z[iz-1]<z[iz])) throw("Input argument 'z' must be monotonically increasing.");

  if (!mxIsNumeric(prhs[8])) throw("Input argument 'ug' must be numeric.");
  if (mxIsSparse(prhs[8])) throw("Input argument 'ug' must not be sparse.");
  if (mxIsEmpty(prhs[8])) throw("Input argument 'ug' must not be empty.");
  const unsigned int nugdim=mxGetNumberOfDimensions(prhs[8]);
  const size_t* const ugdim=mxGetDimensions(prhs[8]);
  const unsigned int nugComp=ugdim[0];

  const unsigned int nGreenDim=((nugdim>4)?(nugdim-4):1);
  unsigned int* const greenDim=new(nothrow) unsigned int[nGreenDim];
  unsigned int nGrSet=1;
  greenDim[0]=1;
  if (nugdim>4)
  {
    for (unsigned int iDim=4; iDim<nugdim; iDim++)
    {
      greenDim[iDim-4]=ugdim[iDim];
      nGrSet=nGrSet*ugdim[iDim];
    }
  }

  if (probDim==2)
  {
    if (!(probAxi) && (!((nugComp==1)|(nugComp==4)|(nugComp==9)))) throw("The first dimension of input argument 'ug' for a 2D problem must be 1 (out-of-plane), 4 (in-plane) or 9 for a 2.5D problem.");
    if (probAxi  && (!(nugComp==5))) throw("The first dimension of input argument 'ug' must be 4 for an axisymmetric problem .");
  }
  else if (probDim==3) if (!(nugComp==5)) throw("The first dimension of input argument 'ug' must be 5 for a 3D problem.");
  if (!((unsigned)nzs == ((nugdim>1)? ugdim[1]:1))) throw("Input arguments 'ug' and 'zs' are incompatible");
  if (!((unsigned)nr  == ((nugdim>2)? ugdim[2]:1))) throw("Input arguments 'ug' and 'r' are incompatible");
  if (!((unsigned)nz  == ((nugdim>3)? ugdim[3]:1))) throw("Input arguments 'ug' and 'z' are incompatible");
  const double* const ugRe=mxGetPr(prhs[8]);
  const double* const ugIm=mxGetPi(prhs[8]);
  const bool ugCmplx=mxIsComplex(prhs[8]);

  unsigned int ntgComp;
  if (nugComp==1) ntgComp=2;  // 2D, out-of-plane
  if (nugComp==4) ntgComp=6;  // 2D, in-plane
  if (nugComp==5) ntgComp=10; // 3D / axisymmetric
  if (nugComp==9) ntgComp=18; // 2.5D

  const bool sgIn=(nrhs==10);
  if (TmatOut)
  {
    if (sgIn)
    {
      if (!mxIsNumeric(prhs[9])) throw("Input argument 'sg' must be numeric.");
      if (mxIsSparse(prhs[9])) throw("Input argument 'sg' must not be sparse.");
      const unsigned int ntgdim=mxGetNumberOfDimensions(prhs[9]);
      const size_t* const tgdim=mxGetDimensions(prhs[9]);
      if (!(tgdim[0]==(unsigned)ntgComp)) throw("The first dimension of input argument 'sg' has incorrect size.");
      if (!(nugdim==ntgdim)) throw("Matrix dimensions of input arguments 'ug' and 'sg' must agree.");
      for (unsigned int iDim=1; iDim<nugdim; iDim++)
      {
        if (!(ugdim[iDim]==tgdim[iDim])) throw("Matrix dimensions of input arguments 'ug' and 'sg' are incompatible.");
      }
    }
  }
  
  // Periodic problems 
  if (probPeriodic){
    if (!mxIsNumeric(prhs[10])) throw("Input argument 'L' must be numeric.");
    if (mxIsSparse(prhs[10])) throw("Input argument 'L' must not be sparse.");
    if (mxIsComplex(prhs[10])) throw("Input argument 'L' must be real.");
    if (!(mxGetNumberOfElements(prhs[10])==1)) throw("Input argument 'L' must be a scalar.");
    
    if (!mxIsNumeric(prhs[11])) throw("Input argument 'ky' must be numeric.");
    if (mxIsSparse(prhs[11])) throw("Input argument 'ky' must not be sparse.");
    if (mxIsComplex(prhs[11])) throw("Input argument 'ky' must be real.");
    if ((mxGetNumberOfDimensions(prhs[11])>2) ||
      ((mxGetM(prhs[11])>1) && (mxGetN(prhs[11])>1)))
                  throw("Input argument 'ky' must be a scalar or a vector.");
    
    if (!mxIsNumeric(prhs[12])) throw("Input argument 'nmax' must be numeric.");
    if (mxIsSparse(prhs[12])) throw("Input argument 'nmax' must not be sparse.");
    if (mxIsComplex(prhs[12])) throw("Input argument 'nmax' must be real.");
    if (!(mxGetNumberOfElements(prhs[12])==1)) throw("Input argument 'nmax' must be a scalar.");
  }
  
  const double L=(probPeriodic ? mxGetScalar(prhs[10]) : -1.0);
  const double* const ky=(probPeriodic ? mxGetPr(prhs[11]) : 0);
  const unsigned int nWave=(probPeriodic ? mxGetNumberOfElements(prhs[11]) : 0);
  const unsigned int nmax=(probPeriodic ? (unsigned int)mxGetScalar(prhs[12]) : 0);
  
  // Create dummy sg array if T is requested and no sg is passed.
  // 
  size_t* const tgdim=new size_t[nugdim];
  tgdim[0]=(unsigned)ntgComp;
  for (unsigned int iDim=1; iDim<nugdim; iDim++) tgdim[iDim]=ugdim[iDim];
  if (sgIn) for (unsigned int iDim=0; iDim<nugdim; iDim++) tgdim[iDim]=(unsigned int)0;
  mxArray* sgdummy=mxCreateNumericArray(nugdim,tgdim,mxDOUBLE_CLASS, mxREAL);    
  delete [] tgdim;
  
  const double* const tgRe=(TmatOut?(sgIn? mxGetPr(prhs[9]):mxGetPr(sgdummy)):0);
  const double* const tgIm=(((TmatOut)&(sgIn))?mxGetPi(prhs[9]):0);
  const bool tgCmplx=(((TmatOut)&(sgIn))?mxIsComplex(prhs[9]):false);
  
  // Number of degrees of freedom points per collocation point.
  unsigned int nColDof;
  if (nugComp==1) nColDof=1;                 // 2D, out-of-plane
  if (nugComp==4) nColDof=2;                 // 2D, in-plane
  if ((nugComp==5) && !(probAxi)) nColDof=3; // 3D
  if ((nugComp==5) &&   probAxi ) nColDof=2; // Axisymmetric
  if (nugComp==9) nColDof=3;                 // 2.5D

  // Vertical receiver coordinate is passed relative
  const bool zRel=false;

  // COPY VARIABLES TO GENERIC ARRAY OF POINTERS GREENPTR
  // The generic pointer has the same layout in both functions
  // bemmat_mex.cpp and bemxfer_mex.cpp
  const unsigned int nGreenPtr=14;
  const unsigned int GreenFunType=1;
  const void** const greenPtr=new(nothrow) const void*[nGreenPtr];
    if (greenPtr==0) throw("Out of memory.");
  greenPtr[0]=&GreenFunType;
  greenPtr[1]=&nzs;
  greenPtr[2]=zs;
  greenPtr[3]=&nr;
  greenPtr[4]=r;
  greenPtr[5]=&nz;
  greenPtr[6]=z;
  greenPtr[7]=ugRe;
  greenPtr[8]=ugIm;
  greenPtr[9]=tgRe;
  greenPtr[10]=tgIm;
  greenPtr[11]=(double* const)0; // tg0Re
  greenPtr[12]=(double* const)0; // tg0Im
  greenPtr[13]=&zRel;
  bemIntegrate(plhs,probAxi,probPeriodic,probDim,nColDof,TmatOut,Nod,nNod,Elt,nElt,Rec,nRec,
               TypeID,TypeName,TypeKeyOpts,nKeyOpt,nEltType,CollPoints,
               nTotalColl,nCentroidColl,greenPtr,nGrSet,nugComp,ugCmplx,tgCmplx,
               greenDim,nGreenDim,L,ky,nWave,nmax);

  mxDestroyArray(sgdummy);               
  delete [] greenPtr;
  delete [] greenDim;
}

//==============================================================================
void IntegrateFsGreenf(mxArray* plhs[], int nrhs,
                       const mxArray* prhs[], const bool probAxi, const bool& probPeriodic,
                       const unsigned int& probDim, const double* const Nod,
                       const unsigned int& nNod, const double* const Elt,
                       const unsigned int& nElt, const double* const Rec, const unsigned int& nRec,
                       const unsigned int* const TypeID,
                       const char* const TypeName[], const char* const TypeKeyOpts[],
                       const unsigned int* const nKeyOpt, const unsigned int& nEltType,
                       const double* const CollPoints, const unsigned int& nTotalColl,
                       const unsigned int& nCentroidColl,
                       const bool& TmatOut)
//==============================================================================
{
  if (!(nrhs==12)) throw("Wrong number of input arguments.");

  if (!mxIsNumeric(prhs[5])) throw("Input argument 'Cs' must be numeric.");
  if (mxIsSparse(prhs[5])) throw("Input argument 'Cs' must not be sparse.");
  if (mxIsComplex(prhs[5])) throw("Input argument 'Cs' must be real.");
  if (!(mxGetNumberOfElements(prhs[5])==1)) throw("Input argument 'Cs' must be a scalar.");

  if (!mxIsNumeric(prhs[6])) throw("Input argument 'Cp' must be numeric.");
  if (mxIsSparse(prhs[6])) throw("Input argument 'Cp' must not be sparse.");
  if (mxIsComplex(prhs[6])) throw("Input argument 'Cp' must be real.");
  if (!(mxGetNumberOfElements(prhs[6])==1)) throw("Input argument 'Cp' must be a scalar.");

  if (!mxIsNumeric(prhs[7])) throw("Input argument 'Ds' must be numeric.");
  if (mxIsSparse(prhs[7])) throw("Input argument 'Ds' must not be sparse.");
  if (mxIsComplex(prhs[7])) throw("Input argument 'Ds' must be real.");
  if (!(mxGetNumberOfElements(prhs[7])==1)) throw("Input argument 'Ds' must be a scalar.");

  if (!mxIsNumeric(prhs[8])) throw("Input argument 'Dp' must be numeric.");
  if (mxIsSparse(prhs[8])) throw("Input argument 'Dp' must not be sparse.");
  if (mxIsComplex(prhs[8])) throw("Input argument 'Dp' must be real.");
  if (!(mxGetNumberOfElements(prhs[8])==1)) throw("Input argument 'Dp' must be a scalar.");

  if (!mxIsNumeric(prhs[9])) throw("Input argument 'rho' must be numeric.");
  if (mxIsSparse(prhs[9])) throw("Input argument 'rho' must not be sparse.");
  if (mxIsComplex(prhs[9])) throw("Input argument 'rho' must be real.");
  if (!(mxGetNumberOfElements(prhs[9])==1)) throw("Input argument 'rho' must be a scalar.");

  if (!mxIsNumeric(prhs[10])) throw("Input argument 'py' must be numeric.");
  if (mxIsSparse(prhs[10])) throw("Input argument 'py' must not be sparse.");
  if (mxIsComplex(prhs[10])) throw("Input argument 'py' must be real.");
  if ((mxGetNumberOfDimensions(prhs[10])>2) ||
      ((mxGetM(prhs[10])>1) && (mxGetN(prhs[10])>1)))
    throw("Input argument 'py' must be a scalar or a vector.");

  if (!mxIsNumeric(prhs[11])) throw("Input argument 'omega' must be numeric.");
  if (mxIsSparse(prhs[11])) throw("Input argument 'omega' must not be sparse.");
  if (mxIsComplex(prhs[11])) throw("Input argument 'omega' must be real.");
  if ((mxGetNumberOfDimensions(prhs[11])>2) ||
      ((mxGetM(prhs[11])>1) && (mxGetN(prhs[11])>1)))
    throw("Input argument 'omega' must be a scalar or a vector.");

  const unsigned int nugComp=9;
  const unsigned int nColDof=3;

  const double Cs  = mxGetScalar(prhs[5]);
  const double Cp  = mxGetScalar(prhs[6]);
  const double Ds  = mxGetScalar(prhs[7]);
  const double Dp  = mxGetScalar(prhs[8]);
  const double rho = mxGetScalar(prhs[9]);
  const double* const py = mxGetPr(prhs[10]);
  const unsigned int nWave = mxGetNumberOfElements(prhs[10]);
  const double* const omega = mxGetPr(prhs[11]);
  const unsigned int nFreq = mxGetNumberOfElements(prhs[11]);

  const unsigned int nGrSet=nFreq*nWave;
  const bool ugCmplx=true;
  const bool tgCmplx=true;

  // OUTPUT ARGUMENT LAST DIMENSIONS
  const unsigned int nGreenDim=2;
  unsigned int* const greenDim=new(nothrow) unsigned int[nGreenDim];
  greenDim[0]=nWave;
  greenDim[1]=nFreq;

  // COPY VARIABLES TO GENERIC ARRAY OF POINTERS GREENPTR
  const unsigned int nGreenPtr=10;
  const unsigned int GreenFunType=2;
  const void** const greenPtr=new(nothrow) const void*[nGreenPtr];
    if (greenPtr==0) throw("Out of memory.");

  greenPtr[0]=&GreenFunType;
  greenPtr[1]=&Cs;
  greenPtr[2]=&Cp;
  greenPtr[3]=&Ds;
  greenPtr[4]=&Dp;
  greenPtr[5]=&rho;
  greenPtr[6]=&nWave;
  greenPtr[7]=&nFreq;
  greenPtr[8]=py;
  greenPtr[9]=omega;

  const double L=-1.0;
  const double* const ky=0;
  const unsigned int nky=0;
  const unsigned int nmax=0;

  bemIntegrate(plhs,probAxi,probPeriodic,probDim,nColDof,TmatOut,Nod,nNod,Elt,nElt,Rec,nRec,
               TypeID,TypeName,TypeKeyOpts,nKeyOpt,nEltType,CollPoints,
               nTotalColl,nCentroidColl,greenPtr,nGrSet,nugComp,ugCmplx,tgCmplx,
               greenDim,nGreenDim,L,ky,nky,nmax);
  delete [] greenPtr;
  delete [] greenDim;
}


//==============================================================================
void IntegrateFsGreen3d(mxArray* plhs[], int nrhs,
                        const mxArray* prhs[], const bool probAxi, const bool& probPeriodic,
                        const unsigned int& probDim, const double* const Nod,
                        const unsigned int& nNod, const double* const Elt,
                        const unsigned int& nElt, const double* const Rec, const unsigned int& nRec,
                        const unsigned int* const TypeID,
                        const char* const TypeName[], const char* const TypeKeyOpts[],
                        const unsigned int* const nKeyOpt, const unsigned int& nEltType,
                        const double* const CollPoints, const unsigned int& nTotalColl,
                        const unsigned int& nCentroidColl,
                        const bool& TmatOut)
//==============================================================================
{
   // INPUT ARGUMENT PROCESSING
  if (probPeriodic)
  {
    if (!(nrhs==14)) throw("Wrong number of input arguments.");
  }
  else
  {
    if (!(nrhs==11)) throw("Wrong number of input arguments.");
  } 

  if (!mxIsNumeric(prhs[5])) throw("Input argument 'Cs' must be numeric.");
  if (mxIsSparse(prhs[5])) throw("Input argument 'Cs' must not be sparse.");
  if (mxIsComplex(prhs[5])) throw("Input argument 'Cs' must be real.");
  if (!(mxGetNumberOfElements(prhs[5])==1)) throw("Input argument 'Cs' must be a scalar.");

  if (!mxIsNumeric(prhs[6])) throw("Input argument 'Cp' must be numeric.");
  if (mxIsSparse(prhs[6])) throw("Input argument 'Cp' must not be sparse.");
  if (mxIsComplex(prhs[6])) throw("Input argument 'Cp' must be real.");
  if (!(mxGetNumberOfElements(prhs[6])==1)) throw("Input argument 'Cp' must be a scalar.");

  if (!mxIsNumeric(prhs[7])) throw("Input argument 'Ds' must be numeric.");
  if (mxIsSparse(prhs[7])) throw("Input argument 'Ds' must not be sparse.");
  if (mxIsComplex(prhs[7])) throw("Input argument 'Ds' must be real.");
  if (!(mxGetNumberOfElements(prhs[7])==1)) throw("Input argument 'Ds' must be a scalar.");

  if (!mxIsNumeric(prhs[8])) throw("Input argument 'Dp' must be numeric.");
  if (mxIsSparse(prhs[8])) throw("Input argument 'Dp' must not be sparse.");
  if (mxIsComplex(prhs[8])) throw("Input argument 'Dp' must be real.");
  if (!(mxGetNumberOfElements(prhs[8])==1)) throw("Input argument 'Dp' must be a scalar.");

  if (!mxIsNumeric(prhs[9])) throw("Input argument 'rho' must be numeric.");
  if (mxIsSparse(prhs[9])) throw("Input argument 'rho' must not be sparse.");
  if (mxIsComplex(prhs[9])) throw("Input argument 'rho' must be real.");
  if (!(mxGetNumberOfElements(prhs[9])==1)) throw("Input argument 'rho' must be a scalar.");

  if (!mxIsNumeric(prhs[10])) throw("Input argument 'omega' must be numeric.");
  if (mxIsSparse(prhs[10])) throw("Input argument 'omega' must not be sparse.");
  if (mxIsComplex(prhs[10])) throw("Input argument 'omega' must be real.");
  if ((mxGetNumberOfDimensions(prhs[10])>2) ||
      ((mxGetM(prhs[10])>1) && (mxGetN(prhs[10])>1)))
    throw("Input argument 'omega' must be a scalar or a vector.");

  const unsigned int nugComp=5;
  const unsigned int nColDof=(probAxi ? 2 : 3);

  const double Cs  = mxGetScalar(prhs[5]);
  const double Cp  = mxGetScalar(prhs[6]);
  const double Ds  = mxGetScalar(prhs[7]);
  const double Dp  = mxGetScalar(prhs[8]);
  const double rho = mxGetScalar(prhs[9]);
  const unsigned int nFreq = mxGetNumberOfElements(prhs[10]);
  const double* const omega = mxGetPr(prhs[10]);

  const unsigned int nGrSet=nFreq;
  const bool ugCmplx=true;
  const bool tgCmplx=true;

  // OUTPUT ARGUMENT LAST DIMENSIONS
  const unsigned int nGreenDim=1;
  unsigned int* const greenDim=new(nothrow) unsigned int[nGreenDim];
  greenDim[0]=nFreq;

  // COPY VARIABLES TO GENERIC ARRAY OF POINTERS GREENPTR
  const unsigned int nGreenPtr=8;
  const unsigned int GreenFunType=3;
  const void** const greenPtr=new(nothrow) const void*[nGreenPtr];
    if (greenPtr==0) throw("Out of memory.");

  greenPtr[0]=&GreenFunType;
  greenPtr[1]=&Cs;
  greenPtr[2]=&Cp;
  greenPtr[3]=&Ds;
  greenPtr[4]=&Dp;
  greenPtr[5]=&rho;
  greenPtr[6]=&nFreq;
  greenPtr[7]=omega;

  // Periodic problems 
  if (probPeriodic){
    if (!mxIsNumeric(prhs[11])) throw("Input argument 'L' must be numeric.");
    if (mxIsSparse(prhs[11])) throw("Input argument 'L' must not be sparse.");
    if (mxIsComplex(prhs[11])) throw("Input argument 'L' must be real.");
    if (!(mxGetNumberOfElements(prhs[11])==1)) throw("Input argument 'L' must be a scalar.");
    
    if (!mxIsNumeric(prhs[12])) throw("Input argument 'ky' must be numeric.");
    if (mxIsSparse(prhs[12])) throw("Input argument 'ky' must not be sparse.");
    if (mxIsComplex(prhs[12])) throw("Input argument 'ky' must be real.");
    if ((mxGetNumberOfDimensions(prhs[12])>2) ||
      ((mxGetM(prhs[12])>1) && (mxGetN(prhs[12])>1)))
                  throw("Input argument 'ky' must be a scalar or a vector.");
                  
    if (!mxIsNumeric(prhs[13])) throw("Input argument 'nmax' must be numeric.");
    if (mxIsSparse(prhs[13])) throw("Input argument 'nmax' must not be sparse.");
    if (mxIsComplex(prhs[13])) throw("Input argument 'nmax' must be real.");
    if (!(mxGetNumberOfElements(prhs[13])==1)) throw("Input argument 'nmax' must be a scalar.");
  }
  const double L=(probPeriodic ? mxGetScalar(prhs[11]) : -1.0);
  const double* const ky=(probPeriodic ? mxGetPr(prhs[12]) : 0);
  const unsigned int nWave=(probPeriodic ? mxGetNumberOfElements(prhs[12]) : 0);
  const unsigned int nmax=(probPeriodic ? (unsigned int)mxGetScalar(prhs[13]) : 0);
  
  bemIntegrate(plhs,probAxi,probPeriodic,probDim,nColDof,TmatOut,Nod,nNod,Elt,nElt,Rec,nRec,
               TypeID,TypeName,TypeKeyOpts,nKeyOpt,nEltType,CollPoints,
               nTotalColl,nCentroidColl,greenPtr,nGrSet,nugComp,ugCmplx,tgCmplx,
               greenDim,nGreenDim,L,ky,nWave,nmax);
  delete [] greenPtr;
  delete [] greenDim;
}

//==============================================================================
void IntegrateFsGreen3d0(mxArray* plhs[], int nrhs,
                         const mxArray* prhs[], const bool probAxi, const bool& probPeriodic,
                         const unsigned int& probDim, const double* const Nod,
                         const unsigned int& nNod, const double* const Elt,
                         const unsigned int& nElt, const double* const Rec, const unsigned int& nRec,
                         const unsigned int* const TypeID,
                         const char* const TypeName[], const char* const TypeKeyOpts[],
                         const unsigned int* const nKeyOpt, const unsigned int& nEltType,
                         const double* const CollPoints, const unsigned int& nTotalColl,
                         const unsigned int& nCentroidColl,
                         const bool& TmatOut)
/* Initialize Green's function for user defined Green's function ('USER')
 *
 *    greenPtr[0]=&GreenFunType;   Green's function type identifier
 *    greenPtr[1]=&nGrSet;         Number of function sets (nFreq,nTime,nCase)
 *    greenPtr[...]                different definition for different
 *                                 Green's function types
 *
 */
//==============================================================================
{
  // INPUT ARGUMENT PROCESSING
  if (!(nrhs==7)) throw("Wrong number of input arguments.");

  if (!mxIsNumeric(prhs[5])) throw("Input argument 'E' must be numeric.");
  if (mxIsSparse(prhs[5])) throw("Input argument 'E' must not be sparse.");
  if (mxIsComplex(prhs[5])) throw("Input argument 'E' must be real.");
  if (!(mxGetNumberOfElements(prhs[5])==1)) throw("Input argument 'E' must be a scalar.");

  if (!mxIsNumeric(prhs[6])) throw("Input argument 'nu' must be numeric.");
  if (mxIsSparse(prhs[6])) throw("Input argument 'nu' must not be sparse.");
  if (mxIsComplex(prhs[6])) throw("Input argument 'nu' must be real.");
  if (!(mxGetNumberOfElements(prhs[6])==1)) throw("Input argument 'nu' must be a scalar.");

  const unsigned int nugComp=5;
  const unsigned int nColDof=(probAxi ? 2 : 3);
  
  const double E  = mxGetScalar(prhs[5]);
  const double nu  = mxGetScalar(prhs[6]);

  const double mu=0.5*E/(1.0+nu);
  const double M=E*(1.0-nu)/(1.0+nu)/(1.0-2.0*nu);
  const double rho = 1.0;

  const double Cs  = sqrt(mu/rho);
  const double Cp  = sqrt(M/rho);
  const double Ds  = 0.0;
  const double Dp  = 0.0;
  const unsigned int nFreq  = 1;
  double* const omega = new(nothrow) double[nFreq];
  omega[0]=0.0;

  const unsigned int nGrSet=nFreq;
  const bool ugCmplx=false;
  const bool tgCmplx=false;

  // OUTPUT ARGUMENT LAST DIMENSIONS
  const unsigned int nGreenDim=1;
  unsigned int* const greenDim=new(nothrow) unsigned int[nGreenDim];
  greenDim[0]=nFreq;

  // COPY VARIABLES TO GENERIC ARRAY OF POINTERS GREENPTR
  const unsigned int nGreenPtr=8;
  const unsigned int GreenFunType=3;
  const void** const greenPtr=new(nothrow) const void*[nGreenPtr];
    if (greenPtr==0) throw("Out of memory.");

  greenPtr[0]=&GreenFunType;
  greenPtr[1]=&Cs;
  greenPtr[2]=&Cp;
  greenPtr[3]=&Ds;
  greenPtr[4]=&Dp;
  greenPtr[5]=&rho;
  greenPtr[6]=&nFreq;
  greenPtr[7]=omega;

  const double L=-1.0;
  const double* const ky=0;
  const unsigned int nWave=0;
  const unsigned int nmax=0;

  bemIntegrate(plhs,probAxi,probPeriodic,probDim,nColDof,TmatOut,Nod,nNod,Elt,nElt,Rec,nRec,
               TypeID,TypeName,TypeKeyOpts,nKeyOpt,nEltType,CollPoints,
               nTotalColl,nCentroidColl,greenPtr,nGrSet,nugComp,ugCmplx,tgCmplx,
               greenDim,nGreenDim,L,ky,nWave,nmax);
  delete [] greenPtr;
  delete [] greenDim;
  delete [] omega;
}

//==============================================================================
void IntegrateFsGreen2d_inplane(mxArray* plhs[], int nrhs,
                                const mxArray* prhs[], const bool probAxi, const bool& probPeriodic,
                                const unsigned int& probDim, const double* const Nod,
                                const unsigned int& nNod, const double* const Elt,
                                const unsigned int& nElt, const double* const Rec, const unsigned int& nRec,
                                const unsigned int* const TypeID,
                                const char* const TypeName[], const char* const TypeKeyOpts[],
                                const unsigned int* const nKeyOpt, const unsigned int& nEltType,
                                const double* const CollPoints, const unsigned int& nTotalColl,
                                const unsigned int& nCentroidColl,
                                const bool& TmatOut)
/* Initialize Green's function for user defined Green's function ('USER')
 *
 *    greenPtr[0]=&GreenFunType;   Green's function type identifier
 *    greenPtr[1]=&nGrSet;         Number of function sets (nFreq,nTime,nCase)
 *    greenPtr[...]                different definition for different
 *                                 Green's function types
 */
//==============================================================================
{
  // INPUT ARGUMENT PROCESSING
  if (!(nrhs==11)) throw("Wrong number of input arguments.");

  if (!mxIsNumeric(prhs[5])) throw("Input argument 'Cs' must be numeric.");
  if (mxIsSparse(prhs[5])) throw("Input argument 'Cs' must not be sparse.");
  if (mxIsComplex(prhs[5])) throw("Input argument 'Cs' must be real.");
  if (!(mxGetNumberOfElements(prhs[5])==1)) throw("Input argument 'Cs' must be a scalar.");

  if (!mxIsNumeric(prhs[6])) throw("Input argument 'Cp' must be numeric.");
  if (mxIsSparse(prhs[6])) throw("Input argument 'Cp' must not be sparse.");
  if (mxIsComplex(prhs[6])) throw("Input argument 'Cp' must be real.");
  if (!(mxGetNumberOfElements(prhs[6])==1)) throw("Input argument 'Cp' must be a scalar.");

  if (!mxIsNumeric(prhs[7])) throw("Input argument 'Ds' must be numeric.");
  if (mxIsSparse(prhs[7])) throw("Input argument 'Ds' must not be sparse.");
  if (mxIsComplex(prhs[7])) throw("Input argument 'Ds' must be real.");
  if (!(mxGetNumberOfElements(prhs[7])==1)) throw("Input argument 'Ds' must be a scalar.");

  if (!mxIsNumeric(prhs[8])) throw("Input argument 'Dp' must be numeric.");
  if (mxIsSparse(prhs[8])) throw("Input argument 'Dp' must not be sparse.");
  if (mxIsComplex(prhs[8])) throw("Input argument 'Dp' must be real.");
  if (!(mxGetNumberOfElements(prhs[8])==1)) throw("Input argument 'Dp' must be a scalar.");

  if (!mxIsNumeric(prhs[9])) throw("Input argument 'rho' must be numeric.");
  if (mxIsSparse(prhs[9])) throw("Input argument 'rho' must not be sparse.");
  if (mxIsComplex(prhs[9])) throw("Input argument 'rho' must be real.");
  if (!(mxGetNumberOfElements(prhs[9])==1)) throw("Input argument 'rho' must be a scalar.");

  if (!mxIsNumeric(prhs[10])) throw("Input argument 'omega' must be numeric.");
  if (mxIsSparse(prhs[10])) throw("Input argument 'omega' must not be sparse.");
  if (mxIsComplex(prhs[10])) throw("Input argument 'omega' must be real.");
  if ((mxGetNumberOfDimensions(prhs[10])>2) ||
      ((mxGetM(prhs[10])>1) && (mxGetN(prhs[10])>1)))
    throw("Input argument 'omega' must be a scalar or a vector.");

  const unsigned int nugComp=4;
  const unsigned int nColDof=2;

  const double Cs  = mxGetScalar(prhs[5]);
  const double Cp  = mxGetScalar(prhs[6]);
  const double Ds  = mxGetScalar(prhs[7]);
  const double Dp  = mxGetScalar(prhs[8]);
  const double rho = mxGetScalar(prhs[9]);
  const double* const omega = mxGetPr(prhs[10]);
  const unsigned int nFreq = mxGetNumberOfElements(prhs[10]);

  const unsigned int nGrSet=nFreq;
  const bool ugCmplx=true;
  const bool tgCmplx=true;

  // OUTPUT ARGUMENT LAST DIMENSIONS
  const unsigned int nGreenDim=1;
  unsigned int* const greenDim=new(nothrow) unsigned int[nGreenDim];
  greenDim[0]=nFreq;

  // COPY VARIABLES TO GENERIC ARRAY OF POINTERS GREENPTR
  const unsigned int nGreenPtr=8;
  const unsigned int GreenFunType=4;
  const void** const greenPtr=new(nothrow) const void*[nGreenPtr];
    if (greenPtr==0) throw("Out of memory.");

  greenPtr[0]=&GreenFunType;
  greenPtr[1]=&Cs;
  greenPtr[2]=&Cp;
  greenPtr[3]=&Ds;
  greenPtr[4]=&Dp;
  greenPtr[5]=&rho;
  greenPtr[6]=&nFreq;
  greenPtr[7]=omega;

  const double L=-1.0;
  const double* const ky=0;
  const unsigned int nWave=0;
  const unsigned int nmax=0;

  bemIntegrate(plhs,probAxi,probPeriodic,probDim,nColDof,TmatOut,Nod,nNod,Elt,nElt,Rec,nRec,
               TypeID,TypeName,TypeKeyOpts,nKeyOpt,nEltType,CollPoints,
               nTotalColl,nCentroidColl,greenPtr,nGrSet,nugComp,ugCmplx,tgCmplx,
               greenDim,nGreenDim,L,ky,nWave,nmax);
  delete [] greenPtr;
  delete [] greenDim;
}

//==============================================================================
void IntegrateFsGreen2d_inplane0(mxArray* plhs[], int nrhs,
                                 const mxArray* prhs[], const bool probAxi, const bool& probPeriodic,
                                 const unsigned int& probDim, const double* const Nod,
                                 const unsigned int& nNod, const double* const Elt,
                                 const unsigned int& nElt, const double* const Rec, const unsigned int& nRec,
                                 const unsigned int* const TypeID,
                                 const char* const TypeName[], const char* const TypeKeyOpts[],
                                 const unsigned int* const nKeyOpt, const unsigned int& nEltType,
                                 const double* const CollPoints, const unsigned int& nTotalColl,
                                 const unsigned int& nCentroidColl,
                                 const bool& TmatOut)
/* Initialize Green's function for user defined Green's function ('FSGREEN2D_INPLANE0')
 *
 *    greenPtr[0]=&GreenFunType;   Green's function type identifier
 *    greenPtr[1]=&nGrSet;         Number of function sets (nFreq,nTime,nCase)
 *    greenPtr[...]                different definition for different
 *                                 Green's function types
 *
 */
//==============================================================================
{
  // INPUT ARGUMENT PROCESSING
  if (!(nrhs==7)) throw("Wrong number of input arguments.");

  if (!mxIsNumeric(prhs[5])) throw("Input argument 'E' must be numeric.");
  if (mxIsSparse(prhs[5])) throw("Input argument 'E' must not be sparse.");
  if (mxIsComplex(prhs[5])) throw("Input argument 'E' must be real.");
  if (!(mxGetNumberOfElements(prhs[5])==1)) throw("Input argument 'E' must be a scalar.");

  if (!mxIsNumeric(prhs[6])) throw("Input argument 'nu' must be numeric.");
  if (mxIsSparse(prhs[6])) throw("Input argument 'nu' must not be sparse.");
  if (mxIsComplex(prhs[6])) throw("Input argument 'nu' must be real.");
  if (!(mxGetNumberOfElements(prhs[6])==1)) throw("Input argument 'nu' must be a scalar.");

  const unsigned int nugComp=4;
  const unsigned int nColDof=2;

  const double E  = mxGetScalar(prhs[5]);
  const double nu  = mxGetScalar(prhs[6]);

  const double mu=0.5*E/(1.0+nu);
  const double M=E*(1.0-nu)/(1.0+nu)/(1.0-2.0*nu);
  const double rho = 1.0;

  const double Cs  = sqrt(mu/rho);
  const double Cp  = sqrt(M/rho);
  const double Ds  = 0.0;
  const double Dp  = 0.0;
  const unsigned int nFreq  = 1;
  double* const omega = new(nothrow) double[nFreq];
  omega[0]=0.0;

  const unsigned int nGrSet=nFreq;
  const bool ugCmplx=false;
  const bool tgCmplx=false;

  // OUTPUT ARGUMENT LAST DIMENSIONS
  const unsigned int nGreenDim=1;
  unsigned int* const greenDim=new(nothrow) unsigned int[nGreenDim];
  greenDim[0]=nFreq;

  // COPY VARIABLES TO GENERIC ARRAY OF POINTERS GREENPTR
  const unsigned int nGreenPtr=8;
  const unsigned int GreenFunType=4;
  const void** const greenPtr=new(nothrow) const void*[nGreenPtr];
    if (greenPtr==0) throw("Out of memory.");

  greenPtr[0]=&GreenFunType;
  greenPtr[1]=&Cs;
  greenPtr[2]=&Cp;
  greenPtr[3]=&Ds;
  greenPtr[4]=&Dp;
  greenPtr[5]=&rho;
  greenPtr[6]=&nFreq;
  greenPtr[7]=omega;

  const double L=-1.0;
  const double* const ky=0;
  const unsigned int nWave=0;
  const unsigned int nmax=0;

  bemIntegrate(plhs,probAxi,probPeriodic,probDim,nColDof,TmatOut,Nod,nNod,Elt,nElt,Rec,nRec,
               TypeID,TypeName,TypeKeyOpts,nKeyOpt,nEltType,CollPoints,
               nTotalColl,nCentroidColl,greenPtr,nGrSet,nugComp,ugCmplx,tgCmplx,
               greenDim,nGreenDim,L,ky,nWave,nmax);
  delete [] greenPtr;
  delete [] greenDim;
  delete [] omega;
}

//==============================================================================
void IntegrateFsGreen2d_outofplane(mxArray* plhs[], int nrhs,
                                   const mxArray* prhs[], const bool probAxi, const bool& probPeriodic,
                                   const unsigned int& probDim, const double* const Nod,
                                   const unsigned int& nNod, const double* const Elt,
                                   const unsigned int& nElt, const double* const Rec, const unsigned int& nRec,
                                   const unsigned int* const TypeID,
                                   const char* const TypeName[], const char* const TypeKeyOpts[],
                                   const unsigned int* const nKeyOpt, const unsigned int& nEltType,
                                   const double* const CollPoints, const unsigned int& nTotalColl,
                                   const unsigned int& nCentroidColl,
                                   const bool& TmatOut)
/* Initialize Green's function for user defined Green's function ('USER')
 *
 *    greenPtr[0]=&GreenFunType;   Green's function type identifier
 *    greenPtr[1]=&nGrSet;         Number of function sets (nFreq,nTime,nCase)
 *    greenPtr[...]                different definition for different
 *                                 Green's function types
 *
 */
//==============================================================================
{
  // INPUT ARGUMENT PROCESSING
  if (!(nrhs==9)) throw("Wrong number of input arguments.");

  if (!mxIsNumeric(prhs[5])) throw("Input argument 'Cs' must be numeric.");
  if (mxIsSparse(prhs[5])) throw("Input argument 'Cs' must not be sparse.");
  if (mxIsComplex(prhs[5])) throw("Input argument 'Cs' must be real.");
  if (!(mxGetNumberOfElements(prhs[5])==1)) throw("Input argument 'Cs' must be a scalar.");

  if (!mxIsNumeric(prhs[6])) throw("Input argument 'Ds' must be numeric.");
  if (mxIsSparse(prhs[6])) throw("Input argument 'Ds' must not be sparse.");
  if (mxIsComplex(prhs[6])) throw("Input argument 'Ds' must be real.");
  if (!(mxGetNumberOfElements(prhs[6])==1)) throw("Input argument 'Ds' must be a scalar.");

  if (!mxIsNumeric(prhs[7])) throw("Input argument 'rho' must be numeric.");
  if (mxIsSparse(prhs[7])) throw("Input argument 'rho' must not be sparse.");
  if (mxIsComplex(prhs[7])) throw("Input argument 'rho' must be real.");
  if (!(mxGetNumberOfElements(prhs[7])==1)) throw("Input argument 'rho' must be a scalar.");

  if (!mxIsNumeric(prhs[8])) throw("Input argument 'omega' must be numeric.");
  if (mxIsSparse(prhs[8])) throw("Input argument 'omega' must not be sparse.");
  if (mxIsComplex(prhs[8])) throw("Input argument 'omega' must be real.");
  if ((mxGetNumberOfDimensions(prhs[8])>2) ||
      ((mxGetM(prhs[8])>1) && (mxGetN(prhs[8])>1)))
    throw("Input argument 'omega' must be a scalar or a vector.");

  const unsigned int nugComp=1;
  const unsigned int nColDof=1;

  const double Cs  = mxGetScalar(prhs[5]);
  const double Ds  = mxGetScalar(prhs[6]);
  const double rho = mxGetScalar(prhs[7]);
  const double* const omega = mxGetPr(prhs[8]);
  const unsigned int nFreq = mxGetNumberOfElements(prhs[8]);

  const unsigned int nGrSet=nFreq;
  const bool ugCmplx=true;
  const bool tgCmplx=true;

  // OUTPUT ARGUMENT LAST DIMENSIONS
  const unsigned int nGreenDim=1;
  unsigned int* const greenDim=new(nothrow) unsigned int[nGreenDim];
  greenDim[0]=nFreq;

  // COPY VARIABLES TO GENERIC ARRAY OF POINTERS GREENPTR
  const unsigned int nGreenPtr=6;
  const unsigned int GreenFunType=5;
  const void** const greenPtr=new(nothrow) const void*[nGreenPtr];
    if (greenPtr==0) throw("Out of memory.");

  greenPtr[0]=&GreenFunType;
  greenPtr[1]=&Cs;
  greenPtr[2]=&Ds;
  greenPtr[3]=&rho;
  greenPtr[4]=&nFreq;
  greenPtr[5]=omega;
  
  const double L=-1.0;
  const double* const ky=0;
  const unsigned int nWave=0;
  const unsigned int nmax=0;
  
  bemIntegrate(plhs,probAxi,probPeriodic,probDim,nColDof,TmatOut,Nod,nNod,Elt,nElt,Rec,nRec,
               TypeID,TypeName,TypeKeyOpts,nKeyOpt,nEltType,CollPoints,
               nTotalColl,nCentroidColl,greenPtr,nGrSet,nugComp,ugCmplx,tgCmplx,
               greenDim,nGreenDim,L,ky,nWave,nmax);
  delete [] greenPtr;
  delete [] greenDim;
}

//==============================================================================
void IntegrateFsGreen2d_outofplane0(mxArray* plhs[], int nrhs,
                                    const mxArray* prhs[], const bool probAxi, const bool& probPeriodic,
                                    const unsigned int& probDim, const double* const Nod,
                                    const unsigned int& nNod, const double* const Elt,
                                    const unsigned int& nElt, const double* const Rec, const unsigned int& nRec,
                                    const unsigned int* const TypeID,
                                    const char* const TypeName[], const char* const TypeKeyOpts[],
                                    const unsigned int* const nKeyOpt, const unsigned int& nEltType,
                                    const double* const CollPoints, const unsigned int& nTotalColl,
                                    const unsigned int& nCentroidColl,
                                    const bool& TmatOut)
/* Initialize Green's function for user defined Green's function ('IntegrateFsGreen2d_outofplane0')
 *
 *
 *    greenPtr[0]=&GreenFunType;   Green's function type identifier
 *    greenPtr[1]=&nGrSet;         Number of function sets (nFreq,nTime,nCase)
 *    greenPtr[...]                different definition for different
 *                                 Green's function types
 *
 */
//==============================================================================
{
  // INPUT ARGUMENT PROCESSING
  if (!(nrhs==6)) throw("Wrong number of input arguments.");

  if (!mxIsNumeric(prhs[5])) throw("Input argument 'mu' must be numeric.");
  if (mxIsSparse(prhs[5])) throw("Input argument 'mu' must not be sparse.");
  if (mxIsComplex(prhs[5])) throw("Input argument 'mu' must be real.");
  if (!(mxGetNumberOfElements(prhs[5])==1)) throw("Input argument 'mu' must be a scalar.");

  const unsigned int nugComp=1;
  const unsigned int nColDof=1;

  const double mu=mxGetScalar(prhs[5]);
  const double rho = 1.0;
  const double Cs  = sqrt(mu/rho);
  const double Ds  = 0.0;
  const unsigned int nFreq = 1;
  double* const omega = new(nothrow) double[1];
  omega[0]=0.0;

  const unsigned int nGrSet=nFreq;
  const bool ugCmplx=false;
  const bool tgCmplx=false;

  // OUTPUT ARGUMENT LAST DIMENSIONS
  const unsigned int nGreenDim=1;
  unsigned int* const greenDim=new(nothrow) unsigned int[nGreenDim];
  greenDim[0]=nFreq;

  // COPY VARIABLES TO GENERIC ARRAY OF POINTERS GREENPTR
  const unsigned int nGreenPtr=6;
  const unsigned int GreenFunType=5;
  const void** const greenPtr=new(nothrow) const void*[nGreenPtr];
    if (greenPtr==0) throw("Out of memory.");

  greenPtr[0]=&GreenFunType;
  greenPtr[1]=&Cs;
  greenPtr[2]=&Ds;
  greenPtr[3]=&rho;
  greenPtr[4]=&nFreq;
  greenPtr[5]=omega;

  const double L=-1.0;
  const double* const ky=0;
  const unsigned int nWave=0;
  const unsigned int nmax=0;
  
  bemIntegrate(plhs,probAxi,probPeriodic,probDim,nColDof,TmatOut,Nod,nNod,Elt,nElt,Rec,nRec,
               TypeID,TypeName,TypeKeyOpts,nKeyOpt,nEltType,CollPoints,
               nTotalColl,nCentroidColl,greenPtr,nGrSet,nugComp,ugCmplx,tgCmplx,
               greenDim,nGreenDim,L,ky,nWave,nmax);
  delete [] greenPtr;
  delete [] greenDim;
  delete [] omega;
}

//==============================================================================
void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
//==============================================================================
{
  try
  {
    checklicense();

    // INPUT ARGUMENT PROCESSING
    if (nrhs<5) throw("Not enough input arguments.");
    if (nlhs>2) throw("Too many output arguments.");

    if (!mxIsNumeric(prhs[0])) throw("Input argument 'nod' must be numeric.");
    if (mxIsSparse(prhs[0])) throw("Input argument 'nod' must not be sparse.");
    if (mxIsComplex(prhs[0])) throw("Input argument 'nod' must be real.");
    if (!(mxGetN(prhs[0])==4)) throw("Input argument 'nod' should have 4 columns.");
    const unsigned int nNod=mxGetM(prhs[0]);
    const double* const Nod=mxGetPr(prhs[0]);

    if (!mxIsNumeric(prhs[1])) throw("Input argument 'elt' must be numeric.");
    if (mxIsSparse(prhs[1])) throw("Input argument 'elt' must not be sparse.");
    if (mxIsComplex(prhs[1])) throw("Input argument 'elt' must be real.");
    if (mxGetN(prhs[1])<=3) throw("Input argument 'elt' should have at least 3 columns.");
    const double* const Elt=mxGetPr(prhs[1]);
    const unsigned int nElt=mxGetM(prhs[1]);
    const unsigned int maxEltColumn=mxGetN(prhs[1]);

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
    const unsigned int probDim=bemDimension(Elt,nElt,TypeID,TypeName,TypeKeyOpts,nKeyOpt,nEltType);
    const bool probAxi=isAxisym(Elt,nElt,TypeID,TypeName,TypeKeyOpts,nKeyOpt,nEltType);
    const bool probPeriodic=isPeriodic(Elt,nElt,TypeID,TypeName,TypeKeyOpts,nKeyOpt,nEltType);

    if (!mxIsNumeric(prhs[3])) throw("Input argument 'Rec' must be numeric.");
    if (mxIsSparse(prhs[3])) throw("Input argument 'Rec' must not be sparse.");
    if (mxIsComplex(prhs[3])) throw("Input argument 'Rec' must be real.");
    if (!(mxGetN(prhs[3])==3)) throw("Input argument 'Rec' should have 3 columns.");
    const double* const Rec=mxGetPr(prhs[3]);
    const unsigned int nRec=mxGetM(prhs[3]);

    const bool TmatOut=(nlhs>1);

    // COLLOCATION POINTS: NODAL OR CENTROID
    unsigned int* const NodalColl=new(nothrow) unsigned int[nNod];
    if (NodalColl==0) throw("Out of memory.");
    unsigned int* const CentroidColl=new(nothrow) unsigned int[nElt];
    if (CentroidColl==0) throw("Out of memory.");
    unsigned int nCentroidColl;
    unsigned int nNodalColl;
    BemCollPoints(Elt,Nod,TypeID,nKeyOpt,TypeName,TypeKeyOpts,nEltType,
                   nElt,maxEltColumn,nNod,NodalColl,CentroidColl,nNodalColl,nCentroidColl);

    // COLLOCATION POINT COORDINATES
    unsigned int nTotalColl = nNodalColl + nCentroidColl;
    double* const CollPoints=new(nothrow) double[5*nTotalColl];
    if (CollPoints==0) throw("Out of memory.");
    BemCollCoords(Elt,Nod,TypeID,nKeyOpt,TypeName,TypeKeyOpts,nEltType,
                  CentroidColl,NodalColl,CollPoints,nTotalColl,nElt,nNod);


    // INTEGRATE GREEN'S FUNCTION
    if (!mxIsChar(prhs[4])) throw("Input argument 'green' must be a string.");
    const char* const green = mxArrayToString(prhs[4]);

    if (strcasecmp(green,"user")==0)
    {
      IntegrateGreenUser(plhs,nrhs,prhs,probAxi,probPeriodic,probDim,Nod,nNod,Elt,nElt,
                         Rec,nRec,TypeID,TypeName,TypeKeyOpts,nKeyOpt,
                         nEltType,CollPoints,nTotalColl,nCentroidColl,TmatOut);
    }
    else if (strcasecmp(green,"fsgreenf")==0)
    {
      IntegrateFsGreenf(plhs,nrhs,prhs,probAxi,probPeriodic,probDim,Nod,nNod,Elt,nElt,
                        Rec,nRec,TypeID,TypeName,TypeKeyOpts,nKeyOpt,
                        nEltType,CollPoints,nTotalColl,nCentroidColl,TmatOut);
    }
    else if (strcasecmp(green,"fsgreen2d_inplane")==0)
    {
      IntegrateFsGreen2d_inplane(plhs,nrhs,prhs,probAxi,probPeriodic,probDim,Nod,nNod,Elt,nElt,
                                 Rec,nRec,TypeID,TypeName,TypeKeyOpts,nKeyOpt,
                                 nEltType,CollPoints,nTotalColl,nCentroidColl,TmatOut);
    }
    else if (strcasecmp(green,"fsgreen2d_inplane0")==0)
    {
      IntegrateFsGreen2d_inplane0(plhs,nrhs,prhs,probAxi,probPeriodic,probDim,Nod,nNod,Elt,nElt,
                                  Rec,nRec,TypeID,TypeName,TypeKeyOpts,nKeyOpt,
                                  nEltType,CollPoints,nTotalColl,nCentroidColl,TmatOut);
    }
    else if (strcasecmp(green,"fsgreen2d_outofplane")==0)
    {
      IntegrateFsGreen2d_outofplane(plhs,nrhs,prhs,probAxi,probPeriodic,probDim,Nod,nNod,Elt,nElt,
                                    Rec,nRec,TypeID,TypeName,TypeKeyOpts,nKeyOpt,
                                    nEltType,CollPoints,nTotalColl,nCentroidColl,TmatOut);
    }
    else if (strcasecmp(green,"fsgreen2d_outofplane0")==0)
    {
      IntegrateFsGreen2d_outofplane0(plhs,nrhs,prhs,probAxi,probPeriodic,probDim,Nod,nNod,Elt,nElt,
                                     Rec,nRec,TypeID,TypeName,TypeKeyOpts,nKeyOpt,
                                     nEltType,CollPoints,nTotalColl,nCentroidColl,TmatOut);
    }
    else if (strcasecmp(green,"fsgreen3d")==0)
    {
      IntegrateFsGreen3d(plhs,nrhs,prhs,probAxi,probPeriodic,probDim,Nod,nNod,Elt,nElt,
                         Rec,nRec,TypeID,TypeName,TypeKeyOpts,nKeyOpt,
                         nEltType,CollPoints,nTotalColl,nCentroidColl,TmatOut);
    }
    else if (strcasecmp(green,"fsgreen3d0")==0)
    {
      IntegrateFsGreen3d0(plhs,nrhs,prhs,probAxi,probPeriodic,probDim,Nod,nNod,Elt,nElt,
                          Rec,nRec,TypeID,TypeName,TypeKeyOpts,nKeyOpt,
                          nEltType,CollPoints,nTotalColl,nCentroidColl,TmatOut);
    }
    else
    {
      throw("Unknown fundamental solution type for input argument 'green'.");
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
  }
  catch (const char* exception)
  {
    mexErrMsgTxt(exception);
  }
}
