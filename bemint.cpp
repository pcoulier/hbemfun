#include "eltdef.h"
#include "bemcollpoints.h"
#include "gausspw.h"
#include "shapefun.h"
#include "mex.h"
#include <complex>
using namespace std;

void bemint(const double* const Nod, const unsigned int& nNod,
            const double* const Elt, const unsigned int& iElt,const unsigned int& nElt,
            unsigned int* const  TypeID, unsigned int* const  nKeyOpt,
            char* TypeName[], char* TypeKeyOpts[],const unsigned int& nEltType,
            const unsigned int* const eltCollIndex, complex<double>* const OutMat,
            const unsigned int& nDof, const complex<double>* const t,
            const complex<double>* const u, const unsigned int& ntMode, const unsigned int& nuMode,
            const unsigned int& ntSet, const unsigned int& nColDof, bool probAxi)
{
  // Read element properties.
  const unsigned int EltType = (unsigned int)(Elt[nElt+iElt]);
  unsigned int Parent;
  unsigned int nEltNod;
  unsigned int nEltColl;
  unsigned int ShapeTypeN;
  unsigned int ShapeTypeM;
  unsigned int EltDim;
  unsigned int AxiSym;
  unsigned int Periodic;
  unsigned int nGauss0;
  unsigned int nEltDiv;
  unsigned int nGaussSing;
  unsigned int nEltDivSing;
  
  eltdef(EltType,TypeID,TypeName,TypeKeyOpts,nKeyOpt,nEltType,Parent,
          nEltNod,nEltColl,ShapeTypeN,ShapeTypeM,EltDim,AxiSym,Periodic,nGauss0,nEltDiv,
          nGaussSing,nEltDivSing);

  // Number of Gaussian points
  unsigned int nXi;
  unsigned int nGauss;
  if (Parent == 0)
  {
    nGauss=7;
    nXi=7;
  }
  else if (Parent == 1)
  {
    nGauss=7;
    nXi=7;
  }
  else if (Parent == 2)
  {
    nGauss=4;
    nXi=16;
  }

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

  // Determine sample points for the element type.
  double* const xi=new(nothrow) double[2*nXi];
  if (xi==0) throw("Out of memory.");
  double* const H=new(nothrow) double[nXi];
  if (H==0) throw("Out of memory.");

  if (Parent == 0) gausspw1D(1,nGauss,xi,H);
  else if (Parent == 1) gausspwtri(nGauss,xi,H);
  else gausspw2D(1,nGauss,xi,H);

  // Shape functions in the sample points.
  double* const N=new(nothrow) double[nXi*nEltNod];
    if (N==0) throw("Out of memory.");
  double* const M=new(nothrow) double[nXi*nEltColl];
    if (M==0) throw("Out of memory.");
  double* const dN=new(nothrow) double[2*nXi*nEltNod];
    if (dN==0) throw("Out of memory.");
  double* const nat=new(nothrow) double[6*nXi];
  if (nat==0) throw("Out of memory.");
  double* const Jac=new(nothrow) double[nXi];
  if (Jac==0) throw("Out of memory.");

  shapefun(ShapeTypeN,nXi,xi,N);
  shapefun(ShapeTypeM,nXi,xi,M);
  shapederiv(ShapeTypeN,nXi,xi,dN);
  shapenatcoord(dN,nEltNod,nXi,EltNod,nat,EltDim);
  jacobian(nat,nXi,Jac,EltDim);

  for (unsigned int iXi=0; iXi<nXi; iXi++)
  {
     // Account for Radius in axisymmetric analysis
     double xiRadius=0.0;
     if (probAxi) for (unsigned int iEltNod=0; iEltNod<nEltNod; iEltNod++) xiRadius+=6.28318530717959*N[nEltNod*iXi+iEltNod]*EltNod[0*nEltNod+iEltNod];
     else xiRadius=1.0;

     for (unsigned int iuMode=0; iuMode<nuMode; iuMode++)
     {
      // first evaluate displacements for a single displacement mode
      complex<double> uvalx(0.0,0.0);
      complex<double> uvaly(0.0,0.0);
      complex<double> uvalz(0.0,0.0);
      for (unsigned int iEltColl=0; iEltColl<nEltColl; iEltColl++)
      {
        unsigned int iDof=nColDof*eltCollIndex[iEltColl];
        if (nColDof>0) uvalx += M[nEltColl*iXi+iEltColl]*u[iuMode*nDof +iDof+0];
        if (nColDof>1) uvaly += M[nEltColl*iXi+iEltColl]*u[iuMode*nDof +iDof+1];
        if (nColDof>2) uvalz += M[nEltColl*iXi+iEltColl]*u[iuMode*nDof +iDof+2];
      }

      // then multiply the displacement mode with all traction modes for all sets.
      for (unsigned int itMode=0; itMode<ntMode; itMode++)
      {
        for (unsigned int itSet=0; itSet<ntSet; itSet++)
        {
          complex<double> tvalx(0.0,0.0);
          complex<double> tvaly(0.0,0.0);
          complex<double> tvalz(0.0,0.0);
          for (unsigned int iEltColl=0; iEltColl<nEltColl; iEltColl++)
          {
            unsigned int iDof=nColDof*eltCollIndex[iEltColl];
            if (nColDof>0) tvalx += M[nEltColl*iXi+iEltColl]*t[itSet*ntMode*nDof+itMode*nDof+iDof+0];
            if (nColDof>1) tvaly += M[nEltColl*iXi+iEltColl]*t[itSet*ntMode*nDof+itMode*nDof+iDof+1];
            if (nColDof>2) tvalz += M[nEltColl*iXi+iEltColl]*t[itSet*ntMode*nDof+itMode*nDof+iDof+2];
          }
          OutMat[itSet*ntMode*nuMode+itMode*nuMode+iuMode]+=xiRadius*H[iXi]*Jac[iXi]*(tvalx*uvalx+tvaly*uvaly+tvalz*uvalz);
        }
      }
    }
  }
  delete [] EltNod;
  delete [] xi;
  delete [] H;
  delete [] N;
  delete [] M;
  delete [] dN;
  delete [] nat;
  delete [] Jac;
}
