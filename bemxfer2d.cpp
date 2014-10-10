#include "eltdef.h"
#include "bemcollpoints.h"
#include "gausspw.h"
#include "shapefun.h"
#include "bemnormal.h"
#include "greeneval2d.h"
#include "greenrotate2d.h"
#include <new>
#include <math.h>

using namespace std;

inline double sign(const double& a)
{
  if (a==0.0) return 1.0;
  else return (a>0.0 ? 1.0 : -1.0);
}
inline double sqr(const double& a)
{
  return a*a;
}

//======================================================================
// 2D BOUNDARY ELEMENT INTEGRATION
//======================================================================
void bemxfer2d(const double* const Nod,const unsigned int& nNod,
               const double* const Elt,const unsigned int& iElt,
               const unsigned int& nElt, const unsigned int* const EltCollIndex,
               const double* const Rec, const unsigned int nRec,
               bool* const boundaryRec,
               double* const URe, double* const UIm,
               double* const TRe, double* const TIm,
               const bool TmatOut,
               const unsigned int& nDof, const unsigned int& nRecDof,
               const unsigned int* const  TypeID,const unsigned int* const  nKeyOpt,
               const char* const* TypeName, const char* const* TypeKeyOpts,
               const unsigned int& nEltType,
               const void* const* const greenPtr, const unsigned int& nGrSet,
               const unsigned int& nugComp, const bool& ugCmplx, const bool& tgCmplx)
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
  unsigned int nGauss;
  unsigned int nEltDiv;
  unsigned int nGaussSing;
  unsigned int nEltDivSing;
  eltdef(EltType,TypeID,TypeName,TypeKeyOpts,nKeyOpt,nEltType,Parent,nEltNod,
         nEltColl,ShapeTypeN,ShapeTypeM,EltDim,AxiSym,Periodic,nGauss,nEltDiv,nGaussSing,nEltDivSing);
  const unsigned int nXi=nGauss*nEltDiv;

  unsigned int nColDof;                // Number of DOFs per collocation point
  if (nugComp==1) nColDof=1;  // 2D, out-of-plane (y)
  if (nugComp==4) nColDof=2;  // 2D, in-plane (x,z)
  if (nugComp==9) nColDof=3;  // 2.5D (x,y,z)

  unsigned int ntgComp;                // Number of components in green's stresses
  if (nugComp==1) ntgComp=2;  // 2D, out-of-plane (y)
  if (nugComp==4) ntgComp=6;  // 2D, in-plane (x,z)
  if (nugComp==9) ntgComp=18;  // 2.5D (x,y,z)

  int NodIndex;
  unsigned int NodID;
  double* const EltNod =new(nothrow) double[3*nEltNod];
  if (EltNod==0) throw("Out of memory.");

  // Determine coordinates of element nodes (element iElt).
  for (unsigned int iEltNod=0; iEltNod<nEltNod; iEltNod++)
  {
    NodID=(unsigned int)(Elt[(2+iEltNod)*nElt+iElt]);
    BemNodeIndex(Nod,nNod,NodID,NodIndex);
    EltNod[0*nEltNod+iEltNod]=Nod[1*nNod+NodIndex];
    EltNod[1*nEltNod+iEltNod]=Nod[2*nNod+NodIndex];
    EltNod[2*nEltNod+iEltNod]=Nod[3*nNod+NodIndex];
  }

  // Determine sample points for the element type (NumGauss,NumEltDiv).
  double* const xi=new(nothrow) double[nXi];
  if (xi==0) throw("Out of memory.");
  double* const H=new(nothrow) double[nXi];
  if (H==0) throw("Out of memory.");
  gausspw1D(nEltDiv,nGauss,xi,H);

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
  double* const xiCart=new(nothrow) double[2*nXi];
  if (xiCart==0) throw("Out of memory.");
  double* const normal=new(nothrow) double[3*nXi];
  if (normal==0) throw("Out of memory.");

  shapefun(ShapeTypeN,nXi,xi,N);
  shapefun(ShapeTypeM,nXi,xi,M);
  shapederiv(ShapeTypeN,nXi,xi,dN);
  shapenatcoord(dN,nEltNod,nXi,EltNod,nat,EltDim);
  jacobian(nat,nXi,Jac,EltDim);
  bemnormal(nat,nXi,EltDim,normal);

  for (unsigned int icomp=0; icomp<2*nXi; icomp++) xiCart[icomp]=0.0;
  for (unsigned int iXi=0; iXi<nXi; iXi++)
  {
    for (unsigned int iEltNod=0; iEltNod<nEltNod; iEltNod++)
    {
      xiCart[2*iXi+0]+=N[nEltNod*iXi+iEltNod]*EltNod[0*nEltNod+iEltNod];
      xiCart[2*iXi+1]+=N[nEltNod*iXi+iEltNod]*EltNod[2*nEltNod+iEltNod];
    }
  }

  double* const UgrRe=new(nothrow) double[nugComp*nGrSet];
  if (UgrRe==0) throw("Out of memory.");
  double* const UgrIm=new(nothrow) double[nugComp*nGrSet];
  if (UgrIm==0) throw("Out of memory.");
  double* const TgrRe=new(nothrow) double[ntgComp*nGrSet];
  if (TgrRe==0) throw("Out of memory.");
  double* const TgrIm=new(nothrow) double[ntgComp*nGrSet];
  if (TgrIm==0) throw("Out of memory.");
  double* const Tgr0Re=0;
  double* const Tgr0Im=0;

  double* const TXiRe=new(nothrow) double[nugComp*nGrSet];
  if (TXiRe==0) throw("Out of memory.");
  double* const TXiIm=new(nothrow) double[nugComp*nGrSet];
  if (TXiIm==0) throw("Out of memory.");
  double* const TXi0Re=0;
  double* const TXi0Im=0;

  // INITIALIZE INTERPOLATION OF GREEN'S FUNCTION
  unsigned int r1=0;
  unsigned int r2=1;
  bool extrapFlag=false;
  double* const interpr=new(nothrow) double[2];
  if (interpr==0) throw("Out of memory.");
  unsigned int z1=0;
  unsigned int z2=1;
  double* const interpz=new(nothrow) double[2];
  if (interpz==0) throw("Out of memory.");
  unsigned int zs1=0;

  for (unsigned int iComp=0; iComp<nugComp*nGrSet;iComp++)
  {
    TXiRe[iComp]=0.0;
    TXiIm[iComp]=0.0;
  }


  for (unsigned int iRec=0; iRec<nRec; iRec++)
  {
    if (!(boundaryRec[iRec]))  // If receiver not on interface
    {
      for (unsigned int iXi=0; iXi<nXi; iXi++)
      {
        const double Xdiff=xiCart[2*iXi+0]-Rec[0*nRec+iRec];
        const double Zdiff=xiCart[2*iXi+1]-Rec[2*nRec+iRec];
        const double xiR=fabs(Xdiff);
        const double xiZ=Zdiff;
        const double Xsgn=(Xdiff==0 ? 0.0 : sign(Xdiff));
      
        // EVALUATE GREEN'S FUNCTION
        const bool tg0Cmplx=false;
        greeneval2d(greenPtr,nGrSet,nugComp,ntgComp,ugCmplx,tgCmplx,tg0Cmplx,
                    xiR,xiZ,Xsgn,r1,r2,z1,z2,zs1,interpr,interpz,extrapFlag,
                    TmatOut,Rec,nRec,iRec,2,UgrRe,UgrIm,TgrRe,TgrIm,Tgr0Re,
                    Tgr0Im);
        greenrotate2d(normal,iXi,nGrSet,ntgComp,tgCmplx,tg0Cmplx,TgrRe,TgrIm,
                      Tgr0Re,Tgr0Im,TXiRe,TXiIm,TXi0Re,TXi0Im,TmatOut);
      
        // SUM UP RESULTS, FOR ALL ELEMENT COLLOCATION POINTS
        for (unsigned int iEltColl=0; iEltColl<nEltColl; iEltColl++)
        {
          double sumutil=H[iXi]*M[nEltColl*iXi+iEltColl]*Jac[iXi];
          unsigned int rowBeg=nColDof*iRec;
          unsigned int colBeg=nColDof*EltCollIndex[iEltColl];
      
          for (unsigned int iGrSet=0; iGrSet<nGrSet; iGrSet++)
          {
            unsigned int ind0 =nRecDof*nDof*iGrSet;
            if (nugComp==1){
              URe[ind0+nRecDof*(colBeg+0)+rowBeg+0]+=sumutil*UgrRe[nugComp*iGrSet+0];
              if (ugCmplx) UIm[ind0+nRecDof*(colBeg+0)+rowBeg+0]+=sumutil*UgrIm[nugComp*iGrSet+0];
              if (TmatOut)
              {
                TRe[ind0+nRecDof*(colBeg+0)+rowBeg+0]+=sumutil*TXiRe[nugComp*iGrSet+0];
                if (tgCmplx) TIm[ind0+nRecDof*(colBeg+0)+rowBeg+0]+=sumutil*TXiIm[nugComp*iGrSet+0];
              }
            }
            else if (nugComp==4){
              URe[ind0+nRecDof*(colBeg+0)+rowBeg+0]+=sumutil*UgrRe[nugComp*iGrSet+0];
              URe[ind0+nRecDof*(colBeg+1)+rowBeg+0]+=sumutil*UgrRe[nugComp*iGrSet+1];
              URe[ind0+nRecDof*(colBeg+0)+rowBeg+1]+=sumutil*UgrRe[nugComp*iGrSet+2];
              URe[ind0+nRecDof*(colBeg+1)+rowBeg+1]+=sumutil*UgrRe[nugComp*iGrSet+3];
              if (ugCmplx)
              {
                UIm[ind0+nRecDof*(colBeg+0)+rowBeg+0]+=sumutil*UgrIm[nugComp*iGrSet+0];
                UIm[ind0+nRecDof*(colBeg+1)+rowBeg+0]+=sumutil*UgrIm[nugComp*iGrSet+1];
                UIm[ind0+nRecDof*(colBeg+0)+rowBeg+1]+=sumutil*UgrIm[nugComp*iGrSet+2];
                UIm[ind0+nRecDof*(colBeg+1)+rowBeg+1]+=sumutil*UgrIm[nugComp*iGrSet+3];
              }
              if (TmatOut)
              {
                TRe[ind0+nRecDof*(colBeg+0)+rowBeg+0]+=sumutil*TXiRe[nugComp*iGrSet+0];     // txx
                TRe[ind0+nRecDof*(colBeg+1)+rowBeg+0]+=sumutil*TXiRe[nugComp*iGrSet+1];     // txz
                TRe[ind0+nRecDof*(colBeg+0)+rowBeg+1]+=sumutil*TXiRe[nugComp*iGrSet+2];     // tzx
                TRe[ind0+nRecDof*(colBeg+1)+rowBeg+1]+=sumutil*TXiRe[nugComp*iGrSet+3];     // tzz
                if (tgCmplx)
                {
                  TIm[ind0+nRecDof*(colBeg+0)+rowBeg+0]+=sumutil*TXiIm[nugComp*iGrSet+0];   // txx
                  TIm[ind0+nRecDof*(colBeg+1)+rowBeg+0]+=sumutil*TXiIm[nugComp*iGrSet+1];   // txz
                  TIm[ind0+nRecDof*(colBeg+0)+rowBeg+1]+=sumutil*TXiIm[nugComp*iGrSet+2];   // tzx
                  TIm[ind0+nRecDof*(colBeg+1)+rowBeg+1]+=sumutil*TXiIm[nugComp*iGrSet+3];   // tzz
                }
              }
            }
            else if (nugComp==9)
            {
              URe[ind0+nRecDof*(colBeg+0)+rowBeg+0]+=sumutil*UgrRe[nugComp*iGrSet+0];
              URe[ind0+nRecDof*(colBeg+1)+rowBeg+0]+=sumutil*UgrRe[nugComp*iGrSet+1];
              URe[ind0+nRecDof*(colBeg+2)+rowBeg+0]+=sumutil*UgrRe[nugComp*iGrSet+2];
              URe[ind0+nRecDof*(colBeg+0)+rowBeg+1]+=sumutil*UgrRe[nugComp*iGrSet+3];
              URe[ind0+nRecDof*(colBeg+1)+rowBeg+1]+=sumutil*UgrRe[nugComp*iGrSet+4];
              URe[ind0+nRecDof*(colBeg+2)+rowBeg+1]+=sumutil*UgrRe[nugComp*iGrSet+5];
              URe[ind0+nRecDof*(colBeg+0)+rowBeg+2]+=sumutil*UgrRe[nugComp*iGrSet+6];
              URe[ind0+nRecDof*(colBeg+1)+rowBeg+2]+=sumutil*UgrRe[nugComp*iGrSet+7];
              URe[ind0+nRecDof*(colBeg+2)+rowBeg+2]+=sumutil*UgrRe[nugComp*iGrSet+8];
              if (ugCmplx)
              {
                UIm[ind0+nRecDof*(colBeg+0)+rowBeg+0]+=sumutil*UgrIm[nugComp*iGrSet+0];
                UIm[ind0+nRecDof*(colBeg+1)+rowBeg+0]+=sumutil*UgrIm[nugComp*iGrSet+1];
                UIm[ind0+nRecDof*(colBeg+2)+rowBeg+0]+=sumutil*UgrIm[nugComp*iGrSet+2];
                UIm[ind0+nRecDof*(colBeg+0)+rowBeg+1]+=sumutil*UgrIm[nugComp*iGrSet+3];
                UIm[ind0+nRecDof*(colBeg+1)+rowBeg+1]+=sumutil*UgrIm[nugComp*iGrSet+4];
                UIm[ind0+nRecDof*(colBeg+2)+rowBeg+1]+=sumutil*UgrIm[nugComp*iGrSet+5];
                UIm[ind0+nRecDof*(colBeg+0)+rowBeg+2]+=sumutil*UgrIm[nugComp*iGrSet+6];
                UIm[ind0+nRecDof*(colBeg+1)+rowBeg+2]+=sumutil*UgrIm[nugComp*iGrSet+7];
                UIm[ind0+nRecDof*(colBeg+2)+rowBeg+2]+=sumutil*UgrIm[nugComp*iGrSet+8];
              }
              if (TmatOut)
              {
                TRe[ind0+nRecDof*(colBeg+0)+rowBeg+0]+=sumutil*TXiRe[nugComp*iGrSet+0];  // txx
                TRe[ind0+nRecDof*(colBeg+1)+rowBeg+0]+=sumutil*TXiRe[nugComp*iGrSet+1];  // txy
                TRe[ind0+nRecDof*(colBeg+2)+rowBeg+0]+=sumutil*TXiRe[nugComp*iGrSet+2];  // txz
                TRe[ind0+nRecDof*(colBeg+0)+rowBeg+1]+=sumutil*TXiRe[nugComp*iGrSet+3];  // tyx
                TRe[ind0+nRecDof*(colBeg+1)+rowBeg+1]+=sumutil*TXiRe[nugComp*iGrSet+4];  // tyy
                TRe[ind0+nRecDof*(colBeg+2)+rowBeg+1]+=sumutil*TXiRe[nugComp*iGrSet+5];  // tyz
                TRe[ind0+nRecDof*(colBeg+0)+rowBeg+2]+=sumutil*TXiRe[nugComp*iGrSet+6];  // tzx
                TRe[ind0+nRecDof*(colBeg+1)+rowBeg+2]+=sumutil*TXiRe[nugComp*iGrSet+7];  // tzy
                TRe[ind0+nRecDof*(colBeg+2)+rowBeg+2]+=sumutil*TXiRe[nugComp*iGrSet+8];  // tzz
                if (tgCmplx)
                {
                  TIm[ind0+nRecDof*(colBeg+0)+rowBeg+0]+=sumutil*TXiIm[nugComp*iGrSet+0];  // txx
                  TIm[ind0+nRecDof*(colBeg+1)+rowBeg+0]+=sumutil*TXiIm[nugComp*iGrSet+1];  // txy
                  TIm[ind0+nRecDof*(colBeg+2)+rowBeg+0]+=sumutil*TXiIm[nugComp*iGrSet+2];  // txz
                  TIm[ind0+nRecDof*(colBeg+0)+rowBeg+1]+=sumutil*TXiIm[nugComp*iGrSet+3];  // tyx
                  TIm[ind0+nRecDof*(colBeg+1)+rowBeg+1]+=sumutil*TXiIm[nugComp*iGrSet+4];  // tyy
                  TIm[ind0+nRecDof*(colBeg+2)+rowBeg+1]+=sumutil*TXiIm[nugComp*iGrSet+5];  // tyz
                  TIm[ind0+nRecDof*(colBeg+0)+rowBeg+2]+=sumutil*TXiIm[nugComp*iGrSet+6];  // tzx
                  TIm[ind0+nRecDof*(colBeg+1)+rowBeg+2]+=sumutil*TXiIm[nugComp*iGrSet+7];  // tzy
                  TIm[ind0+nRecDof*(colBeg+2)+rowBeg+2]+=sumutil*TXiIm[nugComp*iGrSet+8];  // tzz
                }
              }
            }
          }
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
  delete [] xiCart;
  delete [] normal;
  delete [] interpr;
  delete [] interpz;
  delete [] UgrRe;
  delete [] UgrIm;
  delete [] TgrRe;
  delete [] TgrIm;
  delete [] TXiRe;
  delete [] TXiIm;
}
