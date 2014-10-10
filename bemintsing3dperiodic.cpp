#include "eltdef.h"
#include "shapefun.h"
#include "bemnormal.h"
#include "gausspw.h"
#include "bemcollpoints.h"
#include "greeneval3d.h"
#include "greenrotate3d.h"
#include "mex.h"
#include <new>
#include <math.h>
#include <complex> 

using namespace std;

/******************************************************************************/
const complex<double> i(0.0,1.0);
/******************************************************************************/
inline double sign(const double& a)
{
  if (a==0.0) return 1.0;
  else return (a>0.0 ? 1.0 : -1.0);
}
/******************************************************************************/
inline double sqr(const double& a)
{
  return a*a;
}

//============================================================================//
//  THREE-DIMENSIONAL SINGULAR INTEGRATION
//============================================================================//
void bemintsing3dperiodic(const double* const Nod, const unsigned int& nNod,
                          const double* const Elt, const unsigned int& iElt,
                          const unsigned int& nElt,
                          const unsigned int* const  TypeID, const unsigned int* const  nKeyOpt,
                          const char* const TypeName[], const char* const TypeKeyOpts[], 
                          const unsigned int& nEltType,
                          const double* const Coll,const unsigned int& nColl, const unsigned int& iColl, 
                          const unsigned int* const EltCollIndex, const unsigned int& nDof,
                          const double* const xiSing, const void* const* const greenPtr, 
                          const unsigned int& nGrSet, const bool& ugCmplx, 
                          const bool& tgCmplx, const bool& tg0Cmplx, double* const URe, 
                          double* const UIm, double* const TRe, double* const TIm, 
                          const bool& UmatOut,const bool& TmatOut, const unsigned int nWave)
{
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
         nEltColl,ShapeTypeN,ShapeTypeM,EltDim,AxiSym,Periodic,nGauss,nEltDiv,
         nGaussSing,nEltDivSing);
  const unsigned int nXi=nEltDivSing*nEltDivSing*nGaussSing*nGaussSing;

  int NodIndex;
  unsigned int NodID;
  double* const NodCoord =new(nothrow) double[3*nEltNod];
  if (NodCoord==0) throw("Out of memory.");

  for (unsigned int iEltNod=0; iEltNod<nEltNod; iEltNod++)
  {
    NodID=(unsigned int)(Elt[(2+iEltNod)*nElt+iElt]);
    BemNodeIndex(Nod,nNod,NodID,NodIndex);
    NodCoord[0*nEltNod+iEltNod]=Nod[1*nNod+NodIndex];
    NodCoord[1*nEltNod+iEltNod]=Nod[2*nNod+NodIndex];
    NodCoord[2*nEltNod+iEltNod]=Nod[3*nNod+NodIndex];
  }   
      
  // ELEMENT TRIANGLE DIVISION
  unsigned int nDiv;
  double* const am=new(nothrow) double[8];
  if (am==0) throw("Out of memory.");
  double* const a1=new(nothrow) double[8];
  if (a1==0) throw("Out of memory.");
  double* const a2=new(nothrow) double[8];
  if (a2==0) throw("Out of memory.");
  double* const rhom=new(nothrow) double[8];
  if (rhom==0) throw("Out of memory.");
  double* const rho1=new(nothrow) double[8];
  if (rho1==0) throw("Out of memory.");
  double* const rho2=new(nothrow) double[8];
  if (rho2==0) throw("Out of memory.");
  triangdiv(xiSing,Parent,nDiv,am,a1,a2,rhom,rho1,rho2);

  // H and xi in v1,v2
  double* const v=new(nothrow) double[2*nXi];
  if (v==0) throw("Out of memory.");
  double* const H=new(nothrow) double[nXi];
  if (H==0) throw("Out of memory.");
  gausspw2D(nEltDivSing,nGaussSing,v,H);

  double* const a=new(nothrow) double[nXi];
  if (a==0) throw("Out of memory.");
  double* const rho=new(nothrow) double[nXi];
  if (rho==0) throw("Out of memory.");

  double* const xi=new(nothrow) double[2*nXi];
  if (xi==0) throw("Out of memory.");
  double* const N=new(nothrow) double[nXi*nEltNod];
  if (N==0) throw("Out of memory.");
  double* const M=new(nothrow) double[nXi*nEltColl];
  if (M==0) throw("Out of memory.");
  double* const Mmod=new(nothrow) double[nXi*nEltColl];
  if (Mmod==0) throw("Out of memory.");
  double* const dN=new(nothrow) double[2*nXi*nEltNod];
  if (dN==0) throw("Out of memory.");
  double* const nat=new(nothrow) double[6*nXi];
  if (nat==0) throw("Out of memory.");
  double* const Jac=new(nothrow) double[nXi];
  if (Jac==0) throw("Out of memory.");
  double* const normal=new(nothrow) double[3*nXi];
  if (normal==0) throw("Out of memory.");
  double* const xiCart=new(nothrow) double[3*nXi];
  if (xiCart==0) throw("Out of memory.");
  double* const UgrRe=new(nothrow) double[5*nGrSet];
  if (UgrRe==0) throw("Out of memory.");
  double* const UgrIm=new(nothrow) double[5*nGrSet];
  if (UgrIm==0) throw("Out of memory.");
  double* const TgrRe=new(nothrow) double[10*nGrSet];
  if (TgrRe==0) throw("Out of memory.");
  double* const TgrIm=new(nothrow) double[10*nGrSet];
  if (TgrIm==0) throw("Out of memory.");
  double* const Tgr0Re=new(nothrow) double[10*nGrSet];
  if (Tgr0Re==0) throw("Out of memory.");
  double* const Tgr0Im=new(nothrow) double[10*nGrSet];
  if (Tgr0Im==0) throw("Out of memory.");
  double* const TXi0Re=new(nothrow) double[9*nGrSet];
  if (TXi0Re==0) throw("Out of memory.");
  double* const TXi0Im=new(nothrow) double[9*nGrSet];
  if (TXi0Im==0) throw("Out of memory.");
  double* const UXiRe=new(nothrow) double[9*nGrSet];
  if (UXiRe==0) throw("Out of memory.");
  double* const UXiIm=new(nothrow) double[9*nGrSet];
  if (UXiIm==0) throw("Out of memory.");
  double* const TXiRe=new(nothrow) double[9*nGrSet];
  if (TXiRe==0) throw("Out of memory.");
  double* const TXiIm=new(nothrow) double[9*nGrSet];
  if (TXiIm==0) throw("Out of memory.");
  
  for (unsigned int iComp=0; iComp<9*nGrSet;iComp++)
  {
    TXiRe[iComp]=0.0;
    TXiIm[iComp]=0.0;
    TXi0Re[iComp]=0.0;
    TXi0Im[iComp]=0.0;
  }

  // Initialize interpolation of Green's function
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

  for (unsigned int iDiv=0; iDiv<nDiv; iDiv++) if (rhom[iDiv]>1e-10)
  {
    for (unsigned int iXi=0; iXi<nXi; iXi++)
    {
      a[iXi]  =0.5*(a2[iDiv]-a1[iDiv])*v[nXi+iXi]+0.5*(a2[iDiv]+a1[iDiv]);
      rho[iXi]=0.5*rhom[iDiv]/cos(a[iXi]-am[iDiv])*(1+v[iXi]);
      xi[iXi]=xiSing[0]+rho[iXi]*cos(a[iXi]);
      xi[nXi+iXi]=xiSing[1]+rho[iXi]*sin(a[iXi]);
    }
    
    shapefun(ShapeTypeN,nXi,xi,N);
    shapefun(ShapeTypeM,nXi,xi,M);
    shapederiv(ShapeTypeN,nXi,xi,dN);
    shapenatcoord(dN,nEltNod,nXi,NodCoord,nat,EltDim);
    jacobian(nat,nXi,Jac,EltDim);
    if (TmatOut) bemnormal(nat,nXi,EltDim,normal);

    for (unsigned int icomp=0; icomp<3*nXi; icomp++) xiCart[icomp]=0.0;
    for (unsigned int iXi=0; iXi<nXi; iXi++)
    {
      for (unsigned int iNod=0; iNod<nEltNod; iNod++)
      {
        xiCart[3*iXi+0]+=N[nEltNod*iXi+iNod]*NodCoord[0*nEltNod+iNod];
        xiCart[3*iXi+1]+=N[nEltNod*iXi+iNod]*NodCoord[1*nEltNod+iNod];
        xiCart[3*iXi+2]+=N[nEltNod*iXi+iNod]*NodCoord[2*nEltNod+iNod];
      }
      const double Xdiff=xiCart[3*iXi+0]-Coll[2*nColl+iColl];
      const double Ydiff=xiCart[3*iXi+1]-Coll[3*nColl+iColl];
      const double Zdiff=xiCart[3*iXi+2]-Coll[4*nColl+iColl];
      
      const double xiR=sqrt(Xdiff*Xdiff+Ydiff*Ydiff);
      const double xiTheta=atan2(Ydiff,Xdiff);
      const double xiZ=Zdiff;
    
      if ((xiR==0)&(xiZ==0)) throw("An integration point coincides with the collocation point for singular integration.");
    
      // EVALUATE GREEN'S FUNCTION
      greeneval3d(greenPtr,nGrSet,ugCmplx,tgCmplx,tg0Cmplx,xiR,xiZ,r1,r2,z1,z2,zs1,
                  interpr,interpz,extrapFlag,UmatOut,TmatOut,Coll,nColl,iColl,4,UgrRe,
                  UgrIm,TgrRe,TgrIm,Tgr0Re,Tgr0Im);
    
      // ROTATE GREEN'S FUNCTIONS
      greenrotate3d(normal,iXi,xiTheta,nGrSet,ugCmplx,
                    tgCmplx,tg0Cmplx,UgrRe,UgrIm,TgrRe,TgrIm,
                    Tgr0Re,Tgr0Im,UXiRe,UXiIm,TXiRe,TXiIm,TXi0Re,
                    TXi0Im,UmatOut,TmatOut);
      
      for (unsigned int iWave=0; iWave<nWave; iWave++)                                        
      {                                                                  
        for (unsigned int iEltColl=0; iEltColl<nEltColl; iEltColl++)
        {
          double sumutil=rho[iXi]*H[iXi]*M[nEltColl*iXi+iEltColl]*Jac[iXi]*0.25*(a2[iDiv]-a1[iDiv])*rhom[iDiv]/cos(a[iXi]-am[iDiv]);
          unsigned int rowBeg=3*iColl;
          unsigned int colBeg=3*EltCollIndex[iEltColl];
          for (unsigned int iGrSet=0; iGrSet<nGrSet; iGrSet++)
          {
            const unsigned int ind0 =nDof*nDof*(nGrSet*iWave+iGrSet);
            URe[ind0+nDof*(colBeg+0)+rowBeg+0]+=sumutil*UXiRe[9*iGrSet+0];  // ugxx 
            URe[ind0+nDof*(colBeg+1)+rowBeg+0]+=sumutil*UXiRe[9*iGrSet+1];  // ugxy 
            URe[ind0+nDof*(colBeg+2)+rowBeg+0]+=sumutil*UXiRe[9*iGrSet+2];  // ugxz 
            URe[ind0+nDof*(colBeg+0)+rowBeg+1]+=sumutil*UXiRe[9*iGrSet+3];  // ugyx 
            URe[ind0+nDof*(colBeg+1)+rowBeg+1]+=sumutil*UXiRe[9*iGrSet+4];  // ugyy 
            URe[ind0+nDof*(colBeg+2)+rowBeg+1]+=sumutil*UXiRe[9*iGrSet+5];  // ugyz 
            URe[ind0+nDof*(colBeg+0)+rowBeg+2]+=sumutil*UXiRe[9*iGrSet+6];  // ugzx 
            URe[ind0+nDof*(colBeg+1)+rowBeg+2]+=sumutil*UXiRe[9*iGrSet+7];  // ugzy 
            URe[ind0+nDof*(colBeg+2)+rowBeg+2]+=sumutil*UXiRe[9*iGrSet+8];  // ugzz 
            if (ugCmplx)
            {
              UIm[ind0+nDof*(colBeg+0)+rowBeg+0]+=sumutil*UXiIm[9*iGrSet+0];
              UIm[ind0+nDof*(colBeg+1)+rowBeg+0]+=sumutil*UXiIm[9*iGrSet+1];
              UIm[ind0+nDof*(colBeg+2)+rowBeg+0]+=sumutil*UXiIm[9*iGrSet+2];
              UIm[ind0+nDof*(colBeg+0)+rowBeg+1]+=sumutil*UXiIm[9*iGrSet+3];
              UIm[ind0+nDof*(colBeg+1)+rowBeg+1]+=sumutil*UXiIm[9*iGrSet+4];
              UIm[ind0+nDof*(colBeg+2)+rowBeg+1]+=sumutil*UXiIm[9*iGrSet+5];
              UIm[ind0+nDof*(colBeg+0)+rowBeg+2]+=sumutil*UXiIm[9*iGrSet+6];
              UIm[ind0+nDof*(colBeg+1)+rowBeg+2]+=sumutil*UXiIm[9*iGrSet+7];
              UIm[ind0+nDof*(colBeg+2)+rowBeg+2]+=sumutil*UXiIm[9*iGrSet+8];
            }
            if (TmatOut)
            {
              TRe[ind0+nDof*(colBeg+0)+rowBeg+0]+=sumutil*TXiRe[9*iGrSet+0];  // txx
              TRe[ind0+nDof*(colBeg+1)+rowBeg+0]+=sumutil*TXiRe[9*iGrSet+1];  // txy
              TRe[ind0+nDof*(colBeg+2)+rowBeg+0]+=sumutil*TXiRe[9*iGrSet+2];  // txz
              TRe[ind0+nDof*(colBeg+0)+rowBeg+1]+=sumutil*TXiRe[9*iGrSet+3];  // tyx
              TRe[ind0+nDof*(colBeg+1)+rowBeg+1]+=sumutil*TXiRe[9*iGrSet+4];  // tyy
              TRe[ind0+nDof*(colBeg+2)+rowBeg+1]+=sumutil*TXiRe[9*iGrSet+5];  // tyz
              TRe[ind0+nDof*(colBeg+0)+rowBeg+2]+=sumutil*TXiRe[9*iGrSet+6];  // tzx
              TRe[ind0+nDof*(colBeg+1)+rowBeg+2]+=sumutil*TXiRe[9*iGrSet+7];  // tzy
              TRe[ind0+nDof*(colBeg+2)+rowBeg+2]+=sumutil*TXiRe[9*iGrSet+8];  // tzz
              if (tgCmplx)
              {
                TIm[ind0+nDof*(colBeg+0)+rowBeg+0]+=sumutil*TXiIm[9*iGrSet+0];
                TIm[ind0+nDof*(colBeg+1)+rowBeg+0]+=sumutil*TXiIm[9*iGrSet+1];
                TIm[ind0+nDof*(colBeg+2)+rowBeg+0]+=sumutil*TXiIm[9*iGrSet+2];
                TIm[ind0+nDof*(colBeg+0)+rowBeg+1]+=sumutil*TXiIm[9*iGrSet+3];
                TIm[ind0+nDof*(colBeg+1)+rowBeg+1]+=sumutil*TXiIm[9*iGrSet+4];
                TIm[ind0+nDof*(colBeg+2)+rowBeg+1]+=sumutil*TXiIm[9*iGrSet+5];
                TIm[ind0+nDof*(colBeg+0)+rowBeg+2]+=sumutil*TXiIm[9*iGrSet+6];
                TIm[ind0+nDof*(colBeg+1)+rowBeg+2]+=sumutil*TXiIm[9*iGrSet+7];
                TIm[ind0+nDof*(colBeg+2)+rowBeg+2]+=sumutil*TXiIm[9*iGrSet+8];
              }
              // Account for singular part of Green's function on the diagonal
              TRe[ind0+nDof*(rowBeg+0)+rowBeg+0]-=sumutil*TXi0Re[9*iGrSet+0];  // txx
              TRe[ind0+nDof*(rowBeg+1)+rowBeg+0]-=sumutil*TXi0Re[9*iGrSet+1];  // txy
              TRe[ind0+nDof*(rowBeg+2)+rowBeg+0]-=sumutil*TXi0Re[9*iGrSet+2];  // txz
              TRe[ind0+nDof*(rowBeg+0)+rowBeg+1]-=sumutil*TXi0Re[9*iGrSet+3];  // tyx
              TRe[ind0+nDof*(rowBeg+1)+rowBeg+1]-=sumutil*TXi0Re[9*iGrSet+4];  // tyy
              TRe[ind0+nDof*(rowBeg+2)+rowBeg+1]-=sumutil*TXi0Re[9*iGrSet+5];  // tyz
              TRe[ind0+nDof*(rowBeg+0)+rowBeg+2]-=sumutil*TXi0Re[9*iGrSet+6];  // tzx
              TRe[ind0+nDof*(rowBeg+1)+rowBeg+2]-=sumutil*TXi0Re[9*iGrSet+7];  // tzy
              TRe[ind0+nDof*(rowBeg+2)+rowBeg+2]-=sumutil*TXi0Re[9*iGrSet+8];  // tzz
              if (tgCmplx)
              {
                TIm[ind0+nDof*(rowBeg+0)+rowBeg+0]-=sumutil*TXi0Im[9*iGrSet+0];
                TIm[ind0+nDof*(rowBeg+1)+rowBeg+0]-=sumutil*TXi0Im[9*iGrSet+1];
                TIm[ind0+nDof*(rowBeg+2)+rowBeg+0]-=sumutil*TXi0Im[9*iGrSet+2];
                TIm[ind0+nDof*(rowBeg+0)+rowBeg+1]-=sumutil*TXi0Im[9*iGrSet+3];
                TIm[ind0+nDof*(rowBeg+1)+rowBeg+1]-=sumutil*TXi0Im[9*iGrSet+4];
                TIm[ind0+nDof*(rowBeg+2)+rowBeg+1]-=sumutil*TXi0Im[9*iGrSet+5];
                TIm[ind0+nDof*(rowBeg+0)+rowBeg+2]-=sumutil*TXi0Im[9*iGrSet+6];
                TIm[ind0+nDof*(rowBeg+1)+rowBeg+2]-=sumutil*TXi0Im[9*iGrSet+7];
                TIm[ind0+nDof*(rowBeg+2)+rowBeg+2]-=sumutil*TXi0Im[9*iGrSet+8];
              }
            }
          }
        }
      }
    }
  }
  delete [] NodCoord;
  delete [] am;
  delete [] a1;
  delete [] a2;
  delete [] rhom;
  delete [] rho1;
  delete [] rho2;
  delete [] v;
  delete [] H;
  delete [] a;
  delete [] rho;
  delete [] xi;
  delete [] N;
  delete [] M;
  delete [] Mmod;
  delete [] dN;
  delete [] nat;
  delete [] Jac;
  delete [] normal;
  delete [] xiCart;
  delete [] interpr;
  delete [] interpz;
  delete [] UgrRe;
  delete [] UgrIm;
  delete [] TgrRe;
  delete [] TgrIm;
  delete [] UXiRe;
  delete [] UXiIm;
  delete [] TXiRe;
  delete [] TXiIm;
  delete [] Tgr0Re;
  delete [] Tgr0Im;
  delete [] TXi0Re;
  delete [] TXi0Im;
}


