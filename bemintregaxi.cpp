#include "eltdef.h"
#include "bemcollpoints.h"
#include "gausspw.h"
#include "shapefun.h"
#include "bemnormal.h"
#include "greeneval3d.h"
#include "greenrotate3d.h"
#include <new>
#include <math.h>

#include "mex.h"

using namespace std;

inline double sqr(const double& a)
{
  return a*a;
}
inline double sign(const double& a)
{
  if (a==0.0) return 1.0;
  else return (a>0.0 ? 1.0 : -1.0);
}

//======================================================================
// AXISYMMETRIC REGULAR BOUNDARY ELEMENT INTEGRATION
//======================================================================
void bemintregaxi(const double* const Nod, const unsigned int& nNod,
                  const double* const Elt, const unsigned int& iElt, const unsigned int& nElt,
                  const unsigned int* const  TypeID, const unsigned int* const nKeyOpt,
                  const char* const TypeName[], const char* const TypeKeyOpts[],
                  const unsigned int& nEltType, const double* const Coll,
                  const unsigned int& nColl, const unsigned int* const RegularColl,
                  const unsigned int* const EltCollIndex, const unsigned int& nDof,
                  const void* const* const greenPtr, const unsigned int& nGrSet,
                  const bool& ugCmplx, const bool& tgCmplx,
                  const bool& tg0Cmplx, double* const URe, double* const UIm,
                  double* const TRe, double* const TIm, const bool UmatOut, const bool TmatOut)
{
  // Element properties
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

  int NodIndex;
  unsigned int NodID;
  double* const EltNod =new(nothrow) double[3*nEltNod];
  if (EltNod==0) throw("Out of memory.");

  // Element node coordinates (of element iElt)
  for (unsigned int iEltNod=0; iEltNod<nEltNod; iEltNod++)
  {
    NodID=(unsigned int)(Elt[(2+iEltNod)*nElt+iElt]);
    BemNodeIndex(Nod,nNod,NodID,NodIndex);
    EltNod[0*nEltNod+iEltNod]=Nod[1*nNod+NodIndex];
    EltNod[1*nEltNod+iEltNod]=Nod[2*nNod+NodIndex];
    EltNod[2*nEltNod+iEltNod]=Nod[3*nNod+NodIndex];
  }

  // --- NUMBER OF INTEGRATION POINTS ---
  //   The number of integration points in Azimuthal direction is determined
  //   so that the distance between the integration points is similar to the
  //   distance in the radial direction, with a minimum of nGauss integration
  //   points.
  double xLen=EltNod[0*nEltNod+0]-EltNod[0*nEltNod+nEltNod-1];
  double zLen=EltNod[2*nEltNod+0]-EltNod[2*nEltNod+nEltNod-1];
  double eLen=sqrt(xLen*xLen+zLen*zLen);
  
  double RadiusAverage=0.50*(EltNod[0*nEltNod+0]+EltNod[0*nEltNod+nEltNod-1]); // Average radius
  unsigned int nAzimDiv=(unsigned int)(ceil(3.141592653589793*RadiusAverage/eLen));
  const unsigned int nXi1=nGauss*nEltDiv;         // Number of Radial points.
  unsigned int nXi2=nGauss;                       // Number of Azimuthal points.
  if ((eLen>0)&&(RadiusAverage>0)&&(nAzimDiv>nXi2)) nXi2=nAzimDiv;
  const unsigned int nXi=nXi1*nXi2;              // Total number of Integration points.
  
  // --- GENERATE SURFACE SAMPLING --
  //   For an axisymmetric analysis, the surface is constructed from
  //   a generating line mesh.
  //
  // Integration points and weights for the axisymmetric element
  // --- In Radial direction
  double* const xi1=new(nothrow) double[nXi1];
    if (xi1==0) throw("Out of memory.");
  double* const H1=new(nothrow) double[nXi1];
    if (H1==0) throw("Out of memory.");
  gausspw1D(nEltDiv,nGauss,xi1,H1);
  // --- In Azimuthal direction
  double* const xi2=new(nothrow) double[nXi2];
    if (xi2==0) throw("Out of memory.");
  double* const H2=new(nothrow) double[nXi2];
    if (H2==0) throw("Out of memory.");
  gausspw1D(nXi2,1,xi2,H2);

  // Weights for integration points on the surface
  double* const H=new(nothrow) double[nXi];
    if (H==0) throw("Out of memory.");
  for (unsigned int iXi1=0; iXi1<nXi1; iXi1++)  // radial direction
  {
    for (unsigned int iXi2=0; iXi2<nXi2; iXi2++) // azimuthal direction
    {
      const unsigned int iXi=nXi2*iXi1+iXi2;
      H[iXi]=H1[iXi1]*H2[iXi2];
    }
  }

  // Shape functions in the sample points (Radial sampling).
  double* const N1=new(nothrow) double[nXi1*nEltNod];
    if (N1==0) throw("Out of memory.");
  double* const M1=new(nothrow) double[nXi1*nEltColl];
    if (M1==0) throw("Out of memory.");
  double* const dN1=new(nothrow) double[2*nXi1*nEltNod];
    if (dN1==0) throw("Out of memory.");
  double* const nat1=new(nothrow) double[6*nXi1];
    if (nat1==0) throw("Out of memory.");
  double* const Jac1=new(nothrow) double[nXi1];
    if (Jac1==0) throw("Out of memory.");
  double* const normal1=new(nothrow) double[3*nXi1];
    if (normal1==0) throw("Out of memory.");

  shapefun(ShapeTypeN,nXi1,xi1,N1);
  shapefun(ShapeTypeM,nXi1,xi1,M1);
  shapederiv(ShapeTypeN,nXi1,xi1,dN1);
  shapenatcoord(dN1,nEltNod,nXi1,EltNod,nat1,EltDim);
  jacobian(nat1,nXi1,Jac1,EltDim);
  bemnormal(nat1,nXi1,EltDim,normal1);

  // Nodal coordinates of all integration points on the surface
  double* const xiCart=new(nothrow) double[3*nXi];
    if (xiCart==0) throw("Out of memory.");
  for (unsigned int icomp=0; icomp<3*nXi; icomp++) xiCart[icomp]=0.0;

  // Non-linear mapping for theta (cfr. Book Dominguez)
  //    theta = pi/4*(1 + xi)^2;
  double* const xi2theta=new(nothrow) double[nXi2];
    if (xi2theta==0) throw("Out of memory.");
  for (unsigned int iXi2=0;iXi2<nXi2;iXi2++) xi2theta[iXi2]=0.78539816339745*sqr(1.0+xi2[iXi2]);

  unsigned int* const xiiXi2=new(nothrow) unsigned int[nXi];
    if (xiiXi2==0) throw("Out of memory.");

  // Jacobian
  double* const Jac=new(nothrow) double[nXi];
    if (Jac==0) throw("Out of memory.");
  for (unsigned int iXi1=0; iXi1<nXi1; iXi1++)  // Radial direction
  {
    double xi1r=0.0;
    double xi1z=0.0;
    for (unsigned int iEltNod=0; iEltNod<nEltNod; iEltNod++)
    {
      // Coordinates of points on the generating line
      xi1r+=N1[nEltNod*iXi1+iEltNod]*EltNod[0*nEltNod+iEltNod];
      xi1z+=N1[nEltNod*iXi1+iEltNod]*EltNod[2*nEltNod+iEltNod];
    }
    for (unsigned int iXi2=0; iXi2<nXi2; iXi2++) // Azimuthal direction
    {
      const unsigned int iXi=nXi2*iXi1+iXi2;
      xiiXi2[iXi]=iXi2;
      xiCart[3*iXi+0]=xi1r*cos(xi2theta[iXi2]);
      xiCart[3*iXi+1]=xi1r*sin(xi2theta[iXi2]);
      xiCart[3*iXi+2]=xi1z;
      
      // Integration according to:
      // Jac = r * J1 * pi/2*(1+xi2)
      // where J1 is the Jacobian for the two-dimensional integration.
      Jac[iXi]= Jac1[iXi1]*xi1r*1.57079632679490*(1.0+xi2[iXi2]);
    }
  }

  // Normal for all integration points on the surface
  double* const normal=new(nothrow) double[3*nXi];
  if (normal==0) throw("Out of memory.");

  for (unsigned int iXi1=0; iXi1<nXi1; iXi1++)  // Radial direction
  {
    double xi1nr=normal1[3*iXi1+0];
    double xi1nz=normal1[3*iXi1+2];
    for (unsigned int iXi2=0; iXi2<nXi2; iXi2++) // Azimuthal direction
    {
      const unsigned int iXi=nXi2*iXi1+iXi2;
      normal[3*iXi+0]= xi1nr*cos(xi2theta[iXi2]);
      normal[3*iXi+1]= xi1nr*sin(xi2theta[iXi2]);
      normal[3*iXi+2]= xi1nz;
    }
  }

  // Shape function for all integration points on the surface.
  double* const M=new(nothrow) double[nXi*nEltColl];
    if (M==0) throw("Out of memory.");

  for (unsigned int iXi1=0; iXi1<nXi1; iXi1++)  // Radial direction
  {
    for (unsigned int iXi2=0; iXi2<nXi2; iXi2++) // Azimuthal direction
    {
      const unsigned int iXi=nXi2*iXi1+iXi2;
      for (unsigned int iEltColl=0; iEltColl<nEltColl; iEltColl++)
      {
        M[nEltColl*iXi+iEltColl]=M1[nEltColl*iXi1+iEltColl];
      }
    }
  }

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

  double* const UXiRe=new(nothrow) double[9*nGrSet];
  if (UXiRe==0) throw("Out of memory.");
  double* const UXiIm=new(nothrow) double[9*nGrSet];
  if (UXiIm==0) throw("Out of memory.");
  double* const TXiRe=new(nothrow) double[9*nGrSet];
  if (TXiRe==0) throw("Out of memory.");
  double* const TXiIm=new(nothrow) double[9*nGrSet];
  if (TXiIm==0) throw("Out of memory.");
  double* const TXi0Re=new(nothrow) double[9*nGrSet];
  if (TXi0Re==0) throw("Out of memory.");
  double* const TXi0Im=new(nothrow) double[9*nGrSet];
  if (TXi0Im==0) throw("Out of memory.");

  for (unsigned int iComp=0; iComp<9*nGrSet;iComp++)
  {
    TXiRe[iComp]=0.0;
    TXiIm[iComp]=0.0;
    TXi0Re[iComp]=0.0;
    TXi0Im[iComp]=0.0;
  }

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


  for (unsigned int iColl=0; iColl<nColl; iColl++)
  {
    if (RegularColl[iColl]==1)
    {
      for (unsigned int iXi=0; iXi<nXi; iXi++)
      {
        const double Xdiff=xiCart[3*iXi+0]-Coll[2*nColl+iColl];
        const double Ydiff=xiCart[3*iXi+1]-Coll[3*nColl+iColl];
        const double Zdiff=xiCart[3*iXi+2]-Coll[4*nColl+iColl];
      
        const double xiR=sqrt(Xdiff*Xdiff + Ydiff*Ydiff);
        const double xiTheta=atan2(Ydiff,Xdiff);
        const double xiZ=Zdiff;
      
        // EVALUATE GREEN'S FUNCTION
        greeneval3d(greenPtr,nGrSet,ugCmplx,tgCmplx,tg0Cmplx,xiR,xiZ,r1,r2,z1,z2,zs1,
                    interpr,interpz,extrapFlag,UmatOut,TmatOut,Coll,nColl,iColl,4,UgrRe,
                    UgrIm,TgrRe,TgrIm,Tgr0Re,Tgr0Im);
        greenrotate3d(normal,iXi,xiTheta,nGrSet,ugCmplx,
                      tgCmplx,tg0Cmplx,UgrRe,UgrIm,TgrRe,TgrIm,
                      Tgr0Re,Tgr0Im,UXiRe,UXiIm,TXiRe,TXiIm,TXi0Re,
                      TXi0Im,UmatOut,TmatOut);
        // SUM UP RESULTS, FOR ALL COLLOCATION POINTS
        for (unsigned int iEltColl=0; iEltColl<nEltColl; iEltColl++)
        {
          double sumutil=H[iXi]*M[nEltColl*iXi+iEltColl]*Jac[iXi];
          unsigned int rowBeg=2*iColl;
          unsigned int colBeg=2*EltCollIndex[iEltColl];
          for (unsigned int iGrSet=0; iGrSet<nGrSet; iGrSet++)
          {
            const unsigned int ind0 =nDof*nDof*iGrSet;
            const double theta=xi2theta[xiiXi2[iXi]];
            const double costheta=cos(theta);
            const double sintheta=sin(theta);
      
            // Axisymmetric analysis:
            // [(ugxx cost + ugxy sint)  ugxz]
            // [(ugzx cost + ugzy sint)  ugzz]
            URe[ind0+nDof*(colBeg+0)+rowBeg+0]+=2.0*sumutil*(costheta*UXiRe[9*iGrSet+0]+sintheta*UXiRe[9*iGrSet+1]);
            URe[ind0+nDof*(colBeg+1)+rowBeg+0]+=2.0*sumutil*(UXiRe[9*iGrSet+2]);
            URe[ind0+nDof*(colBeg+0)+rowBeg+1]+=2.0*sumutil*(costheta*UXiRe[9*iGrSet+6]+sintheta*UXiRe[9*iGrSet+7]);
            URe[ind0+nDof*(colBeg+1)+rowBeg+1]+=2.0*sumutil*(UXiRe[9*iGrSet+8]);
            if (ugCmplx)
            {
              UIm[ind0+nDof*(colBeg+0)+rowBeg+0]+=2.0*sumutil*(costheta*UXiIm[9*iGrSet+0]+sintheta*UXiIm[9*iGrSet+1]);
              UIm[ind0+nDof*(colBeg+1)+rowBeg+0]+=2.0*sumutil*(UXiIm[9*iGrSet+2]);
              UIm[ind0+nDof*(colBeg+0)+rowBeg+1]+=2.0*sumutil*(costheta*UXiIm[9*iGrSet+6]+sintheta*UXiIm[9*iGrSet+7]);
              UIm[ind0+nDof*(colBeg+1)+rowBeg+1]+=2.0*sumutil*(UXiIm[9*iGrSet+8]);
            }
            if (TmatOut)
            {
              TRe[ind0+nDof*(colBeg+0)+rowBeg+0]+=2.0*sumutil*(costheta*TXiRe[9*iGrSet+0]+sintheta*TXiRe[9*iGrSet+1]);
              TRe[ind0+nDof*(colBeg+1)+rowBeg+0]+=2.0*sumutil*(TXiRe[9*iGrSet+2]);
              TRe[ind0+nDof*(colBeg+0)+rowBeg+1]+=2.0*sumutil*(costheta*TXiRe[9*iGrSet+6]+sintheta*TXiRe[9*iGrSet+7]);
              TRe[ind0+nDof*(colBeg+1)+rowBeg+1]+=2.0*sumutil*(TXiRe[9*iGrSet+8]);
              if (tgCmplx)
              {
                TIm[ind0+nDof*(colBeg+0)+rowBeg+0]+=2.0*sumutil*(costheta*TXiIm[9*iGrSet+0]+sintheta*TXiIm[9*iGrSet+1]);
                TIm[ind0+nDof*(colBeg+1)+rowBeg+0]+=2.0*sumutil*(TXiIm[9*iGrSet+2]);
                TIm[ind0+nDof*(colBeg+0)+rowBeg+1]+=2.0*sumutil*(costheta*TXiIm[9*iGrSet+6]+sintheta*TXiIm[9*iGrSet+7]);
                TIm[ind0+nDof*(colBeg+1)+rowBeg+1]+=2.0*sumutil*(TXiIm[9*iGrSet+8]);
              }
              /* Account for singular part of Green's function on the
               *  diagonal terms.
               *  The static Green's fundamental solution is projected on
               *  The rigid body displacement. Therefore, the Green's functions
               *  may not be rotated: 
               */
              
              TRe[ind0+nDof*(rowBeg+0)+rowBeg+0]-=2.0*sumutil*TXi0Re[9*iGrSet+0];
              TRe[ind0+nDof*(rowBeg+1)+rowBeg+0]-=2.0*sumutil*TXi0Re[9*iGrSet+2];
              TRe[ind0+nDof*(rowBeg+0)+rowBeg+1]-=2.0*sumutil*TXi0Re[9*iGrSet+6];
              TRe[ind0+nDof*(rowBeg+1)+rowBeg+1]-=2.0*sumutil*TXi0Re[9*iGrSet+8];
              if (tgCmplx)
              {
                TIm[ind0+nDof*(rowBeg+0)+rowBeg+0]-=2.0*sumutil*TXi0Im[9*iGrSet+0];
                TIm[ind0+nDof*(rowBeg+1)+rowBeg+0]-=2.0*sumutil*TXi0Im[9*iGrSet+2];
                TIm[ind0+nDof*(rowBeg+0)+rowBeg+1]-=2.0*sumutil*TXi0Im[9*iGrSet+6];
                TIm[ind0+nDof*(rowBeg+1)+rowBeg+1]-=2.0*sumutil*TXi0Im[9*iGrSet+8];
              }
            }
          }
        }
      }
    }
  }
  delete [] EltNod;
  delete [] xi1;
  delete [] H1;
  delete [] xi2;
  delete [] H2;
  delete [] H;
  delete [] N1;
  delete [] M1;
  delete [] dN1;
  delete [] nat1;
  delete [] Jac1;
  delete [] normal1;
  delete [] xiCart;
  delete [] xi2theta;
  delete [] xiiXi2;
  delete [] Jac;
  delete [] normal;
  delete [] M;
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
