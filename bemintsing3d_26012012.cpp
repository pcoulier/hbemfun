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

//============================================================================//
//  THREE-DIMENSIONAL SINGULAR INTEGRATION
//============================================================================//
void bemintsing3d(const double* const Nod, const int& nNod,
                  const double* const Elt, const int& iElt,
                  const int& nElt,
                  const int* const  TypeID, const int* const  nKeyOpt,
                  const char* const TypeName[], const char* const TypeKeyOpts[], 
                  const int& nEltType,
                  const double* const Coll,const int& nColl, const int& iColl,  const int& iuniqueColl,
                  const int* const EltCollIndex, const int& nDof,
                  const double* const xiSing, const void* const* const greenPtr, 
                  const int& nGrSet, const bool& ugCmplx, 
                  const bool& tgCmplx, const bool& tg0Cmplx, double* const URe, 
                  double* const UIm, double* const TRe, double* const TIm, 
                  const bool& UmatOut,
				  const bool& TmatOut,
				  // const double* const s,
				  const bool& spassed,
				  const int& ms,
            	  const int& ns,
				  // const int* const scolli,
 				  const int* const scompi,
				  // const int* const uniquescolli,
				  // const int* const Nuniquescolli,
				  const int* const nuniquescolli,
				  const int* const uniquescolliind,
				  const int* const scollj,
				  const int* const scompj,
				  // const int* const uniquescollj,
				  // const int* const Nuniquescollj,
				  // const int* const nuniquescollj,
				  // const int* const uniquescolljind,
				  const bool* const InListuniquecollj,
				  const int& nuniquescollicumul,
				  const int* const inddiag,
				  const bool& ondiag,
				  const bool* const blockdiag,
				  const bool* const blocks,
				  const int& NEltCollConsider,
				  const int* const EltParent, const int* const nEltNod, const int* const nEltColl,
				  const int* const EltShapeN, const int* const EltShapeM, const int* const EltDim,
				  const int* const AxiSym, const int* const Periodic, const int* const nGauss,
				  const int* const nEltDiv, const int* const nGaussSing, const int* const nEltDivSing,
				  const double* const NodCoord)
{
  const int EltType = int(Elt[nElt+iElt]);
  /*
  int Parent;
  int nEltNod;
  int nEltColl;
  int ShapeTypeN;
  int ShapeTypeM;
  int EltDim;
  int AxiSym;
  int Periodic;
  int nGauss;
  int nEltDiv;
  int nGaussSing;
  int nEltDivSing;
  eltdef(EltType,TypeID,TypeName,TypeKeyOpts,nKeyOpt,nEltType,Parent,nEltNod,
         nEltColl,ShapeTypeN,ShapeTypeM,EltDim,AxiSym,Periodic,nGauss,nEltDiv,
         nGaussSing,nEltDivSing);
  */		 
  const int nXi=nEltDivSing[iElt]*nEltDivSing[iElt]*nGaussSing[iElt]*nGaussSing[iElt];

  /*
  int NodIndex;
  int NodID;
  double* const NodCoord =new(nothrow) double[3*nEltNod[iElt]];
  if (NodCoord==0) throw("Out of memory.");
 
  
  for (int iEltNod=0; iEltNod<nEltNod[iElt]; iEltNod++)
  {
    NodID=int(Elt[(2+iEltNod)*nElt+iElt]);
    BemNodeIndex(Nod,nNod,NodID,NodIndex);
    NodCoord[0*nEltNod[iElt]+iEltNod]=Nod[1*nNod+NodIndex];
    NodCoord[1*nEltNod[iElt]+iEltNod]=Nod[2*nNod+NodIndex];
    NodCoord[2*nEltNod[iElt]+iEltNod]=Nod[3*nNod+NodIndex];
  }   
  */
	
  // ELEMENT TRIANGLE DIVISION
  int nDiv;
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
  triangdiv(xiSing,EltParent[iElt],nDiv,am,a1,a2,rhom,rho1,rho2);

  // H and xi in v1,v2
  double* const v=new(nothrow) double[2*nXi];
  if (v==0) throw("Out of memory.");
  double* const H=new(nothrow) double[nXi];
  if (H==0) throw("Out of memory.");
  gausspw2D(nEltDivSing[iElt],nGaussSing[iElt],v,H);

  double* const a=new(nothrow) double[nXi];
  if (a==0) throw("Out of memory.");
  double* const rho=new(nothrow) double[nXi];
  if (rho==0) throw("Out of memory.");

  double* const xi=new(nothrow) double[2*nXi];
  if (xi==0) throw("Out of memory.");
  double* const N=new(nothrow) double[nXi*nEltNod[iElt]];
  if (N==0) throw("Out of memory.");
  double* const M=new(nothrow) double[nXi*nEltColl[iElt]];
  if (M==0) throw("Out of memory.");
  double* const Mmod=new(nothrow) double[nXi*nEltColl[iElt]];
  if (Mmod==0) throw("Out of memory.");
  double* const dN=new(nothrow) double[2*nXi*nEltNod[iElt]];
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
  
  for (int iComp=0; iComp<9*nGrSet;iComp++)
  {
    TXiRe[iComp]=0.0;
    TXiIm[iComp]=0.0;
    TXi0Re[iComp]=0.0;
    TXi0Im[iComp]=0.0;
  }

  // Initialize interpolation of Green's function
  int r1=0;
  int r2=1;
  bool extrapFlag=false;
  double* const interpr=new(nothrow) double[2];
  if (interpr==0) throw("Out of memory.");
  int z1=0;
  int z2=1;
  double* const interpz=new(nothrow) double[2];
  if (interpz==0) throw("Out of memory.");
  int zs1=0;

  // mexPrintf("Running bemintsing3d...  \n "); // DEBUG
  
 
  // if (s!=0)
  if (spassed)
  {
  // mexPrintf("ondiag: %s\n",ondiag ? "true" : "false"); // DEBUG
  
   for (int iDiv=0; iDiv<nDiv; iDiv++) if (rhom[iDiv]>1e-10)
  {
    for (int iXi=0; iXi<nXi; iXi++)
    {
      a[iXi]  =0.5*(a2[iDiv]-a1[iDiv])*v[nXi+iXi]+0.5*(a2[iDiv]+a1[iDiv]);
      rho[iXi]=0.5*rhom[iDiv]/cos(a[iXi]-am[iDiv])*(1+v[iXi]);
      xi[iXi]=xiSing[0]+rho[iXi]*cos(a[iXi]);
      xi[nXi+iXi]=xiSing[1]+rho[iXi]*sin(a[iXi]);
    }
    
    shapefun(EltShapeN[iElt],nXi,xi,N);
    shapefun(EltShapeM[iElt],nXi,xi,M);
    shapederiv(EltShapeN[iElt],nXi,xi,dN);
    shapenatcoord(dN,nEltNod[iElt],nXi,NodCoord,nat,EltDim[iElt]);
    jacobian(nat,nXi,Jac,EltDim[iElt]);
    if (TmatOut) bemnormal(nat,nXi,EltDim[iElt],normal);

    for (int icomp=0; icomp<3*nXi; icomp++) xiCart[icomp]=0.0;
    for (int iXi=0; iXi<nXi; iXi++)
    {
      for (int iNod=0; iNod<nEltNod[iElt]; iNod++)
      {
        xiCart[3*iXi+0]+=N[nEltNod[iElt]*iXi+iNod]*NodCoord[0*nEltNod[iElt]+iNod];
        xiCart[3*iXi+1]+=N[nEltNod[iElt]*iXi+iNod]*NodCoord[1*nEltNod[iElt]+iNod];
        xiCart[3*iXi+2]+=N[nEltNod[iElt]*iXi+iNod]*NodCoord[2*nEltNod[iElt]+iNod];
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
      int nEltCollConsider=0;
	  
      for (int iEltColl=0; iEltColl<nEltColl[iElt]; iEltColl++)
      {
	  
	  // mexPrintf("ondiag : %s\n",ondiag ? "true" : "false"); // DEBUG
	  // // for (int i=0; i<Nuniquescolli[0]; i++)
	// {
		// mexPrintf("blockdiag[iuniqueColl][%d]: %s \n",iuniqueColl, blockdiag[iuniqueColl] ? "true": "false");
	// }
	
	  
	  // if (ondiag==true && blockdiag[iColl])
	  if (ondiag==true && blockdiag[iuniqueColl])
	  {
		// mexPrintf(" in loop ... \n");
		
		double sumutil=rho[iXi]*H[iXi]*M[nEltColl[iElt]*iXi+iEltColl]*Jac[iXi]*0.25*(a2[iDiv]-a1[iDiv])*rhom[iDiv]/cos(a[iXi]-am[iDiv]);
		for (int iGrSet=0; iGrSet<nGrSet; iGrSet++)
		{
		const int ind0=ms*ns*iGrSet;
		
		// // // TRe[ind0+inddiag[9*iColl+0]]-=sumutil*TXi0Re[9*iGrSet+0];																			 
		// // // TRe[ind0+inddiag[9*iColl+1]]-=sumutil*TXi0Re[9*iGrSet+1];																			 
		// // // TRe[ind0+inddiag[9*iColl+2]]-=sumutil*TXi0Re[9*iGrSet+2];																			 
		// // // TRe[ind0+inddiag[9*iColl+3]]-=sumutil*TXi0Re[9*iGrSet+3];																			 
		// // // TRe[ind0+inddiag[9*iColl+4]]-=sumutil*TXi0Re[9*iGrSet+4];																			 
		// // // TRe[ind0+inddiag[9*iColl+5]]-=sumutil*TXi0Re[9*iGrSet+5];																			 
		// // // TRe[ind0+inddiag[9*iColl+6]]-=sumutil*TXi0Re[9*iGrSet+6];																			 
		// // // TRe[ind0+inddiag[9*iColl+7]]-=sumutil*TXi0Re[9*iGrSet+7];																			 
		// // // TRe[ind0+inddiag[9*iColl+8]]-=sumutil*TXi0Re[9*iGrSet+8];																			 
		if (TmatOut)
        {
		TRe[ind0+inddiag[9*iuniqueColl+0]]-=sumutil*TXi0Re[9*iGrSet+0];																			 
		TRe[ind0+inddiag[9*iuniqueColl+1]]-=sumutil*TXi0Re[9*iGrSet+1];																			 
		TRe[ind0+inddiag[9*iuniqueColl+2]]-=sumutil*TXi0Re[9*iGrSet+2];																			 
		TRe[ind0+inddiag[9*iuniqueColl+3]]-=sumutil*TXi0Re[9*iGrSet+3];																			 
		TRe[ind0+inddiag[9*iuniqueColl+4]]-=sumutil*TXi0Re[9*iGrSet+4];																			 
		TRe[ind0+inddiag[9*iuniqueColl+5]]-=sumutil*TXi0Re[9*iGrSet+5];																			 
		TRe[ind0+inddiag[9*iuniqueColl+6]]-=sumutil*TXi0Re[9*iGrSet+6];																			 
		TRe[ind0+inddiag[9*iuniqueColl+7]]-=sumutil*TXi0Re[9*iGrSet+7];																			 
		TRe[ind0+inddiag[9*iuniqueColl+8]]-=sumutil*TXi0Re[9*iGrSet+8];																			 
		}
		
		}
	  }
	  
	  
	  else
	  {
		if (InListuniquecollj[EltCollIndex[iEltColl]])
		{
		
		if (blocks[iuniqueColl]==true && ms>1)
		{
		
		// mexPrintf("Adrie Koster rules... \n");
		double sumutil=rho[iXi]*H[iXi]*M[nEltColl[iElt]*iXi+iEltColl]*Jac[iXi]*0.25*(a2[iDiv]-a1[iDiv])*rhom[iDiv]/cos(a[iXi]-am[iDiv]);
		int rowBeg=3*iuniqueColl;
        int colBeg=3*EltCollIndex[iEltColl]-3*(NEltCollConsider+nEltCollConsider);
		
		// if (iXi==0 && iDiv==0)
		// {
		// mexPrintf("rowBeg: %d \n",rowBeg);
		// mexPrintf("colBeg: %d \n",colBeg);
		// }
		
        for (int iGrSet=0; iGrSet<nGrSet; iGrSet++)
        {
          int ind0=ms*ns*iGrSet;
			
          URe[ind0+ms*(colBeg+0)+rowBeg+0]+=sumutil*UXiRe[9*iGrSet+0];
          URe[ind0+ms*(colBeg+1)+rowBeg+0]+=sumutil*UXiRe[9*iGrSet+1];
          URe[ind0+ms*(colBeg+2)+rowBeg+0]+=sumutil*UXiRe[9*iGrSet+2];
          URe[ind0+ms*(colBeg+0)+rowBeg+1]+=sumutil*UXiRe[9*iGrSet+3];
          URe[ind0+ms*(colBeg+1)+rowBeg+1]+=sumutil*UXiRe[9*iGrSet+4];
          URe[ind0+ms*(colBeg+2)+rowBeg+1]+=sumutil*UXiRe[9*iGrSet+5];
          URe[ind0+ms*(colBeg+0)+rowBeg+2]+=sumutil*UXiRe[9*iGrSet+6];
          URe[ind0+ms*(colBeg+1)+rowBeg+2]+=sumutil*UXiRe[9*iGrSet+7];
          URe[ind0+ms*(colBeg+2)+rowBeg+2]+=sumutil*UXiRe[9*iGrSet+8];
          if (ugCmplx)
          {
            UIm[ind0+ms*(colBeg+0)+rowBeg+0]+=sumutil*UXiIm[9*iGrSet+0];
            UIm[ind0+ms*(colBeg+1)+rowBeg+0]+=sumutil*UXiIm[9*iGrSet+1];
            UIm[ind0+ms*(colBeg+2)+rowBeg+0]+=sumutil*UXiIm[9*iGrSet+2];
            UIm[ind0+ms*(colBeg+0)+rowBeg+1]+=sumutil*UXiIm[9*iGrSet+3];
            UIm[ind0+ms*(colBeg+1)+rowBeg+1]+=sumutil*UXiIm[9*iGrSet+4];
            UIm[ind0+ms*(colBeg+2)+rowBeg+1]+=sumutil*UXiIm[9*iGrSet+5];
            UIm[ind0+ms*(colBeg+0)+rowBeg+2]+=sumutil*UXiIm[9*iGrSet+6];
            UIm[ind0+ms*(colBeg+1)+rowBeg+2]+=sumutil*UXiIm[9*iGrSet+7];
            UIm[ind0+ms*(colBeg+2)+rowBeg+2]+=sumutil*UXiIm[9*iGrSet+8];
          }
          if (TmatOut)
          {
            TRe[ind0+ms*(colBeg+0)+rowBeg+0]+=sumutil*TXiRe[9*iGrSet+0];  // txx
            TRe[ind0+ms*(colBeg+1)+rowBeg+0]+=sumutil*TXiRe[9*iGrSet+1];  // txy
            TRe[ind0+ms*(colBeg+2)+rowBeg+0]+=sumutil*TXiRe[9*iGrSet+2];  // txz
            TRe[ind0+ms*(colBeg+0)+rowBeg+1]+=sumutil*TXiRe[9*iGrSet+3];  // tyx
            TRe[ind0+ms*(colBeg+1)+rowBeg+1]+=sumutil*TXiRe[9*iGrSet+4];  // tyy
            TRe[ind0+ms*(colBeg+2)+rowBeg+1]+=sumutil*TXiRe[9*iGrSet+5];  // tyz
            TRe[ind0+ms*(colBeg+0)+rowBeg+2]+=sumutil*TXiRe[9*iGrSet+6];  // tzx
            TRe[ind0+ms*(colBeg+1)+rowBeg+2]+=sumutil*TXiRe[9*iGrSet+7];  // tzy
            TRe[ind0+ms*(colBeg+2)+rowBeg+2]+=sumutil*TXiRe[9*iGrSet+8];  // tzz
            if (tgCmplx)
            {
              TIm[ind0+ms*(colBeg+0)+rowBeg+0]+=sumutil*TXiIm[9*iGrSet+0];
              TIm[ind0+ms*(colBeg+1)+rowBeg+0]+=sumutil*TXiIm[9*iGrSet+1];
              TIm[ind0+ms*(colBeg+2)+rowBeg+0]+=sumutil*TXiIm[9*iGrSet+2];
              TIm[ind0+ms*(colBeg+0)+rowBeg+1]+=sumutil*TXiIm[9*iGrSet+3];
              TIm[ind0+ms*(colBeg+1)+rowBeg+1]+=sumutil*TXiIm[9*iGrSet+4];
              TIm[ind0+ms*(colBeg+2)+rowBeg+1]+=sumutil*TXiIm[9*iGrSet+5];
              TIm[ind0+ms*(colBeg+0)+rowBeg+2]+=sumutil*TXiIm[9*iGrSet+6];
              TIm[ind0+ms*(colBeg+1)+rowBeg+2]+=sumutil*TXiIm[9*iGrSet+7];
              TIm[ind0+ms*(colBeg+2)+rowBeg+2]+=sumutil*TXiIm[9*iGrSet+8];
            }
		  }
		  
		}
		

		}
		
		
		
		else
		{
		
        double sumutil=rho[iXi]*H[iXi]*M[nEltColl[iElt]*iXi+iEltColl]*Jac[iXi]*0.25*(a2[iDiv]-a1[iDiv])*rhom[iDiv]/cos(a[iXi]-am[iDiv]);
			
			// for (int iuniquescolliind=0; iuniquescolliind<nuniquescolli[iColl]; iuniquescolliind++)
			for (int iuniquescolliind=0; iuniquescolliind<nuniquescolli[iuniqueColl]; iuniquescolliind++)
				{
				
				
				if (scollj[uniquescolliind[nuniquescollicumul+iuniquescolliind]] == EltCollIndex[iEltColl])
					{
					
					
					
					
					
					for (int iGrSet=0; iGrSet<nGrSet; iGrSet++)
					{
					const int ind0=ms*ns*iGrSet;
					
					if (ondiag==false)
					{
					
					// // mexPrintf("ind0+uniquescolliind[nuniquescollicumul+iuniquescolliind]: %d\n",ind0+uniquescolliind[nuniquescollicumul+iuniquescolliind]); // DEBUG
					// // mexPrintf("nuniquescollicumul+iuniquescolliind: %d\n",nuniquescollicumul+iuniquescolliind); // DEBUG
					// // mexPrintf("3*scompi[uniquescolliind[nuniquescollicumul+iuniquescolliind]]+scompj[uniquescolliind[nuniquescollicumul+iuniquescolliind]]: %d \n",3*scompi[uniquescolliind[nuniquescollicumul+iuniquescolliind]]+scompj[uniquescolliind[nuniquescollicumul+iuniquescolliind]]);
					// // mexPrintf("uniquescolliind[nuniquescollicumul+iuniquescolliind]: %d \n",uniquescolliind[nuniquescollicumul+iuniquescolliind]);
					if (UmatOut)
					{
					URe[ind0+uniquescolliind[nuniquescollicumul+iuniquescolliind]]+=sumutil*UXiRe[9*iGrSet+3*scompi[uniquescolliind[nuniquescollicumul+iuniquescolliind]]+
																											 scompj[uniquescolliind[nuniquescollicumul+iuniquescolliind]]];
																											 
					if (ugCmplx)
					{	
						UIm[ind0+uniquescolliind[nuniquescollicumul+iuniquescolliind]]+=sumutil*UXiIm[9*iGrSet+3*scompi[uniquescolliind[nuniquescollicumul+iuniquescolliind]]+
																											     scompj[uniquescolliind[nuniquescollicumul+iuniquescolliind]]];
					}
					}
					
					}

					if (TmatOut)
					{
						if (ondiag==false)
						{	
						// mexPrintf("test...\n");
						// mexPrintf("ind0+uniquescolliind[nuniquescollicumul+iuniquescolliind]: %d \n",ind0+uniquescolliind[nuniquescollicumul+iuniquescolliind]);
						// mexPrintf("lala...: %d \n",9*iGrSet+3*scompi[uniquescolliind[nuniquescollicumul+iuniquescolliind]]+scompj[uniquescolliind[nuniquescollicumul+iuniquescolliind]]);
						
						TRe[ind0+uniquescolliind[nuniquescollicumul+iuniquescolliind]]+=sumutil*TXiRe[9*iGrSet+3*scompi[uniquescolliind[nuniquescollicumul+iuniquescolliind]]+
																											 scompj[uniquescolliind[nuniquescollicumul+iuniquescolliind]]];
						// mexPrintf("aftertest...\n");
																											 
						 if (tgCmplx)
						{
							TIm[ind0+uniquescolliind[nuniquescollicumul+iuniquescolliind]]+=sumutil*TXiIm[9*iGrSet+3*scompi[uniquescolliind[nuniquescollicumul+iuniquescolliind]]+
																											 scompj[uniquescolliind[nuniquescollicumul+iuniquescolliind]]];
						}
						
						}
						
						// double temp = 	scollj[uniquescolliind[nuniquescollicumul+iuniquescolliind]];
												// mexPrintf("temp == Coll[nColl+iColl]: %s", temp == Coll[iColl] ? true : false);
						// if ( temp == Coll[iColl])
						// if ( temp == 0 && Coll[iColl]==0)
						// if ( temp == 0 )
						// / if ((double) scollj[uniquescolliind[nuniquescollicumul+iuniquescolliind]] == Coll[nColl+iColl]-1)
												
							// mexPrintf("sing...\n");	
						
						// if (ondiag==true && inddiag[9*iColl+3*scompi[uniquescolliind[nuniquescollicumul+iuniquescolliind]]+scompj[uniquescolliind[nuniquescollicumul+iuniquescolliind]]]!=-1)
						if (ondiag==true && inddiag[9*iuniqueColl+3*scompi[uniquescolliind[nuniquescollicumul+iuniquescolliind]]+scompj[uniquescolliind[nuniquescollicumul+iuniquescolliind]]]!=-1)
						// // // // // // // // // // if (ondiag==true && inddiag[9*iColl+scompi[uniquescolliind[nuniquescollicumul+iuniquescolliind]]+3*scompj[uniquescolliind[nuniquescollicumul+iuniquescolliind]]]!=-1)
						{
						// TRe[ind0+uniquescolliind[nuniquescollicumul+iuniquescolliind]]-=sumutil*TXi0Re[9*iGrSet+3*scompi[uniquescolliind[nuniquescollicumul+iuniquescolliind]]+
																											 // scompj[uniquescolliind[nuniquescollicumul+iuniquescolliind]]];
						
						// TRe[ind0+inddiag[9*iColl+3*scompi[uniquescolliind[nuniquescollicumul+iuniquescolliind]]+scompj[uniquescolliind[nuniquescollicumul+iuniquescolliind]]]]-=sumutil*TXi0Re[9*iGrSet+3*scompi[uniquescolliind[nuniquescollicumul+iuniquescolliind]]+scompj[uniquescolliind[nuniquescollicumul+iuniquescolliind]]];																			 
						TRe[ind0+inddiag[9*iuniqueColl+3*scompi[uniquescolliind[nuniquescollicumul+iuniquescolliind]]+scompj[uniquescolliind[nuniquescollicumul+iuniquescolliind]]]]-=sumutil*TXi0Re[9*iGrSet+3*scompi[uniquescolliind[nuniquescollicumul+iuniquescolliind]]+scompj[uniquescolliind[nuniquescollicumul+iuniquescolliind]]];																			 
						// // // // // // // // // // TRe[ind0+inddiag[9*iColl+scompi[uniquescolliind[nuniquescollicumul+iuniquescolliind]]+3*scompj[uniquescolliind[nuniquescollicumul+iuniquescolliind]]]]-=sumutil*TXi0Re[9*iGrSet+3*scompi[uniquescolliind[nuniquescollicumul+iuniquescolliind]]+scompj[uniquescolliind[nuniquescollicumul+iuniquescolliind]]];																			 
							
						if (tgCmplx)
						{
						
							TIm[ind0+inddiag[9*iColl+3*scompi[uniquescolliind[nuniquescollicumul+iuniquescolliind]]+scompj[uniquescolliind[nuniquescollicumul+iuniquescolliind]]]]-=sumutil*TXi0Im[9*iGrSet+3*scompi[uniquescolliind[nuniquescollicumul+iuniquescolliind]]+scompj[uniquescolliind[nuniquescollicumul+iuniquescolliind]]];
							// // // // // // // // // // TIm[ind0+inddiag[9*iColl+scompi[uniquescolliind[nuniquescollicumul+iuniquescolliind]]+3*scompj[uniquescolliind[nuniquescollicumul+iuniquescolliind]]]]-=sumutil*TXi0Im[9*iGrSet+3*scompi[uniquescolliind[nuniquescollicumul+iuniquescolliind]]+scompj[uniquescolliind[nuniquescollicumul+iuniquescolliind]]];
							// TIm[ind0+uniquescolliind[nuniquescollicumul+iuniquescolliind]]-=sumutil*TXi0Im[9*iGrSet+3*scompi[uniquescolliind[nuniquescollicumul+iuniquescolliind]]+
																											 // scompj[uniquescolliind[nuniquescollicumul+iuniquescolliind]]];
							
						}
											
						// TRe[ind0+inddiag[3*scompi[uniquescolliind[nuniquescollicumul+iuniquescolliind]]+scompj[uniquescolliind[nuniquescollicumul+iuniquescolliind]]]]-=sumutil*TXi0Re[9*iGrSet+3*scompi[uniquescolliind[nuniquescollicumul+iuniquescolliind]]+scompj[uniquescolliind[nuniquescollicumul+iuniquescolliind]]];
						
						}
						
						
						
						
					}
																											 
					}																									
				
					}
				}
				
				
			}	
				
				
				
		
		}
		else
		{
		nEltCollConsider++;
		}
		
		}
	
		
		/*
		double sumutil=rho[iXi]*H[iXi]*M[nEltColl*iXi+iEltColl]*Jac[iXi]*0.25*(a2[iDiv]-a1[iDiv])*rhom[iDiv]/cos(a[iXi]-am[iDiv]);
		int rowBeg=3*iColl;
        int colBeg=3*EltCollIndex[iEltColl];
        for (int iGrSet=0; iGrSet<nGrSet; iGrSet++)
        {
          int ind0 =nDof*nDof*iGrSet;
          URe[ind0+nDof*(colBeg+0)+rowBeg+0]+=sumutil*UXiRe[9*iGrSet+0];
          URe[ind0+nDof*(colBeg+1)+rowBeg+0]+=sumutil*UXiRe[9*iGrSet+1];
          URe[ind0+nDof*(colBeg+2)+rowBeg+0]+=sumutil*UXiRe[9*iGrSet+2];
          URe[ind0+nDof*(colBeg+0)+rowBeg+1]+=sumutil*UXiRe[9*iGrSet+3];
          URe[ind0+nDof*(colBeg+1)+rowBeg+1]+=sumutil*UXiRe[9*iGrSet+4];
          URe[ind0+nDof*(colBeg+2)+rowBeg+1]+=sumutil*UXiRe[9*iGrSet+5];
          URe[ind0+nDof*(colBeg+0)+rowBeg+2]+=sumutil*UXiRe[9*iGrSet+6];
          URe[ind0+nDof*(colBeg+1)+rowBeg+2]+=sumutil*UXiRe[9*iGrSet+7];
          URe[ind0+nDof*(colBeg+2)+rowBeg+2]+=sumutil*UXiRe[9*iGrSet+8];
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
            TRe[ind0+nDof*(rowBeg+0)+rowBeg+0]-=sumutil*TXi0Re[9*iGrSet+0];
            TRe[ind0+nDof*(rowBeg+1)+rowBeg+0]-=sumutil*TXi0Re[9*iGrSet+1];
            TRe[ind0+nDof*(rowBeg+2)+rowBeg+0]-=sumutil*TXi0Re[9*iGrSet+2];
            TRe[ind0+nDof*(rowBeg+0)+rowBeg+1]-=sumutil*TXi0Re[9*iGrSet+3];
            TRe[ind0+nDof*(rowBeg+1)+rowBeg+1]-=sumutil*TXi0Re[9*iGrSet+4];
            TRe[ind0+nDof*(rowBeg+2)+rowBeg+1]-=sumutil*TXi0Re[9*iGrSet+5];
            TRe[ind0+nDof*(rowBeg+0)+rowBeg+2]-=sumutil*TXi0Re[9*iGrSet+6];
            TRe[ind0+nDof*(rowBeg+1)+rowBeg+2]-=sumutil*TXi0Re[9*iGrSet+7];
            TRe[ind0+nDof*(rowBeg+2)+rowBeg+2]-=sumutil*TXi0Re[9*iGrSet+8];
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
		*/
		
      }
    }
  }
  
  
  }
  
  
  else
  {

  for (int iDiv=0; iDiv<nDiv; iDiv++) if (rhom[iDiv]>1e-10)
  {
    for (int iXi=0; iXi<nXi; iXi++)
    {
      a[iXi]  =0.5*(a2[iDiv]-a1[iDiv])*v[nXi+iXi]+0.5*(a2[iDiv]+a1[iDiv]);
      rho[iXi]=0.5*rhom[iDiv]/cos(a[iXi]-am[iDiv])*(1+v[iXi]);
      xi[iXi]=xiSing[0]+rho[iXi]*cos(a[iXi]);
      xi[nXi+iXi]=xiSing[1]+rho[iXi]*sin(a[iXi]);
    }
    
    shapefun(EltShapeN[iElt],nXi,xi,N);
    shapefun(EltShapeM[iElt],nXi,xi,M);
    shapederiv(EltShapeN[iElt],nXi,xi,dN);
    shapenatcoord(dN,nEltNod[iElt],nXi,NodCoord,nat,EltDim[iElt]);
    jacobian(nat,nXi,Jac,EltDim[iElt]);
    if (TmatOut) bemnormal(nat,nXi,EltDim[iElt],normal);

    for (int icomp=0; icomp<3*nXi; icomp++) xiCart[icomp]=0.0;
    for (int iXi=0; iXi<nXi; iXi++)
    {
      for (int iNod=0; iNod<nEltNod[iElt]; iNod++)
      {
        xiCart[3*iXi+0]+=N[nEltNod[iElt]*iXi+iNod]*NodCoord[0*nEltNod[iElt]+iNod];
        xiCart[3*iXi+1]+=N[nEltNod[iElt]*iXi+iNod]*NodCoord[1*nEltNod[iElt]+iNod];
        xiCart[3*iXi+2]+=N[nEltNod[iElt]*iXi+iNod]*NodCoord[2*nEltNod[iElt]+iNod];
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
      
      for (int iEltColl=0; iEltColl<nEltColl[iElt]; iEltColl++)
      {
        double sumutil=rho[iXi]*H[iXi]*M[nEltColl[iElt]*iXi+iEltColl]*Jac[iXi]*0.25*(a2[iDiv]-a1[iDiv])*rhom[iDiv]/cos(a[iXi]-am[iDiv]);
        int rowBeg=3*iColl;
        int colBeg=3*EltCollIndex[iEltColl];
        for (int iGrSet=0; iGrSet<nGrSet; iGrSet++)
        {
          int ind0 =nDof*nDof*iGrSet;
          URe[ind0+nDof*(colBeg+0)+rowBeg+0]+=sumutil*UXiRe[9*iGrSet+0];
          URe[ind0+nDof*(colBeg+1)+rowBeg+0]+=sumutil*UXiRe[9*iGrSet+1];
          URe[ind0+nDof*(colBeg+2)+rowBeg+0]+=sumutil*UXiRe[9*iGrSet+2];
          URe[ind0+nDof*(colBeg+0)+rowBeg+1]+=sumutil*UXiRe[9*iGrSet+3];
          URe[ind0+nDof*(colBeg+1)+rowBeg+1]+=sumutil*UXiRe[9*iGrSet+4];
          URe[ind0+nDof*(colBeg+2)+rowBeg+1]+=sumutil*UXiRe[9*iGrSet+5];
          URe[ind0+nDof*(colBeg+0)+rowBeg+2]+=sumutil*UXiRe[9*iGrSet+6];
          URe[ind0+nDof*(colBeg+1)+rowBeg+2]+=sumutil*UXiRe[9*iGrSet+7];
          URe[ind0+nDof*(colBeg+2)+rowBeg+2]+=sumutil*UXiRe[9*iGrSet+8];
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
            TRe[ind0+nDof*(rowBeg+0)+rowBeg+0]-=sumutil*TXi0Re[9*iGrSet+0];
            TRe[ind0+nDof*(rowBeg+1)+rowBeg+0]-=sumutil*TXi0Re[9*iGrSet+1];
            TRe[ind0+nDof*(rowBeg+2)+rowBeg+0]-=sumutil*TXi0Re[9*iGrSet+2];
            TRe[ind0+nDof*(rowBeg+0)+rowBeg+1]-=sumutil*TXi0Re[9*iGrSet+3];
            TRe[ind0+nDof*(rowBeg+1)+rowBeg+1]-=sumutil*TXi0Re[9*iGrSet+4];
            TRe[ind0+nDof*(rowBeg+2)+rowBeg+1]-=sumutil*TXi0Re[9*iGrSet+5];
            TRe[ind0+nDof*(rowBeg+0)+rowBeg+2]-=sumutil*TXi0Re[9*iGrSet+6];
            TRe[ind0+nDof*(rowBeg+1)+rowBeg+2]-=sumutil*TXi0Re[9*iGrSet+7];
            TRe[ind0+nDof*(rowBeg+2)+rowBeg+2]-=sumutil*TXi0Re[9*iGrSet+8];
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
  
  // delete [] NodCoord;
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


