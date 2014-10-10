#include "eltdef.h"
#include "bemcollpoints.h"
#include "gausspw.h"
#include "shapefun.h"
#include "bemnormal.h"
#include "greeneval3d.h"
#include "greenrotate3d.h"
#include <math.h>
#include <time.h>
#include <new>

#include "mex.h"

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
// THREE-DIMENSIONAL REGULAR BOUNDARY ELEMENT INTEGRATION
//======================================================================
void bemintreg3d(const double* const Nod, const unsigned int& nNod,
                 const double* const Elt, const unsigned int& iElt, const unsigned int& nElt,
                 const unsigned int* const  TypeID, const unsigned int* const nKeyOpt,
                 const char* const TypeName[], const char* const TypeKeyOpts[],
                 const unsigned int& nEltType, const double* const Coll,
                 const unsigned int& nColl, const unsigned int* const RegularColl,
                 const unsigned int* const EltCollIndex, const unsigned int& nDof,
                 const void* const* const greenPtr, const unsigned int& nGrSet,
                 const bool& ugCmplx, const bool& tgCmplx,
                 const bool& tg0Cmplx, double* const URe, double* const UIm,
                 double* const TRe, double* const TIm, const bool UmatOut,const bool TmatOut,
				 const double* const s,
				 const unsigned int& ms,
            	 const unsigned int& ns,
				 // const unsigned int* const scolli,
 				 const unsigned int* const scompi,
				 const unsigned int* const uniquescolli,
				 const unsigned int* const Nuniquescolli,
				 const unsigned int* const nuniquescolli,
				 const unsigned int* const uniquescolliind,
				 const unsigned int* const scollj,
				 const unsigned int* const scompj,
				 // const unsigned int* const uniquescollj,
				 // const unsigned int* const Nuniquescollj,
				 // const unsigned int* const nuniquescollj,
				 // const unsigned int* const uniquescolljind,
				 const bool* const InListuniquecollj,
				 const unsigned int* const inddiag,
				 const bool& ondiag)
{
  // ELEMENT PROPERTIES
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

  // NUMBER OF GAUSSIAN POINTS
  unsigned int nXi;
  if (Parent == 1) nXi=nGauss;
  else if (Parent == 2) nXi=nEltDiv*nEltDiv*nGauss*nGauss;

  int NodIndex;
  unsigned int NodID;
  double* const EltNod =new(nothrow) double[3*nEltNod];
  if (EltNod==0) throw("Out of memory.");

  // DETERMINE COORDINATES OF ELEMENT NODES (OF ELEMENT IELT)
  for (unsigned int iEltNod=0; iEltNod<nEltNod; iEltNod++)
  {
        NodID=(unsigned int)(Elt[(2+iEltNod)*nElt+iElt]);
        BemNodeIndex(Nod,nNod,NodID,NodIndex);
        EltNod[0*nEltNod+iEltNod]=Nod[1*nNod+NodIndex];
        EltNod[1*nEltNod+iEltNod]=Nod[2*nNod+NodIndex];
        EltNod[2*nEltNod+iEltNod]=Nod[3*nNod+NodIndex];
  }

  // DETERMINE SAMPLE POINTS FOR THE ELEMENT TYPE (NumGauss * NumEltDiv)
  double* const xi=new(nothrow) double[2*nXi];
    if (xi==0) throw("Out of memory.");
  double* const H=new(nothrow) double[nXi];
    if (H==0) throw("Out of memory.");

  if (Parent == 1) gausspwtri(nGauss,xi,H);
  else gausspw2D(nEltDiv,nGauss,xi,H);

  // unsigned int nEltCollToConsider=0;
  // if (s!=0)
  // {
     // for (unsigned int iEltColl=0; iEltColl<nEltColl; iEltColl++)
     // {
	  	 // if (InListuniquecollj[EltCollIndex[iEltColl]])
	  	 // {
		  	// nEltCollToConsider++;
		 // }
     // }
  // }
  // else
  // {
   	   // nEltCollToConsider=nEltColl;
  // }
		  											   
  // SHAPE FUNCTIONS IN THE SAMPLE POINTS
  double* const N=new(nothrow) double[nXi*nEltNod];
    if (N==0) throw("Out of memory.");
 double* const M=new(nothrow) double[nXi*nEltColl];
   if (M==0) throw("Out of memory.");
  // double* const M=new(nothrow) double[nXi*nEltCollToConsider];
     // if (M==0) throw("Out of memory.");  
  double* const dN=new(nothrow) double[2*nXi*nEltNod];
    if (dN==0) throw("Out of memory.");
  double* const nat=new(nothrow) double[6*nXi];
    if (nat==0) throw("Out of memory.");
  double* const Jac=new(nothrow) double[nXi];
    if (Jac==0) throw("Out of memory.");
  double* const xiCart=new(nothrow) double[3*nXi];
    if (xiCart==0) throw("Out of memory.");
  double* const normal=new(nothrow) double[3*nXi];
    if (normal==0) throw("Out of memory.");

  shapefun(ShapeTypeN,nXi,xi,N);
  shapefun(ShapeTypeM,nXi,xi,M);
  
  shapederiv(ShapeTypeN,nXi,xi,dN);
  shapenatcoord(dN,nEltNod,nXi,EltNod,nat,EltDim);
  jacobian(nat,nXi,Jac,EltDim);
  if (TmatOut) bemnormal(nat,nXi,EltDim,normal);
  

  // NODAL COORDINATES
  for (unsigned int icomp=0; icomp<3*nXi; icomp++) xiCart[icomp]=0.0;
  for (unsigned int iXi=0; iXi<nXi; iXi++)
  {
    for (unsigned int iEltNod=0; iEltNod<nEltNod; iEltNod++)
    {
      xiCart[3*iXi+0]+=N[nEltNod*iXi+iEltNod]*EltNod[0*nEltNod+iEltNod];
      xiCart[3*iXi+1]+=N[nEltNod*iXi+iEltNod]*EltNod[1*nEltNod+iEltNod];
      xiCart[3*iXi+2]+=N[nEltNod*iXi+iEltNod]*EltNod[2*nEltNod+iEltNod];
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

  
  //  s  not empty
  if (s!=0)
  {	 
	unsigned int nuniquescollicumul=0;
  
	for (unsigned int iuniquescolli=0; iuniquescolli<Nuniquescolli[0]; iuniquescolli++) // bepalen welke nodig
	{
	
	// unsigned int* const inddiag=new(nothrow) unsigned int[9];
       // if (inddiag==0) throw("Out of memory."); 

		// for (unsigned int iscolliOnDiag=0; iscolliOnDiag<ms*ns; iscolliOnDiag++)
		// {
			// if (scolli[iscolliOnDiag]==uniquescolli[iuniquescolli] && scolliOnDiag[iscolliOnDiag]==true)
			// {
				
				// mexPrintf("inddiag ... \n");
				// // mexPrintf("iscolliOnDiag: %d \n",iscolliOnDiag);
				// // mexPrintf("3*scompi[iscolliOnDiag]+scompj[iscolliOnDiag]: %d \n",3*scompi[iscolliOnDiag]+scompj[iscolliOnDiag]);
				// inddiag[3*scompi[iscolliOnDiag]+scompj[iscolliOnDiag]] = iscolliOnDiag;
			// }
		// }

		// for (unsigned int i=0; i<9; i++)	{ mexPrintf("inddiag [%d]: %f \n",i,inddiag[i]);}
		
	   
	// mexPrintf("uniquescolli[iuniquescolli]: %d \n",uniquescolli[iuniquescolli]);
    
	if (RegularColl[uniquescolli[iuniquescolli]]==1)
    {
	// mexPrintf("Is Regular: %d \n",uniquescolli[iuniquescolli]);
	
		for (unsigned int iXi=0; iXi<nXi; iXi++)
		{
	   	  			
        const double Xdiff=xiCart[3*iXi+0]-Coll[2*nColl+uniquescolli[iuniquescolli]];
        const double Ydiff=xiCart[3*iXi+1]-Coll[3*nColl+uniquescolli[iuniquescolli]];
        const double Zdiff=xiCart[3*iXi+2]-Coll[4*nColl+uniquescolli[iuniquescolli]];

        const double xiR=sqrt(Xdiff*Xdiff + Ydiff*Ydiff);
        const double xiTheta=atan2(Ydiff,Xdiff);
        const double xiZ=Zdiff;

        // EVALUATE GREEN'S FUNCTION
        greeneval3d(greenPtr,nGrSet,ugCmplx,tgCmplx,tg0Cmplx,xiR,xiZ,r1,r2,z1,z2,zs1,
                    interpr,interpz,extrapFlag,UmatOut,TmatOut,Coll,nColl,uniquescolli[iuniquescolli],4,UgrRe,
                    UgrIm,TgrRe,TgrIm,Tgr0Re,Tgr0Im);
        greenrotate3d(normal,iXi,xiTheta,nGrSet,ugCmplx,
                      tgCmplx,tg0Cmplx,UgrRe,UgrIm,TgrRe,TgrIm,
                      Tgr0Re,Tgr0Im,UXiRe,UXiIm,TXiRe,TXiIm,TXi0Re,
                      TXi0Im,UmatOut,TmatOut);
					  
	// URe[0]=10; 
        // SUM UP RESULTS, FOR ALL COLLOCATION POINTS        
		for (unsigned int iEltColl=0; iEltColl<nEltColl; iEltColl++)
        {
		 	// mexPrintf("ConsiderColl %d: %s\n",EltCollIndex[iEltColl],InListuniquecollj[EltCollIndex[iEltColl]] ? "true" : "false"); // DEBUG
		 	
			 	
		 	if (InListuniquecollj[EltCollIndex[iEltColl]])
		 	{


		 	
			double sumutil=H[iXi]*M[nEltColl*iXi+iEltColl]*Jac[iXi];
			
				for (unsigned int iuniquescolliind=0; iuniquescolliind<nuniquescolli[iuniquescolli]; iuniquescolliind++)
				{
					
					if (scollj[uniquescolliind[nuniquescollicumul+iuniquescolliind]] == EltCollIndex[iEltColl])
					{
					

					for (unsigned int iGrSet=0; iGrSet<nGrSet; iGrSet++)
					{
					
					
					const unsigned int ind0=ms*ns*iGrSet;
					// mexPrintf("ondiag : %s\n",ondiag ? "true" : "false"); // DEBUG
					if (ondiag==false)
					{
					// URe[ind0+uniquescolliind[nuniquescollicumul+iuniquescolliind]]+=sumutil*UXiRe[9*iGrSet+3*scompi[uniquescolliind[nuniquescollicumul+iuniquescolliind]]+scompj[uniquescolliind[nuniquescollicumul+iuniquescolliind]]]; 
					// unsigned int index=ind0+uniquescolliind[nuniquescollicumul+iuniquescolliind];
					// double value=sumutil*UXiRe[9*iGrSet+3*scompi[uniquescolliind[nuniquescollicumul+iuniquescolliind]]+scompj[uniquescolliind[nuniquescollicumul+iuniquescolliind]]]; 
					// time_t  assign_time= clock(); 
					// URe[index]+=value; 
					// float time1 = (float) (clock() - assign_time) / CLOCKS_PER_SEC; 
					// mexPrintf("time for assign was %f seconds\n", time1);	
					
					if (ugCmplx)
					{
						// UIm[ind0+uniquescolliind[nuniquescollicumul+iuniquescolliind]]+=sumutil*UXiIm[9*iGrSet+3*scompi[uniquescolliind[nuniquescollicumul+iuniquescolliind]]+scompj[uniquescolliind[nuniquescollicumul+iuniquescolliind]]];	
					// ind0+uniquescolliind[nuniquescollicumul+iuniquescolliind];
					// sumutil*UXiIm[9*iGrSet+3*scompi[uniquescolliind[nuniquescollicumul+iuniquescolliind]]+scompj[uniquescolliind[nuniquescollicumul+iuniquescolliind]]];	
					}
					
					}
					/*
					// mexPrintf("On diagonal:  scolli: %d \t scollj: %d \n",uniquescolli[iuniquescolli],scollj[uniquescolliind[nuniquescollicumul+iuniquescolliind]]);
					
					
						if (TmatOut)
						{
							// if (scollj[uniquescolliind[nuniquescollicumul+iuniquescolliind]]==0)
							{
								// mexPrintf("scollj[uniquescolliind[nuniquescollicumul+iuniquescolliind]]: %d \n",scollj[uniquescolliind[nuniquescollicumul+iuniquescolliind]]);
							}
							// if (scollj[uniquescolliind[nuniquescollicumul+iuniquescolliind]]==0)
							// if (uniquescolli[iuniquescolli]==scollj[uniquescolliind[nuniquescollicumul+iuniquescolliind]])
							// if (uniquescolli[iuniquescolli]<4)
							{
								// mexPrintf("On diagonal:  scolli: %d \t scollj: %d \n",uniquescolli[iuniquescolli],scollj[uniquescolliind[nuniquescollicumul+iuniquescolliind]]);
							}
							
							if (ondiag==false)
							{
							TRe[ind0+uniquescolliind[nuniquescollicumul+iuniquescolliind]]+=sumutil*TXiRe[9*iGrSet+3*scompi[uniquescolliind[nuniquescollicumul+iuniquescolliind]]+
																											 scompj[uniquescolliind[nuniquescollicumul+iuniquescolliind]]];
							
							
							if ((double) scollj[uniquescolliind[nuniquescollicumul+iuniquescolliind]] == Coll[nColl+iuniquescolli]-1)
							{
								// mexPrintf("reg...\n");
							}
							
							// TRe[ind0+nDof*(rowBeg+0)+rowBeg+0]-=sumutil*TXi0Re[9*iGrSet+0];  // txx
							
							
							// mexPrintf("TRe[ind0+uniquescolliind[nuniquescollicumul+iuniquescolliind]]: %e \n",TRe[ind0+uniquescolliind[nuniquescollicumul+iuniquescolliind]]);
							// mexPrintf("sumutil: %e \n",sumutil);
							// mexPrintf("TXiRe[9*iGrSet+3*scompi[uniquescolliind[nuniquescollicumul+iuniquescolliind]]+scompj[uniquescolliind[nuniquescollicumul+iuniquescolliind]]]: %e \n",  TXiRe[9*iGrSet+3*scompi[uniquescolliind[nuniquescollicumul+iuniquescolliind]]+scompj[uniquescolliind[nuniquescollicumul+iuniquescolliind]]]);
																											 
							if (tgCmplx)
							{
								TIm[ind0+uniquescolliind[nuniquescollicumul+iuniquescolliind]]+=sumutil*TXiIm[9*iGrSet+3*scompi[uniquescolliind[nuniquescollicumul+iuniquescolliind]]+
																											 scompj[uniquescolliind[nuniquescollicumul+iuniquescolliind]]];	
							}
							}
							
							
							// mexPrintf("inddiag[3*scompi[uniquescolliind[nuniquescollicumul+iuniquescolliind]]+scompj[uniquescolliind[nuniquescollicumul+iuniquescolliind]]] %f \n",inddiag[3*scompi[uniquescolliind[nuniquescollicumul+iuniquescolliind]]+scompj[uniquescolliind[nuniquescollicumul+iuniquescolliind]]]);
							
							if (ondiag==true && inddiag[9*iuniquescolli+3*scompi[uniquescolliind[nuniquescollicumul+iuniquescolliind]]+scompj[uniquescolliind[nuniquescollicumul+iuniquescolliind]]]!=-1)
							{
							// mexPrintf("Calculating ondiag... \n");
							TRe[ind0+inddiag[9*iuniquescolli+3*scompi[uniquescolliind[nuniquescollicumul+iuniquescolliind]]+scompj[uniquescolliind[nuniquescollicumul+iuniquescolliind]]]]-=sumutil*TXi0Re[9*iGrSet+3*scompi[uniquescolliind[nuniquescollicumul+iuniquescolliind]]+scompj[uniquescolliind[nuniquescollicumul+iuniquescolliind]]];
							
							// if (inddiag[9*iuniquescolli+3*scompi[uniquescolliind[nuniquescollicumul+iuniquescolliind]]+scompj[uniquescolliind[nuniquescollicumul+iuniquescolliind]]]==-1)
							// {
								// mexPrintf("scompi: %d \n",scompi[uniquescolliind[nuniquescollicumul+iuniquescolliind]]);
								// mexPrintf("scompj: %d \n",scompj[uniquescolliind[nuniquescollicumul+iuniquescolliind]]);
							// }
							
							// mexPrintf("inddiag[...]: %d \n",inddiag[9*iuniquescolli+3*scompi[uniquescolliind[nuniquescollicumul+iuniquescolliind]]+scompj[uniquescolliind[nuniquescollicumul+iuniquescolliind]]]);
			
							// TRe[ind0+uniquescolliind[nuniquescollicumul+iuniquescolliind]]-=sumutil*TXi0Re[9*iGrSet+3*scompi[uniquescolliind[nuniquescollicumul+iuniquescolliind]]+
																											 // scompj[uniquescolliind[nuniquescollicumul+iuniquescolliind]]];
																										 
							
							if (tgCmplx)
							{
							
								TIm[ind0+inddiag[9*iuniquescolli+3*scompi[uniquescolliind[nuniquescollicumul+iuniquescolliind]]+scompj[uniquescolliind[nuniquescollicumul+iuniquescolliind]]]]-=sumutil*TXi0Im[9*iGrSet+3*scompi[uniquescolliind[nuniquescollicumul+iuniquescolliind]]+scompj[uniquescolliind[nuniquescollicumul+iuniquescolliind]]];
								// TIm[ind0+uniquescolliind[nuniquescollicumul+iuniquescolliind]]-=sumutil*TXi0Im[9*iGrSet+3*scompi[uniquescolliind[nuniquescollicumul+iuniquescolliind]]+
																											 // scompj[uniquescolliind[nuniquescollicumul+iuniquescolliind]]];		
							}
							}
						
						}
						
					*/	
					
					}
					
					
					}
					
					
					
				}
			}
          }
          
        }
      }
  nuniquescollicumul+=nuniquescolli[iuniquescolli];
  // delete [] inddiag;  
  }
  } 
  
   	 
   //  s empty: all collocation points
  else  
  {
  
  for (unsigned int iColl=0; iColl<nColl; iColl++) // bepalen welke nodig
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
		  // mexPrintf("sumutil %f: \n",sumutil); // DEBUG
          unsigned int rowBeg=3*iColl;  // welke rijpositie -> sColi
          unsigned int colBeg=3*EltCollIndex[iEltColl]; // welke colompositie 3*sColj
          for (unsigned int iGrSet=0; iGrSet<nGrSet; iGrSet++)
          {
            const unsigned int ind0 =nDof*nDof*iGrSet;
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
              // Account for singular part of Green's function on the
              // diagonal terms.
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
  
 }

  delete [] EltNod;
  delete [] xi;
  delete [] H;
  delete [] N;
  delete [] M;
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
