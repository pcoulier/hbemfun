#include <string.h>
#include "eltdef.h"
#include "bemcollpoints.h"
#include "bemintreg3d.h"
#include "bemintreg3dnodiag.h"
#include "bemintreg3ddiag.h"
#include "bemintreg3dperiodic.h"
#include "bemintreg2d.h"
#include "bemintregaxi.h"
#include "bemintsing3d.h"
#include "bemintsing3dperiodic.h"
#include "bemintsing2d.h"
#include "bemintsingaxi.h"
#include "bemdimension.h"
#include "bemisaxisym.h"
#include "s2coll.h"
#include "checklicense.h"
#include <math.h>
#include <time.h>
#include <new>
#include "mex.h"

#ifndef __GNUC__
#define strcasecmp _strcmpi
#endif

#ifndef _int64_
typedef long long int int64;
typedef unsigned long long int uint64;
#endif

using namespace std;

//==============================================================================
void bemmat(const bool& probAxi, const bool& probPeriodic, const unsigned int& probDim,
            const unsigned int& nColDof, const bool& UmatOut, const bool& TmatOut,
            const double* const Nod, const unsigned int& nNod,
            const double* const Elt, const unsigned int& nElt,
            const unsigned int* const TypeID,
            const char* const TypeName[], const char* const TypeKeyOpts[],
            const unsigned int* const nKeyOpt,
            const unsigned int& nEltType, const double* const CollPoints,
            const unsigned int& nTotalColl,
            const void* const* const greenPtr, const unsigned int& nGrSet,
            const unsigned int& nugComp, const bool& ugCmplx,
            const bool& tgCmplx, const bool& tg0Cmplx,
            double* const URe, double* const UIm,
            double* const TRe, double* const TIm,
            const double* const s, const unsigned int& ms, const unsigned int& ns,
            const double L, const double* const ky, const unsigned int nWave, const unsigned int nmax,
			const unsigned int* const EltParent, const unsigned int* const nEltNod, const unsigned int* const nEltColl,
			const unsigned int* const EltShapeN, const unsigned int* const EltShapeM, const unsigned int* const EltDim,
			const unsigned int* const AxiSym, const unsigned int* const Periodic, const unsigned int* const nGauss,
			const unsigned int* const nEltDiv, const unsigned int* const nGaussSing, const unsigned int* const nEltDivSing,
			const unsigned int* const ncumulEltCollIndex, const unsigned int* const eltCollIndex,
			const unsigned int* const ncumulSingularColl, const unsigned int* const nSingularColl, const int& NSingularColl, 
			const unsigned int* const RegularColl, 
			const unsigned int* const ncumulEltNod, const double* const EltNod,
			const unsigned int* const RefEltType,  const unsigned int* const ncumulnXi, const unsigned int* const nXi, const double* const xi, const double* const H,
			const unsigned int* const ncumulNshape, const double* const Nshape, const double* const Mshape, const double* const dNshape)
//==============================================================================
{



//int Test=42949672960;
//mexPrintf("Test: %d \n",Test);
//
//int64 Test64=42949672960;
//mexPrintf("Test64: %lld \n",Test64);
//
//uint64 Testu64=42949672960;
//mexPrintf("Testu64: %lld \n",Testu64);




//mexPrintf("nTotalColl: %d \n",nTotalColl);
//mexPrintf("nEltColl: %d \n",nEltColl);
  time_t  start_bemmat = clock();  
  // time_t  start_prepro_bemmat = clock();  
  
  	
  // time_t  start_prepro1_bemmat = clock();  
  
  // DEGREES OF FREEDOM
  unsigned int nDof=nColDof*nTotalColl;

  // INITIALIZE VARIABLES
  
  /*
  int* const RegularColl=new(nothrow) int[2*nTotalColl];
  if (RegularColl==0) throw("Out of memory.");
  int nRegularColl;
  int nSingularColl;
  */
  
  double* const xiSing=new(nothrow) double[2];
  if (xiSing==0) throw("Out of memory.");

   // Hier: welke elementen dragen bij tot de matrix? 
   bool spassed;
   
   if (s!=0)
   {spassed=true;}
   else
   {spassed=false;}
   
   unsigned int* const scolli=new(nothrow) unsigned int[ms*ns];
       if (scolli==0) throw("Out of memory scolli.");
   unsigned int* const scompi=new(nothrow) unsigned int[ms*ns];
       if (scompi==0) throw("Out of memory scompi.");
   // int* const uniquescolli=new(nothrow) int[ms*ns];
   unsigned int* const uniquescolli=new(nothrow) unsigned int[nDof];
       if (uniquescolli==0) throw("Out of memory uniquescolli.");
   unsigned int* const Nuniquescolli=new(nothrow) unsigned int[1];
       if (Nuniquescolli==0) throw("Out of memory Nuniquescolli.");       
   // int* const nuniquescolli=new(nothrow) int[ms*ns];
   unsigned int* const nuniquescolli=new(nothrow) unsigned int[nDof];
       if (nuniquescolli==0) throw("Out of memory nuniquescolli.");       
   unsigned int* const uniquescolliind=new(nothrow) unsigned int[ms*ns];
       if (uniquescolliind==0) throw("Out of memory uniquescolliind.");     
       
   unsigned int* const scollj=new(nothrow) unsigned int[ms*ns];
       if (scollj==0) throw("Out of memory scollj.");
   unsigned int* const scompj=new(nothrow) unsigned int[ms*ns];
       if (scompj==0) throw("Out of memory scompj.");
   // int* const uniquescollj=new(nothrow) int[ms*ns];
   unsigned int* const uniquescollj=new(nothrow) unsigned int[nDof];
       if (uniquescollj==0) throw("Out of memory uniquescollj.");
   unsigned int* const Nuniquescollj=new(nothrow) unsigned int[1];
       if (Nuniquescollj==0) throw("Out of memory Nuniquescollj.");       
   // int* const nuniquescollj=new(nothrow) int[ms*ns];
   unsigned int* const nuniquescollj=new(nothrow) unsigned int[nDof];
       if (nuniquescollj==0) throw("Out of memory nuniquescollj.");       
   unsigned int* const uniquescolljind=new(nothrow) unsigned int[ms*ns];
       if (uniquescolljind==0) throw("Out of memory uniquescolljind.");    
     
   // List of all collocation points
   bool* const InListuniquecollj=new(nothrow) bool[nTotalColl];
       if (InListuniquecollj==0) throw("Out of memory InListuniquecollj."); 
	   
	int* const DeltaInListuniquecollj=new(nothrow) int[nTotalColl];
       if (DeltaInListuniquecollj==0) throw("Out of memory DeltaInListuniquecollj."); 	   
	   
	bool* const scolliOnDiag=new(nothrow) bool[ms*ns];
       if (scolliOnDiag==0) throw("Out of memory scolliOnDiag."); 
       
    
	bool* const scolliOnDiaginddiag=new(nothrow) bool[ms*ns];
       if (scolliOnDiaginddiag==0) throw("Out of memory scolliOnDiaginddiag."); 
    
    Nuniquescolli[0]=0;
	Nuniquescollj[0]=0;
	   
   //  Get collocation points of interest
   if (s !=0)
   {	
		// time_t  start_s2coll = clock();  
		// float time_s2coll = (float) (clock() - start_s2coll) / CLOCKS_PER_SEC; 
 	    s2coll(s,ms,ns,nDof,probDim,scolli,scompi,uniquescolli,Nuniquescolli,nuniquescolli,uniquescolliind,
                                     scollj,scompj,uniquescollj,Nuniquescollj,nuniquescollj,uniquescolljind,
									 scolliOnDiag);
		// mexPrintf("time for s2coll was %f seconds\n", time_s2coll);
		// time_t  start_InList = clock();
		
         for (unsigned int iInListuniquecollj=0; iInListuniquecollj<nTotalColl; iInListuniquecollj++)   
         {
             InListuniquecollj[iInListuniquecollj]=0;
         }
         for (unsigned int iuniquescollj=0; iuniquescollj<Nuniquescollj[0]; iuniquescollj++)
         {
             InListuniquecollj[uniquescollj[iuniquescollj]]=1;
         }


		DeltaInListuniquecollj[0]=0; 
		for (unsigned int iDeltaInListuniquecollj=1; iDeltaInListuniquecollj<nTotalColl; iDeltaInListuniquecollj++)   
         {
             if (InListuniquecollj[iDeltaInListuniquecollj-1]==0)
			 {
				DeltaInListuniquecollj[iDeltaInListuniquecollj]=DeltaInListuniquecollj[iDeltaInListuniquecollj-1]+1;
			 }
			 else
			 {
				DeltaInListuniquecollj[iDeltaInListuniquecollj]=DeltaInListuniquecollj[iDeltaInListuniquecollj-1];
			 }
         }
		 
		 // for (int i=0; i<nTotalColl; i++)
		// {
			// mexPrintf("DeltaInListuniquecollj[%d]: %d \n",i,DeltaInListuniquecollj[i]);
		// }
	
		 
		 
	 	// float time_InList = (float) (clock() - start_InList) / CLOCKS_PER_SEC; 
		// mexPrintf("time for _InList was %f seconds\n", time_InList);
	
   }
   else // All collocation points
   {
       for (unsigned int iInListuniquecollj=0; iInListuniquecollj<nTotalColl; iInListuniquecollj++)   
         {
             InListuniquecollj[iInListuniquecollj]=1;
         }   
   }
	
for (unsigned int i=0; i<ms*ns; i++)
{
	scolliOnDiaginddiag[i]=scolliOnDiag[i];
}
	
	// for (int i=0; i<ms*ns; i++)
	// {
		// mexPrintf("scolli[%d]: %d \n",i,scolli[i]);
	// }
	
	// for (int i=0; i<ms*ns; i++)
	// {
		// mexPrintf("scollj[%d]: %d \n",i,scollj[i]);
	// }
	
	// for (int i=0; i<ms*ns; i++)
	// {
			 // mexPrintf("scolliOnDiag [%d]: %s\n",i,scolliOnDiag [i] ? "true" : "false"); // DEBUG
	// }
	
	// for (int i=0; i<ms*ns; i++)
	// {
			 // mexPrintf("uniquescolliind [%d]: %d \n",i,uniquescolliind[i]); // DEBUG
	// }
	
	// for (int i=0; i<Nuniquescolli[0]; i++)
	// {
			 // mexPrintf("uniquescolli [%d]: %d \n",i,uniquescolli[i]); // DEBUG
	// }


	// for (int i=0; i<Nuniquescollj[0]; i++)
	// {
			 // mexPrintf("uniquescollj [%d]: %d \n",i,uniquescollj[i]); // DEBUG
	// }
	
	bool ondiag=false;
	unsigned int NOnDiag=0;
	
	// while (ondiag ==false && iOnDiag<ms*ns)
	for (unsigned int iOnDiag=0; iOnDiag<ms*ns;iOnDiag++)
	{
		if (scolliOnDiag[iOnDiag]==true)
		{
			NOnDiag++;
			ondiag=true;
			// mexPrintf("ondiag : %s\n",ondiag ? "true" : "false"); // DEBUG
		}
		// else
		// {
			// iOnDiag++;
		// }
	}
	
	 // mexPrintf("ondiag : %s\n",ondiag ? "true" : "false"); // DEBUG
	

	
	// mexPrintf("NOnDiag : %d \n",NOnDiag); // DEBUG
	// mexPrintf("NOnDiagUnique %d \n",NOnDiagUnique);
	
	// for (int i=0; i<ms*ns; i++)
	// {
			 // mexPrintf("scolliOnDiag [%d]: %s\n",i,scolliOnDiag [i] ? "true" : "false"); // DEBUG
	// }
	
	// for (int i=0; i<ms*ns; i++)
	// {
		 // mexPrintf("scolliOnDiag[%d]: %s\n",i,scolliOnDiag [i] ? "true" : "false"); // DEBUG
	// }
	// for (int i=0; i<NOnDiag; i++)
	// {
		 // mexPrintf("scolliOnDiagUnique[%d]: %f \n",i,scolliOnDiagUnique[i]); // DEBUG
	// }
	// time_t  start_sdiag = clock();  
	// float time_sdiag = (float) (clock() - start_sdiag) / CLOCKS_PER_SEC; 
 	    
		
	
	
		
		
	// mexPrintf("time for sdiag was %f seconds\n", time_sdiag);	
		
	// for (int i=0;i<ms*ns;i++)
	// {
		// mexPrintf("s[%d] : %f \n",i,s[i]); // DEBUG
	// }
	
	// for (int i=0;i<NOnDiagUnique*nDof;i++)
	// {
		// mexPrintf("sdiag[%d] : %f \n",i,sdiag[i]); // DEBUG
	// }

	
	// float time_prepro1_bemmat = (float) (clock() - start_prepro1_bemmat) / CLOCKS_PER_SEC; 
	// mexPrintf("time for prepro1_bemmat was %f seconds\n", time_prepro1_bemmat);
	
	// time_t  start_prepro2_bemmat = clock();  
	

	
	// float time_InListDiag = (float) (clock() - start_InListDiag) / CLOCKS_PER_SEC; 
	// mexPrintf("time for _InListDiag was %f seconds\n", time_InListDiag);
	
	
	
	// for (int i=0; i<NOnDiag*nDof; i++)
	// {
		// mexPrintf("sdiagcolli[%d]: %d \n",i,sdiagcolli[i]);
	// }
	
	// for (int i=0; i<NOnDiag*nDof; i++)
	// {
		// mexPrintf("sdiagcollj[%d]: %d \n",i,sdiagcollj[i]);
	// }
	
	// for (int i=0; i<NOnDiag*nDof; i++)
	// {
		// mexPrintf("sdiagcompj[%d]: %d \n",i,sdiagcompj[i]);
	// }
	
	
	// mexPrintf("Nuniquesdiagcolli[0]: %d \n",Nuniquesdiagcolli[0]);
	// mexPrintf("Nuniquesdiagcollj[0]: %d \n",Nuniquesdiagcollj[0]);
	
	// for (int i=0; i<Nuniquesdiagcolli[0]; i++)
	// {
		// mexPrintf("uniquesdiagcolli[%d]: %d \n",i,uniquesdiagcolli[i]);
	// }
	
	// for (int i=0; i<NOnDiagUnique*nDof; i++)
	// {
		// mexPrintf("uniquesdiagcolliind[%d]: %d \n",i,uniquesdiagcolliind[i]);
	// }
	
	// for (int i=0; i<NOnDiagUnique*nDof; i++)
	// {
		// mexPrintf("uniquesdiagcolljind[%d]: %d \n",i,uniquesdiagcolljind[i]);
	// }
	
	// for (int i=0; i<Nuniquesdiagcollj[0]; i++)
	// {
		// mexPrintf("uniquesdiagcollj[%d]: %d \n",i,uniquesdiagcollj[i]);
	// }
	
	
	// mexPrintf("Nuniquescolli[0]: %d \n",Nuniquescolli[0]);
	

	
	
	
	
	
	
	int* const inddiag=new(nothrow) int[probDim*probDim*Nuniquescolli[0]];
		if (inddiag==0) throw("Out of memory inddiag."); 

	for (unsigned int i=0; i<probDim*probDim*Nuniquescolli[0]; i++)
	{
		inddiag[i]=-1;
	}
	
	
	
		
  	for (unsigned int iuniquescolli=0; iuniquescolli<Nuniquescolli[0]; iuniquescolli++) // bepalen welke nodig
	{
		for (unsigned int iscolliOnDiag=0; iscolliOnDiag<ms*ns; iscolliOnDiag++)
		{
			// if (scolli[iscolliOnDiag]==uniquescolli[iuniquescolli] && scolliOnDiag[iscolliOnDiag]==true)
			if (scolli[iscolliOnDiag]==uniquescolli[iuniquescolli] && scolliOnDiaginddiag[iscolliOnDiag]==true)
			{
				inddiag[probDim*probDim*iuniquescolli+3*scompi[iscolliOnDiag]+scompj[iscolliOnDiag]] = iscolliOnDiag;
			}
		}
	}

   bool* const blockdiag=new(nothrow) bool[Nuniquescolli[0]];
      if (blockdiag==0) throw("Out of memory blockdiag."); 

	
	for (unsigned int iuniquescolli=0; iuniquescolli<Nuniquescolli[0];iuniquescolli++)
	{
		blockdiag[iuniquescolli]=true;

		unsigned int iblockdiag=0;
		while (iblockdiag<probDim*probDim && blockdiag[iuniquescolli]==true)
		{
			if (inddiag[probDim*probDim*iuniquescolli+iblockdiag]==-1)
			{blockdiag[iuniquescolli]=false;}
			iblockdiag++;
		}
	}
	
	// for (int i=0; i<Nuniquescolli[0]; i++)
	// {
		// mexPrintf("blockdiag[%d]: %s \n",i, blockdiag[i] ? "true": "false");
	// }
	
	bool* const blocks=new(nothrow) bool[Nuniquescolli[0]];
      if (blocks==0) throw("Out of memory blocks."); 
	
	
	unsigned int nuniquescollicumulver=0;
	for (unsigned int iuniquescolli=0; iuniquescolli<Nuniquescolli[0];iuniquescolli++)
	{
		bool comp0=false;
		bool comp1=false;
		bool comp2=false;
		
		blocks[iuniquescolli]=false;
		
		unsigned int iuniquescolliind=0;
		while (iuniquescolliind<nuniquescolli[iuniquescolli] && blocks[iuniquescolli]==false && (comp0==false || comp1==false || comp2==false))
		{
		if (scompi[uniquescolliind[nuniquescollicumulver+iuniquescolliind]]==0){comp0=true;}
		if (scompi[uniquescolliind[nuniquescollicumulver+iuniquescolliind]]==1){comp1=true;}
		if (scompi[uniquescolliind[nuniquescollicumulver+iuniquescolliind]]==2){comp2=true;}
		
		if (comp0==true && comp1 == true && comp2==true)
		{
		blocks[iuniquescolli]=true;
		}
		
		iuniquescolliind++;
		}
		
		// // mexPrintf("blocks[iuniquescolli]: %s \n", blocks[iuniquescolli]? "true": "false");
		
		// mexPrintf("comp0 %s \n",comp0 ? "true": "false");
		// mexPrintf("comp1 %s \n",comp1 ? "true": "false");
		// mexPrintf("comp2 %s \n",comp2 ? "true": "false");
		
		if ((nuniquescolli[iuniquescolli] % (probDim*probDim)) !=0){blocks[iuniquescolli]=false;}
		// mexPrintf("mod: %d \n",nuniquescolli[iuniquescolli] % (probDim*probDim));
		
		// mexPrintf("blocks[iuniquescolli]: %s \n", blocks[iuniquescolli]? "true": "false");
		
	nuniquescollicumulver+=nuniquescolli[iuniquescolli];	
	}
		
	
	// for (int i=0; i<Nuniquescolli[0]; i++)
	// {
		// mexPrintf("blocks[%d]: %s \n",i, blocks[i] ? "true": "false");
	// }
	
	// // // // // // // // // // for (int iuniquescolli=0; iuniquescolli<Nuniquescolli[0]; iuniquescolli++) // bepalen welke nodig
	// // // // // // // // // // {
		// // // // // // // // // // for (int iscolliOnDiag=0; iscolliOnDiag<ms*ns; iscolliOnDiag++)
		// // // // // // // // // // {
			// // // // // // // // // // // if (scolli[iscolliOnDiag]==uniquescolli[iuniquescolli] && scolliOnDiag[iscolliOnDiag]==true)
			// // // // // // // // // // if (scolli[iscolliOnDiag]==uniquescolli[iuniquescolli] && scolliOnDiaginddiag[iscolliOnDiag]==true)
			// // // // // // // // // // {
				// // // // // // // // // // inddiag[probDim*probDim*iuniquescolli+scompi[iscolliOnDiag]+3*scompj[iscolliOnDiag]] = iscolliOnDiag;
			// // // // // // // // // // }
		// // // // // // // // // // }
	// // // // // // // // // // }
	
	
	
	
	
	


	// for (int i=0; i<probDim*probDim*Nuniquescolli[0]; i++)
	// {
		// mexPrintf("inddiag[%d]: %d \n",i,inddiag[i]);
	// }


	// for (int is=0; is<ms*ns; is++)
	// {
		// if (scolli[is]==scollj[is])
		// {
			// mexPrintf("On diagonal version 0...\n");
		// }
	// }
	
	
	
	
	
	

	// bool ondiag=false;
	
	
	// bool* const ListUniqueColliOnDiag=new(nothrow) bool[Nuniquescolli];
       // if (ListUniqueColliOnDiag==0) throw("Out of memory."); 
   
	// int* const uniquescolliindondiag=new(nothrow) int[ms*ns];
       // if (uniquescolliindondiag==0) throw("Out of memory.");
		 
	
	
	
	// for (int iListUniqueColliOnDiag=0; iListUniqueColliOnDiag<Nuniquescolli; iListUniqueColliOnDiag++)   
    // {
        // ListUniqueColliOnDiag[iListUniqueColliOnDiag]=0;
    // }
			
	
	
	// int nuniquescollicumulver=0;	
	// for (int iuniquescolli=0; iuniquescolli<Nuniquescolli[0]; iuniquescolli++) // bepalen welke nodig
	// {
		// for (int iuniquescolliind=0; iuniquescolliind<nuniquescolli[iuniquescolli]; iuniquescolliind++)
		// {		
			// if (uniquescolli[iuniquescolli] == scollj[uniquescolliind[nuniquescollicumulver+iuniquescolliind]])
			// {
				// ListUniqueColliOnDiag[iuniquescolli]=1;
				
				// mexPrintf("On diagonal version 1...\n");
			// }
		// }
		// nuniquescollicumulver+=nuniquescolli[iuniquescolli];
	// }
	
	
	
	
	 // time_t  start_time_loop = clock();  
	
	
		// mexPrintf("ondiag : %s\n",ondiag ? "true" : "false"); // DEBUG
	

	unsigned int EltType_prev = int(Elt[nElt+0]);
	unsigned int nXi_loc=nXi[EltType_prev-1];
	double* const xi_loc=new(nothrow) double[2*nXi_loc];
	if (xi_loc==0) throw("Out of memory.");
	double* const H_loc=new(nothrow) double[nXi_loc];
	if (H_loc==0) throw("Out of memory.");
	
	double* const N_loc=new(nothrow) double[nXi_loc*nEltNod[RefEltType[EltType_prev-1]]];
	if (N_loc==0) throw("Out of memory.");
		
	double* const M_loc=new(nothrow) double[nXi_loc*nEltNod[RefEltType[EltType_prev-1]]];
	if (M_loc==0) throw("Out of memory.");
		
	double* const dN_loc=new(nothrow) double[2*nXi_loc*nEltNod[RefEltType[EltType_prev-1]]];
	if (dN_loc==0) throw("Out of memory.");
	
	
	for(unsigned int iXi=0; iXi<2*nXi_loc; iXi++)
	{
		xi_loc[iXi] = xi[2*ncumulnXi[EltType_prev-1]+iXi];
	}
	for(unsigned int iH=0; iH<nXi_loc; iH++)
	{
		H_loc[iH] = H[ncumulnXi[EltType_prev-1]+iH];
	}	
		
		
	for(unsigned int iN=0; iN<nXi_loc*nEltNod[RefEltType[EltType_prev-1]]; iN++)
	{
		N_loc[iN] = Nshape[ncumulNshape[EltType_prev-1]+iN];
		M_loc[iN] = Mshape[ncumulNshape[EltType_prev-1]+iN];
	}
	for(unsigned int idN=0; idN<2*nXi_loc*nEltNod[RefEltType[EltType_prev-1]]; idN++)
	{
		dN_loc[idN] = dNshape[2*ncumulNshape[EltType_prev-1]+idN];
	}

	
		
   unsigned int NEltCollConsider=0;
   
	// float time_prepro2_bemmat = (float) (clock() - start_prepro2_bemmat) / CLOCKS_PER_SEC; 
	// mexPrintf("time for prepro2_bemmat was %f seconds\n", time_prepro2_bemmat);
   
	// float time_prepro_bemmat = (float) (clock() - start_prepro_bemmat) / CLOCKS_PER_SEC; 
	// mexPrintf("time for prepro_bemmat was %f seconds\n", time_prepro_bemmat);
		
	// time_t  start_bemmat_loop = clock();  
	
	// float timeTest=0.0;
	// float timeTest_sing=0.0;
		
		
		// /*
  // ELEMENT-BY-ELEMENT INTEGRATION 
  for (unsigned int iElt=0; iElt<nElt; iElt++)  // Enkel loop over nodige elementen
  {
	//mexPrintf("UmatOut: %s \n", UmatOut ? "true": "false");
	
	
    unsigned int EltType = int(Elt[nElt+iElt]);
	/*
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
    int nGauss;
    int nEltDivSing;
	
	eltdef(EltType,TypeID,TypeName,TypeKeyOpts,nKeyOpt,nEltType,EltParent,nEltNod,
         nEltColl,EltShapeN,EltShapeM,EltDim,AxiSym,Periodic,nGauss,nEltDiv,
         nGaussSing,nEltDivSing);
		 
	
    int* const eltCollIndex=new(nothrow) int[nEltColl];
    if (eltCollIndex==0) throw("Out of memory.");
    
	
	BemEltCollIndex(Elt,iElt,nElt,CollPoints,nCentroidColl,nTotalColl,
                    nEltColl,nEltNod,eltCollIndex);
	BemRegularColl(Elt,iElt,nElt,Nod,nNod,CoincNodes,SlavesExist,
                   CollPoints,nCentroidColl,nTotalColl,RegularColl,
                   nRegularColl,nSingularColl,TypeID,nKeyOpt,TypeName,
                   TypeKeyOpts,nEltType);
				   
				   
	*/
 
    unsigned int iEltCollIndex=0;
    bool ConsiderElement=0;

    // Check if element iElt should be considered or not
	// time_t  start_eltconsider = clock();
    // // // while(iEltCollIndex<nEltColl && ConsiderElement==0)  
    // // // {
       // // // if (InListuniquecollj[eltCollIndex[iEltCollIndex]])
       // // // {
          // // // ConsiderElement=1;
       // // // }
       // // // iEltCollIndex++;
    // // // }         
	
	while(iEltCollIndex<nEltColl[iElt] && ConsiderElement==0)  
    {
       if (InListuniquecollj[eltCollIndex[ncumulEltCollIndex[iElt]+iEltCollIndex]])
       {
          ConsiderElement=1;
       }
       iEltCollIndex++;
    }       
	
	
	if(ConsiderElement==0)     
    {
		NEltCollConsider+=nEltColl[iElt];
	}
	// float time_eltconsider = (float) (clock() - start_eltconsider) / CLOCKS_PER_SEC; 
	// if (iElt==0)	 
	// {mexPrintf("time for nElt*eltconsider was %f seconds\n", nElt*time_eltconsider);}
	
	
	
    
    if(ConsiderElement==1)     // Element iElt should be considered
    {
	
	//
	unsigned int* const eltCollIndex_loc=new(nothrow) unsigned int[nEltColl[iElt]];
	if (eltCollIndex_loc==0) throw("Out of memory.");
	
	for(unsigned int iEltCollIndex=0; iEltCollIndex<nEltColl[iElt]; iEltCollIndex++)
	{
			eltCollIndex_loc[iEltCollIndex] = eltCollIndex[ncumulEltCollIndex[iElt]+iEltCollIndex];
					
	}
	
	//
	unsigned int* const RegularColl_loc=new(nothrow) unsigned int[2*nTotalColl];
		if (RegularColl_loc==0) throw("Out of memory.");

		
	// mexPrintf("iElt: %u \n",iElt);			
	// // // // // for(unsigned int i=0; i<2*nTotalColl; i++)
	// // // // // {
		   // // // // // RegularColl_loc[i]=RegularColl[(uint64)(2*nTotalColl*iElt+i)];
		   // // // // // // mexPrintf("RegularColl_loc[i]: %u \n",RegularColl_loc[i]);
	// // // // // }
	
	for(unsigned int i=0; i<nTotalColl; i++)
	{
		   RegularColl_loc[i]=1;
	}

	for(unsigned int iSingular=0; iSingular < nSingularColl[iElt]; iSingular++)
	{
		RegularColl_loc[RegularColl[(uint64)(ncumulSingularColl[iElt]+iSingular)]]=0;
		RegularColl_loc[nTotalColl+RegularColl[(uint64)(ncumulSingularColl[iElt]+iSingular)]]=RegularColl[(uint64)(NSingularColl+ncumulSingularColl[iElt]+iSingular)];
	}
	
	// for(unsigned int i=0; i<2*nTotalColl; i++)
	// {
			   // mexPrintf("RegularColl_loc[%u]: %u \n",i,RegularColl_loc[i]);
	// }
		
		
	//
	double* const EltNod_loc=new(nothrow) double[3*nEltNod[iElt]];
	if (EltNod_loc==0) throw("Out of memory.");


	for(unsigned int i=0; i<3*nEltNod[iElt]; i++)
	{
		EltNod_loc[i] = EltNod[3*ncumulEltNod[iElt]+i];
	}
	
	//
	
	if (EltType != EltType_prev)
	{
	nXi_loc=nXi[EltType-1];
	
	delete [] xi_loc;
	double* const xi_loc=new(nothrow) double[2*nXi_loc];
	if (xi_loc==0) throw("Out of memory.");
	
	delete [] H_loc;	
	double* const H_loc=new(nothrow) double[nXi_loc];
	if (H_loc==0) throw("Out of memory.");
	
	delete [] N_loc;
	double* const N_loc=new(nothrow) double[nXi_loc*nEltNod[RefEltType[EltType-1]]];
	if (N_loc==0) throw("Out of memory.");
		
	delete [] M_loc;
	double* const M_loc=new(nothrow) double[nXi_loc*nEltNod[RefEltType[EltType-1]]];
	if (M_loc==0) throw("Out of memory.");
		
	delete [] dN_loc;
	double* const dN_loc=new(nothrow) double[2*nXi_loc*nEltNod[RefEltType[EltType-1]]];
	if (dN_loc==0) throw("Out of memory.");
	

	
		
	for(unsigned int iXi=0; iXi<2*nXi_loc; iXi++)
	{
		xi_loc[iXi] = xi[2*ncumulnXi[EltType-1]+iXi];
	}
	for(unsigned int iH=0; iH<nXi_loc; iH++)
	{
		H_loc[iH] = H[ncumulnXi[EltType-1]+iH];
	}
	

		
	for(unsigned int iN=0; iN<nXi_loc*nEltNod[RefEltType[EltType-1]]; iN++)
	{
		N_loc[iN] = Nshape[ncumulNshape[EltType-1]+iN];
		M_loc[iN] = Mshape[ncumulNshape[EltType-1]+iN];
	}
	for(unsigned int idN=0; idN<2*nXi_loc*nEltNod[RefEltType[EltType-1]]; idN++)
	{
		dN_loc[idN] = dNshape[2*ncumulNshape[EltType-1]+idN];
	}


		
	}	
		
	// for(int i=0; i<2*nXi[EltType]; i++) {
                     // mexPrintf("xi_loc[%d]: %f\n",i,xi_loc[i]); // DEBUG
	// }
	// for(int i=0; i<nXi[EltType]; i++) {
                     // mexPrintf("H_loc[%d]: %f\n",i,H_loc[i]); // DEBUG
	// }
	
		
	// for(int i=0; i<3*nEltNod[iElt]; i++) {
                     // mexPrintf("EltNod_loc[%d]: %f\n",i,EltNod_loc[i]); // DEBUG
	// }
		
	// // BemRegularColl(Elt,iElt,nElt,Nod,nNod,CoincNodes,SlavesExist,
                   // // CollPoints,nCentroidColl,nTotalColl,RegularColl,
                   // // nRegularColl,nSingularColl,TypeID,nKeyOpt,TypeName,
                   // // TypeKeyOpts,nEltType);
	
	
	 //mexPrintf("iElt %d \n",iElt);

	 // time_t  time_perelem = clock();  
	
    // REGULAR INTEGRATION
    if (probDim==3)
    {
      if (probPeriodic){
        bemintreg3dperiodic(Nod,nNod,Elt,iElt,nElt,TypeID,nKeyOpt,TypeName,TypeKeyOpts,
                            nEltType,CollPoints,nTotalColl,RegularColl_loc,eltCollIndex_loc,nDof,greenPtr,
                            nGrSet,ugCmplx,tgCmplx,tg0Cmplx,URe,UIm,TRe,TIm,UmatOut,TmatOut,L,ky,nWave,nmax);
      }
      else
      {
        // bemintreg3d(Nod,nNod,Elt,iElt,nElt,TypeID,nKeyOpt,TypeName,TypeKeyOpts,
                    // nEltType,CollPoints,nTotalColl,RegularColl,eltCollIndex,nDof,greenPtr,
                    // nGrSet,ugCmplx,tgCmplx,tg0Cmplx,URe,UIm,TRe,TIm,TmatOut,
					// s,ms,ns,
					// scolli,scompi,uniquescolli,Nuniquescolli,nuniquescolli,uniquescolliind,
					// scollj,scompj,uniquescollj,Nuniquescollj,nuniquescollj,uniquescolljind,InListuniquecollj);

  // time_t  start_bemmat_elt = clock();  
	/*				
	bemintreg3dnodiag(Nod,nNod,Elt,iElt,nElt,TypeID,nKeyOpt,TypeName,TypeKeyOpts,
                    nEltType,CollPoints,nTotalColl,RegularColl,eltCollIndex,nDof,greenPtr,
                    nGrSet,ugCmplx,tgCmplx,tg0Cmplx,URe,UIm,TRe,TIm,UmatOut,TmatOut,
					spassed,ms,ns,
					scompi,uniquescolli,Nuniquescolli,nuniquescolli,uniquescolliind,
					scollj,scompj,InListuniquecollj,
					NULL,0,blocks,NEltCollConsider);
	*/
	// for(int i=0; i<nXi_loc*nEltNod[iElt]; i++) {
                     // mexPrintf("N[%d]: %f \n",i,N_loc[i]); // DEBUG
	// }
	// for(int i=0; i<2*nXi_loc*nEltNod[iElt]; i++) {
                     // mexPrintf("dN[%d]: %f \n",i,dN_loc[i]); // DEBUG
	// }
	
	// 26/01/2012
	// bemintreg3dnodiag(Nod,nNod,Elt,iElt,nElt,TypeID,nKeyOpt,TypeName,TypeKeyOpts,
                    // nEltType,CollPoints,nTotalColl,RegularColl_loc,eltCollIndex_loc,nDof,greenPtr,
                    // nGrSet,ugCmplx,tgCmplx,tg0Cmplx,URe,UIm,TRe,TIm,UmatOut,TmatOut,
					// spassed,ms,ns,
					// scompi,uniquescolli,Nuniquescolli,nuniquescolli,uniquescolliind,
					// scollj,scompj,InListuniquecollj,
					// NULL,0,blocks,NEltCollConsider,
					// EltParent,nEltNod,nEltColl,EltShapeN,EltShapeM,EltDim,AxiSym,Periodic,nGauss,nEltDiv,nGaussSing,nEltDivSing,
					// EltNod_loc,
					// nXi_loc,xi_loc,H_loc,N_loc,M_loc,dN_loc);
		
        // mexPrintf("bemintreg3dnodiag not used ... \n");
		bemintreg3dnodiag(
//					Nod,nNod,
					Elt,iElt,nElt,
//					TypeID,nKeyOpt,TypeName,TypeKeyOpts,
//                    nEltType,
					CollPoints,nTotalColl,RegularColl_loc,eltCollIndex_loc,nDof,greenPtr,
                    nGrSet,ugCmplx,tgCmplx,tg0Cmplx,URe,UIm,TRe,TIm,UmatOut,TmatOut,
					spassed,ms,ns,
					scompi,uniquescolli,Nuniquescolli,nuniquescolli,uniquescolliind,
					scollj,scompj,InListuniquecollj,DeltaInListuniquecollj,
//					NULL,0,
					blocks,NEltCollConsider,
//					EltParent,
					nEltNod,nEltColl,
//					EltShapeN,EltShapeM,
					EltDim,
//					AxiSym,Periodic,nGauss,nEltDiv,nGaussSing,nEltDivSing,
					EltNod_loc,
					nXi_loc,
//					xi_loc,
					H_loc,N_loc,M_loc,dN_loc);
        
	//mexPrintf("Regular %d \n",iElt);
// float time_bemmat_elt = (float) (clock() - start_bemmat_elt) / CLOCKS_PER_SEC; 
// timeTest+=time_bemmat_elt;

// mexPrintf("time for bemmat_elt was %f seconds\n", time_bemmat_elt);	
// mexPrintf("time for timeTest was %f seconds\n", timeTest);	

      }
    }
    else if (probDim==2)
    {
      if (probAxi)
      {
        bemintregaxi(Nod,nNod,Elt,iElt,nElt,TypeID,nKeyOpt,TypeName,TypeKeyOpts,
                     nEltType,CollPoints,nTotalColl,RegularColl_loc,eltCollIndex_loc,nDof,greenPtr,
                     nGrSet,ugCmplx,tgCmplx,tg0Cmplx,URe,UIm,TRe,TIm,UmatOut,TmatOut);
      }
      else
      {
        bemintreg2d(Nod,nNod,Elt,iElt,nElt,TypeID,nKeyOpt,TypeName,
                    TypeKeyOpts,nEltType,CollPoints,nTotalColl,RegularColl_loc,
                    eltCollIndex_loc,nDof,greenPtr,nGrSet,nugComp,ugCmplx,
                    tgCmplx,tg0Cmplx,URe,UIm,TRe,TIm,TmatOut);
      }
    }

	

    // SINGULAR INTEGRATION
	
    
    double* const eltNodXi=new(nothrow) double[2*nEltNod[iElt]];
    if (eltNodXi==0) throw("Out of memory.");
    eltnoddef(EltType,TypeID,TypeName,nEltType,eltNodXi);

	unsigned int nuniquescollicumul=0;
	
	
	//  s  not empty
	if (spassed)
	{
		for (unsigned int iuniquescolli=0; iuniquescolli<Nuniquescolli[0]; iuniquescolli++) // bepalen welke nodig
		{		
		
		// int* const inddiag=new(nothrow) int[9];
			// if (inddiag==0) throw("Out of memory."); 

		// for (int iscolliOnDiag=0; iscolliOnDiag<ms*ns; iscolliOnDiag++)
		// {
			// if (scolli[iscolliOnDiag]==uniquescolli[iuniquescolli] && scolliOnDiag[iscolliOnDiag]==true)
			// {
				// inddiag[3*scompi[iscolliOnDiag]+scompj[iscolliOnDiag]] = iscolliOnDiag;
			// }
		// }
		
		
		
			if (RegularColl_loc[uniquescolli[iuniquescolli]]==0)
			{
			  	// mexPrintf("Singular integration required...: \t [iuniquescolli] : %d \n",uniquescolli[iuniquescolli]);
				
				
				if (probDim==3)
				{
					if ((CollPoints[uniquescolli[iuniquescolli]]==1) && (EltParent[iElt]==1)) // Triangle element centroid;
					{
						xiSing[0]=3.333333333333333e-01;
						xiSing[1]=3.333333333333333e-01;
					}
					else if (CollPoints[uniquescolli[iuniquescolli]]==1)  // Quadrilateral element centroid or line element;
					{
						xiSing[0]=0.0;
						xiSing[1]=0.0;
					}
					else if (CollPoints[uniquescolli[iuniquescolli]]==2)
					{
						int iEltNod=RegularColl_loc[nTotalColl+uniquescolli[iuniquescolli]];
						xiSing[0]=eltNodXi[0*nEltNod[iElt]+iEltNod];
						xiSing[1]=eltNodXi[1*nEltNod[iElt]+iEltNod];
					}
          
					if (probPeriodic)
					{		
						bemintsing3dperiodic(Nod,nNod,Elt,iElt,nElt,TypeID,nKeyOpt,TypeName,
											TypeKeyOpts,nEltType,CollPoints,nTotalColl,uniquescolli[iuniquescolli],
											eltCollIndex_loc,nDof,xiSing,greenPtr,nGrSet,ugCmplx,
											tgCmplx,tg0Cmplx,URe,UIm,TRe,TIm,UmatOut,TmatOut,nWave);
					}
					else 
					{
						// bemintsing3d(Nod,nNod,Elt,iElt,nElt,TypeID,nKeyOpt,TypeName,
									// TypeKeyOpts,nEltType,CollPoints,nTotalColl,uniquescolli[iuniquescolli],
									// eltCollIndex,nDof,xiSing,greenPtr,nGrSet,ugCmplx,
									// tgCmplx,tg0Cmplx,URe,UIm,TRe,TIm,TmatOut,
									// s,ms,ns,
									// scolli,scompi,uniquescolli,Nuniquescolli,nuniquescolli,uniquescolliind,
									// scollj,scompj,uniquescollj,Nuniquescollj,nuniquescollj,uniquescolljind,InListuniquecollj,nuniquescollicumul);
						// mexPrintf("sing3d...\n"); // DEBUG	
						
						// time_t  start_bemmat_elt_sing = clock();  
						
						// 26/01/2012
						// bemintsing3d(Nod,nNod,Elt,iElt,nElt,TypeID,nKeyOpt,TypeName,
									// TypeKeyOpts,nEltType,CollPoints,nTotalColl,uniquescolli[iuniquescolli],iuniquescolli,
									// eltCollIndex_loc,nDof,xiSing,greenPtr,nGrSet,ugCmplx,
									// tgCmplx,tg0Cmplx,URe,UIm,TRe,TIm,UmatOut,TmatOut,
									// spassed,ms,ns,
									// scompi,nuniquescolli,uniquescolliind,
									// scollj,scompj,InListuniquecollj,nuniquescollicumul,
									// inddiag,0,blockdiag,blocks,NEltCollConsider,
									// EltParent,nEltNod,nEltColl,EltShapeN,EltShapeM,EltDim,AxiSym,Periodic,nGauss,nEltDiv,nGaussSing,nEltDivSing,
									// EltNod_loc);
						// mexPrintf("bemintsing3d not used ... \n");			
						bemintsing3d(
									// Nod,nNod,
									Elt,iElt,nElt,
									// TypeID,nKeyOpt,TypeName,TypeKeyOpts,nEltType,
									CollPoints,nTotalColl,uniquescolli[iuniquescolli],iuniquescolli,
									eltCollIndex_loc,nDof,xiSing,greenPtr,nGrSet,ugCmplx,
									tgCmplx,tg0Cmplx,URe,UIm,TRe,TIm,UmatOut,TmatOut,
									spassed,ms,ns,
									scompi,nuniquescolli,uniquescolliind,
									scollj,scompj,InListuniquecollj,DeltaInListuniquecollj,nuniquescollicumul,
									inddiag,0,blockdiag,blocks,NEltCollConsider,
									EltParent,nEltNod,nEltColl,EltShapeN,EltShapeM,EltDim,
									// AxiSym,Periodic,nGauss,nEltDiv,
									nGaussSing,nEltDivSing,
									EltNod_loc);
						
						// float time_bemmat_elt_sing = (float) (clock() - start_bemmat_elt_sing) / CLOCKS_PER_SEC; 
						// timeTest_sing+=time_bemmat_elt_sing;

						// mexPrintf("time for bemmat_elt_sing was %f seconds\n", time_bemmat_elt_sing);	
						// mexPrintf("time for timeTest_sing was %f seconds\n", timeTest_sing);	

									
					}
          
				}
				else if (probDim==2)
				{
					if (probAxi)
					{
						bemintsingaxi(Nod,nNod,Elt,iElt,nElt,TypeID,nKeyOpt,TypeName,
									TypeKeyOpts,nEltType,CollPoints,nTotalColl,uniquescolli[iuniquescolli],
									eltCollIndex_loc,nDof,greenPtr,nGrSet,ugCmplx,
									tgCmplx,tg0Cmplx,URe,UIm,TRe,TIm,UmatOut,TmatOut);
					}
					else
					{
						bemintsing2d(Nod,nNod,Elt,iElt,nElt,TypeID,nKeyOpt,TypeName,TypeKeyOpts,
									nEltType,CollPoints,nTotalColl,uniquescolli[iuniquescolli],eltCollIndex_loc,
									nDof,greenPtr,nGrSet,nugComp,ugCmplx,tgCmplx,
									tg0Cmplx,URe,UIm,TRe,TIm,UmatOut,TmatOut);
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
	// int* const inddiag=new(nothrow) int[9];
			// if (inddiag==0) throw("Out of memory."); 
	
    for (unsigned int iColl=0; iColl<nTotalColl; iColl++)
    {
      if (RegularColl_loc[iColl]==0)
      {
		
        
		if (probDim==3)
        {
          if ((CollPoints[iColl]==1) && (EltParent[iElt]==1)) // Triangle element centroid;
          {
            xiSing[0]=3.333333333333333e-01;
            xiSing[1]=3.333333333333333e-01;
          }
          else if (CollPoints[iColl]==1)  // Quadrilateral element centroid or line element;
          {
            xiSing[0]=0.0;
            xiSing[1]=0.0;
          }
          else if (CollPoints[iColl]==2)
          {
            unsigned int iEltNod=RegularColl_loc[nTotalColl+iColl];
            xiSing[0]=eltNodXi[0*nEltNod[iElt]+iEltNod];
            xiSing[1]=eltNodXi[1*nEltNod[iElt]+iEltNod];
          }
          
          if (probPeriodic)
          {
            bemintsing3dperiodic(Nod,nNod,Elt,iElt,nElt,TypeID,nKeyOpt,TypeName,
                                 TypeKeyOpts,nEltType,CollPoints,nTotalColl,iColl,
                                 eltCollIndex_loc,nDof,xiSing,greenPtr,nGrSet,ugCmplx,
                                 tgCmplx,tg0Cmplx,URe,UIm,TRe,TIm,UmatOut,TmatOut,nWave);
          }
          else 
          {
            // bemintsing3d(Nod,nNod,Elt,iElt,nElt,TypeID,nKeyOpt,TypeName,
                        // TypeKeyOpts,nEltType,CollPoints,nTotalColl,iColl,
                        // eltCollIndex,nDof,xiSing,greenPtr,nGrSet,ugCmplx,
                        // tgCmplx,tg0Cmplx,URe,UIm,TRe,TIm,TmatOut,
						// s,ms,ns,
						// scolli,scompi,uniquescolli,Nuniquescolli,nuniquescolli,uniquescolliind,
						// scollj,scompj,uniquescollj,Nuniquescollj,nuniquescollj,uniquescolljind,InListuniquecollj,nuniquescollicumul);

		// mexPrintf("Test ... \n");
		
		// time_t  start_bemmat_elt_sing = clock();  
		
		// 26/01/2012
		// bemintsing3d(Nod,nNod,Elt,iElt,nElt,TypeID,nKeyOpt,TypeName,
                        // TypeKeyOpts,nEltType,CollPoints,nTotalColl,iColl,iColl,
                        // eltCollIndex_loc,nDof,xiSing,greenPtr,nGrSet,ugCmplx,
                        // tgCmplx,tg0Cmplx,URe,UIm,TRe,TIm,UmatOut,TmatOut,
						// spassed,ms,ns,
						// scompi,nuniquescolli,uniquescolliind,
						// scollj,scompj,InListuniquecollj,nuniquescollicumul,
						// inddiag,0,blockdiag,blocks,NEltCollConsider,
						// EltParent,nEltNod,nEltColl,EltShapeN,EltShapeM,EltDim,AxiSym,Periodic,nGauss,nEltDiv,nGaussSing,nEltDivSing,
						// EltNod_loc);

						
						
			bemintsing3d(
						// Nod,nNod,
						Elt,iElt,nElt,
						// TypeID,nKeyOpt,TypeName,TypeKeyOpts,nEltType,
						CollPoints,nTotalColl,iColl,iColl,
                        eltCollIndex_loc,nDof,xiSing,greenPtr,nGrSet,ugCmplx,
                        tgCmplx,tg0Cmplx,URe,UIm,TRe,TIm,UmatOut,TmatOut,
						spassed,ms,ns,
						scompi,nuniquescolli,uniquescolliind,
						scollj,scompj,InListuniquecollj,DeltaInListuniquecollj,nuniquescollicumul,
						inddiag,0,blockdiag,blocks,NEltCollConsider,
						EltParent,nEltNod,nEltColl,EltShapeN,EltShapeM,EltDim,
						// AxiSym,Periodic,nGauss,nEltDiv,
						nGaussSing,nEltDivSing,
						EltNod_loc);
						
						
			
			// float time_bemmat_elt_sing = (float) (clock() - start_bemmat_elt_sing) / CLOCKS_PER_SEC; 
			// timeTest_sing+=time_bemmat_elt_sing;

			// mexPrintf("time for bemmat_elt_sing was %f seconds\n", time_bemmat_elt_sing);	
			// mexPrintf("time for timeTest_sing was %f seconds\n", timeTest_sing);	

						
		  }
          
        }
        else if (probDim==2)
        {
          if (probAxi)
          {
            bemintsingaxi(Nod,nNod,Elt,iElt,nElt,TypeID,nKeyOpt,TypeName,
                          TypeKeyOpts,nEltType,CollPoints,nTotalColl,iColl,
                          eltCollIndex_loc,nDof,greenPtr,nGrSet,ugCmplx,
                          tgCmplx,tg0Cmplx,URe,UIm,TRe,TIm,UmatOut,TmatOut);
          }
          else
          {
            bemintsing2d(Nod,nNod,Elt,iElt,nElt,TypeID,nKeyOpt,TypeName,TypeKeyOpts,
                         nEltType,CollPoints,nTotalColl,iColl,eltCollIndex_loc,
                         nDof,greenPtr,nGrSet,nugComp,ugCmplx,tgCmplx,
                         tg0Cmplx,URe,UIm,TRe,TIm,UmatOut,TmatOut);
          }
        }
      }
    }
	
	// delete [] inddiag;  
	
	}
	
	
	
	
    // delete [] eltCollIndex;
	delete [] eltCollIndex_loc;
	delete [] RegularColl_loc;
	delete [] EltNod_loc;
    delete [] eltNodXi;  
	// delete [] xi_loc;
	// delete [] H_loc;
	// float time5 = (float) (clock() - time_perelem) / CLOCKS_PER_SEC;
	// mexPrintf("time for time_perelem was %f seconds\n", time5);		
    }
	
	EltType_prev=EltType;
	
  }
  
  
  delete [] xi_loc;
  delete [] H_loc;
  delete [] N_loc;
  delete [] M_loc;
  delete [] dN_loc;
  // */
  
  // 

// float time2 = (float) (clock() - start_time_loop) / CLOCKS_PER_SEC; 
// mexPrintf("time for loop1 was %f seconds\n", time2);	

	
	
	
	
	 // time_t  start_time_loop_diag = clock();  
	float timeTest_3ddiag=0.0;
	float timeTest_diag_sing=0.0;
	
	
	if (s!=0 && TmatOut==true && ondiag==true)
	{
		// mexPrintf(" In loopdiag ...\n"); // DEBUG
		
			double* const scolliOnDiagUnique=new(nothrow) double[NOnDiag];
			if (scolliOnDiagUnique==0) throw("Out of memory."); 
	
			for (unsigned int iscolliOnDiagUnique=0; iscolliOnDiagUnique<NOnDiag; iscolliOnDiagUnique++)
			{
				scolliOnDiagUnique[iscolliOnDiagUnique]=-1;
			}
	
			unsigned int NOnDiagUnique=0;
			
				// // // // // 
				for (unsigned int iOnDiag=0; iOnDiag<ms*ns;iOnDiag++)
				{
					if (scolliOnDiag[iOnDiag]==true)
					{
						bool key=false;
						unsigned int i=0;
			
						while (key==false && i<NOnDiagUnique)
						{
							if (((unsigned int) (s[iOnDiag]-1) % nDof)==scolliOnDiagUnique[i])
							{
								scolliOnDiag[iOnDiag]=false;
								key=true;
							}
							i++;
						}
						if(key==false)
						{	
								scolliOnDiagUnique[NOnDiagUnique]=((unsigned int) (s[iOnDiag]-1) % nDof);
								NOnDiagUnique++;
						}			
					}
				}
				// // // // // 
		
			
			
			// mexPrintf("NOnDiag : %d \n",NOnDiag); // DEBUG
			// mexPrintf("NOnDiagUnique : %d \n",NOnDiagUnique); // DEBUG
			// mexPrintf("nTotalColl : %d \n",nTotalColl); // DEBUG
			
			double* const sdiag=new(nothrow) double[NOnDiagUnique*nDof];
			// // // // // double* const sdiag=new(nothrow) double[NOnDiag*nDof];
			if (sdiag==0) throw("Out of memory."); 
		
			unsigned int OnDiagcumul=0;
	
	
			for (unsigned int iOnDiag=0; iOnDiag<ms*ns; iOnDiag++)
			{
				if (scolliOnDiag[iOnDiag]==true)
				{	
					uint64 aa =  (uint64) (s[iOnDiag]);
					uint64 bb = (uint64) (aa-1) % nDof;
					// mexPrintf("Mod: %lld \n",bb);
					for (unsigned int iscolliOnDiag=0;iscolliOnDiag<nDof;iscolliOnDiag++)
					{
						// sdiag[nDof*iOnDiag+iscolliOnDiag]=(int (s[iOnDiag]-1) % nDof)+iscolliOnDiag;
						// // sdiag[nDof*OnDiagcumul+iscolliOnDiag]=(int (s[iOnDiag]-1) % nDof)+iscolliOnDiag*nDof+1;
						// // // // // sdiag[NOnDiagUnique*iscolliOnDiag+OnDiagcumul]=(int (s[iOnDiag]-1) % nDof)+iscolliOnDiag*nDof+1;
						// // // // // sdiag[NOnDiagUnique*iscolliOnDiag+OnDiagcumul]=(uint64 ( ((uint64) (s[iOnDiag]-1))) % nDof)+iscolliOnDiag*nDof+1;
						// uint64 aa =  (uint64) (s[iOnDiag]);
						// bb = (aa-1) % nDof;
						uint64 cc= (uint64) iscolliOnDiag*nDof;
						sdiag[NOnDiagUnique*iscolliOnDiag+OnDiagcumul]=(uint64)((bb)+(cc)+1);
						// // // sdiag[NOnDiagUnique*iscolliOnDiag+OnDiagcumul]=(uint64)((bb)+iscolliOnDiag*nDof+1);
					// // // // // sdiag[NOnDiag*iscolliOnDiag+OnDiagcumul]=(int (s[iOnDiag]-1) % nDof)+iscolliOnDiag*nDof+1;
					
						// if (iscolliOnDiag==0)
						// {
						// mexPrintf("sdiag[%d] : %f \n",NOnDiagUnique*iscolliOnDiag+OnDiagcumul,sdiag[NOnDiagUnique*iscolliOnDiag+OnDiagcumul]); // DEBUG
						// }
						// if (iscolliOnDiag==(nDof-1))
						// {
						// mexPrintf("cc: %lld \n",cc);
						// mexPrintf("sdiag[%d] : %f \n",NOnDiagUnique*iscolliOnDiag+OnDiagcumul,sdiag[NOnDiagUnique*iscolliOnDiag+OnDiagcumul]); // DEBUG
						// }
					}
				OnDiagcumul++;	
				}	
			}
	
	

	
	
	// // // // // mexPrintf(" Check 1...\n"); // DEBUG
	
	
			unsigned int* const sdiagcolli=new(nothrow) unsigned int[NOnDiagUnique*nDof];
				if (sdiagcolli==0) throw("Out of memory sdiagcolli.");
			unsigned int* const sdiagcompi=new(nothrow) unsigned int[NOnDiagUnique*nDof];
				if (sdiagcompi==0) throw("Out of memory sdiagcompi.");
			// int* const uniquesdiagcolli=new(nothrow) int[NOnDiagUnique*nDof];
			unsigned int* const uniquesdiagcolli=new(nothrow) unsigned int[nDof];
				if (uniquesdiagcolli==0) throw("Out of memory uniquesdiagcolli.");
			unsigned int* const Nuniquesdiagcolli=new(nothrow) unsigned int[1];
				if (Nuniquesdiagcolli==0) throw("Out of memory Nuniquesdiagcolli.");       
			// int* const nuniquesdiagcolli=new(nothrow) int[NOnDiagUnique*nDof];
			unsigned int* const nuniquesdiagcolli=new(nothrow) unsigned int[nDof];
				if (nuniquesdiagcolli==0) throw("Out of memory nuniquesdiagcolli.");       
			unsigned int* const uniquesdiagcolliind=new(nothrow) unsigned int[NOnDiagUnique*nDof];
				if (uniquesdiagcolliind==0) throw("Out of memory uniquesdiagcolliind.");     
       
			unsigned int* const sdiagcollj=new(nothrow) unsigned int[NOnDiagUnique*nDof];
				if (sdiagcollj==0) throw("Out of memory scollj.");
			unsigned int* const sdiagcompj=new(nothrow) unsigned int[NOnDiagUnique*nDof];
				if (sdiagcompj==0) throw("Out of memory sdiagcompj.");
			// int* const uniquesdiagcollj=new(nothrow) int[NOnDiagUnique*nDof];
			unsigned int* const uniquesdiagcollj=new(nothrow) unsigned int[nDof];
				if (uniquesdiagcollj==0) throw("Out of memory uniquesdiagcollj.");
			unsigned int* const Nuniquesdiagcollj=new(nothrow) unsigned int[1];
				if (Nuniquesdiagcollj==0) throw("Out of memory Nuniquesdiagcollj.");       
			// int* const nuniquesdiagcollj=new(nothrow) int[NOnDiagUnique*nDof];
			unsigned int* const nuniquesdiagcollj=new(nothrow) unsigned int[nDof];
				if (nuniquesdiagcollj==0) throw("Out of memory nuniquesdiagcollj.");       
			unsigned int* const uniquesdiagcolljind=new(nothrow) unsigned int[NOnDiagUnique*nDof];
				if (uniquesdiagcolljind==0) throw("Out of memory uniquesdiagcolljind.");    

			bool* const InListuniquediagcollj=new(nothrow) bool[nTotalColl];
				if (InListuniquediagcollj==0) throw("Out of memory InListuniquediagcollj."); 
				
			int* const DeltaInListuniquediagcollj=new(nothrow) int[nTotalColl];
				if (DeltaInListuniquediagcollj==0) throw("Out of memory DeltaInListuniquediagcollj."); 

			bool* const sdiagcolliOnDiag=new(nothrow) bool[NOnDiagUnique*nDof];
				if (sdiagcolliOnDiag==0) throw("Out of memory sdiagcolliOnDiag."); 
				

		// for (int i=0;i<NOnDiagUnique*nDof;i++)
		// {
			// sdiagcolli[i]=0;
			// sdiagcompi[i]=0;
			// sdiagcollj[i]=0;
			// sdiagcompj[i]=0;
		// }	
				
				
			// unsigned int a=1;
	
   // time_t  start_s2coll_diag = clock();  	
			s2coll(sdiag,1,NOnDiagUnique*nDof,nDof,probDim,sdiagcolli,sdiagcompi,uniquesdiagcolli,Nuniquesdiagcolli,nuniquesdiagcolli,uniquesdiagcolliind,
										      sdiagcollj,sdiagcompj,uniquesdiagcollj,Nuniquesdiagcollj,nuniquesdiagcollj,uniquesdiagcolljind,
											  sdiagcolliOnDiag);
   // // // // s2coll(sdiag,NOnDiagUnique,nDof,nDof,probDim,sdiagcolli,sdiagcompi,uniquesdiagcolli,Nuniquesdiagcolli,nuniquesdiagcolli,uniquesdiagcolliind,
										      // // // // sdiagcollj,sdiagcompj,uniquesdiagcollj,Nuniquesdiagcollj,nuniquesdiagcollj,uniquesdiagcolljind,
											  // // // // sdiagcolliOnDiag);											  
   // float time_s2coll_diag = (float) (clock() - start_s2coll_diag) / CLOCKS_PER_SEC; 
   // mexPrintf("time for s2coll_diag was %f seconds\n", time_s2coll_diag);
	
	 // time_t  start_InListDiag = clock();  
	// for (int i=0; i<NOnDiagUnique*nDof; i++)
	// {
		// mexPrintf("sdiagcolli[%d]: %u \n",i,sdiagcolli[i]);
	// // }	 
	// for (int i=0; i<441; i++)
	// {
		// mexPrintf("uniquesdiagcolli[%d]: %u \n",i,uniquesdiagcolli[i]);
	// }	
	// for (int i=0; i<441; i++)
	// {
		// mexPrintf("uniquesdiagcollj[%d]: %u \n",i,uniquesdiagcollj[i]);
	// }	
		
	// for (int i=0;i<NOnDiagUnique*nDof;i++)
	// {
		// if (sdiagcolli[i]>10000)
		// {
		// mexPrintf("sdiag[%d] : %f \n",i,sdiag[i]); // DEBUG
		// }
	// }
		
		// mexPrintf("Nuniquescolli : %u \n",Nuniquescolli[0]); // DEBUG
		// mexPrintf("Nuniquesdiagcolli : %u \n",Nuniquesdiagcolli[0]); // DEBUG
		
		
		
			for (unsigned int iInListuniquediagcollj=0; iInListuniquediagcollj<nTotalColl; iInListuniquediagcollj++)   
			{
				
				// mexPrintf("iInListuniquediagcollj : %d \n",iInListuniquediagcollj); // DEBUG
				InListuniquediagcollj[iInListuniquediagcollj]=0;
			}
			for (unsigned int iuniquesdiagcollj=0; iuniquesdiagcollj<Nuniquesdiagcollj[0]; iuniquesdiagcollj++)
			{
				// mexPrintf("uniquesdiagcollj[iuniquesdiagcollj] : %d \n",uniquesdiagcollj[iuniquesdiagcollj]);
				InListuniquediagcollj[uniquesdiagcollj[iuniquesdiagcollj]]=1;
			}
			
			DeltaInListuniquediagcollj[0]=0; 
			for (unsigned int iDeltaInListuniquediagcollj=1; iDeltaInListuniquediagcollj<nTotalColl; iDeltaInListuniquediagcollj++)   
			{
			// mexPrintf("iDeltaInListuniquediagcollj : %d \n",iDeltaInListuniquediagcollj); // DEBUG
				if (InListuniquediagcollj[iDeltaInListuniquediagcollj-1]==0)
				{
					DeltaInListuniquediagcollj[iDeltaInListuniquediagcollj]=DeltaInListuniquediagcollj[iDeltaInListuniquediagcollj-1]+1;
				}
				else
				{
					DeltaInListuniquediagcollj[iDeltaInListuniquediagcollj]=DeltaInListuniquediagcollj[iDeltaInListuniquediagcollj-1];
				}
			}
	
			
			
			
			// mexPrintf(" Check 3...\n"); // DEBUG
			
	// // // // //
	
		// // // // // for (unsigned int iOnDiag=0; iOnDiag<ms*ns;iOnDiag++)
		// // // // // {
			// // // // // if (scolliOnDiag[iOnDiag]==true)
			// // // // // {
				// // // // // bool key=false;
				// // // // // unsigned int i=0;
				
				// // // // // while (key==false && i<NOnDiagUnique)
				// // // // // {
					// // // // // if (((unsigned int) (s[iOnDiag]-1) % nDof)==scolliOnDiagUnique[i])
					// // // // // {
						// // // // // scolliOnDiag[iOnDiag]=false;
						// // // // // key=true;
					// // // // // }
					// // // // // i++;
				// // // // // }
				// // // // // if(key==false)
				// // // // // {	
						// // // // // scolliOnDiagUnique[NOnDiagUnique]=((unsigned int) (s[iOnDiag]-1) % nDof);
						// // // // // NOnDiagUnique++;
				// // // // // }			
			// // // // // }
		// // // // // }
		
	// // // // //
		
		

		
	// mexPrintf(" Check 4...\n"); // DEBUG		
		
		
		EltType_prev = (unsigned int)(Elt[nElt+0]);
		nXi_loc=nXi[EltType_prev-1];
		
		double* const xi_loc=new(nothrow) double[2*nXi_loc];
		if (xi_loc==0) throw("Out of memory.");
		
		double* const H_loc=new(nothrow) double[nXi_loc];
		if (H_loc==0) throw("Out of memory.");
		
		for(unsigned int iXi=0; iXi<2*nXi_loc; iXi++)
		{
			xi_loc[iXi] = xi[2*ncumulnXi[EltType_prev-1]+iXi];
		}
		
		for(unsigned int iH=0; iH<nXi_loc; iH++)
		{
			H_loc[iH] = H[ncumulnXi[EltType_prev-1]+iH];
		}
		
		double* const N_loc=new(nothrow) double[nXi_loc*nEltNod[RefEltType[EltType_prev-1]]];
		if (N_loc==0) throw("Out of memory.");
		
		double* const M_loc=new(nothrow) double[nXi_loc*nEltNod[RefEltType[EltType_prev-1]]];
		if (M_loc==0) throw("Out of memory.");
		
		double* const dN_loc=new(nothrow) double[2*nXi_loc*nEltNod[RefEltType[EltType_prev-1]]];
		if (dN_loc==0) throw("Out of memory.");
		
		
		for(unsigned int iN=0; iN<nXi_loc*nEltNod[RefEltType[EltType_prev-1]]; iN++)
		{
			N_loc[iN] = Nshape[ncumulNshape[EltType_prev-1]+iN];
			M_loc[iN] = Mshape[ncumulNshape[EltType_prev-1]+iN];
		}
		for(unsigned int idN=0; idN<2*nXi_loc*nEltNod[RefEltType[EltType_prev-1]]; idN++)
		{
			dN_loc[idN] = dNshape[2*ncumulNshape[EltType_prev-1]+idN];
		}
	
		NEltCollConsider=0;
		
		// mexPrintf(" Check 5...\n"); // DEBUG

		for (unsigned int iElt=0; iElt<nElt; iElt++)  // Enkel loop over nodige elementen
		{
			// mexPrintf("iElt : %d \n",iElt); // DEBUG
			//
			unsigned int* const eltCollIndex_loc=new(nothrow) unsigned int[nEltColl[iElt]];
			if (eltCollIndex_loc==0) throw("Out of memory.");
	
			for(unsigned int iEltCollIndex=0; iEltCollIndex<nEltColl[iElt]; iEltCollIndex++)
			{
				eltCollIndex_loc[iEltCollIndex] = eltCollIndex[ncumulEltCollIndex[iElt]+iEltCollIndex];
					
			}
	
			//
			unsigned int* const RegularColl_loc=new(nothrow) unsigned int[2*nTotalColl];
			if (RegularColl_loc==0) throw("Out of memory.");
		
			// // // // // for(unsigned int i=0; i<2*nTotalColl; i++)
			// // // // // {
				// // // // // RegularColl_loc[i]=RegularColl[2*nTotalColl*iElt+i];
			// // // // // }		
			
			for(unsigned int i=0; i<nTotalColl; i++)
			{
				RegularColl_loc[i]=1;
			}

			for(unsigned int iSingular=0; iSingular < nSingularColl[iElt]; iSingular++)
			{
				RegularColl_loc[RegularColl[(uint64)(ncumulSingularColl[iElt]+iSingular)]]=0;
				RegularColl_loc[nTotalColl+RegularColl[(uint64)(ncumulSingularColl[iElt]+iSingular)]]=RegularColl[(uint64)(NSingularColl+ncumulSingularColl[iElt]+iSingular)];
			}
			
			//
			double* const EltNod_loc=new(nothrow) double[3*nEltNod[iElt]];
			if (EltNod_loc==0) throw("Out of memory.");


			for(unsigned int i=0; i<3*nEltNod[iElt]; i++)
			{
				EltNod_loc[i] = EltNod[3*ncumulEltNod[iElt]+i];
			}
		
			unsigned int EltType = (unsigned int)(Elt[nElt+iElt]);
			
			
			//
			if (EltType != EltType_prev)
			{
			nXi_loc=nXi[EltType-1];
	
			delete [] xi_loc;
			double* const xi_loc=new(nothrow) double[2*nXi_loc];
			if (xi_loc==0) throw("Out of memory.");
	
			delete [] H_loc;	
			double* const H_loc=new(nothrow) double[nXi_loc];
			if (H_loc==0) throw("Out of memory.");
		
			delete [] N_loc;	
			double* const N_loc=new(nothrow) double[nXi_loc*nEltNod[RefEltType[EltType-1]]];
			if (N_loc==0) throw("Out of memory.");
		
			delete [] M_loc;
			double* const M_loc=new(nothrow) double[nXi_loc*nEltNod[RefEltType[EltType-1]]];
			if (M_loc==0) throw("Out of memory.");
		
			delete [] dN_loc;
			double* const dN_loc=new(nothrow) double[2*nXi_loc*nEltNod[RefEltType[EltType-1]]];
			if (dN_loc==0) throw("Out of memory.");
	

			for(unsigned int iXi=0; iXi<2*nXi_loc; iXi++)
			{
				xi_loc[iXi] = xi[2*ncumulnXi[EltType-1]+iXi];
			}
			for(unsigned int iH=0; iH<nXi_loc; iH++)
			{
				H_loc[iH] = H[ncumulnXi[EltType-1]+iH];
			}
			
			for(unsigned int iN=0; iN<nXi_loc*nEltNod[RefEltType[EltType-1]]; iN++)
			{
				N_loc[iN] = Nshape[ncumulNshape[EltType-1]+iN];
				M_loc[iN] = Mshape[ncumulNshape[EltType-1]+iN];
			}
			for(unsigned int idN=0; idN<2*nXi_loc*nEltNod[RefEltType[EltType-1]]; idN++)
			{
				dN_loc[idN] = dNshape[2*ncumulNshape[EltType-1]+idN];
			}

			
			
			}		
	
			
			
			/*
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
			eltdef(EltType,TypeID,TypeName,TypeKeyOpts,nKeyOpt,nEltType,EltParent,nEltNod,
				nEltColl,EltShapeN,EltShapeM,EltDim,AxiSym,Periodic,nGauss,nEltDiv,
				nGaussSing,nEltDivSing);

			int* const eltCollIndex=new(nothrow) int[nEltColl];
			if (eltCollIndex==0) throw("Out of memory.");
			BemEltCollIndex(Elt,iElt,nElt,CollPoints,nCentroidColl,nTotalColl,
							nEltColl,nEltNod,eltCollIndex);
			BemRegularColl(Elt,iElt,nElt,Nod,nNod,CoincNodes,SlavesExist,
						CollPoints,nCentroidColl,nTotalColl,RegularColl,
						nRegularColl,nSingularColl,TypeID,nKeyOpt,TypeName,
						TypeKeyOpts,nEltType);
			*/

			
			
			if (probDim==3)
			{
			if (probPeriodic){
				bemintreg3dperiodic(Nod,nNod,Elt,iElt,nElt,TypeID,nKeyOpt,TypeName,TypeKeyOpts,
									nEltType,CollPoints,nTotalColl,RegularColl_loc,eltCollIndex_loc,nDof,greenPtr,
									nGrSet,ugCmplx,tgCmplx,tg0Cmplx,URe,UIm,TRe,TIm,UmatOut,TmatOut,L,ky,nWave,nmax);
			}
			else
			{
			
				// for (int i=0; i<Nuniquescolli[0]; i++)
				// {
					// mexPrintf("blockdiag[%d]: %s \n",i, blockdiag[i] ? "true": "false");
				// }
				
				time_t  start_bemmat_elt_3ddiag = clock();  
				
				bemintreg3ddiag(Nod,nNod,Elt,iElt,nElt,TypeID,nKeyOpt,TypeName,TypeKeyOpts,
                    nEltType,CollPoints,nTotalColl,RegularColl_loc,eltCollIndex_loc,nDof,greenPtr,
                    nGrSet,ugCmplx,tgCmplx,tg0Cmplx,URe,UIm,TRe,TIm,UmatOut,TmatOut,
					sdiag,1.0,NOnDiagUnique*nDof,
					sdiagcompi,uniquesdiagcolli,Nuniquesdiagcolli,nuniquesdiagcolli,uniquesdiagcolliind,
					sdiagcollj,sdiagcompj,InListuniquediagcollj,
					inddiag,ondiag,blockdiag,
					EltParent,nEltNod,nEltColl,EltShapeN,EltShapeM,EltDim,AxiSym,Periodic,nGauss,nEltDiv,nGaussSing,nEltDivSing,
					EltNod_loc,
					nXi_loc,xi_loc,H_loc,
					N_loc,M_loc,dN_loc);
					
				float time_bemmat_elt_3ddiag = (float) (clock() - start_bemmat_elt_3ddiag) / CLOCKS_PER_SEC; 
				timeTest_3ddiag+=time_bemmat_elt_3ddiag;

				//mexPrintf("time for bemmat_elt_3ddiag was %f seconds\n", time_bemmat_elt_3ddiag);	
				//mexPrintf("time for timeTest_3ddiag was %f seconds\n", timeTest_3ddiag);
					
					
			}
			}
			else if (probDim==2)
			{
			if (probAxi)
			{
				bemintregaxi(Nod,nNod,Elt,iElt,nElt,TypeID,nKeyOpt,TypeName,TypeKeyOpts,
							nEltType,CollPoints,nTotalColl,RegularColl_loc,eltCollIndex_loc,nDof,greenPtr,
							nGrSet,ugCmplx,tgCmplx,tg0Cmplx,URe,UIm,TRe,TIm,UmatOut,TmatOut);
			}
			else
			{
				bemintreg2d(Nod,nNod,Elt,iElt,nElt,TypeID,nKeyOpt,TypeName,
							TypeKeyOpts,nEltType,CollPoints,nTotalColl,RegularColl_loc,
							eltCollIndex_loc,nDof,greenPtr,nGrSet,nugComp,ugCmplx,
							tgCmplx,tg0Cmplx,URe,UIm,TRe,TIm,TmatOut);
			}
			}
						
						

						
		double* const eltNodXi=new(nothrow) double[2*nEltNod[iElt]];
		if (eltNodXi==0) throw("Out of memory.");
		eltnoddef(EltType,TypeID,TypeName,nEltType,eltNodXi);

		// mexPrintf("iElt : %d \n",iElt); // DEBUG		
		
		unsigned int nuniquesdiagcollicumul=0;
		for (unsigned int iuniquesdiagcolli=0; iuniquesdiagcolli<Nuniquesdiagcolli[0]; iuniquesdiagcolli++) // bepalen welke nodig
		{
		// mexPrintf("uniquesdiagcolli[iuniquesdiagcolli]] : %d \n",uniquesdiagcolli[iuniquesdiagcolli]); // DEBUG		
		// mexPrintf("RegularColl_loc[uniquesdiagcolli[iuniquesdiagcolli]] : %d \n",RegularColl_loc[uniquesdiagcolli[iuniquesdiagcolli]]); // DEBUG		
		if (RegularColl_loc[uniquesdiagcolli[iuniquesdiagcolli]]==0)
		{
				// mexPrintf("CollPoints[uniquesdiagcolli[iuniquesdiagcolli]] : %d \n",CollPoints[uniquesdiagcolli[iuniquesdiagcolli]]); // DEBUG		
				// mexPrintf("EltParent[iElt] : %d \n",EltParent[iElt]); // DEBUG		
				if (probDim==3)
				{
					if ((CollPoints[uniquesdiagcolli[iuniquesdiagcolli]]==1) && (EltParent[iElt]==1)) // Triangle element centroid;
					{
						xiSing[0]=3.333333333333333e-01;
						xiSing[1]=3.333333333333333e-01;
					}
					else if (CollPoints[uniquesdiagcolli[iuniquesdiagcolli]]==1)  // Quadrilateral element centroid or line element;
					{
						xiSing[0]=0.0;
						xiSing[1]=0.0;
					}
					else if (CollPoints[uniquesdiagcolli[iuniquesdiagcolli]]==2)
					{
						unsigned int iEltNod=RegularColl_loc[nTotalColl+uniquesdiagcolli[iuniquesdiagcolli]];
						xiSing[0]=eltNodXi[0*nEltNod[iElt]+iEltNod];
						xiSing[1]=eltNodXi[1*nEltNod[iElt]+iEltNod];
					}
					
					if (probPeriodic)
					{		
						bemintsing3dperiodic(Nod,nNod,Elt,iElt,nElt,TypeID,nKeyOpt,TypeName,
											TypeKeyOpts,nEltType,CollPoints,nTotalColl,uniquesdiagcolli[iuniquesdiagcolli],
											eltCollIndex_loc,nDof,xiSing,greenPtr,nGrSet,ugCmplx,
											tgCmplx,tg0Cmplx,URe,UIm,TRe,TIm,UmatOut,TmatOut,nWave);
					}
					else 
					{			
						time_t  start_bemmat_elt_diag_sing = clock();  
						
						// 26/01/2012
						// bemintsing3d(Nod,nNod,Elt,iElt,nElt,TypeID,nKeyOpt,TypeName,
									// TypeKeyOpts,nEltType,CollPoints,nTotalColl,uniquesdiagcolli[iuniquesdiagcolli],iuniquesdiagcolli,
									// eltCollIndex_loc,nDof,xiSing,greenPtr,nGrSet,ugCmplx,
									// tgCmplx,tg0Cmplx,URe,UIm,TRe,TIm,UmatOut,TmatOut,
									// sdiag,1.0,NOnDiagUnique*nDof,
									// sdiagcompi,nuniquesdiagcolli,uniquesdiagcolliind,
									// sdiagcollj,sdiagcompj,InListuniquediagcollj,nuniquesdiagcollicumul,
									// inddiag,ondiag,blockdiag,blocks,NEltCollConsider,
									// EltParent,nEltNod,nEltColl,EltShapeN,EltShapeM,EltDim,AxiSym,Periodic,nGauss,nEltDiv,nGaussSing,nEltDivSing,
									// EltNod_loc);
									
						bemintsing3d(
									// Nod,nNod,
									Elt,iElt,nElt,
									// TypeID,nKeyOpt,TypeName,TypeKeyOpts,nEltType,
									CollPoints,nTotalColl,uniquesdiagcolli[iuniquesdiagcolli],iuniquesdiagcolli,
									eltCollIndex_loc,nDof,xiSing,greenPtr,nGrSet,ugCmplx,
									tgCmplx,tg0Cmplx,URe,UIm,TRe,TIm,UmatOut,TmatOut,
									sdiag,1.0,NOnDiagUnique*nDof,
									sdiagcompi,nuniquesdiagcolli,uniquesdiagcolliind,
									sdiagcollj,sdiagcompj,InListuniquediagcollj,DeltaInListuniquediagcollj,nuniquesdiagcollicumul,
									inddiag,ondiag,blockdiag,blocks,NEltCollConsider,
									EltParent,nEltNod,nEltColl,EltShapeN,EltShapeM,EltDim,
									// AxiSym,Periodic,nGauss,nEltDiv,
									nGaussSing,nEltDivSing,
									EltNod_loc);
									
						float time_bemmat_elt_diag_sing = (float) (clock() - start_bemmat_elt_diag_sing) / CLOCKS_PER_SEC; 
						timeTest_diag_sing+=time_bemmat_elt_diag_sing;

						// mexPrintf("time for bemmat_elt_diag_sing was %f seconds\n", time_bemmat_elt_diag_sing);	
						// mexPrintf("time for timeTest_diag_sing was %f seconds\n", timeTest_diag_sing);
						
					}
		  
				}
				else if (probDim==2)
				{
					if (probAxi)
					{
						bemintsingaxi(Nod,nNod,Elt,iElt,nElt,TypeID,nKeyOpt,TypeName,
									TypeKeyOpts,nEltType,CollPoints,nTotalColl,uniquesdiagcolli[iuniquesdiagcolli],
									eltCollIndex_loc,nDof,greenPtr,nGrSet,ugCmplx,
									tgCmplx,tg0Cmplx,URe,UIm,TRe,TIm,UmatOut,TmatOut);
					}
					else
					{
						bemintsing2d(Nod,nNod,Elt,iElt,nElt,TypeID,nKeyOpt,TypeName,TypeKeyOpts,
									nEltType,CollPoints,nTotalColl,uniquesdiagcolli[iuniquesdiagcolli],eltCollIndex_loc,
									nDof,greenPtr,nGrSet,nugComp,ugCmplx,tgCmplx,
									tg0Cmplx,URe,UIm,TRe,TIm,UmatOut,TmatOut);
					}
				}
	
	
			}
			
			nuniquesdiagcollicumul+=nuniquesdiagcolli[iuniquesdiagcolli];
			// delete [] inddiag;  
		

		
			}

		// delete [] eltCollIndex;
		delete [] eltCollIndex_loc;
		delete [] RegularColl_loc;
		delete [] EltNod_loc;
		delete [] eltNodXi;   
		// delete [] xi_loc;
		// delete [] H_loc;
		
						
		}
	
			// // // // // mexPrintf(" Check 6a...\n"); // DEBUG
			
		delete [] xi_loc;
		delete [] H_loc;
		delete [] N_loc;
		delete [] M_loc;
		delete [] dN_loc;

		// // // // // mexPrintf(" Check 6b...\n"); // DEBUG
	
	delete [] sdiag;  
	delete [] sdiagcolli;
	// // // // // mexPrintf(" Check 6c...\n"); // DEBUG
	delete [] sdiagcompi;
	// // // // // mexPrintf(" Check 6d...\n"); // DEBUG
	
	
	delete [] uniquesdiagcolli;  
	delete [] Nuniquesdiagcolli;

	
	delete [] nuniquesdiagcolli;
	delete [] uniquesdiagcolliind;
	delete [] sdiagcollj;
	delete [] sdiagcompj;
	delete [] uniquesdiagcollj;  
	
	
	
	delete [] Nuniquesdiagcollj;
	delete [] nuniquesdiagcollj;
	delete [] uniquesdiagcolljind;
	delete [] sdiagcolliOnDiag;
	// // // // // // // // // // delete [] scolliOnDiaginddiag;
	delete [] InListuniquediagcollj;
	delete [] scolliOnDiagUnique;
	delete [] DeltaInListuniquediagcollj;
	
		// mexPrintf(" Check 6...\n"); // DEBUG
	
	}

	   	// mexPrintf(" Check 7...\n"); // DEBUG
	 
  
  //delete [] MatDim;
  delete [] xiSing;
  // delete [] RegularColl;
  delete [] scolli;
  delete [] scompi;
  delete [] uniquescolli;  
  delete [] Nuniquescolli;
  delete [] nuniquescolli;
  delete [] uniquescolliind;
  delete [] scollj;
  delete [] scompj;
  delete [] uniquescollj;  
  delete [] Nuniquescollj;
  delete [] nuniquescollj;
  delete [] uniquescolljind;
  delete [] InListuniquecollj;
  delete [] scolliOnDiag;  
  delete [] inddiag;  
  delete [] blockdiag;
  delete [] blocks;


  delete [] DeltaInListuniquecollj;
  delete [] scolliOnDiaginddiag;
  
  // delete [] uniquescolliondiag;  


 	  // float time3 = (float) (clock() - start_time_loop_diag) / CLOCKS_PER_SEC; 
// mexPrintf("time for loopdiag was %f seconds\n", time3);			  

// float time_bemmat_loop = (float) (clock() - start_bemmat_loop) / CLOCKS_PER_SEC; 
// mexPrintf("time for bemmat_loop was %f seconds\n", time_bemmat_loop);	

// float time_bemmat = (float) (clock() - start_bemmat) / CLOCKS_PER_SEC; 
// mexPrintf("time for bemmat was %f seconds\n", time_bemmat);	

// mexPrintf("time for timeTest was %f seconds\n", timeTest);	
// mexPrintf("time for timeTest_sing was %f seconds\n", timeTest_sing);	

}
