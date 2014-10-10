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
			
			double* const sdiag=new(nothrow) double[NOnDiagUnique*nDof];
			if (sdiag==0) throw("Out of memory."); 
		
			unsigned int OnDiagcumul=0;
	
	
			for (unsigned int iOnDiag=0; iOnDiag<ms*ns; iOnDiag++)
			{
				if (scolliOnDiag[iOnDiag]==true)
				{
					// mexPrintf("Mod: %d \n",(int (s[iOnDiag]-1) % nDof));
					for (unsigned int iscolliOnDiag=0;iscolliOnDiag<nDof;iscolliOnDiag++)
					{
						// sdiag[nDof*iOnDiag+iscolliOnDiag]=(int (s[iOnDiag]-1) % nDof)+iscolliOnDiag;
						// // sdiag[nDof*OnDiagcumul+iscolliOnDiag]=(int (s[iOnDiag]-1) % nDof)+iscolliOnDiag*nDof+1;
						sdiag[NOnDiagUnique*iscolliOnDiag+OnDiagcumul]=(int (s[iOnDiag]-1) % nDof)+iscolliOnDiag*nDof+1;
					}
				OnDiagcumul++;	
				}	
			}
	
	
	mexPrintf(" Check 1...\n"); // DEBUG
	
	
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
    
	
   // time_t  start_s2coll_diag = clock();  	
			s2coll(sdiag,1.0,NOnDiagUnique*nDof,nDof,probDim,sdiagcolli,sdiagcompi,uniquesdiagcolli,Nuniquesdiagcolli,nuniquesdiagcolli,uniquesdiagcolliind,
										      sdiagcollj,sdiagcompj,uniquesdiagcollj,Nuniquesdiagcollj,nuniquesdiagcollj,uniquesdiagcolljind,
											  sdiagcolliOnDiag);
   // // // // s2coll(sdiag,NOnDiagUnique,nDof,nDof,probDim,sdiagcolli,sdiagcompi,uniquesdiagcolli,Nuniquesdiagcolli,nuniquesdiagcolli,uniquesdiagcolliind,
										      // // // // sdiagcollj,sdiagcompj,uniquesdiagcollj,Nuniquesdiagcollj,nuniquesdiagcollj,uniquesdiagcolljind,
											  // // // // sdiagcolliOnDiag);											  
   // float time_s2coll_diag = (float) (clock() - start_s2coll_diag) / CLOCKS_PER_SEC; 
   // mexPrintf("time for s2coll_diag was %f seconds\n", time_s2coll_diag);
	
	 // time_t  start_InListDiag = clock();  	
		mexPrintf(" Check 2...\n"); // DEBUG
	
			for (unsigned int iInListuniquediagcollj=0; iInListuniquediagcollj<nTotalColl; iInListuniquediagcollj++)   
			{
				InListuniquediagcollj[iInListuniquediagcollj]=0;
			}
			for (unsigned int iuniquesdiagcollj=0; iuniquesdiagcollj<Nuniquesdiagcollj[0]; iuniquesdiagcollj++)
			{
				InListuniquediagcollj[uniquesdiagcollj[iuniquesdiagcollj]]=1;
			}
			
			DeltaInListuniquediagcollj[0]=0; 
			for (unsigned int iDeltaInListuniquediagcollj=1; iDeltaInListuniquediagcollj<nTotalColl; iDeltaInListuniquediagcollj++)   
			{
				if (InListuniquediagcollj[iDeltaInListuniquediagcollj-1]==0)
				{
					DeltaInListuniquediagcollj[iDeltaInListuniquediagcollj]=DeltaInListuniquediagcollj[iDeltaInListuniquediagcollj-1]+1;
				}
				else
				{
					DeltaInListuniquediagcollj[iDeltaInListuniquediagcollj]=DeltaInListuniquediagcollj[iDeltaInListuniquediagcollj-1];
				}
			}
	
			
			
			
			mexPrintf(" Check 3...\n"); // DEBUG
			
	
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
		
		

		
	mexPrintf(" Check 4...\n"); // DEBUG		
		
		
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
		
		mexPrintf(" Check 5...\n"); // DEBUG

		for (unsigned int iElt=0; iElt<nElt; iElt++)  // Enkel loop over nodige elementen
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
		
			for(unsigned int i=0; i<2*nTotalColl; i++)
			{
				RegularColl_loc[i]=RegularColl[2*nTotalColl*iElt+i];
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

		unsigned int nuniquesdiagcollicumul=0;
		for (unsigned int iuniquesdiagcolli=0; iuniquesdiagcolli<Nuniquesdiagcolli[0]; iuniquesdiagcolli++) // bepalen welke nodig
		{
		if (RegularColl_loc[uniquesdiagcolli[iuniquesdiagcolli]]==0)
		{
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
						int iEltNod=RegularColl_loc[nTotalColl+uniquesdiagcolli[iuniquesdiagcolli]];
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

						mexPrintf("time for bemmat_elt_diag_sing was %f seconds\n", time_bemmat_elt_diag_sing);	
						mexPrintf("time for timeTest_diag_sing was %f seconds\n", timeTest_diag_sing);
						
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
	
			mexPrintf(" Check 6a...\n"); // DEBUG
			
		delete [] xi_loc;
		delete [] H_loc;
		delete [] N_loc;
		delete [] M_loc;
		delete [] dN_loc;

		mexPrintf(" Check 6b...\n"); // DEBUG
	
	delete [] sdiag;  
	delete [] sdiagcolli;
	mexPrintf(" Check 6c...\n"); // DEBUG
	delete [] sdiagcompi;
	mexPrintf(" Check 6d...\n"); // DEBUG
	
	
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
	// // // // // delete [] scolliOnDiaginddiag;
	delete [] InListuniquediagcollj;
	delete [] scolliOnDiagUnique;
	delete [] DeltaInListuniquediagcollj;
	
		mexPrintf(" Check 6...\n"); // DEBUG
	
	}