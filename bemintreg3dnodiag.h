#ifndef _BEMINTREG3DNODIAG_
#define _BEMINTREG3DNODIAG_

#ifndef _int64_
typedef long long int int64;
typedef unsigned long long int uint64;
#endif

void bemintreg3dnodiag(
				 // const double* const Nod, const int& nNod,
                 const double* const Elt, const unsigned int& iElt, const unsigned int& nElt,
                 // const int* const  TypeID, const int* const nKeyOpt,
                 // const char* const TypeName[], const char* const TypeKeyOpts[],
                 // const int& nEltType, 
				 const double* const Coll,
                 const unsigned int& nColl, const unsigned int* const RegularColl,
                 const unsigned int* const EltCollIndex, const unsigned int& nDof,
                 const void* const* const greenPtr, const unsigned int& nGrSet,
                 const bool& ugCmplx, const bool& tgCmplx,
                 const bool& tg0Cmplx, double* const URe, double* const UIm,
                 double* const TRe, double* const TIm, const bool UmatOut, const bool TmatOut,
				 // const double* const s,
				 const bool& spassed,
				 const unsigned int& ms,
            	 const unsigned int& ns,
				 // const int* const scolli,
 				 const unsigned int* const scompi,
				 const unsigned int* const uniquescolli,
				 const unsigned int* const Nuniquescolli,
				 const unsigned int* const nuniquescolli,
				 const unsigned int* const uniquescolliind,
				 const unsigned int* const scollj,
				 const unsigned int* const scompj,
				 // const int* const uniquescollj,
				 // const int* const Nuniquescollj,
				 // const int* const nuniquescollj,
				 // const int* const uniquescolljind,
				 const bool* const InListuniquecollj,
				 const int* const DeltaInListuniquecollj,
				 // const int* const inddiag,
				 // const bool& ondiag,
				 const bool* const blocks,
				 const unsigned int& NEltCollConsider,
				 // const int* const EltParent, 
				 const unsigned int* const nEltNod, const unsigned int* const nEltColl,
				 // const int* const EltShapeN, const int* const EltShapeM, 
				 const unsigned int* const EltDim,
				 // const int* const AxiSym, const int* const Periodic, const int* const nGauss,
				 // const int* const nEltDiv, const int* const nGaussSing, const int* const nEltDivSing,
				 const double* const EltNod,
				 const unsigned int& nXi, 
				 // const double* const xi, 
				 const double* const H,
				 const double* const N, const double* const M, const double* const dN);
#endif
