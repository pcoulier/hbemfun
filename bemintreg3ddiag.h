#ifndef _BEMINTREG3DDIAG_
#define _BEMINTREG3DDIAG_

#ifndef _int64_
typedef long long int int64;
typedef unsigned long long int uint64;
#endif

void bemintreg3ddiag(const double* const Nod, const unsigned int& nNod, 
                 const double* const Elt, const unsigned int& iElt, const unsigned int& nElt,
                 const unsigned int* const  TypeID, const unsigned int* const nKeyOpt, 
                 const char* const TypeName[], const char* const TypeKeyOpts[], 
                 const unsigned int& nEltType, const double* const Coll, 
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
				 const int* const inddiag,
				 const bool& ondiag,
				 const bool* const blockdiag,
				 const unsigned int* const EltParent, const unsigned int* const nEltNod, const unsigned int* const nEltColl,
				 const unsigned int* const EltShapeN, const unsigned int* const EltShapeM, const unsigned int* const EltDim,
				 const unsigned int* const AxiSym, const unsigned int* const Periodic, const unsigned int* const nGauss,
				 const unsigned int* const nEltDiv, const unsigned int* const nGaussSing, const unsigned int* const nEltDivSing,
				 const double* const EltNod,
				 const unsigned int& nXi, const double* const xi, const double* const H,
				 const double* const N, const double* const M, const double* const dN);
#endif
