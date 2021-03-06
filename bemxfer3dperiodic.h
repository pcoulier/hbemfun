#ifndef _BEMXFER3DPERIODIC_
#define _BEMXFER3DPERIODIC_
void bemxfer3dperiodic(const double* const Nod,const unsigned int& nNod,
                       const double* const Elt,const unsigned int& iElt,
                       const unsigned int& nElt, const unsigned int* const EltCollIndex,
                       const double* const Rec, const unsigned int nRec,
                       bool* const boundaryRec,
                       double* const URe, double* const UIm,
                       double* const TRe, double* const TIm,
                       const bool UmatOut,
					   const bool TmatOut,
                       const unsigned int& nDof, const unsigned int& nRecDof,
                       const unsigned int* const  TypeID, const unsigned int* const  nKeyOpt,
                       const char* const* TypeName, const char* const* TypeKeyOpts,
                       const unsigned int& nEltType,
                       const void* const* const greenPtr, const unsigned int& nGrSet, 
                       const bool& ugCmplx, const bool& tgCmplx,
                       const double L, const double* const ky, const unsigned int nWave, 
                       const unsigned int nmax);
#endif
