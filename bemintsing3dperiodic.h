#ifndef _BEMINTSING3DPERIODIC_
#define _BEMINTSING3DPERIODIC_
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
                          const bool& UmatOut,const bool& TmatOut, const unsigned int nWave);
#endif
