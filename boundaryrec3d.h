#ifndef _BOUNDARYREC3D_
#define _BOUNDARYREC3D_
void boundaryRec3d(const double* const Nod, const unsigned int& nNod,
                   const double* const Elt, const unsigned int& nElt, const unsigned int& iElt,
                   const unsigned int* const TypeID,
                   const char* const TypeName[], const char* const TypeKeyOpts[],
                   const unsigned int* const nKeyOpt,
                   const unsigned int& nEltType, 
                   const double* const CollPoints,
                   const unsigned int& nTotalColl, const unsigned int& nCentroidColl,
                   const double* const Rec, 
                   const unsigned int& nRec, const unsigned int& nRecDof,
                   bool* const boundaryRec,
                   double* const TRe,const bool TmatOut,
                   const unsigned int& nDof,
                   const unsigned int& nGrSet);
#endif
