#ifndef _BEMNODEINDEX_
#define _BEMNODEINDEX_
void BemNodeIndex(const double* const Nod, const unsigned int& nNod,
                  const unsigned int& NodeID, int& index);
#endif

#ifndef _BEMCOLLPOINTS_
#define _BEMCOLLPOINTS_
void BemCollPoints(const double* const Elt, const double* const Nod,
                   unsigned int* const  TypeID, unsigned int* const  nKeyOpt,
                   const char* const TypeName[], const char* const TypeKeyOpts[], 
                   const unsigned int& nEltType,  
                   const unsigned int& nElt, const unsigned int& maxEltCol, const unsigned int& nNod,
                   unsigned int* const NodalColl,unsigned int* const CentroidColl,
                   unsigned int& nNodalColl, unsigned int& nCentroidColl);
#endif

#ifndef _BEMCOLLCOORDS_
#define _BEMCOLLCOORDS_
void BemCollCoords(const double* const Elt, const double* const Nod,
                   unsigned int* const  TypeID, unsigned int* const  nKeyOpt,
                   const char* const TypeName[],const char* const TypeKeyOpts[],const unsigned int& nEltType,
                   const unsigned int* const CentroidColl,
                   const unsigned int* const NodalColl,double* const CollPoints,
                   const unsigned int& nTotalColl, const unsigned int& nElt,const unsigned int& nNod);
#endif

#ifndef _BEMCOINCNODES_
#define _BEMCOINCNODES_
void BemCoincNodes(const double* const Nod, const unsigned int& nNod,
                   double* const CoincNod, bool& SlavesExist);
#endif


#ifndef _BEMELTCOLLINDEX_
#define _BEMELTCOLLINDEX_
void BemEltCollIndex(const double* const Elt, const unsigned int& iElt, const unsigned int& nElt,
                     const double* const CollPoints, const unsigned int& nCentroidColl,
                     const unsigned int& nTotalColl, const unsigned int& nEltColl,
                     const unsigned int& nEltNod, unsigned int* const eltCollIndex);
#endif


#ifndef _BEMREGULARCOLL_
#define _BEMREGULARCOLL_
void BemRegularColl(const double* const Elt,const unsigned int& iElt,
                    const unsigned int& nElt, const double* const Nod,
                    const unsigned int& nNod, const double* const CoincNod,
                    const bool& SlavesExist, const double* const CollPoints,
                    const unsigned int& nCentroidColl, const unsigned int& nTotalColl,
                    unsigned int* const RegularColl,unsigned int& nRegularColl,unsigned int& nSingularColl,
                    const unsigned int* const  TypeID, const unsigned int* const nKeyOpt,
                    const char* const TypeName[], const char* const TypeKeyOpts[],
                    const unsigned int& nEltType);
#endif
