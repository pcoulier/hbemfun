#ifndef _ELTDEF_
#define _ELTDEF_
void eltdef(const unsigned int& EltType, const unsigned int* const  TypeID,
            const char* const TypeName[], const char* const TypeKeyOpts[],
            const unsigned int* const  nKeyOpt, const unsigned int& nEltType, unsigned int& Parent, 
            unsigned int& nNod, unsigned int& nColl, unsigned int& ShapeTypeN, unsigned int& ShapeTypeM,
            unsigned int& EltDim, unsigned int& AxiSym, unsigned int& Periodic, unsigned int& nGauss, unsigned int& nEltDiv,
            unsigned int& nGaussSing, unsigned int& nEltDivSing);
#endif

#ifndef _ELTNODDEF_
#define _ELTNODDEF_
void eltnoddef(const unsigned int& EltType,const unsigned int* const  TypeID,
               const char* const TypeName[], const unsigned int& nEltType, 
               double* const eltNodXi);
#endif
