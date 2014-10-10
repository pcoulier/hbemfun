#ifndef _BEMINT_
#define _BEMINT_
#include <complex>
void bemint(const double* const Nod, const unsigned int& nNod,
            const double* const Elt, const unsigned int& iElt,const unsigned int& nElt,
            unsigned int* const  TypeID, unsigned int* const  nKeyOpt,
            char* TypeName[], char* TypeKeyOpts[],const unsigned int& nEltType,
            const unsigned int* const eltCollIndex, std::complex<double>* const OutMat,
            const unsigned int& nDof, const std::complex<double>* const t,
            const std::complex<double>* const u, const unsigned int& ntMode, const unsigned int& nuMode,
            const unsigned int& ntSet, const unsigned int& nColDof, bool probAxi);
#endif
