#ifndef _S2COLL_
#define _S2COLL_

void s2coll(const double* const s,
            const unsigned int& ms,
            const unsigned int& ns,
            const unsigned int& nDof,
            const unsigned int& probDim,
            unsigned int* scolli,
            unsigned int* scompi,
            unsigned int* uniquescolli,
            unsigned int* Nuniquescolli,
            unsigned int* nuniquescolli,
            unsigned int* uniquescolliind,
            unsigned int* scollj,
            unsigned int* scompj,
            unsigned int* uniquescollj,
            unsigned int* Nuniquescollj,
            unsigned int* nuniquescollj,
            unsigned int* uniquescolljind,
			bool* scolliOnDiag);
#endif
