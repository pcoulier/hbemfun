#ifndef _BEMISPERIODIC_
#define _BEMISPERIODIC_

bool isPeriodic(const double* const Elt, const unsigned int& nElt,
                unsigned int* const  TypeID, char* TypeName[], char* TypeKeyOpts[],
                unsigned int* const nKeyOpt, const unsigned int& nEltType);
#endif
