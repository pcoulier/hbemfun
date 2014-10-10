#ifndef _BEMISAXISYM_
#define _BEMISAXISYM_

bool isAxisym(const double* const Elt, const unsigned int& nElt,
              unsigned int* const  TypeID, char* TypeName[], char* TypeKeyOpts[],
              unsigned int* const nKeyOpt, const unsigned int& nEltType);
#endif
