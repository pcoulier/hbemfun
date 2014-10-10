#ifndef _BEMDIMENSION_
#define _BEMDIMENSION_

int bemDimension(const double* const Elt, const unsigned int& nElt,
                 unsigned int* const  TypeID, char* TypeName[], char* TypeKeyOpts[],
                 unsigned int* const nKeyOpt, const unsigned int& nEltType);
#endif
