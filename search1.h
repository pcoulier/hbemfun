#ifndef _SEARCH1_
#define _SEARCH1_
void search1(const double xi, const double* const x,
             const unsigned int& xend, unsigned int& x1, unsigned int& x2,
             double* const interpval, bool& extrapFlag);
void search1alt(const double xi, const double* const x, const unsigned int& xend,
                unsigned int& x1, unsigned int& x2, double* const interpval);
#endif

#ifndef _SEARCHCLOSEST_
#define _SEARCHCLOSEST_
void searchClosest(const double xi, const double* const x,
                   const unsigned int& xend, unsigned int& xind);

#endif
