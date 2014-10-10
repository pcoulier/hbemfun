#ifndef _GAUSSPW1D_NODIV_
#define _GAUSSPW1D_NODIV_
void gausspw1D_nodiv(const int& nGauss, double* const xi,double* const H);
#endif

#ifndef _GAUSSPW1D_
#define _GAUSSPW1D_
void gausspw1D(const int& nEltDiv, const int& nGauss, double* const xi,double* const H);
#endif

#ifndef _GAUSSPW2D_
#define _GAUSSPW2D_
void gausspw2D(const int& nEltDiv, const int& nGauss, double* const xi,
               double* const H);
#endif

#ifndef _GAUSSPWTRI_
#define _GAUSSPWTRI_
void gausspwtri(const int& nGauss, double* const xi ,double* const H);
#endif
