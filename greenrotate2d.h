#ifndef _GREENROTATE2D_
#define _GREENROTATE2D_
void greenrotate2d(const double* const normal,const unsigned int& iXi, 
                   const unsigned int& nGrSet,const unsigned int& ntgComp,
                   const bool& tgCmplx,const bool& tg0Cmplx,
                   const double* const TgrRe, const double* const TgrIm,
                   const double* const Tgr0Re, const double* const Tgr0Im,
                   double* const TXiRe, double* const TXiIm, double* const TXi0Re,
                   double* const TXi0Im, const bool& TmatOut);
#endif
