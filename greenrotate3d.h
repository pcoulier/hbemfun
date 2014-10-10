#ifndef _GREENROTATE3D_
#define _GREENROTATE3D_
void greenrotate3d(const double* const normal,const unsigned int& iXi, 
                   const double& xiTheta,const unsigned int& nGrSet,const bool& ugCmplx,
                   const bool& tgCmplx,const bool& tg0Cmplx,
                   const double* const UgrRe,const double* const UgrIm,
                   const double* const TgrRe, const double* const TgrIm,
                   const double* const Tgr0Re, const double* const Tgr0Im,
                   double* const UXiRe, double* const UXiIm,
                   double* const TXiRe, double* const TXiIm, double* const TXi0Re,
                   double* const TXi0Im, const bool& UmatOut,const bool& TmatOut);
#endif
