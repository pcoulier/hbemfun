#include "search1.h"
#include <complex>
#include "math.h"
#include "fsgreenf.h"
#include "fsgreen2d_inplane.h"
#include "fsgreen2d_outofplane.h"

#ifndef __GNUC__
#define isnan(x) ((x) != (x))
#endif

using namespace std;
//==============================================================================
void greeneval2d(const void* const* const greenPtr, const unsigned int& nGrSet,
                 const unsigned int& nugComp, const unsigned int& ntgComp,
                 const bool& ugCmplx, const bool& tgCmplx, const bool& tg0Cmplx,
                 const double& xiR, const double& xiZ, const double& Xsgn, unsigned int& r1,
                 unsigned int& r2, unsigned int& z1, unsigned int& z2, unsigned int& zs1, double* const interpr,
                 double* const interpz, bool& extrapFlag, const bool& TmatOut,
                 const double* const Coll, const unsigned int& nColl, const unsigned int& iColl,
                 const unsigned int& zPos, double* const UgrRe, double* const UgrIm,
                 double* const TgrRe, double* const TgrIm, double* const Tgr0Re,
                 double* const Tgr0Im)
//==============================================================================
{
  const unsigned int GreenFunType=*((const unsigned int*)greenPtr[0]);
  const bool calcTg0=Tgr0Re!=0;

  if (GreenFunType==1)
  {
    // RESOLVE GREEN'S FUNCTION POINTER ARRAY
    const unsigned int nzs=*((const unsigned int*)greenPtr[1]);
    const double* const zs =(const double* const)greenPtr[2];
    const unsigned int nr=*((const unsigned int*)greenPtr[3]);
    const double* const r =(const double* const)greenPtr[4];
    const unsigned int nz=*((const unsigned int*)greenPtr[5]);
    const double* const z =(const double* const)greenPtr[6];
    const double* const ugRe =(const double* const)greenPtr[7];
    const double* const ugIm =(const double* const)greenPtr[8];
    const double* const tgRe =(const double* const)greenPtr[9];
    const double* const tgIm =(const double* const)greenPtr[10];
    const double* const tg0Re =(const double* const)greenPtr[11];
    const double* const tg0Im =(const double* const)greenPtr[12];
    const bool zRel=*((const bool*)greenPtr[13]);

    const unsigned int rend=nr-1;
    const unsigned int zend=nz-1;
    const unsigned int zsend=nzs-1;
    searchClosest(Coll[zPos*nColl+iColl],zs,zsend,zs1);

    // Use actual source to receiver distance if zRel=false
    //const double xiZabs=(zRel?xiZ:xiZ+zs[zs1]);
    const double xiZabs=(zRel?xiZ:xiZ+Coll[zPos*nColl+iColl]);
    search1(xiR,r,rend,r1,r2,interpr,extrapFlag);
      if (extrapFlag && (r1==rend)) throw("Range of input argument 'x' insufficient.");

    search1(xiZabs,z,zend,z1,z2,interpz,extrapFlag);
      if (extrapFlag) throw("Range of input argument 'z' insufficient.");

    const double fac11 = interpr[0]*interpz[0];
    const double fac12 = interpr[0]*interpz[1];
    const double fac21 = interpr[1]*interpz[0];
    const double fac22 = interpr[1]*interpz[1];

    for (unsigned int iGrSet=0; iGrSet<nGrSet; iGrSet++)
    {
      // ind   = nugComp*(izs+nzs*(ir+nr*(iz+nz*iGrSet)));
      const unsigned int ind11=zs1+nzs*(r1+nr*(z1+nz*iGrSet));
      const unsigned int ind12=zs1+nzs*(r1+nr*(z2+nz*iGrSet));
      const unsigned int ind21=zs1+nzs*(r2+nr*(z1+nz*iGrSet));
      const unsigned int ind22=zs1+nzs*(r2+nr*(z2+nz*iGrSet));

      if (nugComp==1)
      {
        UgrRe[nugComp*iGrSet+0]=ugRe[nugComp*ind11+0]*fac11+ugRe[nugComp*ind12+0]*fac12+ugRe[nugComp*ind21+0]*fac21+ugRe[nugComp*ind22+0]*fac22;
        if (ugCmplx)
        {
          UgrIm[nugComp*iGrSet+0]=ugIm[nugComp*ind11+0]*fac11+ugIm[nugComp*ind12+0]*fac12+ugIm[nugComp*ind21+0]*fac21+ugIm[nugComp*ind22+0]*fac22;
        }
        if (TmatOut)
        {
          TgrRe[ntgComp*iGrSet+0]=tgRe[ntgComp*ind11+0]*fac11+tgRe[ntgComp*ind12+0]*fac12+tgRe[ntgComp*ind21+0]*fac21+tgRe[ntgComp*ind22+0]*fac22;
          TgrRe[ntgComp*iGrSet+1]=tgRe[ntgComp*ind11+1]*fac11+tgRe[ntgComp*ind12+1]*fac12+tgRe[ntgComp*ind21+1]*fac21+tgRe[ntgComp*ind22+1]*fac22;
          if (tgCmplx)
          {
            TgrIm[ntgComp*iGrSet+0]=tgIm[ntgComp*ind11+0]*fac11+tgIm[ntgComp*ind12+0]*fac12+tgIm[ntgComp*ind21+0]*fac21+tgIm[ntgComp*ind22+0]*fac22;
            TgrIm[ntgComp*iGrSet+1]=tgIm[ntgComp*ind11+1]*fac11+tgIm[ntgComp*ind12+1]*fac12+tgIm[ntgComp*ind21+1]*fac21+tgIm[ntgComp*ind22+1]*fac22;
          }
          if (calcTg0)
          {
            Tgr0Re[ntgComp*iGrSet+0]=tg0Re[ntgComp*ind11+0]*fac11+tg0Re[ntgComp*ind12+0]*fac12+tg0Re[ntgComp*ind21+0]*fac21+tg0Re[ntgComp*ind22+0]*fac22;
            Tgr0Re[ntgComp*iGrSet+1]=tg0Re[ntgComp*ind11+1]*fac11+tg0Re[ntgComp*ind12+1]*fac12+tg0Re[ntgComp*ind21+1]*fac21+tg0Re[ntgComp*ind22+1]*fac22;
            if (tg0Cmplx)
            {
              Tgr0Im[ntgComp*iGrSet+0]=tg0Im[ntgComp*ind11+0]*fac11+tg0Im[ntgComp*ind12+0]*fac12+tg0Im[ntgComp*ind21+0]*fac21+tg0Im[ntgComp*ind22+0]*fac22;
              Tgr0Im[ntgComp*iGrSet+1]=tg0Im[ntgComp*ind11+1]*fac11+tg0Im[ntgComp*ind12+1]*fac12+tg0Im[ntgComp*ind21+1]*fac21+tg0Im[ntgComp*ind22+1]*fac22;
            }
          }
        }
      }
      else if (nugComp==4)
      {
        UgrRe[nugComp*iGrSet+0]=ugRe[nugComp*ind11+0]*fac11+ugRe[nugComp*ind12+0]*fac12+ugRe[nugComp*ind21+0]*fac21+ugRe[nugComp*ind22+0]*fac22;
        UgrRe[nugComp*iGrSet+1]=ugRe[nugComp*ind11+2]*fac11+ugRe[nugComp*ind12+2]*fac12+ugRe[nugComp*ind21+2]*fac21+ugRe[nugComp*ind22+2]*fac22;
        UgrRe[nugComp*iGrSet+2]=ugRe[nugComp*ind11+1]*fac11+ugRe[nugComp*ind12+1]*fac12+ugRe[nugComp*ind21+1]*fac21+ugRe[nugComp*ind22+1]*fac22;
        UgrRe[nugComp*iGrSet+3]=ugRe[nugComp*ind11+3]*fac11+ugRe[nugComp*ind12+3]*fac12+ugRe[nugComp*ind21+3]*fac21+ugRe[nugComp*ind22+3]*fac22;
        if (ugCmplx)
        {
          UgrIm[nugComp*iGrSet+0]=ugIm[nugComp*ind11+0]*fac11+ugIm[nugComp*ind12+0]*fac12+ugIm[nugComp*ind21+0]*fac21+ugIm[nugComp*ind22+0]*fac22;
          UgrIm[nugComp*iGrSet+1]=ugIm[nugComp*ind11+2]*fac11+ugIm[nugComp*ind12+2]*fac12+ugIm[nugComp*ind21+2]*fac21+ugIm[nugComp*ind22+2]*fac22;
          UgrIm[nugComp*iGrSet+2]=ugIm[nugComp*ind11+1]*fac11+ugIm[nugComp*ind12+1]*fac12+ugIm[nugComp*ind21+1]*fac21+ugIm[nugComp*ind22+1]*fac22;
          UgrIm[nugComp*iGrSet+3]=ugIm[nugComp*ind11+3]*fac11+ugIm[nugComp*ind12+3]*fac12+ugIm[nugComp*ind21+3]*fac21+ugIm[nugComp*ind22+3]*fac22;
        }
        if (TmatOut)
        {
          TgrRe[ntgComp*iGrSet+0]=tgRe[ntgComp*ind11+0]*fac11+tgRe[ntgComp*ind12+0]*fac12+tgRe[ntgComp*ind21+0]*fac21+tgRe[ntgComp*ind22+0]*fac22;
          TgrRe[ntgComp*iGrSet+1]=tgRe[ntgComp*ind11+2]*fac11+tgRe[ntgComp*ind12+2]*fac12+tgRe[ntgComp*ind21+2]*fac21+tgRe[ntgComp*ind22+2]*fac22;
          TgrRe[ntgComp*iGrSet+2]=tgRe[ntgComp*ind11+4]*fac11+tgRe[ntgComp*ind12+4]*fac12+tgRe[ntgComp*ind21+4]*fac21+tgRe[ntgComp*ind22+4]*fac22;
          TgrRe[ntgComp*iGrSet+3]=tgRe[ntgComp*ind11+1]*fac11+tgRe[ntgComp*ind12+1]*fac12+tgRe[ntgComp*ind21+1]*fac21+tgRe[ntgComp*ind22+1]*fac22;
          TgrRe[ntgComp*iGrSet+4]=tgRe[ntgComp*ind11+3]*fac11+tgRe[ntgComp*ind12+3]*fac12+tgRe[ntgComp*ind21+3]*fac21+tgRe[ntgComp*ind22+3]*fac22;
          TgrRe[ntgComp*iGrSet+5]=tgRe[ntgComp*ind11+5]*fac11+tgRe[ntgComp*ind12+5]*fac12+tgRe[ntgComp*ind21+5]*fac21+tgRe[ntgComp*ind22+5]*fac22;
          if (tgCmplx)
          {
            TgrIm[ntgComp*iGrSet+0]=tgIm[ntgComp*ind11+0]*fac11+tgIm[ntgComp*ind12+0]*fac12+tgIm[ntgComp*ind21+0]*fac21+tgIm[ntgComp*ind22+0]*fac22;
            TgrIm[ntgComp*iGrSet+1]=tgIm[ntgComp*ind11+2]*fac11+tgIm[ntgComp*ind12+2]*fac12+tgIm[ntgComp*ind21+2]*fac21+tgIm[ntgComp*ind22+2]*fac22;
            TgrIm[ntgComp*iGrSet+2]=tgIm[ntgComp*ind11+4]*fac11+tgIm[ntgComp*ind12+4]*fac12+tgIm[ntgComp*ind21+4]*fac21+tgIm[ntgComp*ind22+4]*fac22;
            TgrIm[ntgComp*iGrSet+3]=tgIm[ntgComp*ind11+1]*fac11+tgIm[ntgComp*ind12+1]*fac12+tgIm[ntgComp*ind21+1]*fac21+tgIm[ntgComp*ind22+1]*fac22;
            TgrIm[ntgComp*iGrSet+4]=tgIm[ntgComp*ind11+3]*fac11+tgIm[ntgComp*ind12+3]*fac12+tgIm[ntgComp*ind21+3]*fac21+tgIm[ntgComp*ind22+3]*fac22;
            TgrIm[ntgComp*iGrSet+5]=tgIm[ntgComp*ind11+5]*fac11+tgIm[ntgComp*ind12+5]*fac12+tgIm[ntgComp*ind21+5]*fac21+tgIm[ntgComp*ind22+5]*fac22;
          }
          if (calcTg0)
          {
            Tgr0Re[ntgComp*iGrSet+0]=tg0Re[ntgComp*ind11+0]*fac11+tg0Re[ntgComp*ind12+0]*fac12+tg0Re[ntgComp*ind21+0]*fac21+tg0Re[ntgComp*ind22+0]*fac22;
            Tgr0Re[ntgComp*iGrSet+1]=tg0Re[ntgComp*ind11+2]*fac11+tg0Re[ntgComp*ind12+2]*fac12+tg0Re[ntgComp*ind21+2]*fac21+tg0Re[ntgComp*ind22+2]*fac22;
            Tgr0Re[ntgComp*iGrSet+2]=tg0Re[ntgComp*ind11+4]*fac11+tg0Re[ntgComp*ind12+4]*fac12+tg0Re[ntgComp*ind21+4]*fac21+tg0Re[ntgComp*ind22+4]*fac22;
            Tgr0Re[ntgComp*iGrSet+3]=tg0Re[ntgComp*ind11+1]*fac11+tg0Re[ntgComp*ind12+1]*fac12+tg0Re[ntgComp*ind21+1]*fac21+tg0Re[ntgComp*ind22+1]*fac22;
            Tgr0Re[ntgComp*iGrSet+4]=tg0Re[ntgComp*ind11+3]*fac11+tg0Re[ntgComp*ind12+3]*fac12+tg0Re[ntgComp*ind21+3]*fac21+tg0Re[ntgComp*ind22+3]*fac22;
            Tgr0Re[ntgComp*iGrSet+5]=tg0Re[ntgComp*ind11+5]*fac11+tg0Re[ntgComp*ind12+5]*fac12+tg0Re[ntgComp*ind21+5]*fac21+tg0Re[ntgComp*ind22+5]*fac22;
            if (tg0Cmplx)
            {
              Tgr0Im[ntgComp*iGrSet+0]=tg0Im[ntgComp*ind11+0]*fac11+tg0Im[ntgComp*ind12+0]*fac12+tg0Im[ntgComp*ind21+0]*fac21+tg0Im[ntgComp*ind22+0]*fac22;
              Tgr0Im[ntgComp*iGrSet+1]=tg0Im[ntgComp*ind11+2]*fac11+tg0Im[ntgComp*ind12+2]*fac12+tg0Im[ntgComp*ind21+2]*fac21+tg0Im[ntgComp*ind22+2]*fac22;
              Tgr0Im[ntgComp*iGrSet+2]=tg0Im[ntgComp*ind11+4]*fac11+tg0Im[ntgComp*ind12+4]*fac12+tg0Im[ntgComp*ind21+4]*fac21+tg0Im[ntgComp*ind22+4]*fac22;
              Tgr0Im[ntgComp*iGrSet+3]=tg0Im[ntgComp*ind11+1]*fac11+tg0Im[ntgComp*ind12+1]*fac12+tg0Im[ntgComp*ind21+1]*fac21+tg0Im[ntgComp*ind22+1]*fac22;
              Tgr0Im[ntgComp*iGrSet+4]=tg0Im[ntgComp*ind11+3]*fac11+tg0Im[ntgComp*ind12+3]*fac12+tg0Im[ntgComp*ind21+3]*fac21+tg0Im[ntgComp*ind22+3]*fac22;
              Tgr0Im[ntgComp*iGrSet+5]=tg0Im[ntgComp*ind11+5]*fac11+tg0Im[ntgComp*ind12+5]*fac12+tg0Im[ntgComp*ind21+5]*fac21+tg0Im[ntgComp*ind22+5]*fac22;
            }
          }
        }
      }
      else if (nugComp==9)
      {
        UgrRe[nugComp*iGrSet+0]=ugRe[nugComp*ind11+0]*fac11+ugRe[nugComp*ind12+0]*fac12+ugRe[nugComp*ind21+0]*fac21+ugRe[nugComp*ind22+0]*fac22; // ugxx
        UgrRe[nugComp*iGrSet+1]=ugRe[nugComp*ind11+3]*fac11+ugRe[nugComp*ind12+3]*fac12+ugRe[nugComp*ind21+3]*fac21+ugRe[nugComp*ind22+3]*fac22; // ugxy
        UgrRe[nugComp*iGrSet+2]=ugRe[nugComp*ind11+6]*fac11+ugRe[nugComp*ind12+6]*fac12+ugRe[nugComp*ind21+6]*fac21+ugRe[nugComp*ind22+6]*fac22; // ugxz
        UgrRe[nugComp*iGrSet+3]=ugRe[nugComp*ind11+1]*fac11+ugRe[nugComp*ind12+1]*fac12+ugRe[nugComp*ind21+1]*fac21+ugRe[nugComp*ind22+1]*fac22; // ugyx
        UgrRe[nugComp*iGrSet+4]=ugRe[nugComp*ind11+4]*fac11+ugRe[nugComp*ind12+4]*fac12+ugRe[nugComp*ind21+4]*fac21+ugRe[nugComp*ind22+4]*fac22; // ugyy
        UgrRe[nugComp*iGrSet+5]=ugRe[nugComp*ind11+7]*fac11+ugRe[nugComp*ind12+7]*fac12+ugRe[nugComp*ind21+7]*fac21+ugRe[nugComp*ind22+7]*fac22; // ugyz
        UgrRe[nugComp*iGrSet+6]=ugRe[nugComp*ind11+2]*fac11+ugRe[nugComp*ind12+2]*fac12+ugRe[nugComp*ind21+2]*fac21+ugRe[nugComp*ind22+2]*fac22; // ugzx
        UgrRe[nugComp*iGrSet+7]=ugRe[nugComp*ind11+5]*fac11+ugRe[nugComp*ind12+5]*fac12+ugRe[nugComp*ind21+5]*fac21+ugRe[nugComp*ind22+5]*fac22; // ugzy
        UgrRe[nugComp*iGrSet+8]=ugRe[nugComp*ind11+8]*fac11+ugRe[nugComp*ind12+8]*fac12+ugRe[nugComp*ind21+8]*fac21+ugRe[nugComp*ind22+8]*fac22; // ugzz
        if (ugCmplx)
        {
          UgrIm[nugComp*iGrSet+0]=ugIm[nugComp*ind11+0]*fac11+ugIm[nugComp*ind12+0]*fac12+ugIm[nugComp*ind21+0]*fac21+ugIm[nugComp*ind22+0]*fac22; // ugxx
          UgrIm[nugComp*iGrSet+1]=ugIm[nugComp*ind11+3]*fac11+ugIm[nugComp*ind12+3]*fac12+ugIm[nugComp*ind21+3]*fac21+ugIm[nugComp*ind22+3]*fac22; // ugxy
          UgrIm[nugComp*iGrSet+2]=ugIm[nugComp*ind11+6]*fac11+ugIm[nugComp*ind12+6]*fac12+ugIm[nugComp*ind21+6]*fac21+ugIm[nugComp*ind22+6]*fac22; // ugxz
          UgrIm[nugComp*iGrSet+3]=ugIm[nugComp*ind11+1]*fac11+ugIm[nugComp*ind12+1]*fac12+ugIm[nugComp*ind21+1]*fac21+ugIm[nugComp*ind22+1]*fac22; // ugyx
          UgrIm[nugComp*iGrSet+4]=ugIm[nugComp*ind11+4]*fac11+ugIm[nugComp*ind12+4]*fac12+ugIm[nugComp*ind21+4]*fac21+ugIm[nugComp*ind22+4]*fac22; // ugyy
          UgrIm[nugComp*iGrSet+5]=ugIm[nugComp*ind11+7]*fac11+ugIm[nugComp*ind12+7]*fac12+ugIm[nugComp*ind21+7]*fac21+ugIm[nugComp*ind22+7]*fac22; // ugyz
          UgrIm[nugComp*iGrSet+6]=ugIm[nugComp*ind11+2]*fac11+ugIm[nugComp*ind12+2]*fac12+ugIm[nugComp*ind21+2]*fac21+ugIm[nugComp*ind22+2]*fac22; // ugzx
          UgrIm[nugComp*iGrSet+7]=ugIm[nugComp*ind11+5]*fac11+ugIm[nugComp*ind12+5]*fac12+ugIm[nugComp*ind21+5]*fac21+ugIm[nugComp*ind22+5]*fac22; // ugzy
          UgrIm[nugComp*iGrSet+8]=ugIm[nugComp*ind11+8]*fac11+ugIm[nugComp*ind12+8]*fac12+ugIm[nugComp*ind21+8]*fac21+ugIm[nugComp*ind22+8]*fac22; // ugzz
        }
        if (TmatOut)
        {
          TgrRe[ntgComp*iGrSet+ 0]=tgRe[ntgComp*ind11+ 0]*fac11+tgRe[ntgComp*ind12+ 0]*fac12+tgRe[ntgComp*ind21+ 0]*fac21+tgRe[ntgComp*ind22+ 0]*fac22; // sgxxx
          TgrRe[ntgComp*iGrSet+ 1]=tgRe[ntgComp*ind11+ 3]*fac11+tgRe[ntgComp*ind12+ 3]*fac12+tgRe[ntgComp*ind21+ 3]*fac21+tgRe[ntgComp*ind22+ 3]*fac22; // sgxyy
          TgrRe[ntgComp*iGrSet+ 2]=tgRe[ntgComp*ind11+ 6]*fac11+tgRe[ntgComp*ind12+ 6]*fac12+tgRe[ntgComp*ind21+ 6]*fac21+tgRe[ntgComp*ind22+ 6]*fac22; // sgxzz
          TgrRe[ntgComp*iGrSet+ 3]=tgRe[ntgComp*ind11+ 9]*fac11+tgRe[ntgComp*ind12+ 9]*fac12+tgRe[ntgComp*ind21+ 9]*fac21+tgRe[ntgComp*ind22+ 9]*fac22; // sgxxy
          TgrRe[ntgComp*iGrSet+ 4]=tgRe[ntgComp*ind11+12]*fac11+tgRe[ntgComp*ind12+12]*fac12+tgRe[ntgComp*ind21+12]*fac21+tgRe[ntgComp*ind22+12]*fac22; // sgxyz
          TgrRe[ntgComp*iGrSet+ 5]=tgRe[ntgComp*ind11+15]*fac11+tgRe[ntgComp*ind12+15]*fac12+tgRe[ntgComp*ind21+15]*fac21+tgRe[ntgComp*ind22+15]*fac22; // sgxzx
          TgrRe[ntgComp*iGrSet+ 6]=tgRe[ntgComp*ind11+ 1]*fac11+tgRe[ntgComp*ind12+ 1]*fac12+tgRe[ntgComp*ind21+ 1]*fac21+tgRe[ntgComp*ind22+ 1]*fac22; // sgyxx
          TgrRe[ntgComp*iGrSet+ 7]=tgRe[ntgComp*ind11+ 4]*fac11+tgRe[ntgComp*ind12+ 4]*fac12+tgRe[ntgComp*ind21+ 4]*fac21+tgRe[ntgComp*ind22+ 4]*fac22; // sgyyy
          TgrRe[ntgComp*iGrSet+ 8]=tgRe[ntgComp*ind11+ 7]*fac11+tgRe[ntgComp*ind12+ 7]*fac12+tgRe[ntgComp*ind21+ 7]*fac21+tgRe[ntgComp*ind22+ 7]*fac22; // sgyzz
          TgrRe[ntgComp*iGrSet+ 9]=tgRe[ntgComp*ind11+10]*fac11+tgRe[ntgComp*ind12+10]*fac12+tgRe[ntgComp*ind21+10]*fac21+tgRe[ntgComp*ind22+10]*fac22; // sgyxy
          TgrRe[ntgComp*iGrSet+10]=tgRe[ntgComp*ind11+13]*fac11+tgRe[ntgComp*ind12+13]*fac12+tgRe[ntgComp*ind21+13]*fac21+tgRe[ntgComp*ind22+13]*fac22; // sgyyz
          TgrRe[ntgComp*iGrSet+11]=tgRe[ntgComp*ind11+16]*fac11+tgRe[ntgComp*ind12+16]*fac12+tgRe[ntgComp*ind21+16]*fac21+tgRe[ntgComp*ind22+16]*fac22; // sgyzx
          TgrRe[ntgComp*iGrSet+12]=tgRe[ntgComp*ind11+ 2]*fac11+tgRe[ntgComp*ind12+ 2]*fac12+tgRe[ntgComp*ind21+ 2]*fac21+tgRe[ntgComp*ind22+ 2]*fac22; // sgzxx
          TgrRe[ntgComp*iGrSet+13]=tgRe[ntgComp*ind11+ 5]*fac11+tgRe[ntgComp*ind12+ 5]*fac12+tgRe[ntgComp*ind21+ 5]*fac21+tgRe[ntgComp*ind22+ 5]*fac22; // sgzyy
          TgrRe[ntgComp*iGrSet+14]=tgRe[ntgComp*ind11+ 8]*fac11+tgRe[ntgComp*ind12+ 8]*fac12+tgRe[ntgComp*ind21+ 8]*fac21+tgRe[ntgComp*ind22+ 8]*fac22; // sgzzz
          TgrRe[ntgComp*iGrSet+15]=tgRe[ntgComp*ind11+11]*fac11+tgRe[ntgComp*ind12+11]*fac12+tgRe[ntgComp*ind21+11]*fac21+tgRe[ntgComp*ind22+11]*fac22; // sgzxy
          TgrRe[ntgComp*iGrSet+16]=tgRe[ntgComp*ind11+14]*fac11+tgRe[ntgComp*ind12+14]*fac12+tgRe[ntgComp*ind21+14]*fac21+tgRe[ntgComp*ind22+14]*fac22; // sgzyz
          TgrRe[ntgComp*iGrSet+17]=tgRe[ntgComp*ind11+17]*fac11+tgRe[ntgComp*ind12+17]*fac12+tgRe[ntgComp*ind21+17]*fac21+tgRe[ntgComp*ind22+17]*fac22; // sgzzx
          if (tgCmplx)
          {
            TgrIm[ntgComp*iGrSet+ 0]=tgIm[ntgComp*ind11+ 0]*fac11+tgIm[ntgComp*ind12+ 0]*fac12+tgIm[ntgComp*ind21+ 0]*fac21+tgIm[ntgComp*ind22+ 0]*fac22; // sgxxx
            TgrIm[ntgComp*iGrSet+ 1]=tgIm[ntgComp*ind11+ 3]*fac11+tgIm[ntgComp*ind12+ 3]*fac12+tgIm[ntgComp*ind21+ 3]*fac21+tgIm[ntgComp*ind22+ 3]*fac22; // sgxyy
            TgrIm[ntgComp*iGrSet+ 2]=tgIm[ntgComp*ind11+ 6]*fac11+tgIm[ntgComp*ind12+ 6]*fac12+tgIm[ntgComp*ind21+ 6]*fac21+tgIm[ntgComp*ind22+ 6]*fac22; // sgxzz
            TgrIm[ntgComp*iGrSet+ 3]=tgIm[ntgComp*ind11+ 9]*fac11+tgIm[ntgComp*ind12+ 9]*fac12+tgIm[ntgComp*ind21+ 9]*fac21+tgIm[ntgComp*ind22+ 9]*fac22; // sgxxy
            TgrIm[ntgComp*iGrSet+ 4]=tgIm[ntgComp*ind11+12]*fac11+tgIm[ntgComp*ind12+12]*fac12+tgIm[ntgComp*ind21+12]*fac21+tgIm[ntgComp*ind22+12]*fac22; // sgxyz
            TgrIm[ntgComp*iGrSet+ 5]=tgIm[ntgComp*ind11+15]*fac11+tgIm[ntgComp*ind12+15]*fac12+tgIm[ntgComp*ind21+15]*fac21+tgIm[ntgComp*ind22+15]*fac22; // sgxzx
            TgrIm[ntgComp*iGrSet+ 6]=tgIm[ntgComp*ind11+ 1]*fac11+tgIm[ntgComp*ind12+ 1]*fac12+tgIm[ntgComp*ind21+ 1]*fac21+tgIm[ntgComp*ind22+ 1]*fac22; // sgyxx
            TgrIm[ntgComp*iGrSet+ 7]=tgIm[ntgComp*ind11+ 4]*fac11+tgIm[ntgComp*ind12+ 4]*fac12+tgIm[ntgComp*ind21+ 4]*fac21+tgIm[ntgComp*ind22+ 4]*fac22; // sgyyy
            TgrIm[ntgComp*iGrSet+ 8]=tgIm[ntgComp*ind11+ 7]*fac11+tgIm[ntgComp*ind12+ 7]*fac12+tgIm[ntgComp*ind21+ 7]*fac21+tgIm[ntgComp*ind22+ 7]*fac22; // sgyzz
            TgrIm[ntgComp*iGrSet+ 9]=tgIm[ntgComp*ind11+10]*fac11+tgIm[ntgComp*ind12+10]*fac12+tgIm[ntgComp*ind21+10]*fac21+tgIm[ntgComp*ind22+10]*fac22; // sgyxy
            TgrIm[ntgComp*iGrSet+10]=tgIm[ntgComp*ind11+13]*fac11+tgIm[ntgComp*ind12+13]*fac12+tgIm[ntgComp*ind21+13]*fac21+tgIm[ntgComp*ind22+13]*fac22; // sgyyz
            TgrIm[ntgComp*iGrSet+11]=tgIm[ntgComp*ind11+16]*fac11+tgIm[ntgComp*ind12+16]*fac12+tgIm[ntgComp*ind21+16]*fac21+tgIm[ntgComp*ind22+16]*fac22; // sgyzx
            TgrIm[ntgComp*iGrSet+12]=tgIm[ntgComp*ind11+ 2]*fac11+tgIm[ntgComp*ind12+ 2]*fac12+tgIm[ntgComp*ind21+ 2]*fac21+tgIm[ntgComp*ind22+ 2]*fac22; // sgzxx
            TgrIm[ntgComp*iGrSet+13]=tgIm[ntgComp*ind11+ 5]*fac11+tgIm[ntgComp*ind12+ 5]*fac12+tgIm[ntgComp*ind21+ 5]*fac21+tgIm[ntgComp*ind22+ 5]*fac22; // sgzyy
            TgrIm[ntgComp*iGrSet+14]=tgIm[ntgComp*ind11+ 8]*fac11+tgIm[ntgComp*ind12+ 8]*fac12+tgIm[ntgComp*ind21+ 8]*fac21+tgIm[ntgComp*ind22+ 8]*fac22; // sgzzz
            TgrIm[ntgComp*iGrSet+15]=tgIm[ntgComp*ind11+11]*fac11+tgIm[ntgComp*ind12+11]*fac12+tgIm[ntgComp*ind21+11]*fac21+tgIm[ntgComp*ind22+11]*fac22; // sgzxy
            TgrIm[ntgComp*iGrSet+16]=tgIm[ntgComp*ind11+14]*fac11+tgIm[ntgComp*ind12+14]*fac12+tgIm[ntgComp*ind21+14]*fac21+tgIm[ntgComp*ind22+14]*fac22; // sgzyz
            TgrIm[ntgComp*iGrSet+17]=tgIm[ntgComp*ind11+17]*fac11+tgIm[ntgComp*ind12+17]*fac12+tgIm[ntgComp*ind21+17]*fac21+tgIm[ntgComp*ind22+17]*fac22; // sgzzx
          }
          if (calcTg0)
          {
            Tgr0Re[ntgComp*iGrSet+ 0]=tg0Re[ntgComp*ind11+ 0]*fac11+tg0Re[ntgComp*ind12+ 0]*fac12+tg0Re[ntgComp*ind21+ 0]*fac21+tg0Re[ntgComp*ind22+ 0]*fac22; // sgxxx
            Tgr0Re[ntgComp*iGrSet+ 1]=tg0Re[ntgComp*ind11+ 3]*fac11+tg0Re[ntgComp*ind12+ 3]*fac12+tg0Re[ntgComp*ind21+ 3]*fac21+tg0Re[ntgComp*ind22+ 3]*fac22; // sgxyy
            Tgr0Re[ntgComp*iGrSet+ 2]=tg0Re[ntgComp*ind11+ 6]*fac11+tg0Re[ntgComp*ind12+ 6]*fac12+tg0Re[ntgComp*ind21+ 6]*fac21+tg0Re[ntgComp*ind22+ 6]*fac22; // sgxzz
            Tgr0Re[ntgComp*iGrSet+ 3]=tg0Re[ntgComp*ind11+ 9]*fac11+tg0Re[ntgComp*ind12+ 9]*fac12+tg0Re[ntgComp*ind21+ 9]*fac21+tg0Re[ntgComp*ind22+ 9]*fac22; // sgxxy
            Tgr0Re[ntgComp*iGrSet+ 4]=tg0Re[ntgComp*ind11+12]*fac11+tg0Re[ntgComp*ind12+12]*fac12+tg0Re[ntgComp*ind21+12]*fac21+tg0Re[ntgComp*ind22+12]*fac22; // sgxyz
            Tgr0Re[ntgComp*iGrSet+ 5]=tg0Re[ntgComp*ind11+15]*fac11+tg0Re[ntgComp*ind12+15]*fac12+tg0Re[ntgComp*ind21+15]*fac21+tg0Re[ntgComp*ind22+15]*fac22; // sgxzx
            Tgr0Re[ntgComp*iGrSet+ 6]=tg0Re[ntgComp*ind11+ 1]*fac11+tg0Re[ntgComp*ind12+ 1]*fac12+tg0Re[ntgComp*ind21+ 1]*fac21+tg0Re[ntgComp*ind22+ 1]*fac22; // sgyxx
            Tgr0Re[ntgComp*iGrSet+ 7]=tg0Re[ntgComp*ind11+ 4]*fac11+tg0Re[ntgComp*ind12+ 4]*fac12+tg0Re[ntgComp*ind21+ 4]*fac21+tg0Re[ntgComp*ind22+ 4]*fac22; // sgyyy
            Tgr0Re[ntgComp*iGrSet+ 8]=tg0Re[ntgComp*ind11+ 7]*fac11+tg0Re[ntgComp*ind12+ 7]*fac12+tg0Re[ntgComp*ind21+ 7]*fac21+tg0Re[ntgComp*ind22+ 7]*fac22; // sgyzz
            Tgr0Re[ntgComp*iGrSet+ 9]=tg0Re[ntgComp*ind11+10]*fac11+tg0Re[ntgComp*ind12+10]*fac12+tg0Re[ntgComp*ind21+10]*fac21+tg0Re[ntgComp*ind22+10]*fac22; // sgyxy
            Tgr0Re[ntgComp*iGrSet+10]=tg0Re[ntgComp*ind11+13]*fac11+tg0Re[ntgComp*ind12+13]*fac12+tg0Re[ntgComp*ind21+13]*fac21+tg0Re[ntgComp*ind22+13]*fac22; // sgyyz
            Tgr0Re[ntgComp*iGrSet+11]=tg0Re[ntgComp*ind11+16]*fac11+tg0Re[ntgComp*ind12+16]*fac12+tg0Re[ntgComp*ind21+16]*fac21+tg0Re[ntgComp*ind22+16]*fac22; // sgyzx
            Tgr0Re[ntgComp*iGrSet+12]=tg0Re[ntgComp*ind11+ 2]*fac11+tg0Re[ntgComp*ind12+ 2]*fac12+tg0Re[ntgComp*ind21+ 2]*fac21+tg0Re[ntgComp*ind22+ 2]*fac22; // sgzxx
            Tgr0Re[ntgComp*iGrSet+13]=tg0Re[ntgComp*ind11+ 5]*fac11+tg0Re[ntgComp*ind12+ 5]*fac12+tg0Re[ntgComp*ind21+ 5]*fac21+tg0Re[ntgComp*ind22+ 5]*fac22; // sgzyy
            Tgr0Re[ntgComp*iGrSet+14]=tg0Re[ntgComp*ind11+ 8]*fac11+tg0Re[ntgComp*ind12+ 8]*fac12+tg0Re[ntgComp*ind21+ 8]*fac21+tg0Re[ntgComp*ind22+ 8]*fac22; // sgzzz
            Tgr0Re[ntgComp*iGrSet+15]=tg0Re[ntgComp*ind11+11]*fac11+tg0Re[ntgComp*ind12+11]*fac12+tg0Re[ntgComp*ind21+11]*fac21+tg0Re[ntgComp*ind22+11]*fac22; // sgzxy
            Tgr0Re[ntgComp*iGrSet+16]=tg0Re[ntgComp*ind11+14]*fac11+tg0Re[ntgComp*ind12+14]*fac12+tg0Re[ntgComp*ind21+14]*fac21+tg0Re[ntgComp*ind22+14]*fac22; // sgzyz
            Tgr0Re[ntgComp*iGrSet+17]=tg0Re[ntgComp*ind11+17]*fac11+tg0Re[ntgComp*ind12+17]*fac12+tg0Re[ntgComp*ind21+17]*fac21+tg0Re[ntgComp*ind22+17]*fac22; // sgzzx
            if (tg0Cmplx)
            {
              Tgr0Im[ntgComp*iGrSet+ 0]=tg0Im[ntgComp*ind11+ 0]*fac11+tg0Im[ntgComp*ind12+ 0]*fac12+tg0Im[ntgComp*ind21+ 0]*fac21+tg0Im[ntgComp*ind22+ 0]*fac22; // sgxxx
              Tgr0Im[ntgComp*iGrSet+ 1]=tg0Im[ntgComp*ind11+ 3]*fac11+tg0Im[ntgComp*ind12+ 3]*fac12+tg0Im[ntgComp*ind21+ 3]*fac21+tg0Im[ntgComp*ind22+ 3]*fac22; // sgxyy
              Tgr0Im[ntgComp*iGrSet+ 2]=tg0Im[ntgComp*ind11+ 6]*fac11+tg0Im[ntgComp*ind12+ 6]*fac12+tg0Im[ntgComp*ind21+ 6]*fac21+tg0Im[ntgComp*ind22+ 6]*fac22; // sgxzz
              Tgr0Im[ntgComp*iGrSet+ 3]=tg0Im[ntgComp*ind11+ 9]*fac11+tg0Im[ntgComp*ind12+ 9]*fac12+tg0Im[ntgComp*ind21+ 9]*fac21+tg0Im[ntgComp*ind22+ 9]*fac22; // sgxxy
              Tgr0Im[ntgComp*iGrSet+ 4]=tg0Im[ntgComp*ind11+12]*fac11+tg0Im[ntgComp*ind12+12]*fac12+tg0Im[ntgComp*ind21+12]*fac21+tg0Im[ntgComp*ind22+12]*fac22; // sgxyz
              Tgr0Im[ntgComp*iGrSet+ 5]=tg0Im[ntgComp*ind11+15]*fac11+tg0Im[ntgComp*ind12+15]*fac12+tg0Im[ntgComp*ind21+15]*fac21+tg0Im[ntgComp*ind22+15]*fac22; // sgxzx
              Tgr0Im[ntgComp*iGrSet+ 6]=tg0Im[ntgComp*ind11+ 1]*fac11+tg0Im[ntgComp*ind12+ 1]*fac12+tg0Im[ntgComp*ind21+ 1]*fac21+tg0Im[ntgComp*ind22+ 1]*fac22; // sgyxx
              Tgr0Im[ntgComp*iGrSet+ 7]=tg0Im[ntgComp*ind11+ 4]*fac11+tg0Im[ntgComp*ind12+ 4]*fac12+tg0Im[ntgComp*ind21+ 4]*fac21+tg0Im[ntgComp*ind22+ 4]*fac22; // sgyyy
              Tgr0Im[ntgComp*iGrSet+ 8]=tg0Im[ntgComp*ind11+ 7]*fac11+tg0Im[ntgComp*ind12+ 7]*fac12+tg0Im[ntgComp*ind21+ 7]*fac21+tg0Im[ntgComp*ind22+ 7]*fac22; // sgyzz
              Tgr0Im[ntgComp*iGrSet+ 9]=tg0Im[ntgComp*ind11+10]*fac11+tg0Im[ntgComp*ind12+10]*fac12+tg0Im[ntgComp*ind21+10]*fac21+tg0Im[ntgComp*ind22+10]*fac22; // sgyxy
              Tgr0Im[ntgComp*iGrSet+10]=tg0Im[ntgComp*ind11+13]*fac11+tg0Im[ntgComp*ind12+13]*fac12+tg0Im[ntgComp*ind21+13]*fac21+tg0Im[ntgComp*ind22+13]*fac22; // sgyyz
              Tgr0Im[ntgComp*iGrSet+11]=tg0Im[ntgComp*ind11+16]*fac11+tg0Im[ntgComp*ind12+16]*fac12+tg0Im[ntgComp*ind21+16]*fac21+tg0Im[ntgComp*ind22+16]*fac22; // sgyzx
              Tgr0Im[ntgComp*iGrSet+12]=tg0Im[ntgComp*ind11+ 2]*fac11+tg0Im[ntgComp*ind12+ 2]*fac12+tg0Im[ntgComp*ind21+ 2]*fac21+tg0Im[ntgComp*ind22+ 2]*fac22; // sgzxx
              Tgr0Im[ntgComp*iGrSet+13]=tg0Im[ntgComp*ind11+ 5]*fac11+tg0Im[ntgComp*ind12+ 5]*fac12+tg0Im[ntgComp*ind21+ 5]*fac21+tg0Im[ntgComp*ind22+ 5]*fac22; // sgzyy
              Tgr0Im[ntgComp*iGrSet+14]=tg0Im[ntgComp*ind11+ 8]*fac11+tg0Im[ntgComp*ind12+ 8]*fac12+tg0Im[ntgComp*ind21+ 8]*fac21+tg0Im[ntgComp*ind22+ 8]*fac22; // sgzzz
              Tgr0Im[ntgComp*iGrSet+15]=tg0Im[ntgComp*ind11+11]*fac11+tg0Im[ntgComp*ind12+11]*fac12+tg0Im[ntgComp*ind21+11]*fac21+tg0Im[ntgComp*ind22+11]*fac22; // sgzxy
              Tgr0Im[ntgComp*iGrSet+16]=tg0Im[ntgComp*ind11+14]*fac11+tg0Im[ntgComp*ind12+14]*fac12+tg0Im[ntgComp*ind21+14]*fac21+tg0Im[ntgComp*ind22+14]*fac22; // sgzyz
              Tgr0Im[ntgComp*iGrSet+17]=tg0Im[ntgComp*ind11+17]*fac11+tg0Im[ntgComp*ind12+17]*fac12+tg0Im[ntgComp*ind21+17]*fac21+tg0Im[ntgComp*ind22+17]*fac22; // sgzzx
            }
          }
        }
      }
    }
    
    // CHECK FOR NAN VALUES IN Ug,Tg and Tg0
    for (unsigned int iElt=0; iElt<nugComp*nGrSet; iElt++) if (isnan(UgrRe[iElt])) throw("Green's function has a NaN value at the integration point.");
    if (ugCmplx) for (unsigned int iElt=0; iElt<nugComp*nGrSet; iElt++) if (isnan(UgrIm[iElt])) throw("Green's function has a NaN value at the integration point.");
    if (TmatOut)
    {
      for (unsigned int iElt=0; iElt<ntgComp*nGrSet; iElt++) if (isnan(TgrRe[iElt])) throw("Green's function has a NaN value at the integration point.");
      if (tgCmplx) for (unsigned int iElt=0; iElt<ntgComp*nGrSet; iElt++) if (isnan(TgrIm[iElt])) throw("Green's function has a NaN value at the integration point.");
      if (calcTg0)
      {
        for (unsigned int iElt=0; iElt<ntgComp*nGrSet; iElt++) if (isnan(Tgr0Re[iElt])) throw("Green's function has a NaN value at the integration point.");
        if (tg0Cmplx) for (unsigned int iElt=0; iElt<ntgComp*nGrSet; iElt++) if (isnan(Tgr0Im[iElt])) throw("Green's function has a NaN value at the integration point.");
      }
    }
    
  }
  else if (GreenFunType==2) // FSGREENF (2.5D full-space solution)
  {
    // RESOLVE GREEN'S FUNCTION POINTER ARRAY
    const double Cs=*((const double*)greenPtr[1]);
    const double Cp=*((const double*)greenPtr[2]);
    const double Ds=*((const double*)greenPtr[3]);
    const double Dp=*((const double*)greenPtr[4]);
    const double rho=*((const double*)greenPtr[5]);
    const unsigned int nWave=*((const unsigned int*)greenPtr[6]);
    const unsigned int nFreq=*((const unsigned int*)greenPtr[7]);
    const double* const py=(const double* const)greenPtr[8];
    const double* const omega=(const double* const)greenPtr[9];

    // CHANGE WAVENUMBER SIGN (ug(-ky) and tg(-ky) should be integrated)
    double* const  minPy = new(nothrow) double[nWave];
    for (unsigned int iWave=0; iWave<nWave; iWave++) minPy[iWave]=-py[iWave];

    // EVALUATE ANALYTICAL SOLUTION
    const unsigned int nxRec=1;
    const unsigned int nzRec=1;
    complex<double>* const Ug = new(nothrow) complex<double>[3*3*nxRec*nWave*nzRec*nFreq];
    if (Ug==0) throw("Out of memory.");
    complex<double>* Sg = 0;
    complex<double>* Sg0 = 0;
    if (TmatOut)
    {
      Sg = new(nothrow) complex<double>[3*6*nxRec*nWave*nzRec*nFreq];
      if (Sg==0) throw("Out of memory.");
    }
    if (calcTg0)
    {
      Sg0 = new(nothrow) complex<double>[3*6*nxRec*nzRec];  // 2D static solution
      if (Sg0==0) throw("Out of memory.");
      const double py0=0.0;
      const double omega0=0.0;
      fsgreenf(Cs,Cp,Ds,Dp,rho,&xiR,&py0,&xiZ,&omega0,nxRec,1,nzRec,1,Ug,Sg0,false,calcTg0);
    }
    fsgreenf(Cs,Cp,Ds,Dp,rho,&xiR,minPy,&xiZ,omega,nxRec,nWave,nzRec,nFreq,Ug,Sg,true,TmatOut);

    // COPY RESULTS
    for (unsigned int iWave=0; iWave<nWave; iWave++)
    {
      for (unsigned int iFreq=0; iFreq<nFreq; iFreq++)
      {
        unsigned int iGrSet=nWave*iFreq+iWave;
        for (unsigned int iComp=0; iComp<9; iComp++)
        {
          UgrRe[9*iGrSet+iComp]=real(Ug[9*iGrSet+iComp]);
          UgrIm[9*iGrSet+iComp]=imag(Ug[9*iGrSet+iComp]);
        }
        if (TmatOut)
        {
          for (unsigned int iComp=0; iComp<18; iComp++)
          {
            TgrRe[18*iGrSet+iComp]=real(Sg[18*iGrSet+iComp]);
            TgrIm[18*iGrSet+iComp]=imag(Sg[18*iGrSet+iComp]);
            if (calcTg0)
            {
              Tgr0Re[18*iGrSet+iComp]=real(Sg0[iComp]);
              if (tg0Cmplx)
              {
                Tgr0Im[18*iGrSet+iComp]=imag(Sg0[iComp]);
              }
            }
          }
        }
      }
    }
    delete [] Ug;
    delete [] Sg;
    delete [] Sg0;
    delete [] minPy;
  }
  else if (GreenFunType==4) // FSGREEN2D_inplane (2D in-plane full-space solution)
  {
    // RESOLVE GREEN'S FUNCTION POINTER ARRAY
    const double Cs=*((const double*)greenPtr[1]);
    const double Cp=*((const double*)greenPtr[2]);
    const double Ds=*((const double*)greenPtr[3]);
    const double Dp=*((const double*)greenPtr[4]);
    const double rho=*((const double*)greenPtr[5]);
    const unsigned int nFreq=*((const unsigned int*)greenPtr[6]);
    const double* const omega=(const double* const)greenPtr[7];

    // EVALUATE ANALYTICAL SOLUTION
    const unsigned int nxRec=1;
    const unsigned int nzRec=1;
    complex<double>* const Ug = new(nothrow) complex<double>[2*2*nxRec*nzRec*nFreq];
    if (Ug==0) throw("Out of memory.");
    complex<double>* Sg = 0;
    complex<double>* Sg0 = 0;
    if (TmatOut)
    {
      Sg = new(nothrow) complex<double>[2*3*nxRec*nzRec*nFreq];
      if (Sg==0) throw("Out of memory.");
    }
    if (calcTg0)
    {
      Sg0 = new(nothrow) complex<double>[2*3*nxRec*nzRec];
      if (Sg0==0) throw("Out of memory.");
      const double omega0=0.0;
      const unsigned int nFreq0=1;
      fsgreen2d_inplane(Cs,Cp,Ds,Dp,rho,&xiR,&xiZ,&omega0,nxRec,nzRec,nFreq0,Ug,Sg0,false,calcTg0);
    }
    fsgreen2d_inplane(Cs,Cp,Ds,Dp,rho,&xiR,&xiZ,omega,nxRec,nzRec,nFreq,Ug,Sg,true,TmatOut);

    // COPY RESULTS
    for (unsigned int iFreq=0; iFreq<nFreq; iFreq++)
    {
      unsigned int iGrSet=iFreq;
      for (unsigned int iComp=0; iComp<4; iComp++)
      {
        UgrRe[4*iGrSet+iComp]=real(Ug[4*iGrSet+iComp]);
        UgrIm[4*iGrSet+iComp]=imag(Ug[4*iGrSet+iComp]);
      }
      if (TmatOut)
      {
        for (unsigned int iComp=0; iComp<6; iComp++)
        {
          TgrRe[6*iGrSet+iComp]=real(Sg[6*iGrSet+iComp]);
          TgrIm[6*iGrSet+iComp]=imag(Sg[6*iGrSet+iComp]);
          if (calcTg0)
          {
            Tgr0Re[6*iGrSet+iComp]=real(Sg0[iComp]);
            if (tg0Cmplx)
            {
              Tgr0Im[6*iGrSet+iComp]=imag(Sg0[iComp]);
            }
          }
        }
      }
    }
    delete [] Ug;
    delete [] Sg;
    delete [] Sg0;
  }
  else if (GreenFunType==5) // FSGREEN2D_outofplane (2D in-plane full-space solution)
  {
    // RESOLVE GREEN'S FUNCTION POINTER ARRAY
    const double Cs=*((const double*)greenPtr[1]);
    const double Ds=*((const double*)greenPtr[2]);
    const double rho=*((const double*)greenPtr[3]);
    const unsigned int nFreq=*((const unsigned int*)greenPtr[4]);
    const double* const omega=(const double* const)greenPtr[5];

    // EVALUATE ANALYTICAL SOLUTION
    const unsigned int nxRec=1;
    const unsigned int nzRec=1;
    complex<double>* const Ug = new(nothrow) complex<double>[1*1*nxRec*nzRec*nFreq];
    if (Ug==0) throw("Out of memory.");
    complex<double>* Sg = 0;
    complex<double>* Sg0 = 0;
    if (TmatOut)
    {
      Sg = new(nothrow) complex<double>[1*2*nxRec*nzRec*nFreq];
      if (Sg==0) throw("Out of memory.");
    }
    if (calcTg0)
    {
      Sg0 = new(nothrow) complex<double>[1*2*nxRec*nzRec];
      if (Sg0==0) throw("Out of memory.");
      const double omega0=0.0;
      const unsigned int nFreq0=1;
      fsgreen2d_outofplane(Cs,Ds,rho,&xiR,&xiZ,&omega0,nxRec,nzRec,nFreq0,Ug,Sg0,false,calcTg0);
    }
    fsgreen2d_outofplane(Cs,Ds,rho,&xiR,&xiZ,omega,nxRec,nzRec,nFreq,Ug,Sg,true,TmatOut);

    // COPY RESULTS
    for (unsigned int iFreq=0; iFreq<nFreq; iFreq++)
    {
      unsigned int iGrSet=iFreq;
      for (unsigned int iComp=0; iComp<1; iComp++)
      {
        UgrRe[1*iGrSet+iComp]=real(Ug[1*iGrSet+iComp]);
        UgrIm[1*iGrSet+iComp]=imag(Ug[1*iGrSet+iComp]);
      }
      if (TmatOut)
      {
        for (unsigned int iComp=0; iComp<2; iComp++)
        {
          TgrRe[2*iGrSet+iComp]=real(Sg[2*iGrSet+iComp]);
          TgrIm[2*iGrSet+iComp]=imag(Sg[2*iGrSet+iComp]);
          if (calcTg0)
          {
            Tgr0Re[2*iGrSet+iComp]=real(Sg0[iComp]);
            if (tg0Cmplx)
            {
              Tgr0Im[2*iGrSet+iComp]=imag(Sg0[iComp]);
            }
          }
        }
      }
    }
    delete [] Ug;
    delete [] Sg;
    delete [] Sg0;
  }
  else
  {
    throw("Unknown Green's function type in subroutine greeneval2d.");
  }

  // ACCOUNT FOR SYMMETRY/ANTIMETRY OF X
  for (unsigned int iGrSet=0; iGrSet<nGrSet; iGrSet++)
  {
    if (nugComp==4)
    { 
      UgrRe[4*iGrSet+1]=Xsgn*UgrRe[4*iGrSet+1];
      UgrRe[4*iGrSet+2]=Xsgn*UgrRe[4*iGrSet+2];
      if (ugCmplx)
      {
        UgrIm[4*iGrSet+1]=Xsgn*UgrIm[4*iGrSet+1];
        UgrIm[4*iGrSet+2]=Xsgn*UgrIm[4*iGrSet+2];
      }
      
     
    }
    else if (nugComp==9)
    {
      UgrRe[9*iGrSet+1]=Xsgn*UgrRe[9*iGrSet+1];
      UgrRe[9*iGrSet+2]=Xsgn*UgrRe[9*iGrSet+2];
      UgrRe[9*iGrSet+3]=Xsgn*UgrRe[9*iGrSet+3];
      UgrRe[9*iGrSet+6]=Xsgn*UgrRe[9*iGrSet+6];
      if (ugCmplx)
      {
        UgrIm[9*iGrSet+1]=Xsgn*UgrIm[9*iGrSet+1];
        UgrIm[9*iGrSet+2]=Xsgn*UgrIm[9*iGrSet+2];
        UgrIm[9*iGrSet+3]=Xsgn*UgrIm[9*iGrSet+3];
        UgrIm[9*iGrSet+6]=Xsgn*UgrIm[9*iGrSet+6];
      }
    }
    if (TmatOut)
    {
      // 2D out-of-plane
      if (ntgComp==2)
      {
        TgrRe[2*iGrSet+0]=Xsgn*TgrRe[2*iGrSet+0];
        if (tgCmplx)
        {
          TgrIm[2*iGrSet+0]=Xsgn*TgrIm[2*iGrSet+0];
        }
        if (calcTg0)
        {
          Tgr0Re[2*iGrSet+0]=Xsgn*Tgr0Re[2*iGrSet+0];
          if (tg0Cmplx)
          {
            Tgr0Im[2*iGrSet+0]=Xsgn*Tgr0Im[2*iGrSet+0];
          }
        }
      }
      // 2D in-plane
      else if (ntgComp==6)
      {
        TgrRe[6*iGrSet+0]=Xsgn*TgrRe[6*iGrSet+0];
        TgrRe[6*iGrSet+1]=Xsgn*TgrRe[6*iGrSet+1];
        TgrRe[6*iGrSet+5]=Xsgn*TgrRe[6*iGrSet+5];
        if (tgCmplx)
        {
          TgrIm[6*iGrSet+0]=Xsgn*TgrIm[6*iGrSet+0];
          TgrIm[6*iGrSet+1]=Xsgn*TgrIm[6*iGrSet+1];
          TgrIm[6*iGrSet+5]=Xsgn*TgrIm[6*iGrSet+5];
        }

        if (calcTg0)
        {
          Tgr0Re[6*iGrSet+0]=Xsgn*Tgr0Re[6*iGrSet+0];
          Tgr0Re[6*iGrSet+1]=Xsgn*Tgr0Re[6*iGrSet+1];
          Tgr0Re[6*iGrSet+5]=Xsgn*Tgr0Re[6*iGrSet+5];
          if (tg0Cmplx)
          {
            Tgr0Im[6*iGrSet+0]=Xsgn*Tgr0Im[6*iGrSet+0];
            Tgr0Im[6*iGrSet+1]=Xsgn*Tgr0Im[6*iGrSet+1];
            Tgr0Im[6*iGrSet+5]=Xsgn*Tgr0Im[6*iGrSet+5];
          }
        }
      }
      // 2.5D
      else if (ntgComp==18)
      {
        TgrRe[18*iGrSet+ 0]=Xsgn*TgrRe[18*iGrSet+ 0];
        TgrRe[18*iGrSet+ 1]=Xsgn*TgrRe[18*iGrSet+ 1];
        TgrRe[18*iGrSet+ 2]=Xsgn*TgrRe[18*iGrSet+ 2];
        TgrRe[18*iGrSet+ 4]=Xsgn*TgrRe[18*iGrSet+ 4];
        TgrRe[18*iGrSet+ 9]=Xsgn*TgrRe[18*iGrSet+ 9];
        TgrRe[18*iGrSet+11]=Xsgn*TgrRe[18*iGrSet+11];
        TgrRe[18*iGrSet+15]=Xsgn*TgrRe[18*iGrSet+15];
        TgrRe[18*iGrSet+17]=Xsgn*TgrRe[18*iGrSet+17];
        if (tgCmplx)
        {
          TgrIm[18*iGrSet+ 0]=Xsgn*TgrIm[18*iGrSet+ 0];
          TgrIm[18*iGrSet+ 1]=Xsgn*TgrIm[18*iGrSet+ 1];
          TgrIm[18*iGrSet+ 2]=Xsgn*TgrIm[18*iGrSet+ 2];
          TgrIm[18*iGrSet+ 4]=Xsgn*TgrIm[18*iGrSet+ 4];
          TgrIm[18*iGrSet+ 9]=Xsgn*TgrIm[18*iGrSet+ 9];
          TgrIm[18*iGrSet+11]=Xsgn*TgrIm[18*iGrSet+11];
          TgrIm[18*iGrSet+15]=Xsgn*TgrIm[18*iGrSet+15];
          TgrIm[18*iGrSet+17]=Xsgn*TgrIm[18*iGrSet+17];
        }
        if (calcTg0)
        {
          Tgr0Re[18*iGrSet+ 0]=Xsgn*Tgr0Re[18*iGrSet+ 0];
          Tgr0Re[18*iGrSet+ 1]=Xsgn*Tgr0Re[18*iGrSet+ 1];
          Tgr0Re[18*iGrSet+ 2]=Xsgn*Tgr0Re[18*iGrSet+ 2];
          Tgr0Re[18*iGrSet+ 4]=Xsgn*Tgr0Re[18*iGrSet+ 4];
          Tgr0Re[18*iGrSet+ 9]=Xsgn*Tgr0Re[18*iGrSet+ 9];
          Tgr0Re[18*iGrSet+11]=Xsgn*Tgr0Re[18*iGrSet+11];
          Tgr0Re[18*iGrSet+15]=Xsgn*Tgr0Re[18*iGrSet+15];
          Tgr0Re[18*iGrSet+17]=Xsgn*Tgr0Re[18*iGrSet+17];
          if (tg0Cmplx)
          {
            Tgr0Im[18*iGrSet+ 0]=Xsgn*Tgr0Im[18*iGrSet+ 0];
            Tgr0Im[18*iGrSet+ 1]=Xsgn*Tgr0Im[18*iGrSet+ 1];
            Tgr0Im[18*iGrSet+ 2]=Xsgn*Tgr0Im[18*iGrSet+ 2];
            Tgr0Im[18*iGrSet+ 4]=Xsgn*Tgr0Im[18*iGrSet+ 4];
            Tgr0Im[18*iGrSet+ 9]=Xsgn*Tgr0Im[18*iGrSet+ 9];
            Tgr0Im[18*iGrSet+11]=Xsgn*Tgr0Im[18*iGrSet+11];
            Tgr0Im[18*iGrSet+15]=Xsgn*Tgr0Im[18*iGrSet+15];
            Tgr0Im[18*iGrSet+17]=Xsgn*Tgr0Im[18*iGrSet+17];
          }
        }
      }
    }
  }
}
