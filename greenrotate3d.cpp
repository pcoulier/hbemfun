#include <math.h>
#include "mex.h"
//==============================================================================
void greenrotate3d(const double* const normal,const unsigned int& iXi,
                   const double& xiTheta,const unsigned int& nGrSet,const bool& ugCmplx,
                   const bool& tgCmplx,const bool& tg0Cmplx,
                   const double* const UgrRe,const double* const UgrIm,
                   const double* const TgrRe, const double* const TgrIm,
                   const double* const Tgr0Re, const double* const Tgr0Im,
                   double* const UXiRe, double* const UXiIm,
                   double* const TXiRe, double* const TXiIm, double* const TXi0Re,
                   double* const TXi0Im, const bool& UmatOut,const bool& TmatOut)
/*
 * Rotate the 3D Green's displacement and traction functions.
 */
//==============================================================================
{
  const bool calcTg0=Tgr0Re!=0;

  const double cosxi=cos(xiTheta);
  const double sinxi=sin(xiTheta);
  const double cosxi2=cosxi*cosxi;
  const double sinxi2=sinxi*sinxi;
  const double cosxi3=cosxi2*cosxi;
  const double sinxi3=sinxi2*sinxi;

  // ug[0] = ugxr
  // ug[1] = ugxz
  // ug[2] = ugyt
  // ug[3] = ugzr
  // ug[4] = ugzz

  if (UmatOut)
  {	
  // if (iXi==0)
	// { 
	// mexPrintf("test...\n");
	// }

  for (unsigned int iGrSet=0; iGrSet<nGrSet; iGrSet++)
  {
    UXiRe[9*iGrSet+0]=cosxi2*UgrRe[5*iGrSet+0]+sinxi2*UgrRe[5*iGrSet+2];// ugxx
    UXiRe[9*iGrSet+1]=cosxi*sinxi*(UgrRe[5*iGrSet+0]-UgrRe[5*iGrSet+2]);// ugxy
    UXiRe[9*iGrSet+2]=cosxi*UgrRe[5*iGrSet+1];                          // ugxz
    UXiRe[9*iGrSet+3]=UXiRe[9*iGrSet+1];                                // ugyx
    UXiRe[9*iGrSet+4]=sinxi2*UgrRe[5*iGrSet+0]+cosxi2*UgrRe[5*iGrSet+2];// ugyy
    UXiRe[9*iGrSet+5]=sinxi*UgrRe[5*iGrSet+1];                          // ugyz
    UXiRe[9*iGrSet+6]=cosxi*UgrRe[5*iGrSet+3];                          // ugzx
    UXiRe[9*iGrSet+7]=sinxi*UgrRe[5*iGrSet+3];                          // ugzy
    UXiRe[9*iGrSet+8]=UgrRe[5*iGrSet+4];                                // ugzz
    if (ugCmplx)
    {
      UXiIm[9*iGrSet+0]=cosxi2*UgrIm[5*iGrSet+0]+sinxi2*UgrIm[5*iGrSet+2];// ugxx
      UXiIm[9*iGrSet+1]=cosxi*sinxi*(UgrIm[5*iGrSet+0]-UgrIm[5*iGrSet+2]);// ugxy
      UXiIm[9*iGrSet+2]=cosxi*UgrIm[5*iGrSet+1];                          // ugxz
      UXiIm[9*iGrSet+3]=UXiIm[9*iGrSet+1];                                // ugyx
      UXiIm[9*iGrSet+4]=sinxi2*UgrIm[5*iGrSet+0]+cosxi2*UgrIm[5*iGrSet+2];// ugyy
      UXiIm[9*iGrSet+5]=sinxi*UgrIm[5*iGrSet+1];                          // ugyz
      UXiIm[9*iGrSet+6]=cosxi*UgrIm[5*iGrSet+3];                          // ugzx
      UXiIm[9*iGrSet+7]=sinxi*UgrIm[5*iGrSet+3];                          // ugzy
      UXiIm[9*iGrSet+8]=UgrIm[5*iGrSet+4];                                // ugzz
    }
  }
 }
	
	
  // ROTATE THE GREEN'S TRACTION VECTOR
  if (TmatOut)
  {
     for (unsigned int iGrSet=0; iGrSet<nGrSet; iGrSet++)
     {
      const double tgxxxRe=cosxi3*TgrRe[10*iGrSet+0]+cosxi*sinxi2*TgrRe[10*iGrSet+1]+2.0*cosxi*sinxi2*TgrRe[10*iGrSet+4];
      const double tgxyyRe=cosxi*sinxi2*TgrRe[10*iGrSet+0]+cosxi3*TgrRe[10*iGrSet+1]-2.0*cosxi*sinxi2*TgrRe[10*iGrSet+4];
      const double tgxzzRe=cosxi*TgrRe[10*iGrSet+2];
      const double tgxxyRe=cosxi2*sinxi*TgrRe[10*iGrSet+0]-cosxi2*sinxi*TgrRe[10*iGrSet+1]-cosxi2*sinxi*TgrRe[10*iGrSet+4]+sinxi3*TgrRe[10*iGrSet+4];
      const double tgxyzRe=sinxi*cosxi*TgrRe[10*iGrSet+3]-cosxi*sinxi*TgrRe[10*iGrSet+5];
      const double tgxzxRe=cosxi2*TgrRe[10*iGrSet+3]+sinxi2*TgrRe[10*iGrSet+5];

      const double tgyxxRe=sinxi*cosxi2*TgrRe[10*iGrSet+0]-2.0*cosxi2*sinxi*TgrRe[10*iGrSet+4]+sinxi3*TgrRe[10*iGrSet+1];
      const double tgyyyRe=sinxi3*TgrRe[10*iGrSet+0]+2.0*cosxi2*sinxi*TgrRe[10*iGrSet+4]+sinxi*cosxi2*TgrRe[10*iGrSet+1];
      const double tgyzzRe=sinxi*TgrRe[10*iGrSet+2];
      const double tgyxyRe=sinxi2*cosxi*TgrRe[10*iGrSet+0]-cosxi*sinxi2*TgrRe[10*iGrSet+4]+cosxi3*TgrRe[10*iGrSet+4]-sinxi2*cosxi*TgrRe[10*iGrSet+1];
      const double tgyyzRe=sinxi2*TgrRe[10*iGrSet+3]+cosxi2*TgrRe[10*iGrSet+5];
      const double tgyzxRe=sinxi*cosxi*TgrRe[10*iGrSet+3]-cosxi*sinxi*TgrRe[10*iGrSet+5];

      const double tgzxxRe=cosxi2*TgrRe[10*iGrSet+6]+sinxi2*TgrRe[10*iGrSet+7];
      const double tgzyyRe=sinxi2*TgrRe[10*iGrSet+6]+cosxi2*TgrRe[10*iGrSet+7];
      const double tgzzzRe=TgrRe[10*iGrSet+8];
      const double tgzxyRe=cosxi*sinxi*(TgrRe[10*iGrSet+6]-TgrRe[10*iGrSet+7]);
      const double tgzyzRe=sinxi*TgrRe[10*iGrSet+9];
      const double tgzzxRe=cosxi*TgrRe[10*iGrSet+9];

      // Project traction vector on element normal
      TXiRe[9*iGrSet+0]=tgxxxRe*normal[3*iXi+0]+tgxxyRe*normal[3*iXi+1]+tgxzxRe*normal[3*iXi+2];  // txx
      TXiRe[9*iGrSet+1]=tgxxyRe*normal[3*iXi+0]+tgxyyRe*normal[3*iXi+1]+tgxyzRe*normal[3*iXi+2];  // txy
      TXiRe[9*iGrSet+2]=tgxzxRe*normal[3*iXi+0]+tgxyzRe*normal[3*iXi+1]+tgxzzRe*normal[3*iXi+2];  // txz
      TXiRe[9*iGrSet+3]=tgyxxRe*normal[3*iXi+0]+tgyxyRe*normal[3*iXi+1]+tgyzxRe*normal[3*iXi+2];  // tyx
      TXiRe[9*iGrSet+4]=tgyxyRe*normal[3*iXi+0]+tgyyyRe*normal[3*iXi+1]+tgyyzRe*normal[3*iXi+2];  // tyy
      TXiRe[9*iGrSet+5]=tgyzxRe*normal[3*iXi+0]+tgyyzRe*normal[3*iXi+1]+tgyzzRe*normal[3*iXi+2];  // tyz
      TXiRe[9*iGrSet+6]=tgzxxRe*normal[3*iXi+0]+tgzxyRe*normal[3*iXi+1]+tgzzxRe*normal[3*iXi+2];  // tzx
      TXiRe[9*iGrSet+7]=tgzxyRe*normal[3*iXi+0]+tgzyyRe*normal[3*iXi+1]+tgzyzRe*normal[3*iXi+2];  // tzy
      TXiRe[9*iGrSet+8]=tgzzxRe*normal[3*iXi+0]+tgzyzRe*normal[3*iXi+1]+tgzzzRe*normal[3*iXi+2];  // tzz

      if (tgCmplx)
      {
        const double tgxxxIm=cosxi3*TgrIm[10*iGrSet+0]+cosxi*sinxi2*TgrIm[10*iGrSet+1]+2.0*cosxi*sinxi2*TgrIm[10*iGrSet+4];
        const double tgxyyIm=cosxi*sinxi2*TgrIm[10*iGrSet+0]+cosxi3*TgrIm[10*iGrSet+1]-2.0*cosxi*sinxi2*TgrIm[10*iGrSet+4];
        const double tgxzzIm=cosxi*TgrIm[10*iGrSet+2];
        const double tgxxyIm=cosxi2*sinxi*TgrIm[10*iGrSet+0]-cosxi2*sinxi*TgrIm[10*iGrSet+1]-cosxi2*sinxi*TgrIm[10*iGrSet+4]+sinxi3*TgrIm[10*iGrSet+4];
        const double tgxyzIm=sinxi*cosxi*TgrIm[10*iGrSet+3]-cosxi*sinxi*TgrIm[10*iGrSet+5];
        const double tgxzxIm=cosxi2*TgrIm[10*iGrSet+3]+sinxi2*TgrIm[10*iGrSet+5];

        const double tgyxxIm=sinxi*cosxi2*TgrIm[10*iGrSet+0]-2.0*cosxi2*sinxi*TgrIm[10*iGrSet+4]+sinxi3*TgrIm[10*iGrSet+1];
        const double tgyyyIm=sinxi3*TgrIm[10*iGrSet+0]+2.0*cosxi2*sinxi*TgrIm[10*iGrSet+4]+sinxi*cosxi2*TgrIm[10*iGrSet+1];
        const double tgyzzIm=sinxi*TgrIm[10*iGrSet+2];
        const double tgyxyIm=sinxi2*cosxi*TgrIm[10*iGrSet+0]-cosxi*sinxi2*TgrIm[10*iGrSet+4]+cosxi3*TgrIm[10*iGrSet+4]-sinxi2*cosxi*TgrIm[10*iGrSet+1];
        const double tgyyzIm=sinxi2*TgrIm[10*iGrSet+3]+cosxi2*TgrIm[10*iGrSet+5];
        const double tgyzxIm=sinxi*cosxi*TgrIm[10*iGrSet+3]-cosxi*sinxi*TgrIm[10*iGrSet+5];

        const double tgzxxIm=cosxi2*TgrIm[10*iGrSet+6]+sinxi2*TgrIm[10*iGrSet+7];
        const double tgzyyIm=sinxi2*TgrIm[10*iGrSet+6]+cosxi2*TgrIm[10*iGrSet+7];
        const double tgzzzIm=TgrIm[10*iGrSet+8];
        const double tgzxyIm=cosxi*sinxi*(TgrIm[10*iGrSet+6]-TgrIm[10*iGrSet+7]);
        const double tgzyzIm=sinxi*TgrIm[10*iGrSet+9];
        const double tgzzxIm=cosxi*TgrIm[10*iGrSet+9];

        // Project traction vector on element normal
        TXiIm[9*iGrSet+0]=tgxxxIm*normal[3*iXi+0]+tgxxyIm*normal[3*iXi+1]+tgxzxIm*normal[3*iXi+2];
        TXiIm[9*iGrSet+1]=tgxxyIm*normal[3*iXi+0]+tgxyyIm*normal[3*iXi+1]+tgxyzIm*normal[3*iXi+2];
        TXiIm[9*iGrSet+2]=tgxzxIm*normal[3*iXi+0]+tgxyzIm*normal[3*iXi+1]+tgxzzIm*normal[3*iXi+2];
        TXiIm[9*iGrSet+3]=tgyxxIm*normal[3*iXi+0]+tgyxyIm*normal[3*iXi+1]+tgyzxIm*normal[3*iXi+2];
        TXiIm[9*iGrSet+4]=tgyxyIm*normal[3*iXi+0]+tgyyyIm*normal[3*iXi+1]+tgyyzIm*normal[3*iXi+2];
        TXiIm[9*iGrSet+5]=tgyzxIm*normal[3*iXi+0]+tgyyzIm*normal[3*iXi+1]+tgyzzIm*normal[3*iXi+2];
        TXiIm[9*iGrSet+6]=tgzxxIm*normal[3*iXi+0]+tgzxyIm*normal[3*iXi+1]+tgzzxIm*normal[3*iXi+2];
        TXiIm[9*iGrSet+7]=tgzxyIm*normal[3*iXi+0]+tgzyyIm*normal[3*iXi+1]+tgzyzIm*normal[3*iXi+2];
        TXiIm[9*iGrSet+8]=tgzzxIm*normal[3*iXi+0]+tgzyzIm*normal[3*iXi+1]+tgzzzIm*normal[3*iXi+2];
      }


      // SINGULAR PART OF GREEN'S FUNCTION
      if (calcTg0)
      {
        const double tgxxx0Re=cosxi3*Tgr0Re[10*iGrSet+0]+cosxi*sinxi2*Tgr0Re[10*iGrSet+1]+2.0*cosxi*sinxi2*Tgr0Re[10*iGrSet+4];
        const double tgxyy0Re=cosxi*sinxi2*Tgr0Re[10*iGrSet+0]+cosxi3*Tgr0Re[10*iGrSet+1]-2.0*cosxi*sinxi2*Tgr0Re[10*iGrSet+4];
        const double tgxzz0Re=cosxi*Tgr0Re[10*iGrSet+2];
        const double tgxxy0Re=cosxi2*sinxi*Tgr0Re[10*iGrSet+0]-cosxi2*sinxi*Tgr0Re[10*iGrSet+1]-cosxi2*sinxi*Tgr0Re[10*iGrSet+4]+sinxi3*Tgr0Re[10*iGrSet+4];
        const double tgxyz0Re=sinxi*cosxi*Tgr0Re[10*iGrSet+3]-cosxi*sinxi*Tgr0Re[10*iGrSet+5];
        const double tgxzx0Re=cosxi2*Tgr0Re[10*iGrSet+3]+sinxi2*Tgr0Re[10*iGrSet+5];

        const double tgyxx0Re=sinxi*cosxi2*Tgr0Re[10*iGrSet+0]-2.0*cosxi2*sinxi*Tgr0Re[10*iGrSet+4]+sinxi3*Tgr0Re[10*iGrSet+1];
        const double tgyyy0Re=sinxi3*Tgr0Re[10*iGrSet+0]+2.0*cosxi2*sinxi*Tgr0Re[10*iGrSet+4]+sinxi*cosxi2*Tgr0Re[10*iGrSet+1];
        const double tgyzz0Re=sinxi*Tgr0Re[10*iGrSet+2];
        const double tgyxy0Re=sinxi2*cosxi*Tgr0Re[10*iGrSet+0]-cosxi*sinxi2*Tgr0Re[10*iGrSet+4]+cosxi3*Tgr0Re[10*iGrSet+4]-sinxi2*cosxi*Tgr0Re[10*iGrSet+1];
        const double tgyyz0Re=sinxi2*Tgr0Re[10*iGrSet+3]+cosxi2*Tgr0Re[10*iGrSet+5];
        const double tgyzx0Re=sinxi*cosxi*Tgr0Re[10*iGrSet+3]-cosxi*sinxi*Tgr0Re[10*iGrSet+5];

        const double tgzxx0Re=cosxi2*Tgr0Re[10*iGrSet+6]+sinxi2*Tgr0Re[10*iGrSet+7];
        const double tgzyy0Re=sinxi2*Tgr0Re[10*iGrSet+6]+cosxi2*Tgr0Re[10*iGrSet+7];
        const double tgzzz0Re=Tgr0Re[10*iGrSet+8];
        const double tgzxy0Re=cosxi*sinxi*(Tgr0Re[10*iGrSet+6]-Tgr0Re[10*iGrSet+7]);
        const double tgzyz0Re=sinxi*Tgr0Re[10*iGrSet+9];
        const double tgzzx0Re=cosxi*Tgr0Re[10*iGrSet+9];

        TXi0Re[9*iGrSet+0]=  tgxxx0Re*normal[3*iXi+0]+tgxxy0Re*normal[3*iXi+1]+tgxzx0Re*normal[3*iXi+2];  // txx
        TXi0Re[9*iGrSet+1]=  tgxxy0Re*normal[3*iXi+0]+tgxyy0Re*normal[3*iXi+1]+tgxyz0Re*normal[3*iXi+2];  // txy
        TXi0Re[9*iGrSet+2]=  tgxzx0Re*normal[3*iXi+0]+tgxyz0Re*normal[3*iXi+1]+tgxzz0Re*normal[3*iXi+2];  // txz
        TXi0Re[9*iGrSet+3]=  tgyxx0Re*normal[3*iXi+0]+tgyxy0Re*normal[3*iXi+1]+tgyzx0Re*normal[3*iXi+2];  // tyx
        TXi0Re[9*iGrSet+4]=  tgyxy0Re*normal[3*iXi+0]+tgyyy0Re*normal[3*iXi+1]+tgyyz0Re*normal[3*iXi+2];  // tyy
        TXi0Re[9*iGrSet+5]=  tgyzx0Re*normal[3*iXi+0]+tgyyz0Re*normal[3*iXi+1]+tgyzz0Re*normal[3*iXi+2];  // tyz
        TXi0Re[9*iGrSet+6]=  tgzxx0Re*normal[3*iXi+0]+tgzxy0Re*normal[3*iXi+1]+tgzzx0Re*normal[3*iXi+2];  // tzx
        TXi0Re[9*iGrSet+7]=  tgzxy0Re*normal[3*iXi+0]+tgzyy0Re*normal[3*iXi+1]+tgzyz0Re*normal[3*iXi+2];  // tzy
        TXi0Re[9*iGrSet+8]=  tgzzx0Re*normal[3*iXi+0]+tgzyz0Re*normal[3*iXi+1]+tgzzz0Re*normal[3*iXi+2];  // tzz
        if (tg0Cmplx)
        {
          const double tgxxx0Im=cosxi3*Tgr0Im[10*iGrSet+0]+cosxi*sinxi2*Tgr0Im[10*iGrSet+1]+2.0*cosxi*sinxi2*Tgr0Im[10*iGrSet+4];
          const double tgxyy0Im=cosxi*sinxi2*Tgr0Im[10*iGrSet+0]+cosxi3*Tgr0Im[10*iGrSet+1]-2.0*cosxi*sinxi2*Tgr0Im[10*iGrSet+4];
          const double tgxzz0Im=cosxi*Tgr0Im[10*iGrSet+2];
          const double tgxxy0Im=cosxi2*sinxi*Tgr0Im[10*iGrSet+0]-cosxi2*sinxi*Tgr0Im[10*iGrSet+1]-cosxi2*sinxi*Tgr0Im[10*iGrSet+4]+sinxi3*Tgr0Im[10*iGrSet+4];
          const double tgxyz0Im=sinxi*cosxi*Tgr0Im[10*iGrSet+3]-cosxi*sinxi*Tgr0Im[10*iGrSet+5];
          const double tgxzx0Im=cosxi2*Tgr0Im[10*iGrSet+3]+sinxi2*Tgr0Im[10*iGrSet+5];

          const double tgyxx0Im=sinxi*cosxi2*Tgr0Im[10*iGrSet+0]-2.0*cosxi2*sinxi*Tgr0Im[10*iGrSet+4]+sinxi3*Tgr0Im[10*iGrSet+1];
          const double tgyyy0Im=sinxi3*Tgr0Im[10*iGrSet+0]+2.0*cosxi2*sinxi*Tgr0Im[10*iGrSet+4]+sinxi*cosxi2*Tgr0Im[10*iGrSet+1];
          const double tgyzz0Im=sinxi*Tgr0Im[10*iGrSet+2];
          const double tgyxy0Im=sinxi2*cosxi*Tgr0Im[10*iGrSet+0]-cosxi*sinxi2*Tgr0Im[10*iGrSet+4]+cosxi3*Tgr0Im[10*iGrSet+4]-sinxi2*cosxi*Tgr0Im[10*iGrSet+1];
          const double tgyyz0Im=sinxi2*Tgr0Im[10*iGrSet+3]+cosxi2*Tgr0Im[10*iGrSet+5];
          const double tgyzx0Im=sinxi*cosxi*Tgr0Im[10*iGrSet+3]-cosxi*sinxi*Tgr0Im[10*iGrSet+5];

          const double tgzxx0Im=cosxi2*Tgr0Im[10*iGrSet+6]+sinxi2*Tgr0Im[10*iGrSet+7];
          const double tgzyy0Im=sinxi2*Tgr0Im[10*iGrSet+6]+cosxi2*Tgr0Im[10*iGrSet+7];
          const double tgzzz0Im=Tgr0Im[10*iGrSet+8];
          const double tgzxy0Im=cosxi*sinxi*(Tgr0Im[10*iGrSet+6]-Tgr0Im[10*iGrSet+7]);
          const double tgzyz0Im=sinxi*Tgr0Im[10*iGrSet+9];
          const double tgzzx0Im=cosxi*Tgr0Im[10*iGrSet+9];

          // project singular part
          TXi0Im[9*iGrSet+0]=tgxxx0Im*normal[3*iXi+0]+tgxxy0Im*normal[3*iXi+1]+tgxzx0Im*normal[3*iXi+2];
          TXi0Im[9*iGrSet+1]=tgxxy0Im*normal[3*iXi+0]+tgxyy0Im*normal[3*iXi+1]+tgxyz0Im*normal[3*iXi+2];
          TXi0Im[9*iGrSet+2]=tgxzx0Im*normal[3*iXi+0]+tgxyz0Im*normal[3*iXi+1]+tgxzz0Im*normal[3*iXi+2];
          TXi0Im[9*iGrSet+3]=tgyxx0Im*normal[3*iXi+0]+tgyxy0Im*normal[3*iXi+1]+tgyzx0Im*normal[3*iXi+2];
          TXi0Im[9*iGrSet+4]=tgyxy0Im*normal[3*iXi+0]+tgyyy0Im*normal[3*iXi+1]+tgyyz0Im*normal[3*iXi+2];
          TXi0Im[9*iGrSet+5]=tgyzx0Im*normal[3*iXi+0]+tgyyz0Im*normal[3*iXi+1]+tgyzz0Im*normal[3*iXi+2];
          TXi0Im[9*iGrSet+6]=tgzxx0Im*normal[3*iXi+0]+tgzxy0Im*normal[3*iXi+1]+tgzzx0Im*normal[3*iXi+2];
          TXi0Im[9*iGrSet+7]=tgzxy0Im*normal[3*iXi+0]+tgzyy0Im*normal[3*iXi+1]+tgzyz0Im*normal[3*iXi+2];
          TXi0Im[9*iGrSet+8]=tgzzx0Im*normal[3*iXi+0]+tgzyz0Im*normal[3*iXi+1]+tgzzz0Im*normal[3*iXi+2];
        }
        else
        {
          TXi0Im[9*iGrSet+0]=0.0;
          TXi0Im[9*iGrSet+1]=0.0;
          TXi0Im[9*iGrSet+2]=0.0;
          TXi0Im[9*iGrSet+3]=0.0;
          TXi0Im[9*iGrSet+4]=0.0;
          TXi0Im[9*iGrSet+5]=0.0;
          TXi0Im[9*iGrSet+6]=0.0;
          TXi0Im[9*iGrSet+7]=0.0;
          TXi0Im[9*iGrSet+8]=0.0;
        }
      }
    }
  }
}
