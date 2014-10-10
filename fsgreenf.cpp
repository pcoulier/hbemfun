/* fsgreenf.cpp
 *
 */

#include <math.h>
#include <complex>
#include "besselh.h"
using namespace std;

/******************************************************************************/
const complex<double> i(0.0,1.0);
/******************************************************************************/
inline complex<double> sqr(const complex<double>& a)
{
  return a*a;
}
/******************************************************************************/
inline double sqr(const double& a)
{
  return a*a;
}
/******************************************************************************/
inline double sign(const double& a)
{
  if (a==0.0) return 0.0;
  else return (a>0.0 ? 1.0 : -1.0);
}
/******************************************************************************/

void fsgreenf(const double Cs, const double Cp,
              const double Ds, const double Dp, const double rho,
              const double* const x,
              const double* const py,
              const double* const z,
              const double* const omega,
              const int nxRec, const int nWave,
              const int nzRec, const int nFreq,
              complex<double>* const Ug, complex<double>* const Sg,
              const bool calcUg, const bool calcSg)
{
  for (int iFreq=0; iFreq<nFreq; iFreq++)
  {
    const complex<double> mu=rho*sqr(Cs)*(1.0+sign(omega[iFreq])*2.0*i*Ds);
    const complex<double> M=rho*sqr(Cp)*(1.0+sign(omega[iFreq])*2.0*i*Dp);

    const complex<double> lambda=M-2.0*mu;
    const complex<double> alpha=sqrt(M/rho);
    const complex<double> beta=sqrt(mu/rho);
    const complex<double> kp=omega[iFreq]/alpha;
    const complex<double> ks=omega[iFreq]/beta;

    const double pi=3.141592653589793;

    for (int iWave=0;iWave<nWave;iWave++)
    {
      double ky;
      if (omega[iFreq]==0) ky=py[iWave]; else ky=omega[iFreq]*py[iWave];

      complex<double> ka=sqrt(sqr(kp)-sqr(ky));
      complex<double> kb=sqrt(sqr(ks)-sqr(ky));
      if (imag(ka)>0) ka=-ka;
      if (imag(kb)>0) kb=-kb;

      for (int ixRec=0;ixRec<nxRec;ixRec++)
      {
        for (int izRec=0;izRec<nzRec;izRec++)
        {
          const double r=sqrt(sqr(x[ixRec])+sqr(z[izRec]));
          const double gx=x[ixRec]/r;
          const double gz=z[izRec]/r;
          int ind=(ixRec+nxRec*(iWave+nWave*(izRec+nzRec*iFreq)));
          if (omega[iFreq]==0)
          {
            if (ky==0)  // 2D static solution
            {
              complex<double> nu=(M-2.0*mu)/(2.0*(M-mu));
              if (calcUg)
              {
                complex<double> Au= 1.0/(8.0*pi*mu*(1.0-nu));
                complex<double> nu3=3.0-4.0*nu;
                const double logr=log(r);
                Ug[9*ind+0]=Au*(gx*gx-nu3*logr);   //ugxx
                Ug[9*ind+1]=0.0;                   //ugxy
                Ug[9*ind+2]=Au*(gz*gx);            //ugxz
                Ug[9*ind+3]=0.0;                   //ugyx
                Ug[9*ind+4]=-1.0/(2.0*pi*mu)*logr; //ugyy
                Ug[9*ind+5]=0.0;                   //ugyz
                Ug[9*ind+6]=Au*(gz*gx);            //ugzx
                Ug[9*ind+7]=0.0;                   //ugzy
                Ug[9*ind+8]=Au*(gz*gz-nu3*logr);   //ugzz
              }
              if (calcSg)
              {
                complex<double> As=-1.0/(4.0*pi*(1.0-nu)*r);
                complex<double> nu2=1.0-2.0*nu;
                Sg[18*ind+0]=As*(2.0*sqr(gx)*gx+nu2*gx);                //sgxxx
                Sg[18*ind+2]=As*(2.0*gx*sqr(gz)-nu2*gx);                //sgxzz
                Sg[18*ind+3]=0.0;                                       //sgxxy
                Sg[18*ind+4]=0.0;                                       //sgxyz
                Sg[18*ind+5]=As*(2.0*sqr(gx)*gz+nu2*gz);                //sgxzx
                Sg[18*ind+6]=0.0;                                       //sgyxx
                Sg[18*ind+7]=0.0;                                       //sgyyy
                Sg[18*ind+8]=0.0;                                       //sgyzz
                Sg[18*ind+9]=-gx/(2.0*pi*r);                            //sgyxy
                Sg[18*ind+10]=-gz/(2.0*pi*r);                           //sgyyz
                Sg[18*ind+11]=0.0;                                      //sgyzx
                Sg[18*ind+12]=As*(2.0*gz*sqr(gx)-nu2*gz);               //sgzxx
                Sg[18*ind+14]=As*(2.0*sqr(gz)*gz+nu2*gz);               //sgzzz
                Sg[18*ind+15]=0.0;                                      //sgzxy
                Sg[18*ind+16]=0.0;                                      //sgzyz
                Sg[18*ind+17]=As*(2.0*sqr(gz)*gx+nu2*gx);               //sgzzx

                // Plane strain condition: syy=nu*(sxx+szz)
                Sg[18*ind+1]=nu*(Sg[18*ind+0]+Sg[18*ind+2]);    //sgxyy
                Sg[18*ind+13]=nu*(Sg[18*ind+12]+Sg[18*ind+14]); //sgzyy
              }
            }
            else  // 2.5D static solution
            {
              // Hankel functions
              
              double kySign=sign(ky);
              ky=abs(ky);
              
              double zr=0;
              double zi=-ky*r;
              double fnu=0.0;
              int HankelM=2;
              int HankelN=3;
              int kode=1;
              double* const cyr = new(nothrow) double[3];
              if (cyr==0) throw("Out of memory.");
              double* const cyi = new(nothrow) double[3];
              if (cyi==0) throw("Out of memory.");
              int nz;
              int ierr;
              zbesh_(&zr,&zi,&fnu,&kode,&HankelM,&HankelN,cyr,cyi,&nz,&ierr);
              complex<double> H0=cyr[0]+i*cyi[0];
              complex<double> H1=cyr[1]+i*cyi[1];
              complex<double> H2=cyr[2]+i*cyi[2];
              delete [] cyr;
              delete [] cyi;

              if (calcUg)
              {
                Ug[9*ind+0]=-i/(8.0*M*mu)*((M+mu)*H0+gx*gx*i*sqrt(ky*ky)*(-M+mu)*r*H1);     //ugxx
                Ug[9*ind+1]=kySign*(gx*ky*(M-mu)*r*H0)/(8.0*M*mu);                          //ugxy
                Ug[9*ind+2]=(-gx*gz*sqrt(ky*ky)*(M-mu)*r*H1)/(8.0*M*mu);                    //ugxz
                Ug[9*ind+3]=kySign*(gx*ky*(M-mu)*r*H0)/(8.0*M*mu);                          //ugyx
                Ug[9*ind+4]=-i/(8.0*M*mu)*(2.0*M*H0+(i*sqrt(ky*ky)*(M-mu)*r*H1));           //ugyy
                Ug[9*ind+5]=kySign*(gz*ky*(M-mu)*r*H0)/(8.0*M*mu);                          //ugyz
                Ug[9*ind+6]=(-gx*gz*sqrt(ky*ky)*(M-mu)*r*H1)/(8.0*M*mu);                    //ugzx
                Ug[9*ind+7]=kySign*(gz*ky*(M-mu)*r*H0)/(8.0*M*mu);                          //ugzy
                Ug[9*ind+8]=-(i*((M+mu)*H0 + gz*gz*i*sqrt(ky*ky)*(-M+mu)*r*H1))/(8.0*M*mu); //ugzz
              }
              if (calcSg)
              {
                Sg[18*ind+0]=(i*sqr(ky)*x[ixRec]*(((M-mu)*x[ixRec]*x[ixRec]*H0)/r+(((2.0*M-mu)*x[ixRec]*x[ixRec]+mu*z[izRec]*z[izRec])*H1)/(i*sqrt(ky*ky)*r*r*r)))/(4.0*M);                              //  sxxx
                Sg[18*ind+1]=-(i*x[ixRec]*(sqr(ky)*(M-mu)*H0+(i*sqrt(sqr(ky))*(M-2.0*mu)*H1)/r))/(4.0*M);                                                                                                //  sxyy
                Sg[18*ind+2]=(i*x[ixRec]*((sqr(ky)*(M-mu)*z[izRec]*z[izRec]*H0)/r+(i*sqrt(sqr(ky))*(-2.0*M*z[izRec]*z[izRec]+mu*(x[ixRec]*x[ixRec]+3.0*z[izRec]*z[izRec]))*H1)/(r*r*r)))/(4.0*M);        //  sxzz
                Sg[18*ind+3]=kySign*(ky*(-mu*H0 + (i*sqrt(sqr(ky))*(M-mu)*x[ixRec]*x[ixRec]*H1)/r))/(4.0*M);                                                                                             //  sxxy
                Sg[18*ind+4]=kySign*(ky*i*sqrt(sqr(ky))*(M-mu)*x[ixRec]*z[izRec]*H1)/(4.0*M*r);                                                                                                          //  sxyz
                Sg[18*ind+5]=-(i*z[izRec]*((i*sqrt(sqr(ky))*mu*H1)/r+(sqr(ky)*(M-mu)*x[ixRec]*x[ixRec]*H2)/r/r))/(4.0*M);                                                                                //  sxzx
                Sg[18*ind+6]=kySign*(ky*(mu*H0+(i*sqrt(sqr(ky))*(M-mu)*x[ixRec]*x[ixRec]*H1)/r))/(4.0*M);                                                                                                //  syxx
                Sg[18*ind+7]=kySign*(ky*((-3.0*M+2.0*mu)*H0 + i*sqrt(sqr(ky))*(-M+mu)*r*H1))/(4.0*M);                                                                                                    //  syyy
                Sg[18*ind+8]=kySign*(ky*(mu*H0+(i*sqrt(sqr(ky))*(M-mu)*z[izRec]*z[izRec]*H1)/r))/(4.0*M);                                                                                                //  syzz
                Sg[18*ind+9]=-0.25*i*x[ixRec]*((sqr(ky)*(M-mu)*H0)/M+(i*sqrt(sqr(ky))*H1)/r);                                                                                                            //  syxy
                Sg[18*ind+10]=-0.25*i*z[izRec]*((sqr(ky)*(M-mu)*H0)/M+(i*sqrt(sqr(ky))*H1)/r);                                                                                                           //  syyz
                Sg[18*ind+11]=kySign*(ky*i*sqrt(sqr(ky))*(M-mu)*x[ixRec]*z[izRec]*H1)/(4.0*M*r);                                                                                                         //  syzx
                Sg[18*ind+12]=(i*z[izRec]*((sqr(ky)*(M-mu)*x[ixRec]*x[ixRec]*H0)/r + (i*sqrt(sqr(ky))*(-2.0*M*x[ixRec]*x[ixRec] + mu*(3.0*x[ixRec]*x[ixRec]+z[izRec]*z[izRec]))*H1)/r/r/r))/(4.0*M);     //  szxx
                Sg[18*ind+13]=-(i*z[izRec]*(sqr(ky)*(M-mu)*H0 + (i*sqrt(sqr(ky))*(M-2.0*mu)*H1)/r))/(4.0*M);                                                                                             //  szyy
                Sg[18*ind+14]=(i*sqr(ky)*z[izRec]*(((M-mu)*z[izRec]*z[izRec]*H0)/r/r + ((mu*x[ixRec]*x[ixRec]+(2.0*M-mu)*z[izRec]*z[izRec])*H1)/(i*sqrt(sqr(ky))*r*r*r)))/(4.0*M);                       //  szzz
                Sg[18*ind+15]=kySign*(ky*i*sqrt(sqr(ky))*(M-mu)*x[ixRec]*z[izRec]*H1)/(4.0*M*r);                                                                                                         //  szxy
                Sg[18*ind+16]=kySign*(ky*(-mu*H0+(i*sqrt(sqr(ky))*(M-mu)*z[izRec]*z[izRec]*H1)/r))/(4.0*M);                                                                                              //  szyz
                Sg[18*ind+17]=-(i*x[ixRec]*((i*sqrt(sqr(ky))*mu*H1)/r + (sqr(ky)*(M-mu)*z[izRec]*z[izRec]*H2)/r/r))/(4.0*M);                                                                             //  szzx
              }
            }
          }
          else  // Dynamic solution
          {
            // Hankel functions
            double fnu=0.0;
            int HankelM=2;
            int HankelN=4;
            int kode=1;
            double* const cyr = new(nothrow) double[4];
            if (cyr==0) throw("Out of memory.");
            double* const cyi = new(nothrow) double[4];
            if (cyi==0) throw("Out of memory.");
            int nz;
            int ierr;

            double zar=real(ka*r);
            double zai=imag(ka*r);
            zbesh_(&zar,&zai,&fnu,&kode,&HankelM,&HankelN,cyr,cyi,&nz,&ierr);
            complex<double> H0a=cyr[0]+i*cyi[0];
            complex<double> H1a=cyr[1]+i*cyi[1];
            complex<double> H2a=cyr[2]+i*cyi[2];
            complex<double> H3a=cyr[3]+i*cyi[3];
            double zbr=real(kb*r);
            double zbi=imag(kb*r);
            zbesh_(&zbr,&zbi,&fnu,&kode,&HankelM,&HankelN,cyr,cyi,&nz,&ierr);
            complex<double> H0b=cyr[0]+i*cyi[0];
            complex<double> H1b=cyr[1]+i*cyi[1];
            complex<double> H2b=cyr[2]+i*cyi[2];
            complex<double> H3b=cyr[3]+i*cyi[3];
            delete [] cyr;
            delete [] cyi;

            complex<double> A = 1.0/(4.0*i*rho*sqr(omega[iFreq]));
            complex<double> B0= H0b - H0a;
            complex<double> B1= kb*H1b - ka*H1a;
            complex<double> B2= kb*kb*H2b - ka*ka*H2a;
            complex<double> B3= kb*kb*kb*H3b - ka*ka*ka*H3a;

            if (calcUg)
            { 
              Ug[9*ind+0]=A*(sqr(ks)*H0b-1.0/r*B1+gx*gx*B2); //ugxx
              Ug[9*ind+1]=i*ky*gx*A*B1;                      //ugxy
              Ug[9*ind+2]=gx*gz*A*B2;                        //ugxz
              Ug[9*ind+3]=i*ky*gx*A*B1;                      //ugyx
              Ug[9*ind+4]=A*(sqr(ks)*H0b-ky*ky*B0);          //ugyy
              Ug[9*ind+5]=i*ky*gz*A*B1;                      //ugyz
              Ug[9*ind+6]=gx*gz*A*B2;                        //ugzx
              Ug[9*ind+7]=i*ky*gz*A*B1;                      //ugzy
              Ug[9*ind+8]=A*(sqr(ks)*H0b-1.0/r*B1+gz*gz*B2); //ugzz
            }
            if (calcSg)
            {
              complex <double> egxvol = gx*A*(-ks*ks*kb*H1b+ky*ky*B1+4.0/r*B2-B3);
              complex <double> egxxx  = gx*A*((2.0/r*B2-ks*ks*kb*H1b)+1.0/r*B2 -gx*gx*B3);
              complex <double> egxzz  = gx*A*(1.0/r*B2 -gz*gz*B3);
              complex <double> egxyy  = gx*ky*ky*A*B1;
              complex <double> egxzx  = A*((1.0/r*B2-1.0/2.0*ks*ks*kb*H1b)*gz - gx*gz*gx*B3);
              complex <double> egxxy  = i*ky*A*((1.0/r*B1-1.0/2.0*ks*ks*H0b) -gx*gx*B2);
              complex <double> egxyz  =-i*ky*A*gz*gx*B2;
              complex <double> egzvol = gz*A*(-ks*ks*kb*H1b+ky*ky*B1+4.0/r*B2-B3);
              complex <double> egzxx  = gz*A*( 1.0/r*B2-gx*gx*B3);
              complex <double> egzzz  = gz*A*((2.0/r*B2-ks*ks*kb*H1b)+1.0/r*B2-gz*gz*B3);
              complex <double> egzyy  = gz*ky*ky*A*B1;
              complex <double> egzzx  = A*((1.0/r*B2 -1.0/2.0*ks*ks*kb*H1b)*gx-gx*gz*gz*B3);
              complex <double> egzxy  =-i*ky*A*gx*gz*B2;
              complex <double> egzyz  = i*ky*A*((1.0/r*B1-0.5*ks*ks*H0b)-gz*gz*B2);
              complex <double> egyvol =i*ky*A*(-ks*ks*H0b+ky*ky*B0+2.0/r*B1-B2);
              complex <double> egyxx  =i*ky*A*(1.0/r*B1-gx*gx*B2);
              complex <double> egyzz  =i*ky*A*(1.0/r*B1-gz*gz*B2);
              complex <double> egyyy  =i*ky*A*(-ks*ks*H0b+ky*ky*B0);
              complex <double> egyzx  =-i*ky*gx*gz*A*B2;
              complex <double> egyxy  =gx*A*(-0.5*ks*ks*kb*H1b+ky*ky*B1);
              complex <double> egyyz  =gz*A*(-0.5*ks*ks*kb*H1b+ky*ky*B1);

              Sg[18*ind+0]=lambda*egxvol+2.0*mu*egxxx;             //  sxxx
              Sg[18*ind+1]=lambda*egxvol+2.0*mu*egxyy;             //  sxyy
              Sg[18*ind+2]=lambda*egxvol+2.0*mu*egxzz;             //  sxzz
              Sg[18*ind+3]=2.0*mu*egxxy;                           //  sxxy
              Sg[18*ind+4]=2.0*mu*egxyz;                           //  sxyz
              Sg[18*ind+5]=2.0*mu*egxzx;                           //  sxzx
              Sg[18*ind+6]=lambda*egyvol+2.0*mu*egyxx;             //  syxx
              Sg[18*ind+7]=lambda*egyvol+2.0*mu*egyyy;             //  syyy
              Sg[18*ind+8]=lambda*egyvol+2.0*mu*egyzz;             //  syzz
              Sg[18*ind+9]=2.0*mu*egyxy;                           //  syxy
              Sg[18*ind+10]=2.0*mu*egyyz;                          //  syyz
              Sg[18*ind+11]=2.0*mu*egyzx;                          //  syzx
              Sg[18*ind+12]=lambda*egzvol+2.0*mu*egzxx;            //  szxx
              Sg[18*ind+13]=lambda*egzvol+2.0*mu*egzyy;            //  szyy
              Sg[18*ind+14]=lambda*egzvol+2.0*mu*egzzz;            //  szzz
              Sg[18*ind+15]=2.0*mu*egzxy;                          //  szxy
              Sg[18*ind+16]=2.0*mu*egzyz;                          //  szyz
              Sg[18*ind+17]=2.0*mu*egzzx;                          //  szzx
            }
          }
        }
      }
    }
  }
}
