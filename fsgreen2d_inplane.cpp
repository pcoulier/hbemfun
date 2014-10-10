/* fsgreen2d_inplane.cpp
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
void fsgreen2d_inplane(const double Cs, const double Cp,
                       const double Ds, const double Dp, const double rho,
                       const double* const x, const double* const z,
                       const double* const omega, const int nxRec,
                       const int nzRec, const int nFreq,
                       complex<double>* const Ug, complex<double>* const Sg,
                       const bool calcUg, const bool calcSg)
{
  for (int iFreq=0; iFreq<nFreq; iFreq++)
  {
    const complex<double> mu=rho*sqr(Cs)*(1.0+sign(omega[iFreq])*2.0*i*Ds);
    const complex<double> M=rho*sqr(Cp)*(1.0+sign(omega[iFreq])*2.0*i*Dp);
    const complex<double> lambda=M-2.0*mu;
    complex<double> nu=(M-2.0*mu)/(2.0*(M-mu));
    const double pi=3.141592653589793;

    for (int ixRec=0;ixRec<nxRec;ixRec++)
    {
      for (int izRec=0;izRec<nzRec;izRec++)
      {
        const double r=sqrt(sqr(x[ixRec])+sqr(z[izRec]));
        const double gx=x[ixRec]/r;
        const double gz=z[izRec]/r;
        int ind=(ixRec+nxRec*(izRec+nzRec*iFreq));
        
        // Static solution
        if (omega[iFreq]==0)
        {
          if (calcUg)
          {
            const double logr=log(r);
            complex<double> Au= 1.0/(8.0*pi*mu*(1.0-nu));
            complex<double> nu3=3.0-4.0*nu;
            Ug[4*ind+0]=Au*(gx*gx-nu3*logr);          //ugxx
            Ug[4*ind+1]=Au*(gz*gx);                   //ugxz
            Ug[4*ind+2]=Au*(gz*gx);                   //ugzx
            Ug[4*ind+3]=Au*(gz*gz-nu3*logr);          //ugzz
          }
          if (calcSg)
          {
            complex<double> As=-1.0/(4.0*pi*(1.0-nu)*r);
            complex<double> nu2=1.0-2.0*nu;
            Sg[0+6*ind]=As*(2.0*sqr(gx)*gx+nu2*gx);   //sgxxx
            Sg[1+6*ind]=As*(2.0*gx*sqr(gz)-nu2*gx);   //sgxzz
            Sg[2+6*ind]=As*(2.0*sqr(gx)*gz+nu2*gz);   //sgxzx
            Sg[3+6*ind]=As*(2.0*gz*sqr(gx)-nu2*gz);   //sgzxx
            Sg[4+6*ind]=As*(2.0*sqr(gz)*gz+nu2*gz);   //sgzzz
            Sg[5+6*ind]=As*(2.0*sqr(gz)*gx+nu2*gx);   //sgzzx
          }
        }
        // Dynamic solution
        else
        {
          const complex<double> Cpc=sqrt(M/rho);
          const complex<double> Csc=sqrt(mu/rho);
          complex<double> kp=omega[iFreq]/Cpc;
          complex<double> ks=omega[iFreq]/Csc;
          const complex<double> beta=Csc/Cpc;
          if (imag(kp)>0) kp=-kp;
          if (imag(ks)>0) ks=-ks;

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

          double zar=real(kp*r);
          double zai=imag(kp*r);
          zbesh_(&zar,&zai,&fnu,&kode,&HankelM,&HankelN,cyr,cyi,&nz,&ierr);
          complex<double> H0a=cyr[0]+i*cyi[0];
          complex<double> H1a=cyr[1]+i*cyi[1];
          complex<double> H2a=cyr[2]+i*cyi[2];
          complex<double> H3a=cyr[3]+i*cyi[3];
          double zbr=real(ks*r);
          double zbi=imag(ks*r);
          zbesh_(&zbr,&zbi,&fnu,&kode,&HankelM,&HankelN,cyr,cyi,&nz,&ierr);
          complex<double> H0b=cyr[0]+i*cyi[0];
          complex<double> H1b=cyr[1]+i*cyi[1];
          complex<double> H2b=cyr[2]+i*cyi[2];
          complex<double> H3b=cyr[3]+i*cyi[3];
          delete [] cyr;
          delete [] cyi;

          complex<double> A = 1.0/(4.0*i*rho*sqr(omega[iFreq]));
          complex<double> B2= ks*ks*H2b - kp*kp*H2a;
          if (calcUg)
          {
            complex<double> B1= ks*H1b - kp*H1a;
            Ug[4*ind+0]=A*(sqr(ks)*H0b-1.0/r*B1+gx*gx*B2);  //ugxx
            Ug[4*ind+1]=gx*gz*A*B2;                         //ugxz
            Ug[4*ind+2]=gx*gz*A*B2;                         //ugzx
            Ug[4*ind+3]=A*(sqr(ks)*H0b-1.0/r*B1+gz*gz*B2);  //ugzz
          }
          if (calcSg)
          {
            complex<double> B3= ks*ks*ks*H3b - kp*kp*kp*H3a;
            complex <double> egxvol = gx*A*(-ks*ks*ks*H1b+4.0/r*B2-B3);
            complex <double> egxxx  = gx*A*((2.0/r*B2-ks*ks*ks*H1b)+1.0/r*B2 -gx*gx*B3);
            complex <double> egxzz  = gx*A*(1.0/r*B2 -gz*gz*B3);
            complex <double> egxzx  = A*((1.0/r*B2-1.0/2.0*ks*ks*ks*H1b)*gz - gx*gz*gx*B3);
            complex <double> egzvol = gz*A*(-ks*ks*ks*H1b+4.0/r*B2-B3);
            complex <double> egzxx  = gz*A*( 1.0/r*B2-gx*gx*B3);
            complex <double> egzzz  = gz*A*((2.0/r*B2-ks*ks*ks*H1b)+1.0/r*B2-gz*gz*B3);
            complex <double> egzzx  = A*((1.0/r*B2 -1.0/2.0*ks*ks*ks*H1b)*gx-gx*gz*gz*B3);

            Sg[6*ind+0]=lambda*egxvol+2.0*mu*egxxx;         //sgxxx
            Sg[6*ind+1]=lambda*egxvol+2.0*mu*egxzz;         //sgxzz
            Sg[6*ind+2]=2.0*mu*egxzx;                       //sgxzx
            Sg[6*ind+3]=lambda*egzvol+2.0*mu*egzxx;         //sgzxx
            Sg[6*ind+4]=lambda*egzvol+2.0*mu*egzzz;         //sgzzz
            Sg[6*ind+5]=2.0*mu*egzzx;                       //sgzzx
          }
        }
      }
    }
  }
}
