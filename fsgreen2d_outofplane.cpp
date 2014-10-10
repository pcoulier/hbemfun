/* fsgreen2d_outofplane.cpp
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
void fsgreen2d_outofplane(const double Cs, const double Ds, const double rho,
                          const double* const x, const double* const z,
                          const double* const omega, const int nxRec, 
                          const int nzRec, const int nFreq, 
                          complex<double>* const Ug, complex<double>* const Sg,
                          const bool calcUg, const bool calcSg)
{
  for (int iFreq=0; iFreq<nFreq; iFreq++)
  {
    const complex<double> mu=rho*sqr(Cs)*(1.0+sign(omega[iFreq])*2.0*i*Ds);
    const double pi=3.141592653589793;

    for (int ixRec=0;ixRec<nxRec;ixRec++)
    {
      for (int izRec=0;izRec<nzRec;izRec++)
      {
        const double r=sqrt(sqr(x[ixRec])+sqr(z[izRec]));
        const double gx=x[ixRec]/r;
        const double gz=z[izRec]/r;
        int ind=(ixRec+nxRec*(izRec+nzRec*iFreq));
        if (omega[iFreq]==0)
        {
          if (calcUg)
          {
            Ug[0+ind]=-1.0/(2.0*pi*mu)*log(r);                          //ugyy
          }
          if (calcSg)
          {
            Sg[0+2*ind]=-gx/(2.0*pi*r);                                 //sgyxy
            Sg[1+2*ind]=-gz/(2.0*pi*r);                                 //sgyyz
          }
        }
        else  // Dynamic solution
        {
          const complex<double> Csc=sqrt(mu/rho);
          complex<double> ks=omega[iFreq]/Csc;
          if (imag(ks)>0) ks=-ks;
    
          // Hankel functions
          double fnu=0.0;
          int HankelM=2;
          int HankelN=2;
          int kode=1;
          double* const cyr = new(nothrow) double[2];
          if (cyr==0) throw("Out of memory.");
          double* const cyi = new(nothrow) double[2];
          if (cyi==0) throw("Out of memory.");
          int nz;
          int ierr;

          double zsr=real(ks*r);
          double zsi=imag(ks*r);
          zbesh_(&zsr,&zsi,&fnu,&kode,&HankelM,&HankelN,cyr,cyi,&nz,&ierr);
          complex<double> H0b=cyr[0]+i*cyi[0];
          complex<double> H1b=cyr[1]+i*cyi[1];
          delete [] cyr;
          delete [] cyi;

          int ind=(ixRec+nxRec*(izRec+nzRec*iFreq));
          if (calcUg)
          {
            Ug[0+ind]=1.0/(4.0*i*rho*sqr(omega[iFreq]))*(sqr(ks)*H0b); //ugyy
          }  
          if (calcSg)
          {
            complex<double> A = 1.0/(4.0*i*rho*sqr(omega[iFreq]));
            complex <double> Sutil=2.0*mu*A*(-0.5*ks*ks*ks*H1b);
            Sg[0+2*ind]=Sutil*gx;                                        //sgyxy
            Sg[1+2*ind]=Sutil*gz;                                        //sgyyz
          }
        }
      }
    }
  }
}
