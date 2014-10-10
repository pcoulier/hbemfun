#include <math.h>
#include <complex>
#include "mex.h"
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

void fsgreen3d(const double Cs, const double Cp,
               const double Ds, const double Dp, const double rho,
               const double* const r,
               const double* const z,
               const double* const omega, const int& nrRec,
               const int& nzRec, const int& nFreq,
               complex<double>* const Ug, complex<double>* const Sg,
               const bool calcUg, const bool calcSg)
{
  const double pi=3.141592653589793;
  for (int iFreq=0; iFreq<nFreq; iFreq++)
  {
    // CONSTANTS
    const complex<double> mu=rho*sqr(Cs)*(1.0+sign(omega[iFreq])*2.0*i*Ds);
    const complex<double> M=rho*sqr(Cp)*(1.0+sign(omega[iFreq])*2.0*i*Dp);
    const complex<double> nu=(M-2.0*mu)/(2.0*(M-mu));
    const complex<double> lambda=M-2.0*mu;
    const complex<double> Csc=sqrt(mu/rho);
    const complex<double> Cpc=sqrt(M/rho);

    for (int irRec=0;irRec<nrRec;irRec++)
    {
      for (int izRec=0;izRec<nzRec;izRec++)
      {
        const int ind=(irRec+nrRec*(izRec+nzRec*iFreq));

        if (omega[iFreq]==0.0) // STATIC GREEN'S FUNCTIONS
        {
          const double R=sqrt(sqr(r[irRec])+sqr(z[izRec]));
          const double rr= r[irRec]/R;
          const double rz= z[izRec]/R;
          if (calcUg)
          {
		 
			// mexPrintf("test...\n");
	
            const complex<double> facu=1.0/(16.0*pi*mu*(1.0-nu)*R);
            Ug[5*ind+0]=facu*(rr*rr+(3.0-4.0*nu));  // ugxr
            Ug[5*ind+1]=facu*(rr*rz);               // ugxz
            Ug[5*ind+2]=facu*(3.0-4.0*nu);          // ugyt
            Ug[5*ind+3]=facu*(rz*rr);               // ugzr
            Ug[5*ind+4]=facu*(rz*rz+(3.0-4.0*nu));  // ugzz
          }
          if (calcSg)
          {
            const complex<double> facs= -1.0/(8.0*pi*(1.0-nu))/sqr(R);
            Sg[10*ind+0]= facs*(3.0*rr*rr*rr+(1.0-2.0*nu)*rr);  // sgxrr
            Sg[10*ind+1]=-facs*(1.0-2.0*nu)*rr;                 // sgxtt
            Sg[10*ind+2]= facs*(3.0*rr*rz*rz-(1.0-2.0*nu)*rr);  // sgxzz
            Sg[10*ind+3]= facs*(3.0*rr*rr*rz+(1.0-2.0*nu)*rz);  // sgxzr
            Sg[10*ind+4]= facs*(1.0-2.0*nu)*rr;                 // sgyrt
            Sg[10*ind+5]= facs*(1.0-2.0*nu)*rz;                 // sgytz
            Sg[10*ind+6]= facs*(3.0*rz*rr*rr-(1.0-2.0*nu)*rz);  // sgzrr
            Sg[10*ind+7]=-facs*(1.0-2.0*nu)*rz;                 // sgztt
            Sg[10*ind+8]= facs*(3.0*rz*rz*rz+(1.0-2.0*nu)*rz);  // sgzzz
            Sg[10*ind+9]= facs*(3.0*rz*rz*rr+(1.0-2.0*nu)*rr);  // sgzzr
          }
        }
        else
        {
          const double R=sqrt(sqr(r[irRec])+sqr(z[izRec]));
          const complex<double> a=Csc/Cpc;
          const complex<double> gr= r[irRec]/R;
          const complex<double> gy= 0.0;
          const complex<double> gz= z[izRec]/R;
          const complex<double> Op=omega[iFreq]*R/Cpc;
          const complex<double> Os=omega[iFreq]*R/Csc;
          const complex<double> Psi=exp(-i*Op)*a*a*(i/Op+1.0/sqr(Op))+exp(-i*Os)*(1.0-i/Os-1.0/sqr(Os));
          const complex<double> Chi=exp(-i*Op)*a*a*(1.0-3.0*i/Op-3.0/sqr(Op))-exp(-i*Os)*(1.0-3.0*i/Os-3.0/sqr(Os));
          const complex<double> fac=1.0/(4.0*pi*mu*R);
          if (calcUg)
          {
            Ug[5*ind+0]= fac*(Psi+Chi*gr*gr);  // ugxr
            Ug[5*ind+1]= fac*Chi*gr*gz;        // ugxz
            Ug[5*ind+2]= fac*Psi;              // ugyt
            Ug[5*ind+3]= fac*Chi*gz*gr;        // ugzr
            Ug[5*ind+4]= fac*(Psi+Chi*gz*gz);  // ugzz
          }
          if (calcSg)
          {
            const complex<double> expS=exp(-i*omega[iFreq]*R/Csc);
            const complex<double> expP=exp(-i*omega[iFreq]*R/Cpc);
            const complex<double> DpsiDr=-i*expS*omega[iFreq]/Csc+2.0*i*Csc*expS/omega[iFreq]/sqr(R)-1.0/R*expS+sqr(Csc)/(R*R*R)*(2.0*expS/sqr(omega[iFreq])+expP*(-2.0/sqr(omega[iFreq])-2.0*i*R/omega[iFreq]/Cpc+sqr(R)/sqr(Cpc)));
            const complex<double> DchiDr=1.0/Csc/sqr(omega[iFreq])/(R*R*R)*(-6.0*i*sqr(Csc)*expS*omega[iFreq]*R+3.0*Csc*expS*sqr(omega[iFreq])*sqr(R)+i*expS*(omega[iFreq]*R)*(omega[iFreq]*R)*(omega[iFreq]*R)+Csc*Csc*Csc*(-6.0*expS+1.0/(Cpc*Cpc*Cpc)*expP*(6.0*Cpc*Cpc*Cpc+6.0*i*Cpc*Cpc*omega[iFreq]*R-3.0*Cpc*sqr(omega[iFreq]*R)-i*omega[iFreq]*R*omega[iFreq]*R*omega[iFreq]*R)));

            const complex<double> DurrDr0=fac*(gr*((DpsiDr-Psi/R)+(DchiDr-3.0*Chi/R)*gr*gr)+Chi/R*2.0*gr);
            const complex<double> DurrDz0=fac*(gz*((DpsiDr-Psi/R)+(DchiDr-3.0*Chi/R)*gr*gr));
            const complex<double> DuyyDr0=fac*(gr*((DpsiDr-Psi/R)+(DchiDr-3.0*Chi/R)*gy*gy));
            const complex<double> DuyyDz0=fac*(gz*((DpsiDr-Psi/R)+(DchiDr-3.0*Chi/R)*gy*gy));
            const complex<double> DuzzDr0=fac*(gr*((DpsiDr-Psi/R)+(DchiDr-3.0*Chi/R)*gz*gz));
            const complex<double> DuzzDz0=fac*(gz*((DpsiDr-Psi/R)+(DchiDr-3.0*Chi/R)*gz*gz)+Chi/R*2.0*gz);
            const complex<double> DuryDy0=fac*(gy*((DchiDr-3.0*Chi/R)*gr*gy)+Chi/R*gr);
            const complex<double> DuyzDy0=fac*(gy*((DchiDr-3.0*Chi/R)*gy*gz)+Chi/R*gz);
            const complex<double> DuzrDr0=fac*(gr*((DchiDr-3.0*Chi/R)*gz*gr)+Chi/R*gz);
            const complex<double> DuzrDz0=fac*(gz*((DchiDr-3.0*Chi/R)*gz*gr)+Chi/R*gr);

            const complex<double> exxx=DurrDr0;
            const complex<double> exyy=DuryDy0;
            const complex<double> exzz=DuzrDz0;
            const complex<double> exzx=0.5*(DuzrDr0+DurrDz0);
            const complex<double> eyxy=0.5*(DuryDy0+DuyyDr0);
            const complex<double> eyyz=0.5*(DuyyDz0+DuyzDy0);
            const complex<double> ezxx=DuzrDr0;
            const complex<double> ezyy=DuyzDy0;
            const complex<double> ezzz=DuzzDz0;
            const complex<double> ezzx=0.5*(DuzzDr0+DuzrDz0);
            const complex<double> exvol=exxx+exyy+exzz;
            const complex<double> ezvol=ezxx+ezyy+ezzz;

            Sg[10*ind+0]=lambda*exvol+2.0*mu*exxx;   // sgxrr
            Sg[10*ind+1]=lambda*exvol+2.0*mu*exyy;   // sgxtt
            Sg[10*ind+2]=lambda*exvol+2.0*mu*exzz;   // sgxzz
            Sg[10*ind+3]=2.0*mu*exzx;                // sgxzr
            Sg[10*ind+4]=2.0*mu*eyxy;                // sgyrt
            Sg[10*ind+5]=2.0*mu*eyyz;                // sgytz
            Sg[10*ind+6]=lambda*ezvol+2.0*mu*ezxx;   // sgzrr
            Sg[10*ind+7]=lambda*ezvol+2.0*mu*ezyy;   // sgztt
            Sg[10*ind+8]=lambda*ezvol+2.0*mu*ezzz;   // sgzzz
            Sg[10*ind+9]=2.0*mu*ezzx;                // sgzzr
          }
        }
      }
    }
  }
}
