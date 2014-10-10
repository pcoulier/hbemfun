#include <math.h>
#include <complex>
using namespace std;


/******************************************************************************/
inline double sqr(const double& a)
{
  return a*a;
}

/******************************************************************************/
void fsgreen3dt(const double& Cs, const double& Cp, const double& rho,
                const int& ftyp, const double& delt, const double* const r, 
                const double* const z, const double* const t, const int& nrRec,
                const int& nzRec, const int& nTime, double* const Ug, 
                double* const Sg, const bool calcUg, const bool calcSg)
{
  const double pi=3.141592653589793;
  const double mu=rho*Cs*Cs;
  const double beta=Cs/Cp;
  
  for (int iTime=0; iTime<nTime; iTime++){
    for (int irRec=0; irRec<nrRec; irRec++){
      for (int izRec=0; izRec<nzRec; izRec++){

        const int ind=(irRec+nrRec*(izRec+nzRec*iTime));
        const double R=sqrt(r[irRec]*r[irRec]+z[izRec]*z[izRec]);

        // Retarded time
        const double trp=t[iTime]-R/Cp;
        const double trs=t[iTime]-R/Cs;
        double fs;
        double fsd;
        double fp;
        double fpd;
        double A;
        if (ftyp==0){
          fs=0.0;
          fsd=0.0;
          fp=0.0;
          fpd=0.0;
          if ( (t[iTime]<(R/Cs)) && (t[iTime]>(R/Cp)) ){
            A=t[iTime]*Cs/R;
          }
          else{
            A=0.0;
          }
        }
        else if (ftyp==1){
          fs=(double)(trs>(-delt))-(double)(trs>0.0);
          fsd=0.0;
          fp=(double)(trp>(-delt))-(double)(trp>0.0);
          fpd=0.0;
          const double ll=max(-delt,t[iTime]-R/Cs);
          const double ul=min(0.0,t[iTime]-beta*R/Cs);
          if (ll<ul){
            A=sqr(Cs/R)*(0.5*(sqr(ul)-sqr(ll))-t[iTime]*(ul-ll)); 
          }
          else{
            A=0.0;
          }
        }
        else if (ftyp==2){
          fs= (1.0+trs/delt)*((trs>-delt)-(trs>0.0))+(1.0-trs/delt)*((trs>0.0)-(trs>delt));
          fsd=delt*((trs>-delt)-(trs>0.0))-1.0/delt*((trs>0.0)-(trs>delt));
          fp= (1.0+trp/delt)*((trp>-delt)-(trp>0.0))+(1.0-trp/delt)*((trp>0.0)-(trp>delt));
          fpd=delt*((trp>-delt)-(trp>0.0))-1.0/delt*((trp>0.0)-(trp>delt));
          const double ll=max(-delt,t[iTime]-R/Cs);        
          const double ul=min(delt,t[iTime]-beta*R/Cs);
          if (ll<ul){
            if (ul<-delt) A=0.0;
            else if (ul<0.0){
              A=sqr(Cs/R)*(-t[iTime]*(ul-ll)+0.5*(1.0-t[iTime]/delt)*(sqr(ul)-sqr(ll))+1.0/3.0/delt*(ul*ul*ul-ll*ll*ll));
            }
            else if (ll<0.0){
              A=sqr(Cs/R)*(-t[iTime]*(-ll)+0.5*(1.0-t[iTime]/delt)*(-sqr(ll))+1.0/3.0/delt*(-ll*ll*ll))+sqr(Cs/R)*(t[iTime]*(ul)-0.5*(1.0+t[iTime]/delt)*( sqr(ul))+1.0/3.0/delt*( ul*ul*ul));
            }
            else{
              A=sqr(Cs/R)*( t[iTime]*(ul-ll)-0.5*(1.0+t[iTime]/delt)*(sqr(ul)-sqr(ll))+1.0/3.0/delt*(ul*ul*ul-ll*ll*ll));
            }
          }
          else{
            A=0.0;
          }
        }
        else throw("undefined function type ftyp.");
        
        const double B=-3.0*A-(fs-beta*beta*fp);
        const double Cpar=-5.0*B+(fs+Cs*fsd)-beta*beta*(fp+Cp*fpd);
        const double Fs=2.0*B-(fs+Cs*fsd);
        const double Fp=2.0*B-(1.0-2.0*beta*beta)*(fp+Cp*fpd);

        const double rr=r[irRec]/R;
        const double rz=z[izRec]/R;
  
        // DISPLACEMENTS
        if (calcUg){
        const double facu=1.0/(4.0*pi*mu*R);
          Ug[5*ind+0]=facu*(B*rr*rr+(A+fs));                        // ugxr
          Ug[5*ind+1]=facu*B*rr*rz;                                 // ugxz
          Ug[5*ind+2]=facu*(A+fs);                                  // ugyt
          Ug[5*ind+3]=facu*B*rz*rr;                                 // ugzr
          Ug[5*ind+4]=facu*(B*rz*rz+(A+fs));                        // ugzz
        }

        // STRESSES
        if (calcSg){
          const double facs=1.0/(4.0*pi*R*R);
          Sg[10*ind+0]=facs*(2.0*Cpar*rr*rr*rr+rr*Fp+2.0*rr*Fs);  // sgxrr
          Sg[10*ind+1]=facs*(rr*Fp);                              // sgxtt
          Sg[10*ind+2]=facs*(2.0*Cpar*rz*rz*rr+rr*Fp);            // sgxzz
          Sg[10*ind+3]=facs*(2.0*Cpar*rz*rr*rr+rz*Fs);            // sgxzr
          Sg[10*ind+4]=facs*(rr*Fs);                              // sgyrt
          Sg[10*ind+5]=facs*(rz*Fs);                              // sgytz
          Sg[10*ind+6]=facs*(2.0*Cpar*rr*rr*rz+rz*Fp);            // sgzrr
          Sg[10*ind+7]=facs*(rz*Fp);                              // sgztt
          Sg[10*ind+8]=facs*(2.0*Cpar*rz*rz*rz+rz*Fp+2.0*rz*Fs);  // sgzzz
          Sg[10*ind+9]=facs*(2.0*Cpar*rz*rr*rz+rr*Fs);            // sgzzr
        }
      }
    }
  }
}
