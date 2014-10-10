#include <complex>
#include <math.h>
using namespace std;

void bemnormal(const double* const a, const unsigned int& nXi, const unsigned int& EltDim, double* const normal)
{
  if (EltDim==2)
  {
    for (unsigned int iXi=0; iXi<nXi; iXi++)
    {
      double nx=a[6*iXi+1]*a[6*iXi+5]-a[6*iXi+2]*a[6*iXi+4];
      double ny=a[6*iXi+2]*a[6*iXi+3]-a[6*iXi+0]*a[6*iXi+5];
      double nz=a[6*iXi+0]*a[6*iXi+4]-a[6*iXi+1]*a[6*iXi+3];
      double nLen=sqrt(nx*nx+ny*ny+nz*nz);
      //if (nLen==0) nLen=1;
      normal[3*iXi+0]= nx/nLen;
      normal[3*iXi+1]= ny/nLen;
      normal[3*iXi+2]= nz/nLen;
    }
  }
  else if (EltDim==1)
  {
  	for (unsigned int iXi=0; iXi<nXi; iXi++)
    {
      double nx=-a[2*iXi+1];
      double nz= a[2*iXi+0];

      double nLen=sqrt(nx*nx+nz*nz);
      if (nLen==0) nLen=1;
      normal[3*iXi+0]= nx/nLen;
      normal[3*iXi+1]= 0.0;
      normal[3*iXi+2]= nz/nLen;
    }
  }
}
