#include <new>
#include <math.h>
using namespace std;

inline double sqr(const double& a)
{
  return a*a;
}

void shapefun(const unsigned int& ShapeType, const unsigned int& nXi,
              const double* const xi ,double* const N)
{
  if (ShapeType==1)       // Constant shape function
  {
  	for (unsigned int iXi=0; iXi<nXi; iXi++)
  	{
  	  N[iXi]=1.0;
  	}
  }
  else if  (ShapeType==2) // 3-node triangular shape function
  {
  	for (unsigned int iXi=0; iXi<nXi; iXi++)
  	{
      N[3*iXi+0]= 1.0-xi[iXi]-xi[nXi+iXi];
      N[3*iXi+1]= xi[iXi];
      N[3*iXi+2]= xi[nXi+iXi];
    }
  }
  else if  (ShapeType==3) // 6-node triangular shape function
  {
  	for (unsigned int iXi=0; iXi<nXi; iXi++)
  	{
      N[6*iXi+0]=(1.0-xi[iXi]-xi[nXi+iXi])*(1.0-2.0*xi[iXi]-2.0*xi[nXi+iXi]);
      N[6*iXi+1]=-xi[iXi]*(1.0-2.0*xi[iXi]);
      N[6*iXi+2]=-xi[nXi+iXi]*(1.0-2.0*xi[nXi+iXi]);
      N[6*iXi+3]=4.0*xi[iXi]*(1.0-xi[iXi]-xi[nXi+iXi]);
      N[6*iXi+4]=4.0*xi[iXi]*xi[nXi+iXi];
      N[6*iXi+5]=4.0*xi[nXi+iXi]*(1.0-xi[iXi]-xi[nXi+iXi]);
    }
  }
  else if  (ShapeType==4) // 4-node quadrilateral shape function
  {
  	for (unsigned int iXi=0; iXi<nXi; iXi++)
  	{
  	  N[4*iXi+0]=0.25*(1-xi[iXi])*(1-xi[nXi+iXi]);
      N[4*iXi+1]=0.25*(1+xi[iXi])*(1-xi[nXi+iXi]);
      N[4*iXi+2]=0.25*(1+xi[iXi])*(1+xi[nXi+iXi]);
      N[4*iXi+3]=0.25*(1-xi[iXi])*(1+xi[nXi+iXi]);
    }
  }
  else if  (ShapeType==5) // 8-node quadrilateral shape function
  {
  	for (unsigned int iXi=0; iXi<nXi; iXi++)
  	{
  	  N[8*iXi+0]=0.25*(1-xi[iXi])*(1-xi[nXi+iXi])*(-1-xi[iXi]-xi[nXi+iXi]);
      N[8*iXi+1]=0.25*(1+xi[iXi])*(1-xi[nXi+iXi])*(-1+xi[iXi]-xi[nXi+iXi]);
      N[8*iXi+2]=0.25*(1+xi[iXi])*(1+xi[nXi+iXi])*(-1+xi[iXi]+xi[nXi+iXi]);
      N[8*iXi+3]=0.25*(1-xi[iXi])*(1+xi[nXi+iXi])*(-1-xi[iXi]+xi[nXi+iXi]);
      N[8*iXi+4]=0.50*(1-xi[nXi+iXi])*(1-xi[iXi]*xi[iXi]);
      N[8*iXi+5]=0.50*(1+xi[iXi])*(1-xi[nXi+iXi]*xi[nXi+iXi]);
      N[8*iXi+6]=0.50*(1+xi[nXi+iXi])*(1-xi[iXi]*xi[iXi]);
      N[8*iXi+7]=0.50*(1-xi[iXi])*(1-xi[nXi+iXi]*xi[nXi+iXi]);
    }
  }
  else if  (ShapeType==6) // 9-node quadrilateral shape function
  {
  	for (unsigned int iXi=0; iXi<nXi; iXi++)
  	{
  	  N[9*iXi+0]=0.25*xi[iXi]*(xi[iXi]-1.0)*xi[nXi+iXi]*(xi[nXi+iXi]-1.0);
      N[9*iXi+1]=0.25*xi[iXi]*(xi[iXi]+1.0)*xi[nXi+iXi]*(xi[nXi+iXi]-1.0);
      N[9*iXi+2]=0.25*xi[iXi]*(xi[iXi]+1.0)*xi[nXi+iXi]*(xi[nXi+iXi]+1.0);
      N[9*iXi+3]=0.25*xi[iXi]*(xi[iXi]-1.0)*xi[nXi+iXi]*(xi[nXi+iXi]+1.0);
      N[9*iXi+4]=0.50*(1.0-xi[iXi]*xi[iXi])*xi[nXi+iXi]*(xi[nXi+iXi]-1.0);
      N[9*iXi+5]=0.50*xi[iXi]*(xi[iXi]+1.0)*(1.0-xi[nXi+iXi]*xi[nXi+iXi]);
      N[9*iXi+6]=0.50*(1.0-xi[iXi]*xi[iXi])*xi[nXi+iXi]*(xi[nXi+iXi]+1.0);
      N[9*iXi+7]=0.50*xi[iXi]*(xi[iXi]-1.0)*(1.0-xi[nXi+iXi]*xi[nXi+iXi]);
      N[9*iXi+8]=(1.0-xi[iXi]*xi[iXi])*(1.0-xi[nXi+iXi]*xi[nXi+iXi]);
    }
  }
  else if  (ShapeType==7) // 1D linear shape function
  {
  	for (unsigned int iXi=0; iXi<nXi; iXi++)
  	{
  	  N[2*iXi+0]=0.50*(1.0-xi[iXi]);
      N[2*iXi+1]=0.50*(1.0+xi[iXi]);
    }
  }
  else if  (ShapeType==8) // 1D quadratic shape function
  {
  	for (unsigned int iXi=0; iXi<nXi; iXi++)
  	{
  	  N[3*iXi+0]=0.50*xi[iXi]*(xi[iXi]-1.0);
      N[3*iXi+1]=1.0-sqr(xi[iXi]);
      N[3*iXi+2]=0.50*xi[iXi]*(xi[iXi]+1.0);
    }
  }
  else if  (ShapeType==9) // 1D cubic shape function
  {
  	for (unsigned int iXi=0; iXi<nXi; iXi++)
  	{
  	  N[4*iXi+0]=-0.0625*(1.0-xi[iXi])*(1.0-9.0*xi[iXi]*xi[iXi]);
      N[4*iXi+1]= 0.5625*(1.0-3.0*xi[iXi])*(1.0-xi[iXi]*xi[iXi]);
      N[4*iXi+2]= 0.5625*(1.0+3.0*xi[iXi])*(1.0-xi[iXi]*xi[iXi]);
      N[4*iXi+3]=-0.0625*(1.0+xi[iXi])*(1.0-9.0*xi[iXi]*xi[iXi]);
    }
  }
  else throw("Unknown shape function type.");
}



void shapederiv(const unsigned int& ShapeType, const unsigned int& nXi,
                const double* const xi ,double* const dN)
{
	if (ShapeType==1)       // Constant shape function
  {
  	for (unsigned int iXi=0; iXi<nXi; iXi++)
  	{
  	  dN[iXi]=0.0;
  	}
  }
	else if (ShapeType==2)  // 3-node triangular shape function
  {
  	for (unsigned int iXi=0; iXi<nXi; iXi++)
  	{
  	  dN[6*iXi+0]=-1.0;
      dN[6*iXi+1]= 1.0;
      dN[6*iXi+2]= 0.0;
      dN[6*iXi+3]=-1.0;
      dN[6*iXi+4]= 0.0;
      dN[6*iXi+5]= 1.0;
  	}
  }
	else if (ShapeType==3)  // 6-node triangular shape function
  {
  	for (unsigned int iXi=0; iXi<nXi; iXi++)
  	{
      dN[12*iXi+0] =-3.0+4.0*xi[iXi]+4.0*xi[nXi+iXi];
      dN[12*iXi+1] =-1.0+4.0*xi[iXi];
      dN[12*iXi+2] = 0.0;
      dN[12*iXi+3] = 4.0-8.0*xi[iXi]-4.0*xi[nXi+iXi];
      dN[12*iXi+4] = 4.0*xi[nXi+iXi];
      dN[12*iXi+5] =-4.0*xi[nXi+iXi];
      dN[12*iXi+6] =-3.0+4.0*xi[iXi]+4.0*xi[nXi+iXi];
      dN[12*iXi+7] = 0.0;
      dN[12*iXi+8] =-1.0+4.0*xi[nXi+iXi];
      dN[12*iXi+9] =-4.0*xi[iXi];
      dN[12*iXi+10]= 4.0*xi[iXi];
      dN[12*iXi+11]= 4.0-4.0*xi[iXi]-8.0*xi[nXi+iXi];
  	}
  }
	else if (ShapeType==4)  // 4-node quadrilateral shape function
  {
  	for (unsigned int iXi=0; iXi<nXi; iXi++)
  	{
      dN[8*iXi+0]=-0.25*(1.0-xi[nXi+iXi]);
      dN[8*iXi+1]= 0.25*(1.0-xi[nXi+iXi]);
      dN[8*iXi+2]= 0.25*(1.0+xi[nXi+iXi]);
      dN[8*iXi+3]=-0.25*(1.0+xi[nXi+iXi]);
      dN[8*iXi+4]=-0.25*(1.0-xi[iXi]);
      dN[8*iXi+5]=-0.25*(1.0+xi[iXi]);
      dN[8*iXi+6]= 0.25*(1.0+xi[iXi]);
      dN[8*iXi+7]= 0.25*(1.0-xi[iXi]);
  	}
  }
  else if (ShapeType==5)  // 8-node quadrilateral shape function
  {
  	for (unsigned int iXi=0; iXi<nXi; iXi++)
  	{
      dN[16*iXi+0]= 0.25*(1.0-xi[nXi+iXi])*(2.0*xi[iXi]+xi[nXi+iXi]);
      dN[16*iXi+1]= 0.25*(1.0-xi[nXi+iXi])*(2.0*xi[iXi]-xi[nXi+iXi]);
      dN[16*iXi+2]= 0.25*(1.0+xi[nXi+iXi])*(2.0*xi[iXi]+xi[nXi+iXi]);
      dN[16*iXi+3]= 0.25*(1.0+xi[nXi+iXi])*(2.0*xi[iXi]-xi[nXi+iXi]);
      dN[16*iXi+4]=-xi[iXi]*(1.0-xi[nXi+iXi]);
      dN[16*iXi+5]= 0.50*(1.0-xi[nXi+iXi]*xi[nXi+iXi]);
      dN[16*iXi+6]=-xi[iXi]*(1.0+xi[nXi+iXi]);
      dN[16*iXi+7]=-0.50*(1.0-xi[nXi+iXi]*xi[nXi+iXi]);
      dN[16*iXi+8]= 0.25*(1.0-xi[iXi])*( xi[iXi]+2.0*xi[nXi+iXi]);
      dN[16*iXi+9]= 0.25*(1.0+xi[iXi])*(-xi[iXi]+2.0*xi[nXi+iXi]);
      dN[16*iXi+10]= 0.25*(1.0+xi[iXi])*( xi[iXi]+2.0*xi[nXi+iXi]);
      dN[16*iXi+11]= 0.25*(1.0-xi[iXi])*(-xi[iXi]+2.0*xi[nXi+iXi]);
      dN[16*iXi+12]=-0.50*(1.0-xi[iXi]*xi[iXi]);
      dN[16*iXi+13]=-xi[nXi+iXi]*(1.0+xi[iXi]);
      dN[16*iXi+14]= 0.50*(1.0-xi[iXi]*xi[iXi]);
      dN[16*iXi+15]=-xi[nXi+iXi]*(1.0-xi[iXi]);
  	}
  }
  else if (ShapeType==6)  // 9-node quadrilateral shape function
  {
  	for (unsigned int iXi=0; iXi<nXi; iXi++)
  	{
      dN[18*iXi+0] = 0.25*xi[nXi+iXi]*(xi[nXi+iXi]-1.0)*(2.0*xi[iXi]-1.0);
      dN[18*iXi+1] = 0.25*xi[nXi+iXi]*(xi[nXi+iXi]-1.0)*(2.0*xi[iXi]+1.0);
      dN[18*iXi+2] = 0.25*xi[nXi+iXi]*(xi[nXi+iXi]+1.0)*(2.0*xi[iXi]+1.0);
      dN[18*iXi+3] = 0.25*xi[nXi+iXi]*(xi[nXi+iXi]+1.0)*(2.0*xi[iXi]-1.0);
      dN[18*iXi+4] =-xi[iXi]*xi[nXi+iXi]*(xi[nXi+iXi]-1.0);
      dN[18*iXi+5] =-0.50*(xi[nXi+iXi]-1.0)*(xi[nXi+iXi]+1.0)*(2.0*xi[iXi]+1.0);
      dN[18*iXi+6] =-xi[iXi]*xi[nXi+iXi]*(xi[nXi+iXi]+1.0);
      dN[18*iXi+7] =-0.50*(xi[nXi+iXi]-1.0)*(xi[nXi+iXi]+1.0)*(2.0*xi[iXi]-1.0);
      dN[18*iXi+8] = 2.0*xi[iXi]*(xi[nXi+iXi]-1.0)*(1.0+xi[nXi+iXi]);
      dN[18*iXi+9] = 0.25*xi[iXi]*(2.0*xi[nXi+iXi]-1.0)*(xi[iXi]-1.0);
      dN[18*iXi+10]= 0.25*xi[iXi]*(2.0*xi[nXi+iXi]-1.0)*(xi[iXi]+1.0);
      dN[18*iXi+11]= 0.25*xi[iXi]*(2.0*xi[nXi+iXi]+1.0)*(xi[iXi]+1.0);
      dN[18*iXi+12]= 0.25*xi[iXi]*(2.0*xi[nXi+iXi]+1.0)*(xi[iXi]-1.0);
      dN[18*iXi+13]=-0.50*(xi[iXi]-1.0)*(xi[iXi]+1.0)*(2.0*xi[nXi+iXi]-1.0);
      dN[18*iXi+14]=-xi[iXi]*(xi[iXi]+1.0)*xi[nXi+iXi];
      dN[18*iXi+15]=-0.50*(xi[iXi]-1.0)*(xi[iXi]+1.0)*(2.0*xi[nXi+iXi]+1.0);
      dN[18*iXi+16]=-xi[iXi]*(xi[iXi]-1.0)*xi[nXi+iXi];
      dN[18*iXi+17]= 2.0*(xi[iXi]-1.0)*(1.0+xi[iXi])*xi[nXi+iXi];
  	}
  }
  else if (ShapeType==7)  // 1D 2-node linear shape function
  {
  	for (unsigned int iXi=0; iXi<nXi; iXi++)
  	{
      dN[2*iXi+0] = -0.5;
      dN[2*iXi+1] =  0.5;
  	}
  }
  else if (ShapeType==8)  // 1D 3-node quadratic shape function
  {
  	for (unsigned int iXi=0; iXi<nXi; iXi++)
  	{
      dN[3*iXi+0] = xi[iXi]-0.5;
      dN[3*iXi+1] =-2.0*xi[iXi];
      dN[3*iXi+2] = xi[iXi]+0.5;
  	}
  }
  else if (ShapeType==9)  // 1D 4-node cubic shape function
  {
  	for (unsigned int iXi=0; iXi<nXi; iXi++)
  	{
      dN[4*iXi+0] =  0.0625-1.6875*sqr(xi[iXi])+1.125*xi[iXi];
      dN[4*iXi+1] = -1.6875+5.0625*sqr(xi[iXi])-1.125*xi[iXi];
      dN[4*iXi+2] =  1.6875-5.0625*sqr(xi[iXi])-1.125*xi[iXi];
      dN[4*iXi+3] = -0.0625+1.6875*sqr(xi[iXi])+1.125*xi[iXi];
  	}
  }
  else throw("Unknown shape function type.");
}

void shapenatcoord(const double* const dN, const unsigned int& nNod,
                   const unsigned int& nXi, const double* const NodCoord,
                   double* const nat, const unsigned int& EltDim)
/*
 * Compute the element natural basis vector.
 */
{
  for (unsigned int i=0; i< EltDim*3*nXi; i++) nat[i]=0.0;
  if (EltDim==2)      // 3D problem
  {
	  for (unsigned int iXi=0; iXi<nXi; iXi++)
	  {
	    for (unsigned int iNod=0; iNod<nNod; iNod++)
	    {
	      nat[6*iXi+0]+=dN[2*nNod*iXi+iNod]*NodCoord[0*nNod+iNod];
	      nat[6*iXi+1]+=dN[2*nNod*iXi+iNod]*NodCoord[1*nNod+iNod];
	      nat[6*iXi+2]+=dN[2*nNod*iXi+iNod]*NodCoord[2*nNod+iNod];
	      nat[6*iXi+3]+=dN[2*nNod*iXi+nNod+iNod]*NodCoord[0*nNod+iNod];
	      nat[6*iXi+4]+=dN[2*nNod*iXi+nNod+iNod]*NodCoord[1*nNod+iNod];
	      nat[6*iXi+5]+=dN[2*nNod*iXi+nNod+iNod]*NodCoord[2*nNod+iNod];
	    }
    }
  }
  else if (EltDim==1) // 2D problem
  {
    for (unsigned int iXi=0; iXi<nXi; iXi++)
  	{
  		for (unsigned int iNod=0; iNod<nNod; iNod++)
  		{
  			nat[2*iXi+0]+=dN[nNod*iXi+iNod]*NodCoord[0*nNod+iNod];
	      nat[2*iXi+1]+=dN[nNod*iXi+iNod]*NodCoord[2*nNod+iNod];
  	  }
  	}
  }
}

void jacobian(const double* const a, const unsigned int& nXi, double* const Jac,
                                                              const unsigned int& EltDim)
/*JACOBIAN Compute the element Jacobian.
 *  a         Element natural basis vector.
 *  nXi       Number of (Gaussian) points.
 *  EltDim    Element dimension 1 or 2 for 1D and 2D elements,
 *            respectively.
 *  Jac       Element Jacobian in the Gaussian points.
 */
{
	if (EltDim==2)
	{
	  for (unsigned int iXi=0; iXi<nXi; iXi++)
	  {
	     Jac[iXi]=sqrt((a[6*iXi+0]*a[6*iXi+0]+a[6*iXi+1]*a[6*iXi+1]+a[6*iXi+2]*a[6*iXi+2])*(a[6*iXi+3]*a[6*iXi+3]+a[6*iXi+4]*a[6*iXi+4]+a[6*iXi+5]*a[6*iXi+5])-sqr((a[6*iXi+0]*a[6*iXi+3]+a[6*iXi+1]*a[6*iXi+4]+a[6*iXi+2]*a[6*iXi+5])));
	  }
  }
  else if (EltDim==1)
  {
  	for (unsigned int iXi=0; iXi<nXi; iXi++)
	  {
	     Jac[iXi]=sqrt(a[2*iXi+0]*a[2*iXi+0]+a[2*iXi+1]*a[2*iXi+1]);
	  }
  }
}


void triangdiv(const double* const xiSing, const unsigned int& Parent, unsigned int& nDiv,
               double* const am, double* const a1, double* const a2,
               double* const rhom, double* const rho1, double* const rho2)
/*
 * Compute the element triangle division for singular integration.
 *
 *  angles are in the interval [-pi,pi]
 */
{
	double* const corner=new(nothrow) double[8];
	if (corner==0) throw("Out of memory. Type HELP MEMORY for your options.");

  if (Parent==1)      // Triangular parent element
	{
	  nDiv=3;
	  am[0]=-1.57079632679490;
	  am[1]= 0.78539816339745;
	  am[2]= 3.14159265358979;

	  corner[0]=0.0-xiSing[0];
	  corner[1]=1.0-xiSing[0];
	  corner[2]=0.0-xiSing[0];
	  corner[3]=0.0-xiSing[1];
	  corner[4]=0.0-xiSing[1];
	  corner[5]=1.0-xiSing[1];
  }
  else if (Parent==2) // Square parent element
  {
	  nDiv=4;
	  am[0]=-1.57079632679490;
	  am[1]= 0.00000000000000;
	  am[2]= 1.57079632679490;
	  am[3]= 3.14159265358979;

    corner[0]=-1.0-xiSing[0];
    corner[1]= 1.0-xiSing[0];
    corner[2]= 1.0-xiSing[0];
    corner[3]=-1.0-xiSing[0];
    corner[4]=-1.0-xiSing[1];
    corner[5]=-1.0-xiSing[1];
    corner[6]= 1.0-xiSing[1];
    corner[7]= 1.0-xiSing[1];
  }


	for (unsigned int iDiv=0; iDiv<nDiv; iDiv++)
	{
		rho1[iDiv]=sqrt(corner[iDiv]*corner[iDiv] + corner[nDiv+iDiv]*corner[nDiv+iDiv]);
		a1[iDiv]=atan2(corner[nDiv+iDiv],corner[iDiv]);
	}
	for (unsigned int iDiv=0; iDiv<nDiv-1; iDiv++)
	{
		a2[iDiv]=a1[iDiv+1];
		rho2[iDiv]=rho1[iDiv+1];
	}
	a2[nDiv-1]=a1[0];
	rho2[nDiv-1]=rho1[0];

	for (unsigned int iDiv=0; iDiv<nDiv; iDiv++)
	{
		   if ((am[iDiv]+1e-10)<a1[iDiv]) a1[iDiv]=a1[iDiv]-6.28318530717959;
	     if ((a2[iDiv]+1e-10)<am[iDiv]) a2[iDiv]=a2[iDiv]+6.28318530717959;
     rhom[iDiv]=rho1[iDiv]*cos(a1[iDiv]-am[iDiv]);
	}
	delete [] corner;
}
