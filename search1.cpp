#include <math.h>
#include "mex.h"

void search1(const double xi, const double* const x, const unsigned int& indEnd,
             unsigned int& ind1, unsigned int& ind2, double* const interpval, bool& extrapFlag)
/* SEARCH1 performs a 1D table lookup.  XI is the ordinate of the search
 * value, while X is a strictly monotonically increasing vector with
 * indEnd+1 entries. The function returns the indices ind1 and ind2 where XI is
 * located in the table X, together with the linear interpolation
 * coefficients INTERPVAL[0] and INTERPVAL[1] as:
 *         XI=INTERPVAL[0]*X[ind1]+INTERPVAL[1]*X[ind2]
 *
 * A good guess for the values ind1 increases performance.
 */
{
	// limit initial guess
	if (ind1>indEnd-1) ind1=indEnd-1;
	if (indEnd==0)              // only one element in x: ind1=ind2=0
  {
  	ind1=0;
  	ind2=0;
  	interpval[0]=1.0;
  	interpval[1]=0.0;
  }
  else if (xi<=x[0])         // lower extrapolation
  {
  	ind1=0;
  	ind2=1;
  	interpval[0]=(x[ind2]-xi)/(x[ind2]-x[ind1]);
  	interpval[1]=1.0-interpval[0];
  }
  else if (xi>=x[indEnd])    // upper extrapolation
  {
  	ind1=indEnd;
  	ind2=indEnd;
  	interpval[0]=0.0;
  	interpval[1]=1.0;
  }
  else if ((x[ind1]<=xi) && (xi<=x[ind1+1]))  // initial guess correct
  {
  	ind2=ind1+1;
    interpval[0]=(x[ind2]-xi)/(x[ind2]-x[ind1]);
    interpval[1]=1.0-interpval[0];
  }
  else if (xi<x[ind1])        // downward search
	{
  	while (xi<x[ind1]) ind1--;
	  ind2=ind1+1;
    interpval[0]=(x[ind2]-xi)/(x[ind2]-x[ind1]);
    interpval[1]=1.0-interpval[0];
  }
  else if (x[ind1+1]<xi )        // upward search
	{
  	while (x[ind1+1]<xi) ind1++;
	  ind2=ind1+1;
    interpval[0]=(x[ind2]-xi)/(x[ind2]-x[ind1]);
    interpval[1]=1.0-interpval[0];
  }
  else throw("Unexpected error in search1");
  
  if ((x[0]-xi>0)) extrapFlag=true;
  else if ((xi-x[indEnd]>0)) extrapFlag=true;
  else extrapFlag=false;
//     mexPrintf("x[0]: %f \n",x[0]);
//     mexPrintf("x[indEnd]: %f \n",x[indEnd]);
//     mexPrintf("xi: %f \n",xi);
}

void searchClosest(const double xi, const double* const x, const unsigned int& xend,
                                    unsigned int& xind)
/* SEARCHCLOSEST performs a 1D table lookup, where the closest value is
 * returned.  XI is the ordinate of the search value, while X is a 
 * strictly monotonically increasing vector with XEND entries. The 
 * function returns the index XIND of the point in X closest to XI.
 *
 * A good guess for the value XIND increases performance.
 *
 */                                    
{
	// limit initial guess
	if (xind > xend) xind=xend;
	
	// lookup closest point
  if (xend==0)             xind=0;            // only one element in x: xind=0
  else if (xi <= x[0])     xind=0;            // lower extrapolation
  else if (xi >= x[xend])  xind=xend;         // upper extrapolation
  else if (xi == x[xind])  xind=xind;         // initial guess correct
  else if (xi <  x[xind])                      // downward search
  {
  	while (xi < x[xind]) xind--;
  	if (((xi-x[xind])/(x[xind+1]-x[xind]))>0.50) xind+=1;
  }   
  else if (x[xind] < xi)                     // upward search
  {
  	while (x[xind] < xi) xind++;
  	if (((x[xind]-xi)/(x[xind]-x[xind-1]))>0.50) xind-=1;
  }
  else throw("unexpected error in searchClosest");
}
