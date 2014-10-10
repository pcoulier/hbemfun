/*BEMTIMECONV   Time convolution of Green's function.
 *
 *   u = BEMTIMECONV(t,ug,tBem,delt,type) computes the convolution of the 
 *   Green's function ug and the boundary element shape function psi. The 
 *   convolution is performed over the last dimension of ugt. The convolution
 *   is defined as:
 *
 *       /+inf
 *   u = |       ug(tau) * psi(tBem-tau) d tau
 *       /-inf
 *
 *   t     Time sampling of the Green's function (1 * nTime) or (nTime * 1).
 *         The vector should be monotonically increasing.
 *   ug    Green's functions (... * nTime). This is usually a Green's
 *         displacement or traction. The convolution is performed over
 *         the last dimension of ug. The function is assumed to be zero
 *         outside the sampling interval.
 *   tBem  Time sampling (1 * nTimeBem) or (nTimeBem * 1) for which the
 *         convolution is evaluated.
 *   delt  Time step that is used in the definition of the shape function
 *   type  Shape function type.
 *         1: Constant shape function from -delt to 0. Used for displacements.
 *         2: Triangular shape function from -delt to delt. Used for tractions
 *         3: Modified triangular shape function from -delt to 0. Used for 
 *         tractions
 *   u     Convolution (... * nTimeBem) evaluated at times tBem.
 */

/* $Make: mex -O -output bemtimeconv bemtimeconv_mex.cpp search1.cpp$*/

#include "mex.h"
#include <complex>
#include "search1.h"
#include "checklicense.h"

using namespace std;
void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{
  try
  {
    checklicense();
    
    if (nrhs<5) throw("Not enough input arguments.");
    if (nrhs>5) throw("Too many input arguments.");
    if (nlhs>1) throw("Too many output arguments.");

    if (!mxIsNumeric(prhs[0])) throw("Input argument 't' must be numeric.");
    if (mxIsSparse(prhs[0])) throw("Input argument 't' must not be sparse.");
    if (mxIsComplex(prhs[0])) throw("Input argument 't' must be real.");
    const double* const t=mxGetPr(prhs[0]);
    const unsigned int nTime=mxGetNumberOfElements(prhs[0]);
    for (unsigned int it=1; it<nTime; it++) if (!(t[it-1]<t[it])) throw("Input argument 't' must be monotonically increasing.");
    const double tMin=t[0];
    const double tMax=t[nTime-1];

    if (!mxIsNumeric(prhs[1])) throw("Input argument 'ug' must be numeric.");
    if (mxIsSparse(prhs[1])) throw("Input argument 'ug' must not be sparse.");
    if (mxIsComplex(prhs[1])) throw("Input argument 'ug' must be real.");
    const double* const ugt=mxGetPr(prhs[1]);
    const unsigned int nugDims=mxGetNumberOfDimensions(prhs[1]);
    const size_t* const ugDim=mxGetDimensions(prhs[1]);
    unsigned int nugComp=1;
    if ((nugDims==2)&&(ugDim[1]==1))
    {
      if (!(mxGetNumberOfElements(prhs[1])==mxGetNumberOfElements(prhs[0]))) throw("Number of elements in 'ug' must be equal to number of elements in t.");
    }
    else if (!(ugDim[nugDims-1]==(unsigned)nTime)) throw("Last dimension of 'ug' must be equal to number of elements in t.");
    else for (unsigned int iDim=0; iDim<nugDims-1; iDim++) nugComp=nugComp*ugDim[iDim];

    if (!mxIsNumeric(prhs[2])) throw("Input argument 'tBem' must be numeric.");
    if (mxIsSparse(prhs[2])) throw("Input argument 'tBem' must not be sparse.");
    if (mxIsComplex(prhs[2])) throw("Input argument 'tBem' must be real.");
    const double* const tBem=mxGetPr(prhs[2]);
    const unsigned int nTimeBem=mxGetNumberOfElements(prhs[2]);

    if (!mxIsNumeric(prhs[3])) throw("Input argument 'delt' must be numeric.");
    if (mxIsSparse(prhs[3])) throw("Input argument 'delt' must not be sparse.");
    if (mxIsComplex(prhs[3])) throw("Input argument 'delt' must be real.");
    const double* const deltutil=mxGetPr(prhs[3]);
    const double delt = deltutil[0];
    if (delt==0) throw("Input argument 'delt' must be zero.");

    if (!mxIsNumeric(prhs[4])) throw("Input argument 'type' must be numeric.");
    if (mxIsSparse(prhs[4])) throw("Input argument 'type' must not be sparse.");
    if (mxIsComplex(prhs[4])) throw("Input argument 'type' must be real.");
    const double* const typeutil=mxGetPr(prhs[4]);
    const unsigned int type = (unsigned int)(typeutil[0]);

    // OUTPUT ARGUMENTS
    size_t* const outDims = new(nothrow) size_t[nugDims];
    if (outDims==0) throw("Out of memory.");
    if ((nugDims==2)&&(ugDim[1]==1))
    {
      outDims[0]=1;
      outDims[1]=nTimeBem;
    }
    else
    {
      for (unsigned int iDim=0; iDim<nugDims-1; iDim++) outDims[iDim]=ugDim[iDim];
      outDims[nugDims-1]=nTimeBem;
    }

    plhs[0]=mxCreateNumericArray(nugDims,outDims,mxDOUBLE_CLASS,mxREAL);
    double* const u=mxGetPr(plhs[0]);

    // Initialise interpolation parameters
    unsigned int t00=0;
    unsigned int t01=0;
    unsigned int t10=0;
    unsigned int t11=0;
    unsigned int t20=0;
    unsigned int t21=0;
    double* const ival0 = new(nothrow) double[2];
    double* const ival1 = new(nothrow) double[2];
    double* const ival2 = new(nothrow) double[2];

    for (unsigned int iTime=0; iTime<nTimeBem; iTime++)
    {
      double t0=tBem[iTime]-delt;
      double t1=tBem[iTime];
      double t2=tBem[iTime]+delt;
      t0=max(t0,tMin);
      t1=max(t1,tMin);
      t2=max(t2,tMin);
      t0=min(t0,tMax);
      t1=min(t1,tMax);
      t2=min(t2,tMax);

      bool extrapFlag=false;
      search1(t0,t,nTime-1,t00,t01,ival0,extrapFlag);
      search1(t1,t,nTime-1,t10,t11,ival1,extrapFlag);
      search1(t2,t,nTime-1,t20,t21,ival2,extrapFlag);

      if (type==1)  // Constant shape function
      {
      unsigned int InnerPoints= t20-t10;
      for (unsigned int iComp=0; iComp<nugComp; iComp++)
      {
          unsigned int iPos=nugComp*iTime+iComp;
          u[iPos]= 0.5*(t[t11]-t1)*(ival1[0]*ugt[nugComp*t10+iComp]+ival1[1]*ugt[nugComp*t11+iComp]);
          u[iPos]+=0.5*(t[t11]-t1)*ugt[nugComp*t11+iComp];
          for (unsigned int tInt=0; tInt<InnerPoints-1; tInt++)
          {
            u[iPos]+=0.5*(t[t11+tInt+1]-t[t11+tInt+0])*ugt[nugComp*(t11+tInt+0)+iComp];
            u[iPos]+=0.5*(t[t11+tInt+1]-t[t11+tInt+0])*ugt[nugComp*(t11+tInt+1)+iComp];
          }
          u[iPos]+=0.5*(t2-t[t20])* ugt[nugComp*t20+iComp];
          u[iPos]+=0.5*(t2-t[t20])*(ival2[0]*ugt[nugComp*t20+iComp]+ival2[1]*ugt[nugComp*t21+iComp]);
        }
      }
      else if (type==2)
      {
        unsigned int ntInt=t10-t01;
        for (unsigned int iComp=0; iComp<nugComp; iComp++)
        {
          unsigned int iPos=nugComp*iTime+iComp;
          u[iPos]+=0.5*(t[t01]-t0)*((t[t01]-t0)/delt)*ugt[nugComp*t01+iComp];
          for (unsigned int tInt=0; tInt<ntInt; tInt++)
          {
            u[iPos]+=0.5*(t[t01+tInt+1]-t[t01+tInt+0])*((t[t01+tInt+0]-t0)/delt)*ugt[nugComp*(t01+tInt+0)+iComp];
            u[iPos]+=0.5*(t[t01+tInt+1]-t[t01+tInt+0])*((t[t01+tInt+1]-t0)/delt)*ugt[nugComp*(t01+tInt+1)+iComp];
          }
          u[iPos]+=0.5*(t1-t[t10])*((t[t10]-t0)/delt)*ugt[nugComp*t10+iComp];
          u[iPos]+=0.5*(t1-t[t10])*1.0*(ival1[0]*ugt[nugComp*t10+iComp]+ival1[1]*ugt[nugComp*t11+iComp]);
        }
        
        ntInt=t20-t11;
        for (unsigned int iComp=0; iComp<nugComp; iComp++)
        {
          //tau=(1-(t-t1)/delt)
          unsigned int iPos=nugComp*iTime+iComp;
          u[iPos]+= 0.5*(t[t11]-t1)*(1.0-0.0)*(ival1[0]*ugt[nugComp*t10+iComp]+ival1[1]*ugt[nugComp*t11+iComp]);
          u[iPos]+=0.5*(t[t11]-t1)*(1.0-(t[t11]-t1)/delt)*ugt[nugComp*t11+iComp];
          for (unsigned int tInt=0; tInt<ntInt; tInt++)
          {
            u[iPos]+=0.5*(t[t11+tInt+1]-t[t11+tInt+0])*(1.0-(t[t11+tInt+0]-t1)/delt)*ugt[nugComp*(t11+tInt+0)+iComp];
            u[iPos]+=0.5*(t[t11+tInt+1]-t[t11+tInt+0])*(1.0-(t[t11+tInt+1]-t1)/delt)*ugt[nugComp*(t11+tInt+1)+iComp];
          }
          u[iPos]+=0.5*(t2-t[t20])*(1.0-(t[t20]-t1)/delt)*ugt[nugComp*t20+iComp];
        }
      }
      else if (type==3) // Modified shape function
      {
        unsigned int InnerPoints= t20-t10;
        for (unsigned int iComp=0; iComp<nugComp; iComp++)
        {
          unsigned int iPos=nugComp*iTime+iComp;
          u[iPos]= 0.5*(t[t11]-t1)*(1.0-0.0)*(ival1[0]*ugt[nugComp*t10+iComp]+ival1[1]*ugt[nugComp*t11+iComp]);
          u[iPos]+=0.5*(t[t11]-t1)*(1.0-(t[t11]-t1)/delt)*ugt[nugComp*t11+iComp];
          for (unsigned int tInt=0; tInt<InnerPoints-1; tInt++)
          {
            u[iPos]+=0.5*(t[t11+tInt+1]-t[t11+tInt+0])*(1.0-(t[t11+tInt+0]-t1)/delt)*ugt[nugComp*(t11+tInt+0)+iComp];
            u[iPos]+=0.5*(t[t11+tInt+1]-t[t11+tInt+0])*(1.0-(t[t11+tInt+1]-t1)/delt)*ugt[nugComp*(t11+tInt+1)+iComp];
          }
          u[iPos]+=0.5*(t2-t[t20])*(1.0-(t[t20]-t1)/delt)*ugt[nugComp*t20+iComp];
        }
      }
      else throw("Unknown shape function type.");
    }
    delete [] outDims;
    delete [] ival0;
    delete [] ival1;
    delete [] ival2;
  }
  catch (const char* exception)
  {
    mexErrMsgTxt(exception);
  }
}
